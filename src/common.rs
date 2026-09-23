use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;

/// Get k-mer length from the first line of an index file that has a non-empty
/// k-mer. All alleles of a variant share the same k-mer length, so the first
/// non-empty k-mer of any allele is taken. Returns an error if no line in the
/// index yields k > 0.
pub fn k_from_index(index: &str) -> Result<i64, io::Error> {
    let file = File::open(index)?;
    let reader = BufReader::new(file);
    for line_result in reader.lines() {
        let line = line_result?;
        let split_line: Vec<&str> = line.split(',').collect();
        if split_line.len() < 8 {
            continue;
        }
        for kmers_by_allele in split_line[7].split('|') {
            for kmer in kmers_by_allele.split(';') {
                if !kmer.is_empty() {
                    return Ok(kmer.len() as i64);
                }
            }
        }
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        "no k-mer length > 0 found in index",
    ))
}

/// Read chromosome lengths from a `.fai` file.
pub fn read_fai(fasta_path: &str) -> Result<HashMap<String, i64>, Box<dyn std::error::Error>> {
    let mut fai_path = fasta_path.to_string();
    fai_path.push_str(".fai");
    let file = File::open(fai_path)?;
    let reader = BufReader::new(file);
    let mut chrom_lengths: HashMap<String, i64> = HashMap::new();
    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split('\t').collect();
        let chrom_name = fields[0].to_string();
        let chrom_length = fields[1].parse::<i64>();
        chrom_lengths.insert(chrom_name, chrom_length?);
    }
    Ok(chrom_lengths)
}

/// Uppercase sequence and convert non-ATGC characters to N.
pub fn stand_seq(input: &str) -> String {
    input
        .chars()
        .map(|c| match c {
            'A' | 'T' | 'C' | 'G' => c,
            'a' => 'A',
            't' => 'T',
            'c' => 'C',
            'g' => 'G',
            _ => 'N',
        })
        .collect()
}

/// Slide a window of k bases through a sequence, keeping the lexicographically
/// smaller of each k-mer and its reverse complement. K-mers containing
/// non-ATGC characters are skipped.
pub fn get_canonical_kmers(sequence: &str, k: usize) -> Vec<String> {
    let mut canonical_kmers = Vec::new();
    if k == 0 || k > sequence.len() {
        return canonical_kmers;
    }
    for i in 0..=(sequence.len() - k) {
        let kmer_slice = &sequence[i..i + k];
        if kmer_slice.chars().all(|c| "ATGC".contains(c)) {
            let reverse_complement = reverse_complement(kmer_slice);
            if kmer_slice < reverse_complement.as_str() {
                canonical_kmers.push(kmer_slice.to_string());
            } else {
                canonical_kmers.push(reverse_complement.to_string());
            }
        }
    }
    canonical_kmers
}

/// Compute the reverse complement of a DNA sequence.
pub fn reverse_complement(dna_sequence: &str) -> String {
    dna_sequence
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => base,
        })
        .collect()
}

/// Pack a canonical ATGC k-mer into an integer (2 bits per base) for cheap
/// hashing. Returns `None` if the k-mer is empty (the empty-allele
/// pseudo-entry ""), contains non-ATGC characters, or is longer than 31 bases.
pub fn pack_kmer(kmer: &str) -> Option<u64> {
    if kmer.is_empty() || kmer.len() > 31 {
        return None;
    }
    let mut packed: u64 = 0;
    for c in kmer.bytes() {
        let two_bit = match c {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        packed = (packed << 2) | two_bit;
    }
    Some(packed)
}

/// Unpack a 2-bit-coded k-mer of length `k` back into its ATGC string.
pub fn unpack_kmer(packed: u64, k: usize) -> String {
    const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    let mut kmer = String::with_capacity(k);
    for i in (0..2 * k).step_by(2).rev() {
        let two_bit = ((packed >> i) & 0b11) as usize;
        kmer.push(BASES[two_bit]);
    }
    kmer
}

/// Key for hashing k-mers: packed `u64` for k <= 31, heap `String` fallback
/// for k > 31 or non-ATGC k-mers (e.g. the empty-allele pseudo-entry "").
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum KmerKey {
    Packed(u64),
    Str(String),
}

impl KmerKey {
    pub fn from_kmer(kmer: &str) -> KmerKey {
        match pack_kmer(kmer) {
            Some(packed) => KmerKey::Packed(packed),
            None => KmerKey::Str(kmer.to_string()),
        }
    }

    /// Convert back to the k-mer string, given k for packed keys.
    pub fn to_kmer(&self, k: usize) -> String {
        match self {
            KmerKey::Packed(packed) => unpack_kmer(*packed, k),
            KmerKey::Str(kmer) => kmer.clone(),
        }
    }
}

/// Slide a window of k bases through a 2-bit-coded sequence, keeping the
/// packed form of the lexicographically smaller of each k-mer and its reverse
/// complement. K-mers containing non-ATGC characters are skipped. The
/// sequence must be uppercase ATGC/N (as produced by `stand_seq`); non-ATGC
/// characters (e.g. N) are encoded as 4 so they can be detected and skipped.
pub fn get_canonical_kmers_packed(sequence: &str, k: usize) -> Vec<u64> {
    let mut packed_kmers = Vec::new();
    if k == 0 || k > 31 || k > sequence.len() {
        return packed_kmers;
    }
    // Encode sequence as 2-bit codes, with 4 marking non-ATGC characters.
    let codes: Vec<u8> = sequence
        .bytes()
        .map(|c| match c {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        })
        .collect();
    let kmer_mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let mut forward: u64 = 0;
    let mut reverse: u64 = 0;
    let mut run_without_bad_char = 0;
    for i in 0..codes.len() {
        let code = codes[i];
        if code == 4 {
            forward = 0;
            reverse = 0;
            run_without_bad_char = 0;
            continue;
        }
        forward = ((forward << 2) | code as u64) & kmer_mask;
        reverse = (reverse >> 2) | ((3 - code as u64) << (2 * (k - 1)));
        run_without_bad_char += 1;
        if run_without_bad_char == k {
            // Reverse was computed in lock-step with forward, so it holds the
            // reverse complement of the k-mer ending at i.
            if forward < reverse {
                packed_kmers.push(forward);
            } else {
                packed_kmers.push(reverse);
            }
            run_without_bad_char -= 1;
        }
    }
    packed_kmers
}

/// Read one field of an index file, splitting on `|` within each line.
pub fn read_index_field(index: &str, column: usize) -> Result<Vec<Vec<String>>, io::Error> {
    let file = File::open(index)?;
    let reader = BufReader::new(file);
    let mut result = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        let field = fields[column];
        let field_vec: Vec<String> = field.split('|').map(|s| s.to_owned()).collect();
        result.push(field_vec);
    }
    Ok(result)
}

/// Write a Vec<String> to a file, one string per line.
pub fn write_strings(strings: Vec<String>, output: &str) -> io::Result<()> {
    let mut file = File::create(output)?;
    for string in strings.iter() {
        writeln!(file, "{}", string)?;
    }
    Ok(())
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    #[test]
    fn test_stand_seq() {
        let test_seq = "ATGWATGAaAGCCcCCNC";
        let result = stand_seq(test_seq);
        let expected = "ATGNATGAAAGCCCCCNC";
        assert_eq!(result, expected);
    }

    #[test]
    fn test_reverse_complement() {
        let test_seq = "ATGCCAGTTAACA";
        let result = reverse_complement(test_seq);
        let expected = "TGTTAACTGGCAT";
        assert_eq!(result, expected);
    }

    #[test]
    fn test_canonical_kmers() {
        let test_seq = "ATGCCAGTTAACA";
        let result = get_canonical_kmers(test_seq, 10);
        let expected = vec!["ATGCCAGTTA", "TGCCAGTTAA", "GCCAGTTAAC", "CCAGTTAACA"];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_pack_kmer_round_trip() {
        for k in [1, 2, 5, 31] {
            let kmer = "ACGTACGTACGTACGTACGTACGTACGTACGT"[..k].to_string();
            let packed = pack_kmer(&kmer).unwrap();
            assert_eq!(unpack_kmer(packed, k), kmer);
        }
        assert_eq!(pack_kmer(""), None);
        assert_eq!(pack_kmer("ATGCN"), None);
        assert_eq!(pack_kmer("ACGTACGTACGTACGTACGTACGTACGTACGT"), None); // 32-mers
    }

    #[test]
    fn test_unpack_kmer() {
        // "AAAC" packs to 1 (A=0,C=1)
        assert_eq!(unpack_kmer(1, 4), "AAAC");
        assert_eq!(unpack_kmer(0, 3), "AAA");
        assert_eq!(unpack_kmer(0b111111, 3), "TTT");
    }

    #[test]
    fn test_packed_canonical_matches_string_canonical() {
        // The packed rolling canonicalizer must agree with the string one,
        // including reverse-complement selection and N-skipping.
        let seq = stand_seq("ATGCCAGTTAACAacgtNNNNtgtcaggctat");
        for k in [1, 4, 7, 31] {
            let string_kmers = get_canonical_kmers(&seq, k);
            let packed_kmers = get_canonical_kmers_packed(&seq, k);
            let packed_as_strings: Vec<String> =
                packed_kmers.iter().map(|p| unpack_kmer(*p, k)).collect();
            assert_eq!(string_kmers, packed_as_strings, "k = {}", k);
        }
    }

    #[test]
    fn test_packed_canonical_skips_k_at_non_atgc_run() {
        // k-mer overlapping an N must be skipped entirely.
        let seq = stand_seq("ACGTACGTNACGTACGT");
        let packed = get_canonical_kmers_packed(&seq, 5);
        let expected: Vec<String> = get_canonical_kmers(&seq, 5);
        assert_eq!(packed.len(), expected.len());
        for (p, e) in packed.iter().zip(expected.iter()) {
            assert_eq!(&unpack_kmer(*p, 5), e);
        }
    }

    #[test]
    fn test_kmer_key_round_trip() {
        assert_eq!(KmerKey::from_kmer("AAAC"), KmerKey::Packed(1));
        assert_eq!(KmerKey::from_kmer("").to_kmer(0), "");
        assert_eq!(KmerKey::from_kmer("ATGCN").to_kmer(0), "ATGCN");
        assert_eq!(
            KmerKey::from_kmer("ACGTACGTACGTACGTACGTACGTACGTACGT").to_kmer(0),
            "ACGTACGTACGTACGTACGTACGTACGTACGT"
        );
        assert_eq!(KmerKey::from_kmer("AAAC").to_kmer(4), "AAAC");
    }

    #[test]
    fn test_k_from_index_skips_empty_kmers_in_first_line() {
        let dir = tempfile::tempdir().unwrap();
        let index_path = dir.path().join("index.txt");
        // First line: empty-allele pseudo-entry with no k-mers (field 8 = "|").
        // Second line: one allele with three 5-mers.
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,AAAAA,REF|ALT,AAAAA|ACAAA,0|1,|\n",
                "1,1,200,CCCCC,REF|ALT,CCCCC|CCGGG,1|1,ACACA;CCCCT;TTTAA|GGGGG\n",
            ),
        )
        .unwrap();
        let k = k_from_index(index_path.to_str().unwrap()).unwrap();
        assert_eq!(k, 5);
    }

    #[test]
    fn test_k_from_index_uses_second_allele_when_first_is_empty() {
        let dir = tempfile::tempdir().unwrap();
        let index_path = dir.path().join("index.txt");
        // Allele 0 has no k-mers (pseudo-entry ""), allele 1 has 5-mers.
        std::fs::write(
            &index_path,
            "0,1,100,AAAAA,REF|ALT,AAAAA|ACAAA,0|1,|GGGGG\n",
        )
        .unwrap();
        let k = k_from_index(index_path.to_str().unwrap()).unwrap();
        assert_eq!(k, 5);
    }

    #[test]
    fn test_k_from_index_skips_short_lines() {
        let dir = tempfile::tempdir().unwrap();
        let index_path = dir.path().join("index.txt");
        std::fs::write(
            &index_path,
            concat!("not an index line\n", "0,1,100,CCCC,REF,CCCC,1,ACGTA\n"),
        )
        .unwrap();
        let k = k_from_index(index_path.to_str().unwrap()).unwrap();
        assert_eq!(k, 5);
    }

    #[test]
    fn test_k_from_index_errors_when_no_kmer_found() {
        let dir = tempfile::tempdir().unwrap();
        let index_path = dir.path().join("index.txt");
        std::fs::write(
            &index_path,
            concat!("0,1,100,AAAA,REF|ALT,AAAA|ACAA,0|1,|\n"),
        )
        .unwrap();
        let result = k_from_index(index_path.to_str().unwrap());
        assert!(result.is_err());
    }
}
