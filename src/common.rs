use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;

/// Get k-mer length from the first line of an index file.
pub fn k_from_index(index: &str) -> Result<i64, io::Error> {
    let file = File::open(index)?;
    let mut reader = BufReader::new(file);
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;
    let split_line: Vec<&str> = first_line.split(',').collect();
    let kmers = split_line[7];
    let kmers_list: Vec<&str> = kmers.split('|').collect();
    let kmers_by_allele: Vec<Vec<&str>> =
        kmers_list.iter().map(|s| s.split(';').collect()).collect();
    let first_inner_vec = &kmers_by_allele[0];
    let first_element = first_inner_vec[0];
    Ok(first_element.len() as i64)
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
}
