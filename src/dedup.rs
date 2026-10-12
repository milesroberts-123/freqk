use std::collections::HashMap;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use std::io::{BufWriter, Write};

use crate::common;

/// Pack a canonical ATGC k-mer into an integer (2 bits per base) for cheap
/// hashing, and key for counting k-mers in a hash map: packed `u64` for
/// k <= 31 (the normal case); heap `String` fallback for longer k-mers so
/// behavior is unchanged for any k. Both live in `common`, shared with `count`.
#[cfg(test)]
use crate::common::pack_kmer;
use crate::common::KmerKey;

/// Parse the k-mer field (column 7) of one index line into a flat list of
/// k-mers, one entry per k-mer occurrence (including the empty-allele
/// pseudo-entry "", matching the original counting behavior).
fn kmers_from_line(line: &str) -> Vec<&str> {
    let fields: Vec<&str> = line.split(',').collect();
    fields[7].split('|').flat_map(|s| s.split(';')).collect()
}

/// Count occurrences of every canonical reference k-mer, over all
/// chromosomes of the reference genome.
///
/// One streaming pass: each chromosome is fetched in full and its canonical
/// k-mers of length `k` are counted (k-mers overlapping non-ATGC bases are
/// skipped, as everywhere else). The VCF plays no part: ref-dedup decides
/// from reference occurrence counts alone, using the per-allele thresholds
/// documented on [`ref_dedup`].
pub fn ref_counts(fasta_path: &str, k: i64) -> Result<HashMap<String, u32>, io::Error> {
    log::info!("Counting canonical reference k-mers...");
    let faidx = rust_htslib::faidx::Reader::from_path(fasta_path)
        .map_err(|e| io::Error::other(format!("Error opening FASTA index: {}", e)))?;
    let chrom_lengths = common::read_fai(fasta_path)
        .map_err(|e| io::Error::other(format!("Error reading chromosome lengths: {}", e)))?;
    let mut counts: HashMap<String, u32> = HashMap::new();
    for (chrom, length) in chrom_lengths.iter() {
        log::debug!("Counting k-mers of CHROM: {} (length {})", chrom, length);
        let seq = common::fetch_fasta(&faidx, chrom, 0, *length as u64)
            .map_err(|e| io::Error::other(format!("Fetching {} failed: {}", chrom, e)))?;
        for ref_kmer in common::get_canonical_kmers(&seq, k.max(0) as usize) {
            *counts.entry(ref_kmer).or_insert(0) += 1;
        }
    }
    log::info!(
        "{} unique canonical reference k-mers, of which {} occur more than once",
        counts.len(),
        counts.values().filter(|&&c| c > 1).count()
    );
    Ok(counts)
}

/// Run ref-dedup end to end: count reference k-mer occurrences, then rewrite
/// the index dropping k-mers that occur outside their own variant's locus
/// (see [`remove_ref_kmers`]).
pub fn ref_dedup(
    index: &str,
    fasta_path: &str,
    output: &str,
    min_alleles: usize,
) -> Result<(), io::Error> {
    let k = common::k_from_index(index)?;
    log::info!("k is: {:?}", k);
    let counts = ref_counts(fasta_path, k)?;
    remove_ref_kmers(index, output, &counts, min_alleles)
}

/// Maximum reference occurrences a kept k-mer may have, given its allele
/// position within its index row: the REF allele (group 0) is allowed one
/// occurrence, its own variant's locus; an ALT-allele k-mer can never match
/// at its own variant (its allele base differs from the reference base
/// there), so every reference occurrence is a stray and must be dropped.
fn max_kept_occurrences(group_index: usize) -> u32 {
    if group_index == 0 {
        1
    } else {
        0
    }
}

/// Count the alleles of a deduplicated variant that have at least one
/// non-empty k-mer.
fn count_alleles_with_kmers<S: AsRef<str>>(dedup_kmers: &[Vec<S>]) -> usize {
    dedup_kmers
        .iter()
        .filter(|inner_vec| inner_vec.iter().any(|s| !s.as_ref().is_empty()))
        .count()
}

/// Remove k-mers found in the reference from an index and write a new index.
///
/// DROP RULE: every index k-mer spans its variant, so a REF-allele k-mer
/// always occurs exactly once in the reference — at its own variant — while
/// an ALT-allele k-mer never does. A k-mer therefore has an occurrence
/// outside its own variant's locus iff its reference count exceeds its
/// allele's threshold (REF: 1, ALT: 0); such k-mers are dropped.
///
/// Streams the index in a single pass instead of loading it into memory: each
/// line is parsed, its k-mers filtered against the reference counts, and the
/// line rewritten immediately. Memory is proportional to one index line, not
/// index size. Rows are dropped entirely when fewer than `min_alleles` alleles
/// retain at least one k-mer (0 keeps every row).
pub fn remove_ref_kmers(
    index: &str,
    output: &str,
    ref_counts: &HashMap<String, u32>,
    min_alleles: usize,
) -> Result<(), io::Error> {
    log::info!("Writing new index, removing k-mers with excess reference occurrences...");
    let mut buffered_file = BufWriter::new(File::create(output)?);
    let file = File::open(index)?;
    let reader = BufReader::new(file);
    let mut dropped = 0;
    let mut written = 0;
    for (line_number, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 8 {
            log::warn!(
                "Skipping malformed index line {} (expected 8 comma-separated fields): {}",
                line_number + 1,
                line
            );
            continue;
        }

        let kmers_by_allele: Vec<Vec<&str>> = fields[7]
            .split('|')
            .map(|s| s.split(';').collect())
            .collect();
        let dedup_kmers: Vec<Vec<&str>> = kmers_by_allele
            .iter()
            .enumerate()
            .map(|(group_index, inner_vec)| {
                let max_kept = max_kept_occurrences(group_index);
                inner_vec
                    .iter()
                    .filter(|s| match ref_counts.get(**s) {
                        None => true,
                        Some(&count) => count <= max_kept,
                    })
                    .cloned()
                    .collect()
            })
            .collect();

        if min_alleles > 0 && count_alleles_with_kmers(&dedup_kmers) < min_alleles {
            dropped += 1;
            continue;
        }
        written += 1;

        let num_kmers_per_allele = dedup_kmers
            .iter()
            .map(|inner_vec| {
                if *inner_vec == vec![""] {
                    "0".to_string()
                } else {
                    inner_vec.len().to_string()
                }
            })
            .collect::<Vec<String>>()
            .join("|");
        let dedup_string = dedup_kmers
            .iter()
            .map(|inner_vec| inner_vec.join(";"))
            .collect::<Vec<String>>()
            .join("|");
        let parts = [
            fields[0],
            fields[1],
            fields[2],
            fields[3],
            fields[4],
            fields[5],
            &num_kmers_per_allele,
            &dedup_string,
        ];
        writeln!(buffered_file, "{}", parts.join(","))?;
    }
    if min_alleles > 0 && dropped > 0 {
        log::info!(
            "Dropped {} index rows with fewer than {} alleles retaining k-mers",
            dropped,
            min_alleles
        );
    }
    log::info!("Wrote {} index rows", written);
    Ok(())
}

/// Remove k-mers shared across variants from an index and write a new index.
///
/// Streams the index file twice instead of loading it into memory: the first
/// pass counts k-mer occurrences in a `HashMap<KmerKey, u32>`, the second pass
/// re-reads each line, drops duplicated k-mers, and rewrites the line.
/// Memory is proportional to the number of unique k-mers, not index size.
/// Rows are dropped entirely when fewer than `min_alleles` alleles retain at
/// least one k-mer (0 keeps every row).
pub fn find_dup_kmers_across_var(
    index: &str,
    output: &str,
    min_alleles: usize,
) -> Result<(), io::Error> {
    log::info!("First pass: counting allele-specific k-mers...");
    let mut counts: HashMap<KmerKey, u32> = HashMap::new();
    {
        let file = File::open(index)?;
        let reader = BufReader::new(file);
        for line_result in reader.lines() {
            let line = line_result?;
            for kmer in kmers_from_line(&line) {
                *counts.entry(KmerKey::from_kmer(kmer)).or_insert(0) += 1;
            }
        }
    }

    log::info!("Removing k-mers found only once from the count table...");
    counts.retain(|_, count| *count > 1);

    log::info!("Second pass: removing k-mers found more than once and writing new index...");
    let mut buffered_file = BufWriter::new(File::create(output)?);
    let file = File::open(index)?;
    let reader = BufReader::new(file);
    let mut dropped = 0;

    for (line_number, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 8 {
            log::warn!(
                "Skipping malformed index line {} (expected 8 comma-separated fields): {}",
                line_number + 1,
                line
            );
            continue;
        }

        let kmers_by_allele: Vec<Vec<&str>> = fields[7]
            .split('|')
            .map(|s| s.split(';').collect())
            .collect();
        let dedup_kmers: Vec<Vec<&str>> = kmers_by_allele
            .iter()
            .map(|inner_vec| {
                inner_vec
                    .iter()
                    .filter(|s| !counts.contains_key(&KmerKey::from_kmer(s)))
                    .cloned()
                    .collect()
            })
            .collect();

        if min_alleles > 0 && count_alleles_with_kmers(&dedup_kmers) < min_alleles {
            dropped += 1;
            continue;
        }

        let num_kmers_per_allele = dedup_kmers
            .iter()
            .map(|inner_vec| {
                if *inner_vec == vec![""] {
                    "0".to_string()
                } else {
                    inner_vec.len().to_string()
                }
            })
            .collect::<Vec<String>>()
            .join("|");

        let dedup_string = dedup_kmers
            .iter()
            .map(|inner_vec| inner_vec.join(";"))
            .collect::<Vec<String>>()
            .join("|");

        let parts = [
            fields[0],
            fields[1],
            fields[2],
            fields[3],
            fields[4],
            fields[5],
            &num_kmers_per_allele,
            &dedup_string,
        ];
        writeln!(buffered_file, "{}", parts.join(","))?;
    }

    if min_alleles > 0 && dropped > 0 {
        log::info!(
            "Dropped {} index rows with fewer than {} alleles retaining k-mers",
            dropped,
            min_alleles
        );
    }
    log::info!("Writing successful! :D");
    Ok(())
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    #[test]
    fn test_pack_kmer() {
        // "AAAC" packs to 0b000001 = 1
        assert_eq!(pack_kmer("AAAC"), Some(1));
        // Full 31-mer is the packing limit: 0b0011011010... (2 bits per base)
        let packed_31 = "ACGTACGTACGTACGTACGTACGTACGTACC"
            .bytes()
            .fold(0u64, |acc, c| {
                (acc << 2)
                    | match c {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        _ => 3,
                    }
            });
        assert_eq!(
            pack_kmer("ACGTACGTACGTACGTACGTACGTACGTACC"),
            Some(packed_31)
        );
        // Empty allele pseudo-entry and non-ATGC fall back to None
        assert_eq!(pack_kmer(""), None);
        assert_eq!(pack_kmer("ATGCN"), None);
        // Longer than 31 falls back to None
        assert_eq!(pack_kmer("ACGTACGTACGTACGTACGTACGTACGTACGT"), None);
    }

    #[test]
    fn test_kmer_key_roundtrip() {
        assert_eq!(KmerKey::from_kmer("AAAC"), KmerKey::Packed(1));
        assert_eq!(KmerKey::from_kmer(""), KmerKey::Str(String::new()));
        assert_eq!(
            KmerKey::from_kmer("ACGTACGTACGTACGTACGTACGTACGTACGT"),
            KmerKey::Str("ACGTACGTACGTACGTACGTACGTACGTACGT".to_string())
        );
    }

    #[test]
    fn test_kmers_from_line() {
        let line = "0,1,100,SEQ,REF|ALT,SEQ|ALT,2|1,AAA;CCC|GGG";
        assert_eq!(kmers_from_line(line), vec!["AAA", "CCC", "GGG"]);
        // Empty allele pseudo-entry is preserved for counting
        let line_empty = "0,1,100,SEQ,REF|ALT,SEQ|ALT,0|0,|";
        assert_eq!(kmers_from_line(line_empty), vec!["", ""]);
    }

    #[test]
    fn test_count_alleles_with_kmers() {
        assert_eq!(count_alleles_with_kmers(&[vec!["AAA"], vec!["GGG"]]), 2);
        assert_eq!(count_alleles_with_kmers(&[vec![], vec!["GGG"]]), 1);
        assert_eq!(count_alleles_with_kmers(&[vec![""], vec![]]), 0);
    }

    #[test]
    fn test_find_dup_kmers_across_var_end_to_end() {
        let dir = tempfile::tempdir().expect("tempdir");
        let index_path = dir.path().join("index.txt");
        let output_path = dir.path().join("dedup.txt");
        // Variant 1: AAA appears twice within the line (duplicate) -> removed
        // everywhere; CCC and TTT appear once and stay.
        // Variant 2 has an empty allele (pseudo-entry "" in both fields).
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,AAAA,REF|ALT,AAAA|ACAA,2|2,AAA;CCC|AAA;TTT\n",
                "1,1,200,CCCC,REF|ALT,CCCC|CCGG,0|1,|GGG\n",
            ),
        )
        .expect("write index");
        find_dup_kmers_across_var(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            0,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        let expected = concat!(
            "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,CCC|TTT\n",
            // GGG appears once -> stays; empty allele still yields 0.
            "1,1,200,CCCC,REF|ALT,CCCC|CCGG,0|1,|GGG\n",
        );
        assert_eq!(result, expected);
    }

    #[test]
    fn test_find_dup_kmers_across_var_min_alleles_drops_rows() {
        let dir = tempfile::tempdir().expect("tempdir");
        let index_path = dir.path().join("index.txt");
        let output_path = dir.path().join("dedup.txt");
        // Variant 1 keeps 1 k-mer in each allele after dedup (CCC, TTT unique).
        // Variant 2 loses its only k-mer (AAA appears three times -> removed),
        // leaving both alleles empty.
        // Variant 3 has a pre-existing empty allele and GGG (unique) in allele 1.
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,CCC|TTT\n",
                "1,1,200,CCCC,REF|ALT,CCCC|CCGG,1|1,AAA|AAA\n",
                "2,1,300,TTTT,REF|ALT,TTTT|TTGG,1|0,GGG|\n",
            ),
        )
        .expect("write index");
        // min_alleles = 1: variant 2 is dropped (0 alleles with k-mers),
        // variants 1 and 3 stay.
        find_dup_kmers_across_var(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            1,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        let expected = concat!(
            "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,CCC|TTT\n",
            // GGG is unique -> retained; the empty allele yields 0.
            "2,1,300,TTTT,REF|ALT,TTTT|TTGG,1|0,GGG|\n",
        );
        assert_eq!(result, expected);
    }

    #[test]
    fn test_find_dup_kmers_across_var_min_alleles_two() {
        let dir = tempfile::tempdir().expect("tempdir");
        let index_path = dir.path().join("index.txt");
        let output_path = dir.path().join("dedup.txt");
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,CCC|TTT\n",
                "1,1,200,CCCC,REF|ALT,CCCC|CCGG,0|1,|GGG\n",
            ),
        )
        .expect("write index");
        // min_alleles = 2: variant 2 has only one allele with a k-mer -> dropped.
        find_dup_kmers_across_var(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            2,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        let expected = "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,CCC|TTT\n";
        assert_eq!(result, expected);
    }

    #[test]
    fn test_remove_ref_kmers_min_alleles() {
        let dir = tempfile::tempdir().expect("tempdir");
        let index_path = dir.path().join("index.txt");
        let output_path = dir.path().join("dedup.txt");
        // Variant row 1: REF allele k-mer AAA occurs twice in the reference —
        // one occurrence is its own locus, so a second one is a stray ->
        // dropped; TTT never occurs -> kept. Row 2: GGG occurs once (its own
        // locus) -> kept for the REF allele; CCC occurs twice -> dropped for
        // the ALT allele.
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,AAA|TTT\n",
                "1,1,200,CCCC,REF|ALT,CCCC|CCGG,1|1,GGG|CCC\n",
            ),
        )
        .expect("write index");
        let mut ref_counts = HashMap::new();
        ref_counts.insert("AAA".to_string(), 2);
        ref_counts.insert("GGG".to_string(), 1);
        ref_counts.insert("CCC".to_string(), 2);

        // min_alleles = 0 keeps all rows.
        remove_ref_kmers(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            &ref_counts,
            0,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        let expected = concat!(
            "0,1,100,AAAA,REF|ALT,AAAA|ACAA,0|1,|TTT\n",
            "1,1,200,CCCC,REF|ALT,CCCC|CCGG,1|0,GGG|\n",
        );
        assert_eq!(result, expected);

        // min_alleles = 1 keeps both rows (each retains one allele).
        remove_ref_kmers(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            &ref_counts,
            1,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        assert_eq!(result, expected);

        // min_alleles = 2 drops both rows (one allele each).
        remove_ref_kmers(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            &ref_counts,
            2,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        assert_eq!(result, "");
    }

    #[test]
    fn test_remove_ref_kmers_two_tier_thresholds() {
        let dir = tempfile::tempdir().expect("tempdir");
        let index_path = dir.path().join("index.txt");
        let output_path = dir.path().join("dedup.txt");
        // The same single-occurrence k-mer is kept in a REF allele but dropped
        // in an ALT allele; any counted k-mer is dropped from both.
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,AAAA,REF|ALT2,AAAA|ACAA|AGAA,1|1|1,AAA|TTT|AAA\n",
                "1,1,200,CCCC,REF|ALT,CCCC|CCGG,1|1,AAA|AAA\n",
            ),
        )
        .expect("write index");
        let mut ref_counts = HashMap::new();
        ref_counts.insert("AAA".to_string(), 1);
        remove_ref_kmers(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            &ref_counts,
            0,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        // Row 0: REF group AAA (count 1, own locus) and ALT group TTT
        // (count 0, no occurrence) are both kept; ALT group AAA (count 1,
        // a stray) is dropped.
        let expected = concat!(
            "0,1,100,AAAA,REF|ALT2,AAAA|ACAA|AGAA,1|1|0,AAA|TTT|\n",
            "1,1,200,CCCC,REF|ALT,CCCC|CCGG,1|0,AAA|\n",
        );
        assert_eq!(result, expected);
    }

    #[test]
    fn test_ref_counts() {
        let dir = tempfile::tempdir().expect("tempdir");
        let fasta_path = dir.path().join("ref.fasta");
        // k = 5: ACGTA occurs on chromosome a twice (own canonical form);
        // GGAAA occurs once; chromosome c is shorter than k and contributes
        // nothing.
        std::fs::write(
            &fasta_path,
            concat!(
                ">a\n",
                "ACGTAACGTAGGGG\n",
                ">b\n",
                "GGAAA\n",
                ">c\n",
                "AC\n",
            ),
        )
        .expect("write fasta");
        // .fai: name, length, offset, linebases, linewidth. Offsets account
        // for each ">name\n" header line (3 bytes) before its sequence.
        std::fs::write(
            dir.path().join("ref.fasta.fai"),
            "a\t14\t3\t14\t15\nb\t5\t21\t5\t6\nc\t2\t30\t2\t3\n",
        )
        .expect("write fai");
        let counts = ref_counts(fasta_path.to_str().expect("path"), 5).expect("counts");
        assert_eq!(counts.get("ACGTA"), Some(&2));
        assert_eq!(counts.get("GGAAA"), Some(&1));
        // "GGGG" plus shorter k-mers are not counted; only length-5 k-mers,
        // of which ACGTA (at offsets 0 and 5) is the only duplicate.
        assert_eq!(counts.get("AAAAA"), None);
        let duplicated = counts.values().filter(|&&c| c > 1).count();
        assert_eq!(duplicated, 1);
    }
}
