use std::fs::File;
use std::io::{prelude::*, BufReader};

/// Filter an index file, streaming one pass, writing rows that pass to the
/// output. Row-level filtering only (the whole variant is dropped); the index
/// format is unchanged, so the output can be used anywhere an index is.
///
/// Criteria (all optional, combined with AND):
/// - `min_alleles`: keep rows where at least this many alleles have at least
///   one allele-specific k-mer (0 keeps all; same semantics as the dedup
///   `-m` flag).
/// - `min_kmers_per_allele`: keep rows where every allele has at least this
///   many allele-specific k-mers. An allele with no k-mers (the "" pseudo-
///   entry) counts as 0, so any positive threshold drops the row. 0 keeps all.
///
/// Rows with fewer than 8 comma-separated fields are malformed and fail the
/// run, since they would corrupt downstream parsing.
pub fn filter_index(
    index: &str,
    output: &str,
    min_alleles: usize,
    min_kmers_per_allele: usize,
) -> Result<(), std::io::Error> {
    log::info!(
        "Filtering index: INDEX: {} OUTPUT: {} MIN_ALLELES: {} MIN_KMERS_PER_ALLELE: {}",
        index,
        output,
        min_alleles,
        min_kmers_per_allele
    );
    // Write to a temp sibling and rename on success, so a malformed index
    // leaves no partial output file behind.
    let tmp_output = format!("{}.tmp", output);
    let buffered_file = File::create(&tmp_output)?;
    let mut buffered_writer = std::io::BufWriter::new(buffered_file);
    let file = File::open(index)?;
    let reader = BufReader::new(file);
    let mut kept = 0;
    let mut dropped = 0;
    for (line_number, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 8 {
            let _ = std::fs::remove_file(&tmp_output);
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Malformed index line {} (expected 8 comma-separated fields): {}",
                    line_number + 1,
                    line
                ),
            ));
        }

        // Field 6 (0-based) is the number of k-mers per allele, e.g. "3|0|2".
        // The pseudo-entry "0" for an empty allele is handled like any other
        // number.
        let mut kmers_per_allele: Vec<usize> = Vec::new();
        let mut parse_error: Option<std::io::Error> = None;
        for s in fields[6].split('|') {
            match s.trim().parse::<usize>() {
                Ok(c) => kmers_per_allele.push(c),
                Err(_) => {
                    parse_error = Some(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "Invalid k-mer count '{}' on index line {} (expected a non-negative integer)",
                            s,
                            line_number + 1
                        ),
                    ));
                    break;
                }
            }
        }
        if let Some(e) = parse_error {
            let _ = std::fs::remove_file(&tmp_output);
            return Err(e);
        }

        // Both filters on the parsed counts.
        let alleles_with_kmers = kmers_per_allele.iter().filter(|&&c| c > 0).count();
        let min_alleles_ok = min_alleles == 0 || alleles_with_kmers >= min_alleles;
        let min_kmers_ok = min_kmers_per_allele == 0
            || kmers_per_allele.iter().all(|&c| c >= min_kmers_per_allele);

        if min_alleles_ok && min_kmers_ok {
            writeln!(buffered_writer, "{}", line)?;
            kept += 1;
        } else {
            dropped += 1;
            log::debug!(
                "Dropped index line {} (alleles with k-mers: {}, k-mers per allele: {:?})",
                line_number + 1,
                alleles_with_kmers,
                kmers_per_allele
            );
        }
    }
    log::info!("Kept {} index rows, dropped {} rows", kept, dropped);
    // Flush before the rename so the destination is complete.
    buffered_writer.flush()?;
    std::fs::rename(&tmp_output, output)?;
    Ok(())
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    fn write_temp(dir: &std::path::Path, name: &str, content: &str) -> String {
        let path = dir.join(name);
        std::fs::write(&path, content).unwrap();
        path.to_str().unwrap().to_string()
    }

    /// 3 lines: line 1 has alleles 2|3 (both non-zero), line 2 has 1|0 (one
    /// empty allele), line 3 has 0|0 (no k-mers anywhere).
    const INDEX: &str = concat!(
        "0,1,100,AAAAA,REF|ALT,AAAAA|ACAAA,2|3,AAAAA;CAAAC|GGGGT;TTTTA;CCCAG\n",
        "1,1,200,CCCCC,REF|ALT,CCCCC|CCGGG,1|0,CCCCC|\n",
        "2,1,300,GGGGG,REF|ALT,GGGGG|GGGGT,0|0,|\n",
    );

    fn run_filter(dir: &std::path::Path, min_alleles: usize, min_kmers: usize) -> Vec<String> {
        let index = write_temp(dir, "index.txt", INDEX);
        let out = dir.join("filtered.txt");
        filter_index(&index, out.to_str().unwrap(), min_alleles, min_kmers).unwrap();
        std::fs::read_to_string(&out)
            .unwrap()
            .lines()
            .map(|l| l.to_string())
            .collect()
    }

    #[test]
    fn test_no_flags_keeps_all() {
        let dir = tempfile::tempdir().unwrap();
        let lines = run_filter(dir.path(), 0, 0);
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0], INDEX.lines().next().unwrap());
    }

    #[test]
    fn test_min_alleles_drops_rows_with_fewer() {
        let dir = tempfile::tempdir().unwrap();
        // Line 1: 2 alleles with k-mers. Line 2: 1. Line 3: 0.
        let lines = run_filter(dir.path(), 2, 0);
        assert_eq!(lines.len(), 1);
        assert!(lines[0].starts_with("0,"));
    }

    #[test]
    fn test_min_kmers_per_allele_drops_rows_with_weak_alleles() {
        let dir = tempfile::tempdir().unwrap();
        // Line 1: 2|3 (all >= 2). Line 2: 1|0 (fails). Line 3: 0|0 (fails).
        let lines = run_filter(dir.path(), 0, 2);
        assert_eq!(lines.len(), 1);
        assert!(lines[0].starts_with("0,"));
    }

    #[test]
    fn test_min_kmers_per_allele_one_drops_empty_alleles() {
        let dir = tempfile::tempdir().unwrap();
        // Only line 1 has every allele >= 1.
        let lines = run_filter(dir.path(), 0, 1);
        assert_eq!(lines.len(), 1);
        assert!(lines[0].starts_with("0,"));
    }

    #[test]
    fn test_both_filters_combined() {
        let dir = tempfile::tempdir().unwrap();
        let lines = run_filter(dir.path(), 2, 1);
        // Line 1 passes both; lines 2 and 3 fail min_alleles.
        assert_eq!(lines.len(), 1);
        assert!(lines[0].starts_with("0,"));
    }

    #[test]
    fn test_malformed_line_is_err() {
        let dir = tempfile::tempdir().unwrap();
        let index = write_temp(dir.path(), "index.txt", "short line\n");
        let out = dir.path().join("filtered.txt");
        let result = filter_index(&index, out.to_str().unwrap(), 0, 0);
        assert!(result.is_err());
        assert!(!out.exists());
    }

    #[test]
    fn test_non_numeric_kmer_count_is_err() {
        let dir = tempfile::tempdir().unwrap();
        let index = write_temp(
            dir.path(),
            "index.txt",
            "0,1,100,AAAAA,REF|ALT,AAAAA|ACAAA,abc|3,AAAAA|GGGGG\n",
        );
        let out = dir.path().join("filtered.txt");
        let result = filter_index(&index, out.to_str().unwrap(), 0, 0);
        assert!(result.is_err());
    }

    #[test]
    fn test_empty_index_yields_empty_output() {
        let dir = tempfile::tempdir().unwrap();
        let index = write_temp(dir.path(), "index.txt", "");
        let out = dir.path().join("filtered.txt");
        filter_index(&index, out.to_str().unwrap(), 0, 0).unwrap();
        let content = std::fs::read_to_string(&out).unwrap();
        assert_eq!(content, "");
    }
}
