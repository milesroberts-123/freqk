use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use crate::common;

mod freq_from_hetmers;

/// Read a k-mer count table (two tab-separated columns) into (sequence, count)
/// pairs. The whole table is parsed before any format check runs, so format
/// errors are reported against the full table. Unreadable files, unparseable
/// counts, and empty k-mer strings are hard errors; lines with the wrong
/// column count are skipped with a warning.
fn load_kmers(input: &str, minimum: usize) -> Result<Vec<(String, usize)>, std::io::Error> {
    log::info!("Loading k-mer count file {}...", input);
    let file = File::open(input)?;
    let reader = BufReader::new(file);
    let mut kmers = Vec::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() != 2 {
            log::warn!(
                "Skipping line {} that does not have two tab-separated columns:",
                line_number + 1
            );
            log::warn!("{}", line);
            continue;
        }

        let seq = parts[0];
        if seq.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Empty k-mer on line {} of {} (hetmers cannot process empty k-mers)",
                    line_number + 1,
                    input
                ),
            ));
        }
        let count: usize = parts[1].parse().map_err(|e| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Invalid count value on line {} of {}: {}",
                    line_number + 1,
                    input,
                    e
                ),
            )
        })?;

        if count >= minimum {
            kmers.push((seq.to_string(), count));
        }
    }

    Ok(kmers)
}

/// Check that k-mers are lexicographically sorted (first 1000 elements).
fn check_sort(seqs: &[String]) -> bool {
    let limit = seqs.len().min(1000);
    let seqs_sub = &seqs[..limit];

    let mut sorted_seqs = seqs_sub.to_vec();
    sorted_seqs.sort();

    let result = seqs_sub == sorted_seqs;

    log::info!("Input sorted: {}", result);
    result
}

/// Check that only ATGC are in the alphabet (first 1000 elements).
fn check_letters(seqs: &[String]) -> bool {
    let limit = seqs.len().min(1000);
    let seqs_sub = &seqs[..limit];

    let result = seqs_sub
        .iter()
        .all(|seq| seq.chars().all(|c| matches!(c, 'A' | 'T' | 'G' | 'C')));
    log::info!("Only ATGC: {}", result);
    result
}

/// Run all input checks, returning an error naming the first offending k-mer
/// on failure.
fn all_checks(seqs: &[String]) -> Result<(), std::io::Error> {
    log::info!("Checking input format...");
    if !check_sort(seqs) {
        // Find the first out-of-order pair for a helpful message.
        let offender = seqs
            .windows(2)
            .find(|w| w[0] > w[1])
            .map(|w| format!("'{}' before '{}'", w[0], w[1]))
            .unwrap_or_else(|| "unknown position".to_string());
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "Input k-mers are not lexicographically sorted ({})",
                offender
            ),
        ));
    }

    if !check_letters(seqs) {
        let offender = seqs
            .iter()
            .find(|seq| !seq.chars().all(|c| matches!(c, 'A' | 'T' | 'G' | 'C')))
            .cloned()
            .unwrap_or_default();
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "Input k-mers contain characters other than ATGC ('{}')",
                offender
            ),
        ));
    }

    Ok(())
}

/// Remove the central base from each k-mer.
fn extract_border(seqs: &[String]) -> Result<Vec<String>, std::io::Error> {
    if seqs.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "No k-mers retained after applying the minimum count filter",
        ));
    }
    let k = seqs[0].len();
    log::info!("k is {}", k);
    let k_half = k / 2;

    Ok(seqs
        .iter()
        .map(|s| format!("{}{}", &s[..k_half], &s[k_half + 1..]))
        .collect())
}

/// Reverse complement each sequence.
fn rev_comp(seqs: &[String]) -> Vec<String> {
    log::debug!("Reverse complementing...");
    seqs.iter().map(|s| common::reverse_complement(s)).collect()
}

/// Hash each sequence with SHA-256, keeping the first 8 bytes.
fn hash_seqs(seqs: &[String]) -> Vec<u64> {
    log::debug!("Hashing...");
    let hash_fn = |s: &String| {
        let mut hasher = Sha256::new();
        hasher.update(s.as_bytes());
        u64::from_be_bytes(hasher.finalize()[..8].try_into().unwrap())
    };

    seqs.iter().map(hash_fn).collect()
}

/// Take the elementwise minimum of two hash vectors.
fn min_hash(hash1: Vec<u64>, hash2: Vec<u64>) -> Vec<u64> {
    log::debug!("Getting the minimum hash...");
    hash1
        .iter()
        .zip(hash2.iter())
        .map(|(x, y)| *x.min(y))
        .collect()
}

/// Group identical hashes, mapping each hash to the indices where it occurs.
/// A BTreeMap keeps the output files in ascending-hash order, deterministic
/// across runs.
fn group_hashes(hashes: Vec<u64>) -> BTreeMap<u64, Vec<usize>> {
    log::debug!("Grouping unique hashes into a dictionary...");
    let mut d: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    for (i, num) in hashes.iter().enumerate() {
        d.entry(*num).or_default().push(i);
    }

    d
}

/// Keep only hash groups with exactly `alleles` members.
fn filter_groups(input: BTreeMap<u64, Vec<usize>>, alleles: usize) -> BTreeMap<u64, Vec<usize>> {
    log::debug!("Filtering hash groups by number of alleles...");
    input
        .into_iter()
        .filter(|(_, v)| v.len() == alleles)
        .collect()
}

/// Extract hetmer sequences, counts, and hashes from hash groups.
fn extract_hetmers(
    hashdict: BTreeMap<u64, Vec<usize>>,
    seqs: Vec<String>,
    counts: Vec<usize>,
) -> (Vec<String>, Vec<String>, Vec<u64>) {
    log::debug!("Extracting counts and sequences...");

    let hetmer_seqs: Vec<String> = hashdict
        .values()
        .map(|indices| {
            indices
                .iter()
                .map(|&i| seqs[i].clone())
                .collect::<Vec<String>>()
                .join(",")
        })
        .collect();

    let hetmer_counts: Vec<String> = hashdict
        .values()
        .map(|indices| {
            indices
                .iter()
                .map(|&i| counts[i].to_string())
                .collect::<Vec<String>>()
                .join(",")
        })
        .collect();

    (hetmer_seqs, hetmer_counts, hashdict.into_keys().collect())
}

/// Write a vector of strings to a file, one line per element.
fn write_file(output: &[String], prefix: &str, suffix: &str) -> std::io::Result<()> {
    log::info!("Saving results to {}_{}...", prefix, suffix);
    let mut file = File::create(format!("{}_{}", prefix, suffix))?;
    if !output.is_empty() {
        writeln!(file, "{}", output.join("\n"))?;
    }
    Ok(())
}

/// Find hetmers in a k-mer count table and write results to output files.
/// Input is a table of k-mer counts from any counter (kmc, jellyfish, or
/// `freqk count -c`), with two tab-separated columns: k-mer, count.
#[allow(clippy::too_many_arguments)]
pub fn kmers_to_hetmers(
    input: &str,
    output: &str,
    minimum: usize,
    alleles: usize,
    pool: i32,
    coverage: f64,
    alpha: f64,
    beta: f64,
    sigma: f64,
) -> Result<(), std::io::Error> {
    // load k-mers (the whole table is parsed before format checks run)
    let kmers = load_kmers(input, minimum)?;
    let seqs: Vec<String> = kmers.iter().map(|(seq, _)| seq.clone()).collect();
    let counts: Vec<usize> = kmers.iter().map(|(_, count)| *count).collect();

    // input checks
    all_checks(&seqs)?;

    // remove central base from each k-mer
    let borders = extract_border(&seqs)?;

    // reverse complement borders
    let revborders = rev_comp(&borders);

    // get hash of borders
    let hashbord = hash_seqs(&borders);
    let hashrevbord = hash_seqs(&revborders);

    // compare forward and reverse hash and take the min
    let min_hashes = min_hash(hashbord, hashrevbord);

    // group the hashes into a dictionary
    let grouped_hashes = group_hashes(min_hashes);

    // remove hashes that had only one or more than 2 k-mers per group
    let filtered_groups = filter_groups(grouped_hashes, alleles);

    // extract sequences for each hash group
    let hetmers = extract_hetmers(filtered_groups, seqs, counts);

    // empirical frequencies
    let empirical_frequencies = freq_from_hetmers::counts_to_frequencies(&hetmers.1);

    // bayesian allele states
    let bayes_states =
        freq_from_hetmers::counts_to_bayes_state(&hetmers.1, pool, coverage, minimum, alpha, beta)?;

    // check for hetmers with weirdly high coverage
    let check_these_hetmers =
        freq_from_hetmers::high_cov_hetmers(&hetmers.1, sigma, pool, coverage);

    // write output files
    write_file(&hetmers.0, output, "seqs.csv")?;
    write_file(&hetmers.1, output, "counts.csv")?;
    write_file(
        &hetmers
            .2
            .iter()
            .map(|num| num.to_string())
            .collect::<Vec<String>>(),
        output,
        "hashes.csv",
    )?;
    write_file(
        &bayes_states
            .into_iter()
            .map(|s| s.to_string())
            .collect::<Vec<String>>(),
        output,
        "bayes_states.csv",
    )?;
    write_file(&empirical_frequencies, output, "empirical_freqs.csv")?;
    write_file(&check_these_hetmers, output, "bad_hetmers.csv")?;
    Ok(())
}

// test functions
#[cfg(test)]
mod unit_tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn odd_k_borders() {
        let test_vec = vec![
            "ATGCA".to_string(),
            "TTGAT".to_string(),
            "GGATA".to_string(),
        ];
        let result = extract_border(&test_vec).unwrap();
        let expected = vec!["ATCA".to_string(), "TTAT".to_string(), "GGTA".to_string()];
        assert_eq!(result, expected);
    }

    #[test]
    fn even_k_borders() {
        let test_vec = vec![
            "ATGCAT".to_string(),
            "TTGATC".to_string(),
            "GGATAA".to_string(),
        ];
        let result = extract_border(&test_vec).unwrap();
        let expected = vec![
            "ATGAT".to_string(),
            "TTGTC".to_string(),
            "GGAAA".to_string(),
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn empty_input_border_is_err() {
        let result = extract_border(&[]);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("No k-mers retained"));
    }

    #[test]
    fn short_rev_comp() {
        let test_vec = vec![
            "ATGCAT".to_string(),
            "TTGATC".to_string(),
            "GGATAA".to_string(),
        ];
        let result = rev_comp(&test_vec);
        let expected = vec![
            "ATGCAT".to_string(),
            "GATCAA".to_string(),
            "TTATCC".to_string(),
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn small_min_hash() {
        let test_vec_1 = vec![59888, 1, 100];
        let test_vec_2 = vec![59887, 5000, 101];
        let result = min_hash(test_vec_1, test_vec_2);
        let expected = vec![59887, 1, 100];
        assert_eq!(result, expected);
    }

    #[test]
    fn seqs_from_hashmap() {
        let mut input = BTreeMap::new();
        input.insert(9875, vec![0, 4]);
        input.insert(1111, vec![1, 2]);
        input.insert(2222, vec![3, 5]);

        let seqs = vec![
            "TCGTC".to_string(),
            "AATAA".to_string(),
            "AAGAA".to_string(),
            "GATGA".to_string(),
            "TCATC".to_string(),
            "GAAGA".to_string(),
        ];
        let counts = vec![1, 10, 9, 2, 6, 100];

        let result = extract_hetmers(input, seqs, counts);

        let actual: HashSet<_> = result.0.into_iter().collect();
        let expected: HashSet<_> = vec![
            "TCGTC,TCATC".to_string(),
            "GATGA,GAAGA".to_string(),
            "AATAA,AAGAA".to_string(),
        ]
        .into_iter()
        .collect();

        assert_eq!(actual, expected);
    }

    #[test]
    fn counts_from_hashmap() {
        let mut input = BTreeMap::new();
        input.insert(9875, vec![0, 4]);
        input.insert(1111, vec![1, 2]);
        input.insert(2222, vec![3, 5]);

        let seqs = vec![
            "TCGTC".to_string(),
            "AATAA".to_string(),
            "AAGAA".to_string(),
            "GATGA".to_string(),
            "TCATC".to_string(),
            "GAAGA".to_string(),
        ];
        let counts = vec![1, 10, 9, 2, 6, 100];

        let result = extract_hetmers(input, seqs, counts);

        let actual: HashSet<_> = result.1.into_iter().collect();
        let expected: HashSet<_> = vec!["1,6".to_string(), "2,100".to_string(), "10,9".to_string()]
            .into_iter()
            .collect();

        assert_eq!(actual, expected);
    }

    #[test]
    fn hashes_from_hashmap() {
        let mut input = BTreeMap::new();
        input.insert(9875, vec![0, 4]);
        input.insert(1111, vec![1, 2]);
        input.insert(2222, vec![3, 5]);

        let seqs = vec![
            "TCGTC".to_string(),
            "AATAA".to_string(),
            "AAGAA".to_string(),
            "GATGA".to_string(),
            "TCATC".to_string(),
            "GAAGA".to_string(),
        ];
        let counts = vec![1, 10, 9, 2, 6, 100];

        let result = extract_hetmers(input, seqs, counts);

        let actual: HashSet<_> = result.2.into_iter().collect();
        let expected: HashSet<_> = vec![9875, 2222, 1111].into_iter().collect();

        assert_eq!(actual, expected);
    }

    #[test]
    fn two_alleles() {
        let mut input = BTreeMap::new();
        input.insert(1, vec![0, 1]);
        input.insert(2, vec![2, 3, 4]); // should be filtered out
        input.insert(3, vec![2, 3, 4, 5]); // should be filtered out
        input.insert(4, vec![2]); // should be filtered out
        input.insert(5, vec![5, 6]);

        let alleles = 2;
        let result = filter_groups(input, alleles);

        let mut expected = BTreeMap::new();
        expected.insert(1, vec![0, 1]);
        expected.insert(5, vec![5, 6]);

        assert_eq!(result, expected);
    }

    #[test]
    fn three_alleles() {
        let mut input = BTreeMap::new();
        input.insert(1, vec![0, 1]);
        input.insert(2, vec![2, 3, 4]);
        input.insert(3, vec![2, 3, 4, 5]);
        input.insert(4, vec![2]);
        input.insert(5, vec![5, 6]);

        let alleles = 3;
        let result = filter_groups(input, alleles);

        let mut expected = BTreeMap::new();
        expected.insert(2, vec![2, 3, 4]);

        assert_eq!(result, expected);
    }

    #[test]
    fn four_alleles() {
        let mut input = BTreeMap::new();
        input.insert(1, vec![0, 1]);
        input.insert(2, vec![2, 3, 4]);
        input.insert(3, vec![2, 3, 4, 5]);
        input.insert(4, vec![2]);
        input.insert(5, vec![5, 6]);

        let alleles = 4;
        let result = filter_groups(input, alleles);

        let mut expected = BTreeMap::new();
        expected.insert(3, vec![2, 3, 4, 5]);

        assert_eq!(result, expected);
    }

    #[test]
    fn filter_groups_no_matches() {
        let mut input = BTreeMap::new();
        input.insert(1, vec![0]);
        input.insert(2, vec![1, 2, 3]);

        let result = filter_groups(input, 2);
        let expected: BTreeMap<u64, Vec<usize>> = BTreeMap::new();

        assert_eq!(result, expected);
    }

    #[test]
    fn all_checks_accepts_sorted_atgc() {
        let seqs = vec!["AAAAA".to_string(), "CCCCC".to_string()];
        assert!(all_checks(&seqs).is_ok());
    }

    #[test]
    fn all_checks_rejects_unsorted_with_offender() {
        let seqs = vec!["GGGGG".to_string(), "CCCCC".to_string()];
        let result = all_checks(&seqs);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("'GGGGG' before 'CCCCC'"));
    }

    #[test]
    fn all_checks_rejects_non_atgc_with_offender() {
        let seqs = vec!["AAAAN".to_string()];
        let result = all_checks(&seqs);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("AAAAN"));
    }

    #[test]
    fn write_file_empty_vec_creates_empty_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out");
        write_file(&[], path.to_str().unwrap(), "test.csv").unwrap();
        let content = std::fs::read_to_string(dir.path().join("out_test.csv")).unwrap();
        assert_eq!(content, "");
    }
}
