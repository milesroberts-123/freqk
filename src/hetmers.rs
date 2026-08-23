use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use crate::common;

mod freq_from_hetmers;

/// Read a k-mer count table (two tab-separated columns) into (sequence, count) pairs.
fn load_kmers(input: &str, minimum: usize) -> Vec<(String, usize)> {
    log::info!("Loading k-mer count file {}...", input);
    let file = File::open(input).expect("Unable to open file");
    let reader = BufReader::new(file);
    let mut kmers = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() != 2 {
            log::warn!("Skipping line that does not have two tab-separated columns:");
            log::warn!("{}", line);
            continue;
        }

        let seq = parts[0].to_string();
        let count: usize = parts[1].parse().expect("Invalid count value");

        if count >= minimum {
            kmers.push((seq, count));
        }
    }

    kmers
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

/// Run all input checks, panicking on failure.
fn all_checks(seqs: &[String]) {
    log::info!("Checking input format...");
    if !check_sort(seqs) {
        panic!("Input k-mers are not lexicographically sorted");
    }

    if !check_letters(seqs) {
        panic!("Input k-mers contain characters other than ATGC")
    }
}

/// Remove the central base from each k-mer.
fn extract_border(seqs: &[String]) -> Vec<String> {
    let k = seqs[0].len();
    log::info!("k is {}", k);
    let k_half = k / 2;

    seqs.iter()
        .map(|s| format!("{}{}", &s[..k_half], &s[k_half + 1..]))
        .collect()
}

/// Reverse complement each sequence.
fn rev_comp(seqs: &[String]) -> Vec<String> {
    log::debug!("Reverse complementing...");
    seqs.iter().map(|s| common::reverse_complement(s)).collect()
}

/// Hash each sequence with SHA-256, keeping the first 8 bytes.
fn hash_seqs(seqs: Vec<String>) -> Vec<u64> {
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
fn group_hashes(hashes: Vec<u64>) -> HashMap<u64, Vec<usize>> {
    log::debug!("Grouping unique hashes into a dictionary...");
    let mut d: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, num) in hashes.iter().enumerate() {
        d.entry(*num).or_default().push(i);
    }

    d
}

/// Keep only hash groups with exactly `alleles` members.
fn filter_groups(input: HashMap<u64, Vec<usize>>, alleles: usize) -> HashMap<u64, Vec<usize>> {
    log::debug!("Filtering hash groups by number of alleles...");
    input
        .into_iter()
        .filter(|(_, v)| v.len() == alleles)
        .collect()
}

/// Extract hetmer sequences, counts, and hashes from hash groups.
fn extract_hetmers(
    hashdict: HashMap<u64, Vec<usize>>,
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
fn write_file(output: &[String], prefix: &str, suffix: &str) {
    log::info!("Saving results to {}_{}...", prefix, suffix);
    let mut file = File::create(format!("{}_{}", prefix, suffix)).expect("Unable to create file");
    writeln!(file, "{}", output.join("\n")).expect("Unable to write to file");
}

/// Find hetmers in a k-mer count table and write results to output files.
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
) {
    // load k-mers
    let kmers = load_kmers(input, minimum);
    let seqs: Vec<String> = kmers.iter().map(|(seq, _)| seq.clone()).collect();
    let counts: Vec<usize> = kmers.iter().map(|(_, count)| *count).collect();

    // input checks
    all_checks(&seqs);

    // remove central base from each k-mer
    let borders = extract_border(&seqs);

    // reverse complement borders
    let revborders = rev_comp(&borders);

    // get hash of borders
    let hashbord = hash_seqs(borders);
    let hashrevbord = hash_seqs(revborders);

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
        freq_from_hetmers::counts_to_bayes_state(&hetmers.1, pool, coverage, minimum, alpha, beta);

    // check for hetmers with weirdly high coverage
    let check_these_hetmers =
        freq_from_hetmers::high_cov_hetmers(&hetmers.1, sigma, pool, coverage);

    // write output files
    write_file(&hetmers.0, output, "seqs.csv");
    write_file(&hetmers.1, output, "counts.csv");
    write_file(
        &hetmers
            .2
            .iter()
            .map(|num| num.to_string())
            .collect::<Vec<String>>(),
        output,
        "hashes.csv",
    );
    write_file(
        &bayes_states
            .into_iter()
            .map(|s| s.to_string())
            .collect::<Vec<String>>(),
        output,
        "bayes_states.csv",
    );
    write_file(&empirical_frequencies, output, "empirical_freqs.csv");
    write_file(&check_these_hetmers, output, "bad_hetmers.csv");
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
        let result = extract_border(&test_vec);
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
        let result = extract_border(&test_vec);
        let expected = vec![
            "ATGAT".to_string(),
            "TTGTC".to_string(),
            "GGAAA".to_string(),
        ];
        assert_eq!(result, expected);
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
        let mut input = HashMap::new();
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
        let mut input = HashMap::new();
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
        let mut input = HashMap::new();
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
        let mut input = HashMap::new();
        input.insert(1, vec![0, 1]);
        input.insert(2, vec![2, 3, 4]); // should be filtered out
        input.insert(3, vec![2, 3, 4, 5]); // should be filtered out
        input.insert(4, vec![2]); // should be filtered out
        input.insert(5, vec![5, 6]);

        let alleles = 2;
        let result = filter_groups(input, alleles);

        let mut expected = HashMap::new();
        expected.insert(1, vec![0, 1]);
        expected.insert(5, vec![5, 6]);

        assert_eq!(result, expected);
    }

    #[test]
    fn three_alleles() {
        let mut input = HashMap::new();
        input.insert(1, vec![0, 1]);
        input.insert(2, vec![2, 3, 4]);
        input.insert(3, vec![2, 3, 4, 5]);
        input.insert(4, vec![2]);
        input.insert(5, vec![5, 6]);

        let alleles = 3;
        let result = filter_groups(input, alleles);

        let mut expected = HashMap::new();
        expected.insert(2, vec![2, 3, 4]);

        assert_eq!(result, expected);
    }

    #[test]
    fn four_alleles() {
        let mut input = HashMap::new();
        input.insert(1, vec![0, 1]);
        input.insert(2, vec![2, 3, 4]);
        input.insert(3, vec![2, 3, 4, 5]);
        input.insert(4, vec![2]);
        input.insert(5, vec![5, 6]);

        let alleles = 4;
        let result = filter_groups(input, alleles);

        let mut expected = HashMap::new();
        expected.insert(3, vec![2, 3, 4, 5]);

        assert_eq!(result, expected);
    }

    #[test]
    fn filter_groups_no_matches() {
        let mut input = HashMap::new();
        input.insert(1, vec![0]);
        input.insert(2, vec![1, 2, 3]);

        let result = filter_groups(input, 2);
        let expected: HashMap<u64, Vec<usize>> = HashMap::new();

        assert_eq!(result, expected);
    }
}
