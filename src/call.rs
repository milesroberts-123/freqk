use crate::common;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};

/// Write a 2D array of allele frequencies to a file, one `|`-joined line per site.
fn write_allele_freqs(counts: Vec<Vec<f32>>, output_file: &str) {
    let mut output_str = Vec::new();
    for one_site in &counts {
        let one_site_str = one_site
            .iter()
            .map(|num| num.to_string())
            .collect::<Vec<String>>()
            .join("|");
        output_str.push(one_site_str);
    }
    let _ = common::write_strings(output_str, output_file);
}

/// Convert a 2D array of strings to u32, dropping unparseable entries.
fn vec_vec_string_to_vec_vec_u32(string_vec_vec: Vec<Vec<String>>) -> Vec<Vec<u32>> {
    string_vec_vec
        .into_iter()
        .map(|inner_vec_strings| {
            inner_vec_strings
                .into_iter()
                .filter_map(|s| s.parse::<u32>().ok())
                .collect()
        })
        .collect()
}

/// Divide two 2D arrays elementwise, returning NaN where the divisor is zero.
fn elementwise_division_2d(
    vec_a: &[Vec<u32>],
    vec_b: &[Vec<u32>],
) -> Result<Vec<Vec<f32>>, io::Error> {
    if vec_a.len() != vec_b.len() {
        log::error!("Vectors of different lengths");
    }

    let result: Vec<Vec<f32>> = vec_a
        .iter()
        .zip(vec_b.iter())
        .map(|(inner_a, inner_b)| {
            if inner_a.len() != inner_b.len() {
                log::error!("Inner vectors of different lengths");
            }
            inner_a
                .iter()
                .zip(inner_b.iter())
                .map(|(&a, &b)| {
                    if b == 0 {
                        log::warn!("Division by zero detected because there are zero allele-specific k-mers for a particular variant (frequency estimate will be NaN).");
                        f32::NAN
                    } else {
                        (a as f32)/(b as f32)
                    }
                })
                .collect()
        })
        .collect();

    Ok(result)
}

/// Normalize each row of counts by its sum to get allele frequencies.
fn normalized_counts_to_allele_freq(norm_counts: Vec<Vec<f32>>) -> Vec<Vec<f32>> {
    let mut sums: Vec<f32> = Vec::new();
    for inner_vec in &norm_counts {
        sums.push(inner_vec.iter().sum());
    }
    let mut result = Vec::new();
    for (norm_count, sum) in norm_counts.iter().zip(sums.iter()) {
        let mut allele_freqs = Vec::new();
        for count in norm_count {
            allele_freqs.push(count / sum)
        }
        result.push(allele_freqs);
    }
    result
}

/// Convert counts by allele into allele frequencies and write them to a file.
pub fn call_from_counts(index: &str, counts: &str, output: &str) -> Result<(), io::Error> {
    log::info!("Reading index for number of unique k-mers per allele...");
    let num_uniq_kmers_per_allele = common::read_index_field(index, 6);
    log::info!("Converting data to u32...");
    let num_uniq_kmers_per_allele_u32 = vec_vec_string_to_vec_vec_u32(num_uniq_kmers_per_allele?);
    log::info!("Parsing counts by allele...");
    let file = File::open(counts)?;
    let reader = BufReader::new(file);
    let mut kmer_counts_per_allele = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<String> = line.split('|').map(|s| s.to_owned()).collect();
        log::debug!("Parsed counts by allele: {:?}", fields);
        kmer_counts_per_allele.push(fields);
    }
    log::info!("Converting data to u32...");
    let kmer_counts_per_allele_u32 = vec_vec_string_to_vec_vec_u32(kmer_counts_per_allele);
    log::info!("Normalizing k-mer counts by number of allele-specific k-mers...");
    let counts_per_kmer =
        elementwise_division_2d(&kmer_counts_per_allele_u32, &num_uniq_kmers_per_allele_u32);
    log::info!("Converting normalized counts to allele frequencies...");
    let allele_freq = normalized_counts_to_allele_freq(counts_per_kmer?);
    log::info!("Writing allele frequency estimates to {}", output);
    write_allele_freqs(allele_freq, output);
    Ok(())
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    #[test]
    fn test_elementwise_division_2d() {
        let test_kmer_counts_per_allele = vec![vec![55, 11], vec![22, 10]];
        let test_num_uniq_kmers_per_allele = vec![vec![10, 4], vec![2, 20]];
        let expected = vec![vec![5.5, 2.75], vec![11.0, 0.5]];
        let result = elementwise_division_2d(
            &test_kmer_counts_per_allele,
            &test_num_uniq_kmers_per_allele,
        );
        assert_eq!(result.unwrap(), expected);
    }
}
