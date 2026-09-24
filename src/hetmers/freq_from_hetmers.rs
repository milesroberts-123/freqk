/// Parse comma-joined count pairs into (minor, major) count tuples.
/// Pairs that do not have exactly two parseable counts become `None`.
fn parse_count_pairs(count_pairs: &[String]) -> Vec<Option<(f64, f64)>> {
    count_pairs
        .iter()
        .map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                if let (Ok(num1), Ok(num2)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                    return Some((num1.min(num2), num1.max(num2)));
                }
            }
            None
        })
        .collect()
}

/// Compute empirical minor-allele frequencies from hetmer count pairs.
/// The output is parallel to the input (and to seqs.csv/counts.csv): rows
/// with unparseable counts or a zero total are written as NA instead of
/// being dropped, so row i of every output file describes the same hetmer.
pub fn counts_to_frequencies(count_pairs: &[String]) -> Vec<String> {
    log::info!("Calculating frequencies...");
    parse_count_pairs(count_pairs)
        .into_iter()
        .map(|pair| match pair {
            Some((min_num, max_num)) => {
                let sum = min_num + max_num;
                if sum != 0.0 {
                    (min_num / sum).to_string()
                } else {
                    "NA".to_string()
                }
            }
            None => "NA".to_string(),
        })
        .collect()
}

/// Tag hetmers with really high total coverage (potentially due to paralogous sequences).
pub fn high_cov_hetmers(count_pairs: &[String], sigma: f64, n: i32, cov: f64) -> Vec<String> {
    log::info!("Checking for questionable hetmers...");
    let stderr = ((n as f64) * cov).sqrt();
    parse_count_pairs(count_pairs)
        .into_iter()
        .filter_map(|pair| {
            let (num1, num2) = pair?;
            if num1 + num2 > sigma * stderr {
                Some("1".to_string())
            } else {
                Some("0".to_string())
            }
        })
        .collect()
}

/// Compute the factorial of a non-negative integer.
fn factorial(x: usize) -> f64 {
    (1..=x).map(|i| i as f64).product()
}

/// Compute the truncation constant for a Poisson distribution.
pub fn truncation_constant(c: usize, lambda: f64) -> f64 {
    let sum: f64 = (0..c)
        .map(|x| {
            let numerator = (-lambda).exp() * lambda.powi(x as i32);
            let denominator = factorial(x); // == gamma(x + 1)
            numerator / denominator
        })
        .sum();

    1.0 - sum
}

/// Compute the posterior mode of the minor k-mer count.
pub fn posterior_min_kmer_count(
    x: f64,
    z: f64,
    n: i32,
    cov: f64,
    c: usize,
    alpha: f64,
    beta: f64,
) -> Result<usize, std::io::Error> {
    let mut likelihood_times_prior = Vec::new();
    let mut total_probability = 0.0;
    let max_minor_count = n / 2; // minor allele can't have frequency above 1/2 by definition

    for i in 1..max_minor_count {
        let p = i as f64 / n as f64;
        let lambda_x = (i as f64) * cov;
        let lambda_y = ((n - i) as f64) * cov;

        let tx = truncation_constant(c, lambda_x);
        let ty = truncation_constant(c, lambda_y);

        let likelihood = p.powf(x + alpha - 1.0) * (1.0 - p).powf((z - x) + beta - 1.0);
        let likelihood_truncated = likelihood / (tx * ty);
        likelihood_times_prior.push(likelihood_truncated);
        total_probability += likelihood;
    }

    let posterior: Vec<f64> = likelihood_times_prior
        .iter()
        .map(|val| val / total_probability)
        .collect();

    let max_index = posterior
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(idx, _)| idx)
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Cannot compute the Bayesian allele state for a hetmer with pool size {} (the minor-allele frequency grid is empty for pools smaller than 2)",
                    n
                ),
            )
        })?;

    Ok(max_index + c)
}

/// Calculate the posterior distribution for allele count of each hetmer.
/// The output is parallel to the input (and to seqs.csv/counts.csv): rows
/// with unparseable counts are written as 0.
pub fn counts_to_bayes_state(
    count_pairs: &[String],
    n: i32,
    cov: f64,
    c: usize,
    alpha: f64,
    beta: f64,
) -> Result<Vec<usize>, std::io::Error> {
    log::info!("Calculating posterior...");
    parse_count_pairs(count_pairs)
        .into_iter()
        .map(|pair| match pair {
            Some((x, max_num)) => {
                let z = x + max_num;
                posterior_min_kmer_count(x, z, n, cov, c, alpha, beta)
            }
            None => Ok(0),
        })
        .collect()
}

// test functions
#[cfg(test)]
mod unit_tests {
    use super::*;

    #[test]
    fn a_few_freqs() {
        let count_pairs = vec![
            "5,120".to_string(),
            "1,9".to_string(),
            "20,140".to_string(),
            "22,22".to_string(),
        ];
        let result = counts_to_frequencies(&count_pairs);
        let expected = vec![
            "0.04".to_string(),
            "0.1".to_string(),
            "0.125".to_string(),
            "0.5".to_string(),
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn freqs_parallel_with_na_for_bad_pairs() {
        // Malformed and zero-total pairs become NA in place; the output stays
        // parallel to seqs.csv/counts.csv.
        let count_pairs = vec![
            "5,120".to_string(),
            "garbage".to_string(),
            "0,0".to_string(),
            "1,3".to_string(),
        ];
        let result = counts_to_frequencies(&count_pairs);
        assert_eq!(
            result,
            vec![
                "0.04".to_string(),
                "NA".to_string(),
                "NA".to_string(),
                "0.25".to_string()
            ]
        );
        assert_eq!(result.len(), count_pairs.len());
    }

    #[test]
    fn bayes_state_rejects_small_pool() {
        // Pool < 2 has an empty minor-frequency grid; must error, not panic.
        let count_pairs = vec!["5,120".to_string()];
        let result = counts_to_bayes_state(&count_pairs, 1, 50.0, 2, 0.05, 0.05);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("pool size 1"));
    }

    #[test]
    fn bayes_state_none_pairs_are_zero() {
        let count_pairs = vec!["garbage".to_string(), "5,120".to_string()];
        let result = counts_to_bayes_state(&count_pairs, 20, 50.0, 2, 0.05, 0.05).unwrap();
        assert_eq!(result[0], 0);
        assert!(result[1] > 0);
        assert_eq!(result.len(), 2);
    }

    // helper function to test equality of floating point numbers
    fn round_to_decimals(x: f64, decimals: u32) -> f64 {
        let factor = 10f64.powi(decimals as i32);
        (x * factor).round() / factor
    }

    #[test]
    fn small_truncation() {
        let result = round_to_decimals(truncation_constant(5, 10.0), 7);
        let expected = 0.9707473;
        assert_eq!(result, expected);
    }
}
