// functions to check that input is properly formated
/// Check that k-mers are lexicographically sorted (first 1000 elements).
pub fn check_sort(seqs: &[String]) -> bool {
    let limit = seqs.len().min(1000);
    let seqs_sub = &seqs[..limit];

    let mut sorted_seqs = seqs_sub.to_vec();
    sorted_seqs.sort();

    let result = seqs_sub == sorted_seqs;

    println!("Input sorted: {}", result);
    result
}

/// Check that only ATGC are in the alphabet (first 1000 elements).
pub fn check_letters(seqs: &[String]) -> bool {
    let limit = seqs.len().min(1000);
    let seqs_sub = &seqs[..limit];

    let result = seqs_sub
        .iter()
        .all(|seq| seq.chars().all(|c| matches!(c, 'A' | 'T' | 'G' | 'C')));
    println!("Only ATGC: {}", result);
    result
}

/// Run all input checks, panicking on failure.
pub fn all_checks(seqs: &[String]) {
    println!("Checking input format...");
    if !check_sort(seqs) {
        panic!(":(");
    }

    if !check_letters(seqs) {
        panic!(":(")
    }
}
