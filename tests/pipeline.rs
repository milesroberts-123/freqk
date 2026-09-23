//! End-to-end pipeline test: runs the README quick-start commands against the
//! tracked test fixtures and asserts the canonical outputs, locking the
//! `calls.txt` md5 = ea5a558fc94766653a81d604755f04e4 contract.
//!
//! Skipped automatically when the gzip fixtures are absent (e.g. a fresh clone
//! without the force-added fastq files).

use std::process::Command;

const CALLS_MD5: &str = "ea5a558fc94766653a81d604755f04e4";
const COUNTS_BY_KMER_MD5: &str = "241b92bec1473f440842b9f129d5451e";

fn fixture(path: &str) -> String {
    format!("{}/tests/{}", env!("CARGO_MANIFEST_DIR"), path)
}

fn fastq_fixtures_available() -> bool {
    [
        "trimmed_paired_R1_1_0.fastq.gz",
        "trimmed_paired_R2_1_0.fastq.gz",
        "trimmed_unpaired_R1_1_0.fastq.gz",
        "trimmed_unpaired_R2_1_0.fastq.gz",
    ]
    .iter()
    .all(|f| std::path::Path::new(&fixture(f)).exists())
}

fn run(args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_freqk"))
        .args(args)
        .status()
        .expect("Could not run freqk binary");
    assert!(
        status.success(),
        "freqk {} failed with status {:?}",
        args.join(" "),
        status
    );
}

fn md5(path: &std::path::Path) -> String {
    let output = Command::new("md5sum")
        .arg(path)
        .output()
        .expect("md5sum failed (is it installed?)");
    let stdout = String::from_utf8_lossy(&output.stdout);
    stdout
        .split_whitespace()
        .next()
        .expect("md5sum produced no output")
        .to_string()
}

#[test]
fn quick_start_pipeline_produces_canonical_output() {
    if !fastq_fixtures_available() {
        eprintln!("Skipping: fastq fixtures not present (gitignored, force-added in repo)");
        return;
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let index = dir.path().join("index.txt");
    let var_index = dir.path().join("var_index.txt");
    let ref_index = dir.path().join("ref_index.txt");
    let counts_by_allele = dir.path().join("counts_by_allele.txt");
    let counts_by_kmer = dir.path().join("counts_by_kmer.txt");
    let calls = dir.path().join("calls.txt");

    run(&[
        "index",
        "-f",
        &fixture("test.fasta"),
        "--vcf",
        &fixture("test.vcf.gz"),
        "-o",
        index.to_str().unwrap(),
        "-k",
        "31",
    ]);
    run(&[
        "var-dedup",
        "--index",
        index.to_str().unwrap(),
        "-o",
        var_index.to_str().unwrap(),
    ]);
    run(&[
        "ref-dedup",
        "-i",
        var_index.to_str().unwrap(),
        "-o",
        ref_index.to_str().unwrap(),
        "-f",
        &fixture("test.fasta"),
        "--vcf",
        &fixture("test.vcf.gz"),
    ]);
    run(&[
        "count",
        "-i",
        ref_index.to_str().unwrap(),
        "-r",
        &fixture("trimmed_paired_R1_1_0.fastq.gz"),
        "-r",
        &fixture("trimmed_paired_R2_1_0.fastq.gz"),
        "-r",
        &fixture("trimmed_unpaired_R1_1_0.fastq.gz"),
        "-r",
        &fixture("trimmed_unpaired_R2_1_0.fastq.gz"),
        "-n",
        "4",
        "-f",
        counts_by_allele.to_str().unwrap(),
        "-c",
        counts_by_kmer.to_str().unwrap(),
    ]);
    run(&[
        "call",
        "-i",
        ref_index.to_str().unwrap(),
        "-c",
        counts_by_allele.to_str().unwrap(),
        "-o",
        calls.to_str().unwrap(),
    ]);

    assert_eq!(
        md5(&calls),
        CALLS_MD5,
        "calls.txt changed; check whether the difference is intended and update the contract md5"
    );
    assert_eq!(
        md5(&counts_by_kmer),
        COUNTS_BY_KMER_MD5,
        "counts_by_kmer.txt changed; check whether the difference is intended and update the contract md5"
    );
}
