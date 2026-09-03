# AGENTS.md

freqk: Rust CLI estimating allele frequencies in pool-seq data from allele-specific
k-mer counts. Single binary, clap-derive subcommands: index, var-dedup, ref-dedup,
count, call, hetmers.

## Commands
- `cargo build` — debug binary at `target/debug/freqk`
- `cargo test` — all tests are inline `#[cfg(test)] mod unit_tests` blocks at the
  bottom of each `src/*.rs` (there is no Rust integration suite in `tests/`;
  `tests/` holds data fixtures only)
- `cargo fmt` before committing — CI runs super-linter with rustfmt checks for
  Rust 2021/2024 (`.github/workflows/linter.yml`)

## End-to-end verification (README quick start)
```bash
target/debug/freqk index -f tests/test.fasta --vcf tests/test.vcf.gz -o index.txt -k 31
target/debug/freqk var-dedup --index index.txt -o var_index.txt
target/debug/freqk ref-dedup -i var_index.txt -o ref_index.txt -f tests/test.fasta --vcf tests/test.vcf.gz
target/debug/freqk count -i ref_index.txt -r tests/test.fastq.gz -n 4 -f counts_by_allele.txt -c counts_by_kmer.txt
target/debug/freqk call -i ref_index.txt -c counts_by_allele.txt -o calls.txt
```
- Success check: `md5sum calls.txt` == `ea5a558fc94766653a81d604755f04e4`
- Gotcha: `tests/test.fastq.gz` is NOT in git (`.gitignore` ignores `*.fastq.gz`);
  only `test.fasta[.fai]` and `test.vcf.gz[.tbi]` are tracked. Generate/obtain a
  fastq to run the `count` step.
- Intermediate outputs (`*.txt`) are gitignored; don't commit them.

## Input requirements
- FASTA must be `samtools faidx`-indexed (`.fai` alongside it).
- VCF must be sorted, bgzipped, tabix-indexed, and normalized so variants at the
  same position are one multiallelic record (`bcftools norm -m +any`).
- Chromosome names must match between FASTA and VCF.
- `-k` is only needed for `index`; downstream commands read k from the index.

## Architecture
- `main.rs` wires clap subcommands to workflow fns in sibling modules:
  `index.rs` (build allele-specific k-mer index), `dedup.rs` (var-dedup across
  variants; ref-dedup vs reference), `count.rs` (k-mer counting over fastq with
  manual `nthreads` parallelism), `call.rs` (normalize counts → allele freqs),
  `hetmers.rs` + `hetmers/freq_from_hetmers.rs`, `common.rs` (shared IO/parsing).
- Formats: index is CSV with pipe-separated fields inside (alleles, k-mers,
  counts); `call` output is one line per variant, pipe-separated allele freqs.

## Conventions / gotchas
- Work happens on `dev`; changes reach `main` via merge PRs (matches git history).
- `CONTRIBUTING.md` is a running to-do list; consult it and check off items done.
- Don't swap the `fastq` crate for paraseq — known conflict with rust-htslib
  (documented in CONTRIBUTING).
- Every subcommand has `-v` verbosity flags (`-vvv` = debug); README asks bug
  reports to include `-vvv` output.
- `index` handles edge cases already: uppercases, converts non-ATGC to N, skips
  non-ATGC k-mers, skips IUPAC-mismatch and chromosome-tip variants. Extend the
  existing unit tests for these paths instead of restructuring.