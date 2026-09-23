# Contributing

## To-do

- [ ] add likelihood ratio test from HAWK paper?

- [ ] add simulate subcommand to simulate variants similar to mutatrix

- [x] add file does not exist errors for call subcommand: `common::ensure_readable`
  pre-flight checks all subcommand inputs and exits 1 with the missing path;
  dedup/call errors are propagated in main instead of being discarded with
  `let _ =`, and count/call output-write failures exit 1.

- [x] bug for variants within k of chromosome start? Why doesn't test data work anymore?

- [x] index: convert non-ATGC to N, and lowercase to uppercase

- [x] index: skip k-mers with non-ATGC chars

- [x] debug ref-dedup

- [x] include reference lengths in decision of what sequence to extract: start --- k --- REF length --- k --- end

- [x] how are N's handeled?

- [x] how are non ATGC-handeled?

- [x] index, if variant is within k bp of previous and next variant, also check the allele lengths, if the total length is > k bp, then you should still be able to extract k-mers?

- [x] debug ref-dedup: some variants with 0 reference allele-specific k-mers end up with a count of 1 at the end of ref-dedup, but don't actually have a reference-allele-specific k-mer in the index -> will probably lead to some allele frequencies being 0 instead of nan

- [x] add typical counting speed: release build, 18184-read fixture with a
  20933-k-mer index scans in ~6 s with 4 threads; see the memory-reduction
  entry below for a benchmark harness recipe.

- [x] reduce `count` memory: k-mers stored as packed u64 keys (k <= 31) instead of
  Strings; hashset shared across threads with Arc; counts map keyed by packed k-mer.
  Measured on a 6.2M k-mer index: 996 MB -> 501 MB with 1 thread, 2.49 GB -> 505 MB
  with 4 threads; counts byte-identical before/after.

- [x] stream ref-dedup in a single pass (like var-dedup): remove_ref_kmers no
  longer loads the whole index via read_index_kmers; each line is parsed,
  filtered against the reference hashset, and written immediately. RSS is
  flat (~250 MB) from 340 to 100k variants; output byte-identical.

- [x] add unit tests

- [x] add integration test (tests/pipeline.rs): runs the README quick-start
  against the tracked fixtures and asserts the canonical calls.txt and
  counts_by_kmer.txt md5s; runs in CI.

- [x] add cargo test + cargo clippy to CI (linter.yml test job)

- [x] write counts_by_kmer table in deterministic sorted order (rows sorted
  lexicographically by k-mer, stable across runs and thread counts); hetmers
  outputs still unordered, tracked in the cleanup item below.

- [ ] add methods to structs

- [x] cleanup hetmers subcommand: group_hashes/extract_hetmers now use a
  BTreeMap, so the six output files (seqs, counts, hashes, bad_hetmers,
  bayes_states, empirical_freqs) are in ascending-hash order and deterministic
  across runs. load_kmers rejects empty k-mer strings and unreadable files
  with clear errors; hetmers inputs are pre-flighted with ensure_readable.
  Remaining polish (separate to-dos): make load_kmers parse the whole table
  before running input checks, unify hetmers CLI list style (space-separated)
  with count (comma-separated), methods on structs.

- [ ] add q-mers?

- [x] debug ref-dedup: variants near chromosome ends need two regions extracted, region prior to the variant and the region between the variant and chromsome tip

- [x] hash info for variants so that you skip duplicates?

- [x] drop variants at tip of chromosome

- [x] trim beginning of sequence so that variants within k bp of chromosome start can be included

- [x] for ref-dedup, what to do if there are chromosomes in fasta that have no variants in vcf?

- [x] ref-dedup, grab chromosome end if there's no next record in vcf

- [x] for ref-dedup, put k-mers into hashset at end of each loop?

- [x] add verbose flag with clap, improve logging: https://rust-cli.github.io/book/tutorial/output.html

- [x] index step, for variants within k bp, only get k-mers that overlap one variant

- [x] skip invariant sites, sites where REF and ALT are the same or there is no ALT information

- [x] get k-mer length from index in ref-dedup

- [x] check that REF allele matches fasta file

- [x] add input checkers, check that vcf is sorted for freqk index

- [x] add filter or warning for when 0 or only a few k-mers tag a variant (if not many k-mers tag a variant and the variant is rare in the pool, then coverage needs to be absurdly high in order to accurately estimate allele frequency)

- [x] parallelize for loop over reads

- [x] remove leftover code

- [x] organize code into modules

- [x] use paraseq instead of fastqrs - paraseq has a conflict that prevents it from working with rusthtslib, also seems to do similarly to fastqrs in benchmarks

- [x] convert counts into allele frequencies

- [x] skip over variants at chromosome ends and overlapping variants

- [x] check that reference sequence in vcf matches fasta

- [x] var-dedup: add more print messages

- [x] bug: var-dedup, variant at end of one chromosome and start of another chromosome are labeled as overlapping and skipped

- [x] bug: ref-dedup, lots of variants are skipped

- [x] deduplicate: remove any allele-specific k-mers found elsewhere in the reference genome

- [x] dedeuplicate index: remove any k-mers that are shared across variants

- [x] determine k-mer length from index

- [x] when indexing, count number of k-mers for each allele

- [x] for count subcommand: load in k-mers from index as hash set, loop over input reads with a fastx iterator, slide window to get k-mers, only count k-mers if they're in the hash set
