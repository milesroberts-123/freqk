# Contributing

## To-do

- [ ] add likelihood ratio test from HAWK paper?

- [ ] add simulate subcommand to simulate variants similar to mutatrix

- [ ] add file does not exist errors for call subcommand

- [x] bug for variants within k of chromosome start? Why doesn't test data work anymore?

- [x] index: convert non-ATGC to N, and lowercase to uppercase

- [x] index: skip k-mers with non-ATGC chars

- [x] debug ref-dedup

- [x] include reference lengths in decision of what sequence to extract: start --- k --- REF length --- k --- end

- [x] how are N's handeled?

- [x] how are non ATGC-handeled?

- [x] index, if variant is within k bp of previous and next variant, also check the allele lengths, if the total length is > k bp, then you should still be able to extract k-mers?

- [x] debug ref-dedup: some variants with 0 reference allele-specific k-mers end up with a count of 1 at the end of ref-dedup, but don't actually have a reference-allele-specific k-mer in the index -> will probably lead to some allele frequencies being 0 instead of nan

- [ ] add typical counting speed

- [x] add unit tests

- [ ] add methods to structs

- [ ] cleanup hetmers subcommand

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
