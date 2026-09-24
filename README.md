# freqk

Estimate frequencies of known variants in pool-seq data from k-mer counts

Author: Miles Roberts

Contact: leave a github issue

**Note:** Before leaving an issue, re-run the problematic command(s) with `-vvv` to get full debug information for that command and post the debug information with your issue.

## Installation

Grab rust binary from release page

## Inputs

1. FASTA file of reference sequence/path:

* chromosome names that match VCF file

* indexed with `samtools faidx`

2. VCF file of variants: 

* REF alleles correspond to sequences in FASTA file 

* chromosome names that match FASTA file

* sorted

* bgzipped `bgzip myfile.vcf`

* normalized (variants at same position are represented as one multiallelic record) `bcftools norm -m +any`

* Index of vcf file `tabix myfile.vcf.gz`

3. Pooled DNA sequencing reads. Read pairing information does not matter and reads can be passed as a comma-separated list of fastq files (gzipped or not), like so:

```
# for gzip compressed fastq files
-r r1.fastq.gz,r2.fastq.gz,u.fastq.gz

# for uncompressed fastq files
-r r1.fastq,r2.fastq,u.fastq
```

## Outputs

Estimates of allele frequencies for variants in the vcf file. These are output in a text file with pipe separators

```
0.01|0.99
1|0
0.75|0.25
0.3|0.5|0.2
```

### Reading outputs into R for further analysis

If the VCF file includes variable number of alleles per site (i.e. not only bi-alleleic sites), then the number of entries per line varies. If you want to load this file in R. You can do something like the following:

Merge the index with the allele frequncy estimates

## Quick start

### short version

```bash
# download repo
git clone https://github.com/milesroberts-123/freqk.git
cd freqk
# try freqk on small test dataset
freqk index -f tests/test.fasta --vcf tests/test.vcf.gz -o index.txt -k 31
freqk var-dedup --index index.txt --output var_index.txt
freqk ref-dedup -i var_index.txt -o ref_index.txt -f tests/test.fasta --vcf tests/test.vcf.gz
freqk count -i ref_index.txt -r tests/trimmed_paired_R1_1_0.fastq.gz,tests/trimmed_paired_R2_1_0.fastq.gz,tests/trimmed_unpaired_R1_1_0.fastq.gz,tests/trimmed_unpaired_R2_1_0.fastq.gz -n 4 -f counts_by_allele.txt -c counts_by_kmer.txt
freqk call -i ref_index.txt -c counts_by_allele.txt -o calls.txt
```

If the example was successful, the `md5sum` of the final calls.txt file should be `27c8de92b97963229c7cc8e4c3569ad6`

### step-by-step breakdown

1. Index the panel of reference variants

`freqk index -f tests/test.fasta --vcf tests/test.vcf.gz -o index.txt -k 31`

2. (Optional but recommended) Deduplicate index

This step can be done in 2 stages. You can do one, both, or neither stages in any order, but generally doing both before proceeding to the counting step (3) is strongly recommended for the most rigorous results.

The faster of the two deduplication steps is `var-dedup`. This simply scans through the index and removes any putatively allele-specific k-mers that are actually found in multiple alleles of different variants in the index. For example:

`freqk var-dedup --index index.txt --output var_index.txt`

The slower deduplication step is `ref-dedup`. This step removes any allele-specific k-mers that are found elsewhere in the reference sequence. 

`freqk ref-dedup -i var_index.txt -o ref_index.txt -f tests/test.fasta --vcf tests/test.vcf.gz`

3. Count the indexed k-mers in the pool-seq reads (multithreaded). 

For example, counting indexed k-mers with four threads (`-n 4`) looks like this:

`freqk count -i ref_index.txt -r tests/trimmed_paired_R1_1_0.fastq.gz,tests/trimmed_paired_R2_1_0.fastq.gz,tests/trimmed_unpaired_R1_1_0.fastq.gz,tests/trimmed_unpaired_R2_1_0.fastq.gz -n 4 -f counts_by_allele.txt -c counts_by_kmer.txt`

4. Normalize allele-specific k-mer counts into allele frequencies

This step just divides the counts of allele-specific k-mers in the reads by the number of allele-specific k-mers in the index.

`freqk call -i ref_index.txt -c counts_by_allele.txt -o calls.txt`

## Command line interfaces

All commands have a verbosity flag. Only errors are output by default, but adding `-v` will make warnings print, `-vv` means info will also print, and `-vvv` means debug data will also print.

### help

```bash
$ freqk help
Usage: freqk [OPTIONS] <COMMAND>

Commands:
  index      Get k-mers specific to each allele of each variant
  var-dedup  Deduplicate index of k-mers shared across variants
  count      Count k-mers by allele
  call       Convert counts by allele into allele frequencies
  ref-dedup  Deduplicate index of reference k-mers
  filter     Filter index rows by allele-specific k-mer content
  hetmers    Find het-mers in a k-mer count table (e.g. from kmc or jellyfish)
  help       Print this message or the help of the given subcommand(s)

Options:
  -v, --verbose...  Increase logging verbosity
  -q, --quiet...    Decrease logging verbosity
  -h, --help        Print help
  -V, --version     Print version
```

### index

```bash
$ freqk index -h
Get k-mers specific to each allele of each variant

Usage: freqk index [OPTIONS] --fasta <FASTA> --vcf <VCF> --output <OUTPUT> --kmer <KMER>

Options:
  -f, --fasta <FASTA>    fasta file of reference genome
      --vcf <VCF>        vcf file of variations between reference and other genomes
  -o, --output <OUTPUT>  name of the index file to be output
  -k, --kmer <KMER>      kmer length for building the index
  -v, --verbose...       Increase logging verbosity
  -q, --quiet...         Decrease logging verbosity
  -h, --help             Print help
```

### filter

Filter index rows by allele-specific k-mer content. Works on any index
(`var-dedup` or `ref-dedup` output), so thresholds can be tuned without
re-running the expensive `index` step. Rows are dropped whole: a frequency
estimate is only meaningful when all alleles are tagged.

```bash
$ freqk filter -h
Filter index rows by allele-specific k-mer content

Usage: freqk filter [OPTIONS] --index <INDEX> --output <OUTPUT>

Options:
  -i, --index <INDEX>
          path to index file
  -o, --output <OUTPUT>
          path to filtered index file
  -m, --min-alleles <MIN_ALLELES>
          Keep index rows only if at least this many alleles have >= 1 allele-specific k-mer (0 keeps all rows) [default: 0]
  -k, --min-kmers-per-allele <MIN_KMERS_PER_ALLELE>
          Keep index rows only if every allele has >= this many allele-specific k-mers (0 keeps all rows) [default: 0]
  -v, --verbose...
          Increase logging verbosity
  -q, --quiet...
          Decrease logging verbosity
  -h, --help
          Print help
```

Example: keep variants where both alleles are tagged by at least 3
allele-specific k-mers:

```bash
freqk filter -i ref_index.txt -o filtered_index.txt -k 3
```

### hetmers

Find *het-mers*: pairs of k-mers whose "borders" (the k-mer with its central
base removed) are identical, meaning the pair differs at exactly one base —
a heterozygous SNV. hetmers discovers these sites **de novo** from k-mer
counts alone, with no VCF or reference genome, so it is a parallel branch of
the workflow rather than a stage in the index → count → call pipeline (which
estimates frequencies at *known* variants).

The input is a k-mer count table with two tab-separated columns (k-mer,
count), produced by any k-mer counter such as [kmc](https://github.com/refresh-bio/KMC)
or [jellyfish](https://github.com/gmarcais/Jellyfish) (or `freqk count -c`,
which is already sorted as hetmers requires). Counts below the minimum are
filtered out; the table must be lexicographically sorted by k-mer.

For each hetmer (group of exactly `-l/--alleles` k-mers sharing one border),
six output files are written (`<prefix>_*.csv`, all in ascending-hash order):
`seqs` (comma-joined k-mers), `counts` (their counts), `hashes` (the border
hash), `empirical_freqs` (minor-allele count fraction; `NA` when counts are
unparseable or sum to zero), `bayes_states` (posterior mode of the minor
k-mer count, given pool size, coverage, and the alpha/beta prior), and
`bad_hetmers` (1 = total coverage above `sigma * sqrt(pool * coverage)`,
suggesting paralogous sequence).

```bash
$ freqk hetmers -h
Find het-mers in a k-mer count table (e.g. from kmc or jellyfish)

Usage: freqk hetmers [OPTIONS] --inputs <INPUTS>... --outputs <OUTPUTS>... --minimums <MINIMUMS>... --coverages <COVERAGES>... --pools <POOLS>... --alphas <ALPHAS>... --betas <BETAS>... --sigmas <SIGMAS>...

Options:
  -i, --inputs <INPUTS>...        comma-separated list of k-mer count tables (two tab-separated columns: k-mer, count)
  -o, --outputs <OUTPUTS>...      comma-separated list of output file prefixes
  -m, --minimums <MINIMUMS>...    comma-separated list of minimum k-mer counts
  -l, --alleles <ALLELES>         number of alleles in each hetmer [default: 2]
  -c, --coverages <COVERAGES>...  comma-separated list of mean k-mer coverages
  -p, --pools <POOLS>...          comma-separated list of pool sizes
  -a, --alphas <ALPHAS>...        comma-separated list of alpha shape parameters
  -b, --betas <BETAS>...          comma-separated list of beta shape parameters
  -s, --sigmas <SIGMAS>...        comma-separated list of sigma thresholds
  -v, --verbose...
          Increase logging verbosity
  -q, --quiet...
          Decrease logging verbosity
  -h, --help
          Print help
```

Example: find het-mers in a k-mer count table from a diploid pool of 20
samples at 50x mean coverage:

```bash
freqk hetmers -i kmer_counts.txt -o hetmers -m 2 -c 50 -p 20 -a 0.05 -b 0.05 -s 0.9
```
