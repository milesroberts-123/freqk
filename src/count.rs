use fastq::{parse_path, Record};
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use std::sync::Arc;

use crate::common::{
    self, increment, new_sharded_counts, prefill_from_index, shard_of, ShardedCounts,
};

/// Count indexed k-mers in one reads file, incrementing the shared sharded
/// count map in place. Every worker thread looks up its k-mers in the same
/// pre-filled map, so there are no per-thread maps and no merge step.
/// Progress is logged every `print_frequency` reads (per worker thread).
pub fn count_target_kmers_in_reads(
    reads: &str,
    k: usize,
    nthreads: usize,
    print_frequency: usize,
    sharded_counts: Arc<ShardedCounts>,
) {
    parse_path(Some(reads), |parser| {
        parser
            .parallel_each::<(), (), _>(nthreads, move |record_sets| {
                let mut num_records = 0;
                // Memo of the last incremented key: long identical-k-mer runs
                // (poly-A tails, repeated adaptors) would otherwise re-lock
                // the same shard for every read.
                let mut last_kmer: Option<u64> = None;
                for record_set in record_sets {
                    for record in record_set.iter() {
                        if num_records % print_frequency == 0 {
                            log::info!("Reads processed: {}", num_records);
                        }
                        num_records += 1;
                        let read_kmers = common::get_canonical_kmers_packed(
                            std::str::from_utf8(record.seq())
                                .expect("Invalid UTF-8 in read sequence"),
                            k,
                        );
                        for read_kmer in read_kmers {
                            if last_kmer == Some(read_kmer) {
                                continue;
                            }
                            increment(&sharded_counts, read_kmer);
                            last_kmer = Some(read_kmer);
                        }
                    }
                }
            })
            .expect("Invalid fastq file");
    })
    .expect("Invalid compression");
}

/// Count indexed k-mers across several reads files into one shared map.
/// All files increment the same sharded tables, so accumulation is free.
pub fn count_target_kmers_in_reads_files(
    index: &str,
    reads_files: &[String],
    k: usize,
    nthreads: usize,
    print_frequency: usize,
) -> Arc<ShardedCounts> {
    let sharded_counts = new_sharded_counts();
    log::info!("Pre-filling shared count map from index...");
    if let Err(e) = prefill_from_index(index, &sharded_counts) {
        log::error!("Loading index at {} failed: {}", index, e);
        std::process::exit(1);
    }
    for reads in reads_files {
        log::info!("Counting k-mers in reads file: {}", reads);
        count_target_kmers_in_reads(reads, k, nthreads, print_frequency, sharded_counts.clone());
        log::info!(
            "Accumulated counts: {} unique k-mers (after {})",
            common::sharded_len(&sharded_counts),
            reads
        );
    }
    sharded_counts
}

/// Write k-mer counts to a file. All non-zero entries are collected into one
/// table and sorted numerically by packed key; since 2-bit codes A<C<G<T are
/// monotone, numeric key order equals lexicographic k-mer order, so the
/// output is the same fully sorted table as before without materializing any
/// k-mer strings for the sort. Zero-count k-mers (pre-filled, never observed)
/// are skipped. Shards are hash-scattered, so a global sort is required for
/// ordering (per-shard order alone is not global).
pub fn write_kmers(kmer_counts: &ShardedCounts, k: usize, output: &str) -> io::Result<()> {
    let mut file = File::create(output)?;
    let mut rows: Vec<(u64, usize)> = Vec::new();
    for shard in kmer_counts.iter() {
        let shard = shard.lock().unwrap();
        rows.extend(
            shard
                .iter()
                .filter(|&(_, &count)| count > 0)
                .map(|(&key, &count)| (key, count)),
        );
    }
    rows.sort_unstable();
    for (key, count) in rows {
        writeln!(file, "{}\t{}", common::unpack_kmer(key, k), count)?;
    }
    Ok(())
}

/// Sum k-mer counts into totals per allele, one line per index entry.
pub fn combine_counts_by_allele(
    index: &str,
    counts: &ShardedCounts,
) -> Result<Vec<String>, io::Error> {
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    let mut result = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 8 {
            continue;
        }
        let kmers = fields[7];
        let kmers_list: Vec<&str> = kmers.split('|').collect();

        let kmers_by_allele: Vec<Vec<&str>> =
            kmers_list.iter().map(|s| s.split(';').collect()).collect();

        let mut counts_by_allele = Vec::new();

        for allele in kmers_by_allele {
            let mut total_allele_count: usize = 0;
            for kmer in allele {
                if let Some(packed) = common::pack_kmer(kmer) {
                    let shard = counts[shard_of(packed)].lock().unwrap();
                    if let Some(kmer_count) = shard.get(&packed) {
                        total_allele_count += *kmer_count;
                    }
                }
            }
            counts_by_allele.push(total_allele_count);
        }

        let counts_by_allele_str = counts_by_allele
            .iter()
            .map(|num| num.to_string())
            .collect::<Vec<String>>()
            .join("|");

        result.push(counts_by_allele_str);
    }
    Ok(result)
}

/// Count indexed k-mers in reads and write per-allele and per-k-mer counts.
/// When `count_output` is None, the per-k-mer count table (and its sort) is
/// skipped entirely.
pub fn count_workflow(
    index: &str,
    reads_files: &[String],
    nthreads: usize,
    print_frequency: usize,
    freq_output: &str,
    count_output: Option<&str>,
) {
    log::info!("Loading index at {} ...", index);
    let k = common::k_from_index(index);
    log::info!("k is: {:?}", k);
    log::debug!("Counting indexed k-mers in reads...");
    let k = k.expect("Cannot parse kmer length from index.") as usize;
    let kmer_counts =
        count_target_kmers_in_reads_files(index, reads_files, k, nthreads, print_frequency);
    match count_output {
        Some(count_output) => {
            if let Err(e) = write_kmers(&kmer_counts, k, count_output) {
                log::error!("Writing k-mer counts to {} failed: {}", count_output, e);
                std::process::exit(1);
            }
        }
        None => log::debug!("No -c output given, skipping per-k-mer count table"),
    }
    log::debug!("Combining k-mer counts by allele...");
    let counts_by_allele = combine_counts_by_allele(index, &kmer_counts);
    log::debug!("Writing counts by allele...");
    match counts_by_allele {
        Ok(strings) => {
            if let Err(e) = common::write_strings(strings, freq_output) {
                log::error!("Writing counts by allele to {} failed: {}", freq_output, e);
                std::process::exit(1);
            }
        }
        Err(e) => {
            log::error!("Combining counts by allele failed: {}", e);
            std::process::exit(1);
        }
    }
    log::debug!("Successfully wrote counts by allele!");
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    /// Write a tiny index and return its path string.
    fn tiny_index(dir: &std::path::Path, content: &str) -> String {
        let path = dir.join("index.txt");
        std::fs::write(&path, content).unwrap();
        path.to_str().unwrap().to_string()
    }

    const INDEX_ROWS: &str = concat!(
        "0,1,100,AAAAA,REF|ALT,AAAAA|ACAAA,0|1,|\n",
        "1,1,200,CCCCC,REF|ALT,CCCCC|CCGGG,1|1,AAACA;CCCCT|GGGGG\n",
    );

    #[test]
    fn test_shard_of_spreads_keys() {
        // Consecutive packed k-mers must not all land on one shard.
        let shards: std::collections::HashSet<usize> = (0..1000).map(|i| shard_of(i)).collect();
        assert!(shards.len() > 100);
    }

    #[test]
    fn test_sharded_order_equals_lexicographic() {
        // write_kmers must output globally lexicographically sorted k-mers:
        // a single numeric sort over hash-scattered shards is equivalent to
        // string sort, because 2-bit codes are monotone (A<C<G<T). Only
        // k-mers pre-filled from the index are counted; unknown keys are
        // dropped by increment().
        let kmers = ["TTTTT", "AAACA", "GCGAT", "ACGTA", "GGGGG", "CCCCC"];
        let dir = tempfile::tempdir().unwrap();
        let index = tiny_index(
            dir.path(),
            concat!(
                "0,1,100,AAAAA,REF|ALT,AAAAA|ACAAA,0|1,",
                "TTTTT;GCGAT|ACGTA;CCCCC;AAACA;GGGGG\n",
            ),
        );
        let sharded = new_sharded_counts();
        prefill_from_index(&index, &sharded).unwrap();
        for kmer in kmers {
            increment(&sharded, common::pack_kmer(kmer).unwrap());
        }
        let out = dir.path().join("counts.txt");
        write_kmers(&sharded, 5, out.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&out).unwrap();
        let got: Vec<&str> = content
            .lines()
            .map(|l| l.split('\t').next().unwrap())
            .collect();
        let mut sorted = kmers.to_vec();
        sorted.sort();
        assert_eq!(got, sorted);
    }

    #[test]
    fn test_prefill_and_increment() {
        let dir = tempfile::tempdir().unwrap();
        let index = tiny_index(dir.path(), INDEX_ROWS);
        let sharded = new_sharded_counts();
        prefill_from_index(&index, &sharded).unwrap();
        // Every ATGC k-mer of the index exists at count 0.
        for kmer in ["AAACA", "CCCCT", "GGGGG"] {
            let packed = common::pack_kmer(kmer).unwrap();
            let shard = sharded[shard_of(packed)].lock().unwrap();
            assert_eq!(shard.get(&packed), Some(&0), "{}", kmer);
        }
        // Increment flips observed counts; unknown keys stay absent.
        increment(&sharded, common::pack_kmer("AAACA").unwrap());
        increment(&sharded, common::pack_kmer("AAACA").unwrap());
        let packed = common::pack_kmer("AAACA").unwrap();
        let shard = sharded[shard_of(packed)].lock().unwrap();
        assert_eq!(shard.get(&packed), Some(&2));
    }

    #[test]
    fn test_write_kmers_skips_zero_and_sorts() {
        let dir = tempfile::tempdir().unwrap();
        let index = tiny_index(dir.path(), INDEX_ROWS);
        let sharded = new_sharded_counts();
        prefill_from_index(&index, &sharded).unwrap();
        increment(&sharded, common::pack_kmer("GGGGG").unwrap());
        increment(&sharded, common::pack_kmer("AAACA").unwrap());
        let out = dir.path().join("counts.txt");
        write_kmers(&sharded, 5, out.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&out).unwrap();
        assert_eq!(content, "AAACA\t1\nGGGGG\t1\n");
    }

    #[test]
    fn test_combine_counts_by_allele() {
        let dir = tempfile::tempdir().unwrap();
        let index = tiny_index(dir.path(), INDEX_ROWS);
        let sharded = new_sharded_counts();
        prefill_from_index(&index, &sharded).unwrap();
        increment(&sharded, common::pack_kmer("AAACA").unwrap());
        increment(&sharded, common::pack_kmer("GGGGG").unwrap());
        let counts = combine_counts_by_allele(&index, &sharded).unwrap();
        // Line 1: empty allele 0, empty allele 1 -> "0|0".
        // Line 2: allele 0 = AAACA + CCCCT, allele 1 = GGGGG -> "1|1".
        assert_eq!(counts, vec!["0|0", "1|1"]);
    }

    #[test]
    fn test_kmer_key_still_available_for_other_commands() {
        // KmerKey remains used by dedup/hetmers; make sure string fallback
        // keys still resolve distinctly.
        use crate::common::KmerKey;
        let a = KmerKey::from_kmer("ACGTACGTACGTACGTACGTACGTACGTACGT");
        let b = KmerKey::from_kmer("ACGTACGTACGTACGTACGTACGTACGTACGA");
        assert_ne!(a, b);
    }
}
