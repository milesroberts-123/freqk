use fastq::{parse_path, Record};
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};

use crate::common;

/// Load all allele-specific k-mers from an index file into a single hashset.
pub fn build_kmer_hashset(index: &str) -> Result<HashSet<String>, io::Error> {
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    let mut kmers_hashset: HashSet<String> = HashSet::new();

    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        let kmers = fields[7];
        let kmers_list: Vec<&str> = kmers.split('|').collect();

        let kmers_by_allele: Vec<Vec<&str>> =
            kmers_list.iter().map(|s| s.split(';').collect()).collect();

        let kmers_all: Vec<&str> = kmers_by_allele.into_iter().flatten().collect();

        kmers_hashset.extend(kmers_all.iter().map(|s| s.to_string()));
    }

    Ok(kmers_hashset)
}

fn merge_hashmaps(vec_of_maps: Vec<HashMap<String, usize>>) -> HashMap<String, usize> {
    let mut merged_map: HashMap<String, usize> = HashMap::new();

    for map in vec_of_maps {
        for (key, value) in map {
            *merged_map.entry(key).or_insert(0) += value;
        }
    }
    merged_map
}

/// Count indexed k-mers in reads, in parallel across `nthreads` threads.
pub fn count_target_kmers_in_reads(
    index: &str,
    reads: &str,
    k: i64,
    nthreads: usize,
) -> HashMap<String, usize> {
    let kmers_hashset = build_kmer_hashset(index).expect("Error loading index");
    let k = k as usize;
    let merged_counts: HashMap<String, usize> = parse_path(Some(reads), |parser| {
        let results: Vec<HashMap<String, usize>> = parser
            .parallel_each(nthreads, move |record_sets| {
                let mut kmer_counts: HashMap<String, usize> = HashMap::new();
                let kmers_hashset = kmers_hashset.clone();
                let mut num_records = 0;
                for record_set in record_sets {
                    for record in record_set.iter() {
                        if num_records % 10000 == 0 {
                            log::info!("Reads processed: {}", num_records);
                        }
                        num_records += 1;
                        let read_kmers = common::get_canonical_kmers(
                            std::str::from_utf8(record.seq())
                                .expect("Invalid UTF-8 in read sequence"),
                            k,
                        );
                        for read_kmer in &read_kmers {
                            if kmers_hashset.contains(read_kmer) {
                                let count = kmer_counts.entry(read_kmer.to_string()).or_insert(0);
                                *count += 1;
                            }
                        }
                    }
                }
                kmer_counts
            })
            .expect("Invalid fastq file");
        log::debug!("Merging hashmaps...");
        merge_hashmaps(results)
    })
    .expect("Invalid compression");
    merged_counts
}

/// Write k-mer counts to a file.
pub fn write_kmers(kmer_counts: HashMap<String, usize>, output: &str) -> io::Result<()> {
    let mut file = File::create(output)?;

    for (key, value) in kmer_counts.iter() {
        writeln!(file, "{}\t{}", key, value)?;
    }

    Ok(())
}

/// Sum k-mer counts into totals per allele, one line per index entry.
pub fn combine_counts_by_allele(
    index: &str,
    counts: HashMap<String, usize>,
) -> Result<Vec<String>, io::Error> {
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    let mut result = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        let kmers = fields[7];
        let kmers_list: Vec<&str> = kmers.split('|').collect();

        let kmers_by_allele: Vec<Vec<&str>> =
            kmers_list.iter().map(|s| s.split(';').collect()).collect();

        let mut counts_by_allele = Vec::new();

        for allele in kmers_by_allele {
            let mut total_allele_count: usize = 0;
            for kmer in allele {
                if let Some(kmer_count) = counts.get(kmer) {
                    total_allele_count += *kmer_count;
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
pub fn count_workflow(
    index: &str,
    reads: &str,
    nthreads: usize,
    freq_output: &str,
    count_output: &str,
) {
    log::info!("Loading index at {} into a hashset...", index);
    let k = common::k_from_index(index);
    log::info!("k is: {:?}", k);
    log::debug!("Counting indexed k-mers in reads...");
    let kmer_counts = count_target_kmers_in_reads(
        index,
        reads,
        k.expect("Cannot parse kmer length from index."),
        nthreads,
    );
    let _ = write_kmers(kmer_counts.clone(), count_output);
    log::debug!("Combining k-mer counts by allele...");
    let counts_by_allele = combine_counts_by_allele(index, kmer_counts);
    log::debug!("Writing counts by allele...");
    let _ = common::write_strings(
        counts_by_allele.expect("Error writing counts by allele"),
        freq_output,
    );
    log::debug!("Successfully wrote counts by allele!");
}
