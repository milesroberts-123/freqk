use bio::bio_types::genome::AbstractLocus;
use bio::io::fasta::IndexedReader;
use rust_htslib::bcf::Read;
use rust_htslib::bcf::Reader;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use std::io::{BufWriter, Write};

use crate::common;

/// Read the k-mer field (column 7) of an index file into a 3D array:
/// one entry per index line, one inner vec per allele, one string per k-mer.
fn read_index_kmers(index: &str) -> Result<Vec<Vec<Vec<String>>>, io::Error> {
    let file = File::open(index)?;
    let reader = BufReader::new(file);
    let mut data = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        let kmers = fields[7];
        let kmers_list: Vec<&str> = kmers.split('|').collect();
        let kmers_by_allele: Vec<Vec<String>> = kmers_list
            .iter()
            .map(|s| s.split(';').map(|x| x.to_string()).collect())
            .collect();
        data.push(kmers_by_allele);
    }
    Ok(data)
}

/// Pack a canonical ATGC k-mer into an integer (2 bits per base) for cheap
/// hashing, and key for counting k-mers in a hash map: packed `u64` for
/// k <= 31 (the normal case); heap `String` fallback for longer k-mers so
/// behavior is unchanged for any k. Both live in `common`, shared with `count`.
#[cfg(test)]
use crate::common::pack_kmer;
use crate::common::KmerKey;

/// Parse the k-mer field (column 7) of one index line into a flat list of
/// k-mers, one entry per k-mer occurrence (including the empty-allele
/// pseudo-entry "", matching the original counting behavior).
fn kmers_from_line(line: &str) -> Vec<&str> {
    let fields: Vec<&str> = line.split(',').collect();
    fields[7].split('|').flat_map(|s| s.split(';')).collect()
}

/// Build a hashset of all k-mers in non-variable reference regions.
pub fn reference_hashset(index: &str, fasta_path: &str, vcf_path: &str) -> HashSet<String> {
    log::info!("Reading k-mer length from index...");
    let k = common::k_from_index(index).expect("Error reading k-mer length from index.");
    log::info!("k is: {:?}", k);
    log::info!("Build hashset of reference k-mers...");
    let mut vcf_reader = Reader::from_path(vcf_path).expect("Error opening file.");
    let mut faidx = IndexedReader::from_file(&fasta_path.to_string()).unwrap();
    let chrom_lengths = common::read_fai(fasta_path);
    log::info!("Chromosome lengths:");
    log::info!("{:?}", chrom_lengths);
    let mut start = 1;
    let mut ref_kmers_hashset = HashSet::new();
    let mut chrom_visited: HashSet<String> = HashSet::new();
    let mut vcf_iterator = vcf_reader.records().peekable();
    while let Some(record_result) = vcf_iterator.next() {
        let record = record_result.expect("Failure reading record");
        let pos = record.pos() - 1;
        let chrom = record.contig();
        chrom_visited.insert(chrom.into());
        let end = pos - k;
        log::debug!("Processing record CHROM: {} POS: {}", chrom, pos);
        log::debug!("Current region start: {}", start);
        log::debug!("Extracting allele sequences from VCF...");
        let mut alleles = String::new();
        for allele in record.alleles() {
            for c in allele {
                alleles.push(char::from(*c))
            }
            alleles.push(' ')
        }
        let alleles_list: Vec<&str> = alleles.split_whitespace().collect();
        log::debug!("Alleles: {:?}", alleles_list);
        let ref_allele_len = alleles_list[0].len() as i64;
        log::debug!("Reference allele length: {}", ref_allele_len);
        if pos <= 1 {
            log::warn!("Skipping current variant (CHROM: {} POS: {}) because its at the beginning of the chromosome", chrom, pos);
            start = 1 + k + (ref_allele_len - 1);
            continue;
        }
        if (pos - start) < k {
            log::warn!(
                "Current variant (CHROM: {} POS: {}) within k bp of region start ({}), so skipping",
                chrom,
                pos,
                start
            );
            start = pos + k + (ref_allele_len - 1);
            continue;
        }
        if let Some(next_ref) = vcf_iterator.peek() {
            let next_result = next_ref.as_ref().unwrap();
            let pos_next = next_result.pos() - 1;
            let chrom_next = next_result.contig();
            log::debug!("Next record is CHROM: {} POS: {}", chrom_next, pos_next);
            if chrom != chrom_next {
                log::info!(
                    "Chromosome boundry reached between CHROM: {} POS: {} and CHROM: {} POS: {}",
                    chrom,
                    pos,
                    chrom_next,
                    pos_next
                );
                if let Some(chrom_end) = chrom_lengths
                    .as_ref()
                    .expect("Error reading chromosome lengths")
                    .get(chrom)
                {
                    log::debug!(
                        "First, extracting sequence before current variant: {} {} - {}",
                        chrom,
                        start,
                        end
                    );
                    faidx
                        .fetch(chrom, start.try_into().unwrap(), end.try_into().unwrap())
                        .expect("Could not fetch interval");
                    log::debug!("Reading sequence...");
                    let mut seq = Vec::new();
                    faidx.read(&mut seq).expect("Could not read interval");
                    let seq_string =
                        String::from_utf8(seq.to_vec()).expect("Invalid UTF-8 sequence");
                    log::debug!("Extract canonical k-mers...");
                    let ref_kmers: Vec<String> =
                        common::get_canonical_kmers(&seq_string, k as usize);
                    log::debug!("Putting k-mers into hashset...");
                    for ref_kmer in ref_kmers {
                        ref_kmers_hashset.insert(ref_kmer);
                    }
                    log::debug!(
                        "Second, extracting sequence from POS: {} to end of {} at : {:?}",
                        pos,
                        chrom,
                        chrom_end
                    );
                    start = pos + k + (ref_allele_len - 1);
                    if start >= *chrom_end {
                        log::debug!("Start exceeds chrom end, skipping");
                        start = 1;
                        continue;
                    } else if start < *chrom_end {
                        log::debug!("start within k bp of chrom end, extracting");
                        faidx
                            .fetch(chrom, start.try_into().unwrap(), *chrom_end as u64)
                            .expect("Could not fetch interval");
                        start = 1;
                    }
                } else {
                    log::error!("Error getting length of chromosome");
                    panic!();
                }
            } else {
                log::debug!("Extracting sequence: {}:{}-{}", chrom, start, end);
                faidx
                    .fetch(chrom, start.try_into().unwrap(), end.try_into().unwrap())
                    .expect("Could not fetch interval");
                start = pos + k + (ref_allele_len - 1);
            }
        } else {
            log::debug!("No next record, so end of VCF reached. Grab remainder of chromosome");
            if let Some(chrom_end) = chrom_lengths
                .as_ref()
                .expect("Error reading chromosome lengths")
                .get(chrom)
            {
                log::debug!(
                    "First, extracting sequence before current variant: {} {} - {}",
                    chrom,
                    start,
                    end
                );
                faidx
                    .fetch(chrom, start.try_into().unwrap(), end.try_into().unwrap())
                    .expect("Could not fetch interval");
                log::debug!("Reading sequence...");
                let mut seq = Vec::new();
                faidx.read(&mut seq).expect("Could not read interval");
                let seq_string = String::from_utf8(seq.to_vec()).expect("Invalid UTF-8 sequence");
                log::debug!("Extract canonical k-mers...");
                let ref_kmers: Vec<String> = common::get_canonical_kmers(&seq_string, k as usize);
                log::debug!("Putting k-mers into hashset...");
                for ref_kmer in ref_kmers {
                    ref_kmers_hashset.insert(ref_kmer);
                }
                log::debug!(
                    "Second, attempting to extract sequence from POS: {} to end of {} at : {:?}",
                    pos,
                    chrom,
                    chrom_end
                );
                start = pos + k + (ref_allele_len - 1);
                if start >= *chrom_end {
                    log::debug!("POS within k bp of chrom end, breaking loop");
                    break;
                } else if start < *chrom_end {
                    log::debug!("POS not within k bp of chrom end, extracting");
                    faidx
                        .fetch(chrom, (pos + k).try_into().unwrap(), *chrom_end as u64)
                        .expect("Could not fetch interval");
                    start = 1;
                }
            } else {
                log::error!("Error getting length of chromosome.");
            }
        }
        log::debug!("Reading sequence...");
        let mut seq = Vec::new();
        faidx.read(&mut seq).expect("Could not read interval");
        let seq_string = String::from_utf8(seq.to_vec()).expect("Invalid UTF-8 sequence");
        log::debug!("Extract canonical k-mers...");
        let ref_kmers: Vec<String> = common::get_canonical_kmers(&seq_string, k as usize);
        log::debug!("Putting k-mers into hashset...");
        for ref_kmer in ref_kmers {
            ref_kmers_hashset.insert(ref_kmer);
        }
    }
    let binding = chrom_lengths
        .as_ref()
        .expect("Error unpacking chromosome lengths");
    let unvisted_chroms: Vec<String> = binding
        .keys()
        .filter(|x| !chrom_visited.contains(*x))
        .cloned()
        .collect();
    log::warn!("Unvisited chromosomes: {:?}", unvisted_chroms);
    for unvis_chrom in unvisted_chroms {
        if let Some(chrom_end) = chrom_lengths
            .as_ref()
            .expect("Error reading chromosome lengths")
            .get(&unvis_chrom)
        {
            log::debug!(
                "Extracting sequence on {:?} from 1 to {}",
                unvis_chrom,
                chrom_end
            );
            faidx
                .fetch(&unvis_chrom, 1, *chrom_end as u64)
                .expect("Could not fetch interval");
            log::debug!("Reading sequence...");
            let mut seq = Vec::new();
            faidx.read(&mut seq).expect("Could not read interval");
            let seq_string = String::from_utf8(seq.to_vec()).expect("Invalid UTF-8 sequence");
            log::debug!("Extract canonical k-mers...");
            let ref_kmers: Vec<String> = common::get_canonical_kmers(&seq_string, k as usize);
            log::debug!("Putting k-mers into hashset...");
            for ref_kmer in ref_kmers {
                ref_kmers_hashset.insert(ref_kmer);
            }
        } else {
            log::error!("Error getting length of chromosome.");
        }
    }
    ref_kmers_hashset
}

/// Remove k-mers found in the reference from an index and write a new index.
pub fn remove_ref_kmers(
    index: &str,
    output: &str,
    ref_hashset: HashSet<String>,
) -> Result<(), io::Error> {
    log::info!("Opening index...");
    let mut data = read_index_kmers(index)?;
    log::info!("Removing k-mers found in non-variable sequences...");
    for inner_vec in &mut data {
        for inner_inner_vec in inner_vec {
            inner_inner_vec.retain(|s| !ref_hashset.contains(s));
        }
    }
    log::info!("Writing new index...");
    let mut buffered_file = BufWriter::new(File::create(output)?);
    let file = File::open(index)?;
    let reader = BufReader::new(file);
    for (i, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        let dedup_kmers = &data[i];

        let num_kmers_per_allele = dedup_kmers
            .iter()
            .map(|inner_vec| {
                if *inner_vec == vec![""] {
                    "0".to_string()
                } else {
                    inner_vec.len().to_string()
                }
            })
            .collect::<Vec<String>>()
            .join("|");
        let dedup_string = dedup_kmers
            .iter()
            .map(|inner_vec| inner_vec.join(";"))
            .collect::<Vec<String>>()
            .join("|");
        let parts = [
            fields[0],
            fields[1],
            fields[2],
            fields[3],
            fields[4],
            fields[5],
            &num_kmers_per_allele,
            &dedup_string,
        ];
        writeln!(buffered_file, "{}", parts.join(","))?;
    }
    Ok(())
}

/// Remove k-mers shared across variants from an index and write a new index.
///
/// Streams the index file twice instead of loading it into memory: the first
/// pass counts k-mer occurrences in a `HashMap<KmerKey, u32>`, the second pass
/// re-reads each line, drops duplicated k-mers, and rewrites the line.
/// Memory is proportional to the number of unique k-mers, not index size.
pub fn find_dup_kmers_across_var(index: &str, output: &str) -> Result<(), io::Error> {
    log::info!("First pass: counting allele-specific k-mers...");
    let mut counts: HashMap<KmerKey, u32> = HashMap::new();
    {
        let file = File::open(index)?;
        let reader = BufReader::new(file);
        for line_result in reader.lines() {
            let line = line_result?;
            for kmer in kmers_from_line(&line) {
                *counts.entry(KmerKey::from_kmer(kmer)).or_insert(0) += 1;
            }
        }
    }

    log::info!("Removing k-mers found only once from the count table...");
    counts.retain(|_, count| *count > 1);

    log::info!("Second pass: removing k-mers found more than once and writing new index...");
    let mut buffered_file = BufWriter::new(File::create(output)?);
    let file = File::open(index)?;
    let reader = BufReader::new(file);

    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();

        let kmers_by_allele: Vec<Vec<&str>> = fields[7]
            .split('|')
            .map(|s| s.split(';').collect())
            .collect();
        let dedup_kmers: Vec<Vec<&str>> = kmers_by_allele
            .iter()
            .map(|inner_vec| {
                inner_vec
                    .iter()
                    .filter(|s| !counts.contains_key(&KmerKey::from_kmer(s)))
                    .cloned()
                    .collect()
            })
            .collect();

        let num_kmers_per_allele = dedup_kmers
            .iter()
            .map(|inner_vec| {
                if *inner_vec == vec![""] {
                    "0".to_string()
                } else {
                    inner_vec.len().to_string()
                }
            })
            .collect::<Vec<String>>()
            .join("|");

        let dedup_string = dedup_kmers
            .iter()
            .map(|inner_vec| inner_vec.join(";"))
            .collect::<Vec<String>>()
            .join("|");

        let parts = [
            fields[0],
            fields[1],
            fields[2],
            fields[3],
            fields[4],
            fields[5],
            &num_kmers_per_allele,
            &dedup_string,
        ];
        writeln!(buffered_file, "{}", parts.join(","))?;
    }

    log::info!("Writing successful! :D");
    Ok(())
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    #[test]
    fn test_pack_kmer() {
        // "AAAC" packs to 0b000001 = 1
        assert_eq!(pack_kmer("AAAC"), Some(1));
        // Full 31-mer is the packing limit: 0b0011011010... (2 bits per base)
        let packed_31 = "ACGTACGTACGTACGTACGTACGTACGTACC"
            .bytes()
            .fold(0u64, |acc, c| {
                (acc << 2)
                    | match c {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        _ => 3,
                    }
            });
        assert_eq!(
            pack_kmer("ACGTACGTACGTACGTACGTACGTACGTACC"),
            Some(packed_31)
        );
        // Empty allele pseudo-entry and non-ATGC fall back to None
        assert_eq!(pack_kmer(""), None);
        assert_eq!(pack_kmer("ATGCN"), None);
        // Longer than 31 falls back to None
        assert_eq!(pack_kmer("ACGTACGTACGTACGTACGTACGTACGTACGT"), None);
    }

    #[test]
    fn test_kmer_key_roundtrip() {
        assert_eq!(KmerKey::from_kmer("AAAC"), KmerKey::Packed(1));
        assert_eq!(KmerKey::from_kmer(""), KmerKey::Str(String::new()));
        assert_eq!(
            KmerKey::from_kmer("ACGTACGTACGTACGTACGTACGTACGTACGT"),
            KmerKey::Str("ACGTACGTACGTACGTACGTACGTACGTACGT".to_string())
        );
    }

    #[test]
    fn test_kmers_from_line() {
        let line = "0,1,100,SEQ,REF|ALT,SEQ|ALT,2|1,AAA;CCC|GGG";
        assert_eq!(kmers_from_line(line), vec!["AAA", "CCC", "GGG"]);
        // Empty allele pseudo-entry is preserved for counting
        let line_empty = "0,1,100,SEQ,REF|ALT,SEQ|ALT,0|0,|";
        assert_eq!(kmers_from_line(line_empty), vec!["", ""]);
    }

    #[test]
    fn test_find_dup_kmers_across_var_end_to_end() {
        let dir = tempfile::tempdir().expect("tempdir");
        let index_path = dir.path().join("index.txt");
        let output_path = dir.path().join("dedup.txt");
        // Variant 1: AAA appears twice within the line (duplicate) -> removed
        // everywhere; CCC and TTT appear once and stay.
        // Variant 2 has an empty allele (pseudo-entry "" in both fields).
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,AAAA,REF|ALT,AAAA|ACAA,2|2,AAA;CCC|AAA;TTT\n",
                "1,1,200,CCCC,REF|ALT,CCCC|CCGG,0|1,|GGG\n",
            ),
        )
        .expect("write index");
        find_dup_kmers_across_var(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        let expected = concat!(
            "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,CCC|TTT\n",
            // GGG appears once -> stays; empty allele still yields 0.
            "1,1,200,CCCC,REF|ALT,CCCC|CCGG,0|1,|GGG\n",
        );
        assert_eq!(result, expected);
    }
}
