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
pub fn find_dup_kmers_across_var(index: &str, output: &str) -> Result<(), io::Error> {
    log::info!("Reading index...");
    let mut data = read_index_kmers(index)?;

    log::info!("First pass: counting allele-specific k-mers...");
    let mut counts: HashMap<String, usize> = HashMap::new();

    for inner_vec in &data {
        for inner_inner_vec in inner_vec {
            for s in inner_inner_vec {
                *counts.entry(s.to_string()).or_insert(0) += 1;
            }
        }
    }

    let dup_kmers: Vec<String> = counts
        .into_iter()
        .filter(|(_key, value)| *value > 1)
        .map(|(key, _value)| key)
        .collect();

    let dup_kmers_hashset: HashSet<String> = dup_kmers.into_iter().collect();

    log::info!("Second pass: removing k-mers found more than once...");
    for inner_vec in &mut data {
        for inner_inner_vec in inner_vec {
            inner_inner_vec.retain(|s| !dup_kmers_hashset.contains(s));
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
            .map(|inner_vec| inner_vec.len().to_string())
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
