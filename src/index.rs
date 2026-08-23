use crate::common;
use bio::bio_types::genome::AbstractLocus;
use bio::io::fasta::IndexedReader;
use rust_htslib::bcf::Read;
use rust_htslib::bcf::Reader;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};

/// Remove k-mers shared across alleles of the same locus.
fn find_dup_kmers(mut data: Vec<Vec<String>>) -> Vec<Vec<String>> {
    let mut counts: HashMap<String, usize> = HashMap::new();
    for inner_vec in &data {
        for s in inner_vec {
            *counts.entry(s.to_string()).or_insert(0) += 1;
        }
    }
    let dup_kmers: Vec<String> = counts
        .into_iter()
        .filter(|(_key, value)| *value > 1)
        .map(|(key, _value)| key)
        .collect();
    let dup_kmers_hashset: HashSet<String> = dup_kmers.into_iter().collect();
    for inner_vec in &mut data {
        inner_vec.retain(|s| !dup_kmers_hashset.contains(s));
    }
    data
}

/// Build an index of allele-specific k-mers from a VCF and indexed FASTA.
pub fn index_workflow(vcf_path: &str, fasta_path: &str, output_path: &str, k: i64) {
    let mut vcf_reader = Reader::from_path(vcf_path).expect("Error opening file.");
    let mut faidx = IndexedReader::from_file(&fasta_path.to_string()).unwrap();
    let chrom_lengths = common::read_fai(fasta_path);
    log::info!("Chromosome lengths:");
    log::info!("{:?}", chrom_lengths);
    let mut buffered_file =
        BufWriter::new(File::create(output_path).expect("Error creating output file"));
    let mut i = 0;
    let mut chrom_prev = String::new();
    let mut pos_prev = 0;
    let mut dup_tracker = HashSet::new();
    let mut vcf_iterator = vcf_reader.records().peekable();
    while let Some(record_result) = vcf_iterator.next() {
        i += 1;
        let record = record_result.expect("Failed to read record!");
        let pos = record.pos() - 1;
        let chrom = record.contig();
        let mut ku = k as usize;
        log::debug!(
            "Extracting allele sequences for CHROM: {} POS: {}...",
            chrom,
            pos
        );
        let mut alleles = String::new();
        for allele in record.alleles() {
            for c in allele {
                alleles.push(char::from(*c))
            }
            alleles.push(' ')
        }
        let alleles_list: Vec<&str> = alleles.split_whitespace().collect();
        log::debug!("Alleles: {:?}", alleles_list);
        if (alleles_list[1..alleles_list.len()]).contains(&alleles_list[0]) {
            log::warn!(
                "Invariant site detected, skipping CHROM: {} POS: {}",
                chrom,
                pos
            );
            continue;
        }
        let ref_allele_len = alleles_list[0].len();
        let mut region_start = pos - k + 1;
        let mut region_end = pos + k + (ref_allele_len as i64);
        log::debug!("Putative variant region: {}-{}", region_start, region_end);
        let chrom_pos_str = pos.to_string() + "-" + chrom;
        log::debug!("Checking for duplicate: {:?}", chrom_pos_str);
        if dup_tracker.contains(&chrom_pos_str) {
            log::warn!(
                "Skipping duplicate variant found at CHROM: {} POS: {}",
                chrom,
                pos
            );
            continue;
        }
        dup_tracker.insert(chrom_pos_str);
        if pos <= 1 {
            log::warn!(
                "Skipping variant at tip of chromosome: CHROM: {}, POS: {}",
                chrom,
                pos
            );
            continue;
        }
        if (pos_prev > pos) && (chrom == chrom_prev) {
            log::error!(
                "VCF is not sorted! CHROM: {} POS: {} occurs before CHROM: {} POS: {}",
                chrom_prev,
                pos_prev,
                chrom,
                pos
            );
            panic!();
        }
        if let Some(num_ref) = chrom_lengths
            .as_ref()
            .expect("Failed to read .fa.fai file!")
            .get(chrom)
        {
            let end = *num_ref;
            if end < k {
                log::warn!(
                    "Chromosome shorter than k, skipping CHROM: {} POS: {}",
                    chrom,
                    pos
                );
                continue;
            }
            if pos >= (end - k - (ref_allele_len as i64)) {
                log::warn!(
                    "Reference allele overlaps with chromosome end. Skipping current variant."
                );
                continue;
            } else if pos >= (end - k) {
                log::debug!(
                    "Variant near chromosome end detected, CHROM: {} POS: {}",
                    chrom,
                    pos
                );
                region_end = end;
                ku = k as usize;
            }
        } else {
            log::warn!(
                "Warning: {} not found in .fa.fai file! Thus, will skip CHROM: {} POS: {}",
                chrom,
                chrom,
                pos
            );
            continue;
        }
        if pos <= k {
            log::debug!(
                "Variant near chromosome start detected, CHROM: {} POS: {}",
                chrom,
                pos
            );
            region_start = 1;
            ku = pos as usize;
        }
        if let Some(next_ref) = vcf_iterator.peek() {
            let next_result = next_ref.as_ref().unwrap();
            let pos_next = next_result.pos() - 1;
            let chrom_next = next_result.contig();
            if (pos_next - pos_prev <= k) && (chrom_next == chrom_prev) {
                log::warn!("Current variant (CHROM: {} POS: {}) within k bp of previous and next variant. Skipping current variant.", chrom, pos);
                pos_prev = pos;
                chrom_prev = chrom.to_string();
                continue;
            }
            if ((pos_next - pos) <= k) && (chrom == chrom_next) {
                log::debug!("Current variant (CHROM: {} POS: {}) within k bp of next variant (CHROM: {} POS: {})", chrom, pos, chrom_next, pos_next);
                if (pos_next - pos) <= ref_allele_len.try_into().unwrap() {
                    log::warn!("Current variant overlaps next variant. Skipping current variant");
                    pos_prev = pos;
                    chrom_prev = chrom.to_string();
                    continue;
                } else {
                    region_end = pos_next;
                }
            }
        }
        if ((pos - pos_prev) <= k) && (chrom == chrom_prev) {
            log::debug!("Current variant (CHROM: {} POS: {}) within k bp previous variant (CHROM: {} POS: {}).", chrom, pos, chrom_prev, pos_prev);
            region_start = pos_prev + 1;
            ku = (pos - pos_prev) as usize;
        }
        if (region_end - region_start) < k {
            log::warn!("Reference region < k bp. Skipping current variant");
            pos_prev = pos;
            chrom_prev = chrom.to_string();
            continue;
        }
        log::debug!(
            "Fetching sequence {}:{}-{} from FASTA...",
            chrom,
            region_start,
            region_end
        );
        faidx
            .fetch(
                chrom,
                region_start.try_into().unwrap(),
                region_end.try_into().unwrap(),
            )
            .expect("Couldn't fetch interval");
        let mut seq = Vec::new();
        faidx.read(&mut seq).expect("Couldn't read the interval");
        let seq_string =
            common::stand_seq(std::str::from_utf8(&seq).expect("Invalid UTF-8 sequence"));
        log::debug!("Sequence length: {}", seq_string.chars().count());
        let mut var_seqs = Vec::new();
        log::debug!(
            "Replace range from start: {} end: {}",
            ku,
            ku + ref_allele_len
        );
        for allele in &alleles_list {
            let mut var_seq = seq_string.clone();
            var_seq.replace_range(ku..ku + ref_allele_len, &common::stand_seq(allele));
            var_seqs.push(var_seq);
        }
        if var_seqs[0] != seq_string {
            // If the only mismatches are positions where the FASTA has N (from IUPAC ambiguity
            // codes like Y, W, M, etc. that stand_seq() converts to N), skip gracefully.
            // k-mers containing N are discarded downstream anyway, so these variants contribute
            // nothing to the index regardless.
            let n_mismatch_only = seq_string
                .chars()
                .zip(var_seqs[0].chars())
                .all(|(f, v)| f == v || f == 'N');
            if n_mismatch_only {
                log::warn!("Skipping variant at CHROM: {} POS: {} — FASTA has N (IUPAC ambiguity code) where VCF REF has '{}'. k-mers with N are skipped downstream anyway.", chrom, pos, alleles_list[0]);
                pos_prev = pos;
                chrom_prev = chrom.to_string();
                continue;
            }
            log::error!("REF allele does not match FASTA at CHROM: {} POS: {}\nFASTA : {}\nVCF   : {}\nREGION: {} \nREF: {}\nWas the FASTA used as the reference to make the VCF?", chrom, pos, seq_string, var_seqs[0], region_start.to_string() + "-" + &region_end.to_string(), alleles_list[0]);
            panic!();
        }
        log::debug!("Extracting canonical k-mers from alleles...");
        let mut kmers_by_allele = Vec::new();
        for var_seq in &var_seqs {
            let kmer_list: Vec<String> = common::get_canonical_kmers(var_seq, k as usize);
            kmers_by_allele.push(kmer_list);
        }
        log::debug!("Dedupping k-mers shared across alleles of the same locus...");
        let kmers_by_allele_no_dups = find_dup_kmers(kmers_by_allele);
        let mut joined_kmers_list: Vec<String> = Vec::new();
        let mut num_kmers_per_allele: Vec<String> = Vec::new();
        for inner_vec in kmers_by_allele_no_dups {
            let inner_joined = inner_vec.join(";");
            joined_kmers_list.push(inner_joined);
            num_kmers_per_allele.push(inner_vec.len().to_string());
        }
        let parts = [
            i.to_string(),
            chrom.to_string(),
            pos.to_string(),
            seq_string.clone(),
            alleles_list.join("|"),
            var_seqs.join("|"),
            num_kmers_per_allele.join("|"),
            joined_kmers_list.join("|"),
        ];
        writeln!(buffered_file, "{}", parts.join(",")).expect("Error writing index");
        chrom_prev = chrom.to_string();
        pos_prev = pos;
    }
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    #[test]
    fn test_dedup_across_alleles() {
        let test_vec = vec![
            vec!["ATGC".to_string(), "AAAA".to_string(), "CCCC".to_string()],
            vec!["CATG".to_string(), "ATGC".to_string(), "GGCG".to_string()],
        ];
        let result = find_dup_kmers(test_vec);
        let expected = vec![
            vec!["AAAA".to_string(), "CCCC".to_string()],
            vec!["CATG".to_string(), "GGCG".to_string()],
        ];
        assert_eq!(result, expected);
    }
}
