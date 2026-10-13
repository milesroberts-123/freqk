use std::collections::HashMap;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

use crate::common;

/// Pack a canonical ATGC k-mer into an integer (2 bits per base) for cheap
/// hashing, and key for counting k-mers in a hash map: packed `u64` for
/// k <= 31 (the normal case); heap `String` fallback for longer k-mers so
/// behavior is unchanged for any k. Both live in `common`, shared with `count`.
use crate::common::pack_kmer;
use crate::common::KmerKey;
use crate::common::ShardedCounts;

/// Number of reference bases assigned per work item in the reference scan.
/// Each item fetches its segment plus the following k-1 bases, so every
/// k-mer owned by the segment (by start position) can be completed.
const REF_SEGMENT_BASES: u64 = 8_000_000;

/// Parse the k-mer field (column 7) of one index line into a flat list of
/// k-mers, one entry per k-mer occurrence (including the empty-allele
/// pseudo-entry "", matching the original counting behavior).
fn kmers_from_line(line: &str) -> Vec<&str> {
    let fields: Vec<&str> = line.split(',').collect();
    fields[7].split('|').flat_map(|s| s.split(';')).collect()
}

/// Exact reference occurrence counts, restricted to the k-mers present in
/// the index. The counting pass only ever needs counts of index k-mers, so
/// memory is proportional to the index, not to the reference genome.
///
/// `Packed` stores `u64`-packed keys (2 bits per base) in the sharded maps
/// shared with `count`; `Strings` is the `KmerKey`/`String` fallback for
/// k > 31. Both count canonical reference k-mers of length exactly `k`,
/// matching the keys the previous genome-wide count map held: an index k-mer
/// of another length, a non-ATGC k-mer, or the empty-allele pseudo-entry
/// yield 0, same as looking them up in that map.
pub enum RefKmerCounts {
    Packed {
        counts: Arc<ShardedCounts>,
        k: usize,
    },
    Strings {
        counts: Mutex<HashMap<KmerKey, u32>>,
        k: usize,
    },
}

impl RefKmerCounts {
    /// Count the canonical length-k k-mers of one fetched reference segment.
    fn count_seq_kmers(&self, seq: &str) {
        match self {
            RefKmerCounts::Packed { counts, k } => {
                for packed in common::get_canonical_kmers_packed(seq, *k) {
                    common::increment(counts, packed);
                }
            }
            RefKmerCounts::Strings { counts, k } => {
                let mut counts = counts.lock().unwrap();
                for kmer in common::get_canonical_kmers(seq, *k) {
                    if let Some(count) = counts.get_mut(&KmerKey::from_kmer(&kmer)) {
                        *count += 1;
                    }
                }
            }
        }
    }

    /// Reference occurrences of one index k-mer; 0 when the k-mer is not a
    /// counted member (different length than the index k, non-ATGC, empty
    /// pseudo-entry) — identical to its lookup in the former full map.
    fn count_of_kmer(&self, kmer: &str) -> usize {
        match self {
            RefKmerCounts::Packed { counts, k } => {
                if kmer.len() != *k {
                    return 0;
                }
                match pack_kmer(kmer) {
                    Some(packed) => {
                        let shard = counts[common::shard_of(packed)].lock().unwrap();
                        *shard.get(&packed).unwrap_or(&0)
                    }
                    None => 0,
                }
            }
            RefKmerCounts::Strings { counts, .. } => counts
                .lock()
                .unwrap()
                .get(&KmerKey::from_kmer(kmer))
                .copied()
                .unwrap_or(0) as usize,
        }
    }

    /// (distinct counted k-mers, of which occur more than once).
    fn stats(&self) -> (usize, usize) {
        match self {
            RefKmerCounts::Packed { counts, .. } => {
                let mut distinct = 0;
                let mut multi = 0;
                for shard in counts.iter() {
                    let shard = shard.lock().unwrap();
                    distinct += shard.len();
                    multi += shard.values().filter(|&&c| c > 1).count();
                }
                (distinct, multi)
            }
            RefKmerCounts::Strings { counts, .. } => {
                let counts = counts.lock().unwrap();
                (counts.len(), counts.values().filter(|&&c| c > 1).count())
            }
        }
    }
}

/// Register every length-k index k-mer as a counted member at count 0, so
/// the reference scan can bump exactly those keys (`increment` ignores
/// unknown keys). Only true length-k k-mers take part: k-mers of other
/// lengths, non-ATGC k-mers, and the empty pseudo-entry are not members,
/// matching the keys of the former genome-wide map.
fn prefill_ref_counts(index: &str, k: usize, packed: bool) -> Result<RefKmerCounts, io::Error> {
    log::info!("Pre-filling count map from index...");
    let counts = if packed {
        RefKmerCounts::Packed {
            counts: common::new_sharded_counts(),
            k,
        }
    } else {
        RefKmerCounts::Strings {
            counts: Mutex::new(HashMap::new()),
            k,
        }
    };
    let file = File::open(index)?;
    let reader = BufReader::new(file);
    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 8 {
            continue;
        }
        for kmer in fields[7].split('|').flat_map(|s| s.split(';')) {
            if kmer.len() != k {
                continue;
            }
            match &counts {
                RefKmerCounts::Packed { counts, .. } => {
                    if let Some(packed_kmer) = pack_kmer(kmer) {
                        let mut shard = counts[common::shard_of(packed_kmer)].lock().unwrap();
                        shard.insert(packed_kmer, 0);
                    }
                }
                RefKmerCounts::Strings { counts, .. } => {
                    counts
                        .lock()
                        .unwrap()
                        .entry(KmerKey::from_kmer(kmer))
                        .or_insert(0);
                }
            }
        }
    }
    Ok(counts)
}

/// Count reference occurrences of every index k-mer: pre-fill the count map
/// from the index, then scan the reference in segments (each segment fetches
/// the following k-1 bases so its k-mers, owned by start position, are all
/// completed exactly once; the next segment starts fresh at the boundary).
/// Scanning is split across `nthreads` workers, each opening its own faidx
/// reader; empty contigs yield the same fetch error as the former
/// whole-chromosome fetch.
fn ref_counts(
    index: &str,
    fasta_path: &str,
    k: i64,
    nthreads: usize,
    segment_bases: u64,
) -> Result<RefKmerCounts, io::Error> {
    let ki = k.max(0) as usize;
    let counts = prefill_ref_counts(index, ki, ki <= 31)?;
    log::info!("Counting reference occurrences of index k-mers...");
    let chrom_lengths = common::read_fai(fasta_path)
        .map_err(|e| io::Error::other(format!("Error reading chromosome lengths: {}", e)))?;

    // Work item: (chromosome, segment start, fetch stop). The fetch extends
    // into the next segment's first k-1 bases; k-mer ownership by start
    // position keeps counts exact.
    let mut items: Vec<(String, u64, u64)> = Vec::new();
    for (chrom, length) in chrom_lengths.iter() {
        if *length <= 0 {
            items.push((chrom.clone(), 0, 0));
            continue;
        }
        let length = *length as u64;
        let mut seg_start: u64 = 0;
        while seg_start < length {
            let seg_end = (seg_start + segment_bases).min(length);
            let stop = (seg_end + ki.saturating_sub(1) as u64).min(length);
            items.push((chrom.clone(), seg_start, stop));
            seg_start = seg_end;
        }
    }

    let cursor = AtomicUsize::new(0);
    let error: Mutex<Option<io::Error>> = Mutex::new(None);
    let nworkers = nthreads.max(1);
    thread::scope(|scope| {
        for _ in 0..nworkers {
            scope.spawn(|| {
                let faidx = match rust_htslib::faidx::Reader::from_path(fasta_path) {
                    Ok(faidx) => faidx,
                    Err(e) => {
                        let mut slot = error.lock().unwrap();
                        if slot.is_none() {
                            *slot = Some(io::Error::other(format!(
                                "Error opening FASTA index: {}",
                                e
                            )));
                        }
                        return;
                    }
                };
                loop {
                    let idx = cursor.fetch_add(1, Ordering::SeqCst);
                    if idx >= items.len() {
                        break;
                    }
                    let (chrom, start, stop) = &items[idx];
                    log::debug!("Scanning {}:{}-{}", chrom, start, stop);
                    match common::fetch_fasta(&faidx, chrom, *start, *stop) {
                        Ok(seq) => counts.count_seq_kmers(&seq),
                        Err(e) => {
                            log::error!("Fetching {} failed: {}", chrom, e);
                            let mut slot = error.lock().unwrap();
                            if slot.is_none() {
                                *slot = Some(io::Error::other(format!(
                                    "Fetching {} failed: {}",
                                    chrom, e
                                )));
                            }
                        }
                    }
                }
            });
        }
    });
    if let Some(e) = error.into_inner().unwrap() {
        return Err(e);
    }

    let (distinct, multi) = counts.stats();
    log::info!(
        "Counted occurrences for {} index k-mers, of which {} occur more than once in the reference",
        distinct,
        multi
    );
    Ok(counts)
}

/// Run ref-dedup end to end: count reference k-mer occurrences, then rewrite
/// the index dropping k-mers that occur outside their own variant's locus
/// (see [`remove_ref_kmers`]).
pub fn ref_dedup(
    index: &str,
    fasta_path: &str,
    output: &str,
    min_alleles: usize,
    nthreads: usize,
) -> Result<(), io::Error> {
    let k = common::k_from_index(index)?;
    log::info!("k is: {:?}", k);
    let counts = ref_counts(index, fasta_path, k, nthreads, REF_SEGMENT_BASES)?;
    remove_ref_kmers(index, output, &counts, min_alleles)
}

/// Maximum reference occurrences a kept k-mer may have, given its allele
/// position within its index row: the REF allele (group 0) is allowed one
/// occurrence, its own variant's locus; an ALT-allele k-mer can never match
/// at its own variant (its allele base differs from the reference base
/// there), so every reference occurrence is a stray and must be dropped.
fn max_kept_occurrences(group_index: usize) -> usize {
    if group_index == 0 {
        1
    } else {
        0
    }
}

/// Count the alleles of a deduplicated variant that have at least one
/// non-empty k-mer.
fn count_alleles_with_kmers<S: AsRef<str>>(dedup_kmers: &[Vec<S>]) -> usize {
    dedup_kmers
        .iter()
        .filter(|inner_vec| inner_vec.iter().any(|s| !s.as_ref().is_empty()))
        .count()
}

/// Remove k-mers found in the reference from an index and write a new index.
///
/// DROP RULE: every index k-mer spans its variant, so a REF-allele k-mer
/// always occurs exactly once in the reference — at its own variant — while
/// an ALT-allele k-mer never does. A k-mer therefore has an occurrence
/// outside its own variant's locus iff its reference count exceeds its
/// allele's threshold (REF: 1, ALT: 0); such k-mers are dropped.
///
/// Streams the index in a single pass instead of loading it into memory: each
/// line is parsed, its k-mers filtered against the reference counts
/// ([`RefKmerCounts`], restricted to index k-mers), and the line rewritten
/// immediately. Memory is proportional to the number of distinct index
/// k-mers, not to the reference genome. Rows are dropped entirely when fewer
/// than `min_alleles` alleles retain at least one k-mer (0 keeps every row).
pub fn remove_ref_kmers(
    index: &str,
    output: &str,
    ref_counts: &RefKmerCounts,
    min_alleles: usize,
) -> Result<(), io::Error> {
    log::info!("Writing new index, removing k-mers with excess reference occurrences...");
    let mut buffered_file = BufWriter::new(File::create(output)?);
    let file = File::open(index)?;
    let reader = BufReader::new(file);
    let mut dropped = 0;
    let mut written = 0;
    for (line_number, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 8 {
            log::warn!(
                "Skipping malformed index line {} (expected 8 comma-separated fields): {}",
                line_number + 1,
                line
            );
            continue;
        }

        let kmers_by_allele: Vec<Vec<&str>> = fields[7]
            .split('|')
            .map(|s| s.split(';').collect())
            .collect();
        let dedup_kmers: Vec<Vec<&str>> = kmers_by_allele
            .iter()
            .enumerate()
            .map(|(group_index, inner_vec)| {
                let max_kept = max_kept_occurrences(group_index);
                inner_vec
                    .iter()
                    .filter(|s| match ref_counts.count_of_kmer(s) {
                        0 => true,
                        count => count <= max_kept,
                    })
                    .cloned()
                    .collect()
            })
            .collect();

        if min_alleles > 0 && count_alleles_with_kmers(&dedup_kmers) < min_alleles {
            dropped += 1;
            continue;
        }
        written += 1;

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
    if min_alleles > 0 && dropped > 0 {
        log::info!(
            "Dropped {} index rows with fewer than {} alleles retaining k-mers",
            dropped,
            min_alleles
        );
    }
    log::info!("Wrote {} index rows", written);
    Ok(())
}

/// Remove k-mers shared across variants from an index and write a new index.
///
/// Streams the index file twice instead of loading it into memory: the first
/// pass counts k-mer occurrences in a `HashMap<KmerKey, u32>`, the second pass
/// re-reads each line, drops duplicated k-mers, and rewrites the line.
/// Memory is proportional to the number of unique k-mers, not index size.
/// Rows are dropped entirely when fewer than `min_alleles` alleles retain at
/// least one k-mer (0 keeps every row).
pub fn find_dup_kmers_across_var(
    index: &str,
    output: &str,
    min_alleles: usize,
) -> Result<(), io::Error> {
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
    let mut dropped = 0;

    for (line_number, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 8 {
            log::warn!(
                "Skipping malformed index line {} (expected 8 comma-separated fields): {}",
                line_number + 1,
                line
            );
            continue;
        }

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

        if min_alleles > 0 && count_alleles_with_kmers(&dedup_kmers) < min_alleles {
            dropped += 1;
            continue;
        }

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

    if min_alleles > 0 && dropped > 0 {
        log::info!(
            "Dropped {} index rows with fewer than {} alleles retaining k-mers",
            dropped,
            min_alleles
        );
    }
    log::info!("Writing successful! :D");
    Ok(())
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    /// Build `RefKmerCounts::Packed` with preset counts for `k`-mers (all
    /// entries must be length-k ATGC k-mers).
    fn packed_ref_counts(k: usize, entries: &[(&str, usize)]) -> RefKmerCounts {
        let sharded = common::new_sharded_counts();
        for (kmer, count) in entries {
            let packed = pack_kmer(kmer).expect("test k-mer must pack");
            let mut shard = sharded[common::shard_of(packed)].lock().unwrap();
            shard.insert(packed, *count);
        }
        RefKmerCounts::Packed { counts: sharded, k }
    }

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
    fn test_count_alleles_with_kmers() {
        assert_eq!(count_alleles_with_kmers(&[vec!["AAA"], vec!["GGG"]]), 2);
        assert_eq!(count_alleles_with_kmers(&[vec![], vec!["GGG"]]), 1);
        assert_eq!(count_alleles_with_kmers(&[vec![""], vec![]]), 0);
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
            0,
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

    #[test]
    fn test_find_dup_kmers_across_var_min_alleles_drops_rows() {
        let dir = tempfile::tempdir().expect("tempdir");
        let index_path = dir.path().join("index.txt");
        let output_path = dir.path().join("dedup.txt");
        // Variant 1 keeps 1 k-mer in each allele after dedup (CCC, TTT unique).
        // Variant 2 loses its only k-mer (AAA appears three times -> removed),
        // leaving both alleles empty.
        // Variant 3 has a pre-existing empty allele and GGG (unique) in allele 1.
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,CCC|TTT\n",
                "1,1,200,CCCC,REF|ALT,CCCC|CCGG,1|1,AAA|AAA\n",
                "2,1,300,TTTT,REF|ALT,TTTT|TTGG,1|0,GGG|\n",
            ),
        )
        .expect("write index");
        // min_alleles = 1: variant 2 is dropped (0 alleles with k-mers),
        // variants 1 and 3 stay.
        find_dup_kmers_across_var(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            1,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        let expected = concat!(
            "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,CCC|TTT\n",
            // GGG is unique -> retained; the empty allele yields 0.
            "2,1,300,TTTT,REF|ALT,TTTT|TTGG,1|0,GGG|\n",
        );
        assert_eq!(result, expected);
    }

    #[test]
    fn test_find_dup_kmers_across_var_min_alleles_two() {
        let dir = tempfile::tempdir().expect("tempdir");
        let index_path = dir.path().join("index.txt");
        let output_path = dir.path().join("dedup.txt");
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,CCC|TTT\n",
                "1,1,200,CCCC,REF|ALT,CCCC|CCGG,0|1,|GGG\n",
            ),
        )
        .expect("write index");
        // min_alleles = 2: variant 2 has only one allele with a k-mer -> dropped.
        find_dup_kmers_across_var(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            2,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        let expected = "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,CCC|TTT\n";
        assert_eq!(result, expected);
    }

    #[test]
    fn test_remove_ref_kmers_min_alleles() {
        let dir = tempfile::tempdir().expect("tempdir");
        let index_path = dir.path().join("index.txt");
        let output_path = dir.path().join("dedup.txt");
        // Variant row 1: REF allele k-mer AAA occurs twice in the reference —
        // one occurrence is its own locus, so a second one is a stray ->
        // dropped; TTT never occurs -> kept. Row 2: GGG occurs once (its own
        // locus) -> kept for the REF allele; CCC occurs twice -> dropped for
        // the ALT allele.
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,AAAA,REF|ALT,AAAA|ACAA,1|1,AAA|TTT\n",
                "1,1,200,CCCC,REF|ALT,CCCC|CCGG,1|1,GGG|CCC\n",
            ),
        )
        .expect("write index");
        let ref_counts = packed_ref_counts(3, &[("AAA", 2), ("GGG", 1), ("CCC", 2)]);

        // min_alleles = 0 keeps all rows.
        remove_ref_kmers(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            &ref_counts,
            0,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        let expected = concat!(
            "0,1,100,AAAA,REF|ALT,AAAA|ACAA,0|1,|TTT\n",
            "1,1,200,CCCC,REF|ALT,CCCC|CCGG,1|0,GGG|\n",
        );
        assert_eq!(result, expected);

        // min_alleles = 1 keeps both rows (each retains one allele).
        remove_ref_kmers(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            &ref_counts,
            1,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        assert_eq!(result, expected);

        // min_alleles = 2 drops both rows (one allele each).
        remove_ref_kmers(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            &ref_counts,
            2,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        assert_eq!(result, "");
    }

    #[test]
    fn test_remove_ref_kmers_two_tier_thresholds() {
        let dir = tempfile::tempdir().expect("tempdir");
        let index_path = dir.path().join("index.txt");
        let output_path = dir.path().join("dedup.txt");
        // The same single-occurrence k-mer is kept in a REF allele but dropped
        // in an ALT allele; any counted k-mer is dropped from both.
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,AAAA,REF|ALT2,AAAA|ACAA|AGAA,1|1|1,AAA|TTT|AAA\n",
                "1,1,200,CCCC,REF|ALT,CCCC|CCGG,1|1,AAA|AAA\n",
            ),
        )
        .expect("write index");
        let ref_counts = packed_ref_counts(3, &[("AAA", 1)]);
        remove_ref_kmers(
            index_path.to_str().expect("path"),
            output_path.to_str().expect("path"),
            &ref_counts,
            0,
        )
        .expect("dedup");
        let result = std::fs::read_to_string(&output_path).expect("read output");
        // Row 0: REF group AAA (count 1, own locus) and ALT group TTT
        // (count 0, no occurrence) are both kept; ALT group AAA (count 1,
        // a stray) is dropped.
        let expected = concat!(
            "0,1,100,AAAA,REF|ALT2,AAAA|ACAA|AGAA,1|1|0,AAA|TTT|\n",
            "1,1,200,CCCC,REF|ALT,CCCC|CCGG,1|0,AAA|\n",
        );
        assert_eq!(result, expected);
    }

    #[test]
    fn test_ref_counts() {
        let dir = tempfile::tempdir().expect("tempdir");
        let fasta_path = dir.path().join("ref.fasta");
        // k = 5. On chromosome a, ACGTA occurs twice (own canonical form);
        // CGTAA (canonical of CGTAA's reverse complement, TTACG being larger)
        // occurs once. Chromosome b: GGAAA. Chromosome d:
        // AAAAA (collides numerically with the 3-mer AAA when packed).
        // Chromosome c is shorter than k and contributes nothing.
        std::fs::write(
            &fasta_path,
            concat!(
                ">a\n",
                "ACGTAACGTAGGGG\n",
                ">b\n",
                "GGAAA\n",
                ">c\n",
                "AC\n",
                ">d\n",
                "AAAAAG\n",
            ),
        )
        .expect("write fasta");
        // .fai: name, length, offset, linebases, linewidth. Offsets account
        // for each ">name\n" header line (3 bytes) before its sequence.
        std::fs::write(
            dir.path().join("ref.fasta.fai"),
            "a\t14\t3\t14\t15\nb\t5\t21\t5\t6\nc\t2\t30\t2\t3\nd\t6\t35\t6\t7\n",
        )
        .expect("write fai");
        // Only length-k index k-mers become counted members. The GTAAC k-mer
        // is in the reference but not in the index, so it must not be
        // tracked; query k-mers of other lengths return 0 even when their
        // packed value collides with an indexed k-mer (AAA vs AAAAA).
        let index_path = dir.path().join("index.txt");
        std::fs::write(
            &index_path,
            concat!(
                "0,1,100,ACGTA,REF|ALT,ACGTA|CGTAA,1|1,ACGTA;GGAAA|CGTAA;ACGTA\n",
                "1,1,200,X,REF|ALT,X|Y,1|1,AAAAA|\n",
                "2,1,300,Z,REF,Z,1,AAA\n",
            ),
        )
        .expect("write index");
        let cases: [(&str, usize); 7] = [
            ("ACGTA", 2),
            ("GGAAA", 1),
            ("CGTAA", 1),
            ("AAAAA", 1),
            ("GTAAC", 0),
            ("TTTTT", 0),
            ("AAA", 0),
        ];
        // Tiny segment sizes force k-mers to complete across segment
        // boundaries; several thread counts exercise the parallel scan.
        for &(nthreads, segment_bases) in &[(1, 8_000_000), (1, 2), (4, 2), (4, 8_000_000)] {
            let counts = ref_counts(
                index_path.to_str().expect("path"),
                fasta_path.to_str().expect("path"),
                5,
                nthreads,
                segment_bases,
            )
            .expect("counts");
            for (kmer, expected) in cases {
                assert_eq!(
                    counts.count_of_kmer(kmer),
                    expected,
                    "kmer {} (threads {}, segment {})",
                    kmer,
                    nthreads,
                    segment_bases
                );
            }
            assert_eq!(
                counts.stats(),
                (4, 1),
                "threads {}, segment {}",
                nthreads,
                segment_bases
            );
        }
    }

    #[test]
    fn test_ref_counts_string_fallback_long_kmers() {
        let dir = tempfile::tempdir().expect("tempdir");
        let fasta_path = dir.path().join("ref.fasta");
        let s1 = "ACGTACGTACGTACGTACGTACGTACGTACGTA"; // 33 bases; canonical (rc starts T)
        let s2 = "GACTGACTGACTGACTGACTGACTGACTGACTG"; // 33 bases; absent from the reference
        std::fs::write(&fasta_path, format!(">a\n{}{}GGGGGGGG\n", s1, s1)).expect("write fasta");
        // Single-line record: 3 header bytes, one 74-base line, linewidth 75.
        std::fs::write(dir.path().join("ref.fasta.fai"), "a\t74\t3\t74\t75\n").expect("write fai");
        std::fs::write(
            dir.path().join("index.txt"),
            format!("0,1,100,S,REF|ALT,S|ALT2,1|1,{}|{}\n", s1, s2),
        )
        .expect("write index");
        for &nthreads in &[1, 2] {
            let counts = ref_counts(
                dir.path().join("index.txt").to_str().unwrap(),
                fasta_path.to_str().expect("path"),
                33,
                nthreads,
                8_000_000,
            )
            .expect("counts");
            assert!(matches!(counts, RefKmerCounts::Strings { .. }));
            assert_eq!(counts.count_of_kmer(s1), 2);
            assert_eq!(counts.count_of_kmer(s2), 0);
            assert_eq!(counts.stats(), (2, 1));
        }
    }
}
