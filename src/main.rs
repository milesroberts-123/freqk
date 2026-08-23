use clap::{Parser, Subcommand};
use itertools::izip;

mod common;
mod index;
mod dedup;
mod count;
mod call;
mod simulate;
mod hetmers;

#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
    #[command(flatten)]
    verbosity: clap_verbosity_flag::Verbosity,
}


#[derive(Subcommand, Debug)]
enum Commands {
    /// Get k-mers specific to each allele of each variant
    Index {
        #[arg(short,long, help = "fasta file of reference genome")]
        fasta: String,
        #[arg(long, help = "vcf file of variations between reference and other genomes")]
        vcf: String,
        #[arg(short,long, help = "name of the index file to be output")]
        output: String,
        #[arg(short,long, help = "kmer length for building the index")]
        kmer: i64,
        #[command(flatten)]
        verbosity: clap_verbosity_flag::Verbosity,
    },
    /// Deduplicate index of k-mers shared across variants
    VarDedup {
        #[arg(short,long, help = "name of index file")]
        index: String,
        #[arg(short,long, help = "name of deduplicated index")]
        output: String,
        #[command(flatten)]
        verbosity: clap_verbosity_flag::Verbosity,
    },
    /// Count k-mers by allele
    Count {
        #[arg(short,long, help = "name of index file")]
        index: String,
        #[arg(short,long, help = "comma-separated list of reads in fastq format (can be gz or not)", value_delimiter = ',', value_name = "FILE")]
        reads: String,
        #[arg(short,long, help = "Number of threads to use (a copy of index will be loaded onto each thread)")]
        nthreads: usize,
        #[arg(short,long, help = "name of output for counts by allele")]
        freq_output: String,
        #[arg(short,long, help = "name of output for raw kmer counts")]
        count_output: String,
        #[command(flatten)]
        verbosity: clap_verbosity_flag::Verbosity,
    },
    /// Convert counts by allele into allele frequencies
    Call {
        #[arg(short,long, help = "name of index file")]
        index: String,
        #[arg(short,long, help = "name of file for counts by allele")]
        counts: String,
        #[arg(short,long, help = "name of output file with allele frequencies")]
        output: String,
        #[command(flatten)]
        verbosity: clap_verbosity_flag::Verbosity,
    },
    /// Deduplicate index of reference k-mers 
    RefDedup {
        #[arg(short,long, help = "path to index file")]
        index: String,
        #[arg(short,long, help = "path to deduped index file")]
        output: String,
        #[arg(short,long, help = "fasta file of reference genome")]
        fasta: String,
        #[arg(long, help = "vcf file of variations between reference and other genomes")]
        vcf: String,
        #[command(flatten)]
        verbosity: clap_verbosity_flag::Verbosity,
    },
    Simulate {
        #[arg(short, long, help = "genome length to simulate")]
        length: usize,
        #[arg(short, long, help = "Mutation rate")]
        snp: f64,
        #[arg(short, long, help = "Zeta distribution shape")]
        shape: f64,
        #[command(flatten)]                                                                                                         verbosity: clap_verbosity_flag::Verbosity,
    },
    /// Count het-mers
    Hetmers {
        /// kmer count table file name
        #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
        inputs: Vec<String>,
        /// prefix for output files
        #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
        outputs: Vec<String>,
        /// minimum k-mer count
        #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
        minimums: Vec<usize>,
        /// number of alleles in each hetmer
        #[arg(short = 'l', long, default_value_t = 2)]
        alleles: usize,
        /// mean k-mer coverage
        #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
        coverages: Vec<f64>,
        /// pool size
        #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
        pools: Vec<i32>,
        /// shape parameter for prior distribution
        #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
        alphas: Vec<f64>,
        /// shape parameter for prior distribution
        #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
        betas: Vec<f64>,
        /// thresholds for determining if k-mer has abnormal copy number
        #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
        sigmas: Vec<f64>,
        #[command(flatten)]
        verbosity: clap_verbosity_flag::Verbosity,
    }
}

// bring it all together!
fn main() {
    let cli = Cli::parse();

    env_logger::Builder::new()
        .filter_level(cli.verbosity.log_level_filter())
        .init();

    match &cli.command {
        Commands::Index { fasta, vcf, output, kmer, verbosity } => {
            verbosity.log_level_filter();
            log::info!("fasta: {}, vcf: {}, k: {}", fasta, vcf, kmer);
            index::index_workflow(vcf, fasta, output, kmer);
        }
        Commands::Count { index, reads, nthreads, freq_output, count_output, verbosity } => {
            verbosity.log_level_filter();
            log::info!("Counting k-mers: INDEX: {}, READS: {:?}, NTHREADS: {}, FREQ OUTPUT: {}, COUNT OUTPUT: {}", index, reads, nthreads, freq_output, count_output);
            count::count_workflow(index, reads, *nthreads, freq_output, count_output);
        }
        Commands::VarDedup { index, output, verbosity } => {
            verbosity.log_level_filter();
            log::info!("Deduplicating index across variants: INDEX: {} OUTPUT: {}", index, output);
            let _ = dedup::find_dup_kmers_across_var(index, output);
        }
        Commands::RefDedup { index, output, fasta, vcf, verbosity } => {
            verbosity.log_level_filter();
            log::info!("Deduplicating index of reference k-mers: INDEX: {} OUTPUT: {} FASTA: {} VCF: {}", index, output, fasta, vcf);
            let ref_hash = dedup::reference_hashset(index, fasta, vcf);
            let _ = dedup::remove_ref_kmers(index, output, ref_hash);
        }
        Commands::Call { index, counts, output, verbosity } => {
            verbosity.log_level_filter();
            log::info!("Converting counts to allele frequencies: {} {} {}", index, counts, output);
            let _ = call::call_from_counts(index, counts, output);
        }
        Commands::Simulate { length, snp, shape, verbosity } => {
            verbosity.log_level_filter();
            log::info!("Simulating a genome with structural variants.");
            let _ = simulate::simulate_workflow(length, snp, shape);
        }
        Commands::Hetmers { inputs, outputs, minimums, alleles, coverages, pools, alphas, betas, sigmas, verbosity } => {
            verbosity.log_level_filter();
            log::info!("Finding hetmers in k-mer counts");
            for (input, output, minimum, coverage, pool, alpha, beta, sigma) in izip!(inputs, outputs, minimums, coverages, pools, alphas, betas, sigmas) {
                hetmers::kmers_to_hetmers(input, output, *minimum, *alleles, *pool, *coverage, *alpha, *beta, *sigma);
            }
        }
    }
}
