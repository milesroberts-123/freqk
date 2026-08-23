use rust_htslib::bcf::{Format, Writer};
use rust_htslib::bcf::header::Header;
use rust_htslib::bcf::record::GenotypeAllele;
use rand::prelude::*;
use rand_distr::Zeta;
use rand::distr::StandardUniform;

pub fn simulate_workflow(length: &usize, snp: &f64, shape: &f64) -> () {
    // sample distribution
    let val: f64 = rand::rng().sample(Zeta::new(1.5).unwrap());
    println!("{}", val);
    // sample uniform distribution
    let val: f32 = rand::rng().sample(StandardUniform);
    println!("f32 from [0, 1): {}", val);

    // Generate random reference genome

    // Generate random mutations
    let mut rng = rand::rng();
    let mut die_range: rand_distr::Iter<_, &mut _, f64> = StandardUniform.sample_iter(&mut rng);
    for i in 1..*length{
        if die_range.next().unwrap() < *snp {
            println!("SNP created at pos {}", i)
        }
    }

    // Distribute mutations across a population
    // to create a VCF

    // Create minimal VCF header with a single contig and a single sample
    let mut header = Header::new();
    let header_contig_line = r#"##contig=<ID=1,length=10>"#;
    header.push_record(header_contig_line.as_bytes());
    let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
    header.push_record(header_gt_line.as_bytes());
    header.push_sample("test_sample".as_bytes());

    // Write uncompressed VCF to stdout with above header and get an empty record
    let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    let mut record = vcf.empty_record();

    // Set chrom and pos to 1 and 7, respectively - note the 0-based positions
    let rid = vcf.header().name2rid(b"1").unwrap();
    record.set_rid(Some(rid));
    record.set_pos(6);

    // Set record genotype to 0|1 - note first allele is always unphased
    let alleles = &[GenotypeAllele::Unphased(0)];
    record.push_genotypes(alleles).unwrap();

    // Write record
    vcf.write(&record).unwrap()
}
