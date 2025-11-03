use crate::kmer::*;

pub fn mapping(args: crate::cmdline::MappingArgs) {
    let mut kmer_set = KmerSet::new(args.k, true);
    kmer_set.add_file_to_kmer_set(&args.reference, args.c);

    let mut total_matched: u32 = 0;
    let mut total_kmers: u32 = 0;
    for file in &args.files {
        let (matched, total) = kmer_set.query_file(file, args.c, args.sample_rate, args.bidirectional, false);
        total_matched += matched;
        total_kmers += total;
    }
    println!("At k={}: estimated overall kmer match rate: {}/{} = {:.4}%", args.k, total_matched, total_kmers, total_matched as f64 / total_kmers as f64 * 100.);
}