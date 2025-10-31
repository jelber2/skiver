use crate::kmer::*;

pub fn mapping(args: crate::cmdline::MappingArgs) {
    let mut kmer_set = KmerSet::new(args.k, true);
    kmer_set.add_file_to_kmer_set(&args.reference, args.c);

    for file in &args.files {
        kmer_set.query_file(file, args.c, args.sample_rate, args.bidirectional);
    }
}