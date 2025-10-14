
use crate::kvmer::*;
use crate::cmdline::AnalyzeArgs;

pub fn analyze(analyze_args: AnalyzeArgs) {
    let mut kvmer_set = KVmerSet::new(analyze_args.k, analyze_args.v, false);
    for file in &analyze_args.files {
        kvmer_set.add_file_to_kvmer_set(file, analyze_args.c);
    }

    //println!("Error rate: {}", kvmer_set.get_stats(analyze_args.threshold));
    let stats = kvmer_set.get_stats(analyze_args.threshold);
    kvmer_set.output_stats(&stats);
}