
use crate::kvmer;
use crate::kvmer::*;
use crate::utils::*;
use crate::inference::*;
use crate::cmdline::AnalyzeArgs;

use std::collections::HashMap;
use log::{info, warn};

pub fn analyze(analyze_args: AnalyzeArgs) {
    
    let mut kvmer_set = KVmerSet::new(analyze_args.k, analyze_args.v, analyze_args.bidirectional);
    for file in &analyze_args.files {
        if is_fastx_file(file) {
            kvmer_set.add_file_to_kvmer_set(file, analyze_args.c);
        } else if is_sketch_file(file) {
            kvmer_set.load(file);
        } else {
            warn!("File format not recognized for file: {}. Skipping.", file);
        }
    }

    if let Some(reference) = &analyze_args.reference {

         let mut reference_kvmer_set = KVmerSet::new(analyze_args.k, analyze_args.v, true);
        reference_kvmer_set.add_file_to_kvmer_set(reference, analyze_args.c);

        let stats = kvmer_set.get_stats_with_reference(analyze_args.threshold, &reference_kvmer_set);
        //kvmer_set.output_stats(&stats);
        let spectrum = error_profile(&stats, analyze_args.bidirectional);
        output_error_spectrum(&spectrum, analyze_args.v);
        
    } else {
        //println!("Error rate: {}", kvmer_set.get_stats(analyze_args.threshold));
        let stats = kvmer_set.get_stats(analyze_args.threshold);
        //kvmer_set.output_stats(&stats);
        let spectrum = error_profile(&stats, analyze_args.bidirectional);
        output_error_spectrum(&spectrum, analyze_args.v);
       
    }

    
    
    /*
    let mut vmer_set = VmerSet::new(analyze_args.v, false);
    for file in &analyze_args.files {
        vmer_set.add_file_first_pass(file, analyze_args.c);
    }
    let relevant_values = vmer_set._get_relevant_values(analyze_args.threshold);
    let mut value_counts: HashMap<u64, u32> = HashMap::new();
    for file in &analyze_args.files {
        vmer_set.add_file_second_pass(file, &mut value_counts, &relevant_values);
    }
    vmer_set.add_value_counts(&value_counts);
    let stats = vmer_set.get_stats(analyze_args.threshold);

    let spectrum = error_profile(&stats, false);
    output_error_spectrum(&spectrum, analyze_args.v);
    */
}