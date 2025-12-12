

use crate::kvmer::*;
use simple_logger::SimpleLogger;
use log::info;
use crate::cmdline::SketchArgs;
use rayon::prelude::*;


pub fn sketch(args: SketchArgs) {
    SimpleLogger::new().with_level(log::LevelFilter::Info).init().unwrap();
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    info!("Processing query files...");
    
    let mut kvmer_set = KVmerSet::new(args.k, args.v, false);
    for file in &args.files {
        kvmer_set.add_file_to_kvmer_set(file, args.c, args.trim_front, args.trim_back);
    }
    info!("Finished processing query files.");

    kvmer_set.dump(&args.output_path);
    info!("Sketch saved to {}", args.output_path);
}