
use crate::kvmer;
use crate::kvmer::*;
use crate::inference::*;
use crate::cmdline::SketchArgs;

use std::collections::HashMap;

pub fn sketch(args: SketchArgs) {
    
    let mut kvmer_set = KVmerSet::new(args.k, args.v, false);
    for file in &args.files {
        kvmer_set.add_file_to_kvmer_set(file, args.c);
    }

    kvmer_set.dump(&args.output_path);
}