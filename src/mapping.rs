use simple_logger::SimpleLogger;
use log::{info, warn};
use needletail::parse_fastx_file;
use serde::{Serialize, Deserialize};

use std::collections::HashMap;
use std::io::Write;

use crate::{seeding::*};

pub fn map(args: crate::cmdline::MapArgs) {
    SimpleLogger::new().with_level(log::LevelFilter::Info).init().unwrap();

    let mut kmer_set = KmerSet::new(args.k, true);
    kmer_set.add_file_to_kmer_set(&args.reference, args.c, args.trim_front, args.trim_back);
    info!("Loaded reference file: {}", args.reference);

    let mut total_matched: u32 = 0;
    let mut total_kmers: u32 = 0;
    info!("Processing query files...");
    for file in &args.files {
        let (matched, total) = kmer_set.query_file(file, args.c, args.lower_bound, args.sample_rate, args.bidirectional, args.print_verbose, args.trim_front, args.trim_back);

        total_matched += matched;
        total_kmers += total;

    }
    info!("Finished processing query files.");
    info!("At k={}: estimated overall kmer match rate: {}/{} = {:.4}%", args.k, total_matched, total_kmers, total_matched as f64 / total_kmers as f64 * 100.);
}


#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct KmerSet {
    pub key_size: u8,
    pub num_kmers: u32,
    pub kmer_map: HashMap<u64, u32>,

    // whether both forward and reverse complement of the reads are included
    bidirectional: bool,
}

/**
 * Utility class for testing read accuracy
 * used only for testing purposes
 */
impl KmerSet {
    pub fn new(key_size: u8, bidirectional: bool) -> Self {
        assert!(key_size <= 32, "Currently, we only support k <= 32.");

        KmerSet {
            key_size,
            num_kmers: 0,
            kmer_map: HashMap::new(),
            bidirectional,
        }
    }


    pub fn add_seed_vector(&mut self, seed_vec: &[u64]) {
        for &kmer in seed_vec {
            let entry = self.kmer_map.entry(kmer).or_insert(0);
            *entry += 1;
        }
        self.num_kmers += seed_vec.len() as u32;
    }

    pub fn query_seed_vector(&self, seed_vec: &[u64]) -> (u32, u32) {
        let mut count: u32 = 0;
        for &kmer in seed_vec {
            if self.kmer_map.contains_key(&kmer) {
                count += 1;
            }
        }
        (count, seed_vec.len() as u32)
    }


    fn extract_markers_masked(&self, string: &[u8], kmer_vec: &mut Vec<u64>, _value_vec: &mut Vec<u64>, c: usize, bidirectional: bool, trim_front: usize, trim_back: usize) {
        let start = std::cmp::min(trim_front, string.len());
        let end = string.len().saturating_sub(trim_back);

        // extract sketched kv-mers from the given sequence string
        #[cfg(any(target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                use crate::avx2_seeding::*;
                unsafe {
                    extract_markers_avx2_masked(&string[start..end], kmer_vec, _value_vec, c, self.key_size as usize, 0 as usize, bidirectional);
                }
            } else {
                fmh_seeds_masked(&string[start..end], kmer_vec, _value_vec, c, self.key_size as usize, 0 as usize, bidirectional);
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            fmh_seeds_masked(&string[start..end], kmer_vec, _value_vec, c, self.key_size as usize, 0 as usize, bidirectional);
        }
    }

    pub fn add_file_to_kmer_set(
        &mut self,
        seq_file: &str,
        c: usize,
        trim_front: usize,
        trim_back: usize,
    ) {
        let reader = parse_fastx_file(&seq_file);
        //println!("Reading file: {}", seq_file);
        if !reader.is_ok() {
            //println!("Not OK Reading file: {}", seq_file);
            warn!("{} is not a valid fasta/fastq file; skipping.", seq_file);
        } else {
            //println!("Reading file: {}", seq_file);
            let mut reader = reader.unwrap();
            while let Some(record) = reader.next() {
                match record {
                    Ok(record) => {
                        let seq = record.seq();
                        let mut kmer_vec: Vec<u64> = Vec::new();
                        let mut _value_vec: Vec<u64> = Vec::new();
                        self.extract_markers_masked(seq.as_ref(), &mut kmer_vec, &mut _value_vec, c, self.bidirectional, trim_front, trim_back);
                        self.add_seed_vector(&kmer_vec);
                    }
                    Err(e) => {
                        warn!("Error reading record: {}", e);
                    }
                }
            }
        }
    }

    /**
     * Estimate the kmer match rate for the given read file.
     * The rate is calculated by the average kmer match rate across all reads,
     * excluding reads that have zero matched kmers.
     */
    pub fn query_file(
        &self,
        seq_file: &str,
        c: usize,
        threshold: u32,
        sample_per_num_read: usize,
        bidirectional: bool,
        print_verbose: bool,
        trim_front: usize,
        trim_back: usize,
    ) -> (u32, u32) {
        let reader = parse_fastx_file(&seq_file);

        //let matched_kmers: Vec<u32> = Vec::new();
        //let total_kmers: Vec<u32> = Vec::new();

        let mut matched_kmers: u32 = 0;
        let mut total_kmers: u32 = 0;

        let mut read_count: usize = 0;
        //println!("Reading file: {}", seq_file);
        if print_verbose {
            println!("total,matched");
        }
        if !reader.is_ok() {
            //println!("Not OK Reading file: {}", seq_file);
            warn!("{} is not a valid fasta/fastq file; skipping.", seq_file);
        } else {
            //println!("Reading file: {}", seq_file);
            let mut reader = reader.unwrap();
            while let Some(record) = reader.next() {
                match record {
                    Ok(record) => {

                        read_count += 1;
                        if read_count % sample_per_num_read != 0 {
                            continue;
                        }

                        let seq = record.seq();
                        let mut kmer_vec: Vec<u64> = Vec::new();
                        let mut _value_vec: Vec<u64> = Vec::new();
                        self.extract_markers_masked(seq.as_ref(), &mut kmer_vec, &mut _value_vec, c, bidirectional, trim_front, trim_back);
                        let (matched, total) = self.query_seed_vector(&kmer_vec);
                        //matched_kmers.push(matched);
                        //total_kmers.push(total);
                        if print_verbose {
                            println!("{},{}", total, matched);
                        }
                        if matched >= threshold {
                            matched_kmers += matched;
                            total_kmers += total;
                        }
                    }
                    Err(e) => {
                        warn!("Error reading record: {}", e);
                    }
                }
            }
            //println!("Total kmers: {}, Matched kmers: {}, Match rate: {:.4}%", total_kmers, matched_kmers, matched_kmers as f64 / total_kmers as f64 * 100.);
        }
        (matched_kmers, total_kmers)
    }

}


