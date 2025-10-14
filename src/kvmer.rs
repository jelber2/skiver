use log::warn;
use needletail::parse_fastx_file;

use std::collections::HashMap;

use crate::{seeding::*, types::*};


pub struct KVmerSet {
    pub key_size: u8,
    pub value_size: u8,
    pub kv_size: u8,
    pub key_value_map: HashMap<u64, HashMap<u64, u32>>,

    // utilities to extract key and value from a kmer hash
    key_mask: u64,
    value_mask: u64,

    // whether both forward and reverse complement of the reads are included
    bidirectional: bool,
}


impl KVmerSet {
    pub fn new(key_size: u8, value_size: u8, bidirectional: bool) -> Self {
        assert!(key_size + value_size <= 32, "Currently, we only support k + v <= 32.");

        let v_mask = (1 << (value_size * 2)) - 1;
        let k_mask = ((1 << (key_size * 2)) - 1) << (value_size * 2);

        KVmerSet {
            key_size,
            value_size,
            kv_size: key_size + value_size,
            key_value_map: HashMap::new(),
            key_mask: k_mask,
            value_mask: v_mask,
            bidirectional,
        }
    }
    

    pub fn _get_neighbors(&self, value: u64) -> HashMap<u64, EditOperation> {
        // get all the values with edit distance 1 from the input value
        
        let mut neighbors: HashMap<u64, EditOperation> = HashMap::new();
        let bases = [0, 1, 2, 3]; // A, C, G, T

        for i in 0..self.value_size {
            let shift = i * 2;
            
            // Substitutions
            for &b in &bases {
                let current_base = (value >> shift) & 0b11;
                if b != current_base {
                    let neighbor = (value & !(0b11 << shift)) | (b << shift);
                    if self.bidirectional {
                        neighbors.insert(neighbor, BASES_TO_SUBSTITUTION_CANONICAL[current_base as usize][b as usize].unwrap());
                    } else {
                        neighbors.insert(neighbor, BASES_TO_SUBSTITUTION[current_base as usize][b as usize].unwrap());
                    }
                    
                }
            }

            // Indels
            for &b in &bases {
                if shift == 0 && b == (value >> shift) & 0b11 {
                    continue; // skip the original base for the first position
                }
                
                let left_part = (value >> (shift + 2)) << ((shift + 2));
                let right_part = value & ((1 << (shift + 2)) - 1);
                let neighbor_insert = left_part | (b << shift) | (right_part >> 2);
                if self.bidirectional {
                    neighbors.entry(neighbor_insert)
                    .and_modify(|op|
                        if *op != BASES_TO_INSERTION_CANONICAL[b as usize].unwrap() {
                            *op = EditOperation::AMBIGUOUS
                        }
                    )
                    .or_insert(BASES_TO_INSERTION_CANONICAL[b as usize].unwrap());
                } else {
                    neighbors.entry(neighbor_insert)
                    .and_modify(|op|
                        if *op != BASES_TO_INSERTION[b as usize].unwrap() {
                            *op = EditOperation::AMBIGUOUS
                        }
                    )
                    .or_insert(BASES_TO_INSERTION[b as usize].unwrap());
                }
                
                
                
                let right_part = value & ((1 << shift) - 1);
                let neighbor_delete = left_part | (right_part << 2) | b;
                let original_base = (value >> shift) & 0b11;
                if self.bidirectional {
                    neighbors.entry(neighbor_delete)
                    .and_modify(|op| 
                        if *op != BASES_TO_DELETION_CANONICAL[original_base as usize].unwrap() {
                            *op = EditOperation::AMBIGUOUS
                        }
                    )
                    .or_insert(BASES_TO_DELETION_CANONICAL[original_base as usize].unwrap());
                } else {
                    neighbors.entry(neighbor_delete)
                    .and_modify(|op| 
                        if *op != BASES_TO_DELETION[original_base as usize].unwrap() {
                            *op = EditOperation::AMBIGUOUS
                        }
                    )
                    .or_insert(BASES_TO_DELETION[original_base as usize].unwrap());
                }
            }
        }

        neighbors

    }



    pub fn to_value_string(&self, kmer: u64) -> String {
        // for debugging: convert a kmer to a string

        let mut s = Vec::with_capacity(self.value_size as usize);
        for i in (0..self.value_size).rev() {
            let shift = i * 2;
            let base = ((kmer >> shift) & 0b11) as usize;
            s.push(crate::types::SEQ_TO_BYTE[base]);
        }
        String::from_utf8(s).unwrap()
    }

    pub fn to_key_string(&self, kmer: u64) -> String {
        // for debugging: convert a kmer to a string

        let mut s = Vec::with_capacity(self.key_size as usize);
        for i in (0..self.key_size).rev() {
            let shift = i * 2;
            let base = ((kmer >> shift) & 0b11) as usize;
            s.push(crate::types::SEQ_TO_BYTE[base]);
        }
        String::from_utf8(s).unwrap()
    }

    pub fn show_neighbors(&self, value: u64) {
        // for debugging: print all the neighbors of a value

        let neighbors = self._get_neighbors(value);
        for (neighbor, op) in neighbors {
            println!("Neighbor: {}, Operation: {:?}", self.to_value_string(neighbor), op);
        }
    }

    pub fn homopolymer_length(&self, key: u64, value: u64) -> u32 {
        let mut longest_homopolymer: u32 = 1;
        let mut current_homopolymer: u32 = 1;
        
        // Find the longest homopolymer at the end of the key
        let mut last_base = key & 0b11;
        for i in 1..self.key_size {
            let shift = i * 2;
            let base = (key >> shift) & 0b11;
            if base == last_base {
                current_homopolymer += 1;
            } else {
                break;
            }
        }
        // Extend the homopolymer into the value
        for i in (0..self.value_size).rev() {
            let shift = i * 2;
            let base = (value >> shift) & 0b11;
            if base == last_base {
                current_homopolymer += 1;
            } else {
                if current_homopolymer > longest_homopolymer {
                    longest_homopolymer = current_homopolymer;
                }
                current_homopolymer = 1;
                last_base = base;
            }
        }

        if current_homopolymer > longest_homopolymer {
            longest_homopolymer = current_homopolymer;
        }

        longest_homopolymer
    }



    pub fn add_seed_vector(&mut self, seed_vec: &[u64]) {
        for &kmer in seed_vec {
            let key = (kmer & self.key_mask) >> (self.value_size * 2);
            let value = kmer & self.value_mask;

            let entry = self.key_value_map.entry(key).or_insert_with(HashMap::new);
            let count = entry.entry(value).or_insert(0);
            *count += 1;
        }
    }


    fn extract_markers_masked(&self, string: &[u8], kmer_vec: &mut Vec<u64>, c: usize) {
        // extract sketched kv-mers from the given sequence string
        #[cfg(any(target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                use crate::avx2_seeding::*;
                unsafe {
                    extract_markers_avx2_masked(string, kmer_vec, c, self.kv_size as usize, self.key_mask as i64, self.bidirectional);
                }
            } else {
                fmh_seeds_masked(string, kmer_vec, c, self.kv_size as usize, self.key_mask, self.bidirectional);
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            fmh_seeds_masked(string, kmer_vec, c, self.kv_size as usize, self.key_mask, self.bidirectional);
        }
    }

    pub fn add_file_to_kvmer_set(
        &mut self,
        seq_file: &str,
        c: usize,
    ) {
        let reader = parse_fastx_file(&seq_file);
        //println!("Reading file: {}", seq_file);
        if !reader.is_ok() {
            //println!("Not OK Reading file: {}", seq_file);
            println!("{} is not a valid fasta/fastq file; skipping.", seq_file);
        } else {
            //println!("Reading file: {}", seq_file);
            let mut reader = reader.unwrap();
            while let Some(record) = reader.next() {
                match record {
                    Ok(record) => {
                        let seq = record.seq();
                        let mut kmer_vec: Vec<u64> = Vec::new();
                        self.extract_markers_masked(seq.as_ref(), &mut kmer_vec, c);
                        self.add_seed_vector(&kmer_vec);
                    }
                    Err(e) => {
                        warn!("Error reading record: {}", e);
                    }
                }
            }
        }
    }

    pub fn containment_index(&self, other: &KVmerSet) -> (f64, f64) {
        // check the key containment index and key-value pair containment index
        // each key/ key-value pair is counted once
        let mut shared_keys = 0;
        let mut shared_key_values = 0;
        let mut total_key_values = 0;

        for (key, value_map) in &self.key_value_map {
            if let Some(other_value_map) = other.key_value_map.get(key) {
                shared_keys += 1;

                for (value, _count) in value_map {
                    if let Some(_other_count) = other_value_map.get(value) {
                        shared_key_values += 1;
                    }
                }
            } 
            total_key_values += value_map.len();
        }

        let key_containment = if self.key_value_map.is_empty() {
            0.0
        } else {
            shared_keys as f64 / self.key_value_map.len() as f64
        };

        let key_value_containment = if total_key_values == 0 {
            0.0
        } else {
            shared_key_values as f64 / total_key_values as f64
        };

        (key_containment, key_value_containment)
    }

    pub fn get_stats(&self, threshold: u32) -> KVmerStats {
        // record the keys and values for output
        let mut keys: Vec<u64> = Vec::new();
        let mut consensus_values: Vec<u64> = Vec::new();

        // count the number of consensus and error kmers
        let mut consensus_counts: Vec<u32> = Vec::new();
        let mut error_counts: Vec<HashMap<EditOperation, u32>> = Vec::new();
        let mut total_counts: Vec<u32> = Vec::new();
        let mut neighbor_counts: Vec<u32> = Vec::new();

        for (key, value_map) in &self.key_value_map {
            let mut max_count = 0;
            let mut sum_count = 0;
            let mut max_value: u64 = 0;

            // find the consensus value
            for (value, count) in value_map {
                sum_count += *count;
                if *count > max_count {
                    max_count = *count;
                    max_value = *value;
                }
            }

            // skip low coverage keys
            if sum_count <= threshold {
                continue;
            }

            // for each non-consensus value, determine if it is a substitution, insertion, or deletion
            // relative to the consensus value
            let mut error_count_map: HashMap<EditOperation, u32> = HashMap::new();
            let neighbors = self._get_neighbors(max_value);
            if neighbors.contains_key(&max_value) {
                // This would confound the X=0 case
                continue;
            }
            let mut num_neighbors = 0;
            //println!("Neighbors of {}: {:?}", max_value, neighbors);
            for (value, count) in value_map {
                if *value != max_value && neighbors.contains_key(value) {
                    let op = neighbors.get(value).unwrap();
                    let entry = error_count_map.entry(*op).or_insert(0);
                    *entry += *count;
                    num_neighbors += count;
                }
            }

            // update the vectors
            keys.push(*key);
            consensus_values.push(max_value);
            consensus_counts.push(max_count);
            error_counts.push(error_count_map);
            total_counts.push(sum_count);
            neighbor_counts.push(num_neighbors);
        }

        KVmerStats {
            k: self.key_size,
            v: self.value_size,
            keys,
            consensus_values,
            consensus_counts,
            total_counts,
            neighbor_counts,
            error_counts,
        }
    /*

        if total_counts.is_empty() {
            0.0
        } else {
            // compute the overall error rate
            let p_0: f64 = consensus_counts.iter().sum::<u32>() as f64 / total_counts.iter().sum::<u32>() as f64;

            let p_1: f64 = error_counts.iter().sum::<u32>() as f64 / total_counts.iter().sum::<u32>() as f64;

            let mut p_1_over_p_0 = Vec::new();
            for i in 0..consensus_counts.len() {
                p_1_over_p_0.push(error_counts[i] as f64 / (consensus_counts[i] as f64 * (self.value_size as f64)));
            }

            //println!("{:?}", p_1_over_p_0);

            // p_1 / p_0 estimator
            //let error_rate = p_1 / p_0 / (self.value_size as f64);

            // p_0 estimator
            //let error_rate = - p_0.ln() / self.value_size as f64;
            println!("Number of consensus kmers: {}", consensus_counts.len());

            //p_1_over_p_0.sort_by(|a, b| a.partial_cmp(b).unwrap());
            //let error_rate = p_1_over_p_0[p_1_over_p_0.len() / 2];// / (self.value_size as f64);

            //let error_rate = p_1_over_p_0.iter().sum::<f64>() / (p_1_over_p_0.len() as f64);

            
            /*
            let mut p_0_vec = Vec::new();
            for i in 0..consensus_counts.len() {
                p_0_vec.push(consensus_counts[i] as f64 / total_counts[i] as f64);
            }
            p_0_vec.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let p_0_mean = p_0_vec.iter().sum::<f64>() / (p_0_vec.len() as f64);

            


            let error_rate = 1.0 - p_0_mean.powf(1.0 / self.value_size as f64);
            */
            
            error_rate
        }
            */

    }

    pub fn output_stats(&self, stats: &KVmerStats) {
        print!("key\tconsensus_value\thomopolymer_length\tconsensus_count\tneighbor_count\ttotal_count\t");
        if self.bidirectional {
            for op in ALL_OPERATIONS_CANONICAL {
                print!("\t{:?}", op);
            }
        } else {
            for op in ALL_OPERATIONS {
                print!("\t{:?}", op);
            }
        }
        print!("\n");


        for i in 0..stats.keys.len() {
            print!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                self.to_key_string(stats.keys[i]),
                self.to_value_string(stats.consensus_values[i]),
                self.homopolymer_length(stats.keys[i], stats.consensus_values[i]),
                stats.consensus_counts[i],
                stats.neighbor_counts[i],
                stats.total_counts[i],
            );
            if self.bidirectional {
                
                for op in ALL_OPERATIONS_CANONICAL {
                    let value = *(stats.error_counts[i]).get(&op).unwrap_or(&0);
                    print!("\t{}", value);
                }
            } else {
                for op in ALL_OPERATIONS {
                    let value = *(stats.error_counts[i]).get(&op).unwrap_or(&0);
                    print!("\t{}", value);
                }
            }
            print!("\n")
        }
    }

}