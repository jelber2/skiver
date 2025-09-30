use log::warn;
use needletail::parse_fastx_file;

use std::collections::{HashMap, HashSet};

use crate::{seeding::*, types::EditOperation};


pub struct KVmerSet {
    pub key_size: u8,
    pub value_size: u8,
    pub key_value_map: HashMap<u64, HashMap<u64, u32>>,

    // utilities to extract key and value from a kmer hash
    key_mask: u64,
    value_mask: u64,
}


impl KVmerSet {
    pub fn new(key_size: u8, value_size: u8) -> Self {
        assert!(key_size + value_size <= 32, "Key size and value size must sum to at most 32");

        let v_mask = (1 << (value_size * 2)) - 1;
        let k_mask = ((1 << (key_size * 2)) - 1) << (value_size * 2);

        KVmerSet {
            key_size,
            value_size,
            key_value_map: HashMap::new(),
            key_mask: k_mask,
            value_mask: v_mask,
        }
    }

    

    fn _get_neighbors(&self, value: u64) -> HashMap<u64, EditOperation> {
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
                    neighbors.insert(neighbor, EditOperation::SUBSTITUTION);
                }
            }

            // Indels
            for &b in &bases {
                let left_part = (value >> (shift + 2)) << ((shift + 2));
                let right_part = value & ((1 << (shift + 2)) - 1);
                let neighbor_insert = left_part | (b << shift) | (right_part >> 2);
                neighbors.entry(neighbor_insert)
                    .and_modify(|op|
                        if *op != EditOperation::INSERTION {
                            *op = EditOperation::AMBIGUOUS
                        }
                    )
                    .or_insert(EditOperation::INSERTION);

                let right_part = value & ((1 << shift) - 1);
                let neighbor_delete = left_part | (right_part << 2) | b;
                neighbors.entry(neighbor_delete)
                    .and_modify(|op| 
                        if *op != EditOperation::DELETION {
                            *op = EditOperation::AMBIGUOUS
                        }
                    )
                    .or_insert(EditOperation::DELETION);
            }
            
        }

        neighbors

    }



    pub fn add_seed_vector(&mut self, seed_vec: &[u64]) {
        for &kmer in seed_vec {
            let key = (kmer & self.key_mask) >> self.value_size;
            let value = kmer & self.value_mask;

            let entry = self.key_value_map.entry(key).or_insert_with(HashMap::new);
            let count = entry.entry(value).or_insert(0);
            *count += 1;
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


    pub fn consistency_index(&self, threshold: u32) -> f64 {
        // for each key, see by how much the values agree
        let mut total_majority = 0;
        let mut total_values = 0;

        for (_key, value_map) in &self.key_value_map {
            let mut max_count = 0;
            let mut sum_count = 0;

            //println!("Value map: {:?}", value_map);

            for (_value, count) in value_map {
                sum_count += *count;
                if *count > max_count {
                    max_count = *count;
                }
            }

            if sum_count <= threshold {
                continue;
            }

            total_majority += max_count;
            total_values += sum_count;
        }

        if total_values == 0 {
            0.0
        } else {
            total_majority as f64 / total_values as f64
        }
    }

    pub fn error_profile(&self, threshold: u32) -> HashMap<EditOperation, f64> {
        // for each key, compute the error rate as 1 - (max_count / total_count)
        let mut error_profile: HashMap<EditOperation, f64> = HashMap::new();
        let mut num_errors: u32 = 0;
        let mut num_consensus: u32 = 0;

        for (_key, value_map) in &self.key_value_map {
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

            num_consensus += max_count;

            // for each non-consensus value, determine if it is a substitution, insertion, or deletion
            // relative to the consensus value
            let neighbors = self._get_neighbors(max_value);
            //println!("Neighbors of {}: {:?}", max_value, neighbors);
            for (value, count) in value_map {
                if *value != max_value {
                    if let Some(op) = neighbors.get(value) {
                        match op {
                            EditOperation::SUBSTITUTION => {
                                error_profile.entry(EditOperation::SUBSTITUTION)
                                    .and_modify(|e| *e += *count as f64)
                                    .or_insert(*count as f64);
                            },
                            EditOperation::INSERTION => {
                                error_profile.entry(EditOperation::INSERTION)
                                    .and_modify(|e| *e += *count as f64)
                                    .or_insert(*count as f64);
                            },
                            EditOperation::DELETION => {
                                error_profile.entry(EditOperation::DELETION)
                                    .and_modify(|e| *e += *count as f64)
                                    .or_insert(*count as f64);
                            },
                            EditOperation::AMBIGUOUS => {},
                        }
                    }
                    
                    num_errors += *count;
                }
            }
        }

        // normalize the error profile
        for (_op, count) in error_profile.iter_mut() {
            *count /= (num_errors + num_consensus) as f64;
        }

        error_profile
    }

    

    
}

pub fn add_file_to_kvmer_set(
    kvmer: &mut KVmerSet,
    seq_file: &str,
    seq_sketch: &SequencesSketch,
    bidirectional: bool,
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
                    let (kmer_vec, kmer_vec_rev) = seq_to_bidirectional_kmer_vec(seq.as_ref(), seq_sketch.k());
                    kvmer.add_seed_vector(&seq_sketch.sketch_kmer_vec(&kmer_vec));
                    if bidirectional {
                        // Also sketch the reverse complement if bidirectional
                        kvmer.add_seed_vector(&seq_sketch.sketch_kmer_vec(&kmer_vec_rev));
                    }
                }
                Err(e) => {
                    warn!("Error reading record: {}", e);
                }
            }
        }
    }
}