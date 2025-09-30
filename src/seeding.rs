// This file contains multiple implementations from sylph (https://github.com/bluenote-1577/sylph). Below is their license.

/*
MIT License

Copyright (c) 2023 Jim Shaw

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

use crate::types::*;
//use needletail::parse_fastx_file;


/**
 * @brief A fast hash function for k-mer represented using 64-bit integer
 * Adopted from https://github.com/bluenote-1577/sylph/blob/main/src/seeding.rs
 */
#[inline]
pub fn mm_hash64(kmer: u64) -> u64 {
    let mut key = kmer;
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    key
}


/**
 * @brief Class to compute the seeds for a given sequence
 */
pub enum SketchingMethod {
    Minimizer { window_size: usize },
    FracMinHash { fraction: f64 },
    MaskedFracMinHash { fraction: f64, mask_size: usize },
}

pub struct SequencesSketch {
    k: usize,
    method: SketchingMethod,
    hash_fn: Option<fn(u64) -> u64>,
}

impl SequencesSketch {
    pub fn new(k: usize, method: SketchingMethod, hash_fn: Option<fn(u64) -> u64>) -> Self {
        if hash_fn.is_none() {
            SequencesSketch { k, method, hash_fn: Some(mm_hash64) }
        } else {
            SequencesSketch { k, method, hash_fn }
        }
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn method(&self) -> String {
        match self.method {
            SketchingMethod::Minimizer { window_size } => format!("Minimizer (window size: {})", window_size),
            SketchingMethod::FracMinHash { fraction } => format!("FracMinHash (fraction: {:.2})", fraction),
            SketchingMethod::MaskedFracMinHash { fraction, mask_size } => format!("MaskedFracMinHash (fraction: {:.2}, mask size: {})", fraction, mask_size),
        }
    }


    fn extract_minimizers(&self, kmer_vec: &[u64], window_size: usize) -> Vec<u64> {
        // [TODO] check the correctness of this function
        // [TODO] add option to use a hash function

        if kmer_vec.len() < window_size {
            return Vec::new();
        }

        let mut minimizers = Vec::new();
        let mut window = Vec::new();

        // initialize the window
        for i in 0..window_size {
            window.push(kmer_vec[i]);
        }
        // find the minimum in the window
        let mut min = window.iter().cloned().min().unwrap();
        minimizers.push(min);

        // slide the window
        for i in 0..=(kmer_vec.len() - window_size - 1) {
            // remove the first element of the window
            let removed = window.remove(0);
            // add the next element to the window
            window.push(kmer_vec[i + window_size]);
            // find the minimum in the window
            if kmer_vec[i + window_size] < min {
                min = kmer_vec[i + window_size];
                minimizers.push(min);
            } else if removed == min {
                // if the removed element was the minimum, find the new minimum
                min = window.iter().cloned().min().unwrap();
                minimizers.push(min);
            }
        }

        minimizers
    }

    fn extract_frac_min_hash(&self, kmer_vec: &[u64], fraction: f64) -> Vec<u64> {
        // [TODO] check the correctness of this function
        let mut seeds = Vec::new();

        let hash_space = 10000 as u64;
        let hash_threshold = (fraction * (hash_space as f64)) as u64;

        for kmer in kmer_vec {
            let mut kmer_hash = *kmer;
            if let Some(hash_fn) = self.hash_fn {
                kmer_hash = hash_fn(kmer_hash);
            }

            if kmer_hash % hash_space < hash_threshold {
                seeds.push(*kmer);
            }
        }

        seeds
    }

    fn extract_masked_frac_min_hash(&self, kmer_vec: &[u64], fraction: f64, mask_size: usize) -> Vec<u64> {
        // FracMinHash, but the last 2*`mask_size` bits are masked out
        // when computing the hash value
        let mut seeds = Vec::new();

        let hash_space = 10000 as u64;
        let hash_threshold = (fraction * (hash_space as f64)) as u64;

        let mask = (1 << (2 * mask_size)) - 1;

        for kmer in kmer_vec {
            let mut kmer_hash = *kmer & !mask;
            if let Some(hash_fn) = self.hash_fn {
                kmer_hash = hash_fn(kmer_hash);
            }

            if kmer_hash % hash_space < hash_threshold {
                seeds.push(*kmer);
            }
        }

        seeds
    }

    pub fn sketch_kmer_vec(&self, kmer_vec: &Vec<u64>) -> Vec<u64> {
        let seeds = match self.method {
            SketchingMethod::Minimizer { window_size } => {
                self.extract_minimizers(kmer_vec, window_size)
            }
            SketchingMethod::FracMinHash { fraction } => {
                self.extract_frac_min_hash(kmer_vec, fraction)
            }
            SketchingMethod::MaskedFracMinHash { fraction, mask_size } => {
                self.extract_masked_frac_min_hash(kmer_vec, fraction, mask_size)
            }
        };

        seeds
    }

}


/**
 * @brief Compute the seeds for a given string
 * Adopted from https://github.com/bluenote-1577/sylph/blob/main/src/seeding.rs
 * 
 * @param string The string to compute the seeds for
 * @param kmer_vec The vector to store the seeds
 * @param k The length of the k-mers
 */
pub fn seq_to_canonical_kmer_vec(
    string: &[u8],
    k: usize
) -> Vec<u64> {
    // If the string is shorter than the k-mer length, return
    if string.len() < k {
        return Vec::new();
    }

    let mut kmer_vec: Vec<u64> = Vec::with_capacity(string.len() - k + 1);

    let mut rolling_kmer_f: u64 = 0;
    let mut rolling_kmer_r: u64 = 0;

    let marker_reverse_shift_dist = 2 * (k - 1);
    let marker_mask = (1 << (2*k)) - 1;
    let len = string.len();

    // Init with the first k-1 nucleotides
    for i in 0..k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        let nuc_r = 3 - nuc_f;

        // update the k-mers
        rolling_kmer_f = ((rolling_kmer_f << 2) & marker_mask) | nuc_f;
        rolling_kmer_r = (rolling_kmer_r >> 2) | (nuc_r << marker_reverse_shift_dist);
    }

    // Iterate through the rest of the string
    for i in k-1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;

        // update the k-mers
        rolling_kmer_f = ((rolling_kmer_f << 2) & marker_mask) | nuc_f;
        rolling_kmer_r = (rolling_kmer_r >> 2) | (nuc_r << marker_reverse_shift_dist);

        // Determine the canonical k-mer
        let canonical_kmer_marker = if rolling_kmer_f < rolling_kmer_r {
            rolling_kmer_f
        } else {
            rolling_kmer_r
        };
        kmer_vec.push(canonical_kmer_marker as u64);
    }

    kmer_vec
}


/**
 * @brief Compute the k-mer vector for a given string
 * @returns two vectors: the k-mer vector and the reverse complement k-mer vector
 * 
 * @param string The string to compute the k-mer vector for
 * @param k The length of the k-mers
 */
pub fn seq_to_bidirectional_kmer_vec(
    string: &[u8],
    k: usize
) -> (Vec<u64>, Vec<u64>) {
    // If the string is shorter than the k-mer length, return
    if string.len() < k {
        return (Vec::new(), Vec::new());
    }

    let mut kmer_vec: Vec<u64> = Vec::with_capacity(string.len() - k + 1);
    let mut kmer_vec_rev: Vec<u64> = Vec::with_capacity(string.len() - k + 1);

    let marker_mask = (1 << (2*k)) - 1;
    let len = string.len();

    // Init with the first k-1 nucleotides
    let mut rolling_kmer_f: u64 = 0;
    let mut rolling_kmer_r: u64 = 0;

    for i in 0..k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        let nuc_r = 3 - nuc_f;

        // update the k-mers
        rolling_kmer_f = ((rolling_kmer_f << 2) & marker_mask) | nuc_f;
        rolling_kmer_r = (rolling_kmer_r >> 2) | (nuc_r << (2 * (k - 1)));
    }

    // Iterate through the rest of the string
    for i in k-1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;

        // update the k-mers
        rolling_kmer_f = ((rolling_kmer_f << 2) & marker_mask) | nuc_f;
        rolling_kmer_r = (rolling_kmer_r >> 2) | (nuc_r << (2 * (k - 1)));

        // push the k-mers to the vector
        kmer_vec.push(rolling_kmer_f);
        kmer_vec_rev.push(rolling_kmer_r);
    }

    // reverse the reverse complement k-mer vector
    kmer_vec_rev.reverse();

    (kmer_vec, kmer_vec_rev)
}

/**
 * @brief Compute the canonical k-mer vector for a given string
 * 
 * @param string The string to compute the k-mer vector for
 * @param k The length of the k-mers
 * @returns A vector of canonical k-mers
 */
pub fn seq_to_unidirectional_kmer_vec(
    string: &[u8],
    k: usize
) -> Vec<u64> {
    // If the string is shorter than the k-mer length, return
    if string.len() < k {
        return Vec::new();
    }

    let mut kmer_vec: Vec<u64> = Vec::with_capacity(string.len() - k + 1);

    let marker_mask = (1 << (2*k)) - 1;
    let len = string.len();

    // Init with the first k-1 nucleotides
    let mut rolling_kmer: u64 = 0;

    for i in 0..k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;

        // update the k-mers
        rolling_kmer = ((rolling_kmer << 2) & marker_mask) | nuc_f;
    }

    // Iterate through the rest of the string
    for i in k-1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;

        // update the k-mers
        rolling_kmer = ((rolling_kmer << 2) & marker_mask) | nuc_f;

        // push the k-mers to the vector
        kmer_vec.push(rolling_kmer);
    }

    kmer_vec
}

/**
 * @brief Convert a k-mer string to a u64 index
 */
pub fn kmer_to_index(kmer: &str) -> u64 {
    // use BYTE_TO_SEQ to convert the k-mer to a u64
    let mut res: u64 = 0;
    for b in kmer.bytes() {
        let nuc = BYTE_TO_SEQ[b as usize] as u64;
        res = (res << 2) | nuc;
    }
    res
}

/**
 * @brief Convert a k-mer index to a string
 * 
 * @param kmer The k-mer index
 * @param k The length of the k-mer
 * @returns The k-mer string
 */
pub fn index_to_kmer(kmer: u64, k: usize) -> String {
    // use SEQ_TO_BYTE to convert the k-mer to a string
    let mut kmer_str = String::new();
    for i in (0..k).rev() {
        let nuc = kmer >> (i * 2) & 3;
        kmer_str.push(SEQ_TO_BYTE[nuc as usize] as char);
    }
    kmer_str
}

/**
 * @brief Print the k-mer out
 */
pub fn print_kmer(kmer: u64, k: usize) {
    // use SEQ_TO_BYTE to convert the k-mer to a string
    let kmer_str = index_to_kmer(kmer, k);
    print!("{}", kmer_str);
}


/**
 * @brief Print the k-mer vector
 * 
 * @param kmer_vec The k-mer vector to print
 * @param k The length of the k-mers
 */
pub fn print_kmer_vec(kmer_vec: &Vec<u64>, k: usize) {
    let kmer_str_vec: Vec<String> = kmer_vec.iter().map(|&kmer| {
        // use SEQ_TO_BYTE to convert the k-mer to a string
        index_to_kmer(kmer, k)
    }).collect();
    println!("{:?}", kmer_str_vec);
}




