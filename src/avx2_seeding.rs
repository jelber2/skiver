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

use std::arch::x86_64::*;
use crate::types::*;

/**
 * A fast hash function for 64-bit kmers using AVX2 instructions.
 */
#[inline]
#[target_feature(enable = "avx2")]
pub unsafe fn mm_hash256_masked(kmer: __m256i, mask: i64) -> __m256i { unsafe {
    // mask the kmer so that only the masked bits are used in the hash
    let mask_vec = _mm256_set_epi64x(mask, mask, mask, mask);
    let masked_kmer = _mm256_and_si256(kmer, mask_vec);

    let mut key = masked_kmer;
    let s1 = _mm256_slli_epi64(key, 21);
    key = _mm256_add_epi64(key, s1);
    
    key = _mm256_xor_si256(key, _mm256_cmpeq_epi64(key, key));

    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 24));
    let s2 = _mm256_slli_epi64(key, 3);
    let s3 = _mm256_slli_epi64(key, 8);

    key = _mm256_add_epi64(key, s2);
    key = _mm256_add_epi64(key, s3);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 14));
    let s4 = _mm256_slli_epi64(key, 2);
    let s5 = _mm256_slli_epi64(key, 4);
    key = _mm256_add_epi64(key, s4);
    key = _mm256_add_epi64(key, s5);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 28));

    let s6 = _mm256_slli_epi64(key, 31);
    key = _mm256_add_epi64(key, s6);

    return key;
}}

/**
 * Shift each 64-bit lane in a __m256i by 2*k - 2 bits to the left.
 * This is used to update the reverse kmer in the rolling hash.
 */
#[target_feature(enable = "avx2")]
pub unsafe fn _shift_mm256_by_k(kmer: __m256i, k: usize) -> __m256i { unsafe {
    // shift left by 2*k - 2
    let shifted = match k {
        
        3 => { _mm256_slli_epi64(kmer, 4) }
        4 => { _mm256_slli_epi64(kmer, 6) }
        5 => { _mm256_slli_epi64(kmer, 8) }
        6 => { _mm256_slli_epi64(kmer, 10) }
        7 => { _mm256_slli_epi64(kmer, 12) }
        8 => { _mm256_slli_epi64(kmer, 14) }
        9 => { _mm256_slli_epi64(kmer, 16) }
        10 => { _mm256_slli_epi64(kmer, 18) }
        11 => { _mm256_slli_epi64(kmer, 20) }
        12 => { _mm256_slli_epi64(kmer, 22) }
        13 => { _mm256_slli_epi64(kmer, 24) }
        14 => { _mm256_slli_epi64(kmer, 26) }
        15 => { _mm256_slli_epi64(kmer, 28) }
        16 => { _mm256_slli_epi64(kmer, 30) }
        17 => { _mm256_slli_epi64(kmer, 32) }
        18 => { _mm256_slli_epi64(kmer, 34) }
        19 => { _mm256_slli_epi64(kmer, 36) }
        20 => { _mm256_slli_epi64(kmer, 38) }
        21 => { _mm256_slli_epi64(kmer, 40) }
        22 => { _mm256_slli_epi64(kmer, 42) }
        23 => { _mm256_slli_epi64(kmer, 44) }
        24 => { _mm256_slli_epi64(kmer, 46) }
        25 => { _mm256_slli_epi64(kmer, 48) }
        26 => { _mm256_slli_epi64(kmer, 50) }
        27 => { _mm256_slli_epi64(kmer, 52) }
        28 => { _mm256_slli_epi64(kmer, 54) }
        29 => { _mm256_slli_epi64(kmer, 56) }
        30 => { _mm256_slli_epi64(kmer, 58) }
        31 => { _mm256_slli_epi64(kmer, 60) }
        32 => { _mm256_slli_epi64(kmer, 62) }
        
        //21 => { _mm256_slli_epi64(kmer, 40) }
        //31 => { _mm256_slli_epi64(kmer, 60) }
        _ => { panic!() }
    };
    return shifted;
}}

/**
 * Extract kmers using FracMinHash from a DNA string using AVX2 instructions.
 * The k-mers are extracted from both the forward and reverse strands.
 */
#[target_feature(enable = "avx2")]
pub unsafe fn extract_markers_avx2_masked(string: &[u8], kmer_vec: &mut Vec<u64>, c: usize, k: usize, mask: i64, bidirectional: bool) { unsafe {
    if string.len() <= k {
        return;
    }

    // divide the string into 4 parts for parallel processing
    let len = (string.len() - k + 1) / 4;
    let string1 = &string[0..len + k - 1];
    let string2 = &string[len..2 * len + k - 1];
    let string3 = &string[2 * len..3 * len + k - 1];
    let string4 = &string[3 * len..4 * len + k - 1];

    /*
    let use_40 = if 2 * (k - 1) == 40 {
        true
    } else if 2 * (k - 1) == 60 {
        false
    } else {
        panic!()
    };
    */

    //const TWO_K_MINUS_2_40: i32 = 40;
    //const TWO_K_MINUS_2_60: i32 = 60;

    let mut rolling_kmer_f_marker = _mm256_set_epi64x(0, 0, 0, 0);
    let mut rolling_kmer_r_marker = _mm256_set_epi64x(0, 0, 0, 0);
    let rev_sub = _mm256_set_epi64x(3, 3, 3, 3);

    // Initialize
    for i in 0..k - 1 {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        let r_nucs = _mm256_sub_epi64(rev_sub, f_nucs);

        rolling_kmer_f_marker = _mm256_slli_epi64(rolling_kmer_f_marker, 2);
        rolling_kmer_f_marker = _mm256_or_si256(rolling_kmer_f_marker, f_nucs);
        if bidirectional {
            rolling_kmer_r_marker = _mm256_srli_epi64(rolling_kmer_r_marker, 2);

            /*
            let shift_nuc_r;
            if use_40 {
                shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_40);
            } else {
                shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_60);
            }
            */
            let shift_nuc_r = _shift_mm256_by_k(r_nucs, k);
            rolling_kmer_r_marker = _mm256_or_si256(rolling_kmer_r_marker, shift_nuc_r);
        }
    }

    let marker_mask = (Kmer::MAX >> (std::mem::size_of::<Kmer>() * 8 - 2 * k)) as i64;
    let rev_marker_mask: i64 = !(0 | (3 << 2 * k - 2));
    //    let rev_marker_mask = i64::from_le_bytes(rev_marker_mask.to_le_bytes());
    //    dbg!(u64::MAX / (c as u64));
    //    dbg!((u64::MAX / (c as u64)) as i64);
    let threshold_marker = u64::MAX / c as u64;

    let mm256_marker_mask = _mm256_set_epi64x(marker_mask, marker_mask, marker_mask, marker_mask);
    let mm256_rev_marker_mask = _mm256_set_epi64x(
        rev_marker_mask,
        rev_marker_mask,
        rev_marker_mask,
        rev_marker_mask,
    );

    // iterate over the string
    for i in k - 1..(len + k - 1) {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        let r_nucs = _mm256_sub_epi64(rev_sub, f_nucs);

        // f_marker = ((f_marker << 2) | f_nuc) & marker_mask
        rolling_kmer_f_marker = _mm256_slli_epi64(rolling_kmer_f_marker, 2);
        rolling_kmer_f_marker = _mm256_or_si256(rolling_kmer_f_marker, f_nucs);
        rolling_kmer_f_marker = _mm256_and_si256(rolling_kmer_f_marker, mm256_marker_mask);

        if bidirectional {
            // r_marker = ((r_marker >> 2) | (r_nuc << (2*(k-1)))) & rev_marker_mask
            rolling_kmer_r_marker = _mm256_srli_epi64(rolling_kmer_r_marker, 2);
            /*
            let shift_nuc_r;
            if use_40 {
                shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_40);
            } else {
                shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_60);
            }
            */
            let shift_nuc_r = _shift_mm256_by_k(r_nucs, k);
            rolling_kmer_r_marker = _mm256_and_si256(rolling_kmer_r_marker, mm256_rev_marker_mask);
            rolling_kmer_r_marker = _mm256_or_si256(rolling_kmer_r_marker, shift_nuc_r);
        }

        

        //let compare_marker = _mm256_cmpgt_epi64(rolling_kmer_r_marker, rolling_kmer_f_marker);

        //let canonical_markers_256 =
        //    _mm256_blendv_epi8(rolling_kmer_r_marker, rolling_kmer_f_marker, compare_marker);

        //        dbg!(rolling_kmer_f_marker,rolling_kmer_r_marker);
        //        dbg!(print_string(u64::from_ne_bytes(_mm256_extract_epi64(rolling_kmer_f_marker,1).to_ne_bytes()), 31));
        
        
        let hash_256_f = mm_hash256_masked(rolling_kmer_f_marker, mask);
        let v1 = _mm256_extract_epi64(hash_256_f, 0) as u64;
        let v2 = _mm256_extract_epi64(hash_256_f, 1) as u64;
        let v3 = _mm256_extract_epi64(hash_256_f, 2) as u64;
        let v4 = _mm256_extract_epi64(hash_256_f, 3) as u64;

        let kmer1 = _mm256_extract_epi64(rolling_kmer_f_marker, 0) as u64;
        let kmer2 = _mm256_extract_epi64(rolling_kmer_f_marker, 1) as u64;
        let kmer3 = _mm256_extract_epi64(rolling_kmer_f_marker, 2) as u64;
        let kmer4 = _mm256_extract_epi64(rolling_kmer_f_marker, 3) as u64;
        //        let threshold_256 = _mm256_cmpgt_epi64(cmp_thresh, hash_256);
        //        let m1 = _mm256_extract_epi64(threshold_256, 0);
        //        let m2 = _mm256_extract_epi64(threshold_256, 1);
        //        let m3 = _mm256_extract_epi64(threshold_256, 2);
        //        let m4 = _mm256_extract_epi64(threshold_256, 3);

        if v1 < threshold_marker {
            kmer_vec.push(kmer1 as u64);
        }
        if v2 < threshold_marker {
            kmer_vec.push(kmer2 as u64);
        }
        if v3 < threshold_marker {
            kmer_vec.push(kmer3 as u64);
        }
        if v4 < threshold_marker {
            kmer_vec.push(kmer4 as u64);
        }

        if bidirectional {
            let hash_256_r = mm_hash256_masked(rolling_kmer_r_marker, mask);
            let vr1 = _mm256_extract_epi64(hash_256_r, 0) as u64;
            let vr2 = _mm256_extract_epi64(hash_256_r, 1) as u64;
            let vr3 = _mm256_extract_epi64(hash_256_r, 2) as u64;
            let vr4 = _mm256_extract_epi64(hash_256_r, 3) as u64;

            let rkmer1 = _mm256_extract_epi64(rolling_kmer_r_marker, 0) as u64;
            let rkmer2 = _mm256_extract_epi64(rolling_kmer_r_marker, 1) as u64;
            let rkmer3 = _mm256_extract_epi64(rolling_kmer_r_marker, 2) as u64;
            let rkmer4 = _mm256_extract_epi64(rolling_kmer_r_marker, 3) as u64;

            if vr1 < threshold_marker {
                kmer_vec.push(rkmer1 as u64);
            }
            if vr2 < threshold_marker {
                kmer_vec.push(rkmer2 as u64);
            }
            if vr3 < threshold_marker {
                kmer_vec.push(rkmer3 as u64);
            }
            if vr4 < threshold_marker {
                kmer_vec.push(rkmer4 as u64);
            }
        }
    }
}}