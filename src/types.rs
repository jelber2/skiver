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

use std::collections::HashMap;

pub type Kmer = u64;

/**
 * A lookup table to convert a byte to a 2-bit sequence.
 * Adopted from https://github.com/bluenote-1577/sylph/blob/main/src/types.rs
 * 
 * A/a -> 0
 * C/c -> 1
 * G/g -> 2
 * T/t -> 3, U/u -> 3
 */
pub const BYTE_TO_SEQ: [u8; 256] = {
    let mut arr = [0u8; 256];

    arr[b'A' as usize] = 0;
    arr[b'C' as usize] = 1;
    arr[b'G' as usize] = 2;
    arr[b'T' as usize] = 3;
    arr[b'U' as usize] = 3;

    arr[b'a' as usize] = 0;
    arr[b'c' as usize] = 1;
    arr[b'g' as usize] = 2;
    arr[b't' as usize] = 3;
    arr[b'u' as usize] = 3;

    arr
};

pub const SEQ_TO_BYTE: [u8; 4] = [b'A', b'C', b'G', b'T'];


/**
 * Definitions for edit operations.
 */
#[derive(Hash, PartialEq, Eq, Debug, Clone, Copy)]
pub enum SubstitutionOperations {
    AC,
    AG,
    AT,
    CA,
    CG,
    CT,
    GA,
    GC,
    GT,
    TA,
    TC,
    TG,
}

#[derive(Hash, PartialEq, Eq, Debug, Clone, Copy)]
pub enum InsertionOperations {
    A,
    C,
    G,
    T,
}

#[derive(Hash, PartialEq, Eq, Debug, Clone, Copy)]
pub enum DeletionOperations {
    A,
    C,
    G,
    T,
}

#[derive(Hash, PartialEq, Eq, Debug, Clone, Copy)]
pub enum EditOperation {
    /* SUBSTITUTION */
    SUBSTITUTION(SubstitutionOperations),

    /* INSERTION */
    INSERTION(InsertionOperations),

    /* DELETION */
    DELETION(DeletionOperations),
    AMBIGUOUS, // when multiple operations can lead to the same neighbor
}

// 2-D array to map (from, to) -> EditOperation
pub const BASES_TO_SUBSTITUTION: [[Option<EditOperation>; 4]; 4] = {
    let mut arr = [[None; 4]; 4];

    arr[0][1] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::AC));
    arr[0][2] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::AG));
    arr[0][3] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::AT));

    arr[1][0] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::CA));
    arr[1][2] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::CG));
    arr[1][3] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::CT));

    arr[2][0] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::GA));
    arr[2][1] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::GC));
    arr[2][3] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::GT));

    arr[3][0] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::TA));
    arr[3][1] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::TC));
    arr[3][2] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::TG));

    arr
};

// Use in the case we include both forward and reverse complements of the reads
pub const BASES_TO_SUBSTITUTION_CANONICAL: [[Option<EditOperation>; 4]; 4] = {
    let mut arr = [[None; 4]; 4];

    arr[0][1] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::AC));
    arr[0][2] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::AG));
    arr[0][3] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::AT));

    arr[1][0] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::AC));
    arr[1][2] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::CG));
    arr[1][3] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::CT));

    arr[2][0] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::AG));
    arr[2][1] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::CG));
    arr[2][3] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::GT));

    arr[3][0] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::AT));
    arr[3][1] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::CT));
    arr[3][2] = Some(EditOperation::SUBSTITUTION(SubstitutionOperations::GT));

    arr
};

pub const BASES_TO_INSERTION: [Option<EditOperation>; 4] = [
    Some(EditOperation::INSERTION(InsertionOperations::A)),
    Some(EditOperation::INSERTION(InsertionOperations::C)),
    Some(EditOperation::INSERTION(InsertionOperations::G)),
    Some(EditOperation::INSERTION(InsertionOperations::T)),
];

pub const BASES_TO_INSERTION_CANONICAL: [Option<EditOperation>; 4] = [
    Some(EditOperation::INSERTION(InsertionOperations::A)),
    Some(EditOperation::INSERTION(InsertionOperations::C)),
    Some(EditOperation::INSERTION(InsertionOperations::C)),
    Some(EditOperation::INSERTION(InsertionOperations::A)),
];

pub const BASES_TO_DELETION: [Option<EditOperation>; 4] = [
    Some(EditOperation::DELETION(DeletionOperations::A)),
    Some(EditOperation::DELETION(DeletionOperations::C)),
    Some(EditOperation::DELETION(DeletionOperations::G)),
    Some(EditOperation::DELETION(DeletionOperations::T)),
];

pub const BASES_TO_DELETION_CANONICAL: [Option<EditOperation>; 4] = [
    Some(EditOperation::DELETION(DeletionOperations::A)),
    Some(EditOperation::DELETION(DeletionOperations::C)),
    Some(EditOperation::DELETION(DeletionOperations::C)),
    Some(EditOperation::DELETION(DeletionOperations::A)),
];

pub const ALL_OPERATIONS: [EditOperation; 21] = [
    EditOperation::SUBSTITUTION(SubstitutionOperations::AC),
    EditOperation::SUBSTITUTION(SubstitutionOperations::AG),
    EditOperation::SUBSTITUTION(SubstitutionOperations::AT),

    EditOperation::SUBSTITUTION(SubstitutionOperations::GA),
    EditOperation::SUBSTITUTION(SubstitutionOperations::GC),
    EditOperation::SUBSTITUTION(SubstitutionOperations::GT),

    EditOperation::SUBSTITUTION(SubstitutionOperations::CA),
    EditOperation::SUBSTITUTION(SubstitutionOperations::CG),
    EditOperation::SUBSTITUTION(SubstitutionOperations::CT),

    EditOperation::SUBSTITUTION(SubstitutionOperations::TA),
    EditOperation::SUBSTITUTION(SubstitutionOperations::TC),
    EditOperation::SUBSTITUTION(SubstitutionOperations::TG),

    EditOperation::INSERTION(InsertionOperations::A),
    EditOperation::INSERTION(InsertionOperations::C),
    EditOperation::INSERTION(InsertionOperations::G),
    EditOperation::INSERTION(InsertionOperations::T),


    EditOperation::DELETION(DeletionOperations::A),
    EditOperation::DELETION(DeletionOperations::C),
    EditOperation::DELETION(DeletionOperations::G),
    EditOperation::DELETION(DeletionOperations::T),

    EditOperation::AMBIGUOUS,
];


pub const ALL_OPERATIONS_CANONICAL: [EditOperation; 11] = [
    EditOperation::SUBSTITUTION(SubstitutionOperations::AC),
    EditOperation::SUBSTITUTION(SubstitutionOperations::AG),
    EditOperation::SUBSTITUTION(SubstitutionOperations::AT),

    EditOperation::SUBSTITUTION(SubstitutionOperations::GC),
    EditOperation::SUBSTITUTION(SubstitutionOperations::GT),

    EditOperation::SUBSTITUTION(SubstitutionOperations::CT),

    EditOperation::INSERTION(InsertionOperations::A),
    EditOperation::INSERTION(InsertionOperations::C),

    EditOperation::DELETION(DeletionOperations::A),
    EditOperation::DELETION(DeletionOperations::C),

    EditOperation::AMBIGUOUS,
];

/**
 * kv-mer statistics for downstream analysis.
 */
pub struct KVmerStats {
    pub k: u8,
    pub v: u8,

    pub keys: Vec<u64>,
    pub consensus_values: Vec<u64>,

    pub consensus_counts: Vec<u32>,
    pub total_counts: Vec<u32>,
    pub neighbor_counts: Vec<u32>,
    pub error_counts: Vec<HashMap<EditOperation, u32>>,

    pub consensus_up_to_v_counts: Vec<Vec<u32>>,
    pub error_up_to_v_counts: Vec<Vec<u32>>,
}


#[derive(Clone)]
pub struct SequenceInfo {
    pub seq: Vec<u8>,
}
