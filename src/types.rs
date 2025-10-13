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


#[derive(Hash, PartialEq, Eq, Debug, Clone, Copy)]
pub enum EditOperation {
    /* SUBSTITUTION */
    AtoC,
    AtoG,
    AtoT,
    CtoA,
    CtoG,
    CtoT,
    GtoA,
    GtoC,
    GtoT,
    TtoA,
    TtoC,
    TtoG,

    /* INSERTION */
    InsertA,
    InsertC,
    InsertG,
    InsertT,

    /* DELETION */
    DeleteA,
    DeleteC,
    DeleteG,
    DeleteT,

    AMBIGUOUS, // when multiple operations can lead to the same neighbor
}

// 2-D array to map (from, to) -> EditOperation
pub const BASES_TO_SUBSTITUTION: [[Option<EditOperation>; 4]; 4] = {
    let mut arr = [[None; 4]; 4];

    arr[0][1] = Some(EditOperation::AtoC);
    arr[0][2] = Some(EditOperation::AtoG);
    arr[0][3] = Some(EditOperation::AtoT);

    arr[1][0] = Some(EditOperation::CtoA);
    arr[1][2] = Some(EditOperation::CtoG);
    arr[1][3] = Some(EditOperation::CtoT);

    arr[2][0] = Some(EditOperation::GtoA);
    arr[2][1] = Some(EditOperation::GtoC);
    arr[2][3] = Some(EditOperation::GtoT);

    arr[3][0] = Some(EditOperation::TtoA);
    arr[3][1] = Some(EditOperation::TtoC);
    arr[3][2] = Some(EditOperation::TtoG);

    arr
};

pub const BASES_TO_INSERTION: [Option<EditOperation>; 4] = [
    Some(EditOperation::InsertA),
    Some(EditOperation::InsertC),
    Some(EditOperation::InsertG),
    Some(EditOperation::InsertT),
];

pub const BASES_TO_DELETION: [Option<EditOperation>; 4] = [
    Some(EditOperation::DeleteA),
    Some(EditOperation::DeleteC),
    Some(EditOperation::DeleteG),
    Some(EditOperation::DeleteT),
];


pub struct KVmerStats {
    pub k: u8,
    pub v: u8,

    pub keys: Vec<u64>,
    pub consensus_values: Vec<u64>,

    pub consensus_counts: Vec<u32>,
    pub error_counts: Vec<u32>,
    pub total_counts: Vec<u32>,
}
