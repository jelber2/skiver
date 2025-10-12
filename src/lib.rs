pub mod kvmer;
pub mod seeding;
pub mod types;
pub mod profile;
pub mod analyze;
pub mod cmdline;

#[cfg(target_arch = "x86_64")]
pub mod avx2_seeding;


#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_kmer_neighbors() {
        let value = 0b11100100; // AAAA

        let value_length = 4;
        let key_length = 5; // arbitrary

        let kvmer = kvmer::KVmerSet::new(key_length, value_length);
        let neighbors = kvmer._get_neighbors(value);
        println!("Original value: {}", kvmer.kmer_to_string(value));
        kvmer.show_neighbors(value);
        let expected_neighbors = vec![
            0b01, // AAAAC
            0b10, // AAAAG
            0b11, // AAAAT
        ];
        for neighbor in &expected_neighbors {
            assert!(neighbors.contains_key(&neighbor));
        }
        assert_eq!(neighbors.len(), expected_neighbors.len());


    }
}