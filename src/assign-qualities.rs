//! Assign support-derived Phred-33 qualities using an existing Skiver sketch.
//!
//! Integration: this module is declared in `src/lib.rs` via `#[path]` and
//! invoked as a subcommand from `src/main.rs`.
//!
//! Scores are heuristic agreement scores, not calibrated error probabilities.
//! Tests below have not been executed.

use anyhow::{bail, ensure, Context, Result};
use clap::Args;
use crate::kvmer::KVmerSet;
use crate::seeding::fmh_seeds_masked;
use needletail::parse_fastx_file;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Write};
use std::path::PathBuf;

#[derive(Args, Debug)]
#[clap(
    about = "Assign support-derived Phred-33 qualities using an existing Skiver sketch.",
    long_about = "Assign support-derived Phred-33 qualities using an existing Skiver sketch.\n\
                  Output defaults to stdout; '-' also selects stdout.\n\
                  Existing output files are never overwritten.\n\
                  Only FASTA input is accepted.\n\
                  Scores quantify sketch agreement, not calibrated sequencing accuracy."
)]
pub struct AssignQualitiesArgs {
    #[clap(short, long, help_heading = "INPUT", help = "Input Skiver sketch file.")]
    pub sketch: PathBuf,

    #[clap(short, long, help_heading = "INPUT", help = "Input FASTA file (gzip optional).")]
    pub fasta: PathBuf,

    #[clap(short, long, default_value = "-", help_heading = "OUTPUT", help = "Output FASTQ file path, or '-' for stdout.")]
    pub output: String,

    #[clap(long, default_value_t = 2, help_heading = "ALGORITHM", help = "Minimum number of observations required for a key to contribute to per-base quality estimates.")]
    pub min_support: u64,

    #[clap(long, default_value_t = 40, value_parser = clap::value_parser!(u8).range(0..=93), help_heading = "ALGORITHM", help = "Maximum Phred quality to assign (0..=93).")]
    pub max_q: u8,

    #[clap(long, default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=93), help_heading = "ALGORITHM", help = "Fallback Phred quality used when no sketch evidence covers a base (must not exceed --max-q).")]
    pub fallback_q: u8,
}

impl AssignQualitiesArgs {
    fn validate(&self) -> Result<()> {
        ensure!(self.min_support > 0, "--min-support must be positive");
        ensure!(self.max_q <= 93, "--max-q must be between 0 and 93");
        ensure!(
            self.fallback_q <= self.max_q,
            "--fallback-q must not exceed --max-q"
        );
        Ok(())
    }

    fn output_path(&self) -> Option<PathBuf> {
        if self.output == "-" {
            None
        } else {
            Some(PathBuf::from(&self.output))
        }
    }
}

#[derive(Debug)]
struct Profile {
    total: u64,
    counts: Vec<[u64; 4]>,
}

type Profiles = HashMap<u64, Profile>;

fn base_code(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn encoded_base(value: u64, offset: usize, v: usize) -> usize {
    ((value >> (2 * (v - 1 - offset))) & 3) as usize
}

fn build_profiles(sketch: KVmerSet, v: usize) -> Result<Profiles> {
    let mut profiles = HashMap::with_capacity(sketch.key_value_qual_map.len());
    for (key, values) in sketch.key_value_qual_map {
        let mut profile = Profile {
            total: 0,
            counts: vec![[0; 4]; v],
        };
        for (value, observations) in values {
            let n = observations.len() as u64;
            profile.total = profile.total.checked_add(n)
                .context("Observation count overflow")?;
            for j in 0..v {
                let b = encoded_base(value, j, v);
                profile.counts[j][b] = profile.counts[j][b]
                    .checked_add(n).context("Base-support count overflow")?;
            }
        }
        profiles.insert(key, profile);
    }
    Ok(profiles)
}

fn error_probability(total: u64, agreeing: u64) -> f64 {
    // Beta(1,1)-smoothed binary disagreement fraction.
    ((total - agreeing) as f64 + 1.0) / (total as f64 + 2.0)
}

fn quality_byte(error: f64, max_q: u8) -> u8 {
    let q = (-10.0 * error.log10()).round().clamp(0.0, max_q as f64);
    q as u8 + 33
}

fn value_position(start: usize, offset: usize, forward: bool) -> Option<usize> {
    if forward {
        start.checked_add(offset)
    } else {
        start.checked_sub(offset + 1)
    }
}

fn assign(
    sequence: &[u8],
    profiles: &Profiles,
    k: usize,
    v: usize,
    options: &AssignQualitiesArgs,
) -> Result<(Vec<u8>, usize)> {
    let mut error_sums = vec![0.0; sequence.len()];
    let mut context_counts = vec![0u64; sequence.len()];
    let mut keys = Vec::new();
    let mut values = Vec::new();
    let mut infos = Vec::new();

    // Split into clean runs: no context may cross an ambiguous base.
    let mut cursor = 0;
    while cursor < sequence.len() {
        while cursor < sequence.len() && base_code(sequence[cursor]).is_none() {
            cursor += 1;
        }
        let begin = cursor;
        while cursor < sequence.len() && base_code(sequence[cursor]).is_some() {
            cursor += 1;
        }
        if cursor - begin < k + v {
            continue;
        }
        let clean: Vec<u8> = sequence[begin..cursor]
            .iter().map(|b| b.to_ascii_uppercase()).collect();
        ensure!(
            clean.len() <= u32::MAX as usize,
            "ACGT run exceeds Skiver's u32 coordinate range"
        );

        keys.clear();
        values.clear();
        infos.clear();

        // c=1 queries all contexts selected by the existing extractor.
        // Membership in the loaded sketch determines available evidence.
        // Both orientations are queried, even for a forward-only sketch:
        // this permits reverse-oriented target sequences to match evidence.
        fmh_seeds_masked(
            &clean, &mut keys, &mut values, &mut infos, 1, k, v, true,
        );
        ensure!(
            keys.len() == infos.len() && keys.len() == values.len(),
            "Marker extractor returned inconsistent vector lengths"
        );

        for (&key, info) in keys.iter().zip(infos.iter()) {
            let profile = match profiles.get(&key) {
                Some(p) if p.total >= options.min_support => p,
                _ => continue,
            };
            for j in 0..v {
                let local = value_position(
                    info.start_index as usize, j, info.is_forward,
                ).context("Invalid value coordinate from marker extractor")?;
                ensure!(local < clean.len(), "Value coordinate outside sequence");
                let mut b = base_code(clean[local]).unwrap();
                if !info.is_forward {
                    b = 3 - b;
                }
                let position = begin + local;
                error_sums[position] += error_probability(
                    profile.total, profile.counts[j][b],
                );
                context_counts[position] += 1;
            }
        }
    }

    let supported = context_counts.iter().filter(|&&n| n > 0).count();
    let qualities = error_sums.iter().zip(context_counts.iter())
        .map(|(&sum, &count)| {
            if count == 0 {
                options.fallback_q + 33
            } else {
                quality_byte(sum / count as f64, options.max_q)
            }
        }).collect();
    Ok((qualities, supported))
}

pub fn run(args: AssignQualitiesArgs) -> Result<()> {
    args.validate()?;
    let output = args.output_path();

    let reader = BufReader::new(
        File::open(&args.sketch).context("Cannot open sketch")?
    );
    let sketch: KVmerSet = bincode::deserialize_from(reader)
        .context("Cannot decode sketch; use the matching Skiver source version")?;

    let k = sketch.key_size as usize;
    let v = sketch.value_size as usize;
    // The attached scalar extractor shifts by 2*k and 2*v; 32 would shift by 64.
    ensure!(
        (1..32).contains(&k) && (1..32).contains(&v),
        "This implementation requires 1 <= k,v <= 31"
    );
    ensure!(sketch.kv_size as usize == k + v, "Inconsistent sketch dimensions");
    let profiles = build_profiles(sketch, v)?;
    let mut fasta = parse_fastx_file(&args.fasta)
        .context("Cannot open FASTA input")?;

    let sink: Box<dyn Write> = match &output {
        Some(path) => Box::new(
            OpenOptions::new().write(true).create_new(true).open(path)
                .context("Cannot create output (existing files are protected)")?
        ),
        None => Box::new(io::stdout()),
    };
    let mut writer = BufWriter::new(sink);
    let (mut records, mut bases, mut supported) = (0u64, 0u64, 0u64);

    while let Some(record) = fasta.next() {
        let record = record.context("Invalid FASTA record")?;
        if record.qual().is_some() {
            bail!("FASTQ input encountered: this command accepts FASTA only");
        }
        let sequence = record.seq();
        ensure!(!sequence.is_empty(), "Empty FASTA record is not supported");
        ensure!(
            sequence.iter().all(|b| b.is_ascii_graphic()),
            "Sequence contains non-printable or non-ASCII characters"
        );
        let (qualities, covered) = assign(&sequence, &profiles, k, v, &args)?;
        writer.write_all(b"@")?;
        writer.write_all(record.id())?;
        writer.write_all(b"\n")?;
        writer.write_all(&sequence)?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(&qualities)?;
        writer.write_all(b"\n")?;
        records += 1;
        bases += sequence.len() as u64;
        supported += covered as u64;
    }
    writer.flush()?;
    eprintln!("Processed {records} records; {supported}/{bases} bases supported");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn options() -> AssignQualitiesArgs {
        AssignQualitiesArgs {
            sketch: PathBuf::new(),
            fasta: PathBuf::new(),
            output: String::from("-"),
            min_support: 1,
            max_q: 40,
            fallback_q: 0,
        }
    }

    #[test]
    fn phred_encoding_and_smoothing() {
        let opts = options();
        assert_eq!(quality_byte(0.01, opts.max_q), b'5'); // Q20
        assert_eq!(quality_byte(1.0, opts.max_q), b'!');
        assert_eq!(quality_byte(1e-10, opts.max_q), b'I'); // capped Q40
        assert!((error_probability(98, 98) - 0.01).abs() < 1e-12);
    }

    #[test]
    fn coordinates_and_encoding() {
        assert_eq!(value_position(5, 2, true), Some(7));
        assert_eq!(value_position(5, 2, false), Some(2));
        assert_eq!(value_position(0, 0, false), None);
        assert_eq!(encoded_base(0b00011011, 2, 4), 2); // ACGT -> G
    }

    #[test]
    fn forward_reverse_and_ambiguous_mapping() {
        let mut profiles = Profiles::new();
        // Key AA, value C supported 98 times.
        profiles.insert(0, Profile {
            total: 98, counts: vec![[0, 98, 0, 0]],
        });
        let opts = options();
        assert_eq!(assign(b"AAC", &profiles, 2, 1, &opts).unwrap().0, b"!!5");
        // RC(AAC) = GTT: the supported base is now the first base.
        assert_eq!(assign(b"GTT", &profiles, 2, 1, &opts).unwrap().0, b"5!!");
        assert_eq!(assign(b"AANC", &profiles, 2, 1, &opts).unwrap().0, b"!!!!");
        assert_eq!(assign(b"aac", &profiles, 2, 1, &opts).unwrap().0, b"!!5");
    }

    #[test]
    fn unsupported_and_short_sequences() {
        let opts = options();
        assert_eq!(
            assign(b"ACN", &Profiles::new(), 2, 1, &opts).unwrap(),
            (b"!!!".to_vec(), 0)
        );
        assert_eq!(
            assign(b"A", &Profiles::new(), 2, 1, &opts).unwrap().0,
            b"!"
        );
    }
}
