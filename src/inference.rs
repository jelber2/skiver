use crate::types::*;
use crate::constants::*;


pub struct ErrorSpectrum {
    pub substitution_rate: [[(f64, f64); 4]; 4],
    pub insertion_rate: [(f64, f64); 4],
    pub deletion_rate: [(f64, f64); 4],

    pub total_error_rate: (f64, f64, f64, f64),
}

fn linear_regression_no_intercept(x: &Vec<u32>, y: &Vec<u32>) -> (f64, f64) {
    // perform linear regression with model y = k * x with no bias
    // find k and its confidence interval (lower, upper).
    let n = x.len() as f64;

    if n <= 1. {
        return (0., 0.);
    }
    //let sum_x: f64 = x.iter().map(|&xi| xi as f64).sum();
    let sum_xy: f64 = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi as f64 * yi as f64).sum();
    let sum_x2: f64 = x.iter().map(|&xi| (xi * xi) as f64).sum();

    // for variance
    let sum_y2: f64 = y.iter().map(|&yi| (yi * yi) as f64).sum();
    let variance: f64 = (sum_y2 - sum_xy * sum_xy / sum_x2) / (n - 1.);
    let std: f64 = (variance / sum_x2).sqrt();


    let k = sum_xy / sum_x2;
    let range = 1.96 * std;

    (k, range)
}

fn linear_regression(x: &Vec<f64>, y: &Vec<f64>) -> (f64, f64, f64, f64) {
    // perform linear regression with model y = k * x + b
    // Return (k, k_confidence_interval, b, b_confidence_interval)
    let n = x.len() as f64;

    if n <= 2. {
        return (0., 0., 0., 0.);
    }
    let sum_x: f64 = x.iter().map(|&xi| xi as f64).sum();
    let sum_y: f64 = y.iter().map(|&yi| yi as f64).sum();
    let sum_xy: f64 = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi as f64 * yi as f64).sum();
    let sum_x2: f64 = x.iter().map(|&xi| (xi * xi) as f64).sum();

    // for variance
    let sum_y2: f64 = y.iter().map(|&yi| (yi * yi) as f64).sum();
    let variance: f64 = (sum_y2 - (sum_y * sum_y) / n - (sum_xy - sum_x * sum_y / n).powi(2) / (sum_x2 - sum_x * sum_x / n)) / (n - 2.);
    let std: f64 = (variance / (sum_x2 - sum_x * sum_x / n)).sqrt();
    let k = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
    let range = 1.96 * std;

    let b = (sum_y - k * sum_x) / n;
    let b_std: f64 = (variance * sum_x2 / (n * sum_x2 - sum_x * sum_x)).sqrt();
    let b_range = 1.96 * b_std;

    (k, range, b, b_range)
}

fn linear_regression_no_intercept_heteroskedastic(x: &Vec<u32>, y: &Vec<u32>) -> (f64, f64) {
    let n: f64 = x.len() as f64;

    if n <= 1. {
        return (0., 0.);
    }

    let ratio = y.iter().zip(x.iter()).map(|(&yi, &xi)| (yi as f64 / xi as f64)).collect::<Vec<f64>>();

    let k: f64 = ratio.iter().copied().sum::<f64>() / ratio.len() as f64;

    // calculate variance
    let sse: f64 = ratio.iter().map(|r| (r - k) * (r - k)).sum();
    let std: f64 = (sse / (n - 1.)).sqrt();
    let range = 1.96 * std / n.sqrt();

    (k, range)
}


fn error_type_rate(stats: &KVmerStats, error: EditOperation) -> (f64, f64) {
    let mut x: Vec<u32> = Vec::new();
    let mut y: Vec<u32> = Vec::new();

    for (i, error_map) in stats.error_counts.iter().enumerate() {
        let count = error_map.get(&error).cloned().unwrap_or(0);
        x.push(stats.consensus_counts[i] + stats.neighbor_counts[i] - count);
        y.push(count);
    }

    linear_regression_no_intercept_heteroskedastic(&x, &y)
}


fn infer_p_over_v(stats: &KVmerStats, v: u8) -> f64 {
    let x = &stats.consensus_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION) as usize];
    let y = &stats.error_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION) as usize];
    println!("Number of data = {}", x.len());

    linear_regression_no_intercept_heteroskedastic(x, y).0 / v as f64
}



fn total_error_rate(stats: &KVmerStats) -> (f64, f64, f64, f64) {
    let mut x: Vec<f64> = Vec::new();
    let mut y: Vec<f64> = Vec::new();

    for v in MIN_VALUE_FOR_ERROR_ESTIMATION..=stats.v {
        y.push(infer_p_over_v(stats, v as u8));
        println!("v = {}, p/v = {}", v, y.last().unwrap());
        x.push(stats.k as f64 + v as f64);
    }

    linear_regression(&x, &y)
}

fn filter(stats: &KVmerStats) -> Vec<bool> {
    // for each error type, filter out the data points with extremely high error rates
    let mut mask: Vec<bool> = vec![true; stats.consensus_counts.len()];


    mask
}


pub fn error_profile(stats: &KVmerStats, bidirectional: bool) -> ErrorSpectrum {
    


    let mut substitution_rate: [[(f64, f64); 4]; 4] = Default::default();
    for i in 0..4 {
        for j in 0..4 {
            if i != j {
                let op = BASES_TO_SUBSTITUTION[i][j].unwrap();
                substitution_rate[i][j] = error_type_rate(stats, op);
            }
        }
    }

    let mut deletion_rate: [(f64, f64); 4] = Default::default();
    for i in 0..4 {
        let op = BASES_TO_DELETION[i].unwrap();
        deletion_rate[i] = error_type_rate(stats, op);
    }

    let mut insertion_rate: [(f64, f64); 4] = Default::default();
    for i in 0..4 {
        let op = BASES_TO_INSERTION[i].unwrap();
        insertion_rate[i] = error_type_rate(stats, op);
    }

    // find the ambiguous error rates
    let (ambiguous_k, ambiguous_range) = error_type_rate(stats, EditOperation::AMBIGUOUS);

    ErrorSpectrum {
        substitution_rate,
        insertion_rate,
        deletion_rate,
        total_error_rate: total_error_rate(stats),
    }
}

pub fn output_error_spectrum(spectrum: &ErrorSpectrum, v: u8) {
    // Output the error spectrum in a human-readable format
    println!("Substitution Rate:");
    for i in 0..4 {
        for j in 0..4 {
            if i != j {
                let (k, range) = spectrum.substitution_rate[i][j];
                println!("  {} -> {}: {:.3} ± {:.3}", i, j, k / v as f64 * 100., range / v as f64 * 100.);
            }
        }
    }

    println!("Deletion Rate:");
    for i in 0..4 {
        let (k, range) = spectrum.deletion_rate[i];
        println!("  {}: {:.3} ± {:.3}", i, k / v as f64 * 100., range / v as f64 * 100.);
    }

    println!("Insertion Rate:");
    for i in 0..4 {
        let (k, range) = spectrum.insertion_rate[i];
        println!("  {}: {:.3} ± {:.3}", i, k / v as f64 * 100., range / v as f64 * 100.);
    }

    let (k, k_range, b, b_range) = spectrum.total_error_rate;
    println!("Total Error Rate: {:.3}, variance: {:.3}", b, (-k).sqrt());
}