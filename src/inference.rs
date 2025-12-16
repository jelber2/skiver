use std::path::Display;

use crate::types::*;
use crate::constants::*;

use scirs2_optimize::least_squares::{robust_least_squares, HuberLoss};
use scirs2_core::ndarray::{array, Array1};


pub struct ErrorSpectrum {
    pub substitution_rate: [[(f32, (f32, f32)); 4]; 4],
    pub insertion_rate: [(f32, (f32, f32)); 4],
    pub deletion_rate: [(f32, (f32, f32)); 4],

    pub estimated_alpha: f32,
    pub estimated_beta: f32,
    pub error_rate_mean: (f32, (f32, f32)),
    pub error_rate_std: (f32, (f32, f32)),
}

pub struct ErrorAnalyzer {
    pub outlier_threshold: f32,
}

impl ErrorAnalyzer {
    pub fn new(outlier_threshold: f32) -> Self {
        Self {
            outlier_threshold,
        }
    }


    /*
    ========================
    Util functions for estimating the mean and variance of error rates
    ========================
    */

    /**
     * Perform linear regression with model y = k * x
     * return the slope k
     */
    fn linear_regression_no_intercept<T>(&self, x: &Vec<T>, y: &Vec<T>) -> f32 {
        let n = x.len() as f32;

        if n <= 1. {
            return 0.;
        }

        let sum_xy = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi * yi).sum();
        let sum_x2 = x.iter().map(|&xi| xi * xi).sum();
        let k = sum_xy as f32 / sum_x2 as f32;

        k
    }

    /**
     * Perform linear regression with model y = k * x + b
     * return the slope k and intercept b
     */
    fn linear_regression<T>(&self, x: &Vec<T>, y: &Vec<T>) -> (f32, f32) {
        let n = x.len() as f32;

        if n <= 2. {
            return (0., 0.);
        }
        let sum_x: f32 = x.iter().sum() as f32;
        let sum_y: f32 = y.iter().sum() as f32;
        let sum_xy: f32 = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi * yi).sum() as f32;
        let sum_x2: f32 = x.iter().map(|&xi| xi * xi).sum() as f32;

        let k = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
        let b = (sum_y - k * sum_x) / n;

        (k, b)
    }

    fn linear_regression_no_intercept_heteroskedastic<T>(&self, x: &Vec<T>, y: &Vec<T>) -> f32 {
        let n: f32 = x.len() as f32;

        if n <= 1. {
            return 0.;
        }

        let ratio = y.iter()
            .zip(x.iter())
            .filter_map(|(&yi, &xi)| if xi != 0 { Some(yi as f32 / xi as f32) } else { None })
            .collect::<Vec<f32>>();

        let k: f32 = ratio.iter().copied().sum::<f32>() / ratio.len() as f32;

        k
    }


fn linear_regression_no_intercept(x: &Vec<u32>, y: &Vec<u32>) -> (f32, f32) {
    // perform linear regression with model y = k * x with no bias
    // find k and its confidence interval (lower, upper).
    let n = x.len() as f32;

    if n <= 1. {
        return (0., 0.);
    }
    //let sum_x: f32 = x.iter().map(|&xi| xi as f32).sum();
    let sum_xy: f32 = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi as f32 * yi as f32).sum();
    let sum_x2: f32 = x.iter().map(|&xi| (xi * xi) as f32).sum();

    // for variance
    let sum_y2: f32 = y.iter().map(|&yi| (yi * yi) as f32).sum();
    let variance: f32 = (sum_y2 - sum_xy * sum_xy / sum_x2) / (n - 1.);
    let std: f32 = (variance / sum_x2).sqrt();


    let k = sum_xy / sum_x2;
    let range = 1.96 * std;

    (k, range)
}

fn linear_regression(x: &Vec<f32>, y: &Vec<f32>) -> (f32, f32, f32, f32) {
    // perform linear regression with model y = k * x + b
    // Return (k, k_confidence_interval, b, b_confidence_interval)
    let n = x.len() as f32;

    if n <= 2. {
        return (0., 0., 0., 0.);
    }
    let sum_x: f32 = x.iter().map(|&xi| xi as f32).sum();
    let sum_y: f32 = y.iter().map(|&yi| yi as f32).sum();
    let sum_xy: f32 = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi as f32 * yi as f32).sum();
    let sum_x2: f32 = x.iter().map(|&xi| (xi * xi) as f32).sum();

    // for variance
    let sum_y2: f32 = y.iter().map(|&yi| (yi * yi) as f32).sum();
    let variance: f32 = (sum_y2 - (sum_y * sum_y) / n - (sum_xy - sum_x * sum_y / n).powi(2) / (sum_x2 - sum_x * sum_x / n)) / (n - 2.);
    let std: f32 = (variance / (sum_x2 - sum_x * sum_x / n)).sqrt();
    let k = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
    let range = 1.96 * std;

    let b = (sum_y - k * sum_x) / n;
    let b_std: f32 = (variance * sum_x2 / (n * sum_x2 - sum_x * sum_x)).sqrt();
    let b_range = 1.96 * b_std;

    (k, range, b, b_range)
}

fn trimmed_mean(x: &Vec<f32>, trim_fraction: f32) -> f32 {
    let mut x = x.clone();
    x.sort_by(f32::total_cmp);
    let trim_count = (x.len() as f32 * trim_fraction).round() as usize;
    let trimmed = &x[trim_count..x.len()-trim_count];
    trimmed.iter().sum::<f32>() / trimmed.len() as f32
}

fn linear_regression_no_intercept_heteroskedastic(x: &Vec<u32>, y: &Vec<u32>) -> (f32, f32) {
    let n: f32 = x.len() as f32;

    if n <= 1. {
        return (0., 0.);
    }

    let ratio = y.iter()
        .zip(x.iter())
        .filter_map(|(&yi, &xi)| if xi != 0 { Some(yi as f32 / xi as f32) } else { None })
        .collect::<Vec<f32>>();

    let k: f32 = ratio.iter().copied().sum::<f32>() / ratio.len() as f32;
    //let k: f32 = trimmed_mean(&ratio, 0.1);

    // calculate variance
    let sse: f32 = ratio.iter().map(|r| (r - k) * (r - k)).sum();
    let std: f32 = (sse / (n - 1.)).sqrt();
    let range = 1.96 * std / n.sqrt();

    (k, range)
}

fn median_of_means(x: &Vec<u32>, y: &Vec<u32>, num_means: usize) -> f32 {
    let n: usize = x.len();

    if n == 0 {
        return 0.;
    }

    let chunk_size = (n + num_means - 1) / num_means; // ceiling division
    let mut means: Vec<f32> = Vec::new();

    for chunk in x.chunks(chunk_size).zip(y.chunks(chunk_size)) {
        let (x_chunk, y_chunk) = chunk;
        let sum_x: f32 = x_chunk.iter().map(|&xi| xi as f32).sum();
        let sum_y: f32 = y_chunk.iter().map(|&yi| yi as f32).sum();

        if sum_x != 0. {
            means.push(sum_y / sum_x);
        }
    }

    means.sort_by(f32::total_cmp);
    let mid = means.len() / 2;
    if means.len() % 2 == 0 {
        (means[mid - 1] + means[mid]) / 2.
    } else {
        means[mid]
    }

}

fn sum_mean(x: &Vec<u32>, y: &Vec<u32>) -> (f32, f32) {
    let n: f32 = x.len() as f32;

    if n == 0. {
        return (0., 0.);
    }

    let sum_x: f32 = x.iter().map(|&xi| xi as f32).sum();
    let sum_y: f32 = y.iter().map(|&yi| yi as f32).sum();

    (sum_y / sum_x, 0.)
}

fn error_type_rate(stats: &KVmerStats, error: EditOperation) -> (f32, f32) {
    let mut x: Vec<u32> = Vec::new();
    let mut y: Vec<u32> = Vec::new();

    for (i, error_map) in stats.error_counts.iter().enumerate() {
        let count = error_map.get(&error).cloned().unwrap_or(0);
        x.push(stats.consensus_counts[i] + stats.neighbor_counts[i] - count);
        y.push(count);
    }
    println!("Error type {:?}: total consensus = {}, total errors = {}", error, x.iter().sum::<u32>(), y.iter().sum::<u32>());

    sum_mean(&x, &y)
}


fn infer_p(stats: &KVmerStats, v: u8) -> f32 {
    let x = &stats.consensus_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION) as usize];
    let y = &stats.error_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION) as usize];

    //let alpha = 1.;
    // add pseudocounts
    //let x: Vec<u32> = x.iter().map(|&xi| xi + alpha as u32).collect();
    //let y: Vec<u32> = y.iter().map(|&yi| yi + alpha as u32).collect();
    //println!("Number of data = {}", x.len());

    //linear_regression_no_intercept_heteroskedastic(x, y).0
    //linear_regression_no_intercept(&x, &y).0

    let y_sum = y.iter().sum::<u32>() as f32;
    let x_sum = x.iter().sum::<u32>() as f32;
    //y.iter().sum::<u32>() as f32 / x.iter().sum::<u32>() as f32
    y_sum / (x_sum + y_sum)
    //median_of_means(x, y, 10) // / v as f32
}

fn infer_hazard_ratio(stats: &KVmerStats, v: u8) -> f32 {
    if v - MIN_VALUE_FOR_ERROR_ESTIMATION == 0 {
        let x = &stats.total_counts;
        let y = &stats.consensus_up_to_v_counts[0];

        
        y.iter().sum::<u32>() as f32 / x.iter().sum::<u32>() as f32
        //let (slope, _) = linear_regression_no_intercept(x, y);
        //slope
    } else {
        let x = &stats.consensus_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION - 1) as usize];
        let y = &stats.consensus_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION) as usize];

        y.iter().sum::<u32>() as f32 / x.iter().sum::<u32>() as f32
        //let (slope, _) = linear_regression_no_intercept(x, y);
        //slope
    }  
}

fn residual_hazard_ratio_beta_distribution(params: &[f32], data: &[f32]) -> Array1<f32> {
    let n = data.len() / 2;
    let x = &data[0..n];
    let y = &data[n..];

    let mut res = Array1::zeros(n);
    for i in 0..n {
        res[i] = y[i] - (1 - params[0] / (params[0] + params[1] + x[i]));
    }
    res
}

fn fit_hazard_ratio_beta_distribution(hazard_ratios: &Vec<f32>, k: usize) -> (f32, f32) {
    let n = hazard_ratios.len();
    let mut data: Vec<f32> = Vec::with_capacity(n * 2);
    for &hr in hazard_ratios.iter() {
        data.push(hr);
    }
    for &hr in hazard_ratios.iter() {
        data.push(1. - hr);
    }

    let initial_params = array![1.0, 1.0];
    let result = robust_least_squares(
        &residual_hazard_ratio_beta_distribution,
        &initial_params,
        HuberLoss::new(1.0),
        None::<fn(&[f32], &[f32]) -> scirs2_core::ndarray::Array2<f32>>, 
        &data, None
    );

    (result.x[0], result.x[1])

}



fn total_error_rate(stats: &KVmerStats) -> (f32, f32, f32, f32) {
    let mut x: Vec<f32> = Vec::new();
    let mut y: Vec<f32> = Vec::new();

    /*
    for v in MIN_VALUE_FOR_ERROR_ESTIMATION..=stats.v {
        y.push(infer_p(stats, v as u8));
        println!("{},", y.last().unwrap());
        x.push(v as f32);
    }
    */
    for v in (MIN_VALUE_FOR_ERROR_ESTIMATION)..=stats.v {
        y.push(infer_hazard_ratio(stats, v as u8));
        println!("{},", y.last().unwrap());
        x.push(v as f32);
    }

    linear_regression(&x, &y)
}

fn find_hazard_ratio_outliers(stats: &KVmerStats) -> Vec<bool> {

    let res = [true; stats.consensus_counts.len()].to_vec();
    let x: &Vec<f32>;
    let y: &Vec<f32>;
    
    for v in (MIN_VALUE_FOR_ERROR_ESTIMATION)..=stats.v {
        if v - MIN_VALUE_FOR_ERROR_ESTIMATION == 0 {
            x = &stats.total_counts;
            y = &stats.consensus_up_to_v_counts[0];
        } else {
            x = &stats.consensus_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION - 1) as usize];
            y = &stats.consensus_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION) as usize];
        }
        let ratio = x.iter().zip(y.iter())
            .filter_map(|(&xi, &yi)| if xi != 0 { Some(yi as f32 / xi as f32) } else { 1. })
            .collect::<Vec<f32>>();

        // sort the ratios and exclude the ratios that = 1., and find the IQR
        let mut sorted_ratio = ratio.clone();
        sorted_ratio.sort_by(f32::total_cmp);
        // exclude ratios that are exactly 1.
        let filtered_ratio: Vec<f32> = sorted_ratio.into_iter().filter(|&r| r < 1.).collect();
        let n = filtered_ratio.len();
        if n == 0 {
            continue;
        }
        let q1 = filtered_ratio[n / 4];
        let q3 = filtered_ratio[(3 * n) / 4];
        let iqr = q3 - q1;

        // [TODO] Change the threshold to be a parameter
        let lower_bound = q1 - 3.0 * iqr;
        for (i, &r) in ratio.iter().enumerate() {
            if r < lower_bound {
                // mark as outlier
                res[i] = false;
            }
        }
    }
    
    res
}

pub fn error_profile(stats: &KVmerStats, bidirectional: bool) -> ErrorSpectrum {
    


    let mut substitution_rate: [[(f32, f32); 4]; 4] = Default::default();
    for i in 0..4 {
        for j in 0..4 {
            if i != j {
                let op = BASES_TO_SUBSTITUTION[i][j].unwrap();
                substitution_rate[i][j] = error_type_rate(stats, op);
            }
        }
    }

    let mut deletion_rate: [(f32, f32); 4] = Default::default();
    for i in 0..4 {
        let op = BASES_TO_DELETION[i].unwrap();
        deletion_rate[i] = error_type_rate(stats, op);
    }

    let mut insertion_rate: [(f32, f32); 4] = Default::default();
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
                println!("  {} -> {}: {:.5} ± {:.3}", i, j, k / v as f32 * 100., range / v as f32 * 100.);
            }
        }
    }

    println!("Deletion Rate:");
    for i in 0..4 {
        let (k, range) = spectrum.deletion_rate[i];
        println!("  {}: {:.5} ± {:.3}", i, k / v as f32 * 100., range / v as f32 * 100.);
    }

    println!("Insertion Rate:");
    for i in 0..4 {
        let (k, range) = spectrum.insertion_rate[i];
        println!("  {}: {:.5} ± {:.3}", i, k / v as f32 * 100., range / v as f32 * 100.);
    }

    let (k, k_range, b, b_range) = spectrum.total_error_rate;
    println!("Slope: {:.6}, Intercept: {:.6}", k, b);
}