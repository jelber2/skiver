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


pub enum RatioEstimationMethod {
    Slope,
    LinearFit,
    RatioMean,
    SumRatio,
}

pub struct ErrorAnalyzer {
    pub outlier_threshold: f32,
    pub ratio_method: RatioEstimationMethod,

    // bootstrap parameters
    pub num_experiments: u32,
    pub data_per_experiment: u32,
}



impl ErrorAnalyzer {
    pub fn new(outlier_threshold: f32, ratio_method: RatioEstimationMethod, num_experiments: u32, data_per_experiment: u32) -> Self {
        outlier_threshold,
        ratio_method,
        num_experiments,
        data_per_experiment,
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
    fn slope<T>(&self, x: &Vec<T>, y: &Vec<T>, indices: &Vec<u32>) -> f32 {
        let n = x.len() as f32;

        if n <= 1. {
            return 0.;
        }

        let sum_xy = indices.iter().map(|&i| x[i] * y[i]).sum();
        let sum_x2 = indices.iter().map(|&i| x[i] * x[i]).sum();
        let k = sum_xy as f32 / sum_x2 as f32;

        k
    }

    /**
     * Perform linear regression with model y = k * x + b
     * return the slope k and intercept b
     */
    fn linear_fit<T>(&self, x: &Vec<T>, y: &Vec<T>, indices: &Vec<u32>) -> (f32, f32) {
        let n = x.len() as f32;

        if n <= 2. {
            return (0., 0.);
        }
        let sum_x: f32 = indices.iter().map(|&i| x[i]).sum() as f32;
        let sum_y: f32 = indices.iter().map(|&i| y[i]).sum() as f32;
        let sum_xy: f32 = indices.iter().map(|&i| x[i] * y[i]).sum() as f32;
        let sum_x2: f32 = indices.iter().map(|&i| x[i] * x[i]).sum() as f32;

        let k = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
        let b = (sum_y - k * sum_x) / n;

        (k, b)
    }

    /**
     * Calculate the mean of the ratios y/x
     */
    fn ratio_mean<T>(&self, x: &Vec<T>, y: &Vec<T>, indices: &Vec<u32>) -> f32 {
        let n: f32 = x.len() as f32;

        if n <= 1. {
            return 0.;
        }

        let ratio = indices.iter()
            .filter_map(|&i| if x[i] != 0 { Some(y[i] as f32 / x[i] as f32) } else { None })
            .collect::<Vec<f32>>();

        let k: f32 = ratio.iter().sum::<f32>() / ratio.len() as f32;

        k
    }

    /**
     * Calculate the sum of ratios sum(y) / sum(x)
     */
    fn sum_ratio<T>(x: &Vec<T>, y: &Vec<T>, indices: &Vec<u32>) -> f32 {
        let n: f32 = x.len() as f32;

        if n == 0. {
            return 0.;
        }

        let sum_x: f32 = indices.iter().map(|&i| x[i]).sum() as f32;
        let sum_y: f32 = indices.iter().map(|&i| y[i]).sum() as f32;

        sum_y / sum_x
    }

    /**
     * The function that uses different methods to calculate the ratio
     */
    fn calculate_ratio<T>(&self, x: &Vec<T>, y: &Vec<T>, indices: &Vec<u32>) -> f32 {
        match self.ratio_method {
            RatioEstimationMethod::Slope => self.slope(x, y, indices),
            RatioEstimationMethod::LinearFit => {
                let (k, _) = self.linear_fit(x, y, indices);
                k
            },
            RatioEstimationMethod::RatioMean => self.ratio_mean(x, y, indices),
            RatioEstimationMethod::SumRatio => Self::sum_ratio(x, y, indices),
        }
    }
    

    /**
     * Identify outliers based on hazard ratios across different v values,
     * return the indices of inliers.
     */
    pub fn find_hazard_ratio_outliers(stats: &KVmerStats) -> Vec<u32> {

        let mut res = [true; stats.consensus_counts.len()].to_vec();
        let mut x: &Vec<f32>;
        let mut y: &Vec<f32>;
        // [TODO] Return also the number of outliers for each v
        //let mut num_outliers = [0; stats.consensus_counts.len()];
        
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
            let lower_bound = q1 - self.outlier_threshold * iqr;
            for (i, &r) in ratio.iter().enumerate() {
                if r < lower_bound {
                    // mark as outlier
                    res[i] = false;
                }
            }
        }

        let indices = res.iter().enumerate()
                                .filter_map(|(i, &is_inlier)| if is_inlier { Some(i as u32) } else { None })
                                .collect();
        
        info!("Identified {} inliers out of {} data points based on hazard ratios ({}%).", indices.len(), res.len(), (indices.len() as f32 / res.len() as f32) * 100.0);

        indices
    }


    pub fn estimate_error_rate(&self, stats: &KVmerStats, error_type: EditOperation, indices: &Vec<u32>) -> (f32, (f32, f32)) {
        let mut x: Vec<u32> = Vec::new();
        let mut y: Vec<u32> = Vec::new();

        for (i, error_map) in stats.error_counts.iter().enumerate() {
            let count = error_map.get(&error_type).cloned().unwrap_or(0);
            x.push(stats.consensus_counts[i] + stats.neighbor_counts[i] - count);
            y.push(count);
        }

        let k = self.calculate_ratio(&x, &y, indices);

        // bootstrap to estimate the 5-95% confidence interval
        let mut bootstrap_estimates: Vec<f32> = Vec::with_capacity(self.num_experiments as usize);
        for _ in 0..self.num_experiments {
            let mut indices_sample: Vec<u32> = Vec::with_capacity(self.data_per_experiment as usize);
            for _ in 0..self.data_per_experiment {
                let idx = indices[rand::random::<usize>() % indices.len()];
                indices_sample.push(idx);
            }
            let k_sample = self.calculate_ratio(&x, &y, &indices_sample);
            bootstrap_estimates.push(k_sample);
        }

        // calculate the 5-95% confidence interval
        bootstrap_estimates.sort_by(f32::total_cmp);
        let lower = bootstrap_estimates[(self.num_experiments as f32 * 0.05) as usize];
        let upper = bootstrap_estimates[(self.num_experiments as f32 * 0.95) as usize];

        (k, (lower, upper))
    }

    pub fn estimate_hazard_ratio_confidence_interval(&self, stats: &KVmerStats, indices: &Vec<u32>) -> ((f32, f32), (f32, f32)) {
        let mut x: &Vec<u32>;
        let mut y: &Vec<u32>;

        let mut alpha_list: Vec<f32> = Vec::new();
        let mut beta_list: Vec<f32> = Vec::new();

        for _ in 0..self.num_experiments {
            let mut indices_sample: Vec<u32> = Vec::with_capacity(self.data_per_experiment as usize);
            for _ in 0..self.data_per_experiment {
                let idx = indices[rand::random::<usize>() % indices.len()];
                indices_sample.push(idx);
            }

            let mut hazard_ratios: Vec<f32> = Vec::new();

            for v in (MIN_VALUE_FOR_ERROR_ESTIMATION)..=stats.v {
                if v - MIN_VALUE_FOR_ERROR_ESTIMATION == 0 {
                    x = &stats.total_counts;
                    y = &stats.consensus_up_to_v_counts[0];
                } else {
                    x = &stats.consensus_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION - 1) as usize];
                    y = &stats.consensus_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION) as usize];
                }

                let h = self.calculate_ratio(x, y, &indices_sample);
                hazard_ratios.push(h);
            }
            // estimate the parameters of the beta distribution
            let (alpha, beta) = fit_hazard_ratio_beta_distribution(&hazard_ratios, self.data_per_experiment as usize);
            alpha_list.push(alpha);
            beta_list.push(beta);
        }

        let mean = alpha_list.iter().zip(beta_list.iter())
            .map(|(&a, &b)| a / (a + b))
            .collect::<Vec<f32>>();
        mean.sort_by(f32::total_cmp);
        let lower_mean = mean[(self.num_experiments as f32 * 0.05) as usize];
        let upper_mean = mean[(self.num_experiments as f32 * 0.95) as usize];

        let std = alpha_list.iter().zip(beta_list.iter())
            .map(|(&a, &b)| ((a * b) / (((a + b) * (a + b)) * (a + b + 1.0))).sqrt())
            .collect::<Vec<f32>>();
        let mut std_sorted = std.clone();
        std_sorted.sort_by(f32::total_cmp);
        let lower_std = std_sorted[(self.num_experiments as f32 * 0.05) as usize];
        let upper_std = std_sorted[(self.num_experiments as f32 * 0.95) as usize];

        ((lower_mean, upper_mean), (lower_std, upper_std))
    }


    pub fn estimate_hazard_ratio(&self, stats: &KVmerStats, indices: &Vec<u32>) -> (f32, f32, f32, f32) {
        let mut x: &Vec<u32>;
        let mut y: &Vec<u32>;

        let mut hazard_ratios: Vec<f32> = Vec::new();

        for v in (MIN_VALUE_FOR_ERROR_ESTIMATION)..=stats.v {
            if v - MIN_VALUE_FOR_ERROR_ESTIMATION == 0 {
                x = &stats.total_counts;
                y = &stats.consensus_up_to_v_counts[0];
            } else {
                x = &stats.consensus_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION - 1) as usize];
                y = &stats.consensus_up_to_v_counts[(v - MIN_VALUE_FOR_ERROR_ESTIMATION) as usize];
            }

            let h = self.calculate_ratio(x, y, indices);
            hazard_ratios.push(h);
        }
        // estimate the parameters of the beta distribution
        (alpha, beta) = fit_hazard_ratio_beta_distribution(&hazard_ratios, indices.len())
        mean = alpha / (alpha + beta);
        std = (alpha * beta / (((alpha + beta) * (alpha + beta)) * (alpha + beta + 1.0))).sqrt();

        (mean, std, alpha, beta)
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