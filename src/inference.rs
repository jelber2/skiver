use std::hash::Hash;
use std::path::Display;
use std::collections::HashMap;


use crate::types::*;
use crate::constants::*;
use crate::cmdline::AnalyzeArgs;

use scirs2_optimize::least_squares::{least_squares, Method};
use scirs2_stats::regression::{linear_regression, ridge_regression, lasso_regression};
use scirs2_core::ndarray::{array, Array1};
use log::info;
use rand::Rng;


pub struct ErrorSpectrum {
    pub estimated_a: (f32, (f32, f32)),
    pub estimated_b: (f32, (f32, f32)),
    
    pub key_coverage: (f32, (f32, f32)),
    pub estimated_coverage: (f32, (f32, f32)),

    pub snp_rate: HashMap<(EditOperation, u8, u8), f32>,

    pub bidirectional: bool,
}


pub enum RatioEstimationMethod {
    Slope,
    LinearFit,
    RatioMean,
    SumRatio,
}

pub struct ErrorAnalyzer {
    /* 
    pub k: u8,

    pub bidirectional: bool,
    pub exclude_outliers: bool,
    pub outlier_threshold: f32,
    pub ratio_method: RatioEstimationMethod,

    // bootstrap parameters
    pub num_experiments: u32,
    pub bootstrap_sample_rate: f32,
    */
    pub args: AnalyzeArgs,

    pub ratio_method: RatioEstimationMethod,
}



impl ErrorAnalyzer {
    pub fn new(args: AnalyzeArgs) -> Self {
        let method = match args.estimation_method.as_str() {
            "slope" => RatioEstimationMethod::Slope,
            "linear_fit" => RatioEstimationMethod::LinearFit,
            "ratio_mean" => RatioEstimationMethod::RatioMean,
            "sum_ratio" => RatioEstimationMethod::SumRatio,
            _ => {
                panic!("Unknown ratio estimation method: {}. Supported methods are: slope, linear_fit, ratio_mean, sum_ratio.", args.estimation_method);
            }
        };


        Self {
            args,
            ratio_method: method,
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
    fn slope(x: &Vec<u32>, y: &Vec<u32>, indices: &Vec<usize>) -> f32 {
        let n = indices.len() as f32;

        if n <= 1. {
            return 0.;
        }

        let sum_xy: f32 = indices.iter().map(|&i| x[i] * y[i]).sum::<u32>() as f32;
        let sum_x2: f32 = indices.iter().map(|&i| x[i] * x[i]).sum::<u32>() as f32;
        if sum_x2 == 0.0 {
            return 0.;
        }
        let k = sum_xy / sum_x2;

        k
    }

    /**
     * Perform linear regression with model y = k * x + b
     * return the slope k and intercept b
     */
    fn linear_fit(x: &Vec<u32>, y: &Vec<u32>, indices: &Vec<usize>) -> (f32, f32) {
        let n = indices.len() as f32;

        if n <= 2. {
            return (0., 0.);
        }
        let sum_x: f32 = indices.iter().map(|&i| x[i]).sum::<u32>() as f32;
        let sum_y: f32 = indices.iter().map(|&i| y[i]).sum::<u32>() as f32;
        let sum_xy: f32 = indices.iter().map(|&i| x[i] * y[i]).sum::<u32>() as f32;
        let sum_x2: f32 = indices.iter().map(|&i| x[i] * x[i]).sum::<u32>() as f32;

        let denom = n * sum_x2 - sum_x * sum_x;
        if denom == 0.0 {
            return (0., 0.);
        }

        let k = (n * sum_xy - sum_x * sum_y) / denom;
        let b = (sum_y - k * sum_x) / n;

        (k, b)
    }

    fn linear_fit_f32(x: &Vec<f32>, y: &Vec<f32>) -> (f32, f32) {
        let n = x.len() as f32;

        if n <= 2. {
            return (0., 0.);
        }
        let sum_x: f32 = x.iter().sum::<f32>();
        let sum_y: f32 = y.iter().sum::<f32>();
        let sum_xy: f32 = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi * yi).sum::<f32>();
        let sum_x2: f32 = x.iter().map(|&xi| xi * xi).sum::<f32>();

        let denom = n * sum_x2 - sum_x * sum_x;
        if denom == 0.0 {
            return (0., 0.);
        }

        let k = (n * sum_xy - sum_x * sum_y) / denom;
        let b = (sum_y - k * sum_x) / n;

        (k, b)
    }

    fn ridge_fit_f32(x: &Vec<f32>, y: &Vec<f32>, lambda: f32) -> (f32, f32) {
        let n = x.len();
        // Means
        let mut sum_x = x.iter().sum::<f32>();
        let mut sum_y = y.iter().sum::<f32>();

        let mean_x = sum_x / n as f32;
        let mean_y = sum_y / n as f32;

        // Centered sums
        let mut sxx = 0.0f32;
        let mut sxy = 0.0f32;
        for i in 0..n {
            let dx = x[i] - mean_x;
            let dy = y[i] - mean_y;
            sxx += dx * dx;
            sxy += dx * dy;
        }

        // Ridge on slope only
        let denom = sxx + lambda;
        let k = if denom != 0.0 { sxy / denom } else { 0.0 };
        let b = mean_y - k * mean_x;

        (k, b)
    }

    /**
     * Calculate the mean of the ratios y/x
     */
    fn ratio_mean(x: &Vec<u32>, y: &Vec<u32>, indices: &Vec<usize>) -> f32 {
        let n: f32 = indices.len() as f32;

        if n <= 1. {
            return 0.;
        }

        let ratio = indices.iter()
            .filter_map(|&i| if x[i] != 0 { Some(y[i] as f32 / x[i] as f32) } else { None })
            .collect::<Vec<f32>>();

        if ratio.is_empty() {
            return 0.;
        }

        let k: f32 = ratio.iter().sum::<f32>() / ratio.len() as f32;

        k
    }

    /**
     * Calculate the sum of ratios sum(y) / sum(x)
     */
    fn sum_ratio(x: &Vec<u32>, y: &Vec<u32>, indices: &Vec<usize>) -> f32 {
        let n: f32 = indices.len() as f32;

        if n == 0. {
            return 0.;
        }

        let sum_x: f32 = indices.iter().map(|&i| x[i]).sum::<u32>() as f32;
        if sum_x == 0.0 {
            return 0.;
        }
        let sum_y: f32 = indices.iter().map(|&i| y[i]).sum::<u32>() as f32;
        //println!("sum_y: {}, sum_x: {}", sum_y, sum_x);
        sum_y / sum_x
    }

    /**
     * The function that uses different methods to calculate the ratio
     */
    fn calculate_ratio(&self, x: &Vec<u32>, y: &Vec<u32>, indices: &Vec<usize>) -> f32 {
        match self.ratio_method {
            RatioEstimationMethod::Slope => Self::slope(x, y, indices),
            RatioEstimationMethod::LinearFit => {
                let (k, _) = Self::linear_fit(x, y, indices);
                k
            },
            RatioEstimationMethod::RatioMean => Self::ratio_mean(x, y, indices),
            RatioEstimationMethod::SumRatio => Self::sum_ratio(x, y, indices),
        }
    }

    fn sum_indices(&self, x: &Vec<u32>, indices: &Vec<usize>) -> u32 {
        indices.iter().map(|&i| x[i]).sum()
    }

    fn random_subsample_with_replacement(x: &Vec<usize>, n: usize) -> Vec<usize> {
        let mut rng = rand::thread_rng();
        (0..n)
            .map(|_| x[rng.gen_range(0..x.len())])
            .collect()
    }

    fn random_subsample_without_replacement(x: &Vec<usize>, n: usize) -> Vec<usize> {
        use rand::seq::SliceRandom;
        let mut rng = rand::thread_rng();
        let mut x_clone = x.clone();
        x_clone.shuffle(&mut rng);
        x_clone.truncate(n);
        x_clone
    }


    /*
    ========================
    Util functions for estimating the parameters of beta distribution
    ========================
    */
    fn residual_hazard_ratio_beta_distribution(params: &[f64], data: &[f64]) -> Array1<f64> {
        let n = data.len() / 2;
        let x = &data[0..n];
        let y = &data[n..];

        let mut res = Array1::zeros(n);
        for i in 0..n {
            res[i] = y[i] - (params[0] / (params[0] + params[1] + x[i]));
        }
        res
    }

    fn residual_hazard_ratio_beta_distribution_fixed_kappa(params: &[f64], data: &[f64]) -> Array1<f64> {
        let n = data.len() / 2;
        let x = &data[0..n];
        let y = &data[n..];

        let mut res = Array1::zeros(n);
        for i in 0..n {
            res[i] = y[i] - (params[0] / (MIN_KAPPA as f64 + x[i]));
        }
        res
    }

    fn fit_hazard_ratio_beta_distribution(&self, hazard_ratios: &Vec<f32>) -> (f32, f32) {
        let n = hazard_ratios.len();
        let mut vec_data: Vec<f64> = Vec::with_capacity(n * 2);
        for i in 1..=n {
            vec_data.push(i as f64 + self.args.k as f64);
        }
        for &hr in hazard_ratios.iter() {
            vec_data.push(hr as f64);
        }

        // convert the Vec<f64> to an Array1<f64> as required by robust_least_squares
        let data = Array1::from_vec(vec_data);

        let initial_params = array![1.0f64, 1.0f64];
        let result = least_squares(
            &Self::residual_hazard_ratio_beta_distribution,
            &initial_params,
            Method::LevenbergMarquardt,
            None::<fn(&[f64], &[f64]) -> scirs2_core::ndarray::Array2<f64>>, 
            &data, None
        ).expect("robust_least_squares failed");

        if result.x[1] <= MIN_KAPPA as f64 {
            // refit with fixed kappa
            let initial_params_fixed = array![1.0f64];
            //warn!("Estimated beta parameter is too small ({}), probably due to high error rate and low coverage. Refitting with fixed alpha + beta = {}.", result.x[1], MIN_KAPPA);
            //warn!("Consider increasing the coverage or using bidirectional kmers to improve the estimation.");
            let result_fixed = least_squares(
                &Self::residual_hazard_ratio_beta_distribution_fixed_kappa,
                &initial_params_fixed,
                Method::LevenbergMarquardt,
                None::<fn(&[f64], &[f64]) -> scirs2_core::ndarray::Array2<f64>>, 
                &data, None
            ).expect("robust_least_squares failed");

            return (result_fixed.x[0] as f32, MIN_KAPPA - result_fixed.x[0] as f32);
        }

        //println!("data: {:?}", data);

        (result.x[0] as f32, result.x[1] as f32)
    }

    /*
    ========================
    Util functions for estimating the parameters of Weibull distribution
    ========================
    */
    fn fit_hazard_ratio_weibull_distribution(&self, hazard_ratios: &Vec<f32>) -> (f32, f32) {
        // Fit hazard ratio = a * (i + k)^b, or log(hazard ratio) = log(a) + b * log(i + k)
        let x = hazard_ratios.iter().enumerate().
            map(|(i, _)| (i as f32 + self.args.k as f32).ln())
            .collect::<Vec<f32>>();
        let y = hazard_ratios.iter()
            .map(|&hr| if hr > 0.0 { hr.ln() } else { 0.0 })
            .collect::<Vec<f32>>();
        //let (b, log_a) = Self::linear_fit_f32(&x, &y);
        let (b, log_a) = Self::ridge_fit_f32(&x, &y, 1.);
        
        
        let a = log_a.exp();

        (a, b)
    }

    

    /**
     * Identify outliers based on hazard ratios across different v values,
     * return the indices of inliers.
     */
    pub fn find_hazard_ratio_outliers(&self, stats: &KVmerStats) -> Vec<usize> {

        let mut res = vec![true; stats.consensus_counts.len()];
        let mut x: &Vec<u32>;
        let mut y: &Vec<u32>;
        // [TODO] Return also the number of outliers for each v
        //let mut num_outliers = [0; stats.consensus_counts.len()];
        
        for v in 1..=(stats.v - self.args.ignore_last_hazard_ratios as u8) {
            if v - 1 == 0 {
                x = &stats.total_counts;
                y = &stats.consensus_up_to_v_counts[0];
            } else {
                x = &stats.consensus_up_to_v_counts[(v - 1 - 1) as usize];
                y = &stats.consensus_up_to_v_counts[(v - 1) as usize];
            }
            let ratio = x.iter().zip(y.iter())
                .map(|(&xi, &yi)| if xi != 0 { yi as f32 / xi as f32 } else { 1. })
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

            let lower_bound = q1 - self.args.outlier_threshold * iqr;
            for (i, &r) in ratio.iter().enumerate() {
                if r < lower_bound {
                    // mark as outlier
                    res[i] = false;
                }
            }
        }

        let indices: Vec<usize> = res.iter().enumerate()
                                .filter_map(|(i, &is_inlier)| if is_inlier { Some(i) } else { None })
                                .collect();
        
        info!("Identified {} inliers out of {} data points based on hazard ratios ({}%).", indices.len(), res.len(), (indices.len() as f32 / res.len() as f32) * 100.0);

        indices
    }


    /**
     * Estimate the error rate and 5-95% confidence interval using bootstrap
     * for all the error types
     */
    pub fn estimate_error_rate(&self, stats: &KVmerStats, indices: &Vec<usize>) -> HashMap<(EditOperation, u8, u8), f32> {
        // initialize the error count arrays
        let mut error_counts: HashMap<(EditOperation, u8, u8), u32> = HashMap::new();

        indices.iter().for_each(|&i| {
            for (op, count_map) in stats.error_counts[i].iter() {
                let count = error_counts.entry(*op).or_insert(0);
                *count += *count_map;
            }
        });
        
        // calculate the mean for each error type using the full error_counts vector
        let mut estimates: HashMap<(EditOperation, u8, u8), f32> = HashMap::new();
        let mut total_count: u32 = 0;
        for op in if self.args.bidirectional { ALL_OPERATIONS_CANONICAL.iter() } else { ALL_OPERATIONS.iter() } {
            for prev_base in 0..4 {
                for next_base in 0..4 {
                    total_count += error_counts.get(&(*op, prev_base, next_base)).unwrap_or(&0);
                }
            }
        }

        for op in if self.args.bidirectional { ALL_OPERATIONS_CANONICAL.iter() } else { ALL_OPERATIONS.iter() } {
            for prev_base in 0..4 {
                for next_base in 0..4 {
                    let count = error_counts.get(&(*op, prev_base, next_base)).unwrap_or(&0);
                    let rate = if total_count > 0 {
                        *count as f32 / total_count as f32
                    } else {
                        0.0
                    };
                    estimates.insert((*op, prev_base, next_base), rate);
                }
            }
        }
        
        /*
        // bootstrap to estimate the 5-95% confidence interval
        let mut bootstrap_estimates: HashMap<EditOperation, Vec<f32>> = HashMap::new();
        for op in operations.iter() {
            bootstrap_estimates.insert(*op, Vec::new());
        }

        for _ in 0..self.num_experiments {
            let indices_sample = Self::random_subsample_with_replacement(indices, indices.len() as usize);
            
            let mut error_count: HashMap<EditOperation, u32> = HashMap::new();
            for op in operations.iter() {
                let sum = self.sum_indices(&error_counts[op], &indices_sample);
                error_count.insert(*op, sum);
            }
            // normalize the error counts so that they sum to 1
            let total_count: u32 = error_count.values().sum();
            for op in operations.iter() {
                let rate = if total_count > 0 {
                    error_count[op] as f32 / total_count as f32
                } else {
                    0.0
                };
                bootstrap_estimates.get_mut(op).unwrap().push(rate);
            }
        }

        // calculate the 5-95% confidence interval
        let mut result: Vec<(f32, (f32, f32))> = Vec::new();
        for op in operations.iter() {
            let mut estimates_op = bootstrap_estimates[op].clone();
            estimates_op.sort_by(f32::total_cmp);
            let n = estimates_op.len();
            let mean = estimates[op];
            let lower = estimates_op[(n as f32 * 0.05) as usize];
            let upper = estimates_op[(n as f32 * 0.95) as usize];
            result.push((mean, (lower, upper)));
        }
        */

        estimates
    }

    /**
     * Returns ((lower_a, upper_a), (lower_b, upper_b), hazard_ratio_list)
     */
    pub fn estimate_hazard_ratio_confidence_interval(&self, stats: &KVmerStats, indices: &Vec<usize>) -> ((f32, f32), (f32, f32), Vec<((f32, f32))>) {
        let mut x: &Vec<u32>;
        let mut y: &Vec<u32>;

        // record the estimated a and b
        let mut a_list: Vec<f32> = Vec::new();
        let mut b_list: Vec<f32> = Vec::new();

        // record hazard ratios for each v
        let mut hazard_ratio_list: Vec<Vec<f32>> = Vec::new();
        for v in 1..=(stats.v - self.args.ignore_last_hazard_ratios as u8) {
            hazard_ratio_list.push(Vec::new());
        }

        for _ in 0..self.args.num_experiments {
            let indices_sample = Self::random_subsample_with_replacement(indices, indices.len() as usize);

            let mut hazard_ratios: Vec<f32> = Vec::new();

            for v in 1..=(stats.v - self.args.ignore_last_hazard_ratios as u8) {
                if v - 1 == 0 {
                    x = &stats.total_counts;
                    y = &stats.consensus_up_to_v_counts[0];
                } else {
                    x = &stats.consensus_up_to_v_counts[(v - 1 - 1) as usize];
                    y = &stats.consensus_up_to_v_counts[(v - 1) as usize];
                }

                let h = self.calculate_ratio(x, y, &indices_sample);
                hazard_ratios.push(1. - h);
                hazard_ratio_list[(v - 1) as usize].push(1. - h);
            }
            // estimate the parameters of the beta distribution
            //let (alpha, beta) = self.fit_hazard_ratio_beta_distribution(&hazard_ratios, (indices.len() as f32 * self.bootstrap_sample_rate) as usize);

            let (a, b) = self.fit_hazard_ratio_weibull_distribution(&hazard_ratios);
            a_list.push(a);
            b_list.push(b);
        }

        a_list.sort_by(f32::total_cmp);
        let lower_a = a_list[(self.args.num_experiments as f32 * 0.05) as usize];
        let upper_a = a_list[(self.args.num_experiments as f32 * 0.95) as usize];

        b_list.sort_by(f32::total_cmp);
        let lower_b = b_list[(self.args.num_experiments as f32 * 0.05) as usize];
        let upper_b = b_list[(self.args.num_experiments as f32 * 0.95) as usize];

        let mut hazard_ratio_range_list: Vec<((f32, f32))> = Vec::new();
        for v in 0..hazard_ratio_list.len() {
            hazard_ratio_list[v].sort_by(f32::total_cmp);
            let h_lower = hazard_ratio_list[v][(self.args.num_experiments as f32 * 0.05) as usize];
            let h_upper = hazard_ratio_list[v][(self.args.num_experiments as f32 * 0.95) as usize];
            hazard_ratio_range_list.push((h_lower, h_upper));
        }

        ((lower_a, upper_a), (lower_b, upper_b), hazard_ratio_range_list)
        /* 
        let mut mean = alpha_list.iter().zip(beta_list.iter())
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
        */
    }


    pub fn estimate_hazard_ratio(&self, stats: &KVmerStats, indices: &Vec<usize>) -> (f32, f32, Vec<f32>) {
        let mut x: &Vec<u32>;
        let mut y: &Vec<u32>;

        let mut hazard_ratios: Vec<f32> = Vec::new();

        for v in 1..=(stats.v - self.args.ignore_last_hazard_ratios as u8) {
            if v - 1 == 0 {
                x = &stats.total_counts;
                y = &stats.consensus_up_to_v_counts[0];
            } else {
                x = &stats.consensus_up_to_v_counts[(v - 1 - 1) as usize];
                y = &stats.consensus_up_to_v_counts[(v - 1) as usize];
            }

            /*
            println!("v={}, x_sum={}, y_sum={}", v, 
                    x.iter().sum::<u32>(), // self.sum_indices(x, indices), 
                    y.iter().sum::<u32>(), //self.sum_indices(y, indices));
            );
            */

            let h = self.calculate_ratio(x, y, indices);
            hazard_ratios.push(1. - h);
        }
        /*
        for &h in hazard_ratios.iter() {
            println!("{},", h);
        }
        */
        /* 
        // estimate the parameters of the beta distribution
        let (alpha, beta) = self.fit_hazard_ratio_beta_distribution(&hazard_ratios, indices.len());
        let mean = alpha / (alpha + beta);
        let std = (alpha * beta / (((alpha + beta) * (alpha + beta)) * (alpha + beta + 1.0))).sqrt();

        (mean, std, alpha, beta)
        */

        let (a, b) = self.fit_hazard_ratio_weibull_distribution(&hazard_ratios);
        //println!("Weibull parameters: alpha = {}, beta = {}", a, b);
        (a, b, hazard_ratios)

        
    }

    pub fn key_coverage(&self, stats: &KVmerStats, indices: &Vec<usize>) -> (f32, (f32, f32)) {
        let mut coverages: Vec<u32> = indices.iter().map(|&i| stats.total_counts[i]).collect();
        coverages.sort_unstable();
        let n = coverages.len();
        if n == 0 {
            return (0., (0., 0.));
        }

        let median_coverage = if n % 2 == 0 {
            (coverages[n / 2 - 1] + coverages[n / 2]) as f32 / 2.0
        } else {
            coverages[n / 2] as f32
        };
        let coverage_ci_lower = coverages[(n as f32 * 0.05) as usize] as f32;
        let coverage_ci_upper = coverages[(n as f32 * 0.95) as usize] as f32;
        
        (median_coverage, (coverage_ci_lower, coverage_ci_upper))
    }

    pub fn estimate_true_coverage(&self, estimated_a: f32, estimated_b: f32, key_coverage: (f32, (f32, f32))) -> (f32, (f32, f32)) {
        // estimate survival rate at k
        let mut survival_rate: f32 = 1.0;
        for i in 1..=self.args.k {
            survival_rate *= 1.0 - estimated_a * ((i as f32).powf(estimated_b));
        }

        if survival_rate <= 0.0 {
            return (0., (0., 0.));
        }
        let mut estimated_coverage = key_coverage.0 / survival_rate;
        let mut estimated_coverage_ci_lower = (key_coverage.1).0 / survival_rate;
        let mut estimated_coverage_ci_upper = (key_coverage.1).1 / survival_rate;

        if !self.args.bidirectional {
            estimated_coverage *= 2.0;
            estimated_coverage_ci_lower *= 2.0;
            estimated_coverage_ci_upper *= 2.0;
        }

        (estimated_coverage, (estimated_coverage_ci_lower, estimated_coverage_ci_upper))
    }


    pub fn analyze(&self, stats: &KVmerStats) -> ErrorSpectrum {
        // exclude the hazard ratio outliers
        let indices = if !self.args.use_all {
            self.find_hazard_ratio_outliers(stats)
        } else {
            (0..stats.consensus_counts.len()).collect()
        };

        // estimate SNP rates
        let error_rates = self.estimate_error_rate(stats, &indices);

        // estimate hazard ratio parameters
        let (a, b, hazard_ratio) = self.estimate_hazard_ratio(stats, &indices);
        let (a_ci, b_ci, hazard_ratio_ci) = self.estimate_hazard_ratio_confidence_interval(stats, &indices);

        if let Some(hazard_ratio_output) = &self.args.hazard_ratio {
            use std::fs::File;
            use std::io::{BufWriter, Write};

            let file = File::create(hazard_ratio_output).expect("Could not create hazard ratio output file.");
            let mut writer = BufWriter::new(file);

            writeln!(writer, "v,hazard_ratio,5th_percentile,95th_percentile").expect("Could not write to hazard ratio output file.");
            for v in 0..hazard_ratio.len() {
                writeln!(writer, "{},{:.6},{:.6},{:.6}", 
                    v + 1, 
                    hazard_ratio[v], 
                    hazard_ratio_ci[v].0, 
                    hazard_ratio_ci[v].1
                ).expect("Could not write to hazard ratio output file.");
            }
        }

        // estimate key coverage
        let key_coverage = self.key_coverage(stats, &indices);
        let estimated_coverage = self.estimate_true_coverage(a, b, key_coverage);

        ErrorSpectrum {
            estimated_a: (a, a_ci),
            estimated_b: (b, b_ci),

            key_coverage: key_coverage,
            estimated_coverage: estimated_coverage,

            snp_rate: error_rates,
            bidirectional: self.args.bidirectional,
        }
    }
}


/**
 * Format the spectrum into a line in a csv file
 */
pub fn spectrum_to_str(spectrum: &ErrorSpectrum, bidirectional: bool) -> String {
    let mut result = String::new();

    if bidirectional != spectrum.bidirectional {
        panic!("The bidirectional flag does not match the spectrum data.");
    }

    // hazard ratio parameters a and b
    result.push_str(&format!("{:.6},{:.6}~{:.6},", spectrum.estimated_a.0, (spectrum.estimated_a.1).0, (spectrum.estimated_a.1).1));
    result.push_str(&format!("{:.6},{:.6}~{:.6},", spectrum.estimated_b.0, (spectrum.estimated_b.1).0, (spectrum.estimated_b.1).1));

    // key coverage and estimated true coverage
    result.push_str(&format!("{:.6},{:.6}~{:.6},", spectrum.key_coverage.0, (spectrum.key_coverage.1).0, (spectrum.key_coverage.1).1));
    result.push_str(&format!("{:.6},{:.6}~{:.6},", spectrum.estimated_coverage.0, (spectrum.estimated_coverage.1).0, (spectrum.estimated_coverage.1).1));

    // SNP rates
    for op in if bidirectional { ALL_OPERATIONS_CANONICAL.iter() } else { ALL_OPERATIONS.iter() } {
        for prev_base in 0..4 {
            for next_base in 0..4 {
                let rate = spectrum.snp_rate.get(&(*op, prev_base, next_base)).unwrap_or(&0.0);
                result.push_str(&format!("{:.6},", rate));
            }
        }
    }

    // remove the last comma
    result.pop();

    result
}


pub fn header_str(bidirectional: bool) -> String {
    let mut result = String::new();

    result.push_str("a,a_5-95th_percentile,");
    result.push_str("b,b_5-95th_percentile,");

    result.push_str("key_median_coverage,key_coverage_5-95th_percentile,");
    result.push_str("true_median_coverage,true_coverage_5-95th_percentile,");
    for op in if bidirectional { ALL_OPERATIONS_CANONICAL.iter() } else { ALL_OPERATIONS.iter() } {
        for prev_base in 0..4 {
            for next_base in 0..4 {
                result.push_str(&sbs96_str(&(*op, prev_base, next_base)));
                result.push(',');
            }
        }
    }

    result
}