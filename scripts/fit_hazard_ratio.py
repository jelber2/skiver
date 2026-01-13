import numpy as np

def fit_weibull_distribution_linear(v_values, p00_stats):

    # fit 1 - a * (v ** b) to p00_stats
    # Or, fit a linear regression to ln(1 - p00_stats) vs. ln(v)

    #y = np.array([np.log(p) if 0 < p < 1 else 1e-9 for p in p00_stats])
    p00_stats = [min(max(p, 1e-4), 1 - 1e-4) for p in p00_stats]

    # complementary log-log
    y = np.log(- np.log(1 - np.array(p00_stats)))
    x = np.log(np.array(v_values))

    # Use ridge regression to fit the line
    from sklearn.linear_model import Lasso, Ridge, RANSACRegressor, HuberRegressor

    #model = Ridge(alpha=1)
    #model = Lasso(alpha=0.001)
    model = HuberRegressor(alpha=10)
    x_reshaped = x.reshape(-1, 1)
    model.fit(x_reshaped, y)
    #b = model.coef_[0]
    #a = model.intercept_
    #print(model.estimator_.coef_, model.estimator_.intercept_)

    beta = model.coef_[0] + 1
    lam = np.exp(model.intercept_) / beta

    #print(f"Fitted Weibull distribution parameters: lambda={lam}, k={beta}")
    print("({}, {}),".format(lam, beta))

    # Calculate the p-value of the coefficient

    # Residuals
    residuals = y - model.predict(x_reshaped)
    residual_sum_of_squares = np.sum(residuals**2)
    degrees_of_freedom = len(y) - 2  # Number of data points minus number of parameters (a and b)
    variance = residual_sum_of_squares / degrees_of_freedom

    # Standard error of the coefficient
    x_mean = np.mean(x)
    standard_error_b = np.sqrt(variance / np.sum((x - x_mean)**2))

    # t-statistic for the coefficient
    #t_stat_b = b / standard_error_b

    # Two-tailed p-value
    #p_value_b = 2 * (1 - stats.t.cdf(np.abs(t_stat_b), df=degrees_of_freedom))
    #print(f"p-value for coefficient b: {p_value_b}")

    # a = log(-log(q))
    #q = np.exp(a)
    #print(f"Estimated q from a: {q}")

    # Plot the fit
    
    import matplotlib.pyplot as plt
    plt.scatter(x, y, label="Data")
    plt.plot(x, model.predict(x_reshaped), color='red', label="Fit")
    plt.xlabel("ln(v)")
    plt.ylabel("ln(1 - p00_stats)")
    plt.legend()
    plt.show()

    #print(f"Fitted Weibull distribution parameters: a={np.exp(a)}, b={b}")
    return lam, beta  # a, b

def fit_weibull_distribution_curve_fit(self, v_values, p00_stats):
    from scipy.optimize import curve_fit

    def weibull_func(v, a, b):
        return a * (v ** b)
    
    # Initial guess for a and b
    initial_guess = [0.0, -1.0]
    params, covariance = curve_fit(weibull_func, v_values, p00_stats, p0=initial_guess)
    a, b = params

    

    # Plot the fit
    import matplotlib.pyplot as plt
    plt.scatter(v_values, p00_stats, label="Data")
    plt.plot(v_values, weibull_func(v_values, *params), color='red', label="Fit")
    plt.xlabel("v")
    plt.ylabel("p00_stats")
    plt.legend()
    plt.show()
    print(f"Fitted Weibull distribution parameters: a={a}, b={b}")

def neg_loglike(params, t, r, d):
    # params are unconstrained: alpha=log(lambda), beta=log(k)
    alpha, beta = params
    lam = np.exp(alpha)
    k = np.exp(beta)

    # Δ_t = (t+1)^k - t^k  (t can be numpy array)
    Dt = (t + 1.0) ** k - (t * 1.0) ** k

    # h(t) = 1 - exp(-lam * Dt) computed stably
    x = -lam * Dt
    h = -np.expm1(x)          # 1 - exp(x) with x negative => hazard
    one_minus_h = np.exp(x)   # exp(-lam*Dt)

    # avoid log(0) if any h is exactly 0 or 1 numerically
    eps = 1e-15
    h = np.clip(h, eps, 1 - eps)
    one_minus_h = np.clip(one_minus_h, eps, 1 - eps)

    ll = np.sum(d * np.log(h) + (r - d) * np.log(one_minus_h))
    return -ll

def fit_weibull_distribution_mle(t, r, d):
    from scipy.optimize import minimize

    # initial guess
    x0 = np.array([np.log(1e-4), np.log(1.0)])

    res = minimize(neg_loglike, x0, args=(t, r, d), method="L-BFGS-B")
    alpha_hat, beta_hat = res.x
    lam_hat, k_hat = np.exp(alpha_hat), np.exp(beta_hat)
    print(f"Fitted Weibull distribution parameters (MLE): lambda={lam_hat}, k={k_hat}")

    # Plot the fit
    import matplotlib.pyplot as plt
    v_values = np.array(t)
    hazard_rates = np.array(d) / np.array(r)
    plt.scatter(v_values, hazard_rates, label="Data")
    t_fit = np.linspace(0, max(v_values), 100)
    Dt_fit = 1 - np.exp(-lam_hat * ((t_fit + 1) ** k_hat - t_fit ** k_hat))
    plt.plot(t_fit, Dt_fit, color='red', label="Fit")
    plt.xlabel("v")
    plt.ylabel("Hazard Rate")
    plt.legend()
    plt.show()


    return lam_hat, k_hat

if __name__ == "__main__":
    import pandas as pd
    from scipy import stats

    #hazard_ratio_df = pd.read_csv("hazard_ratio.csv")
    #hazard_ratio_df = pd.read_csv("../kv-mer-test/output/coverage_dependence_homogeneous_l5/Ecoli_K12_MG1655_depth_100_id_90_exp_1_hazard_ratio.csv")
    #hazard_ratio_df = pd.read_csv("../kv-mer-test/output/coverage_dependence_homogeneous_l5/Ecoli_K12_MG1655_depth_100_id_90_exp_1_hazard_ratio.csv")

    #hazard_ratio_df = pd.read_csv("../kv-mer-test/output/human/HG002_bi_kvmer_hazard_ratio.csv")

    file_list = [
        "../kv-mer-test/output/zymo/ERR3152366_bi_kvmer_hazard_ratio.csv",
        "../kv-mer-test/output/zymo/SRR13128014_bi_kvmer_hazard_ratio.csv",
        "../kv-mer-test/output/zymo/ERR2935851_bi_kvmer_hazard_ratio.csv",
        "../kv-mer-test/output/zymo/SRR7498042_bi_kvmer_hazard_ratio.csv"
    ]

    file_list = [
        "../kv-mer-test/output/human/HG002_R941_bi_kvmer_hazard_ratio.csv",
    ]
    #hazard_ratio_df = pd.read_csv("../kv-mer-test/output/human/HG002_hifi_bi_kvmer_hazard_ratio.csv")
    #hazard_ratio_df = pd.read_csv("/home/ubuntu/kv-mer-test/output/coverage_dependence/Ecoli_K12_MG1655_depth_10_id_100_exp_1_ref_hazard_ratio.csv")
    #hazard_ratio_df = pd.read_csv("/home/ubuntu/kv-mer-test/output/coverage_dependence/Ecoli_K12_MG1655_depth_10_id_90_exp_4_hazard_ratio.csv")

    for file in file_list:
        #print(f"Processing file: {file}")
        hazard_ratio_df = pd.read_csv(file)
        v_values = hazard_ratio_df["t"].values
        p00_stats = hazard_ratio_df["hazard_ratio"].values
        print(p00_stats)
        num_candidates = hazard_ratio_df["num_candidates"].values
        num_survival = hazard_ratio_df["num_survival"].values
        num_failure = num_candidates - num_survival

        #v_values = v_values[1:]
        #p00_stats = p00_stats[1:]

        fit_weibull_distribution_linear(v_values, p00_stats)
        #fit_weibull_distribution_mle(v_values, num_candidates, num_failure)
