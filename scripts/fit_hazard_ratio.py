import numpy as np

def fit_weibull_distribution_linear(v_values, p00_stats):

    # fit 1 - a * (v ** b) to p00_stats
    # Or, fit a linear regression to ln(1 - p00_stats) vs. ln(v)

    #y = np.array([np.log(p) if 0 < p < 1 else 1e-9 for p in p00_stats])

    # complementary log-log
    y = np.log(- np.log(1 - np.array(p00_stats)))
    x = np.log(np.array(v_values))

    # Use ridge regression to fit the line
    from sklearn.linear_model import Lasso, Ridge

    model = Ridge(alpha=1)
    x_reshaped = x.reshape(-1, 1)
    model.fit(x_reshaped, y)
    b = model.coef_[0]
    a = model.intercept_

    # a = log(-log(q))
    q = np.exp(a)
    print(f"Estimated q from a: {q}")

    # Plot the fit
    import matplotlib.pyplot as plt
    plt.scatter(x, y, label="Data")
    plt.plot(x, model.predict(x_reshaped), color='red', label="Fit")
    plt.xlabel("ln(v)")
    plt.ylabel("ln(1 - p00_stats)")
    plt.legend()
    plt.show()

    print(f"Fitted Weibull distribution parameters: a={np.exp(a)}, b={b}")
    return np.exp(a), b  # a, b

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

if __name__ == "__main__":
    import pandas as pd

    #hazard_ratio_df = pd.read_csv("HG002_hazard_ratio.csv")
    #hazard_ratio_df = pd.read_csv("../kv-mer-test/output/coverage_dependence_homogeneous_l5/Ecoli_K12_MG1655_depth_100_id_90_exp_1_hazard_ratio.csv")
    #hazard_ratio_df = pd.read_csv("../kv-mer-test/output/coverage_dependence_homogeneous_l5/Ecoli_K12_MG1655_depth_100_id_90_exp_1_hazard_ratio.csv")

    #hazard_ratio_df = pd.read_csv("../kv-mer-test/output/human/HG002_bi_kvmer_hazard_ratio.csv")
    hazard_ratio_df = pd.read_csv("../kv-mer-test/output/zymo/ERR3152366_bi_kvmer_hazard_ratio.csv")
    hazard_ratio_df = pd.read_csv("../kv-mer-test/output/zymo/SRR7498042_bi_kvmer_hazard_ratio.csv")

    v_values = hazard_ratio_df["t"].values
    p00_stats = hazard_ratio_df["hazard_ratio"].values

    #v_values = v_values[1:]
    #p00_stats = p00_stats[1:]

    fit_weibull_distribution_linear(v_values, p00_stats)
