import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns


ALPHA_BETA_SUM_CONSTRAINT = 20

class KVMerReport:
    def __init__(self, report_data_df):
        self.report_data_df = pd.read_csv(report_data_df)
    
    def calculate_lambda_stats(self, filter=None):
        lambda_stats = []
        v_values = []
        v = 1
        data_filter = filter if filter is not None else self.report_data_df["total_count"] > 0

        while f"consensus_count_up_to_v{v}" in self.report_data_df.columns:
            v_values.append(v)
            consensus_count = self.report_data_df[data_filter][f"consensus_count_up_to_v{v}"].sum()
            error_count = self.report_data_df[data_filter][f"error_count_up_to_v{v}"].sum()
            lambda_stats.append(float(error_count / consensus_count) if consensus_count > 0 else 0)
            v += 1
        return v_values, lambda_stats
    

    def calculate_p00_stats(self, even_weighted=False, filter=None):
        p00_stats = []
        v_values = []
        v = 1
        data_filter = filter if filter is not None else self.report_data_df["total_count"] > 0

        print(f"{len(self.report_data_df[data_filter])} entries after filtering.")

        while f"consensus_count_up_to_v{v}" in self.report_data_df.columns:
            v_values.append(v)
            if v == 1:
                p0_count = self.report_data_df[data_filter][f"total_count"]
                p00_count = self.report_data_df[data_filter][f"consensus_count_up_to_v{v}"]
                #p00_stats.append(float(p00_count / p0_count) if p0_count > 0 else 0)
                #print(f"v={v}, p0_count={p0_count}, p00_count={p00_count}")
            else:
                p0_count = self.report_data_df[data_filter][f"consensus_count_up_to_v{v-1}"]
                p00_count = self.report_data_df[data_filter][f"consensus_count_up_to_v{v}"]
                #p00_stats.append(float(p00_count / p0_count) if p0_count > 0 else 0)
                #print(f"v={v}, p0_count={p0_count}, p00_count={p00_count}")
            if even_weighted:
                ratios = p00_count / p0_count.replace(0, np.nan)
                mean_ratio = ratios.mean()
                p00_stats.append(float(mean_ratio) if not np.isnan(mean_ratio) else 0)
            else:
                total_p0 = p0_count.sum()
                total_p00 = p00_count.sum()
                p00_stats.append(float(total_p00 / total_p0) if total_p0 > 0 else 0)
            v += 1
        return v_values, p00_stats
    
    def analyze_with_plot(self, filter=None):
        # Plot lambda vs. v and p00 vs. v
        v_values, lambda_stats = self.calculate_lambda_stats(filter=filter)
        v_values, p00_stats = self.calculate_p00_stats(filter=filter)
        plt.figure(figsize=(12, 6))
        
        plt.subplot(1, 2, 1)
        x = np.array(v_values)
        k = len(self.report_data_df["key"].iloc[0])

        # Lambda regression and plot its fit on subplot 1
        y_lambda = np.array(lambda_stats)
        if x.size > 1 and y_lambda.size == x.size:
            lr_lambda = stats.linregress(x[1:], y_lambda[1:])
            y_lambda_fit = lr_lambda.intercept + lr_lambda.slope * x
            plt.plot(x, y_lambda_fit, color='red', linestyle='--',
                 label=f"fit: y={lr_lambda.slope:.4f}x+{lr_lambda.intercept:.4f}, R²={lr_lambda.rvalue**2:.4f}")
            plt.plot(x, y_lambda, marker='o')
            plt.legend()
            plt.title('Lambda vs. v')
            plt.xlabel('v')
            plt.ylabel('Lambda')

        # P00 regression: compute and plot its fit on subplot 2
        y_p00 = np.array(p00_stats)
        if x.size > 1 and y_p00.size == x.size:
            """
            # Fit with beta distribution
            fit = self.fit_beta_distribution(x + k, y_p00)
            if fit[0] + fit[1] < ALPHA_BETA_SUM_CONSTRAINT:
                fit = self.fit_beta_distribution_constrain(x + k, y_p00, alpha_beta_sum=ALPHA_BETA_SUM_CONSTRAINT)
            y_p00_fit = 1 - fit[0] / (fit[0] + fit[1] + x + k)
            """

            # Fit with Weibull distribution
            fit = self.fit_weibull_distribution(x + k, y_p00)
            y_p00_fit = 1 - fit[0] * ((x + k) ** fit[1])

           
            plt.subplot(1, 2, 2)
            plt.plot(x + k, y_p00_fit, color='red', linestyle='--',
                 label=f"fit: alpha={fit[0]:.2f}, beta={fit[1]:.2f}")
            plt.plot(x + k, y_p00, marker='o')
            plt.legend()

        
        

            plt.title('P00 vs. v')
            plt.xlabel('v')
            plt.ylabel('P00')
            #plt.ylim(0.8, 1.0)
            plt.tight_layout()
            plt.savefig("lambda_p00_analysis.png", transparent=True)
            plt.show()
    
    def _find_mutation_outlier(self):
        fields = ["AC","AG","AT","GA","GC","GT","CA","CG","CT","TA","TC","TG",
                  "_A","_C","_G","_T",
                  "A_","C_","G_","T_"]

        outliers = self.report_data_df["total_count"] < 0

        for field in fields:
            self.report_data_df["ratio"] = self.report_data_df[field] / self.report_data_df["total_count"]

            nonzero_ratios = self.report_data_df.loc[self.report_data_df["ratio"] > 0, "ratio"]

            q1 = nonzero_ratios.quantile(0.25)
            q3 = nonzero_ratios.quantile(0.75)
            iqr = q3 - q1

            
            upper_bound = q3 + 1.5 * iqr

            print(f"Field: {field}, Upper Bound: {upper_bound}")

            num_outliers = ((self.report_data_df["ratio"] > upper_bound)).sum()
            print(f"Number of outliers in {field}: {num_outliers}")

            outliers = outliers | (self.report_data_df["ratio"] > upper_bound)
        
        print(f"Total number of outliers: {outliers.sum()}")
        
        return outliers
    
    def _find_consensus_outliers(self):
        v = 1
        outliers = self.report_data_df["total_count"] < 0

        while f"consensus_count_up_to_v{v}" in self.report_data_df.columns:
            if v == 1:
                self.report_data_df["ratio"] = self.report_data_df[f"consensus_count_up_to_v{v}"] / self.report_data_df["total_count"]
            else:
                self.report_data_df["ratio"] = self.report_data_df[f"consensus_count_up_to_v{v}"] / self.report_data_df[f"consensus_count_up_to_v{v-1}"]
            
            nonone_ratios = self.report_data_df.loc[self.report_data_df["ratio"] < 1, "ratio"]
            q1 = nonone_ratios.quantile(0.25)
            q3 = nonone_ratios.quantile(0.75)
            iqr = q3 - q1
            lower_bound = q1 - 1.5 * iqr
            print(f"v: {v}, Lower Bound: {lower_bound}")
            num_outliers = ((self.report_data_df["ratio"] < lower_bound)).sum()
            print(f"Number of outliers in consensus ratio up to v{v}: {num_outliers}")

            outliers = outliers | (self.report_data_df["ratio"] < lower_bound)
            v += 1

        print(f"Total number of outliers: {outliers.sum()}")
        return outliers
    
    def _find_majority_outliers(self):
        outliers = self.report_data_df["total_count"] < 0
        self.report_data_df["ratio"] = self.report_data_df["consensus_count"] / self.report_data_df["total_count"]

        nonone_ratios = self.report_data_df.loc[self.report_data_df["ratio"] <= 1, "ratio"]
        q1 = nonone_ratios.quantile(0.25)
        q3 = nonone_ratios.quantile(0.75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr

        # Plot consensus distribution
        plt.figure(figsize=(8, 6))
        sns.histplot(self.report_data_df["ratio"], bins=50, kde=True)

        # Add a vertical line for the lower bound
        plt.axvline(x=lower_bound, color='red', linestyle='--', label='Lower Bound')
        plt.legend()
        plt.title('Majority Consensus Ratio Distribution')
        plt.xlabel('Majority Consensus Ratio')
        plt.ylabel('Frequency')
        plt.show()

        
        print(f"Majority Lower Bound: {lower_bound}")
        num_outliers = ((self.report_data_df["ratio"] < lower_bound)).sum()
        print(f"Number of outliers in majority consensus ratio: {num_outliers}")
        outliers = outliers | (self.report_data_df["ratio"] < lower_bound)
        print(f"Total number of outliers: {outliers.sum()}")
        return outliers
    
    def plot_consensus_distribution(self, v):
        if f"consensus_count_up_to_v{v}" in self.report_data_df.columns:
            if v == 1:
                self.report_data_df["ratio"] = self.report_data_df[f"consensus_count_up_to_v{v}"] / self.report_data_df["total_count"]
            else:
                self.report_data_df["ratio"] = self.report_data_df[f"consensus_count_up_to_v{v}"] / self.report_data_df[f"consensus_count_up_to_v{v-1}"]
            
            plt.figure(figsize=(8, 6))
            sns.histplot(self.report_data_df[self.report_data_df["ratio"] < 1]["ratio"], bins=50, kde=True)
            plt.title(f'Hazard ratio at v={v}')
            plt.xlabel('Hazard Ratio')
            plt.ylabel('Frequency')
            plt.show()
    

    def plot_mutation_spectrum(self, filter=None):
        substitution_fields = ["AC","AG","AT","GA","GC","GT","CA","CG","CT","TA","TC","TG"]
        insertion_fields = ["_A","_C","_G","_T"]
        deletion_fields = ["A_","C_","G_","T_"]

        filt = filter if filter is not None else self.report_data_df["total_count"] > 0

        mutation_spectrum = {}

        for field_set, title in zip([substitution_fields, insertion_fields, deletion_fields], ["Substitutions", "Insertions", "Deletions"]):
            for field in field_set:
                y_sum = self.report_data_df[filt][field].sum()
                x_sum = self.report_data_df[filt]["consensus_count"].sum() +  self.report_data_df[filt]["neighbor_count"].sum() - y_sum

                ratio = y_sum / x_sum if x_sum > 0 else 0
                mutation_spectrum[field] = ratio

                print(f"{field}: {y_sum} / {x_sum} = {ratio:.6f}")
        
        # Normalize so that the total sums to 1
        total = sum(mutation_spectrum.values())
        for key in mutation_spectrum:
            mutation_spectrum[key] /= total
        
        # Plot mutation spectrum
        # On the left plot, a heat map with [A, C, G, T, _] on both axes for substitutions, insertions, and deletions, and on the right plot, a bar plot showing the substitution, insertion, and deletion rates.
        plt.figure(figsize=(14, 6))
        plt.subplot(1, 2, 1)
        bases = ['A', 'C', 'G', 'T', '_']
        matrix = np.zeros((5, 5))
        for field, ratio in mutation_spectrum.items():
            from_base = field[0]
            to_base = field[1]
            i = bases.index(from_base)
            j = bases.index(to_base)
            matrix[i, j] = ratio
            
        sns.heatmap(matrix, xticklabels=bases, yticklabels=bases, annot=True, fmt=".4f", cmap="Blues")
        plt.title('Error/Mutation Spectrum Heatmap')
        plt.xlabel('Observed Base')
        plt.ylabel('Reference Base')
        plt.subplot(1, 2, 2)
        mutation_types = ['Substitutions', 'Insertions', 'Deletions']
        mutation_rates = [sum(mutation_spectrum[field] for field in substitution_fields),
                          sum(mutation_spectrum[field] for field in insertion_fields),
                          sum(mutation_spectrum[field] for field in deletion_fields)]
        sns.barplot(x=mutation_types, y=mutation_rates)
        plt.title('Overall Error/Mutation Spectrum')
        plt.ylabel('Proportion')
        plt.tight_layout()
        plt.savefig("mutation_spectrum.png", transparent=True)
        plt.show()  
        
    

    def fit_beta_distribution(self, v_values, p00_stats):
        # fit 1 - alpha / (alpha + beta + v) to p00_stats
        from scipy.optimize import curve_fit
        def beta_func(v, alpha, beta):
            return 1 - alpha / (alpha + beta + v)
        
        popt, pcov = curve_fit(beta_func, v_values, p00_stats, bounds=(0, [1000., 1000.]))
        alpha, beta = popt

        print(f"Fitted beta distribution parameters: alpha={alpha}, beta={beta}")
        print(f"Mean: {alpha / (alpha + beta)}, Variance: {(alpha * beta) / ((alpha + beta)**2 * (alpha + beta + 1))}")


        return popt  # alpha, beta
    
    def fit_weibull_distribution(self, v_values, p00_stats):

        # fit 1 - a * (v ** b) to p00_stats
        # Or, fit a linear regression to ln(1 - p00_stats) vs. ln(v)

        y = np.array([np.log(1 - p) if p < 1 else -np.inf for p in p00_stats])
        x = np.log(np.array(v_values))

        # Use ridge regression to fit the line
        from sklearn.linear_model import Ridge

        model = Ridge(alpha=1.0)
        x_reshaped = x.reshape(-1, 1)
        model.fit(x_reshaped, y)
        b = model.coef_[0]
        a = model.intercept_

        print(f"Fitted Weibull distribution parameters: a={np.exp(a)}, b={b}")
        return np.exp(a), b  # a, b



    def plot_coverage_distribution(self):
        print("Median coverage:", self.report_data_df[self.report_data_df["total_count"] < 200]["total_count"].median())
        print("95th percentile coverage:", self.report_data_df[self.report_data_df["total_count"] < 200]["total_count"].quantile(0.95))

        plt.figure(figsize=(8, 6))
        sns.histplot(self.report_data_df[self.report_data_df["total_count"] < 200]["total_count"], bins=50, kde=True)
        plt.title('Coverage Distribution')
        plt.xlabel('Total Count (Coverage)')
        plt.ylabel('Frequency')
        plt.show()

    def fit_beta_distribution_constrain(self, v_values, p00_stats, alpha_beta_sum=10):
        # fit 1 - alpha / (alpha + beta_fixed + v) to p00_stats
        from scipy.optimize import curve_fit
        def beta_func(v, alpha):
            return 1 - alpha / (alpha_beta_sum + v)
        
        popt, pcov = curve_fit(beta_func, v_values, p00_stats, bounds=(0, [1000.]))
        alpha = popt[0]

        print(f"Fitted beta distribution parameters with beta fixed at {alpha_beta_sum - alpha}: alpha={alpha}")
        print(f"Mean: {alpha / (alpha + (alpha_beta_sum - alpha))}, Variance: {(alpha * (alpha_beta_sum - alpha)) / ((alpha + (alpha_beta_sum - alpha))**2 * (alpha + (alpha_beta_sum - alpha) + 1))}")


        return alpha, alpha_beta_sum - alpha  # alpha, beta
    

    def estimation_with_different_coverage(self, filter=None):
        coverage_range = (min(self.report_data_df["total_count"]), max(self.report_data_df["total_count"]))
        coverage_sort = sorted(self.report_data_df["total_count"].unique())
        if len(coverage_sort) > 500:
            coverage_sort = coverage_sort[:500]
        mean_res = []
        variance_res = []

        for coverage in coverage_sort:
            if filter is not None:
                combined_filter = filter & (self.report_data_df["total_count"] >= coverage)
            else:
                combined_filter = self.report_data_df["total_count"] >= coverage
            v_values, p00_stats = self.calculate_p00_stats(even_weighted=False, filter=combined_filter)
            fit = self.fit_beta_distribution(v_values, p00_stats)

            if fit[1] < ALPHA_BETA_SUM_CONSTRAINT:
                fit = self.fit_beta_distribution_constrain(v_values, p00_stats, alpha_beta_sum=ALPHA_BETA_SUM_CONSTRAINT)

            mean_res.append(fit[0] / (fit[0] + fit[1]))
            variance_res.append(np.sqrt((fit[0] * fit[1]) / ((fit[0] + fit[1])**2 * (fit[0] + fit[1] + 1))))

        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.plot(coverage_sort, mean_res)
        plt.title('Mean vs. Minimum Coverage')
        plt.xlabel('Minimum Coverage')
        plt.ylabel('Mean')

        plt.subplot(1, 2, 2)
        plt.plot(coverage_sort, variance_res)
        plt.title('Variance vs. Minimum Coverage')
        plt.xlabel('Minimum Coverage')
        plt.ylabel('Variance')
        plt.tight_layout()
        plt.show()



if __name__ == "__main__":
    #report = KVMerReport("./ERR3152366_ref.csv")
    #report = KVMerReport("./ERR3152366_ref.csv")
    #report = KVMerReport("./ERR2935851.csv")
    #report = KVMerReport("./SRR7415629.csv")
    #report = KVMerReport("./HG002.csv")
    report = KVMerReport("./test_90.csv")
    #report = KVMerReport("/home/ubuntu/kv-mer-test/output/multiple_alleles/two_strain_output.csv")
    #report = KVMerReport("/home/ubuntu/kv-mer-test/output/multiple_alleles/K12_MG1655_output.csv")
    #report = KVMerReport("/home/ubuntu/kv-mer-test/output/multiple_alleles/O157_H7_output.csv")
    #report = KVMerReport("./Ecoli_K12_MG1655_id_96.csv")

    #report.plot_consensus_distribution(v=1)
    report.plot_coverage_distribution()
    #filt = ~report._find_consensus_outliers()
    #report.estimation_with_different_coverage(filter=filt)

    #filt = (report.report_data_df["total_count"] >= 11)
    filt = (report.report_data_df["total_count"] >= 5) # & (report.report_data_df["homopolymer_length"] == 2)
    #filt = report.report_data_df["total_count"] > 5
    #v_values, lambda_stats = report.calculate_lambda_stats(filter=filt)
    #lambda_regression = report._linear_regression(v_values, lambda_stats)
    #print("Lambda Stats:", lambda_stats)
    #print("Lambda Regression:", lambda_regression)

    #v_values, p00_stats = report.calculate_p00_stats(filter=filt)
    #p00_regression = report._linear_regression(v_values, p00_stats)
    #print("P00 Stats:", p00_stats)
    report.plot_mutation_spectrum(filter=filt)
    #print("P00 Regression:", p00_regression)

    report.analyze_with_plot(filter=filt)
    

    