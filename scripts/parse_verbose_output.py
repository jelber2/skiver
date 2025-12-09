import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

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
            lambda_stats.append(error_count / consensus_count if consensus_count > 0 else 0)
            v += 1
        return v_values, lambda_stats
    

    def calculate_p00_stats(self, filter=None):
        p00_stats = []
        v_values = []
        v = 1
        data_filter = filter if filter is not None else self.report_data_df["total_count"] > 0

        while f"consensus_count_up_to_v{v}" in self.report_data_df.columns:
            v_values.append(v)
            if v == 1:
                p0_count = self.report_data_df[data_filter][f"total_count"].sum()
                p00_count = self.report_data_df[data_filter][f"consensus_count_up_to_v{v}"].sum()
                p00_stats.append(p00_count / p0_count if p0_count > 0 else 0)
            else:
                p0_count = self.report_data_df[data_filter][f"consensus_count_up_to_v{v-1}"].sum()
                p00_count = self.report_data_df[data_filter][f"consensus_count_up_to_v{v}"].sum()
                p00_stats.append(p00_count / p0_count if p0_count > 0 else 0)
            v += 1
        return v_values, p00_stats
    
    def analyze_with_plot(self, filter=None):
        # Plot lambda vs. v and p00 vs. v
        v_values, lambda_stats = self.calculate_lambda_stats(filter=filter)
        v_values, p00_stats = self.calculate_p00_stats(filter=filter)
        plt.figure(figsize=(12, 6))
        
        plt.subplot(1, 2, 1)
        x = np.array(v_values)

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
            lr_p00 = stats.linregress(x, y_p00)
            y_p00_fit = lr_p00.intercept + lr_p00.slope * x
            plt.subplot(1, 2, 2)
            plt.plot(x, y_p00_fit, color='red', linestyle='--',
                 label=f"fit: y={lr_p00.slope:.4f}x+{lr_p00.intercept:.4f}, R²={lr_p00.rvalue**2:.4f}")
            plt.plot(x, y_p00, marker='o')
            plt.legend()

        
        

            plt.title('P00 vs. v')
            plt.xlabel('v')
            plt.ylabel('P00')
            plt.ylim(0.8, 1.0)
            plt.tight_layout()
            plt.show()
    
    def _find_mutation_outlier(self):
        fields = ["SUBSTITUTION(AC)","SUBSTITUTION(AG)","SUBSTITUTION(AT)","SUBSTITUTION(GA)","SUBSTITUTION(GC)","SUBSTITUTION(GT)","SUBSTITUTION(CA)","SUBSTITUTION(CG)","SUBSTITUTION(CT)","SUBSTITUTION(TA)","SUBSTITUTION(TC)","SUBSTITUTION(TG)","INSERTION(A)","INSERTION(C)","INSERTION(G)","INSERTION(T)","DELETION(A)","DELETION(C)","DELETION(G)","DELETION(T)"]

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








if __name__ == "__main__":
    #report = KVMerReport("./ERR3152366_trim_ref.csv")
    report = KVMerReport("./ERR3152366.csv")
    #report = KVMerReport("./ERR2935851_trim_ref.csv")
    #report = KVMerReport("./SRR7415629_ref.csv")
    
    #filt = (report.report_data_df["homopolymer_length"] > 0) & (report.report_data_df["consensus_count"] >= 5)
    filt = ~report._find_consensus_outliers() & (report.report_data_df["consensus_count"] >= 5)
    #filt = report.report_data_df["total_count"] > 5
    v_values, lambda_stats = report.calculate_lambda_stats(filter=filt)
    #lambda_regression = report._linear_regression(v_values, lambda_stats)
    print("Lambda Stats:", lambda_stats)
    #print("Lambda Regression:", lambda_regression)

    v_values, p00_stats = report.calculate_p00_stats(filter=filt)
    #p00_regression = report._linear_regression(v_values, p00_stats)
    print("P00 Stats:", p00_stats)
    #print("P00 Regression:", p00_regression)

    report.analyze_with_plot(filter=filt)
    

    