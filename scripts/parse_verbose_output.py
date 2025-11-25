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


if __name__ == "__main__":
    report = KVMerReport("./ERR3152366_trim_ref.csv")
    #report = KVMerReport("./ERR3152366_ref.csv")
    #report = KVMerReport("./ERR2935851_trim_ref.csv")
    #report = KVMerReport("./SRR7415629_ref.csv")
    
    filt = report.report_data_df["homopolymer_length"] > 0
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
    

    