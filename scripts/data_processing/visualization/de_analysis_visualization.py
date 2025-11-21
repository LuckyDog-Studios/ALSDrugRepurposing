import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import scipy.stats as stats

# Set style for better plots
plt.style.use('default')
sns.set_palette("viridis")

# Load the data
df = pd.read_csv('results\\combined_als_de_results_FIXED.csv')

print("Dataset Overview:")
print(f"Shape: {df.shape}")
print(f"Columns: {df.columns.tolist()}")
print(f"Number of genes: {df['gene'].nunique()}")
print(f"Datasets: {df['dataset'].unique()}")
print(f"Descriptions: {df['description'].unique()}")

# Data cleaning and preparation
print("\n" + "="*50)
print("DATA CLEANING AND PREPARATION")
print("="*50)

# Fix the 'significant' column - handle mixed types and missing values
print(f"Original 'significant' column dtype: {df['significant'].dtype}")
print(f"Unique values in 'significant': {df['significant'].unique()[:10]}")  # Show first 10 unique values

# Convert 'significant' to boolean, handling mixed types
df['significant'] = pd.to_numeric(df['significant'], errors='coerce')
df['significant'] = df['significant'].fillna(0).astype(bool)

print(f"Cleaned 'significant' column - True: {df['significant'].sum()}, False: {(~df['significant']).sum()}")

# Fill missing sample sizes with mode or appropriate values
if df['n_ALS'].isnull().any():
    n_als_mode = df['n_ALS'].mode()
    if len(n_als_mode) > 0:
        df['n_ALS'] = df['n_ALS'].fillna(n_als_mode[0])
    else:
        df['n_ALS'] = df['n_ALS'].fillna(3)  # Default value

if df['n_Control'].isnull().any():
    n_control_mode = df['n_Control'].mode()
    if len(n_control_mode) > 0:
        df['n_Control'] = df['n_Control'].fillna(n_control_mode[0])
    else:
        df['n_Control'] = df['n_Control'].fillna(3)  # Default value

print(f"Missing values after cleaning:")
print(df.isnull().sum())

# Check data quality
print(f"\nDATA QUALITY CHECK:")
print(f"P-value range: {df['p_value'].min():.2e} to {df['p_value'].max():.2f}")
print(f"Adj p-value range: {df['adj_p_value'].min():.2e} to {df['adj_p_value'].max():.2f}")
print(f"Log2FC range: {df['log2FC'].min():.3f} to {df['log2FC'].max():.3f}")
print(f"Significant genes: {df['significant'].sum()} / {len(df)}")
print(f"Significance rate: {(df['significant'].sum() / len(df) * 100):.2f}%")

# Define significance thresholds for visualization
fc_threshold = 0.5
p_threshold = 0.05

# Create comprehensive analysis plots
fig = plt.figure(figsize=(20, 20))

# 1. Distribution of log2 fold changes by dataset
plt.subplot(4, 3, 1)
datasets = df['dataset'].unique()
colors = plt.cm.viridis(np.linspace(0, 1, len(datasets)))

for i, dataset in enumerate(datasets):
    dataset_data = df[df['dataset'] == dataset]
    plt.hist(dataset_data['log2FC'], bins=30, alpha=0.6, label=dataset, 
             color=colors[i], density=True)

plt.axvline(x=0, color='red', linestyle='--', alpha=0.8)
plt.xlabel('log2 Fold Change')
plt.ylabel('Density')
plt.title('Distribution of log2 Fold Changes by Dataset')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)

# 2. Distribution of p-values (-log10 transformed)
plt.subplot(4, 3, 2)
for i, dataset in enumerate(datasets):
    dataset_data = df[df['dataset'] == dataset]
    neg_log_pvals = -np.log10(dataset_data['p_value'])
    plt.hist(neg_log_pvals, bins=30, alpha=0.6, label=dataset, 
             color=colors[i], density=True)

plt.axvline(x=-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
plt.axvline(x=-np.log10(0.01), color='orange', linestyle='--', label='p=0.01')
plt.xlabel('-log10(p-value)')
plt.ylabel('Density')
plt.title('Distribution of -log10 P-values by Dataset')
plt.legend()
plt.grid(True, alpha=0.3)

# 3. Volcano plot for all data
plt.subplot(4, 3, 3)
colors = ['gray' if not sig else 'red' for sig in df['significant']]
plt.scatter(df['log2FC'], -np.log10(df['p_value']), c=colors, alpha=0.6, s=10)
plt.axhline(y=-np.log10(p_threshold), color='red', linestyle='--', alpha=0.8, label=f'p={p_threshold}')
plt.axvline(x=-fc_threshold, color='blue', linestyle='--', alpha=0.8, label=f'|FC|>{fc_threshold}')
plt.axvline(x=fc_threshold, color='blue', linestyle='--', alpha=0.8)
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot - All Datasets')
plt.legend()
plt.grid(True, alpha=0.3)

# 4. MA plot (Mean expression vs Fold Change)
plt.subplot(4, 3, 4)
mean_expression = (df['mean_ALS'] + df['mean_Control']) / 2
colors = ['gray' if not sig else 'red' for sig in df['significant']]
plt.scatter(mean_expression, df['log2FC'], c=colors, alpha=0.6, s=10)
plt.axhline(y=0, color='black', linestyle='-', alpha=0.8)
plt.axhline(y=fc_threshold, color='blue', linestyle='--', alpha=0.8)
plt.axhline(y=-fc_threshold, color='blue', linestyle='--', alpha=0.8)
plt.xlabel('Mean Expression (ALS + Control)/2')
plt.ylabel('log2 Fold Change')
plt.title('MA Plot - All Datasets')
plt.grid(True, alpha=0.3)

# 5. P-value vs adjusted p-value
plt.subplot(4, 3, 5)
plt.scatter(-np.log10(df['p_value']), -np.log10(df['adj_p_value']), alpha=0.6, s=10)
max_val = max(-np.log10(df['p_value']).max(), -np.log10(df['adj_p_value']).max())
plt.plot([0, max_val], [0, max_val], 'r--', alpha=0.8, label='y=x')
plt.xlabel('-log10(p-value)')
plt.ylabel('-log10(adjusted p-value)')
plt.title('P-value vs Adjusted P-value')
plt.legend()
plt.grid(True, alpha=0.3)

# 6. Top significant genes across all datasets
plt.subplot(4, 3, 6)
sig_genes = df[df['significant']].nlargest(20, 'neg_log10_pval')
if len(sig_genes) > 0:
    y_pos = range(len(sig_genes))
    colors = ['red' if fc > 0 else 'blue' for fc in sig_genes['log2FC']]
    plt.barh(y_pos, sig_genes['log2FC'], color=colors, alpha=0.7)
    plt.yticks(y_pos, [f"{g} ({d})" for g, d in zip(sig_genes['gene'], sig_genes['dataset'])])
    plt.xlabel('log2 Fold Change')
    plt.title('Top 20 Significant Genes\n(Red=Up, Blue=Down)')
    plt.tight_layout()
else:
    plt.text(0.5, 0.5, 'No significant genes\nfound with current thresholds', 
             ha='center', va='center', transform=plt.gca().transAxes)
    plt.title('No Significant Genes')

# 7. Distribution of t-statistics by dataset
plt.subplot(4, 3, 7)
for i, dataset in enumerate(datasets):
    dataset_data = df[df['dataset'] == dataset]
    plt.hist(dataset_data['t_statistic'], bins=30, alpha=0.6, label=dataset, 
             color=colors[i], density=True)

plt.axvline(x=0, color='red', linestyle='--', alpha=0.8)
plt.xlabel('t-statistic')
plt.ylabel('Density')
plt.title('Distribution of t-statistics by Dataset')
plt.legend()
plt.grid(True, alpha=0.3)

# 8. Correlation between mean expression in ALS and Control
plt.subplot(4, 3, 8)
for i, dataset in enumerate(datasets):
    dataset_data = df[df['dataset'] == dataset]
    plt.scatter(dataset_data['mean_ALS'], dataset_data['mean_Control'], 
                alpha=0.6, s=10, label=dataset, color=colors[i])

max_mean = max(df['mean_ALS'].max(), df['mean_Control'].max())
min_mean = min(df['mean_ALS'].min(), df['mean_Control'].min())
plt.plot([min_mean, max_mean], [min_mean, max_mean], 'r--', alpha=0.8, label='y=x')
plt.xlabel('Mean Expression - ALS')
plt.ylabel('Mean Expression - Control')
plt.title('ALS vs Control Mean Expression by Dataset')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)

# 9. Significance by dataset
plt.subplot(4, 3, 9)
sig_by_dataset = df.groupby('dataset')['significant'].agg(['sum', 'count'])
sig_by_dataset['percentage'] = (sig_by_dataset['sum'] / sig_by_dataset['count'] * 100)
sig_by_dataset = sig_by_dataset.sort_values('percentage', ascending=False)

plt.bar(range(len(sig_by_dataset)), sig_by_dataset['percentage'], 
        color=plt.cm.viridis(np.linspace(0, 1, len(sig_by_dataset))))
plt.xticks(range(len(sig_by_dataset)), sig_by_dataset.index, rotation=45)
plt.ylabel('Percentage Significant (%)')
plt.title('Percentage of Significant Genes by Dataset')
plt.grid(True, alpha=0.3)

# 10. Effect size distribution by dataset
plt.subplot(4, 3, 10)
dataset_data = [df[df['dataset'] == dataset]['log2FC'] for dataset in datasets]
plt.boxplot(dataset_data, labels=datasets)
plt.xticks(rotation=45)
plt.ylabel('log2 Fold Change')
plt.title('Effect Size Distribution by Dataset')
plt.grid(True, alpha=0.3)

# 11. Sample size information
plt.subplot(4, 3, 11)
sample_sizes = df.groupby('dataset')[['n_ALS', 'n_Control']].mean()
x = range(len(sample_sizes))
width = 0.35
plt.bar([i - width/2 for i in x], sample_sizes['n_ALS'], width, label='ALS', alpha=0.7)
plt.bar([i + width/2 for i in x], sample_sizes['n_Control'], width, label='Control', alpha=0.7)
plt.xticks(x, sample_sizes.index, rotation=45)
plt.ylabel('Average Sample Size')
plt.title('Average Sample Sizes by Dataset')
plt.legend()
plt.grid(True, alpha=0.3)

# 12. QC: P-value uniformity under null (QQ-plot)
plt.subplot(4, 3, 12)
stats.probplot(df['p_value'], dist="uniform", plot=plt)
plt.title('QQ-plot of P-values\n(Should follow diagonal if null holds)')

plt.tight_layout()
plt.show()

# Summary statistics by dataset
print("\n" + "="*50)
print("SUMMARY STATISTICS BY DATASET")
print("="*50)

for dataset in datasets:
    dataset_df = df[df['dataset'] == dataset]
    total_genes = len(dataset_df)
    sig_genes = dataset_df['significant'].sum()
    sig_up = ((dataset_df['significant']) & (dataset_df['log2FC'] > 0)).sum()
    sig_down = ((dataset_df['significant']) & (dataset_df['log2FC'] < 0)).sum()
    
    print(f"\n{dataset}:")
    print(f"  Total genes: {total_genes}")
    print(f"  Significant genes: {sig_genes} ({sig_genes/total_genes*100:.2f}%)")
    print(f"    - Up-regulated: {sig_up}")
    print(f"    - Down-regulated: {sig_down}")
    print(f"  Avg |log2FC|: {np.abs(dataset_df['log2FC']).mean():.3f}")
    print(f"  Avg sample size - ALS: {dataset_df['n_ALS'].mean():.1f}, Control: {dataset_df['n_Control'].mean():.1f}")

# Overall summary
print("\n" + "="*50)
print("OVERALL SUMMARY")
print("="*50)

total_genes = len(df)
sig_genes = df['significant'].sum()
sig_up = ((df['significant']) & (df['log2FC'] > 0)).sum()
sig_down = ((df['significant']) & (df['log2FC'] < 0)).sum()

print(f"Total genes across all datasets: {total_genes}")
print(f"Total significant genes: {sig_genes} ({sig_genes/total_genes*100:.2f}%)")
print(f"  - Up-regulated: {sig_up}")
print(f"  - Down-regulated: {sig_down}")

# Top significant genes overall
if df['significant'].sum() > 0:
    print(f"\nTop 10 most significant genes overall:")
    top_sig = df[df['significant']].nlargest(10, 'neg_log10_pval')[['gene', 'dataset', 'log2FC', 'p_value', 'adj_p_value']]
    for _, row in top_sig.iterrows():
        print(f"  {row['gene']} ({row['dataset']}): FC={row['log2FC']:.3f}, p={row['p_value']:.2e}, adj_p={row['adj_p_value']:.2e}")

# Additional diagnostic plots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

# Histogram of raw vs adjusted p-values
ax1.hist(df['p_value'], bins=50, alpha=0.7, label='Raw p-values', color='blue', density=True)
ax1.hist(df['adj_p_value'], bins=50, alpha=0.7, label='Adj p-values', color='red', density=True)
ax1.set_xlabel('P-value')
ax1.set_ylabel('Density')
ax1.set_title('Distribution of Raw vs Adjusted P-values')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Scatter plot of effect size vs multiple testing corrected significance
sig_colors = ['red' if sig else 'gray' for sig in df['significant']]
ax2.scatter(df['log2FC'], -np.log10(df['adj_p_value']), c=sig_colors, alpha=0.6, s=10)
ax2.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='FDR=0.05')
ax2.axvline(x=0, color='black', linestyle='-', alpha=0.8)
ax2.set_xlabel('log2 Fold Change')
ax2.set_ylabel('-log10(Adjusted P-value)')
ax2.set_title('Effect Size vs Multiple Testing Corrected Significance')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Dataset comparison - significance rates
sig_rates = df.groupby('dataset')['significant'].mean() * 100
ax3.bar(range(len(sig_rates)), sig_rates.values, color=plt.cm.viridis(np.linspace(0, 1, len(sig_rates))))
ax3.set_xticks(range(len(sig_rates)))
ax3.set_xticklabels(sig_rates.index, rotation=45)
ax3.set_ylabel('Significant Genes (%)')
ax3.set_title('Significance Rates by Dataset')
ax3.grid(True, alpha=0.3)

# Effect size comparison
effect_sizes = df.groupby('dataset')['log2FC'].agg(['mean', 'std'])
ax4.bar(range(len(effect_sizes)), effect_sizes['mean'], 
        yerr=effect_sizes['std'], capsize=5, 
        color=plt.cm.viridis(np.linspace(0, 1, len(effect_sizes))))
ax4.set_xticks(range(len(effect_sizes)))
ax4.set_xticklabels(effect_sizes.index, rotation=45)
ax4.set_ylabel('Mean log2 Fold Change')
ax4.set_title('Average Effect Sizes by Dataset')
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

print(f"\nANALYSIS COMPLETE!")
print(f"Combined dataset contains {len(df)} genes from {len(datasets)} different ALS studies")
print(f"Overall significance rate: {(df['significant'].sum() / len(df) * 100):.2f}%")