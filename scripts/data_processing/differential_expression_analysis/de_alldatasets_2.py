import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

class DifferentialExpressionAnalyzer:
    def __init__(self):
        self.results = {}
        self.combined_results = None
        
    def simulate_dataset(self, dataset_id, description, n_genes=5000, n_samples=20):
        """
        Simulate gene expression data for demonstration purposes.
        """
        print(f"Simulating data for {dataset_id}: {description}")
        
        # Simulate gene names
        genes = [f'GENE_{i:05d}' for i in range(n_genes)]
        
        # Split samples into ALS and control groups
        n_als = n_samples // 2
        n_control = n_samples - n_als
        
        # Create sample names
        als_samples = [f'{dataset_id}_ALS_{i}' for i in range(n_als)]
        control_samples = [f'{dataset_id}_CTRL_{i}' for i in range(n_control)]
        all_samples = als_samples + control_samples
        
        # Create expression matrix with some differentially expressed genes
        np.random.seed(hash(dataset_id) % 10000)
        
        expression_data = np.random.normal(8, 2, (n_genes, n_samples))
        
        # Introduce differential expression for some genes (5% of genes)
        n_de_genes = int(n_genes * 0.05)
        de_indices = np.random.choice(n_genes, n_de_genes, replace=False)
        
        for idx in de_indices:
            # ALS samples get higher or lower expression
            fold_change = np.random.choice([-1.5, 1.5])
            expression_data[idx, :n_als] *= fold_change
        
        # Create DataFrame
        df = pd.DataFrame(expression_data, index=genes, columns=all_samples)
        
        # Add group labels
        groups = ['ALS'] * n_als + ['Control'] * n_control
        
        return df, groups, genes
    
    def perform_de_analysis(self, expression_data, groups, dataset_id):
        """
        Perform differential expression analysis using t-tests.
        FIXED: Proper handling of log2 fold change calculation
        """
        print(f"Performing DE analysis for {dataset_id}...")
        
        als_mask = np.array(groups) == 'ALS'
        control_mask = np.array(groups) == 'Control'
        
        results = []
        
        for gene in expression_data.index:
            als_expression = expression_data.loc[gene, als_mask]
            control_expression = expression_data.loc[gene, control_mask]
            
            # T-test
            t_stat, p_value = stats.ttest_ind(als_expression, control_expression)
            
            # Calculate fold change (log2) with proper handling of zeros
            mean_als = np.mean(als_expression)
            mean_control = np.mean(control_expression)
            
            # FIX: Add pseudocount to avoid division by zero and log(0)
            pseudocount = 0.1  # Small pseudocount to handle zeros
            log2_fold_change = np.log2((mean_als + pseudocount) / (mean_control + pseudocount))
            
            results.append({
                'gene': gene,
                'log2_fold_change': log2_fold_change,
                'p_value': p_value,
                'mean_als': mean_als,
                'mean_control': mean_control,
                't_statistic': t_stat
            })
        
        results_df = pd.DataFrame(results)
        
        # Handle NaN p-values
        results_df['p_value'] = results_df['p_value'].fillna(1.0)
        
        # Calculate adjusted p-values (FDR)
        results_df['adj_p_value'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
        
        # Calculate -log10(p-value) for visualization
        results_df['neg_log10_pvalue'] = -np.log10(results_df['p_value'])
        
        # Sort by significance
        results_df = results_df.sort_values('p_value')
        
        return results_df
    
    def analyze_dataset(self, dataset_id, description):
        """
        Complete analysis pipeline for a single dataset.
        """
        print(f"\n{'='*60}")
        print(f"Analyzing: {dataset_id}")
        print(f"Description: {description}")
        print(f"{'='*60}")
        
        # Load/simulate data
        expression_data, groups, genes = self.simulate_dataset(dataset_id, description)
        
        # Perform DE analysis
        de_results = self.perform_de_analysis(expression_data, groups, dataset_id)
        
        # Store results
        self.results[dataset_id] = {
            'description': description,
            'expression_data': expression_data,
            'groups': groups,
            'de_results': de_results,
            'significant_genes': de_results[de_results['adj_p_value'] < 0.05]
        }
        
        # Print summary
        n_significant = len(self.results[dataset_id]['significant_genes'])
        print(f"Found {n_significant} significantly differentially expressed genes (FDR < 0.05)")
        
        # Check for missing values
        missing_fc = de_results['log2_fold_change'].isna().sum()
        missing_pval = de_results['p_value'].isna().sum()
        
        if missing_fc > 0:
            print(f"WARNING: {missing_fc} genes have missing fold change values")
        if missing_pval > 0:
            print(f"WARNING: {missing_pval} genes have missing p-values")
        
        return de_results
    
    def save_individual_datasets(self):
        """
        Save individual dataset results to separate Excel files.
        """
        print(f"\n{'='*60}")
        print("SAVING INDIVIDUAL DATASET RESULTS")
        print(f"{'='*60}")
        
        for dataset_id, result in self.results.items():
            filename = f"{dataset_id}_Differential_Expression_Results.xlsx"
            
            # Create Excel file with multiple sheets
            with pd.ExcelWriter(filename) as writer:
                # Main DE results
                result['de_results'].to_excel(writer, sheet_name='DE_Results', index=False)
                
                # Significant genes only
                sig_genes = result['significant_genes']
                sig_genes.to_excel(writer, sheet_name='Significant_Genes', index=False)
                
                # Upregulated genes
                upregulated = sig_genes[sig_genes['log2_fold_change'] > 0]
                upregulated.to_excel(writer, sheet_name='Upregulated_Genes', index=False)
                
                # Downregulated genes
                downregulated = sig_genes[sig_genes['log2_fold_change'] < 0]
                downregulated.to_excel(writer, sheet_name='Downregulated_Genes', index=False)
            
            print(f"Saved {dataset_id} results to {filename}")
            
            # Print quick summary
            n_up = len(upregulated)
            n_down = len(downregulated)
            print(f"  - Significant genes: {len(sig_genes)} (↑{n_up}, ↓{n_down})")
    
    def combine_results(self):
        """
        Combine results from all datasets using meta-analysis approach.
        """
        print(f"\n{'='*60}")
        print("COMBINING RESULTS FROM ALL DATASETS")
        print(f"{'='*60}")
        
        all_results = []
        
        for dataset_id, result in self.results.items():
            de_df = result['de_results'].copy()
            de_df['dataset'] = dataset_id
            de_df['description'] = result['description']
            all_results.append(de_df)
        
        # Combine all results
        combined_df = pd.concat(all_results, ignore_index=True)
        
        # Check for any missing values in combined data
        missing_summary = combined_df.isnull().sum()
        if missing_summary.any():
            print("Missing values in combined data:")
            print(missing_summary[missing_summary > 0])
        
        # Meta-analysis: Combine p-values using Fisher's method
        meta_results = []
        
        # Get unique genes across all datasets
        all_genes = combined_df['gene'].unique()
        
        for gene in all_genes:
            gene_data = combined_df[combined_df['gene'] == gene]
            
            if len(gene_data) > 1:
                # Fisher's method for combining p-values
                valid_pvals = gene_data['p_value'].dropna()
                if len(valid_pvals) > 0:
                    chi_square = -2 * np.sum(np.log(valid_pvals))
                    combined_p = stats.chi2.sf(chi_square, 2 * len(valid_pvals))
                else:
                    combined_p = 1.0
                
                # Average fold change (handle missing values)
                avg_log2_fc = gene_data['log2_fold_change'].mean()
                
                # Count number of datasets where significant
                n_sig = sum(gene_data['adj_p_value'] < 0.05)
                
            else:
                combined_p = gene_data['p_value'].iloc[0] if not pd.isna(gene_data['p_value'].iloc[0]) else 1.0
                avg_log2_fc = gene_data['log2_fold_change'].iloc[0] if not pd.isna(gene_data['log2_fold_change'].iloc[0]) else 0.0
                n_sig = 1 if gene_data['adj_p_value'].iloc[0] < 0.05 else 0
            
            meta_results.append({
                'gene': gene,
                'combined_p_value': combined_p,
                'avg_log2_fold_change': avg_log2_fc,
                'n_datasets_significant': n_sig,
                'n_datasets_tested': len(gene_data)
            })
        
        meta_df = pd.DataFrame(meta_results)
        
        # Adjust p-values for multiple testing
        meta_df['combined_adj_p_value'] = multipletests(
            meta_df['combined_p_value'], method='fdr_bh'
        )[1]
        
        # Sort by combined significance
        meta_df = meta_df.sort_values('combined_p_value')
        
        self.combined_results = {
            'individual_results': combined_df,
            'meta_analysis': meta_df,
            'significant_in_multiple': meta_df[meta_df['n_datasets_significant'] >= 2]
        }
        
        print(f"Meta-analysis completed:")
        print(f"- Total genes analyzed: {len(meta_df)}")
        print(f"- Genes significant in multiple datasets: {len(self.combined_results['significant_in_multiple'])}")
        
        return self.combined_results
    
    def save_combined_results(self):
        """
        Save combined results to Excel file.
        """
        print(f"\n{'='*60}")
        print("SAVING COMBINED RESULTS")
        print(f"{'='*60}")
        
        if self.combined_results is None:
            print("No combined results to save. Run combine_results() first.")
            return
        
        filename = "ALS_Combined_Differential_Expression_Results.xlsx"
        
        with pd.ExcelWriter(filename) as writer:
            # All individual results
            self.combined_results['individual_results'].to_excel(
                writer, sheet_name='All_Individual_Results', index=False
            )
            
            # Meta-analysis results
            self.combined_results['meta_analysis'].to_excel(
                writer, sheet_name='Meta_Analysis', index=False
            )
            
            # Genes significant in multiple datasets
            self.combined_results['significant_in_multiple'].to_excel(
                writer, sheet_name='Multi_Dataset_Significant', index=False
            )
            
            # Top 100 most significant genes
            top_100 = self.combined_results['meta_analysis'].head(100)
            top_100.to_excel(writer, sheet_name='Top_100_Genes', index=False)
        
        print(f"Combined results saved to {filename}")
    
    def create_detailed_summary(self):
        """
        Create detailed summary of all analyses.
        """
        print(f"\n{'='*60}")
        print("DETAILED ANALYSIS SUMMARY")
        print(f"{'='*60}")
        
        summary_data = []
        
        for dataset_id, result in self.results.items():
            de_df = result['de_results']
            sig_genes = result['significant_genes']
            
            n_total = len(de_df)
            n_sig = len(sig_genes)
            n_up = len(sig_genes[sig_genes['log2_fold_change'] > 0])
            n_down = len(sig_genes[sig_genes['log2_fold_change'] < 0])
            
            # Check data quality
            missing_fc = de_df['log2_fold_change'].isna().sum()
            missing_pval = de_df['p_value'].isna().sum()
            
            summary_data.append({
                'Dataset': dataset_id,
                'Description': result['description'],
                'Total_Genes': n_total,
                'Significant_Genes': n_sig,
                'Upregulated': n_up,
                'Downregulated': n_down,
                'Missing_FC': missing_fc,
                'Missing_Pval': missing_pval,
                'Percent_Significant': (n_sig / n_total) * 100
            })
        
        summary_df = pd.DataFrame(summary_data)
        print(summary_df.to_string(index=False))
        
        # Save summary to CSV
        summary_df.to_csv('ALS_DE_Analysis_Summary.csv', index=False)
        print(f"\nDetailed summary saved to 'ALS_DE_Analysis_Summary.csv'")
        
        return summary_df
    
    def create_summary_plot(self):
        """
        Create summary visualization of results.
        """
        if self.combined_results is None:
            print("Please run combine_results() first.")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot 1: Volcano plot for meta-analysis
        meta_df = self.combined_results['meta_analysis']
        sig_mask = meta_df['combined_adj_p_value'] < 0.05
        fc_threshold = 0.5
        
        axes[0,0].scatter(
            meta_df[~sig_mask]['avg_log2_fold_change'],
            -np.log10(meta_df[~sig_mask]['combined_p_value']),
            alpha=0.5, c='gray', label='Not significant'
        )
        axes[0,0].scatter(
            meta_df[sig_mask]['avg_log2_fold_change'],
            -np.log10(meta_df[sig_mask]['combined_p_value']),
            alpha=0.7, c='red', label='FDR < 0.05'
        )
        axes[0,0].axhline(-np.log10(0.05), linestyle='--', color='red', alpha=0.5)
        axes[0,0].axvline(-fc_threshold, linestyle='--', color='blue', alpha=0.5)
        axes[0,0].axvline(fc_threshold, linestyle='--', color='blue', alpha=0.5)
        axes[0,0].set_xlabel('Average log2(Fold Change)')
        axes[0,0].set_ylabel('-log10(p-value)')
        axes[0,0].set_title('Meta-Analysis: Volcano Plot')
        axes[0,0].legend()
        
        # Plot 2: Number of significant genes per dataset
        sig_counts = []
        dataset_names = []
        for dataset_id, result in self.results.items():
            sig_counts.append(len(result['significant_genes']))
            dataset_names.append(dataset_id)
        
        axes[0,1].bar(dataset_names, sig_counts, color='skyblue')
        axes[0,1].set_title('Significant Genes per Dataset')
        axes[0,1].set_ylabel('Number of Significant Genes (FDR < 0.05)')
        axes[0,1].tick_params(axis='x', rotation=45)
        
        # Plot 3: Distribution of fold changes
        all_fold_changes = []
        for dataset_id, result in self.results.items():
            all_fold_changes.extend(result['de_results']['log2_fold_change'].values)
        
        axes[1,0].hist(all_fold_changes, bins=50, alpha=0.7, color='lightgreen')
        axes[1,0].axvline(0, color='black', linestyle='-', alpha=0.5)
        axes[1,0].set_xlabel('log2(Fold Change)')
        axes[1,0].set_ylabel('Frequency')
        axes[1,0].set_title('Distribution of Fold Changes Across All Datasets')
        
        # Plot 4: Genes significant in multiple datasets
        multi_sig_counts = self.combined_results['meta_analysis']['n_datasets_significant'].value_counts().sort_index()
        axes[1,1].bar(multi_sig_counts.index, multi_sig_counts.values, color='orange')
        axes[1,1].set_xlabel('Number of Datasets Where Significant')
        axes[1,1].set_ylabel('Number of Genes')
        axes[1,1].set_title('Genes Significant in Multiple Datasets')
        
        plt.tight_layout()
        plt.savefig('ALS_DE_Analysis_Summary.png', dpi=300, bbox_inches='tight')
        plt.show()

def main():
    """
    Main function to run the differential expression analysis.
    """
    # Define datasets based on your project
    datasets = {
        'GSE124439': 'Spinal cord - human ALS patients',
        'GSE68605': 'Laser-captured motor neurons', 
        'GSE76253': 'Frontal cortex - ALS vs controls',
        'GSE833': 'Spinal cord - human ALS',
        'GSE56500': 'Muscle tissue from ALS patients',
        'GSE76220': 'Blood gene expression in ALS', 
        'GSE112676': 'iPSC motor neurons from ALS patients',
    }
    
    # Initialize analyzer
    analyzer = DifferentialExpressionAnalyzer()
    
    # Analyze each dataset
    for dataset_id, description in datasets.items():
        analyzer.analyze_dataset(dataset_id, description)
    
    # Save individual dataset results
    analyzer.save_individual_datasets()
    
    # Combine results
    analyzer.combine_results()
    
    # Save combined results
    analyzer.save_combined_results()
    
    # Create detailed summary
    analyzer.create_detailed_summary()
    
    # Create summary plot
    analyzer.create_summary_plot()
    
    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE!")
    print(f"{'='*60}")
    print("Results saved to:")
    print("- Individual dataset files: GSEXXXXX_Differential_Expression_Results.xlsx")
    print("- Combined results: ALS_Combined_Differential_Expression_Results.xlsx")
    print("- Summary statistics: ALS_DE_Analysis_Summary.csv")
    print("- Visualization: ALS_DE_Analysis_Summary.png")

if __name__ == "__main__":
    main()