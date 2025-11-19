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
        In practice, you would load real data from GEO using GEOparse.
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
        np.random.seed(hash(dataset_id) % 10000)  # Seed based on dataset ID
        
        expression_data = np.random.normal(8, 2, (n_genes, n_samples))
        
        # Introduce differential expression for some genes (5% of genes)
        n_de_genes = int(n_genes * 0.05)
        de_indices = np.random.choice(n_genes, n_de_genes, replace=False)
        
        for idx in de_indices:
            # ALS samples get higher or lower expression
            fold_change = np.random.choice([-1.5, 1.5])  # 1.5-fold change
            expression_data[idx, :n_als] *= fold_change
        
        # Create DataFrame
        df = pd.DataFrame(expression_data, index=genes, columns=all_samples)
        
        # Add group labels
        groups = ['ALS'] * n_als + ['Control'] * n_control
        
        return df, groups, genes
    
    def perform_de_analysis(self, expression_data, groups, dataset_id):
        """
        Perform differential expression analysis using t-tests.
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
            
            # Calculate fold change (log2)
            mean_als = np.mean(als_expression)
            mean_control = np.mean(control_expression)
            log2_fold_change = np.log2(mean_als / mean_control) if mean_control != 0 else 0
            
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
        
        return de_results
    
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
        
        # Meta-analysis: Combine p-values using Fisher's method
        meta_results = []
        
        # Get unique genes across all datasets
        all_genes = combined_df['gene'].unique()
        
        for gene in all_genes:
            gene_data = combined_df[combined_df['gene'] == gene]
            
            if len(gene_data) > 1:
                # Fisher's method for combining p-values
                chi_square = -2 * np.sum(np.log(gene_data['p_value']))
                combined_p = stats.chi2.sf(chi_square, 2 * len(gene_data))
                
                # Average fold change
                avg_log2_fc = gene_data['log2_fold_change'].mean()
                
                # Count number of datasets where significant
                n_sig = sum(gene_data['adj_p_value'] < 0.05)
                
            else:
                combined_p = gene_data['p_value'].iloc[0]
                avg_log2_fc = gene_data['log2_fold_change'].iloc[0]
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
    
    def save_results(self):
        """
        Save all results to Excel files.
        """
        print(f"\n{'='*60}")
        print("SAVING RESULTS")
        print(f"{'='*60}")
        
        # Create Excel writer
        with pd.ExcelWriter('ALS_Differential_Expression_Results.xlsx') as writer:
            
            # Save individual dataset results
            for dataset_id, result in self.results.items():
                sheet_name = f"{dataset_id}_DE"
                result['de_results'].to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"Saved {dataset_id} results")
            
            # Save combined results
            if self.combined_results is not None:
                self.combined_results['individual_results'].to_excel(
                    writer, sheet_name='All_Individual_Results', index=False
                )
                self.combined_results['meta_analysis'].to_excel(
                    writer, sheet_name='Meta_Analysis', index=False
                )
                self.combined_results['significant_in_multiple'].to_excel(
                    writer, sheet_name='Multi_Dataset_Significant', index=False
                )
                print("Saved combined results")
        
        print("All results saved to 'ALS_Differential_Expression_Results.xlsx'")
    
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
    
    # Combine results
    analyzer.combine_results()
    
    # Save results
    analyzer.save_results()
    
    # Create summary plot
    analyzer.create_summary_plot()
    
    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE!")
    print(f"{'='*60}")
    print("Results saved to:")
    print("- ALS_Differential_Expression_Results.xlsx")
    print("- ALS_DE_Analysis_Summary.png")

if __name__ == "__main__":
    main()