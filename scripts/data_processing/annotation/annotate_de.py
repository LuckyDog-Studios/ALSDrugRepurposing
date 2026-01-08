import pandas as pd
from pathlib import Path

def extract_gene_symbol(gene_assignment):
    """
    Extract gene symbol from the gene_assignment field.
    """
    if pd.isna(gene_assignment) or gene_assignment == '' or gene_assignment == '---':
        return ''
    
    assignments = str(gene_assignment).split(' /// ')
    gene_symbols = set()
    
    for assignment in assignments:
        fields = [field.strip() for field in assignment.split(' // ')]
        if len(fields) >= 2:
            gene_symbol = fields[1]
            if gene_symbol and gene_symbol != '---':
                gene_symbols.add(gene_symbol)
    
    return ';'.join(sorted(gene_symbols)) if gene_symbols else ''

def annotate_expression_data(expression_file, annotation_file):
    """
    Annotate differential expression data with Gene Symbols.
    """
    
    # Read the expression data
    print(f"Reading expression data from: {expression_file}")
    expr_df = pd.read_csv(expression_file)
    
    # Read the annotation file
    print(f"Reading annotation data from: {annotation_file}")
    annot_df = pd.read_csv(annotation_file, sep='\t', comment='#')
    
    # Extract probe IDs from expression data (first column)
    probe_col = expr_df.columns[0]
    expr_probe_ids = expr_df[probe_col].astype(str).str.replace('"', '').tolist()
    
    # Create mapping from probe ID to gene assignment
    probe_to_gene_assignment = dict(zip(
        annot_df['ID'].astype(str), 
        annot_df['gene_assignment']
    ))
    
    # Map probe IDs to gene symbols
    gene_symbols = []
    for probe_id in expr_probe_ids:
        if probe_id in probe_to_gene_assignment:
            gene_symbols.append(extract_gene_symbol(probe_to_gene_assignment[probe_id]))
        else:
            # Try without suffix
            base_id = probe_id.split('.')[0] if '.' in probe_id else probe_id
            gene_symbols.append(
                extract_gene_symbol(probe_to_gene_assignment.get(base_id, ''))
            )
    
    # Add Gene Symbol column and drop the original probe ID column
    expr_df['Gene_Symbol'] = gene_symbols
    expr_df = expr_df.drop(columns=[probe_col])
    
    # Reorder columns to have Gene_Symbol first
    cols = expr_df.columns.tolist()
    cols.remove('Gene_Symbol')
    expr_df = expr_df[['Gene_Symbol'] + cols]
    
    # Filter out rows without gene symbols (optional - remove if you want to keep all)
    original_count = len(expr_df)
    expr_df = expr_df[expr_df['Gene_Symbol'] != '']
    filtered_count = len(expr_df)
    
    print("\nAnnotation summary:")
    print(f"Total rows: {original_count}")
    print(f"Rows with gene symbols: {filtered_count}")
    print(f"Rows without gene symbols removed: {original_count - filtered_count}")
    
    return expr_df

def main():
    # Hard-coded file paths
    EXPRESSION_FILE = "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE118336\\GSE118336_DE_FUS_Heterozygous_vs_FUS_H517D_Mutant.csv" # Your DE data
    ANNOTATION_FILE = "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\datasets\\annotation\\GPL17586-45144.txt" # The annotation file
    
    # Use command line arguments if provided
    import sys
    if len(sys.argv) > 2:
        EXPRESSION_FILE = sys.argv[1]
        ANNOTATION_FILE = sys.argv[2]
    
    # Check files exist
    for file_path in [EXPRESSION_FILE, ANNOTATION_FILE]:
        if not Path(file_path).exists():
            print(f"Error: File not found: {file_path}")
            print("Usage: python annotate_expression.py <expression_file> <annotation_file>")
            return
    
    # Annotate and save
    annotated_df = annotate_expression_data(EXPRESSION_FILE, ANNOTATION_FILE)
    output_file = "GSE118336_annotated_de_FUS_Heterozygous_vs_FUS_H517D_Mutant.csv"
    annotated_df.to_csv(output_file, index=False)
    print(f"\nSaved to: {output_file}")
    print("\nFirst 5 rows:")
    print(annotated_df.head())

if __name__ == "__main__":
    main()