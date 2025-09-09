import os
import pandas as pd
import numpy as np
from pathlib import Path

def process_single_tsv(filepath):
   
    try:
        # Read TSV file, skipping comment lines that start with '#'
        df = pd.read_csv(filepath, sep='\t', comment='#')
        
        # Filter to keep only protein-coding genes
        df_filtered = df[df['gene_type'] == 'protein_coding'].copy()
        
        # Remove QC rows (genes starting with 'N_')
        df_filtered = df_filtered[~df_filtered['gene_name'].str.startswith('N_')]
        
        # Remove version numbers from Ensembl IDs 
        df_filtered['gene_id_clean'] = df_filtered['gene_id'].str.split('.').str[0]
        
        # Extract only the columns we need: gene_name and tpm_unstranded
        result_df = df_filtered[['gene_name', 'tpm_unstranded']].copy()
        
        return result_df
        
    except Exception as e:
        print(f"Error processing {filepath}: {str(e)}")
        return None

def create_tcga_expression_matrix(root_folder, output_csv='tcga_brca_expression_matrix.csv'):
    
    print("🔍 Searching for TSV files...")
    
    # Dictionary to store all sample data
    sample_data = {}
    tsv_files_found = 0
    failed_files = 0
    
    # Recursively walk through all directories
    for root, dirs, files in os.walk(root_folder):
        for filename in files:
            # Look for RNA-seq TSV files
            if filename.endswith('.tsv') and 'rna_seq' in filename:
                filepath = os.path.join(root, filename)
                
                # Use the parent directory name as sample ID (UUID)
                sample_id = os.path.basename(root)
                
                print(f"Processing: {sample_id}")
                
                # Process the TSV file
                processed_data = process_single_tsv(filepath)
                
                if processed_data is not None and not processed_data.empty:
                    # Set gene_name as index and extract TPM values
                    gene_expression = processed_data.set_index('gene_name')['tpm_unstranded']
                    sample_data[sample_id] = gene_expression
                    tsv_files_found += 1
                else:
                    failed_files += 1
    
    print(f"\n Processing Summary:")
    print(f"   • TSV files found and processed: {tsv_files_found}")
    print(f"   • Failed files: {failed_files}")
    
    if tsv_files_found == 0:
        print("No valid TSV files found!")
        return None
    
    print("\n Merging all samples into expression matrix...")
    
    # Create the unified expression matrix
    # Rows = samples, Columns = genes
    expression_matrix = pd.DataFrame(sample_data).T
    
    # Fill any missing values with 0 (in case some genes are missing in some samples)
    expression_matrix = expression_matrix.fillna(0)
    
    # Sort columns (genes) alphabetically for consistency
    expression_matrix = expression_matrix.sort_index(axis=1)
    
    print(f"\n Saving to: {output_csv}")
    
    # Save to CSV
    expression_matrix.to_csv(output_csv, index_label='sample_id')
    
    print(f"\n Success! Final matrix shape: {expression_matrix.shape}")
    print(f"   • Samples (rows): {expression_matrix.shape[0]}")
    print(f"   • Genes (columns): {expression_matrix.shape[1]}")
    
    # Display some basic statistics
    print(f"\n Data Statistics:")
    print(f"   • Mean expression across all samples: {expression_matrix.values.mean():.2f}")
    print(f"   • Max expression value: {expression_matrix.values.max():.2f}")
    print(f"   • Non-zero values: {(expression_matrix.values > 0).sum():,} / {expression_matrix.size:,}")
    
    return expression_matrix

def verify_output(csv_path):
    """
    Verify the generated CSV file and display basic information.
    
    Args:
        csv_path (str): Path to the generated CSV file
    """
    if not os.path.exists(csv_path):
        print(f" File not found: {csv_path}")
        return
    
    print(f"\n🔍 Verifying output file: {csv_path}")
    
    # Read the CSV and display info
    df = pd.read_csv(csv_path, index_col=0)
    print(f"   • File size: {os.path.getsize(csv_path) / (1024*1024):.1f} MB")
    print(f"   • Matrix shape: {df.shape}")
    print(f"   • Sample IDs (first 5): {list(df.index[:5])}")
    print(f"   • Gene names (first 5): {list(df.columns[:5])}")
    
    return df

# Main execution
if __name__ == "__main__":
    # YOUR CONFIGURATION
    ROOT_FOLDER = "cancer_brca"
    OUTPUT_FILE = "tcga_brca_gene_expression.csv"
    
    print(" TCGA BRCA Gene Expression Matrix Generator")
    print("=" * 50)
    
    # Check if input folder exists
    if not os.path.exists(ROOT_FOLDER):
        print(f" Error: Folder '{ROOT_FOLDER}' not found!")
        print("Please make sure your GDC download folder exists at this path.")
    else:
        # Process the data
        matrix = create_tcga_expression_matrix(ROOT_FOLDER, OUTPUT_FILE)
        
        if matrix is not None:
            # Verify the output
            verify_output(OUTPUT_FILE)
            
            print(f"\n Complete! Your ML-ready dataset is saved as: {OUTPUT_FILE}")
            print("Ready for use in pandas, scikit-learn, XGBoost, and neural networks!")

# Quick usage function for interactive use
def quick_process():
    """Quick function to process your TCGA data"""
    return create_tcga_expression_matrix("/home/dddd/gdc_downloads/", "tcga_brca_expression_matrix.csv")
