# TCGA BRCA Gene Expression Matrix Generator  

This repository provides a Python pipeline to process raw TCGA RNA-seq (BRCA) transcriptome profiling data into a clean, unified gene expression matrix (samples × genes).  

The workflow automatically filters for protein-coding genes, maps Ensembl IDs to HGNC gene symbols, and outputs an ML-ready CSV for downstream applications such as cancer classification, biomarker discovery, or survival analysis.  

---

## Features  
- Recursively scans and processes TCGA RNA-seq TSV files  
- Keeps only protein-coding genes (`gene_type == protein_coding`)  
- Maps Ensembl IDs → HGNC gene names (`gene_name` column)  
- Removes QC rows (e.g., `N_` genes)  
- Cleans Ensembl IDs (removes version suffix like `.15`)  
- Merges all samples into a samples × genes matrix  
- Fills missing values with `0` for consistency  
- Exports to CSV with validation and summary statistics  

---



