#!/bin/bash

# Script to create project directory structure
# Author: Generated script
# Description: Creates a standardized directory structure for analysis project

echo "Creating project directory structure..."

# Create main directories and subdirectories
mkdir -p A_scripts                              # Analysis scripts
mkdir -p B_inputs                               # Input data files
mkdir -p B_inputs/EdgeR_results                 # Data, tables from EdgeR analysis
mkdir -p B_inputs/R_analysis                    # Plots, heatmap, correlation, tables from R analysis
mkdir -p C_intermediate/classification          # Genomic region classification
mkdir -p C_intermediate/filtered                # Filtered CleaveLand results
mkdir -p C_intermediate/function                # Functional annotation mapping
mkdir -p C_intermediate/mapping                 # miRNA-target mapping
mkdir -p C_intermediate/mirna_target            # miRNA-target pairs
mkdir -p D_pdfs_logs/logs                       # To save all generated log files
mkdir -p E_annotation                           # Annotation outputs
mkdir -p F_module_heatmaps                      # Module analysis visualizations
mkdir -p G_final_report                         # Final analysis reports
mkdir -p H_README                               # Documentation

echo "Directory structure created successfully!"
echo ""
echo "Created the following structure:"
echo "├── A_scripts/                              # Analysis scripts"
echo "├── B_inputs/                               # Input data files"
echo "│   ├── EdgeR_results/                      # Data, tables from EdgeR analysis"
echo "│   └── R_analysis/                         # Plots, heatmap, correlation, tables from R analysis"
echo "├── C_intermediate/                         # Intermediate processing files"
echo "│   ├── classification/                     # Genomic region classification"
echo "│   ├── filtered/                           # Filtered CleaveLand results"
echo "│   ├── function/                           # Functional annotation mapping"
echo "│   ├── mapping/                            # miRNA-target mapping"
echo "│   └── mirna_target/                       # miRNA-target pairs"
echo "├── D_pdfs_logs/                            # PDF outputs"
echo "│   └── logs/                               # To save all generated log files"
echo "├── E_annotation/                           # Annotation outputs"
echo "├── F_module_heatmaps/                      # Module analysis visualizations"
echo "├── G_final_report/                         # Final analysis reports"
echo "└── H_README/                               # Documentation"
echo ""
echo "You can now navigate to your directories and start your analysis!"
