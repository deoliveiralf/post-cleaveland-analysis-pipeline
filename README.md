# Post-Cleaveland Analysis Pipeline & Tools

**Complete guideline and tools for miRNA-target post-CleaveLand analysis, table generation, and visualization.**

---

## Overview

This repository provides a modular, reproducible pipeline for analysis of CleaveLand output, annotation enhancement, expression and correlation analysis, and generation of interactive HTML reports and publication-ready plots for plant miRNA research. The tools are designed for researchers working with *Setaria viridis* and similar plant genomes.

**Developed by [Leandro F. de Oliveira](https://github.com/deoliveiralf) | Supported by [LigninLab](https://sites.usp.br/ligninlab/) | University of São Paulo | Year: 2025**

---

## Pipeline Flow

**Steps:**

1. **Filter CleaveLand Results & Extract PDFs**  
   Filter CleaveLand results by score, MFE, p-value, etc. and extract matching PDFs.

2. **Extract miRNA-Target Pairs**  
   Identify miRNA-target pairs from filtered results.

3. **Region Analysis**  
   Annotate cleavage sites and classify by transcript region (5'UTR/CDS/3'UTR).

4. **Rename miRNA IDs**  
   Standardize miRNA IDs using reference FASTA.

5. **Add miRNA Info**  
   Enrich tables with miRNA type and full ID from FASTA.

6. **Rename miRNAs in CPM Table**  
   Map CPM table IDs to standardized names.

7. **Rename miRNA Modules Table**  
   Standardize module table IDs.

8–10. **Expression & Correlation Analysis (in R)**  
   Generate heatmaps and correlation plots for miRNAs and targets.

11. **Merge Tables for Final Report**  
   Integrate all results into interactive HTML, TSV, JSON.

12. **Module Heatmap & Correlation Analysis**  
   Summarize modules, top correlations, and create visual outputs.

13. **Plot Cleavage Sites**  
   Make publication-ready PNG plots showing cleavage locations.

14. **Combine Result Plots**  
   Merge heatmaps and correlation plots into composite figures.

**See the [interactive pipeline HTML](Interactive Post-Cleaveland Analysis Pipeline Overview_v5.html) for a full diagram and details.**

---

## Project Directory Structure

Create the recommended structure with:
```bash
bash create_dirs.sh
```
**Folders:**
- `A_scripts/` - Analysis scripts
- `B_inputs/` - Input data files (CleaveLand outputs, FASTA, GFF3, CPM/TPM tables, etc.)
- `C_intermediate/` - Intermediate files (filtered, mapping, region classification)
- `D_pdfs_logs/` - PDF outputs and logs
- `E_annotation/` - Annotation results
- `F_module_heatmaps/` - Module analysis visualizations
- `G_final_report/` - Final tables and reports
- `H_README/` - Documentation

---

## Quick Start Guide

### Prerequisites

- Python 3.7+
- R (with required packages: `edgeR`, `ComplexHeatmap`, `ggplot2`, etc.)
- Required Python packages: `pandas`, `numpy`, `matplotlib`, `seaborn`
- CleaveLand4 output files
- Genome annotation files (GFF3)
- Expression tables for transcripts and miRNAs
- miRNA identification files (FASTA, list, etc.)

### Essential Steps

1. **Filter CleaveLand results**  
   ```bash
   python 01_filter-cleaveland-results.py full_results.txt --pdf_dir D_pdfs_logs/TPlot.pdfs/ --output_dir C_intermediate/filtered/ --pdf_output_dir D_pdfs_logs/matched_pdfs/
   ```
2. **Extract miRNA-target pairs**  
   ```bash
   python 02_mirna-target-modules.py C_intermediate/filtered/filtered_cleaveland_results.txt > C_intermediate/mirna_target/mirna-target-modules-table.txt
   ```
3. **Annotate cleavage regions**  
   ```bash
   python 03_annotate_and_plot_cleavage_regions.py B_inputs/cleaveland_output.txt C_intermediate/mirna_target/mirna-target-modules-table.txt B_inputs/Sviridis_726_v4.1.gene.gff3 C_intermediate/classification/filtered_annotated_cleavages.csv
   ```
4. **Continue with steps 4–14 as described in the HTML guide.**

### R Scripts

- `08_mirna_heatmap.R` – Generate individual miRNA heatmaps
- `09_target_heatmap.R` – Target gene heatmaps
- `10_mirna_correlation.r` – Correlation analysis

---

## Usage

- **Run each script in order (see the pipeline HTML for command examples and parameter explanations).**
- Use provided bash, Python, and R scripts for each step.
- Outputs will be generated in the corresponding folders.
- Use the merging and integration scripts to produce interactive HTML reports for downstream viewing.

---

## Data Format Requirements

- **CleaveLand outputs:** Standard full_results.txt format
- **FASTA files:** miRNA IDs in the pattern `miR166_Chr09_45260`
- **Expression tables:** Tab-separated values, gene IDs as first column
- **Annotation files:** GFF3 format, gene annotation TSVs
- **Module analysis PNGs:** Naming convention `[miRNA_Chr ID]_[Target ID]_combined_analysis.png`

---

## Visualization

- **Interactive HTML output:**  
  Final tables can be viewed in browser with filtering and data export.
- **PNG plots:**  
  Publication-ready visualizations of heatmaps, correlations, and cleavage sites.

---

## Troubleshooting

- **Missing IDs:** Ensure all IDs are standardized using the renaming scripts.
- **No matches in annotation:** Check gene ID formats and remove variant numbers if necessary.
- **Plotting errors:** Verify input file paths and required columns.

---

## Support & Credits

Developed by [Leandro F. de Oliveira](https://github.com/deoliveiralf)  
Supported by [LigninLab](https://sites.usp.br/ligninlab/), FAPESP, CNPq, Biocel Lab, University of São Paulo.

Contact: `lfdeoliveira@usp.br`

---

## License

Academic use only. Contact the author for licensing details.

---

## References

See the [pipeline HTML overview](Interactive Post-Cleaveland Analysis Pipeline Overview_v5.html) for detailed documentation, command-line examples, and step-by-step instructions.
