#!/usr/bin/env bash
set -euo pipefail

# IMPORTANT:
# First edit all directory variables at the beginning of each R script.
# Run from the repository root.

Rscript scripts/00_preprocessing/SCA_Intestine_count2processdData_Final.R
Rscript scripts/00_preprocessing/SCA_Intestine_dataAnalysis_Final.R
Rscript scripts/00_preprocessing/SCA_Mesentery_combined_Final.R

# These scripts export Seurat counts, metadata, PCA coordinates, gene names,
# and barcodes needed by the Python notebooks.
Rscript scripts/01_seurat_export/RNA_Velocity_Intestine_Final.R
Rscript scripts/01_seurat_export/RNA_Velocity_Mesentery_Final.R
