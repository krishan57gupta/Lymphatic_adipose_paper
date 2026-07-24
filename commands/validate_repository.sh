#!/usr/bin/env bash
set -euo pipefail
required=(
  README.md LICENSE CITATION.cff
  scripts/00_preprocessing/SCA_Intestine_count2processdData_Final.R
  scripts/00_preprocessing/SCA_Intestine_dataAnalysis_Final.R
  scripts/00_preprocessing/SCA_Mesentery_combined_Final.R
  scripts/01_seurat_export/RNA_Velocity_Intestine_Final.R
  scripts/01_seurat_export/RNA_Velocity_Mesentery_Final.R
  scripts/02_velocity_cellrank/Velo_Rank_Intestine_Final.ipynb
  scripts/02_velocity_cellrank/Velo_Rank_Mesentery_Final.ipynb
  scripts/03_mebocost/MEBOCOST_Intestine_Final.ipynb
  scripts/04_motif/bindingMotif_F.ipynb
  scripts/utils/functions.R
)
for f in "${required[@]}"; do
  [[ -f "$f" ]] || { echo "Missing: $f"; exit 1; }
done
echo "Repository structure is complete."
