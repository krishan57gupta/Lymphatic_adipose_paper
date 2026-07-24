if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_packages <- c("Seurat", "SeuratWrappers", "reticulate", "Matrix", "ggplot2",
                   "ggrastr", "ggpubr", "export", "Hmisc", "stringr", "MASS",
                   "scales", "ggbreak", "cowplot", "dplyr", "tidyr",
                   "tidygraph", "ggraph", "patchwork", "ggridges", "extrafont")
bioc_packages <- c("clusterProfiler", "org.Mm.eg.db")

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)

missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc)

message("CellChat may require installation from its official repository/version used by the study.")
sessionInfo()
