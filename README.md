# Lymphatic_switch_paper

Code supporting the single-cell RNA-sequencing, RNA-velocity, Monocle, CellRank,
metabolite-mediated communication, GO enrichment pathways analysis, and transcription-factor motif analyses
for the lymphatic adipose study.

## Data availability

The raw FASTQ files and processed data are deposited under the private GEO
accession **GSE251702**. During peer review, reviewers should use the private
reviewer token supplied in the manuscript. The token is intentionally not
repeated in this public code repository.

## Analysis overview

```text
Private GEO GSE251702
        |
        v
Raw FASTQ files
        |
        |  Cell Ranger 7.1.0, refdata-gex-mm10-2020-A
        v
Cell Ranger output
  - filtered_feature_bc_matrix/
  - possorted_genome_bam.bam
        +-----------------------------------+
        |                                   |
        |                                   | velocyto 0.17.17, (mouse genes.gtf + repeat-mask GTF)
        |                                   v
        |                               loom files
        v                                   |
Seurat preprocessing                        |
  - QC and filtering                        |
  - normalization/integration               |
  - integration with batch correction       |
  - PCA/UMAP/clustering                     |
  - cell-type annotation                    |
  - processed Seurat objects                |
        |                                   |
        | export metadata, counts,          |                
        | PCA, genes and barcodes           |
        +------------------+----------------+
                           |
                           v
                    scVelo + CellRank
                           |
                           v
            RNA velocity and transition results

Processed Seurat expression/metadata
        |
        +--> MEBOCOST
        |
        +--> CellChat / DEG / GO enrichment/ Monocle /manuscript figures
        |
        +--> Statistical analysis

Configured genomic sequence + motifs
        |
        +--> motif scanning
```

## Repository contents

```text
Lymphatic_switch_paper/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── commands/
│   ├── 01_cellranger_count.sh
│   ├── 02_velocyto_run10x.sh
│   └── 03_run_r_pipeline.sh
├── config/
│   └── paths.example.yml
├── environment/
│   ├── environment.yml
│   ├── R_packages.txt
│   └── install_R_packages.R
└── scripts/
    ├── 00_preprocessing/
    │   ├── SCA_Intestine_count2processdData_Final.R
    │   ├── SCA_Intestine_dataAnalysis_Final.R
    │   └── SCA_Mesentery_combined_Final.R
    ├── 01_seurat_export/
    │   ├── RNA_Velocity_Intestine_Final.R
    │   └── RNA_Velocity_Mesentery_Final.R
    ├── 02_velocity_cellrank/
    │   ├── Velo_Rank_Intestine_Final.ipynb
    │   └── Velo_Rank_Mesentery_Final.ipynb
    ├── 03_mebocost/
    │   └── MEBOCOST_Intestine_Final.ipynb
    ├── 04_motif/
    │   └── bindingMotif_F.ipynb
    └── utils/
        └── functions.R
```

## Software requirements

The study used the following principal versions:

- Cell Ranger 7.1.0
- R 4.0.3
- Seurat 4.3.0
- SeuratWrappers 0.3.0
- Python 3.9.18
- Scanpy 1.8.2
- velocyto 0.17.17
- scVelo 0.3.1
- CellRank 2.0.7
- Monocle3 1.3.1
- MEBOCOST 1.2.0
- CellChat 2.0.7
- clusterProfiler 3.18.1
- org.Mm.eg.db 3.12.0
- Biopython 1.85
- pyfaidx 0.9.0.3
- JASPAR 2026 motif resources

The scripts also load plotting and data-manipulation packages listed in
`environment/R_packages.txt`.

## Hardware recommendations

Tested hardware configuration:

- macOS
- 8 CPU cores
- 64 GB RAM
- 1TB HD

RNA-velocity dynamical modeling is configured with `n_jobs=8` in the supplied
notebooks. Larger datasets may require more memory and temporary storage.

## Installation

### Cell Ranger and velocyto

Install Cell Ranger 7.1.0 according to the 10x Genomics instructions and make
the executable available on the command line. Install velocyto 0.17.17 in a
compatible Python environment.

Verify:

```bash
cellranger --version
velocyto --help
```

### Python environment

```bash
conda env create -f environment/environment.yml
conda activate lymphatic-adipose
python -c "import scanpy, scvelo, cellrank; print('Python environment ready')"
```

MEBOCOST may require installation according to the version-specific upstream
instructions if its dependencies are not resolved automatically.

### R environment

Start R 4.0.3 from the repository root:

```r
source("environment/install_R_packages.R")
sessionInfo()
```

Install the exact CellChat and Monocle3 versions used in the study

## Configuration before running

1. Copy the path template:

```bash
cp config/paths.example.yml config/paths.yml
```

2. Replace every placeholder in `config/paths.yml`.

3. Open each R script and notebook and replace its original absolute local
   paths (`mainDir`, `dataDir`, `processedDataDir`, `plotDir`, and related
   variables) with paths on the new system.

The supplied code represents the original analysis and therefore contains
machine-specific absolute paths. The YAML file documents the values that need
to be mapped, but the original scripts do not automatically read YAML.

4. In R scripts that use shared plotting/helper functions, set:

```r
source("scripts/utils/functions.R")
```

or adjust the existing `source()` path to this location.

## Complete execution order

### Step 1 — Obtain the private data

Download the raw FASTQ files from **GSE251702** using the reviewer token stated
in the manuscript. Also download any processed reference files needed to
compare or resume downstream analyses.

### Step 2 — Generate Cell Ranger count matrices

Run `commands/01_cellranger_count.sh` once for every sample after replacing the
placeholders. Cell Ranger produces:

```text
<RUN_ID>/outs/filtered_feature_bc_matrix/
<RUN_ID>/outs/possorted_genome_bam.bam
```

The filtered feature-barcode matrix is the input to the Seurat preprocessing scripts. The Cell Ranger BAM file and the filtered cell-barcode file (`filtered_feature_bc_matrix/barcodes.tsv.gz`) are used as inputs for velocyto.

### Step 3 — Generate velocyto loom files

Run `commands/02_velocyto_run10x.sh` once to process all Cell Ranger output directories listed in `SampleIDs.txt`.

Required inputs:
- Cell Ranger BAM (`outs/possorted_genome_bam.bam`)
- Filtered barcodes (`outs/filtered_feature_bc_matrix/barcodes.tsv.gz`)
- Gene annotation GTF (`<PATH_TO_CELLRANGER_REFERENCE>/genes/genes.gtf`)
- Repeat-mask GTF (`<PATH_TO_REPEAT_MASK_GTF>/mm10_rmsk.gtf`)

Load `velocyto` (v0.17.17) and `samtools` using your HPC environment. The script processes all samples listed in `SampleIDs.txt` and generates loom files for downstream scVelo and CellRank analyses.

### Step 4 — Run intestine Seurat preprocessing

```bash
Rscript scripts/00_preprocessing/SCA_Intestine_count2processdData_Final.R
```

This script reads Cell Ranger count matrices, creates Seurat objects, performs
quality control, normalization, integration, PCA, UMAP, neighbor graph
construction, clustering, annotation, and saves the processed intestine Seurat
object.

The implemented analysis includes retention of cells with 200–2500 detected
genes, exclusion of cells with more than 30% mitochondrial reads, PCA/UMAP on
30 dimensions, and clustering at the configured resolution. Confirm all
parameters against the final manuscript and script before publication.

### Step 5 — Run intestine downstream analyses and figures

```bash
Rscript scripts/00_preprocessing/SCA_Intestine_dataAnalysis_Final.R
```

This script loads the processed intestine object and performs downstream
analyses and figure generation. It uses helper functions from
`scripts/utils/functions.R`; update the `source()` path accordingly.

### Step 6 — Run mesentery processing and analyses

```bash
Rscript scripts/00_preprocessing/SCA_Mesentery_combined_Final.R
```

This combined script performs mesentery preprocessing and downstream analyses
from Cell Ranger count matrices through manuscript outputs.

### Step 7 — Export Seurat information for Python

Run:

```bash
Rscript scripts/01_seurat_export/RNA_Velocity_Intestine_Final.R
Rscript scripts/01_seurat_export/RNA_Velocity_Mesentery_Final.R
```

These bridge scripts export the information required to connect Seurat with
the Python RNA-velocity notebooks, including:

- cell metadata (`metadataF.csv`)
- expression/count information
- PCA coordinates
- gene names (`gene_namesF.csv`)
- cell barcodes

Keep the exported cell barcodes consistent with the loom-file barcode format.
The Python notebooks subset and merge these exported cells with the spliced and
unspliced matrices from velocyto.

### Step 8 — Run scVelo and CellRank

Start Jupyter:

```bash
conda activate lymphatic-adipose
jupyter lab
```

Run all cells in this order:

1. `scripts/02_velocity_cellrank/Velo_Rank_Intestine_Final.ipynb`
2. `scripts/02_velocity_cellrank/Velo_Rank_Mesentery_Final.ipynb`

Before running, update:

- the project/data directories
- the four intestine loom filenames
- the two mesentery loom filenames
- processed-output and plot directories

The notebooks:

1. reconstruct an AnnData object from Seurat-exported counts, metadata, genes,
   barcodes, and PCA coordinates;
2. read and merge the matching loom files;
3. filter and normalize spliced/unspliced counts;
4. calculate moments;
5. recover dynamical parameters;
6. estimate dynamical RNA velocity;
7. construct CellRank velocity and connectivity kernels;
8. calculate transition/fate results and export plots or numerical values.

Important notebook parameters include:

```python
scv.pp.filter_and_normalize(
    adata,
    min_shared_counts=20,
    n_top_genes=2000,
    subset_highly_variable=True
)
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
```

Some CellRank probabilities used in publication-quality R plots were copied
from the computed notebook outputs.

### Step 9 — Run MEBOCOST

Run all cells in:

```text
scripts/03_mebocost/MEBOCOST_Intestine_Final.ipynb
```

This notebook uses the Seurat-derived expression matrix and metadata to create
an AnnData object and perform metabolite-mediated cell–cell communication
analysis with MEBOCOST 1.2.0. Update its input, output, and reference-database
paths before execution.

Numerical communication scores used for customized R plotting should be saved
as source-data CSV files to preserve traceability.

### Step 10 — Run motif analysis

Run all cells in:

```text
scripts/04_motif/bindingMotif_F.ipynb
```

Configure:

- the mouse mm10 reference FASTA
- target promoter or genomic coordinates
- JASPAR 2026 motif inputs
- the output directory

The notebook uses pyfaidx and Biopython motif utilities to retrieve sequence
and scan transcription-factor binding motifs.

## Input–output connection table

| Stage | File | Main input | Main output / next stage |
|---|---|---|---|
| Cell Ranger | `commands/01_cellranger_count.sh` | GEO FASTQ and refdata-gex-mm10-2020-A  | filtered count matrix and BAM |
| velocyto | `commands/02_velocyto_run10x.sh` | Cell Ranger BAM, filtered cell barcodes, genes.gtf and repeat-mask GTF | loom files |
| Intestine preprocessing | `SCA_Intestine_count2processdData_Final.R` | intestine count matrices | processed intestine Seurat RDS |
| Intestine analysis | `SCA_Intestine_dataAnalysis_Final.R` | processed intestine Seurat RDS | downstream results and figures |
| Mesentery pipeline | `SCA_Mesentery_combined_Final.R` | mesentery count matrices | processed object, analyses and figures |
| Intestine export | `RNA_Velocity_Intestine_Final.R` | processed intestine Seurat RDS | metadata/count/PCA/gene/barcode files |
| Mesentery export | `RNA_Velocity_Mesentery_Final.R` | processed mesentery Seurat RDS | metadata/count/PCA/gene/barcode files |
| Intestine velocity/rank | `Velo_Rank_Intestine_Final.ipynb` | intestine exports and four loom files | H5AD, velocity and CellRank results |
| Mesentery velocity/rank | `Velo_Rank_Mesentery_Final.ipynb` | mesentery exports and two loom files | H5AD, velocity and CellRank results |
| Metabolic communication | `MEBOCOST_Intestine_Final.ipynb` | expression and metadata exports | MEBOCOST scores and plots |
| Motif analysis | `bindingMotif_F.ipynb` | FASTA, coordinates and motifs | motif matches and plots |

## Reproducibility notes

- The repository intentionally excludes all biological data.
- Reviewer access is through private GEO accession GSE251702.
- Never place the reviewer token in a public GitHub commit.
- Replace all machine-specific paths before running.
- Preserve sample names and barcode transformations exactly across Seurat and
  loom files.
- Run the scripts in the documented order.
- Record `sessionInfo()` and `conda env export` after final validation.
- Export CellRank and MEBOCOST numerical outputs to source-data tables when
  they are replotted in R.
- Cosmetic figure editing may change labels, fonts, dimensions, or colors, but
  must not alter the underlying values or statistical results.

## License

BSD 3-Clause License. See `LICENSE`.
