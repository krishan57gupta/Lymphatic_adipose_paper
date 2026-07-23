# FOXC2 lymphatic endothelial cell single-cell analysis

This repository contains the analysis code and reproducibility documentation for the intestine and mesentery single-cell RNA-sequencing analyses described in the associated manuscript. It is structured to address the Nature Portfolio **Code and Software Submission Checklist**.

> **Important:** Large data files are not included in this ZIP. Add the two processed Seurat objects, six loom files, source-data tables, and any controlled-access/raw-data links before submission. See `DATA_MANIFEST.tsv` and the placeholder README files inside the data folders.

## Workflow overview

```text
FASTQ files
   |-- Cell Ranger v7.1.0 --> filtered count matrices + BAM files
   |                              |
   |                              |-- Seurat preprocessing/QC/integration/clustering
   |                              |       |-- processed Seurat objects
   |                              |
   |                              |-- velocyto v0.17.17 --> six loom files
   |                                                     |
   |                                                     |-- scVelo + CellRank
   |
   |-- downstream R/Python analyses --> manuscript figures and source data
```

### Study-specific organization

- **Intestine:** two biological datasets/replicate batches, producing four loom files (WT and LOF for each dataset).
- **Mesentery:** one dataset, producing two loom files (WT and LOF).
- `SCA_Intestine_count2processdData_Final.R` starts from Cell Ranger count matrices and generates the processed intestine Seurat object.
- `SCA_Intestine_dataAnalysis_Final.R` starts from that processed object and performs downstream R analyses and figure generation.
- `SCA_Mesentery_combined_Final.R` performs mesentery preprocessing and downstream R analyses.
- Separate notebooks/scripts perform RNA velocity, CellRank, MEBOCOST, and motif scanning.

## Repository layout

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── CODE_AVAILABILITY.md
├── DATA_AVAILABILITY.md
├── DATA_MANIFEST.tsv
├── FIGURE_SCRIPT_MAP.tsv
├── MANUAL_VALUES_LOG.tsv
├── config/
├── environment/
├── data/
├── metadata/
├── processed_data/
├── scripts/
├── source_data/
├── expected_output/
├── figures/
├── supplementary/
└── docs/
```

## System requirements

### Tested software versions

| Component | Version |
|---|---:|
| Cell Ranger | 7.1.0 |
| velocyto | 0.17.17 |
| R | 4.3.0 |
| Python | 3.9.18 |
| GraphPad Prism | 10 |
| Seurat | 4.3.0 |
| SeuratWrappers | 0.4.0 |
| Scanpy | 1.8.2 |
| scVelo | 0.3.1 |
| Monocle3 | 1.3.1 |
| CellRank | 2.0.7 |
| MEBOCOST | 1.2.0 |
| CellChat | 2.0.7 |
| clusterProfiler | 3.18.1 |
| org.Mm.eg.db | 3.18.0 |
| Biopython | 1.85 |
| pyfaidx | 0.9.0.3 |
| JASPAR database | 2026 release |

### Operating system and hardware

Fill in the exact tested systems before submission:

- Operating system: `[Linux distribution/version or macOS version]`
- CPU: `[model and number of cores]`
- RAM: `[GB]`
- Storage required for the complete workflow: `[GB/TB]`
- GPU: not required unless one was used in your actual workflow.

Cell Ranger, velocyto, and full single-cell analyses can require substantially more memory and storage than the reduced demonstration workflow.

## Installation

### Python environment

```bash
conda env create -f environment/environment_velocity.yml
conda activate foxc2_velocity
python -m pip install -r environment/requirements_motif.txt
```

Record the actual installation time on the tested machine in `docs/TESTED_SYSTEMS.md`.

### R environment

```r
install.packages("renv")
renv::restore()
```

The supplied `renv.lock.template` is a starting template. Replace it with a lockfile generated from the exact analysis environment:

```r
renv::snapshot()
writeLines(capture.output(sessionInfo()), "environment/sessionInfo.txt")
```

## Required input files

### Raw and Cell Ranger outputs

For each sample, the Cell Ranger output directory should contain at least:

```text
outs/filtered_feature_bc_matrix/
outs/possorted_genome_bam.bam
```

Raw FASTQ files do not need to be duplicated in GitHub when deposited in a suitable public repository. Add accession numbers and links to `DATA_AVAILABILITY.md`.

### Processed objects and loom files

Place or externally archive:

```text
processed_data/intestine/combined_Intestine_ABC_F.rds
processed_data/intestine/loom/<four intestine loom files>
processed_data/mesentery/combined_Mesentery_ABC_F.rds
processed_data/mesentery/loom/<two mesentery loom files>
```

If these files exceed GitHub's limits, deposit them in Zenodo, Figshare, Dryad, or another appropriate repository and list the persistent links and checksums in `DATA_MANIFEST.tsv`.

## Running the analysis

### 1. Cell Ranger

Edit `config/project_config.sh`, then run the commands in:

```bash
bash scripts/01_preprocessing/01_cellranger/run_cellranger_template.sh
```

The file is a template because sample identifiers, FASTQ paths, and transcriptome reference paths must match your deposited data.

### 2. velocyto

Use:

```bash
bash scripts/01_preprocessing/02_velocyto/run_velocyto_template.sh
```

Expected output: one `.loom` file per sample/dataset combination.

### 3. Intestine preprocessing

```bash
Rscript scripts/01_preprocessing/03_intestine/SCA_Intestine_count2processdData_Final.R
```

Expected output: processed intestine Seurat object plus QC and clustering outputs. Before running, replace the placeholder project root in the script or use the path instructions in `config/README.md`.

### 4. Intestine downstream analysis

```bash
Rscript scripts/02_downstream/intestine/SCA_Intestine_dataAnalysis_Final.R
```

### 5. Mesentery analysis

```bash
Rscript scripts/01_preprocessing/04_mesentery/SCA_Mesentery_combined_Final.R
```

### 6. RNA velocity, CellRank, MEBOCOST, and motif analysis

Open the corresponding notebooks in JupyterLab:

```bash
jupyter lab
```

Execute notebooks in their saved order after updating input/output paths according to `config/README.md`.

## Demonstration

Nature requests a small real or simulated dataset, expected output, and run-time information. The directory `data/demo/` contains instructions and placeholders. Because no small demo dataset was provided in the current upload, this repository includes a **demo specification rather than fabricated biological data**.

Before peer review, add either:

1. a small subset of real cells/genes that can be shared, or
2. a simulated count matrix and metadata that exercise the code without being interpreted biologically.

Document the expected outputs and measured run time in `expected_output/demo/README.md` and `docs/TESTED_SYSTEMS.md`.

## Values transferred between programs

Some Python-derived values (for example, CellRank fate probabilities or RNA-velocity summaries) and MEBOCOST-derived values were recorded and then used in R to create publication-quality plots. This is acceptable only when the transfer is fully traceable.

Use `MANUAL_VALUES_LOG.tsv` to record:

- source notebook and cell/output location;
- exact numerical values or exported table;
- destination R script and figure;
- whether values were copied manually or imported from a file;
- verification performed.

**Best practice:** export these values directly to CSV/TSV from Python and read them in R. Do not rely only on manually typed numbers.

## Figure reproduction

`FIGURE_SCRIPT_MAP.tsv` maps manuscript figures to scripts, inputs, and outputs. Complete all rows before submission. Each quantitative figure should have an accompanying source-data table in `source_data/`.

## Manual figure editing

Document unavoidable cosmetic edits in `docs/FIGURE_EDITING_NOTES.md`. Edits must not change data values, statistical results, point locations, or scientific interpretation.

## Reproducibility checklist

Before submission:

- [ ] Replace all placeholder paths.
- [ ] Add the two processed Seurat objects or persistent download links.
- [ ] Add the six loom files or persistent download links.
- [ ] Add raw-data accession numbers.
- [ ] Add sample metadata and cell annotations.
- [ ] Export Python/MEBOCOST-derived plot values as machine-readable tables.
- [ ] Complete the figure-to-script map.
- [ ] Add source data for every quantitative figure.
- [ ] Generate `renv.lock`, `sessionInfo.txt`, and a final Conda YAML from the actual environments.
- [ ] Record exact operating system, hardware, installation time, and run times.
- [ ] Add a small demo dataset and expected outputs.
- [ ] Run the repository on a clean machine or ask an unfamiliar colleague to test it.
- [ ] Create a tagged release and archive it in a DOI-minting repository.

## License

The repository currently includes a BSD 3-Clause License template. Confirm institutional ownership and licensing requirements before public release.

## Citation

Update `CITATION.cff` with the final manuscript title, author list, repository URL, release version, and DOI.
