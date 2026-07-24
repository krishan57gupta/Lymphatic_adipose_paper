#!/usr/bin/env bash
set -euo pipefail

# Run once per sample after downloading FASTQ files from private GEO accession GSE251702.
# Replace all values in angle brackets.

CELLRANGER="<PATH_TO_CELLRANGER>/cellranger"
TRANSCRIPTOME="<PATH_TO_MM10_REFERENCE>"
FASTQ_DIR="<PATH_TO_SAMPLE_FASTQ_DIRECTORY>"
SAMPLE_PREFIX="<FASTQ_SAMPLE_PREFIX>"
RUN_ID="<UNIQUE_SAMPLE_RUN_ID>"
OUTPUT_ROOT="<PATH_TO_CELLRANGER_OUTPUT_ROOT>"

mkdir -p "${OUTPUT_ROOT}"
cd "${OUTPUT_ROOT}"

"${CELLRANGER}" count \
  --id="${RUN_ID}" \
  --transcriptome="${TRANSCRIPTOME}" \
  --fastqs="${FASTQ_DIR}" \
  --sample="${SAMPLE_PREFIX}" \
  --localcores=8 \
  --localmem=64

# Main downstream input:
# ${OUTPUT_ROOT}/${RUN_ID}/outs/filtered_feature_bc_matrix/
# BAM input for velocyto:
# ${OUTPUT_ROOT}/${RUN_ID}/outs/possorted_genome_bam.bam
