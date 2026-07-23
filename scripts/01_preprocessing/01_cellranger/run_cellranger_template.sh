#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../../../config/project_config.sh" 2>/dev/null || true

# Complete one row per sample. Keep sample identifiers consistent with metadata.
# cellranger count \
#   --id=SAMPLE_ID \
#   --transcriptome="$CELLRANGER_TRANSCRIPTOME" \
#   --fastqs=/path/to/fastqs \
#   --sample=SAMPLE_PREFIX \
#   --localcores=16 \
#   --localmem=64

echo "Template only: add the exact sample-specific Cell Ranger commands used in the study."
