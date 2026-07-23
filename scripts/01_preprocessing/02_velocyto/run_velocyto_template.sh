#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../../../config/project_config.sh" 2>/dev/null || true

# Example for each Cell Ranger sample directory:
# velocyto run10x \
#   -m "$VELOCYTO_MASK_GTF" \
#   /path/to/cellranger/SAMPLE_ID \
#   "$VELOCYTO_GTF"

echo "Template only: add the six exact velocyto commands used to generate four intestine and two mesentery loom files."
