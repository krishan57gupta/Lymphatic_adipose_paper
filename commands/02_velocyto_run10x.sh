#!/usr/bin/env bash
set -euo pipefail

# Run once per Cell Ranger output directory.
# Replace all values in angle brackets.

VELOCYTO="<PATH_TO_VELOCYTO>/velocyto"
CELLRANGER_RUN="<PATH_TO_CELLRANGER_OUTPUT>/<RUN_ID>"
GTF="<PATH_TO_MM10_REFERENCE>/genes/genes.gtf"

"${VELOCYTO}" run10x "${CELLRANGER_RUN}" "${GTF}"

# Expected loom file:
# ${CELLRANGER_RUN}/velocyto/<RUN_ID>.loom
