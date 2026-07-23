#!/usr/bin/env bash
set -euo pipefail
find processed_data source_data metadata -type f ! -name README.md -print0 | sort -z | xargs -0 shasum -a 256 > checksums.sha256
echo "Wrote checksums.sha256"
