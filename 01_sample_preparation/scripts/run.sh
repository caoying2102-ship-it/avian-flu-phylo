#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

python3 "$SCRIPT_DIR/prepare_samples.py" \
  --input-dir "$STAGE_DIR/input" \
  --output-dir "$STAGE_DIR/output"
