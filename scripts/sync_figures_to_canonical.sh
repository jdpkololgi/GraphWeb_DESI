#!/bin/bash
# Mirror finished figures from a run-specific output dir into the canonical,
# conference-agnostic figure root (shared.config_paths.GRAPHWEB_CANONICAL_FIGURE_DIR),
# under a run-named subfolder. Run dirs stay the source of truth during active
# experimentation; this keeps a browsable, centralised index to pick figures
# from later, instead of figures being scattered across
# /pscratch/.../flowjax_inference_outputs/<variant>/.
#
# Usage: scripts/sync_figures_to_canonical.sh <src_dir> <run_name>
set -euo pipefail

SRC="${1:?usage: sync_figures_to_canonical.sh <src_dir> <run_name>}"
RUN_NAME="${2:?usage: sync_figures_to_canonical.sh <src_dir> <run_name>}"
CANONICAL_ROOT="${GRAPHWEB_CANONICAL_FIGURE_DIR:-/pscratch/sd/d/dkololgi/graphweb_desi/figures}"
DEST="${CANONICAL_ROOT}/${RUN_NAME}"

mkdir -p "$DEST"
rsync -av --include='*.png' --include='*.pdf' --include='*.html' --include='*.svg' \
  --include='*.gif' --include='*.mp4' \
  --exclude='*' "$SRC"/ "$DEST"/
echo "Synced figures: $SRC -> $DEST"
