#!/bin/bash
set -euo pipefail
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PATH="/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin:$PATH"
loa_output=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1/page_fullcap_20261026_v2
loa_render=/global/u2/d/dkololgi/.codex/plugins/cache/openai-bundled/visualize/1.0.50/skills/visualize/scripts/render.py
python "$loa_render" "$loa_output/NGC.html" "$loa_output/NGC.preview.html" --force
python "$loa_render" "$loa_output/SGC.html" "$loa_output/SGC.preview.html" --force
python /global/u2/d/dkololgi/GraphWeb_DESI/workflows/presentation/test_loa_fullcap_page_20261026.py
