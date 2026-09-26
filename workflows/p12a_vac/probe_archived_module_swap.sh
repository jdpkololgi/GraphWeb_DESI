#!/bin/bash
# Diagnostic subshell only. Never calls fba_run or changes the shared installation.
# Deliberately continues after swap errors to test archived-script behavior.
module use /global/common/software/desi/perlmutter/desiconda/20240425-2.2.0/modulefiles
module load fiberassign/main
printf 'load_status=%s\n' "$?"
printf 'before=%s\n' "$(command -v fba_run)"
module swap fiberassign/4.0.0
printf 'swap4_status=%s\n' "$?"
printf 'after4=%s\n' "$(command -v fba_run)"
module swap fiberassign/5.0.0
printf 'swap5_status=%s\n' "$?"
printf 'after5=%s\n' "$(command -v fba_run)"
