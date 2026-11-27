#!/bin/bash
# Independent-seed re-run of every small stored point (N<=12 all fillings; 4x4 n<=0.5; 18t/20t n<=0.44)
# and comparison with the stored results. Writes rerun_check/rerun_small.{csv,md}. ~10 min, <1 GB.
# Usage: SEED=12345 ./run_rerun_small.sh     (default seed 20261126; production used seed 0)
set -euo pipefail; cd "$(dirname "$0")"; source env.sh
export OSC_ED_SEED=${SEED:-20261126}
nice -n 19 $PY rerun_check/rerun_small.py
