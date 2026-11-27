#!/bin/bash
# FULL reproduction (heavy; NOT run as part of the quick check). Needs ~3.5 GB free beyond the 2.5 GB reserve
# for the largest points; 4x4 n=1.0 took ~6 h and 18t n=0.556 / 20t n=0.5 ~0.6 / ~2 h on one core.
# Each point waits for free memory, runs at nice 19 with one thread and a 1.2 GB watchdog, and is resumable.
# Writes into ed_solver/runs and ed_multisize/runs (overwrites stored results!) - copy the package first.
set -euo pipefail; cd "$(dirname "$0")"; source env.sh; ./fetch_official.sh
[ "${I_UNDERSTAND_THIS_IS_HEAVY:-no}" = yes ] || { echo "set I_UNDERSTAND_THIS_IS_HEAVY=yes"; exit 1; }
(cd ed_solver && nice -n 19 $PY run_check.py)
(cd ed_multisize && nice -n 19 $PY run_multisize.py)
./run_extrapolate.sh
