#!/bin/bash
# Re-run the locked, pre-registered extrapolation from the stored run files into rerun_check/extrapolation/
# and check the CSVs are byte-identical to ed_multisize/results/. Seconds, tiny memory.
set -euo pipefail; cd "$(dirname "$0")"; source env.sh; ./fetch_official.sh
echo "504a9e21047076de4039610f67c502f727b49b69100a74ecd681de44d8f64564  ed_multisize/PREREGISTRATION.md" | sha256sum -c -
mkdir -p rerun_check/extrapolation
OSC_EXTRAP_OUT=rerun_check/extrapolation $PY ed_multisize/extrapolate.py
for f in extrapolation_by_filling.csv extrapolation_materials.csv extrapolation_inputs.csv; do
  cmp ed_multisize/results/$f rerun_check/extrapolation/$f && echo "IDENTICAL $f"
done
