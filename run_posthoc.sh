#!/bin/bash
# Post-hoc controls (U=0 pipeline, unclamped CIs, Spearman sensitivity). Seconds, tiny memory.
# Does not change the locked extrapolation CSVs.
set -euo pipefail; cd "$(dirname "$0")"; source env.sh; ./fetch_official.sh
echo "504a9e21047076de4039610f67c502f727b49b69100a74ecd681de44d8f64564  ed_multisize/PREREGISTRATION.md" | sha256sum -c -
mkdir -p rerun_check/posthoc
OSC_POSTHOC_OUT=rerun_check/posthoc $PY ed_multisize/posthoc_controls.py
echo "post-hoc written to rerun_check/posthoc/"
