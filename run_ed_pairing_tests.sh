#!/bin/bash
# Cheap tests for the optional, off-by-default ED estimate_pairing(). Seconds, tiny memory.
# (a) default off vs the fetched official estimate_pairing and CLI
# (b) opt-in vs ed_multisize/results/extrapolation_materials.csv
# (c) documented range behaviour (including a real unextrapolated expected value)
# (d) opt-in wrapper labelled as an ED calibration, not a submission
# Artefacts go to a temp dir unless OSC_ED_PAIRING_OUT is set.
set -euo pipefail
cd "$(dirname "$0")"
source env.sh
# Prefer a local official cache so tests can compare against the fetched file
# without leaving osc_instrument.py in the deposited third_party/ folder.
if [ -z "${OSC_OFFICIAL_CACHE:-}" ] && [ -f ../ed_reference_official_cache/osc_instrument.py ]; then
  export OSC_OFFICIAL_CACHE="$(cd .. && pwd)/ed_reference_official_cache"
fi
echo "504a9e21047076de4039610f67c502f727b49b69100a74ecd681de44d8f64564  ed_multisize/PREREGISTRATION.md" | sha256sum -c -
$PY tests/test_ed_pairing.py
