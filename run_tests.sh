#!/bin/bash
# Unit tests: ed_solver (2x2 vs dense Jordan-Wigner, U=0 vs Wick, particle-hole, K vs plain) and
# ed_multisize (18t/20t geometry, [H,T]=[H,P]=0, U=0 exact, K vs plain at U=8). ~2 + ~20 min, <1 GB.
set -euo pipefail; cd "$(dirname "$0")"; source env.sh; ./fetch_official.sh
nice -n 19 $PY ed_solver/tests/test_ed.py | tee rerun_check/unit_tests_ed_solver.txt
nice -n 19 $PY ed_multisize/tests/test_geometry.py | tee rerun_check/unit_tests_ed_multisize.txt
