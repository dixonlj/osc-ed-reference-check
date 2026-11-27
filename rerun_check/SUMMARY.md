# Re-run check from this package (2026-11-26, America/New_York)

Fresh virtualenv from `requirements.txt` only (Python 3.13.5, numpy 2.5.5, scipy 1.18.1), one thread, nice 19.
No large-memory runs were started (shared machine).

The `rerun_check/extrapolation/` and `rerun_check/posthoc/` folders were refreshed on 26 Nov 2026,
about 1:40 PM ET, from a fresh copy of this package, using the documented commands
`PY=python3 ./run_extrapolate.sh` and `PY=python3 ./run_posthoc.sh`. That re-run includes the
completed 18t n=0.667 point. It replaces earlier folders (7:47 AM and 9:20 AM ET) that still held
the without-18t outputs.

| check | script | result |
|---|---|---|
| Official inputs fetched and sha256-verified (rev 346f184d) | `fetch_official.sh` | 3/3 OK |
| Locked extrapolation re-run, with 18t n=0.667 included (PREREGISTRATION.md sha256 verified; 26 Nov 2026, 1:40 PM ET, fresh copy) | `run_extrapolate.sh` | 3/3 CSVs byte-identical to `ed_multisize/results/` (`extrapolation_by_filling.csv`, `extrapolation_materials.csv`, `extrapolation_inputs.csv`) |
| Post-hoc controls re-run, same fresh copy and data set (26 Nov 2026, 1:40 PM ET) | `run_posthoc.sh` | `posthoc_by_filling.csv`, `posthoc_materials.csv` and `posthoc_spearman.json` byte-identical to `ed_multisize/results/` |
| Independent-seed re-run (seed 20261126 vs production seed 0) | `run_rerun_small.sh` | **32/32 match**; max dE0 2.8e-13, max dP_org 5.2e-12 (`rerun_small.md`) |
| ed_solver unit tests | `run_tests.sh` | 35/35 pass (`unit_tests_ed_solver.txt`) |
| ed_multisize unit tests (18t/20t) | `run_tests.sh` | 48/48 pass (`unit_tests_ed_multisize.txt`) |
| Optional ED `estimate_pairing` tests (26 Nov 2026, about 2:41 PM ET) | `run_ed_pairing_tests.sh` | (a) default off matches the fetched official `estimate_pairing` and CLI on 22 materials; (b) opt-in matches `extrapolation_materials.csv` (20 of 22 score 0); (c) range behaviour, including unextrapolated A=0.008; (d) opt-in wrapper labelled as an ED calibration, not a submission. Tests write to a temp dir by default; copies of the CSVs and calibration JSON are in `rerun_check/ed_pairing/`. |

Not re-run here (hours, several GB): 4x4 n >= 0.625, 18t n = 0.556, 18t n = 0.667, 20t n = 0.5.
Their stored results are the original production runs (4x4 rows independently matched to Ku-ku;
18t n=0.667 is the completed retry of 26 Nov 2026, 12:48 PM ET). Command: `run_full_check.sh`.

`fetch_official.sh` now also pins official `instrument/osc_instrument.py` (sha256 `6b774e52…5236`).
The 3/3 fetch row above is the earlier data-file check used for the extrapolation re-run; those three
data files are unchanged. Official files remain fetch-only and are not stored in this package.
