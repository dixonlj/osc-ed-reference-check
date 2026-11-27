# ED reference check for the OSC pairing estimate (Hubbard U/t = 8)

A checked, reproducible exact-diagonalization (ED) solver for the single-band Hubbard model that OSC uses
for its pairing amplitude A_d, plus a planned multi-size study (partly blind) and a post-hoc U=0 control.
The size trend on these small clusters has the same up-then-down shape for free electrons as for U=8,
so it is a finite-cluster shell/shape effect and cannot confirm or rule out pairing.
**This is a validation / reference package, not a replacement score.** It does not change the
canonical downfold or any official file. An optional ED `estimate_pairing()` is included but is
**off by default** (`OSC_ED_PAIRING=1` to enable), so anyone can test the method inside the official
instrument. With the wrapper opt-in off, the official script is run without replacing
`estimate_pairing`, and the output (exit code, stdout, stderr, JSON) is byte-identical to running
the official file directly. With it on, the output is labelled as an ED calibration, not a
submission; 20 of 22 materials score 0, and this is not a proposed replacement.
Status: DRAFT v2.5, not submitted.

Plain-language write-up: `WRITEUP.md` (PDF: `WRITEUP.pdf`; HTML source used to build the PDF: `WRITEUP.html`).

## Contents

| Path | What |
|---|---|
| `WRITEUP.md` | Plain-language write-up (draft v2.5) |
| `WRITEUP.pdf` | Same write-up as a PDF (tables should be readable) |
| `WRITEUP.html` | HTML source used to build the PDF (Chrome, non-wrapping tables) |
| `ed_multisize/PREREGISTRATION_addendum.md` | Dated addendum (26 Nov 2026): committed with-and-without rule for the 18t n=0.667 point, now completed and included, with a disclosure that the file was re-saved after the result. Does not change the locked plan. |
| `ed_solver/` | ED code: `lattice.py` (clusters), `hubbard.py` (real-space sectors, eigsh + Lanczos deflation, pair correlation, Wick), `hubbard_k.py` (crystal-momentum sectors), `lanczos.py` (low-memory Lanczos), `run_point.py`, `run_check.py` |
| `ed_solver/results/ed_check.{md,csv}` | Correctness check vs Ku-ku: **all 27 unique rows match** (29 checks, including 4×4 n=0.375 and 0.5 with both solvers) under a rule fixed before comparing |
| `ed_solver/runs/*.json` | Stored per-point results (E0, degeneracy, P_org, residuals, runtime, peak memory) |
| `ed_multisize/PREREGISTRATION.md` | Locked analysis plan, sha256 `504a9e21047076de4039610f67c502f727b49b69100a74ecd681de44d8f64564` (fixed before any 18/20-site number; unchanged). Fits at n ≥ 0.722 still use only 8–16-site values already public then, so they are planned but not blind. Only n ≤ 0.556 was blind. The 18t n=0.667 point, obtained after other results were visible, is included under the addendum rule; n = 0.611 and 0.667 are reported both with and without it. Local timestamps are not independently verifiable. |
| `ed_multisize/` | Tilted 18- and 20-site clusters (`geometry.py`), resumable runner, orchestrator, `extrapolate.py` (locked 1/N extrapolation), `posthoc_controls.py` (U=0 control, unclamped CIs, Spearman sensitivity; labelled post-hoc) |
| `ed_multisize/results/` | locked `extrapolation*.{md,csv}` **with** the 18t n=0.667 point; post-hoc `posthoc_controls.md`, `posthoc_by_filling.csv`, `u0_control_chart.png`; unit tests |
| `ed_multisize/results_without_18t_nup6/` | The earlier locked and post-hoc results **without** that point (26 Nov 2026, 2:53 AM ET), kept for the committed comparison |
| `inputs/materials_canonical.csv` | 22 materials: canonical (official band) n_band, N(E_F), official A_d and index (derived OSC pipeline values; **CC-BY-NC 4.0**) |
| `config/frozen_config.json` | All frozen settings (model, operator, tolerances, seed 0, memory rules, versions) |
| `LICENSE` | MIT licence for the **code** |
| `CITATIONS.md` | Data credits; values are CC-BY-NC 4.0 |
| `rerun_check/` | Re-run evidence from this package (independent seed, byte-identical extrapolation, unit tests) |
| `third_party/` | Runtime only: filled by `fetch_official.sh`. Official files are **not bundled** (see `third_party/README.md`). Excluded from the manifest. |
| `ed_pairing.py` | Optional, off-by-default ED `estimate_pairing()`. Default: official Index = round(N(E_F)×A_d×300, 3). Opt-in (`OSC_ED_PAIRING=1` or `use_ed=True`): locked, clamped 1/N A_d from this package. Official `osc_instrument.py` is not bundled. |
| `wrap_osc_instrument.py` | Optional wrapper around a *fetched* official instrument (see the file's docstring). Default-off: runs the official script without replacing `estimate_pairing` (byte-identical exit code, stdout, stderr, JSON). Opt-in: writes a labelled ED calibration file, not a submission. Does not redistribute official code. |
| `tests/test_ed_pairing.py`, `run_ed_pairing_tests.sh` | Cheap tests: (a) default off vs the fetched official `estimate_pairing` and CLI; (b) opt-in vs `extrapolation_materials.csv`; (c) documented range behaviour, including a real unextrapolated expected value; (d) opt-in wrapper labels. Artefacts go to a temp dir unless `OSC_ED_PAIRING_OUT` is set. |
| `POST_DRAFT.md` | Draft discussion-board post for Liam to review and post himself. Not posted, not submitted. |
| `SHA256SUMS` | Manifest of every package file (excluding `third_party/`, `__pycache__/` and `SHA256SUMS` itself) |

## Reproduce

```bash
python3.13 -m venv venv && venv/bin/pip install -r requirements.txt
export PY=$PWD/venv/bin/python
sha256sum -c SHA256SUMS
./fetch_official.sh              # official Ku-ku table, VMC comparison, scoring_meta.json, osc_instrument.py (sha256-checked; not bundled)
./run_extrapolate.sh             # seconds: locked extrapolation from stored runs; CSVs must be byte-identical
./run_posthoc.sh                 # seconds: U=0 control, unclamped CIs, Spearman (does not change locked CSVs)
./run_ed_pairing_tests.sh        # seconds: optional estimate_pairing, off by default
SEED=20261126 ./run_rerun_small.sh  # ~10 min, <1 GB: independent-seed re-run vs stored results
# (The checked re-run in rerun_check/ used seed 20261126. Any other seed, e.g. 12345, is also a valid check.)
./run_tests.sh                   # ~20 min, <1 GB: unit tests (ed_solver 35 checks, ed_multisize 48 checks)
# heavy (hours, several GB; overwrites stored runs - use a copy):
# I_UNDERSTAND_THIS_IS_HEAVY=yes ./run_full_check.sh
```

All scripts run single-threaded at nice 19. Production results used seed 0.

## Rules followed
- Nothing is fitted, rescaled or tuned to the Ku-ku table, the official 6×6 table or any published value;
  those are compared only. The analysis plan was hashed before the new-cluster runs; fits at
  n ≥ 0.722 still use only then-already-public 8–16-site numbers (planned, not blind). The 18t
  n=0.667 point is included under a pre-committed rule; results are also reported without it.
- Canonical downfold and band choice are untouched; materials use the official band (`old_*` columns).
- **Code licence: MIT** (`LICENSE`). **Values: CC-BY-NC 4.0** (non-commercial). See `CITATIONS.md`.
- Official OSC files are fetched, not bundled.

## Changes relative to the working copies (method/ed_solver, method/ed_multisize)
Only paths and a seed hook; physics and settings unchanged:
1. `ed_solver/hubbard.py`: `SEED` reads optional env `OSC_ED_SEED` (default 0).
2. `ed_solver/run_check.py`: reference tables read from `third_party/` (after `fetch_official.sh`).
3. `ed_multisize/extrapolate.py`: `scoring_meta.json` from `third_party/`, materials from `inputs/materials_canonical.csv`,
   optional output dir `OSC_EXTRAP_OUT`.
4. `ed_multisize/posthoc_controls.py` is new in v2 (post-hoc; does not change locked outputs).
5. `ed_pairing.py` and `wrap_osc_instrument.py` were added in v2.4–v2.5: optional, off-by-default ED `estimate_pairing()`. With the wrapper opt-in off, official CLI output is byte-identical to the official script. With it on, output is labelled as an ED calibration, not a submission.

## Known limits (details in WRITEUP.md)
Small clusters (8–20 sites); jumpy size dependence; mixed measurement distances; 18t n=0.667 is
included after two watchdog aborts and a later successful retry (reported with and without);
no 18/20-site data above n=0.667, so 1/N fits at n ≥ 0.722 (where most materials sit) still use
only three already-public cluster sizes and cannot confirm or rule out pairing.
H. Q. Lin, PRB 44, 7151 (1991) has 18-site ED versus U and filling but was not checked (paywalled).
