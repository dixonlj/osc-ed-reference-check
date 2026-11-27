# Draft post for the OSC dataset discussion board

**Status:** draft only. Not posted, not submitted. For Liam Dixon to review and post himself.
**Suggested title:** Reference / validation: ED check of OSC's pairing estimate (calibration, not a rescoring)

---

Entrant: Liam Dixon (self-taught builder). This is a calibration / reference package, not a score change or a prize claim.

I built, with AI agents, an exact solver for the same single-band Hubbard model OSC uses (U/t = 8, organizer-normalized d-wave operator). All 27 unique Ku-ku rows match under a pass/fail rule written down before comparing (29 checks, including two 4×4 fillings run with both solvers; largest pairing difference 3.9e-7). I added 18- and 20-site tilted clusters and a 1/N extrapolation whose analysis plan was locked before any of those new-cluster numbers (`PREREGISTRATION.md`, sha256 `504a9e21…`). Fits at n ≥ 0.722 still use only 8–16-site values already public then, so they are planned but not blind. The 18-site n=0.667 point was added late under a pre-committed addendum rule (results are also reported without it). The addendum file was re-saved after that result was known (the 9:37 AM ET version, sha256 `2070a11e…`, was overwritten at 1:28 PM ET; no copy of that earlier file survives).

A U=0 run of the same pipeline has the same up-then-down sign pattern as at U=8. Free electrons have a nonzero d-wave background with no pairing order, so that sign pattern is largely a finite-cluster shell and shape effect. At n ≥ 0.8 the magnitudes differ a lot. The official 6×6 “dome” is neither confirmed nor ruled out: locked intercepts go negative there, but the intervals are wide and the fit is sensitive to dropping the 16-site cluster. Spearman ρ versus the official table is −0.12 (p = 0.66), in line with the organizers’ published figure of about −0.1.

The package includes an optional ED `estimate_pairing()`, off by default (`OSC_ED_PAIRING=1`), so anyone can test the method inside the official instrument. With it on, the locked, clamped A_d is used and the output is labelled as an ED calibration, not a submission; 20 of 22 materials score 0. That is not a proposed replacement: the clusters are small, the uncertainty is large, and the U=0 control shows the clamp is not a pairing measurement.

Code is MIT. OSC and C2DB-derived values are CC-BY-NC 4.0 (non-commercial use only). Data credits: S. Haastrup et al., 2D Mater. 5, 042002 (2018); M. N. Gjerding et al., 2D Mater. 8, 044002 (2021). Reproduce with `SHA256SUMS` and the `run_*.sh` scripts in the package. Package / repo: https://github.com/dixonlj/osc-ed-reference-check. Every number comes from the deposited code.

---
