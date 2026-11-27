# Citations, data credits and licences

## Licences
- **Code** in this package: MIT (see `LICENSE`).
- **Values** derived from OSC / C2DB, and official OSC files fetched by `fetch_official.sh`:
  **CC-BY-NC 4.0** (non-commercial). `inputs/materials_canonical.csv` is a derived table of
  official-pipeline n_band, N(E_F) and A_d for 22 materials (C2DB UIDs only; no C2DB data dump).

Official third-party files (the Ku-ku table, the 4×4 ED-vs-VMC note, `scoring_meta.json`, and
`instrument/osc_instrument.py`, sha256 `6b774e52…5236`) are **not bundled**. `fetch_official.sh`
copies or downloads them from a pinned dataset revision and checks sha256. The first three are
CC-BY-NC 4.0 data files; `osc_instrument.py` is official OSC instrument code. Attribution is below.

## C2DB
- S. Haastrup et al., "The Computational 2D Materials Database: high-throughput modeling and discovery of atomically thin crystals", *2D Materials* **5**, 042002 (2018).
- M. N. Gjerding et al., "Recent progress of the Computational 2D Materials Database (C2DB)", *2D Materials* **8**, 044002 (2021).

## Open Superconductor Challenge
- FINAL-Bench / VIDRAFT, dataset `FINAL-Bench/OSC-Superconductor`, revision `346f184d3eff51ebf87ffc5da20c6ec74e4cd67e`: `scoring_meta.json` (official 6×6 A_d table) and the canonical downfold values used in `inputs/materials_canonical.csv`. CC-BY-NC 4.0.

## Ku-ku exact-diagonalization table
- Team Ku-ku, "Multi-size ED of the U/t=8 Hubbard model" (`validation/multisize_ED_table_kuku.md`, 1 Oct 2026) and `validation/comparison_4x4_ED_vs_VMC.md`, same dataset, CC-BY-NC 4.0, attribution to Ku-ku.
- Their caveats, which this write-up adopts: clusters with different shapes and different R_max "are not a strict finite-size sequence"; on these sizes they consider "only the sign and ranking of the 4×4 values usable for calibration." They also reported that 6×6 VMC A_d and 4×4 ED are "essentially uncorrelated across filling (Spearman ≈ −0.1)."

## Other ED packages
- Tuli99, OSC dataset discussion #59 (opened 20 October 2026; organizer reply 21 October 2026): an independent Hubbard ED package that also reproduces the Ku-ku table. Their published 3×4, N_up = N_dn = 5 value P_org = 0.0568266126268765 agrees with ours (0.056826612626876716) to about 1e-16. Cited as independent corroboration only.

## Hubbard-model literature
- M. Qin, C.-M. Chung, H. Shi, E. Vitali, C. Hubig, U. Schollwöck, S. R. White, S. Zhang, "Absence of Superconductivity in the Pure Two-Dimensional Hubbard Model", *Phys. Rev. X* **10**, 031016 (2020), https://doi.org/10.1103/PhysRevX.10.031016 (open access, CC-BY). Supports only short-range pairing at U/t around 6–8 and hole dopings 0.1 < h < 0.2 (n ≈ 0.8–0.9). Our small-cluster trend is consistent with that in that window; it is not a validation of it.
- B.-X. Zheng, C.-M. Chung, P. Corboz, G. Ehlers, M.-P. Qin, R. M. Noack, H. Shi, S. R. White, S. Zhang, G. K.-L. Chan, "Stripe order in the underdoped region of the two-dimensional Hubbard model", *Science* **358**, 1155 (2017), https://doi.org/10.1126/science.aam7127. Period-8 stripes near 1/8 doping, which 10–16-site clusters cannot hold.
- H. Q. Lin, "Ground-state properties of the two-dimensional Hubbard model", *Phys. Rev. B* **44**, 7151 (1991), https://doi.org/10.1103/PhysRevB.44.7151. Exact diagonalization of 4×4 and 18-site lattices versus U and filling. We did not compare our 18-site energies with this paper; it is paywalled and we did not obtain a copy. It may contain usable U/t = 8 values.
- T. Tohyama, Y. Inoue, K. Tsutsui, S. Maekawa, *Phys. Rev. B* **72**, 045113 (2005). 20-site study at U/t = 5 and 10, so it does not apply at U/t = 8.
