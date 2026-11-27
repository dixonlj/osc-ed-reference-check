# Pre-registered 1/N extrapolation of A_d (Hubbard U/t=8), run 2026-11-26 02:53 EST

PREREGISTRATION.md sha256 `504a9e21047076de4039610f67c502f727b49b69100a74ecd681de44d8f64564` (matches locked value). Script `extrapolate.py` sha256 `ed9d8fda1076165156e971cd300fc73b1d8e0d4c02f265de5b52be8e93c591bc`.
Comparison with the official 6x6 VMC table is report-only; nothing was fitted or rescaled to it or to any published value.

## Implementation choices (not fixed by the pre-registration)

- Inputs: every ok U=8 run json in method/ed_solver/runs and method/ed_multisize/runs; non-ok (aborted/skipped) points are not used.
- When a filling has both a 'k' and a 'plain' run (4x4 n=0.375, 0.5), the 'k' value is used; the two agree to < 1e-11.
- Fillings evaluated: the official 6x6 fillings n = Ne/36 for every Ne in scoring_meta A_by_Ne (comparison grid), and each material's canonical n_band.
- Material inputs: canonical (official band) columns old_n (= n_band, 4 dp as stored), nef, old_A_d, old_index from best_method_package/results/band_sensitivity.csv; the alternative-band columns are not used.
- Secondary discretized filling (pre-reg b.4): Ne = round(36 n_band), made even, clamped 2..70, Ne>36 -> 72-Ne, n = Ne/36.
- Fit residual reported as RMS of (data - fit) over the clusters in the fit; with 3 clusters the fit has 1 degree of freedom.
- Official values compared as used in scoring: max(0, A_by_Ne). Spearman rank correlation over the official fillings where both exist.
- Data status: 18t n=0.6667 (nup=6) aborted twice by the 1.2 GB memory watchdog (2026-11-25 15:21, 2026-11-26 02:42 EST; box memory taken by other work) and is LEFT OUT (flagged, per pre-reg f / Liam 2026-11-26), so 18t covers n<=0.556; 18t n>=0.778 and 20t n>=0.6 skipped by the pre-registered 5 GB budget.
- Pairing Index both ways = round(nef * A_d * 300, 3), comparison only; materials with no official index (blank) keep an ED value but are flagged.

## Cluster coverage (largest computed filling)

| cluster | N | computed n | n_max |
|---|---|---|---|
| 2x4 | 8 | 0.250, 0.500, 0.750, 1.000 | 1.000 |
| sqrt8 | 8 | 0.250, 0.500, 0.750, 1.000 | 1.000 |
| sqrt10 | 10 | 0.200, 0.400, 0.600, 0.800, 1.000 | 1.000 |
| 3x4 | 12 | 0.167, 0.333, 0.500, 0.667, 0.833, 1.000 | 1.000 |
| 4x4 | 16 | 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000 | 1.000 |
| 18t | 18 | 0.111, 0.222, 0.333, 0.444, 0.556 | 0.556 |
| 20t | 20 | 0.100, 0.200, 0.300, 0.400, 0.500 | 0.500 |

## A_d^ED at the official 6x6 fillings

| Ne | n | A_d^ED | raw A (fit) | B | fit RMS | spread | alt A (1/N^1.5) | clusters used (P_N(n)) | 8-site (report only) | official A_d | flag |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.0556 | 0 | -6.26302e-06 | 0.000143689 | 7.19e-06 | 1.92e-05 | -2.36159e-06 | sqrt10=-1.69589e-33;3x4=1.91755e-05;4x4=-6.29271e-33;18t=4.38139e-26;20t=2.57646e-23 | 2x4=0.000480066;sqrt8=1.46585e-33 | 0.00283 | clamped to 0 (negative) |
| 10 | 0.2778 | 0 | -0.00755673 | 0.149583 | 0.00514 | 0.0168 | -0.00429226 | sqrt10=0.0124892;3x4=-0.0043104;4x4=0.00635516;18t=1.15523e-05;20t=0.000232601 | 2x4=0.00151625;sqrt8=0.00513387 | 0.00718 | clamped to 0 (negative) |
| 12 | 0.3333 | 0 | -0.019756 | 0.315008 | 0.00812 | 0.0279 | -0.0131132 | sqrt10=0.02141;3x4=-0.00649436;4x4=-0.00512256;18t=6.4776e-05;20t=0.00205253 | 2x4=0.00022815;sqrt8=0.0154016 | 0.00483 | clamped to 0 (negative) |
| 14 | 0.3889 | 0 | -0.0359936 | 0.685997 | 0.00887 | 0.0372 | -0.0200015 | sqrt10=0.0303309;3x4=0.0302956;4x4=-0.00686343;18t=-0.00121538;20t=0.00853639 | 2x4=-0.00105995;sqrt8=0.0256693 | 0.00262 | clamped to 0 (negative) |
| 16 | 0.4444 | 0 | -0.0201707 | 0.604633 | 0.0209 | 0.0696 | -0.00430201 | sqrt10=0.0167341;3x4=0.0670855;4x4=0.0206062;18t=-0.00249553;20t=0.00967753 | 2x4=-0.00234805;sqrt8=0.0359371 | 0 | clamped to 0 (negative) |
| 18 | 0.5000 | 0 | -0.000881543 | 0.451266 | 0.0403 | 0.109 | 0.0133201 | sqrt10=-0.00249214;3x4=0.103875;4x4=0.0480759;18t=-0.00478003;20t=0.009483 | 2x4=-0.00363615;sqrt8=0.0462048 | 0 | clamped to 0 (negative) |
| 20 | 0.5556 | 0.0530399 | 0.0530399 | -0.352509 | 0.0419 | 0.103 | 0.0486349 | sqrt10=-0.0217184;3x4=0.0812706;4x4=0.0534297;18t=-0.00706454 | 2x4=0.0402488;sqrt8=0.039296 | 0 |  |
| 22 | 0.6111 | 0.215196 | 0.215196 | -2.27033 | 0.0231 | 0.0888 | 0.15854 | sqrt10=-0.029983;3x4=0.0586657;4x4=0.0587835 | 2x4=0.0841338;sqrt8=0.0323871 | 0.01269 |  |
| 24 | 0.6667 | 0.108646 | 0.108646 | -0.980012 | 0.00644 | 0.0378 | 0.083765 | sqrt10=0.00559868;3x4=0.0360609;4x4=0.0433583 | 2x4=0.128019;sqrt8=0.0254782 | 0.04299 |  |
| 26 | 0.7222 | 0 | -0.0107621 | 0.559145 | 0.00507 | 0.022 | 0.00447796 | sqrt10=0.0411804;3x4=0.0429828;4x4=0.0210068 | 2x4=0.171904;sqrt8=0.0185694 | 0.05849 | clamped to 0 (negative) |
| 28 | 0.7778 | 0 | -0.0958195 | 1.73305 | 0.000923 | 0.0648 | -0.0503495 | sqrt10=0.076762;3x4=0.0499047;4x4=0.0119173 | 2x4=0.175957;sqrt8=0.00437292 | 0.04812 | clamped to 0 (negative) |
| 30 | 0.8333 | 0 | -0.0927443 | 1.76125 | 0.00198 | 0.0657 | -0.0464096 | sqrt10=0.0818256;3x4=0.0568266;4x4=0.0160897 | 2x4=0.140178;sqrt8=-0.0171111 | 0.05569 | clamped to 0 (negative) |
| 32 | 0.8889 | 0 | -0.0641408 | 1.29342 | 0.00171 | 0.0488 | -0.0304926 | sqrt10=0.0665439;3x4=0.0412268;4x4=0.0177721 | 2x4=0.104399;sqrt8=-0.0385952 | 0.04726 | clamped to 0 (negative) |
| 34 | 0.9444 | 0 | -0.0548856 | 1.03131 | 0.00385 | 0.0393 | -0.0283528 | sqrt10=0.0512621;3x4=0.0256269;4x4=0.0119846 | 2x4=0.0686198;sqrt8=-0.0600792 | 0.04523 | clamped to 0 (negative) |
| 36 | 1.0000 | 0 | -0.0456305 | 0.769204 | 0.00598 | 0.0298 | -0.026213 | sqrt10=0.0359803;3x4=0.0100271;4x4=0.00619709 | 2x4=0.0328408;sqrt8=-0.0815632 | 0.03844 | clamped to 0 (negative) |

Spearman rank correlation A_d^ED vs official A_d over 15 fillings: -0.123

## Materials (canonical band; comparison only, not a submission)

| osc_id | formula | n_band | n_eff | A_d^ED | index (ED) | official A_d | official index | A_d^ED at Ne_eff/36 | flag |
|---|---|---|---|---|---|---|---|---|---|
| OSC-01188 | Cr2FeTe4 | 0.7306 | 0.7306 | 0 | 0 | 0.05849 | 43.627 | 0 | clamped to 0 (negative) |
| OSC-01196 | Cr2S4 | 0.7235 | 0.7235 | 0 | 0 | 0.05849 | 42.131 | 0 | clamped to 0 (negative) |
| OSC-02321 | C2F2Mn3 | 0.8667 | 0.8667 | 0 | 0 | 0.05569 | 41.276 | 0 | clamped to 0 (negative) |
| OSC-01294 | Mo2Te4 | 1.28 | 0.7200 | 0 | 0 | 0.05849 | 40.385 | 0 | clamped to 0 (negative) |
| OSC-01181 | CoCr2Te4 | 0.8637 | 0.8637 | 0 | 0 | 0.05569 | 37.928 | 0 | clamped to 0 (negative) |
| OSC-01197 | Cr2Sb2Te6 | 0.8498 | 0.8498 | 0 | 0 | 0.05569 | 37.38 | 0 | clamped to 0 (negative) |
| OSC-01238 | H2Mn2Se4 | 1.1575 | 0.8425 | 0 | 0 | 0.05569 | 36.048 | 0 | clamped to 0 (negative) |
| OSC-01193 | Cr2PdSe4 | 1.163 | 0.8370 | 0 | 0 | 0.05569 | 33.418 | 0 | clamped to 0 (negative) |
| OSC-01287 | Mn6Se8 | 0.8446 | 0.8446 | 0 | 0 | 0.05569 | 33.382 | 0 | clamped to 0 (negative) |
| OSC-01387 | H2Mn2NO2 | 0.8501 | 0.8501 | 0 | 0 | 0.05569 | 32.468 | 0 | clamped to 0 (negative) |
| OSC-01257 | I2S2W2 | 1.2821 | 0.7179 | 0 | 0 | 0.05849 | 31.918 | 0 | clamped to 0 (negative) |
| OSC-02297 | AuBrI3Mn | 1.3012 | 0.6988 | 0.0395803 | 21.286 | 0.05849 | 31.456 | 0 |  |
| OSC-01285 | Mn3Se4 | 0.7984 | 0.7984 | 0 | 0 | 0.04812 | 31.387 | 0 | clamped to 0 (negative) |
| OSC-01068 | Al2Mn2S5 | 0.8403 | 0.8403 | 0 | 0 | 0.05569 | 30.367 | 0 | clamped to 0 (negative) |
| OSC-02255 | Ag2Br6 | 0.886 | 0.8860 | 0 | 0 | 0.04726 | 30.315 | 0 | clamped to 0 (negative) |
| OSC-01354 | C2H2Mn3O2 | 0.8422 | 0.8422 | 0 | 0 | 0.05569 | 29.569 | 0 | clamped to 0 (negative) |
| OSC-01227 | Ga2Mn2Te5 | 1.1874 | 0.8126 | 0 | 0 | 0.05569 | 29.477 | 0 | clamped to 0 (negative) |
| OSC-01070 | Al2Mn2Se5 | 0.7419 | 0.7419 | 0 | 0 | 0.05849 | 28.489 | 0 | clamped to 0 (negative) |
| OSC-01194 | Cr2PdTe4 | 1.0642 | 0.9358 | 0 | 0 | 0.04523 | 28.311 | 0 | clamped to 0 (negative) |
| OSC-01220 | Fe5GeTe2 | 1.2097 | 0.7903 | 0 | 0 | 0.04812 | 28.203 | 0 | clamped to 0 (negative) |
| OSC-01253 | Mn2I2S | 0.7121 | 0.7121 | 0.010994 | 12.033 | 0.05849 | 64.018 | 0 |  |
| OSC-02304 | Fe2Br2 | 0.9626 | 0.9626 | 0 | 0 | 0.04523 | 55.373 | 0 | clamped to 0 (negative) |

Spearman rank correlation of Pairing Index (ED vs official) over 22 materials: 0.184

## Output files

- `results/extrapolation_inputs.csv` sha256 `d1936ed9055c6a67e5c6872f742bf29e1707d2f7292fd8e6c2bd3d44b24638fb`
- `results/extrapolation_by_filling.csv` sha256 `770b699837981f89fe84b16be1a583565383d0d43d6939a911a2025645cbb7ba`
- `results/extrapolation_materials.csv` sha256 `4d7bb5f6c8d7800dd1731c8fe4e51f4714f6cc9071c3c88ccbf95e0e9c2bdbee`
