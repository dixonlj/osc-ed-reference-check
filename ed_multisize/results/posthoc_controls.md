# Post-hoc controls on the locked 1/N extrapolation (2026-11-26 13:26 EST)

**These analyses are post-hoc.** They are not in PREREGISTRATION.md and were not used to choose
the locked method. PREREGISTRATION.md sha256 `504a9e21047076de4039610f67c502f727b49b69100a74ecd681de44d8f64564` (unchanged). Script sha256 `7dcaea5743c44cfaaa5bdd56ce56dea72bf98c064a376aa376a4278c262dc93e`.
U=0 values are taken from Ku-ku's published table (canonical ED; 4x4 n=1 is canonical Wick).
Nothing is fitted or rescaled to the official 6x6 table.

## U=0 vs U=8 through the same pipeline (3 clusters: sqrt10, 3x4, 4x4)

| n | U=0 A (1/N) | U=8 A (same 3) | U=8 − U=0 | locked U=8 A (all covering N≥10) | official A_d |
|---|---|---|---|---|---|
| 0.0556 | 0.0000 | 3.301e-06 | 3.301e-06 | -6.263e-06 | 0.00283 |
| 0.2778 | 0.0106 | -0.0063 | -0.0169 | -0.0076 | 0.00718 |
| 0.3333 | 0.0294 | -0.0521 | -0.0815 | -0.0198 | 0.00483 |
| 0.3889 | 0.0558 | -0.0660 | -0.1219 | -0.0360 | 0.00262 |
| 0.4444 | 0.0762 | 0.0354 | -0.0408 | -0.0202 | 0 |
| 0.5000 | 0.0917 | 0.1468 | 0.0551 | -8.815e-04 | 0 |
| 0.5556 | 0.1214 | 0.1907 | 0.0693 | 0.0530 | 0 |
| 0.6111 | 0.1462 | 0.2152 | 0.0690 | 0.0954 | 0.01269 |
| 0.6667 | 0.0659 | 0.1086 | 0.0427 | 0.0609 | 0.04299 |
| 0.7222 | -0.0410 | -0.0108 | 0.0302 | -0.0108 | 0.05849 |
| 0.7778 | -0.1359 | -0.0958 | 0.0401 | -0.0958 | 0.04812 |
| 0.8333 | -0.2335 | -0.0927 | 0.1408 | -0.0927 | 0.05569 |
| 0.8889 | -0.3419 | -0.0641 | 0.2778 | -0.0641 | 0.04726 |
| 0.9444 | -0.4323 | -0.0549 | 0.3774 | -0.0549 | 0.04523 |
| 1.0000 | -0.5227 | -0.0456 | 0.4771 | -0.0456 | 0.03844 |

## Locked U=8 fit, unclamped, with uncertainty (post-hoc)

| n | A (1/N) ± SE | 95% CI (t, df=n−2) | A (1/√N) | leave-one-out (drop 4x4) | A at N=36 | official |
|---|---|---|---|---|---|---|
| 0.0556 | -0.000 ± 0.000 | [-5.769e-05, 4.517e-05] | -1.802e-05 | -4.439e-06 | -2.272e-06 | 0.00283 |
| 0.2778 | -0.008 ± 0.012 | [-0.0444, 0.0293] | -0.0175 | -0.0106 | -0.0034 | 0.00718 |
| 0.3333 | -0.020 ± 0.018 | [-0.0779, 0.0384] | -0.0399 | -0.0164 | -0.0110 | 0.00483 |
| 0.3889 | -0.036 ± 0.020 | [-0.0995, 0.0275] | -0.0843 | -0.0268 | -0.0169 | 0.00262 |
| 0.4444 | -0.020 ± 0.047 | [-0.1695, 0.1292] | -0.0681 | -0.0222 | -0.0034 | 0 |
| 0.5000 | -0.001 ± 0.091 | [-0.2892, 0.2874] | -0.0437 | -0.0148 | 0.0117 | 0 |
| 0.5556 | 0.053 ± 0.131 | [-0.5098, 0.6159] | 0.0666 | 0.0152 | 0.0432 | 0 |
| 0.6111 | 0.095 ± 0.105 | [-0.3547, 0.5455] | 0.1566 | 0.0569 | 0.0690 | 0.01269 |
| 0.6667 | 0.061 ± 0.038 | [-0.1012, 0.2230] | 0.0914 | 0.0426 | 0.0483 | 0.04299 |
| 0.7222 | -0.011 ± 0.028 | [-0.3606, 0.3390] | -0.0566 | 0.0520 | 0.0048 | 0.05849 |
| 0.7778 | -0.096 ± 0.005 | [-0.1596, -0.0321] | -0.2329 | -0.0844 | -0.0477 | 0.04812 |
| 0.8333 | -0.093 ± 0.011 | [-0.2297, 0.0442] | -0.2324 | -0.0682 | -0.0438 | 0.05569 |
| 0.8889 | -0.064 ± 0.009 | [-0.1824, 0.0541] | -0.1656 | -0.0854 | -0.0282 | 0.04726 |
| 0.9444 | -0.055 ± 0.021 | [-0.3206, 0.2108] | -0.1349 | -0.1025 | -0.0262 | 0.04523 |
| 1.0000 | -0.046 ± 0.033 | [-0.4587, 0.3674] | -0.1043 | -0.1197 | -0.0243 | 0.03844 |

## Spearman rank correlation vs official A_d (15 fillings; post-hoc)

- **locked_clamped_Ad_vs_official:** rho = -0.123, p = 0.662 (permutation p = 0.675); Fisher-z 95% CI [-0.60, 0.42]; bootstrap 95% CI [-0.55, 0.32].
- **unclamped_A_vs_official:** rho = -0.484, p = 0.068 (permutation p = 0.068); Fisher-z 95% CI [-0.80, 0.04]; bootstrap 95% CI [-0.86, 0.05].
- **sqrt10_vs_official:** rho = 0.728, p = 0.002 (permutation p = 0.003); Fisher-z 95% CI [0.34, 0.90]; bootstrap 95% CI [0.29, 0.93].
- **4x4_vs_official:** rho = -0.014, p = 0.960 (permutation p = 0.962); Fisher-z 95% CI [-0.52, 0.50]; bootstrap 95% CI [-0.53, 0.61].

Organizers already reported that 6x6 VMC A_d and 4x4 ED are "essentially uncorrelated across filling (Spearman ≈ −0.1)" (`validation/comparison_4x4_ED_vs_VMC.md`). The locked clamped comparison is a multi-size extension of that finding.

## Output files

- `results/posthoc_by_filling.csv` sha256 `7417f62d3da07405785d2b4534189bc217e20d03f80060e95e313e120bd74d0e`
- `results/posthoc_materials.csv` sha256 `a080f89b67cc787ab97de87dc33fbdc81dff09c15dec757ec2688e834cb13482`
