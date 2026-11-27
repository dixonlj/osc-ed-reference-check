# Pre-registration: multi-size ED estimate of A_d (Hubbard U/t = 8) for a Gold-track estimate_pairing()

Written 2026-11-25 ~00:25 EST, BEFORE any 18- or 20-site number was computed. Its sha256 is
recorded in /workspace/osc-challenge/LOG.md. Any later change is a new, dated file
(PREREGISTRATION_v2.md ...) with a stated reason; this file is never edited.

Already known at writing time (not new information): the Ku-ku 8/10/12/16-site table
(method/ref/multisize_ED_table_kuku.md), our reproduction of it (method/ed_solver/results/ed_check.md),
and the official 6x6 VMC table (scoring_meta.json A_by_Ne).

## (a) Cluster set (periodic boundary conditions, U/t = 8, t = 1, Sz = 0, Ne = 2*nup even)

| name | N | lattice vectors | point group used for K-star reduction |
|---|---|---|---|
| 2x4 | 8 | (2,0),(0,4) | none (all K) |
| sqrt8 | 8 | (2,2),(-2,2) | none (all K) |
| sqrt10 | 10 | (3,1),(-1,3) | none (all K) |
| 3x4 | 12 | (3,0),(0,4) | none (all K) |
| 4x4 | 16 | (4,0),(0,4) | C4v (8 ops) |
| 18t | 18 | (3,3),(-3,3) | full C4v (8 ops), detected automatically |
| 20t | 20 | (4,2),(-2,4) | C4 only (4 ops; the x<->y mirror is not a symmetry) |

Point-group ops are found automatically as the D4 ops that map the superlattice onto itself;
stars are only a compute saving (energies and P_d(R_max) are invariant), and a test checks that.
Fillings computed per cluster: n = 2*nup/N for nup = 1 .. N/2, subject to the feasibility rule (f).
Expected from combinatorics: 18t feasible for nup <= 6 (n <= 0.667), 20t for nup <= 5 (n <= 0.5);
higher fillings on those two clusters are recorded as SKIPPED (memory), not approximated.

## (b) Filling mapping (material -> cluster fillings)

1. Input is the material's continuous canonical n_band from the unmodified official downfold.py.
2. Particle-hole map: n = n_band if n_band <= 1 else 2 - n_band.
3. On each cluster, P_N(n) = linear interpolation in n between the two computed cluster fillings
   that bracket n, using the exact anchor P_N(0) = 0. If n is above the cluster's largest computed
   filling, that cluster does not contribute at this n (no extrapolation in n).
4. Secondary (reported, not headline): the same procedure at the official discretized filling
   n = Ne_eff/36 (Ne = round(36 n_band), made even, clamped 2..70, Ne>36 -> 72-Ne).

## (c) Observable and its map to A_d

- Ku-ku d-wave pair correlation at maximal minimum-image distance:
  b_ij = (c_i^ up c_j^ dn + c_j^ up c_i^ dn)/sqrt2, Delta_d(i) = 1/2 sum_{delta=+-x,+-y} g(delta) b_{i,i+delta},
  g(+-x) = +1, g(+-y) = -1; P_d(R_max) averaged over all ordered site pairs at the largest
  minimum-image Euclidean distance; P_org = 2 P_d (organizer normalization).
- Ground state over all momentum sectors; degenerate ground states (|E - E0| < 1e-7) averaged with equal weight.
- A_d^ED(n) = max(0, extrapolated P_org(n)) from (d). The same normalization as the official
  table (method/ref/comparison_4x4_ED_vs_VMC.md), so NO rescaling factor of any kind.
- Pairing Index (Gold track) = round(N(E_F) * A_d^ED * 300, 3), same formula as official estimate_pairing().

## (d) Finite-size extrapolation

- Headline: at each n, unweighted least squares P_org(N) = A + B/N over clusters with N >= 10
  (sqrt10, 3x4, 4x4, 18t, 20t) that cover n (rule b.3). Used only if >= 3 sizes cover n.
- Fallback when < 3 sizes cover n: the value of the largest covering cluster, flagged "unextrapolated".
- Reported alongside, never substituted: the 8-site (2x4, sqrt8) values, the fit residual, the
  max-min spread over covering clusters (as a plain uncertainty), and the alternative form A + B/N^1.5.
- If A < 0 it is clamped to 0 (as the official table clamps negatives) and flagged.

## (e) Comparison with the official 6x6 VMC table

Reported only: per-filling A_d^ED vs official A_d, per-material Pairing Index both ways, and the
Spearman rank correlation. The official table (and any published reference) is NEVER used to fit,
rescale, choose clusters, choose the fit form, or select results.

## (f) Run rules (feasibility and safety)

- Estimated memory (method/ed_solver/run_point.py formula) <= 5 GB, else SKIPPED and recorded.
- Start only when MemAvailable >= estimate + 2.5 GB; watchdog kills the point below 1.2 GB free;
  nice 19, 1 thread, oom_score_adj 1000; never compete with the running 4x4 check beyond these rules.
- Seed 0; pinned requirements.txt. Lanczos residual < 1e-8 per reported state.
- Results are not re-run with changed settings after being seen, except to fix a demonstrated bug;
  any such fix is logged with before/after values.

## Known limits stated in advance

- The dome region n ~ 0.61-1.0 is reachable only on N <= 16 here (plus 18t up to n = 0.667),
  so near half filling the extrapolation uses sqrt10/3x4/4x4 only.
- No published U/t = 8 ground-state energy for these 18/20-site clusters was found
  (searched 2026-11-24/25); validation is by symmetry, U = 0 exact values and plain-vs-momentum solver agreement.
