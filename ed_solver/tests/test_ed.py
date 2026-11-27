#!/usr/bin/env python3
"""Unit tests for the Hubbard ED solver (plain asserts; run: python tests/test_ed.py).

1. 2x2 cluster vs an independent dense Jordan-Wigner construction of H and Delta_d.
2. U = 0: ED equals the closed-shell Wick (Slater determinant) value, and matches Ku-ku's
   U = 0 column where their U = 0 ground state is non-degenerate.
3. Particle-hole symmetry on the bipartite 2x4 cluster:
   E0(N-n, N-n) = E0(n, n) + U (N - 2n), same degeneracy, same P_d(R_max).
4. Momentum-sector solver equals the plain solver (2x4, 3x4, 4x4 low filling).
"""
import itertools, math, os, sys, time
import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from lattice import Cluster, GFORM
from hubbard import Sector, pd_correlation, wick_pd, pair_terms, configs
from hubbard_k import KSector, translation_tables, c4v_stars
import run_check

PASS = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    PASS.append(bool(cond))


# ---------- 1. 2x2 vs dense Jordan-Wigner ----------
def jw_ops(M):
    I = sp.identity(2, format="csr"); Z = sp.diags([1.0, -1.0]).tocsr()
    a = sp.csr_matrix(np.array([[0.0, 1.0], [0.0, 0.0]]))  # |0>=occ? use basis (empty, occ): a|occ>=|empty>
    ops = []
    for m in range(M):
        mats = [Z] * m + [a] + [I] * (M - m - 1)
        op = mats[0]
        for x in mats[1:]:
            op = sp.kron(op, x, format="csr")
        ops.append(op)
    return ops


def test_2x2_dense():
    c = Cluster("2x2"); N = c.N; M = 2 * N
    # mode order: up sites 0..N-1, then down sites 0..N-1 (matches the product-state ordering)
    # local basis (|empty>, |occ>) with a = [[0,1],[0,0]]; Z = diag(1,-1) counts parity of occupied
    ops = jw_ops(M)
    cu = ops[:N]; cd = ops[N:]
    dimF = 2 ** M
    H = sp.csr_matrix((dimF, dimF))
    for i, j in c.hop_bonds():
        for cs in (cu, cd):
            H = H + (-1.0) * (cs[j].T @ cs[i])
    for U in (0.0, 8.0):
        HU = H + U * sum((cu[i].T @ cu[i]) @ (cd[i].T @ cd[i]) for i in range(N))
        nu = sum(cu[i].T @ cu[i] for i in range(N)).diagonal()
        nd = sum(cd[i].T @ cd[i] for i in range(N)).diagonal()
        R = c.r_max()
        def delta(i):
            out = sp.csr_matrix((dimF, dimF))
            for coef, p, q in pair_terms(c, i):
                out = out + coef * (cu[p] @ cd[q])
            return out
        Dl = [delta(i) for i in range(N)]
        for n in (1, 2):
            sel = np.nonzero((np.round(nu) == n) & (np.round(nd) == n))[0]
            Hs = HU[sel][:, sel].toarray()
            w, v = np.linalg.eigh(Hs)
            E0 = w[0]; g = np.nonzero(np.abs(w - E0) < 1e-7)[0]
            # P via full Fock vectors
            Ps = []
            for s in g:
                psi = np.zeros(dimF); psi[sel] = v[:, s]
                tot = cnt = 0
                for i in range(N):
                    for j in c.shell(i, R):
                        tot += (Dl[j] @ psi) @ (Dl[i] @ psi); cnt += 1
                Ps.append(tot / cnt)
            Pjw = float(np.mean(Ps))
            s_ = Sector(c, n, n, U)
            E0e, wE, vecs, how = s_.ground_states()
            Pe = float(np.mean([pd_correlation(s_, X)[1] for X in vecs]))
            spec_e = np.linalg.eigvalsh(s_.sparse_H().toarray())
            check(f"2x2 U={U:g} n_up=n_dn={n}: full spectrum equals Jordan-Wigner dense",
                  np.allclose(np.sort(spec_e), np.sort(w), atol=1e-10), f"max diff {np.max(np.abs(np.sort(spec_e)-np.sort(w))):.1e}")
            check(f"2x2 U={U:g} n={n}: degeneracy and P_d equal JW", len(vecs) == len(g) and abs(Pe - Pjw) < 1e-12,
                  f"deg {len(vecs)}/{len(g)}  P ours {Pe:.12g}  JW {Pjw:.12g}")


# ---------- 2. U = 0 limit ----------
def test_U0():
    ref = run_check.parse_ref()
    cases = [("2x2", 1), ("2x4", 1), ("2x4", 3), ("sqrt10", 5), ("3x4", 3), ("3x4", 5), ("3x4", 6), ("4x4", 5)]
    for name, n in cases:
        c = Cluster(name)
        w = wick_pd(c, n, n)
        if w is None:
            check(f"U=0 {name} n_up={n}: closed shell", False, "open shell"); continue
        msg = f"Wick P_org={2*w:.9g}"
        if c.N <= 12 and math.comb(c.N, n) ** 2 <= 1_000_000:
            s = Sector(c, n, n, 0.0)
            E0, wE, vecs, how = s.ground_states()
            Pe = float(np.mean([pd_correlation(s, X)[1] for X in vecs]))
            check(f"U=0 {name} n_up={n}: ED equals Wick", len(vecs) == 1 and abs(Pe - w) < 1e-10,
                  f"ED P_org={2*Pe:.9g} {msg}")
        k = ref.get((name, n))
        if k and k["deg_U0"] == 1 and k["P_org_U0"] is not None:
            check(f"U=0 {name} n_up={n}: Wick equals Ku-ku U=0 column",
                  abs(2 * w - k["P_org_U0"]) <= 1e-5 * abs(k["P_org_U0"]) + 1e-8,
                  f"{msg}  Ku-ku {k['P_org_U0']}")


# ---------- 3. particle-hole symmetry ----------
def test_particle_hole():
    c = Cluster("2x4"); U = 8.0; R = c.r_max()
    sh0 = set(c.nbr[d][0] for d in GFORM) | {0}
    j = c.shell(0, R)[0]
    shj = set(c.nbr[d][j] for d in GFORM) | {j}
    check("2x4: Delta(0) and Delta(j in R_max shell) share no site (PH test is valid)", not (sh0 & shj))
    for n in (1, 2, 3):
        a = Sector(c, n, n, U); b = Sector(c, c.N - n, c.N - n, U)
        Ea, _, va, _ = a.ground_states(); Eb, _, vb, _ = b.ground_states()
        Pa = np.mean([pd_correlation(a, X)[1] for X in va]); Pb = np.mean([pd_correlation(b, X)[1] for X in vb])
        check(f"PH 2x4 n_up={n} <-> {c.N-n}: E0 shift U(N-2n)", abs(Eb - (Ea + U * (c.N - 2 * n))) < 1e-9,
              f"{Ea:.10f} vs {Eb:.10f}")
        check(f"PH 2x4 n_up={n} <-> {c.N-n}: same degeneracy and P_d(R_max)",
              len(va) == len(vb) and abs(Pa - Pb) < 1e-10, f"deg {len(va)}/{len(vb)} P {Pa:.10g}/{Pb:.10g}")


# ---------- 4. K-sector solver vs plain ----------
def k_solver_value(c, n, U, use_stars):
    tabs = translation_tables(c, configs(c.N, n))
    if use_stars:
        ks, stars = c4v_stars(c)
    else:
        ks = c.k_points(); stars = [(i, [i]) for i in range(len(ks))]
    res = []
    for i, members in stars:
        s = KSector(c, n, n, U, ks[i], tables=tabs)
        E, Xs, ress, _ = s.lowest_lanczos()
        res.append((E, members, s, Xs))
    Emin = min(r[0] for r in res)
    P = d = 0
    for E, members, s, Xs in res:
        if abs(E - Emin) < 1e-7:
            for X in Xs:
                P += len(members) * s.pd(X)[1]; d += len(members)
    return Emin, d, P / d


def test_k_vs_plain():
    for name, n, stars in (("2x4", 2, False), ("2x4", 3, False), ("3x4", 3, False), ("3x4", 4, False),
                           ("4x4", 2, True), ("4x4", 3, True)):
        c = Cluster(name)
        s = Sector(c, n, n, 8.0)
        E0, _, vecs, _ = s.ground_states()
        P = float(np.mean([pd_correlation(s, X)[1] for X in vecs]))
        Ek, dk, Pk = k_solver_value(c, n, 8.0, stars)
        check(f"K-solver vs plain {name} n_up={n}{' (C4v stars)' if stars else ''}",
              abs(E0 - Ek) < 1e-9 and dk == len(vecs) and abs(P - Pk) < 1e-8,
              f"E {E0:.10f}/{Ek:.10f} deg {len(vecs)}/{dk} P {P:.10g}/{Pk:.10g}")


if __name__ == "__main__":
    t0 = time.time()
    test_2x2_dense()
    test_U0()
    test_particle_hole()
    test_k_vs_plain()
    print(f"\n{sum(PASS)}/{len(PASS)} checks passed in {time.time()-t0:.1f}s")
    sys.exit(0 if all(PASS) else 1)
