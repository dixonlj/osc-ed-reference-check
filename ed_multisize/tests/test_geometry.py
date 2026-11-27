#!/usr/bin/env python3
"""Tests for the tilted 18t / 20t clusters (seed 0). Output: results/unit_tests.txt (via caller).
1 geometry and detected point groups; 2 H commutes with every translation and every detected
point-group op (full space); 3 U=0: K runner E0 = exact free-fermion sum, P_d = Wick (closed shells);
4 U=8: K runner (with stars) vs full-space plain Lanczos, E0/deg/P_org; 5 stars vs all-K energies.
No published U/t=8 18/20-site energy was found, so none is compared (stated, not invented)."""
import json, math, os, subprocess, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import geometry  # noqa
from lattice import Cluster  # noqa
from hubbard import Sector, pd_correlation, wick_pd, configs, lookup  # noqa
from hubbard_k import translation_tables  # noqa

TMP = os.path.join(HERE, "runs_test"); os.makedirs(TMP, exist_ok=True)
PASS, FAIL = [], []


def check(name, ok, detail=""):
    (PASS if ok else FAIL).append(name)
    print(("PASS " if ok else "FAIL ") + name + ("  " + detail if detail else ""), flush=True)


def runner(cluster, nup, U, stars=True):
    sub = os.path.join(TMP, "stars" if stars else "allk")
    tag = f"{cluster}_nup{nup}_U{U:g}_k"
    for f in (tag + ".json", tag + ".partial.json"):
        p = os.path.join(sub, f)
        if os.path.exists(p):
            os.rename(p, p + f".prev{int(time.time())}")  # keep old outputs, never delete
    cmd = [sys.executable, os.path.join(ROOT, "run_point_ms.py"), "--cluster", cluster, "--nup", str(nup),
           "--U", str(U), "--runs", sub, "--wait-min", "5"] + ([] if stars else ["--no-stars"])
    subprocess.run(cmd, check=True, capture_output=True)
    return json.load(open(os.path.join(sub, tag + ".json")))


def perm_sign_tables(c, perm, basis):
    idx = np.zeros(len(basis), dtype=np.int64); sg = np.zeros(len(basis))
    newc = np.zeros(len(basis), dtype=np.uint64)
    for n, conf in enumerate(basis):
        occ = [perm[k] for k in range(c.N) if (int(conf) >> k) & 1]
        inv = sum(1 for x in range(len(occ)) for y in range(x + 1, len(occ)) if occ[x] > occ[y])
        sg[n] = -1 if inv % 2 else 1; newc[n] = sum(1 << m for m in occ)
    ii, ok = lookup(basis, newc); assert ok.all(); idx[:] = ii
    return idx, sg


def apply_perm(X, idx, sg):
    Y = np.zeros_like(X); Y[np.ix_(idx, idx)] = sg[:, None] * sg[None, :] * X
    return Y


# 1 geometry
expect_pg = {"18t": ['E', 'C4', 'C2', 'C4^3', 'Mxy', 'Mx', 'My', "Mxy'"], "20t": ['E', 'C4', 'C2', 'C4^3']}
for name, N in (("18t", 18), ("20t", 20)):
    c = Cluster(name)
    check(f"{name} has {N} sites", c.N == N == len(c.sites))
    check(f"{name} 4 distinct neighbours per site", all(len({c.nbr[d][i] for d in c.nbr}) == 4 for i in range(N)))
    check(f"{name} point group {expect_pg[name]}", geometry.point_group(c) == expect_pg[name], str(geometry.point_group(c)))
    ks = c.k_points()
    ok = all(abs(math.remainder(float(k @ a), 2 * math.pi)) < 1e-9 for k in ks for a in (c.a1, c.a2))
    check(f"{name} {len(ks)} K points with K.a in 2piZ", ok and len(ks) == N)
    ks, st = geometry.stars(c)
    check(f"{name} stars partition K", sorted(m for _, ms in st for m in ms) == list(range(N)), str([len(m) for _, m in st]))

# 2 symmetry of H in the full space (U=8, nup=ndn=2 and 3)
rng = np.random.default_rng(0)
for name in ("18t", "20t"):
    c = Cluster(name)
    for nup in (2, 3):
        s = Sector(c, nup, nup, 8.0); basis = configs(c.N, nup)
        X = rng.standard_normal((s.du, s.dd))
        idx, sgn = translation_tables(c, basis)
        err = 0.0
        for g in range(c.N):
            TX = apply_perm(X, idx[g], sgn[g].astype(float))
            err = max(err, np.abs(s.matvec_mat(TX) - apply_perm(s.matvec_mat(X), idx[g], sgn[g].astype(float))).max())
        check(f"{name} nup={nup} [H,T_g]=0 for all {c.N} translations", err < 1e-10, f"max err {err:.1e}")
        err = 0.0
        for op in geometry.point_group(c):
            pi_, ps = perm_sign_tables(c, geometry.site_map(c, op), basis)
            err = max(err, np.abs(s.matvec_mat(apply_perm(X, pi_, ps)) - apply_perm(s.matvec_mat(X), pi_, ps)).max())
        check(f"{name} nup={nup} [H,P]=0 for detected point-group ops", err < 1e-10, f"max err {err:.1e}")

# 3 U = 0 exact
for name in ("18t", "20t"):
    c = Cluster(name); e = np.linalg.eigvalsh(c.single_particle_h())
    for nup in (1, 2, 3, 4):
        r = runner(name, nup, 0.0)
        Eex = 2 * e[:nup].sum()
        check(f"{name} nup={nup} U=0 E0 = free-fermion {Eex:.10f}", abs(r["E0"] - Eex) < 1e-8, f"got {r['E0']:.10f}")
        w = wick_pd(c, nup, nup)
        if w is not None:
            check(f"{name} nup={nup} U=0 P_d = Wick {w:.10f}", abs(r["P_d"] - w) < 1e-8, f"got {r['P_d']:.10f} deg {r['deg']}")
        else:
            print(f"info {name} nup={nup} U=0 open shell (degenerate Fermi level): Wick check n/a, deg={r['deg']}")

# 4 U = 8: momentum runner vs plain full-space solver; 5 stars vs all K
for name in ("18t", "20t"):
    c = Cluster(name); R = c.r_max()
    for nup in (2, 3):
        rk = runner(name, nup, 8.0)
        s = Sector(c, nup, nup, 8.0); E0, w, vecs, how = s.ground_states()
        Pp = float(np.mean([pd_correlation(s, X, R)[1] for X in vecs]))
        check(f"{name} nup={nup} U=8 E0 K vs plain", abs(rk["E0"] - E0) < 2e-8, f"K {rk['E0']:.10f} plain {E0:.10f}")
        check(f"{name} nup={nup} U=8 degeneracy K vs plain", rk["deg"] == len(vecs), f"K {rk['deg']} plain {len(vecs)}")
        check(f"{name} nup={nup} U=8 P_org K vs plain", abs(rk["P_org"] - 2 * Pp) < 1e-8 + 1e-6 * abs(Pp),
              f"K {rk['P_org']:.10f} plain {2 * Pp:.10f}")
        ra = runner(name, nup, 8.0, stars=False)
        ks, st = geometry.stars(c)
        Ek = [sec["E"] for sec in ra["sectors"]]
        spread = max(max(Ek[m] for m in ms) - min(Ek[m] for m in ms) for _, ms in st)
        check(f"{name} nup={nup} U=8 sector energies equal within each star", spread < 1e-8, f"max spread {spread:.1e}")
        check(f"{name} nup={nup} U=8 all-K run = star run", abs(ra["E0"] - rk["E0"]) < 1e-9 and ra["deg"] == rk["deg"]
              and abs(ra["P_org"] - rk["P_org"]) < 1e-8, f"P_org {ra['P_org']:.10f} vs {rk['P_org']:.10f}")
        print(f"info {name} nup={nup} n={2*nup/c.N:.3f} U=8 E0={rk['E0']:.10f} deg={rk['deg']} P_org={rk['P_org']:.8g}")

print("info published 18/20-site U/t=8 energies: none found (Lin PRB 44,7151 (1991) not accessible; "
      "PRB 72,045113 (2005) EPAPS 20-site values are U/t=5,10,20 only) -> no literature comparison made")
print(f"{len(PASS)}/{len(PASS) + len(FAIL)} checks passed")
sys.exit(1 if FAIL else 0)
