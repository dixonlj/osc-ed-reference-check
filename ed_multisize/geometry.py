"""Tilted 18- and 20-site clusters + automatic point-group detection and K-star grouping.

Registers the new clusters into lattice.CLUSTERS at import time (in memory only; the files in
method/ed_solver are not edited, because the running 4x4 check imports them).
  18t: a1=(3,3),  a2=(-3,3)  -> full C4v expected
  20t: a1=(4,2),  a2=(-2,4)  -> C4 only expected (x<->y mirror maps (4,2) to (2,4), not a lattice vector)
"""
import math, os, sys
import numpy as np

ED = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ed_solver")
sys.path.insert(0, os.path.abspath(ED))
import lattice  # noqa: E402

NEW = {"18t": ((3, 3), (-3, 3)), "20t": ((4, 2), (-2, 4))}
for _k, _v in NEW.items():
    lattice.CLUSTERS.setdefault(_k, _v)

# the 8 operations of the square point group D4 (C4v), as integer 2x2 matrices
D4 = {
    "E": ((1, 0), (0, 1)), "C4": ((0, -1), (1, 0)), "C2": ((-1, 0), (0, -1)), "C4^3": ((0, 1), (-1, 0)),
    "Mxy": ((0, 1), (1, 0)), "Mx": ((-1, 0), (0, 1)), "My": ((1, 0), (0, -1)), "Mxy'": ((0, -1), (-1, 0)),
}


def point_group(cluster):
    """Names of D4 ops that map the superlattice onto itself (so they are symmetries of H)."""
    A = np.array([cluster.a1, cluster.a2], dtype=float).T
    Ainv = np.linalg.inv(A)
    ok = []
    for name, M in D4.items():
        M = np.array(M, dtype=float)
        coeff = Ainv @ M @ A
        if np.allclose(coeff, np.round(coeff), atol=1e-9):
            ok.append(name)
    return ok


def site_map(cluster, op):
    """Site permutation induced by a point-group op (rotation about site (0,0))."""
    M = np.array(D4[op], dtype=int)
    return [cluster.index[cluster.reduce(tuple(M @ np.array(r)))] for r in cluster.sites]


def stars(cluster):
    """K points grouped into stars under the detected point group. Returns (ks, [(rep, members)])."""
    ks = cluster.k_points()
    tw = 2 * math.pi

    def key(k):
        return (round((k[0] % tw) / tw * 1e6) % 1000000, round((k[1] % tw) / tw * 1e6) % 1000000)
    keys = {key(k): i for i, k in enumerate(ks)}
    ops = [np.array(D4[o], dtype=float) for o in point_group(cluster)]
    seen, out = set(), []
    for i, k in enumerate(ks):
        if i in seen:
            continue
        members = sorted({keys[key(M @ k)] for M in ops})
        seen.update(members)
        out.append((i, members))
    return ks, out


def describe(name):
    c = lattice.Cluster(name)
    ks, st = stars(c)
    R = c.r_max()
    return dict(cluster=name, N=c.N, a1=tuple(int(x) for x in c.a1), a2=tuple(int(x) for x in c.a2),
                point_group=point_group(c), n_k=len(ks), n_stars=len(st), star_sizes=[len(m) for _, m in st],
                R_max=R, shell_size=len(c.shell(0, R)))


if __name__ == "__main__":
    for n in ["2x4", "sqrt8", "sqrt10", "3x4", "4x4", "18t", "20t"]:
        print(describe(n))
