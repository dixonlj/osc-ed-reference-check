"""Real-space Hubbard ED in fixed (N_up, N_dn) sectors, no spatial symmetry.

Basis |a,b> = (prod_{i in a, ascending} c^dag_{i,up}) (prod_{j in b, ascending} c^dag_{j,dn}) |0>.
A state is stored as a matrix X[a, b] (up index, down index).

H = -t sum_i sum_{d in +x,+y} sum_s (c^dag_{j s} c_{i s} + h.c.) + U sum_i n_{i up} n_{i dn}

Pair operator (Ku-ku, method/ref/multisize_ED_table_kuku.md):
  b_ij      = (c_{i up} c_{j dn} + c_{j up} c_{i dn}) / sqrt(2)
  Delta_d(i) = 1/2 sum_{delta in +-x, +-y} g(delta) b_{i, i+delta},  g(+-x)=+1, g(+-y)=-1
  P_d(R)    = < Delta_d(j)^dag Delta_d(i) > averaged over all ordered pairs with
              minimum-image separation R
  P_org     = 2 * P_d   (organizers' normalization, stated by Ku-ku)
Degenerate ground states (|E - E0| < DEG_TOL) are averaged with equal weight.
"""
import itertools, math
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.sparse.linalg import LinearOperator, eigsh

from lattice import GFORM

import os as _os
SEED = int(_os.environ.get("OSC_ED_SEED", "0"))  # package: default 0 (frozen); override only for independent-seed re-runs
DEG_TOL = 1e-7
DENSE_MAX = 6000


def popcount(x):
    return np.bitwise_count(np.asarray(x, dtype=np.uint64)).astype(np.int64)


def configs(N, n):
    out = [sum(1 << i for i in c) for c in itertools.combinations(range(N), n)]
    return np.array(sorted(out), dtype=np.uint64)


def lookup(basis, vals):
    idx = np.searchsorted(basis, vals)
    idx = np.minimum(idx, len(basis) - 1)
    ok = basis[idx] == vals
    return idx, ok


def hopping_matrix(cluster, basis, t=1.0):
    """Single-species hopping matrix T[new, old] = <new| -t sum c^dag_j c_i |old>."""
    rows, cols, data = [], [], []
    ar = np.arange(len(basis))
    for i, j in cluster.hop_bonds():
        bi, bj = np.uint64(1 << i), np.uint64(1 << j)
        sel = ((basis & bi) != 0) & ((basis & bj) == 0)
        if not sel.any():
            continue
        c = basis[sel]
        new = c ^ bi ^ bj
        lo, hi = min(i, j), max(i, j)
        between = np.uint64(((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
        sign = 1 - 2 * (popcount(c & between) % 2)
        idx, ok = lookup(basis, new)
        assert ok.all()
        rows.append(idx); cols.append(ar[sel]); data.append(-t * sign.astype(float))
    if not rows:
        return sp.csr_matrix((len(basis), len(basis)))
    T = sp.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(len(basis), len(basis))).tocsr()
    T.sum_duplicates()
    return T


def annihilation_matrix(basis_from, basis_to, p):
    """A[a', a] = sign for c_p |a> = sign |a'>, sign = (-1)^{# occupied below p}."""
    bp = np.uint64(1 << p)
    sel = (basis_from & bp) != 0
    a = basis_from[sel]
    new = a ^ bp
    sign = 1 - 2 * (popcount(a & np.uint64((1 << p) - 1)) % 2)
    idx, ok = lookup(basis_to, new)
    assert ok.all()
    return sp.csr_matrix((sign.astype(float), (idx, np.nonzero(sel)[0])),
                         shape=(len(basis_to), len(basis_from)))


def pair_terms(cluster, i):
    """Delta_d(i) as a list of (coef, p_up, q_dn): sum coef * c_{p up} c_{q dn}."""
    terms = []
    pref = 0.5 / math.sqrt(2.0)
    for d, g in GFORM.items():
        j = cluster.nbr[d][i]
        terms.append((pref * g, i, j))
        terms.append((pref * g, j, i))
    return terms


class Sector:
    def __init__(self, cluster, nup, ndn, U, t=1.0):
        self.c = cluster
        self.nup, self.ndn, self.U, self.t = nup, ndn, U, t
        N = cluster.N
        self.bu = configs(N, nup)
        self.bd = configs(N, ndn)
        self.du, self.dd = len(self.bu), len(self.bd)
        self.dim = self.du * self.dd
        self.Tu = hopping_matrix(cluster, self.bu, t)
        self.Td = self.Tu if ndn == nup else hopping_matrix(cluster, self.bd, t)
        self.D = popcount(self.bu[:, None] & self.bd[None, :]).astype(np.float64) if self.dim <= 5e7 else None

    def matvec_mat(self, X):
        out = self.Tu @ X
        out += (self.Td @ X.T).T
        if self.U != 0.0:
            out += self.U * (self.D * X)
        return out

    def sparse_H(self):
        Iu, Id = sp.identity(self.du, format="csr"), sp.identity(self.dd, format="csr")
        H = sp.kron(self.Tu, Id, format="csr") + sp.kron(Iu, self.Td, format="csr")
        if self.U != 0.0:
            H = H + sp.diags(self.U * self.D.ravel())
        return H.tocsr()

    def ground_states(self, k_init=6, tol=1e-12, ncv=None):
        """Return (E0, energies_found, [ground-state matrices], method)."""
        if self.dim <= DENSE_MAX:
            w, v = sla.eigh(self.sparse_H().toarray())
            E0 = w[0]
            sel = np.nonzero(np.abs(w - E0) < DEG_TOL)[0]
            return E0, w[:max(len(sel) + 2, 6)], [v[:, s].reshape(self.du, self.dd) for s in sel], "dense"
        k = min(k_init, self.dim - 2)
        rng = np.random.default_rng(SEED)
        v0 = rng.standard_normal(self.dim)
        op = LinearOperator((self.dim, self.dim), dtype=np.float64,
                            matvec=lambda x: self.matvec_mat(x.reshape(self.du, self.dd)).ravel())
        nc = ncv or max(2 * k + 1, 20)
        w, v = eigsh(op, k=k, which="SA", v0=v0, tol=tol, ncv=min(nc, self.dim - 1), maxiter=100000)
        order = np.argsort(w); w, v = w[order], v[:, order]
        E0 = w[0]
        found = [v[:, s] for s in np.nonzero(np.abs(w - E0) < DEG_TOL)[0]]
        found = deflate_complete(op, E0, found, self.dim, np.float64, rng, tol)
        return E0, w, [x.reshape(self.du, self.dd) for x in found], "eigsh+deflation"

    def residual(self, X, E):
        return float(np.linalg.norm(self.matvec_mat(X) - E * X))


def deflate_complete(op, E0, found, dim, dtype, rng, tol, shift=50.0, max_extra=64, ncv=20):
    """Lanczos with deflation: add shift*|v><v| for every found ground state and look
    for another state with |E - E0| < DEG_TOL, until none remains (Ku-ku procedure)."""
    found = list(found)
    # orthonormalize (Gram-Schmidt, twice)
    def ortho(vs):
        out = []
        for x in vs:
            y = x.copy()
            for _ in range(2):
                for q in out:
                    y -= np.vdot(q, y) * q
            nrm = np.linalg.norm(y)
            if nrm > 1e-8:
                out.append(y / nrm)
        return out
    found = ortho(found)
    for _ in range(max_extra):
        Q = list(found)
        def mv(x, Q=Q):
            y = op.matvec(x)
            for q in Q:
                y = y + shift * np.vdot(q, x) * q
            return y
        dop = LinearOperator((dim, dim), dtype=dtype, matvec=mv)
        v0 = rng.standard_normal(dim)
        if np.iscomplexobj(np.zeros(1, dtype=dtype)):
            v0 = v0 + 1j * rng.standard_normal(dim)
        w, v = eigsh(dop, k=1, which="SA", v0=v0.astype(dtype), tol=tol, ncv=min(ncv, dim - 1), maxiter=100000)
        if abs(w[0] - E0) < DEG_TOL:
            found = ortho(found + [v[:, 0]])
        else:
            break
    return found


def pd_correlation(sector, X, R=None):
    """Return (R, P_d(R)) for one normalized state X (Ku-ku normalization, not x2)."""
    c = sector.c
    R = c.r_max() if R is None else R
    if sector.nup == 0 or sector.ndn == 0:
        return R, 0.0
    bur = configs(c.N, sector.nup - 1)
    bdr = configs(c.N, sector.ndn - 1)
    Au = [annihilation_matrix(sector.bu, bur, p) for p in range(c.N)]
    Ad = Au if (sector.ndn == sector.nup) else [annihilation_matrix(sector.bd, bdr, q) for q in range(c.N)]
    sgn = -1.0 if sector.nup % 2 else 1.0  # c_dn passes all up operators

    def phi(i):
        out = None
        for coef, p, q in pair_terms(c, i):
            term = coef * sgn * (Au[p] @ (Ad[q] @ X.T).T)
            out = term if out is None else out + term
        return out

    phis = [phi(i) for i in range(c.N)]
    tot, cnt = 0.0, 0
    for i in range(c.N):
        for j in c.shell(i, R):
            tot += float(np.real(np.vdot(phis[j], phis[i])))
            cnt += 1
    return R, tot / cnt


def wick_pd(cluster, nup, ndn, R=None, t=1.0):
    """U=0 closed-shell Slater determinant value of P_d(R) (Ku-ku normalization).

    Returns None if the Fermi level is degenerate (open shell)."""
    R = cluster.r_max() if R is None else R
    h = cluster.single_particle_h(t)
    e, v = np.linalg.eigh(h)
    G = {}
    for s, n in (("up", nup), ("dn", ndn)):
        if 0 < n < cluster.N and abs(e[n] - e[n - 1]) < 1e-9:
            return None
        occ = v[:, :n]
        G[s] = occ @ occ.T  # G[x, y] = <c^dag_x c_y>
    tot, cnt = 0.0, 0
    for i in range(cluster.N):
        Ti = pair_terms(cluster, i)
        for j in cluster.shell(i, R):
            Tj = pair_terms(cluster, j)
            val = 0.0
            for ci, p, q in Ti:
                for cj, pp, qq in Tj:
                    val += ci * cj * G["up"][pp, p] * G["dn"][qq, q]
            tot += val; cnt += 1
    return tot / cnt
