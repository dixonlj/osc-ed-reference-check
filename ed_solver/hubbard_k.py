"""Translation-symmetric (crystal momentum K) Hubbard ED, spin-factorized storage.

A K-sector state satisfies T_g|psi> = exp(i K.r_g)|psi> for every torus translation g.
Its amplitudes obey psi(g a, g b) = exp(-i K.r_g) s_g(a) s_g(b) psi(a, b), where s_g is the
fermion reordering sign of translating a configuration. We store X[rho, b] = psi(r_rho, b)
for one representative r_rho of every up-spin translation orbit and every down config b.

Norm: ||psi||^2 = sum_rho w_rho ||X[rho,:]||^2 with w_rho = N / |Stab(r_rho)|. eigsh works on
Z = sqrt(w) X with the Hermitian operator Z -> sqrt(w) P H P (Z / sqrt(w)), where P projects
rows with a non-trivial stabilizer onto the K-consistent subspace (it is the identity on the
others). Spurious eigenvalue 0 lives in ker(P); ground-state energies here are negative.

P_d(R) for a momentum eigenstate is translation invariant, so the all-pairs average equals
(1/m) sum_{j in shell(0)} <Delta(j)^dag Delta(0)>; it is evaluated row by row without ever
building the full real-space vector.
"""
import math
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, eigsh

from hubbard import (SEED, DEG_TOL, configs, lookup, popcount, hopping_matrix,
                     pair_terms, deflate_complete)
from lanczos import lanczos_lowest


def translation_tables(cluster, basis):
    """For every translation g: index of g(c) and sign s_g(c) for each config c."""
    N = cluster.N
    G = cluster.trans  # G[g][k] = image of site k
    nb = len(basis)
    idx = np.zeros((N, nb), dtype=np.int64)
    sgn = np.zeros((N, nb), dtype=np.int8)
    occ = [[k for k in range(N) if (int(c) >> k) & 1] for c in basis]
    for g in range(N):
        img = G[g]
        newc = np.zeros(nb, dtype=np.uint64)
        for n, sites in enumerate(occ):
            mapped = [int(img[k]) for k in sites]
            inv = 0
            for x in range(len(mapped)):
                for y in range(x + 1, len(mapped)):
                    if mapped[x] > mapped[y]:
                        inv += 1
            sgn[g, n] = -1 if inv % 2 else 1
            newc[n] = sum(1 << m for m in mapped)
        ii, ok = lookup(basis, newc)
        assert ok.all()
        idx[g] = ii
    return idx, sgn


class KSector:
    def __init__(self, cluster, nup, ndn, U, K, t=1.0, tables=None):
        assert nup == ndn, "spin-factorized K solver assumes N_up == N_dn"
        self.c, self.nup, self.ndn, self.U, self.t = cluster, nup, ndn, U, t
        self.K = np.array(K, dtype=float)
        N = cluster.N
        self.bu = configs(N, nup)
        self.bd = self.bu
        self.du = self.dd = len(self.bu)
        if tables is None:
            tables = translation_tables(cluster, self.bu)
        self.tidx, self.tsgn = tables
        r = np.array(cluster.sites, dtype=float)
        self.phase = np.exp(-1j * (r @ self.K))  # exp(-i K.r_g)
        # inverse translation for each g
        ginv = np.zeros(N, dtype=np.int64)
        for g in range(N):
            for h in range(N):
                if cluster.trans[g][cluster.trans[h][0]] == 0 and all(
                        cluster.trans[h][cluster.trans[g][k]] == k for k in range(N)):
                    ginv[g] = h; break
        self.ginv = ginv
        # orbits
        rep_of = -np.ones(self.du, dtype=np.int64)
        g_of = np.zeros(self.du, dtype=np.int64)
        reps, stabs = [], []
        for a in range(self.du):
            if rep_of[a] >= 0:
                continue
            reps.append(a)
            stab = []
            for g in range(N):
                b = self.tidx[g, a]
                if b == a:
                    stab.append(g)
                if rep_of[b] < 0:
                    rep_of[b] = len(reps) - 1
                    g_of[b] = g
            stabs.append(stab)
        self.reps = np.array(reps, dtype=np.int64)
        self.stabs = stabs
        self.rep_of, self.g_of = rep_of, g_of
        self.R = len(reps)
        self.w = np.array([N / len(s) for s in stabs])
        self.sqw = np.sqrt(self.w)[:, None]
        # per-translation down-index permutation and sign:  pinv_g[b] = idx(g^-1 b), sig_g[b] = s_g(g^-1 b)
        self.pinv = np.array([self.tidx[ginv[g]] for g in range(N)])
        self.sig = np.array([self.tsgn[g][self.pinv[g]] for g in range(N)], dtype=np.float64)
        # up-hopping, grouped by translation
        Tu = hopping_matrix(cluster, self.bu, t).tocsr()
        groups = {}
        for rho, r0 in enumerate(self.reps):
            row = Tu.getrow(r0)
            for a, val in zip(row.indices, row.data):
                ra, g = rep_of[a], g_of[a]
                coef = val * self.phase[g] * self.tsgn[g, self.reps[ra]]
                groups.setdefault(g, ([], [], []))
                groups[g][0].append(rho); groups[g][1].append(ra); groups[g][2].append(coef)
        self.Mg = {g: sp.csr_matrix((np.array(d, dtype=complex), (np.array(r_), np.array(c_))),
                                    shape=(self.R, self.R)) for g, (r_, c_, d) in groups.items()}
        self.Td = Tu
        self.D = popcount(self.bu[self.reps][:, None] & self.bd[None, :]).astype(np.float64)
        # projector data for rows with nontrivial stabilizer
        self.prow = [(rho, s) for rho, s in enumerate(stabs) if len(s) > 1]
        self.dim = self.R * self.dd

    def project(self, X):
        for rho, S in self.prow:
            r0 = self.reps[rho]
            acc = np.zeros(self.dd, dtype=complex)
            for h in S:
                acc += self.phase[h] * self.tsgn[h, r0] * self.sig[h] * X[rho, self.pinv[h]]
            X[rho] = acc / len(S)
        return X

    def H_X(self, X):
        out = np.zeros_like(X)
        for g, M in self.Mg.items():
            if g == 0 and np.all(self.pinv[0] == np.arange(self.dd)):
                out += M @ X
            else:
                out += M @ (X[:, self.pinv[g]] * self.sig[g])
        out += (self.Td @ X.T).T
        if self.U != 0.0:
            out += self.U * (self.D * X)
        return out

    def op_matvec(self, z):
        X = z.reshape(self.R, self.dd) / self.sqw
        X = self.project(X.copy())
        Y = self.project(self.H_X(X))
        return (Y * self.sqw).ravel()

    def lowest(self, k=2, tol=1e-12, ncv=None, deflate=True):
        rng = np.random.default_rng(SEED)
        v0 = rng.standard_normal(self.dim) + 1j * rng.standard_normal(self.dim)
        v0 = (self.project(v0.reshape(self.R, self.dd) / self.sqw) * self.sqw).ravel()
        op = LinearOperator((self.dim, self.dim), dtype=np.complex128, matvec=self.op_matvec)
        k = min(k, self.dim - 2)
        w, v = eigsh(op, k=k, which="SA", v0=v0, tol=tol, ncv=min(ncv or max(2 * k + 2, 8), self.dim - 1),
                     maxiter=100000)
        order = np.argsort(w); w, v = w[order], v[:, order]
        E0 = w[0]
        found = [v[:, s] for s in np.nonzero(np.abs(w - E0) < DEG_TOL)[0]]
        if deflate:
            found = deflate_complete(op, E0, found, self.dim, np.complex128, rng, tol, ncv=min(ncv or 8, 12))
        return E0, w, [z.reshape(self.R, self.dd) / self.sqw for z in found]


    # ---- transposed layout (Y = X.T, shape dd x R, C-contiguous) for speed ----
    def _prep_T(self):
        if hasattr(self, "_DT"):
            return
        self._DT = np.ascontiguousarray(self.D.T)
        self._sqwT = np.sqrt(self.w)[None, :]
        self._buf = None

    def project_T(self, Y):
        for rho, S in self.prow:
            r0 = self.reps[rho]
            acc = np.zeros(self.dd, dtype=complex)
            for h in S:
                acc += self.phase[h] * self.tsgn[h, r0] * self.sig[h] * Y[self.pinv[h], rho]
            Y[:, rho] = acc / len(S)
        return Y

    def H_Y(self, Y):
        out = self.Td @ Y
        if self.U != 0.0:
            out += self.U * (self._DT * Y)
        if self._buf is None:
            self._buf = np.empty_like(Y)
        buf = self._buf
        for g, M in self.Mg.items():
            np.take(Y, self.pinv[g], axis=0, out=buf)
            buf *= self.sig[g][:, None]
            out += (M @ buf.T).T
        return out

    def op_Z(self, z):
        """Hermitian operator on Z = sqrt(w)-scaled storage in transposed layout."""
        self._prep_T()
        Y = z.reshape(self.dd, self.R) / self._sqwT
        Y = self.project_T(Y)
        out = self.project_T(self.H_Y(Y))
        out *= self._sqwT
        return out.ravel()

    def energy_lanczos(self, log=None):
        """Lowest energy of this sector only (one Lanczos pass, no eigenvector)."""
        self._prep_T()
        rng = np.random.default_rng(SEED)
        v = rng.standard_normal(self.dim) + 1j * rng.standard_normal(self.dim)
        Y = self.project_T(v.reshape(self.dd, self.R) / self._sqwT) * self._sqwT
        E, _, _, nmv = lanczos_lowest(self.op_Z, Y.ravel(), log=log, energy_only=True)
        return E, nmv

    def lowest_lanczos(self, deflate=True, res_tol=1e-8, log=None, max_deg=8):
        """Lowest state(s) of this sector by low-memory Lanczos (+ deflation for degeneracy).

        Returns (E0, [X matrices in R x dd layout], residuals, matvecs)."""
        self._prep_T()
        rng = np.random.default_rng(SEED)
        def start():
            v = rng.standard_normal(self.dim) + 1j * rng.standard_normal(self.dim)
            Y = self.project_T(v.reshape(self.dd, self.R) / self._sqwT) * self._sqwT
            return Y.ravel()
        E0, z, res, nmv = lanczos_lowest(self.op_Z, start(), res_tol=res_tol, log=log)
        found, ress = [z], [res]
        if deflate:
            while len(found) < max_deg:
                E1, z1, r1, n1 = lanczos_lowest(self.op_Z, start(), res_tol=res_tol, Q=found, log=log)
                nmv += n1
                if abs(E1 - E0) < DEG_TOL:
                    for q in found:
                        z1 = z1 - np.vdot(q, z1) * q
                    z1 /= np.linalg.norm(z1)
                    found.append(z1); ress.append(r1)
                else:
                    break
        Xs = [np.ascontiguousarray((zz.reshape(self.dd, self.R) / self._sqwT).T) for zz in found]
        return E0, Xs, ress, nmv

    def residual(self, X, E):
        Xp = self.project(X.copy())
        r = self.project(self.H_X(Xp)) - E * Xp
        return float(np.sqrt(np.sum(self.w[:, None] * np.abs(r) ** 2)))

    def norm(self, X):
        return float(np.sqrt(np.sum(self.w[:, None] * np.abs(X) ** 2)))

    def row(self, X, a):
        """Full real-space row psi(a, :) for an arbitrary up config index a."""
        rho, g = self.rep_of[a], self.g_of[a]
        return self.phase[g] * self.tsgn[g, self.reps[rho]] * self.sig[g] * X[rho, self.pinv[g]]

    def pd(self, X, R=None):
        c = self.c
        R = c.r_max() if R is None else R
        shell = c.shell(0, R)
        N = c.N
        bur = configs(N, self.nup - 1)
        bdr = bur
        # down annihilation maps on the reduced basis: for each q, full index of b' | q and sign
        dmap = []
        for q in range(N):
            bq = np.uint64(1 << q)
            ok = (bdr & bq) == 0
            full, ok2 = lookup(self.bd, bdr | bq)
            s = 1 - 2 * (popcount(bdr & np.uint64((1 << q) - 1)) % 2)
            dmap.append((np.nonzero(ok)[0], full[ok], s[ok].astype(float)))
        up_idx = {}
        sgn_glob = -1.0 if self.nup % 2 else 1.0
        T0 = pair_terms(c, 0)
        Tj = {j: pair_terms(c, j) for j in shell}
        total = {j: 0.0 for j in shell}
        for ap in bur:
            ap = int(ap)
            rows_cache = {}

            def phi_row(terms):
                out = np.zeros(len(bdr), dtype=complex)
                for coef, p, q in terms:
                    if (ap >> p) & 1:
                        continue
                    a_full = ap | (1 << p)
                    if a_full not in up_idx:
                        up_idx[a_full] = int(lookup(self.bu, np.array([a_full], dtype=np.uint64))[0][0])
                    ai = up_idx[a_full]
                    if ai not in rows_cache:
                        rows_cache[ai] = self.row(X, ai)
                    su = -1.0 if bin(ap & ((1 << p) - 1)).count("1") % 2 else 1.0
                    pos, full, s = dmap[q]
                    out[pos] += coef * sgn_glob * su * s * rows_cache[ai][full]
                return out

            f0 = phi_row(T0)
            if not np.any(f0):
                continue
            for j in shell:
                total[j] += float(np.real(np.vdot(phi_row(Tj[j]), f0)))
        return R, sum(total.values()) / len(shell)


def c4v_stars(cluster):
    """Group K points into stars under the square point group (only for LxL axis-aligned clusters)."""
    ks = cluster.k_points()
    tw = 2 * math.pi

    def key(k):
        return (round((k[0] % tw) / tw * 1e6) % 1000000, round((k[1] % tw) / tw * 1e6) % 1000000)
    keys = {key(k): i for i, k in enumerate(ks)}
    ops = [lambda x, y: (x, y), lambda x, y: (-y, x), lambda x, y: (-x, -y), lambda x, y: (y, -x),
           lambda x, y: (y, x), lambda x, y: (-x, y), lambda x, y: (x, -y), lambda x, y: (-y, -x)]
    seen, stars = set(), []
    for i, k in enumerate(ks):
        if i in seen:
            continue
        members = sorted({keys[key(np.array(op(*k)))] for op in ops})
        seen.update(members)
        stars.append((i, members))
    return ks, stars
