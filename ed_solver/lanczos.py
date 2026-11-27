"""Low-memory Lanczos for the lowest eigenpair of a Hermitian matrix-free operator.

Used for large 4x4 momentum sectors where ARPACK's Krylov basis would not fit the memory
budget. Deterministic: the start vector comes from numpy default_rng(seed). Two passes
(first: tridiagonal coefficients; second: rebuild the Ritz vector), then restarts from the
Ritz vector until the residual ||H x - E x|| is below `res_tol`. Optional deflation vectors Q
are handled by adding shift * |q><q| to the operator (the Ku-ku "Lanczos plus deflation").
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal


def _shifted(matvec, Q, shift):
    if not Q:
        return matvec
    def mv(x):
        y = matvec(x)
        for q in Q:
            y += shift * np.vdot(q, x) * q
        return y
    return mv


def lanczos_lowest(matvec, v0, max_iter=400, e_tol=1e-13, res_tol=1e-8, max_restarts=8,
                   Q=(), shift=50.0, project=None, log=None, energy_only=False):
    mv = _shifted(matvec, list(Q), shift)
    x = v0 / np.linalg.norm(v0)
    total_mv = 0
    E = None
    for restart in range(max_restarts + 1):
        # pass 1
        alphas, betas = [], []
        v_prev = None
        v = x.copy()
        beta = 0.0
        e_old = None
        for it in range(max_iter):
            w = mv(v); total_mv += 1
            if project is not None:
                w = project(w)
            a = float(np.real(np.vdot(v, w)))
            w -= a * v
            if v_prev is not None:
                w -= beta * v_prev
            alphas.append(a)
            b = float(np.linalg.norm(w))
            if it % 5 == 4 or b < 1e-12:
                ev = eigh_tridiagonal(np.array(alphas), np.array(betas), select="i", select_range=(0, 0),
                                      eigvals_only=True)[0]
                if e_old is not None and abs(ev - e_old) < e_tol * max(1.0, abs(ev)):
                    break
                e_old = ev
            if b < 1e-12:
                break
            betas.append(b)
            v_prev, v = v, w / b
            beta = b
        m = len(alphas)
        bt = np.array(betas[:m - 1])
        evals, evecs = eigh_tridiagonal(np.array(alphas), bt, select="i", select_range=(0, 0))
        E = float(evals[0]); c = evecs[:, 0]
        if energy_only:
            if log:
                log(f"lanczos energy-only m={m} E={E:.12f} matvecs={total_mv}")
            return E, None, None, total_mv
        # pass 2: rebuild the Ritz vector
        v_prev = None
        v = x.copy()
        y = c[0] * v
        for it in range(m - 1):
            w = mv(v); total_mv += 1
            if project is not None:
                w = project(w)
            w -= alphas[it] * v
            if v_prev is not None:
                w -= betas[it - 1] * v_prev
            w /= betas[it]
            v_prev, v = v, w
            y += c[it + 1] * v
        y /= np.linalg.norm(y)
        r = mv(y); total_mv += 1
        if project is not None:
            r = project(r)
        Ey = float(np.real(np.vdot(y, r)))
        res = float(np.linalg.norm(r - Ey * y))
        if log:
            log(f"lanczos restart={restart} m={m} E={Ey:.12f} res={res:.2e} matvecs={total_mv}")
        x = y
        E = Ey
        if res < res_tol:
            break
    return E, x, res, total_mv
