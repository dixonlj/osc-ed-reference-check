#!/usr/bin/env python3
"""Optional, off-by-default ED version of OSC estimate_pairing().

Default (OSC_ED_PAIRING unset or 0): the published official formula,
    Index = round(N(E_F) * A_d * 300, 3)
using the caller-supplied A_d. This matches official instrument/osc_instrument.py.
This package does **not** ship that official file (fetch-only if you want a copy).

Opt-in (environment OSC_ED_PAIRING=1, or use_ed=True): A_d is replaced by the
locked 1/N extrapolated amplitude from this package, at the material's band
filling, then clamped at 0. The curves come from
ed_multisize/results/extrapolation_inputs.csv (the stored cluster P_org values
that also produced extrapolation_by_filling.csv). The rule is the one in
PREREGISTRATION.md (sha256 504a9e21…4564): n_eff = n_band if n_band <= 1 else
2 - n_band; linear interpolation of P_org in n on each cluster, with P(0)=0 and
no extrapolation past that cluster's highest computed filling; unweighted least
squares P = A + B/N over covering clusters with N >= 10, if at least 3 cover n,
else the largest covering cluster (flagged unextrapolated); A_d = max(0, A).
Nothing is fitted or tuned to the official 6x6 table or any published value.

This is a calibration option, not a proposed replacement score. In opt-in mode
it scores 20 of the 22 challenge materials as 0. The intervals are wide, and a
U=0 control shows that the sign pattern of the size trend is largely a
finite-cluster shell and shape effect.

Licences: this module is MIT. OSC / C2DB-derived values used at runtime are
CC-BY-NC 4.0 (non-commercial use only).
"""
from __future__ import annotations

import csv
import os
from functools import lru_cache

import numpy as np

SCALE = 300.0
LOCKED_PREREG = "504a9e21047076de4039610f67c502f727b49b69100a74ecd681de44d8f64564"
BIG = ("sqrt10", "3x4", "4x4", "18t", "20t")

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_INPUTS = os.path.join(HERE, "ed_multisize", "results", "extrapolation_inputs.csv")


_ENV_ON = {"1", "true", "yes", "on"}
_ENV_OFF = {"", "0", "false", "no", "off"}


def env_opt_in(name="OSC_ED_PAIRING"):
    """Parse OSC_ED_PAIRING. Returns (enabled: bool, warning_or_None).

    Recognised on: 1, true, yes, on (any case, surrounding spaces ignored).
    Recognised off: unset, empty, 0, false, no, off.
    Any other non-empty value is treated as off and returns a one-line warning
    so it is not silently mistaken for "on".
    """
    raw = os.environ.get(name, "")
    key = str(raw).strip().lower()
    if key in _ENV_ON:
        return True, None
    if key in _ENV_OFF:
        return False, None
    warn = (f"[warn] unrecognised {name}={raw!r}; treating as off. "
            "Use 1/true/yes/on to enable, or 0/false/no/off to disable.")
    return False, warn


def _truthy_env(name, default=""):
    """Back-compat: True only for recognised 'on' values. May print a warning."""
    on, warn = env_opt_in(name)
    if warn:
        print(warn, file=__import__("sys").stderr)
    return on


def lookup_material_n_band(mid, materials_path=None):
    """Canonical n_band (old_n) for one of the 22 package materials, or None."""
    path = materials_path or os.path.join(HERE, "inputs", "materials_canonical.csv")
    if not mid or not os.path.isfile(path):
        return None
    for r in csv.DictReader(open(path)):
        if r.get("osc_id") == mid or r.get("material_id") == mid:
            return float(r["old_n"])
    return None


def official_index(nef, A_d):
    """Published OSC Pairing Index: round(N(E_F) * A_d * 300, 3)."""
    if nef is None or A_d is None:
        raise ValueError("N(E_F) or A_d not available for this material.")
    return round(float(nef) * float(A_d) * SCALE, 3)


def n_eff(n_band):
    """Pre-registered filling map. n_band must be in [0, 2]; otherwise ValueError."""
    if n_band is None:
        raise ValueError("n_band is required")
    x = float(n_band)
    if x < 0 or x > 2:
        raise ValueError(f"n_band={n_band} is outside the documented range [0, 2]")
    return x if x <= 1 else 2 - x


def _interp(xy, n):
    """Linear interpolation in n; xy is a sorted list of (n, P) including (0,0).
    Returns None if n is above the highest computed filling (no extrapolation in n)."""
    nmax = xy[-1][0]
    if n > nmax + 1e-12:
        return None
    xs, ys = zip(*xy)
    return float(np.interp(n, xs, ys))


def _fit(Ns, Ps):
    X = np.column_stack([np.ones(len(Ns)), np.array(Ns, float) ** -1.0])
    coef, *_ = np.linalg.lstsq(X, np.array(Ps, float), rcond=None)
    return float(coef[0]), float(coef[1])


def load_curves_from_csv(inputs_path):
    """Cluster P_org(n) curves from a CSV with columns cluster,N,n,P (P(0)=0 is added)."""
    rows = list(csv.DictReader(open(inputs_path)))
    cs = {}
    for r in rows:
        cl, N, n, P = r["cluster"], int(r["N"]), float(r["n"]), float(r["P"])
        cs.setdefault(cl, {"N": N, "xy": [(0.0, 0.0)]})["xy"].append((n, P))
    for c in cs.values():
        c["xy"].sort()
        c["nmax"] = c["xy"][-1][0]
    return cs


def _estimate_from_curves(cs, n_band):
    """Same locked rule as ed_multisize/extrapolate.py:estimate()."""
    n = n_band if n_band <= 1 else 2 - n_band
    big = {}
    for cl, c in cs.items():
        if cl not in BIG:
            continue
        v = _interp(c["xy"], n)
        if v is not None:
            big[cl] = v
    out = dict(n_eff=n, n_used=len(big), A=None, A_d=None, flag="")
    if not big:
        out["flag"] = "no cluster covers n"
        return out
    Ns = [cs[k]["N"] for k in big]
    Ps = list(big.values())
    if len(big) >= 3:
        A, B = _fit(Ns, Ps)
        out.update(A=A, B=B, flag="")
    else:
        largest = max(big, key=lambda k: cs[k]["N"])
        out.update(A=big[largest], B=None,
                   flag=f"unextrapolated (only {len(big)} cluster(s) N>=10 cover n; value of {largest})")
    out["A_d"] = max(0.0, out["A"])
    if out["A"] < 0:
        out["flag"] = (out["flag"] + "; " if out["flag"] else "") + "clamped to 0 (negative)"
    return out


@lru_cache(maxsize=1)
def _package_estimate():
    """The package's own extrapolate.estimate, used so A_d matches the locked CSVs bit for bit."""
    import importlib.util
    path = os.path.join(HERE, "ed_multisize", "extrapolate.py")
    spec = importlib.util.spec_from_file_location("ed_multisize_extrapolate", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    cs = mod.curves(mod.load_points())
    return mod.estimate, cs


def locked_result(n_band, inputs_path=None):
    """Locked 1/N extrapolation at n_band. Returns dict with n_eff, A, A_d, flag, n_used.

    Documented edge cases:
    - n_band in [0, 2]: mapped to n_eff in [0, 1] by particle-hole (n_eff = n_band or 2-n_band).
    - n_band outside [0, 2]: ValueError.
    - If n_eff is above every stored cluster's highest filling: A and A_d are None,
      flag='no cluster covers n'. With this package's curves, that does not happen
      for n_eff in [0, 1] (the 10-, 12- and 16-site clusters all reach n=1).
    - If 1 or 2 clusters with N>=10 cover n: A is the value of the largest covering
      cluster, flagged 'unextrapolated' (pre-registered fallback).
    - If the intercept A is negative: A_d is 0 and the result is flagged
      'clamped to 0 (negative)'.

    By default the cluster curves are the same stored run files that
    extrapolate.py used (so official-filling and material scores match the
    locked CSVs exactly). Pass inputs_path to use a CSV instead (tests).
    """
    n = n_eff(n_band)  # validates [0, 2]
    if inputs_path:
        e = _estimate_from_curves(load_curves_from_csv(inputs_path), float(n_band))
    else:
        estimate, cs = _package_estimate()
        e = estimate(cs, float(n_band))
    return dict(n_band=float(n_band), n_eff=n,
                n_used=e.get("n_used"), A=e.get("A"), A_d=e.get("A_d"),
                flag=e.get("flag") or "", B=e.get("B"))


def locked_A_d(n_band, inputs_path=None):
    """Clamped locked A_d at n_band, or None if no cluster covers n."""
    return locked_result(n_band, inputs_path)["A_d"]


def estimate_pairing(nef, A_d=None, n_band=None, use_ed=None, inputs_path=None):
    """OSC Pairing Index. Official behaviour unless opt-in is on.

    Parameters
    ----------
    nef : float
        N(E_F) in states/eV/atom.
    A_d : float, optional
        Official (or caller-supplied) d-wave amplitude. Required when opt-in is off.
    n_band : float, optional
        Canonical band filling. Required when opt-in is on. If omitted, the
        environment variable OSC_ED_N_BAND is used if set.
    use_ed : bool, optional
        If None (default), read OSC_ED_PAIRING from the environment (1/true/yes/on
        to enable). If False, official behaviour. If True, locked ED A_d.
    inputs_path : str, optional
        Path to extrapolation_inputs.csv (tests / alternate data).

    Returns
    -------
    float
        Pairing Index, rounded to 3 decimals, same scale as the official instrument.
    """
    if use_ed is None:
        use_ed, warn = env_opt_in("OSC_ED_PAIRING")
        if warn:
            print(warn, file=__import__("sys").stderr)
    if not use_ed:
        return official_index(nef, A_d)
    if n_band is None:
        env_n = os.environ.get("OSC_ED_N_BAND")
        if env_n in (None, ""):
            mid = os.environ.get("OSC_ED_MATERIAL_ID")
            n_band = lookup_material_n_band(mid) if mid else None
        else:
            n_band = float(env_n)
        if n_band is None:
            raise ValueError("n_band is required when OSC_ED_PAIRING is on "
                             "(argument, OSC_ED_N_BAND, or a known --material_id / OSC_ED_MATERIAL_ID)")
    res = locked_result(n_band, inputs_path)
    if res["A_d"] is None:
        raise ValueError(f"no locked A_d at n_band={n_band} (n_eff={res['n_eff']}): {res['flag']}")
    return official_index(nef, res["A_d"])


if __name__ == "__main__":
    import argparse, json
    ap = argparse.ArgumentParser(description="Optional ED estimate_pairing (off by default)")
    ap.add_argument("--nef", type=float, required=True)
    ap.add_argument("--A_d", type=float, default=None, help="official / caller A_d (used when opt-in is off)")
    ap.add_argument("--n_band", type=float, default=None, help="canonical filling (required when opt-in is on)")
    ap.add_argument("--ed", action="store_true", help="force opt-in (same as OSC_ED_PAIRING=1)")
    a = ap.parse_args()
    print(json.dumps({
        "index": estimate_pairing(a.nef, a.A_d, n_band=a.n_band, use_ed=True if a.ed else None),
        "use_ed": True if a.ed else _truthy_env("OSC_ED_PAIRING"),
        "locked": locked_result(a.n_band) if (a.ed or _truthy_env("OSC_ED_PAIRING")) and a.n_band is not None else None,
    }, indent=1))
