#!/usr/bin/env python3
"""Hook the official osc_instrument CLI (fetch-only; not shipped in this package).

This file is a wrapper, not a copy of official code. Run ./fetch_official.sh
(pinned dataset revision 346f184d; osc_instrument.py sha256
6b774e529cb8d50cb68a695f5ba6bdbdcda290031833ecc21533d6a1610b5236), or set
OSC_INSTRUMENT / --instrument to a fetched copy.

Default (OSC_ED_PAIRING unset, 0, false, no, or off): the wrapper runs the
official script and does not replace estimate_pairing. Output (exit code,
stdout, stderr, JSON) is byte-identical to running the official file directly,
apart from the --instrument / --n_band flags this wrapper consumes.

Opt-in (OSC_ED_PAIRING=1 / true / yes / on): this is an ED calibration path,
not a submission path. n_band comes from --n_band, or OSC_ED_N_BAND, or
(when --material_id is one of the 22 package materials) automatic lookup in
inputs/materials_canonical.csv. The official row's solver_delta is not used.
The wrapper writes a separately labelled calibration file (default
ed_calibration.json, never the official submission.json default), uses the
ED A_d and a calibration method name, and prints a stderr warning not to
submit it to the leaderboard.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ed_pairing import (  # noqa: E402
    env_opt_in, estimate_pairing, locked_result, lookup_material_n_band,
)

CALIB_METHOD = "ED 1/N extrapolation (calibration only, not a replacement)"
CALIB_DEFAULT_OUT = "ed_calibration.json"
OFFICIAL_DEFAULT_OUT = "submission.json"


def load_official(path):
    spec = importlib.util.spec_from_file_location("osc_instrument_official", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def official_instrument_path(explicit=None):
    if explicit:
        return explicit
    env = os.environ.get("OSC_INSTRUMENT")
    if env:
        return env
    return os.path.join(HERE, "third_party", "osc_instrument.py")


def _calibration_out(path):
    """Don't write a file that looks like a leaderboard submission."""
    if not path or path == OFFICIAL_DEFAULT_OUT or os.path.basename(path) == OFFICIAL_DEFAULT_OUT:
        return CALIB_DEFAULT_OUT
    root, ext = os.path.splitext(path)
    if root.endswith("_ed_calibration") or os.path.basename(path) == CALIB_DEFAULT_OUT:
        return path
    return f"{root}_ed_calibration{ext or '.json'}"


def _calibration_level(A_d):
    if A_d is None:
        return "ED calibration only, not a submission (no locked A_d)"
    if A_d <= 0:
        return "ED calibration only, not a submission (locked intercept clamped to 0)"
    return "ED calibration only, not a submission (locked finite-cluster intercept)"


def _require_official(path):
    if not os.path.isfile(path):
        print(f"[error] official osc_instrument.py not found at {path}. "
              "Run ./fetch_official.sh (the file is not bundled) and rerun.",
              file=sys.stderr)
        sys.exit(2)


def run_official(path, official_argv):
    """Run the official CLI with estimate_pairing left in place."""
    _require_official(path)
    mod = load_official(path)
    sys.argv = [path] + official_argv
    mod.main()


def run_calibration(official_argv, n_band_arg):
    """Opt-in path: labelled ED calibration JSON, not a submission file."""
    print("[warn] OSC_ED_PAIRING is on: this output is an ED calibration result, "
          "not a leaderboard submission. Do not submit it on the OSC Space.",
          file=sys.stderr)

    ap = argparse.ArgumentParser(description="ED calibration wrapper (not a submission)")
    ap.add_argument("--material_id", default="")
    ap.add_argument("--t", type=float)
    ap.add_argument("--U", type=float)
    ap.add_argument("--nef", type=float)
    ap.add_argument("--A_d", type=float)
    ap.add_argument("--method_name", default=CALIB_METHOD)
    ap.add_argument("--author", default="")
    ap.add_argument("--out", default=OFFICIAL_DEFAULT_OUT)
    a = ap.parse_args(official_argv)

    n_band = n_band_arg
    if n_band is None:
        env_n = os.environ.get("OSC_ED_N_BAND")
        if env_n not in (None, ""):
            n_band = float(env_n)
        elif a.material_id:
            n_band = lookup_material_n_band(a.material_id)
    if n_band is None:
        print("[error] n_band is required when OSC_ED_PAIRING is on. Pass --n_band "
              "or OSC_ED_N_BAND, or a --material_id from inputs/materials_canonical.csv "
              "(automatic lookup). The official solver_delta is not used.",
              file=sys.stderr)
        sys.exit(2)

    nef, A_d_official = a.nef, a.A_d
    formula = None
    if a.material_id and (nef is None or A_d_official is None or formula is None):
        # Fill missing nef / official A_d from the package table (no network).
        row = lookup_material_row(a.material_id)
        if row:
            formula = row.get("formula")
            if nef is None:
                nef = float(row["nef"])
            if A_d_official is None and row.get("old_A_d") not in (None, ""):
                A_d_official = float(row["old_A_d"])
            if n_band_arg is None and os.environ.get("OSC_ED_N_BAND") in (None, ""):
                # already set from lookup_material_n_band above
                pass
        elif nef is None:
            print(f"[error] {a.material_id} not in inputs/materials_canonical.csv "
                  "and --nef was not given.", file=sys.stderr)
            sys.exit(2)
    if nef is None:
        print("[error] provide --material_id OR --nef (and --n_band when opt-in is on)",
              file=sys.stderr)
        sys.exit(2)

    res = locked_result(n_band)
    if res["A_d"] is None:
        print(f"[error] no locked A_d at n_band={n_band}: {res['flag']}", file=sys.stderr)
        sys.exit(2)
    A_d_ed = res["A_d"]
    score = estimate_pairing(nef, n_band=n_band, use_ed=True)
    out = _calibration_out(a.out)
    level = _calibration_level(A_d_ed)
    # If the caller set a custom method_name, keep it, but never silently
    # use the official baseline default.
    method = a.method_name
    if method == "OSC-baseline (published A_d)":
        method = CALIB_METHOD
    ut = 8.0
    if a.t and a.U:
        ut = a.U / a.t
    mid = a.material_id or f"custom(nef={nef},n_band={n_band})"
    sub = dict(
        material_id=mid,
        t=a.t, U=a.U, nef=nef, U_over_t=round(ut, 4),
        A_dwave=A_d_ed,
        official_A_dwave=A_d_official,
        n_band=float(n_band), n_eff=res["n_eff"],
        solver_delta=None,
        provisional_score=score,
        level=level,
        method_name=method,
        author=a.author,
        track="calibration",
        calibration_only=True,
        flag=res["flag"] or "",
        formula=formula,
    )
    print(f"[level] {level}  (score {score}; ED A_d={A_d_ed:.10g}; calibration only, not a submission)")
    json.dump(sub, open(out, "w"), indent=1)
    print(json.dumps(sub, indent=1))
    print(f"\n[OK] wrote {out}. This is an ED calibration file, not a leaderboard submission.")


def lookup_material_row(mid):
    import csv
    path = os.path.join(HERE, "inputs", "materials_canonical.csv")
    if not os.path.isfile(path):
        return None
    for r in csv.DictReader(open(path)):
        if r.get("osc_id") == mid or r.get("material_id") == mid:
            return r
    return None


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--instrument", default=None)
    pre.add_argument("--n_band", type=float, default=None)
    a, rest = pre.parse_known_args(argv)
    path = official_instrument_path(a.instrument)
    on, env_warn = env_opt_in()
    if env_warn:
        print(env_warn, file=sys.stderr)
    if not on:
        # Official path: do not replace estimate_pairing; do not set OSC_ED_N_BAND.
        run_official(path, rest)
        return
    run_calibration(rest, a.n_band)


if __name__ == "__main__":
    main()
