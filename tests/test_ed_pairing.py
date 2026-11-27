#!/usr/bin/env python3
"""Tests for the optional, off-by-default ED estimate_pairing().

(a) Default off: scores match the *fetched official* estimate_pairing on the 22
    materials; when the official CLI is available, wrap_osc_instrument.py is
    also byte-identical to it (exit code, stdout, stderr, JSON).
(b) Opt-in on: scores match ed_multisize/results/extrapolation_materials.csv.
(c) Fillings outside the computed / documented range behave as documented,
    including a real expected value for the unextrapolated fallback.
(d) Opt-in wrapper output is labelled as an ED calibration, not a submission.

Test artefacts are written to a temp directory unless OSC_ED_PAIRING_OUT is set.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, PKG)

from ed_pairing import (  # noqa: E402
    estimate_pairing, locked_A_d, locked_result, n_eff, env_opt_in,
    lookup_material_n_band,
)

MATS = os.path.join(PKG, "inputs", "materials_canonical.csv")
EXTRAP = os.path.join(PKG, "ed_multisize", "results", "extrapolation_materials.csv")
FILL = os.path.join(PKG, "ed_multisize", "results", "extrapolation_by_filling.csv")
PY = os.environ.get("PY", sys.executable)


def _outdir():
    env = os.environ.get("OSC_ED_PAIRING_OUT")
    if env:
        os.makedirs(env, exist_ok=True)
        return env
    return tempfile.mkdtemp(prefix="ed_pairing_test_")


OUTDIR = None  # set in main()


def official_instrument_path():
    candidates = [
        os.environ.get("OSC_INSTRUMENT") or "",
        os.path.join(PKG, "third_party", "osc_instrument.py"),
        os.path.join(os.environ.get("OSC_OFFICIAL_CACHE", ""), "osc_instrument.py")
        if os.environ.get("OSC_OFFICIAL_CACHE") else "",
        os.path.join(os.path.dirname(PKG), "ed_reference_official_cache", "osc_instrument.py"),
    ]
    for p in candidates:
        if p and os.path.isfile(p):
            return p
    return None


def load_official_estimate_pairing(path):
    spec = importlib.util.spec_from_file_location("osc_instrument_official_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.estimate_pairing


def wcsv(path, rows, fields):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})
    return path


def test_a_default_off_matches_official():
    inst = official_instrument_path()
    if inst is None:
        print("SKIP (a) official-file comparison: osc_instrument.py not fetched "
              "(run ./fetch_official.sh or set OSC_INSTRUMENT / OSC_OFFICIAL_CACHE)")
        return "skipped"
    official_fn = load_official_estimate_pairing(inst)
    materials = list(csv.DictReader(open(MATS)))
    assert len(materials) == 22
    os.environ.pop("OSC_ED_PAIRING", None)
    official_rows, wrapper_rows = [], []
    for r in materials:
        nef, A_d = float(r["nef"]), float(r["old_A_d"])
        off = official_fn(nef, A_d)
        wrap = estimate_pairing(nef, A_d)
        official_rows.append(dict(osc_id=r["osc_id"], formula=r["formula"], nef=r["nef"],
                                  A_d=r["old_A_d"], index=f"{off:.3f}"))
        wrapper_rows.append(dict(osc_id=r["osc_id"], formula=r["formula"], nef=r["nef"],
                                 A_d=r["old_A_d"], index=f"{wrap:.3f}"))
        assert wrap == off, (r["osc_id"], wrap, off)
        assert estimate_pairing(nef, A_d, n_band=float(r["old_n"]), use_ed=False) == off
    fields = ["osc_id", "formula", "nef", "A_d", "index"]
    p_off = wcsv(os.path.join(OUTDIR, "official_pipeline.csv"), official_rows, fields)
    p_wrap = wcsv(os.path.join(OUTDIR, "wrapper_off.csv"), wrapper_rows, fields)
    assert open(p_off, "rb").read() == open(p_wrap, "rb").read()
    assert sum(float(r["index"]) > 0 for r in official_rows) == 22

    # Full official CLI vs wrap_osc_instrument.py, default-off, all 22 materials.
    wrap_py = os.path.join(PKG, "wrap_osc_instrument.py")
    for r in materials:
        with tempfile.TemporaryDirectory() as td_off, tempfile.TemporaryDirectory() as td_wrap:
            env = os.environ.copy()
            env.pop("OSC_ED_PAIRING", None)
            env["PYTHONPATH"] = PKG + os.pathsep + env.get("PYTHONPATH", "")
            args = ["--nef", r["nef"], "--A_d", r["old_A_d"], "--out", "out.json"]
            off = subprocess.run([PY, inst, *args], cwd=td_off, env=env,
                                 capture_output=True, text=True)
            wrap = subprocess.run(
                [PY, wrap_py, "--instrument", inst, *args],
                cwd=td_wrap, env=env, capture_output=True, text=True)
            assert off.returncode == 0 and wrap.returncode == 0, (r["osc_id"], off.stderr, wrap.stderr)
            assert off.returncode == wrap.returncode
            assert off.stdout == wrap.stdout, (r["osc_id"], off.stdout, wrap.stdout)
            assert off.stderr == wrap.stderr, (r["osc_id"], off.stderr, wrap.stderr)
            assert open(os.path.join(td_off, "out.json"), "rb").read() == open(
                os.path.join(td_wrap, "out.json"), "rb").read()
    return "ok"


def test_b_opt_in_matches_extrapolation_materials():
    official = {r["osc_id"]: r for r in csv.DictReader(open(MATS))}
    expected = list(csv.DictReader(open(EXTRAP)))
    assert len(expected) == 22
    os.environ["OSC_ED_PAIRING"] = "1"
    try:
        rows = []
        n_zero = 0
        for e in expected:
            mid = e["osc_id"]
            o = official[mid]
            nef, n_band = float(o["nef"]), float(o["old_n"])
            exp_index = float(e["index_ED"]) if e["index_ED"] not in ("", None) else 0.0
            exp_Ad = float(e["A_d_ED"])
            got = estimate_pairing(nef, A_d=float(o["old_A_d"]), n_band=n_band)
            got_Ad = locked_A_d(n_band)
            assert f"{got_Ad:.10g}" == f"{exp_Ad:.10g}", (mid, got_Ad, exp_Ad)
            assert got == exp_index, (mid, got, exp_index)
            assert estimate_pairing(nef, n_band=n_band, use_ed=True) == exp_index
            rows.append(dict(osc_id=mid, formula=e["formula"], n_band=e["n_band"],
                             A_d_ED=e["A_d_ED"], index_ED=e["index_ED"],
                             index_wrapper=f"{got:.3f}"))
            if got == 0:
                n_zero += 1
        assert n_zero == 20
        by_id = {r["osc_id"]: r for r in expected}
        assert float(by_id["OSC-02297"]["index_ED"]) == 21.286
        assert float(by_id["OSC-01253"]["index_ED"]) == 12.033
        # automatic n_band from package table
        assert lookup_material_n_band("OSC-02297") == 1.3012
    finally:
        os.environ.pop("OSC_ED_PAIRING", None)
    wcsv(os.path.join(OUTDIR, "wrapper_on.csv"), rows,
         ["osc_id", "formula", "n_band", "A_d_ED", "index_ED", "index_wrapper"])
    return n_zero


def test_c_range_behaviour():
    em = {r["osc_id"]: r for r in csv.DictReader(open(EXTRAP))}
    n_band = 1.3012
    assert abs(n_eff(n_band) - 0.6988) < 1e-9
    got, exp = locked_A_d(n_band), float(em["OSC-02297"]["A_d_ED"])
    assert f"{got:.10g}" == f"{exp:.10g}", (got, exp)

    for r in csv.DictReader(open(FILL)):
        n = int(r["Ne"]) / 36
        got = locked_A_d(n)
        exp = float(r["A_d"])
        assert f"{got:.10g}" == f"{exp:.10g}", (r["Ne"], n, got, exp)

    assert locked_A_d(0.0) == 0.0
    half = locked_result(1.0)
    assert half["A_d"] == 0.0
    assert "clamped" in half["flag"]

    for bad in (-0.1, 2.5, 3):
        try:
            n_eff(bad)
        except ValueError as e:
            assert "outside the documented range" in str(e)
        else:
            raise AssertionError(f"n_band={bad} should have raised")
        try:
            estimate_pairing(1.0, n_band=bad, use_ed=True)
        except ValueError:
            pass
        else:
            raise AssertionError(f"estimate_pairing should have raised for n_band={bad}")

    os.environ.pop("OSC_ED_N_BAND", None)
    os.environ.pop("OSC_ED_MATERIAL_ID", None)
    try:
        estimate_pairing(1.0, A_d=0.05, use_ed=True)
    except ValueError as e:
        assert "n_band is required" in str(e)
    else:
        raise AssertionError("missing n_band should have raised")

    # Unrecognised env value is off, with a warning.
    os.environ["OSC_ED_PAIRING"] = "enabled"
    try:
        on, warn = env_opt_in()
        assert on is False and warn and "unrecognised" in warn
        # library call should use official A_d and still succeed
        assert estimate_pairing(1.0, 0.05) == 15.0
    finally:
        os.environ.pop("OSC_ED_PAIRING", None)

    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "inputs.csv")
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["cluster", "N", "nup", "n", "P", "deg", "method", "file", "sha"])
            w.writeheader()
            # only 18t-like coverage to n=0.5
            w.writerow(dict(cluster="18t", N=18, nup=5, n=0.5, P=0.01, deg=1, method="k", file="x", sha="0"))
        res = locked_result(0.9, inputs_path=p)
        assert res["A_d"] is None and res["flag"] == "no cluster covers n"
        try:
            estimate_pairing(1.0, n_band=0.9, use_ed=True, inputs_path=p)
        except ValueError as e:
            assert "no locked A_d" in str(e)
        else:
            raise AssertionError("expected no locked A_d")
        # one covering cluster: unextrapolated fallback is P interpolated
        # between (0, 0) and (0.5, 0.01) at n=0.4, i.e. 0.008.
        res2 = locked_result(0.4, inputs_path=p)
        assert "unextrapolated" in res2["flag"]
        assert res2["n_used"] == 1
        assert abs(res2["A"] - 0.008) < 1e-15, res2["A"]
        assert res2["A_d"] == res2["A"]


def test_d_opt_in_wrapper_labels():
    inst = official_instrument_path()
    if inst is None:
        print("SKIP (d) wrapper label test: osc_instrument.py not fetched")
        return "skipped"
    wrap_py = os.path.join(PKG, "wrap_osc_instrument.py")
    mats = {r["osc_id"]: r for r in csv.DictReader(open(MATS))}
    em = {r["osc_id"]: r for r in csv.DictReader(open(EXTRAP))}

    def run_on(mid, extra=None, out_name="submission.json"):
        r = mats[mid]
        env = os.environ.copy()
        env["OSC_ED_PAIRING"] = "1"
        env["PYTHONPATH"] = PKG + os.pathsep + env.get("PYTHONPATH", "")
        env.pop("OSC_ED_N_BAND", None)
        args = [PY, wrap_py, "--instrument", inst,
                "--material_id", mid, "--nef", r["nef"], "--A_d", r["old_A_d"],
                "--out", out_name]
        if extra:
            args.extend(extra)
        return subprocess.run(args, cwd=OUTDIR, env=env, capture_output=True, text=True)

    # Auto n_band from the package table; official default out name must be redirected.
    p = run_on("OSC-02297")
    assert p.returncode == 0, p.stderr
    assert "not a leaderboard submission" in p.stderr
    assert "Do not submit" in p.stderr
    cal_path = os.path.join(OUTDIR, "ed_calibration.json")
    assert os.path.isfile(cal_path), (p.stdout, os.listdir(OUTDIR))
    assert not os.path.isfile(os.path.join(OUTDIR, "submission.json"))
    data = json.loads(open(cal_path).read())
    assert data["provisional_score"] == 21.286
    assert abs(float(data["A_dwave"]) - float(em["OSC-02297"]["A_d_ED"])) < 1e-12 or \
        f"{float(data['A_dwave']):.10g}" == f"{float(em['OSC-02297']['A_d_ED']):.10g}"
    assert float(data["A_dwave"]) != float(mats["OSC-02297"]["old_A_d"])
    assert data["method_name"] == "ED 1/N extrapolation (calibration only, not a replacement)"
    assert "OSC-baseline" not in data["method_name"]
    assert "STRONG" not in data["level"]
    assert "calibration only" in data["level"]
    assert data.get("calibration_only") is True
    assert data.get("track") == "calibration"
    assert "not a leaderboard submission" in p.stdout or "not a leaderboard submission" in p.stderr
    assert "submit this score" not in p.stdout.lower()

    # A clamped material is also labelled as calibration, not "off-dome" from the official A_d.
    p2 = run_on("OSC-01188", extra=["--n_band", mats["OSC-01188"]["old_n"]],
                out_name="mine.json")
    assert p2.returncode == 0, p2.stderr
    # explicit --out that isn't submission.json gets an _ed_calibration suffix
    out2 = os.path.join(OUTDIR, "mine_ed_calibration.json")
    assert os.path.isfile(out2), os.listdir(OUTDIR)
    d2 = json.loads(open(out2).read())
    assert d2["provisional_score"] == 0
    assert float(d2["A_dwave"]) == 0.0
    assert "clamped" in d2["level"].lower() or "clamped" in (d2.get("flag") or "").lower()
    assert "STRONG" not in d2["level"]
    assert "OSC-baseline" not in d2["method_name"]
    return "ok"


def main():
    global OUTDIR
    OUTDIR = _outdir()
    a = test_a_default_off_matches_official()
    print(f"PASS (a) default off matches fetched official estimate_pairing and CLI ({a})")
    n_zero = test_b_opt_in_matches_extrapolation_materials()
    print(f"PASS (b) opt-in matches extrapolation_materials.csv (20 of 22 A_d=0; {n_zero} zeros)")
    test_c_range_behaviour()
    print("PASS (c) fillings outside the computed / documented range behave as documented")
    d = test_d_opt_in_wrapper_labels()
    print(f"PASS (d) opt-in wrapper is labelled as ED calibration, not a submission ({d})")
    print("OK all ed_pairing tests")
    print("test artefacts:", OUTDIR)


if __name__ == "__main__":
    main()
