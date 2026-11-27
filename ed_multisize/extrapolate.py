#!/usr/bin/env python3
"""Pre-registered finite-size extrapolation of the Ku-ku d-wave pair correlation (U/t=8).
Implements PREREGISTRATION.md (sha256 must equal LOCKED, checked at start; the file is never edited).

(b) n_eff = n_band if n_band <= 1 else 2 - n_band. Per cluster, P_N(n) = linear interpolation between the
    two computed fillings (status ok) that bracket n, with the exact anchor P_N(0)=0; a cluster whose largest
    computed filling is below n does not contribute (no extrapolation in n).
(c) observable = P_org at R_max (degenerate-GS average) straight from the run json; no rescaling.
(d) headline: unweighted least squares P = A + B/N over covering clusters with N >= 10, only if >= 3 cover n;
    else fallback = largest covering cluster, flagged 'unextrapolated'. A_d^ED = max(0, A), flagged if clamped.
    Reported, never substituted: 8-site values, fit RMS residual, max-min spread, alt fit A + B/N^1.5.
(e) official 6x6 table (scoring_meta.json A_by_Ne) is compared only (never fitted or rescaled to).
Implementation choices not fixed by the pre-registration are listed in CHOICES and written to the MD."""
import csv, glob, hashlib, json, math, os, sys, time
import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
M = os.path.dirname(HERE)
LOCKED = "504a9e21047076de4039610f67c502f727b49b69100a74ecd681de44d8f64564"
PREREG = os.path.join(HERE, "PREREGISTRATION.md")
RUN_DIRS = [os.path.join(M, "ed_solver", "runs"), os.path.join(HERE, "runs")]
META = os.path.join(M, "third_party", "scoring_meta.json")  # fetched by fetch_official.sh (sha256-checked)
MATS = os.path.join(M, "inputs", "materials_canonical.csv")  # canonical columns from best_method_package/results/band_sensitivity.csv
BIG = {"sqrt10", "3x4", "4x4", "18t", "20t"}
SMALL = {"2x4", "sqrt8"}
CHOICES = [
    "Inputs: every ok U=8 run json in method/ed_solver/runs and method/ed_multisize/runs; non-ok (aborted/skipped) points are not used.",
    "When a filling has both a 'k' and a 'plain' run (4x4 n=0.375, 0.5), the 'k' value is used; the two agree to < 1e-11.",
    "Fillings evaluated: the official 6x6 fillings n = Ne/36 for every Ne in scoring_meta A_by_Ne (comparison grid), and each material's canonical n_band.",
    "Material inputs: canonical (official band) columns old_n (= n_band, 4 dp as stored), nef, old_A_d, old_index from best_method_package/results/band_sensitivity.csv; the alternative-band columns are not used.",
    "Secondary discretized filling (pre-reg b.4): Ne = round(36 n_band), made even, clamped 2..70, Ne>36 -> 72-Ne, n = Ne/36.",
    "Fit residual reported as RMS of (data - fit) over the clusters in the fit; with 3 clusters the fit has 1 degree of freedom.",
    "Official values compared as used in scoring: max(0, A_by_Ne). Spearman rank correlation over the official fillings where both exist.",
    "Data status: 18t n=0.6667 (nup=6) aborted twice by the 1.2 GB memory watchdog (2026-11-25 15:21, 2026-11-26 02:42 EST) and was left out of the first locked run (2026-11-26 02:53 EST; 18t then covered n<=0.556). A gated retry completed 2026-11-26 12:48:57 EST (runs/18t_nup6_U8_k.json, status ok) and is INCLUDED here, so 18t covers n<=0.667. 18t n>=0.778 and 20t n>=0.6 remain skipped by the pre-registered 5 GB budget. The withheld (without this point) results are kept for the committed with-and-without comparison.",
    "Pairing Index both ways = round(nef * A_d * 300, 3), comparison only; materials with no official index (blank) keep an ED value but are flagged.",
]


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def load_points():
    pts = {}
    for d in RUN_DIRS:
        for p in glob.glob(os.path.join(d, "*_U8_*.json")):
            r = json.load(open(p))
            if r.get("status") != "ok" or r.get("U") != 8.0 or r.get("cluster") == "2x2":
                continue
            key = (r["cluster"], r["nup"])
            if key in pts and pts[key]["method"] == "k":
                continue
            pts[key] = dict(cluster=r["cluster"], N=r["N"], nup=r["nup"], n=2 * r["nup"] / r["N"], P=r["P_org"],
                            deg=r["deg"], method=r["method"], file=os.path.relpath(p, M), sha=sha(p))
    return pts


def curves(pts):
    cs = {}
    for (cl, nup), v in pts.items():
        cs.setdefault(cl, {"N": v["N"], "xy": [(0.0, 0.0)]})["xy"].append((v["n"], v["P"]))
    for c in cs.values():
        c["xy"].sort()
        c["nmax"] = c["xy"][-1][0]
    return cs


def interp(c, n):
    if n > c["nmax"] + 1e-12:
        return None
    xs, ys = zip(*c["xy"])
    return float(np.interp(n, xs, ys))


def fit(Ns, Ps, power=1.0):
    X = np.column_stack([np.ones(len(Ns)), np.array(Ns, float) ** -power])
    coef, *_ = np.linalg.lstsq(X, np.array(Ps), rcond=None)
    rms = float(np.sqrt(np.mean((X @ coef - np.array(Ps)) ** 2)))
    return float(coef[0]), float(coef[1]), rms


def estimate(cs, n_band):
    n = n_band if n_band <= 1 else 2 - n_band
    big = {cl: interp(c, n) for cl, c in cs.items() if cl in BIG}
    big = {k: v for k, v in big.items() if v is not None}
    small = {cl: interp(c, n) for cl, c in cs.items() if cl in SMALL}
    out = dict(n_eff=n, clusters=";".join(f"{k}={big[k]:.6g}" for k in sorted(big, key=lambda k: cs[k]["N"])),
               P8=";".join(f"{k}={v:.6g}" for k, v in sorted(small.items()) if v is not None), n_used=len(big))
    if not big:
        return dict(out, A=None, A_d=None, flag="no cluster covers n")
    Ns = [cs[k]["N"] for k in big]; Ps = list(big.values())
    out["spread"] = max(Ps) - min(Ps)
    if len(big) >= 3:
        A, B, rms = fit(Ns, Ps); A15, B15, rms15 = fit(Ns, Ps, 1.5)
        out.update(A=A, B=B, fit_rms=rms, A_alt_N15=A15, fit_rms_alt=rms15, flag="")
    else:
        largest = max(big, key=lambda k: cs[k]["N"])
        out.update(A=big[largest], B=None, fit_rms=None, A_alt_N15=None, fit_rms_alt=None,
                   flag=f"unextrapolated (only {len(big)} cluster(s) N>=10 cover n; value of {largest})")
    out["A_d"] = max(0.0, out["A"])
    if out["A"] < 0:
        out["flag"] = (out["flag"] + "; " if out["flag"] else "") + "clamped to 0 (negative)"
    return out


def fmt(v, d=6):
    if v is None or v == "":
        return ""
    return f"{v:.{d}g}" if isinstance(v, float) else str(v)


def main():
    h = sha(PREREG)
    if h != LOCKED:
        sys.exit(f"PREREGISTRATION.md sha256 {h} != locked {LOCKED}; refusing to run")
    pts = load_points(); cs = curves(pts)
    meta = json.load(open(META))
    res = os.environ.get("OSC_EXTRAP_OUT", os.path.join(HERE, "results")); os.makedirs(res, exist_ok=True)
    # 1 official-filling grid
    grid = []
    for ne in sorted(int(k) for k in meta["A_by_Ne"]):
        off_raw = meta["A_by_Ne"][str(ne)][1]
        e = estimate(cs, ne / 36)
        grid.append(dict(Ne=ne, n=round(ne / 36, 6), official_A_d_raw=off_raw, official_A_d=max(0.0, off_raw), **e))
    both = [(g["A_d"], g["official_A_d"]) for g in grid if g["A_d"] is not None]
    rho = spearmanr([a for a, _ in both], [b for _, b in both]).statistic if len(both) > 2 else float("nan")
    # 2 materials
    mats = []
    for r in csv.DictReader(open(MATS)):
        nb, nef = float(r["old_n"]), float(r["nef"])
        e = estimate(cs, nb)
        Ne = int(round(nb * 36)); Ne -= Ne % 2; Ne = max(2, min(70, Ne)); Ne_eff = 72 - Ne if Ne > 36 else Ne
        e2 = estimate(cs, Ne_eff / 36)
        off_idx = r["old_index"]
        flag = e["flag"] + ("" if off_idx else ("; " if e["flag"] else "") + "no official index (excluded officially)")
        mats.append(dict(osc_id=r["osc_id"], formula=r["formula"], n_band=nb, n_eff=round(e["n_eff"], 6), nef=nef,
                         A_d_ED=e["A_d"], A_raw=e["A"], fit_rms=e.get("fit_rms"), spread=e.get("spread"), n_clusters=e["n_used"],
                         index_ED=None if e["A_d"] is None else round(nef * e["A_d"] * 300, 3),
                         A_d_ED_discretized=e2["A_d"], Ne_eff=Ne_eff,
                         official_A_d=r["old_A_d"], official_index=off_idx, flag=flag))
    ok_m = [m for m in mats if m["index_ED"] is not None and m["official_index"]]
    rho_m = spearmanr([m["index_ED"] for m in ok_m], [float(m["official_index"]) for m in ok_m]).statistic if len(ok_m) > 2 else float("nan")
    # write CSVs
    def wcsv(name, rows):
        p = os.path.join(res, name)
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
            for r in rows:
                w.writerow({k: fmt(v, 10) if isinstance(v, float) else ("" if v is None else v) for k, v in r.items()})
        return p
    p_in = wcsv("extrapolation_inputs.csv", sorted(pts.values(), key=lambda v: (v["N"], v["cluster"], v["nup"])))
    p_grid = wcsv("extrapolation_by_filling.csv", grid)
    p_mat = wcsv("extrapolation_materials.csv", mats)
    # MD
    md = os.path.join(res, "extrapolation.md")
    with open(md, "w") as f:
        f.write(f"# Pre-registered 1/N extrapolation of A_d (Hubbard U/t=8), run {time.strftime('%Y-%m-%d %H:%M %Z')}\n\n")
        f.write(f"PREREGISTRATION.md sha256 `{h}` (matches locked value). Script `extrapolate.py` sha256 `{sha(__file__)}`.\n")
        f.write("Comparison with the official 6x6 VMC table is report-only; nothing was fitted or rescaled to it or to any published value.\n\n")
        f.write("## Implementation choices (not fixed by the pre-registration)\n\n" + "".join(f"- {c}\n" for c in CHOICES))
        f.write("\n## Cluster coverage (largest computed filling)\n\n| cluster | N | computed n | n_max |\n|---|---|---|---|\n")
        for cl, c in sorted(cs.items(), key=lambda kv: kv[1]["N"]):
            f.write(f"| {cl} | {c['N']} | {', '.join(f'{x:.3f}' for x, _ in c['xy'][1:])} | {c['nmax']:.3f} |\n")
        f.write("\n## A_d^ED at the official 6x6 fillings\n\n| Ne | n | A_d^ED | raw A (fit) | B | fit RMS | spread | alt A (1/N^1.5) | clusters used (P_N(n)) | 8-site (report only) | official A_d | flag |\n|" + "---|" * 12 + "\n")
        for g in grid:
            f.write(f"| {g['Ne']} | {g['n']:.4f} | {fmt(g['A_d'])} | {fmt(g['A'])} | {fmt(g.get('B'))} | {fmt(g.get('fit_rms'), 3)} | {fmt(g.get('spread'), 3)} | {fmt(g.get('A_alt_N15'))} | {g['clusters']} | {g['P8']} | {g['official_A_d']:.5g} | {g['flag']} |\n")
        f.write(f"\nSpearman rank correlation A_d^ED vs official A_d over {len(both)} fillings: {rho:.3f}\n")
        f.write("\n## Materials (canonical band; comparison only, not a submission)\n\n| osc_id | formula | n_band | n_eff | A_d^ED | index (ED) | official A_d | official index | A_d^ED at Ne_eff/36 | flag |\n|---|---|---|---|---|---|---|---|---|---|\n")
        for m in mats:
            f.write(f"| {m['osc_id']} | {m['formula']} | {m['n_band']} | {m['n_eff']:.4f} | {fmt(m['A_d_ED'])} | {fmt(m['index_ED'])} | {m['official_A_d']} | {m['official_index']} | {fmt(m['A_d_ED_discretized'])} | {m['flag']} |\n")
        f.write(f"\nSpearman rank correlation of Pairing Index (ED vs official) over {len(ok_m)} materials: {rho_m:.3f}\n")
        f.write("\n## Output files\n\n")
        for p in (p_in, p_grid, p_mat):
            f.write(f"- `{os.path.relpath(p, HERE)}` sha256 `{sha(p)}`\n")
    print("wrote", md); print("spearman filling", rho, "materials", rho_m)


if __name__ == "__main__":
    main()
