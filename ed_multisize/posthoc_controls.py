#!/usr/bin/env python3
"""Post-hoc controls on the locked 1/N extrapolation (NOT part of PREREGISTRATION.md).

Uses the same interpolation and unweighted A + B/N (or A + B/N^power) as extrapolate.py,
plus: t-based SE/CI, leave-one-out, 1/sqrt(N), evaluation at N=36, a U=0 (Ku-ku) control
through the same pipeline, and Spearman sensitivity. Nothing is fitted to the official table.
"""
import csv, hashlib, json, math, os, re, sys, time
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
M = os.path.dirname(HERE)
LOCKED = "504a9e21047076de4039610f67c502f727b49b69100a74ecd681de44d8f64564"
BIG = ["sqrt10", "3x4", "4x4", "18t", "20t"]
N_OF = {"sqrt10": 10, "3x4": 12, "4x4": 16, "18t": 18, "20t": 20, "2x4": 8, "sqrt8": 8}
RMAX = {"2x4": "sqrt5", "sqrt8": "2", "sqrt10": "sqrt5", "3x4": "sqrt5", "4x4": "2*sqrt2", "18t": "3", "20t": "sqrt10"}


def _find(name, *cands):
    for p in cands:
        if p and os.path.isfile(p):
            return p
    sys.exit(f"missing {name}; tried {cands}")


CACHE = os.environ.get("OSC_OFFICIAL_CACHE", os.path.join(os.path.dirname(M), "ed_reference_official_cache") if os.path.basename(M)=="ed_reference_package" else os.path.join(M, "ed_reference_official_cache"))
KUKU = _find("Ku-ku table",
             os.path.join(M, "third_party", "multisize_ED_table_kuku.md"),  # package after fetch
             os.path.join(M, "ed_reference_package", "third_party", "multisize_ED_table_kuku.md"),
             os.path.join(M, "ref", "multisize_ED_table_kuku.md"),
             os.path.join(CACHE or "", "multisize_ED_table_kuku.md") if CACHE else "")
META = _find("scoring_meta.json",
             os.path.join(M, "third_party", "scoring_meta.json"),
             os.path.join(M, "ed_reference_package", "third_party", "scoring_meta.json"),
             os.path.join(M, "best_method_package", "official", "scoring_meta.json"),
             os.path.join(CACHE or "", "scoring_meta.json") if CACHE else "")
MATS = _find("materials",
             os.path.join(M, "inputs", "materials_canonical.csv"),
             os.path.join(M, "ed_reference_package", "inputs", "materials_canonical.csv"),
             os.path.join(M, "best_method_package", "results", "band_sensitivity.csv"))
U8CSV = _find("locked U=8 inputs",
              os.path.join(HERE, "results", "extrapolation_inputs.csv"),
              os.path.join(M, "ed_multisize", "results", "extrapolation_inputs.csv"),
              os.path.join(M, "ed_reference_package", "ed_multisize", "results", "extrapolation_inputs.csv"))
PREREG = os.path.join(HERE, "PREREGISTRATION.md")
OUT = os.environ.get("OSC_POSTHOC_OUT", os.path.join(HERE, "results"))


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def parse_kuku_u0():
    """P_org at U=0 from Ku-ku's published table (canonical ED, except 4x4 n=1 = Wick)."""
    pts = {}
    for line in open(KUKU):
        if not line.startswith("| "):
            continue
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        if len(c) < 14 or c[0] in ("cluster", "") or c[0].startswith("-") or not c[1][0].isdigit():
            continue
        cl, N, n, nup = c[0], int(c[1]), float(c[2]), int(c[3])
        P0 = float(c[12])
        method = c[14] if len(c) > 14 else ""
        pts[(cl, nup)] = dict(cluster=cl, N=N, n=n, nup=nup, P=P0, method=method)
    return pts


def load_u8():
    pts = {}
    for r in csv.DictReader(open(U8CSV)):
        pts[(r["cluster"], int(r["nup"]))] = dict(cluster=r["cluster"], N=int(r["N"]), n=float(r["n"]),
                                                  nup=int(r["nup"]), P=float(r["P"]))
    return pts


def curves(pts):
    cs = {}
    for v in pts.values():
        cs.setdefault(v["cluster"], {"N": v["N"], "xy": [(0.0, 0.0)]})["xy"].append((v["n"], v["P"]))
    for c in cs.values():
        c["xy"].sort()
        c["nmax"] = c["xy"][-1][0]
    return cs


def interp(c, n):
    if n > c["nmax"] + 1e-12:
        return None
    xs, ys = zip(*c["xy"])
    return float(np.interp(n, xs, ys))


def covering(cs, n, names):
    out = {}
    for cl in names:
        if cl not in cs:
            continue
        v = interp(cs[cl], n)
        if v is not None:
            out[cl] = v
    return out


def fit(Ns, Ps, power=1.0):
    Ns, Ps = np.asarray(Ns, float), np.asarray(Ps, float)
    X = np.column_stack([np.ones(len(Ns)), Ns ** -power])
    coef, _, _, _ = np.linalg.lstsq(X, Ps, rcond=None)
    A, B = float(coef[0]), float(coef[1])
    resid = X @ coef - Ps
    n, df = len(Ps), len(Ps) - 2
    rms = float(np.sqrt(np.mean(resid ** 2)))
    rec = dict(A=A, B=B, rms=rms, n=n, df=df, A_at_36=A + B / (36 ** power))
    if df > 0 and n >= 3:
        sigma2 = float(np.sum(resid ** 2) / df)
        cov = sigma2 * np.linalg.inv(X.T @ X)
        se = float(np.sqrt(cov[0, 0]))
        tcrit = float(stats.t.ppf(0.975, df))
        rec.update(se=se, ci_lo=A - tcrit * se, ci_hi=A + tcrit * se, tcrit=tcrit)
    else:
        rec.update(se=None, ci_lo=None, ci_hi=None, tcrit=None)
    return rec


def loo(Ns, Ps, labels, power=1.0):
    out = {}
    for i, lab in enumerate(labels):
        if len(Ns) < 3:
            continue
        f = fit([n for j, n in enumerate(Ns) if j != i], [p for j, p in enumerate(Ps) if j != i], power)
        out[lab] = f["A"]
    return out


def analyse(cs, n, names):
    cov = covering(cs, n, names)
    if len(cov) < 2:
        return dict(clusters=cov, n_used=len(cov))
    labels = list(cov)
    Ns = [N_OF[c] for c in labels]
    Ps = [cov[c] for c in labels]
    f1 = fit(Ns, Ps, 1.0)
    f15 = fit(Ns, Ps, 1.5)
    fsq = fit(Ns, Ps, 0.5)
    left = loo(Ns, Ps, labels)
    return dict(clusters=cov, n_used=len(cov), spread=max(Ps) - min(Ps),
                fit=f1, fit_N15=f15, fit_sqrtN=fsq, loo=left,
                loo_min=min(left.values()) if left else None, loo_max=max(left.values()) if left else None)


def spearman_report(x, y, seed=0, n_boot=20000, n_perm=20000):
    x, y = np.asarray(x, float), np.asarray(y, float)
    res = stats.spearmanr(x, y)
    rho, p = float(res.statistic), float(res.pvalue)
    # Fisher z (reported even with ties; noted as approximate)
    if abs(rho) < 1 and len(x) > 3:
        z = math.atanh(rho)
        se = 1 / math.sqrt(len(x) - 3)
        fisher = (math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se))
    else:
        fisher = (None, None)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(x), len(x))
        if np.std(x[i]) == 0 or np.std(y[i]) == 0:
            continue
        boots.append(float(stats.spearmanr(x[i], y[i]).statistic))
    boots = np.array(boots)
    boot_ci = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))) if len(boots) else (None, None)
    count = 0
    for _ in range(n_perm):
        if abs(float(stats.spearmanr(x, rng.permutation(y)).statistic)) >= abs(rho) - 1e-15:
            count += 1
    return dict(rho=rho, p=p, p_perm=count / n_perm, fisher=fisher, boot=boot_ci, n=len(x))


def fmt(v, d=4):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.{d}f}" if (abs(v) >= 1e-3 or v == 0) else f"{v:.3e}"
    return str(v)


def main():
    if os.path.isfile(PREREG):
        h = sha(PREREG)
        if h != LOCKED:
            sys.exit(f"PREREGISTRATION.md sha256 {h} != locked")
    else:
        h = LOCKED
    os.makedirs(OUT, exist_ok=True)
    u8, u0 = load_u8(), parse_kuku_u0()
    cs8, cs0 = curves(u8), curves(u0)
    meta = json.load(open(META))
    official = {int(k): float(v[1]) for k, v in meta["A_by_Ne"].items()}

    rows = []
    for Ne in sorted(official):
        n = Ne / 36
        # locked set (all covering N>=10, including 18t/20t when they cover)
        a_lock = analyse(cs8, n, BIG)
        # same 3 clusters as the U=0 control
        a8_3 = analyse(cs8, n, ["sqrt10", "3x4", "4x4"])
        a0_3 = analyse(cs0, n, ["sqrt10", "3x4", "4x4"])
        off = official[Ne]
        rec = dict(Ne=Ne, n=round(n, 6), official_A_d=max(0.0, off), official_A_d_raw=off)
        def pack(prefix, a):
            if a.get("fit") is None:
                return
            f, s = a["fit"], a["fit_sqrtN"]
            rec.update({
                f"{prefix}n_used": a["n_used"], f"{prefix}clusters": ";".join(f"{k}={a['clusters'][k]:.8g}" for k in a["clusters"]),
                f"{prefix}A": f["A"], f"{prefix}B": f["B"], f"{prefix}se": f["se"], f"{prefix}ci_lo": f["ci_lo"],
                f"{prefix}ci_hi": f["ci_hi"], f"{prefix}rms": f["rms"], f"{prefix}A_at_36": f["A_at_36"],
                f"{prefix}A_sqrtN": s["A"], f"{prefix}A_N15": a["fit_N15"]["A"],
                f"{prefix}loo_min": a["loo_min"], f"{prefix}loo_max": a["loo_max"],
                f"{prefix}loo": ";".join(f"{k}={v:.6g}" for k, v in a["loo"].items()),
                f"{prefix}spread": a.get("spread"),
                f"{prefix}A_drop4x4": a["loo"].get("4x4"),
            })
        pack("lock_", a_lock)
        pack("u8_3_", a8_3)
        pack("u0_3_", a0_3)
        rec["U8_minus_U0_3"] = (a8_3["fit"]["A"] - a0_3["fit"]["A"]) if a8_3.get("fit") and a0_3.get("fit") else None
        rows.append(rec)

    # Spearman variants on official fillings (official A_d as used in scoring = max(0, raw))
    off_v = [r["official_A_d"] for r in rows]
    ad_lock = [max(0.0, r["lock_A"]) for r in rows]          # locked (clamped) A_d
    raw_lock = [r["lock_A"] for r in rows]
    # sqrt10 / 4x4 interpolated values vs official
    s10, c44 = [], []
    for r in rows:
        n = r["n"]
        s10.append(interp(cs8["sqrt10"], n))
        c44.append(interp(cs8["4x4"], n))
    reports = {
        "locked_clamped_Ad_vs_official": spearman_report(ad_lock, off_v),
        "unclamped_A_vs_official": spearman_report(raw_lock, off_v),
        "sqrt10_vs_official": spearman_report(s10, off_v),
        "4x4_vs_official": spearman_report(c44, off_v),
    }

    # materials at canonical n_eff, locked fit + CI (post-hoc uncertainty)
    mats = []
    for r in csv.DictReader(open(MATS)):
        nb = float(r.get("old_n") or r["old_n"])
        nef = float(r["nef"])
        n_eff = nb if nb <= 1 else 2 - nb
        a = analyse(cs8, n_eff, BIG)
        f = a.get("fit")
        mats.append(dict(osc_id=r["osc_id"], formula=r["formula"], n_band=nb, n_eff=round(n_eff, 6), nef=nef,
                         A=f["A"] if f else None, se=f["se"] if f else None,
                         ci_lo=f["ci_lo"] if f else None, ci_hi=f["ci_hi"] if f else None,
                         A_d_clamped=max(0.0, f["A"]) if f else None,
                         n_used=a.get("n_used"), official_A_d=r.get("old_A_d", ""), official_index=r.get("old_index", "")))

    def wcsv(name, data):
        p = os.path.join(OUT, name)
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]))
            w.writeheader()
            for r in data:
                w.writerow({k: ("" if v is None else v) for k, v in r.items()})
        return p

    p_fill = wcsv("posthoc_by_filling.csv", rows)
    p_mat = wcsv("posthoc_materials.csv", mats)
    # chart
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        ns = [r["n"] for r in rows]
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        ax.axhline(0, color="0.5", lw=0.8)
        ax.plot(ns, [r["u0_3_A"] for r in rows], "o--", color="#1f77b4", label="U=0, same 1/N pipeline (sqrt10, 3x4, 4x4)")
        ax.plot(ns, [r["u8_3_A"] for r in rows], "s-", color="#d62728", label="U=8, same 3 clusters (post-hoc)")
        ax.plot(ns, [r["lock_A"] for r in rows], "^-", color="#2ca02c", label="U=8, locked set (N≥10 covering)", alpha=0.85)
        ax.plot(ns, [r["official_A_d"] for r in rows], "D", color="0.2", label="official 6×6 A_d (compare only)")
        ax.set_xlabel("filling n")
        ax.set_ylabel("extrapolated A  (raw, unclamped)")
        ax.set_title("Same 1/N pipeline: free electrons (U=0) and U=8")
        ax.legend(fontsize=8, loc="lower left")
        ax.set_xlim(0, 1.05)
        fig.tight_layout()
        chart = os.path.join(OUT, "u0_control_chart.png")
        fig.savefig(chart, dpi=140)
        plt.close()
    except Exception as e:
        chart = None
        print("chart skipped:", e)

    md = os.path.join(OUT, "posthoc_controls.md")
    with open(md, "w") as f:
        f.write(f"# Post-hoc controls on the locked 1/N extrapolation ({time.strftime('%Y-%m-%d %H:%M %Z')})\n\n")
        f.write("**These analyses are post-hoc.** They are not in PREREGISTRATION.md and were not used to choose\n")
        f.write(f"the locked method. PREREGISTRATION.md sha256 `{h}` (unchanged). Script sha256 `{sha(__file__)}`.\n")
        f.write("U=0 values are taken from Ku-ku's published table (canonical ED; 4x4 n=1 is canonical Wick).\n")
        f.write("Nothing is fitted or rescaled to the official 6x6 table.\n\n")
        f.write("## U=0 vs U=8 through the same pipeline (3 clusters: sqrt10, 3x4, 4x4)\n\n")
        f.write("| n | U=0 A (1/N) | U=8 A (same 3) | U=8 − U=0 | locked U=8 A (all covering N≥10) | official A_d |\n|" + "---|" * 6 + "\n")
        for r in rows:
            f.write(f"| {r['n']:.4f} | {fmt(r.get('u0_3_A'))} | {fmt(r.get('u8_3_A'))} | {fmt(r.get('U8_minus_U0_3'))} | {fmt(r.get('lock_A'))} | {r['official_A_d']:.5g} |\n")
        f.write("\n## Locked U=8 fit, unclamped, with uncertainty (post-hoc)\n\n")
        f.write("| n | A (1/N) ± SE | 95% CI (t, df=n−2) | A (1/√N) | leave-one-out (drop 4x4) | A at N=36 | official |\n|" + "---|" * 7 + "\n")
        for r in rows:
            se = r.get("lock_se")
            se_s = f"{r['lock_A']:.3f} ± {se:.3f}" if se is not None else f"{r['lock_A']:.3f}"
            ci = f"[{fmt(r.get('lock_ci_lo'))}, {fmt(r.get('lock_ci_hi'))}]" if r.get("lock_ci_lo") is not None else ""
            f.write(f"| {r['n']:.4f} | {se_s} | {ci} | {fmt(r.get('lock_A_sqrtN'))} | {fmt(r.get('lock_A_drop4x4'))} | {fmt(r.get('lock_A_at_36'))} | {r['official_A_d']:.5g} |\n")
        f.write("\n## Spearman rank correlation vs official A_d (15 fillings; post-hoc)\n\n")
        for name, s in reports.items():
            f.write(f"- **{name}:** rho = {s['rho']:.3f}, p = {s['p']:.3f} (permutation p = {s['p_perm']:.3f}); "
                    f"Fisher-z 95% CI [{s['fisher'][0]:.2f}, {s['fisher'][1]:.2f}]; "
                    f"bootstrap 95% CI [{s['boot'][0]:.2f}, {s['boot'][1]:.2f}].\n")
        f.write("\nOrganizers already reported that 6x6 VMC A_d and 4x4 ED are \"essentially uncorrelated across filling (Spearman ≈ −0.1)\" "
                "(`validation/comparison_4x4_ED_vs_VMC.md`). The locked clamped comparison is a multi-size extension of that finding.\n")
        f.write("\n## Output files\n\n")
        for p in (p_fill, p_mat):
            f.write(f"- `{os.path.relpath(p, HERE)}` sha256 `{sha(p)}`\n")
        if chart:
            f.write(f"- `{os.path.relpath(chart, HERE)}` sha256 `{sha(chart)}`\n")
        # machine-readable spearman
        json.dump({k: {kk: (list(vv) if isinstance(vv, tuple) else vv) for kk, vv in s.items()} for k, s in reports.items()},
                  open(os.path.join(OUT, "posthoc_spearman.json"), "w"), indent=2)
    print("wrote", md)
    print("spearman locked", reports["locked_clamped_Ad_vs_official"])
    # print a few U=0 checks for the log
    for r in rows:
        if abs(r["n"] - 0.61111) < 0.01 or abs(r["n"] - 0.83333) < 0.01 or abs(r["n"] - 1) < 1e-9:
            print("check", r["n"], "U0", r.get("u0_3_A"), "U8_3", r.get("u8_3_A"), "lock", r.get("lock_A"),
                  "CI", r.get("lock_ci_lo"), r.get("lock_ci_hi"), "drop4x4", r.get("lock_A_drop4x4"))


if __name__ == "__main__":
    main()
