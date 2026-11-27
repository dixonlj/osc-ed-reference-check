#!/usr/bin/env python3
"""Correctness check of the ED solver against Ku-ku's published multi-size ED table.

Runs each point as its own low-priority, single-threaded process (nice 19), skips points
that already have an ok JSON (resumable), and rebuilds results/ed_check.{csv,md} after each.
Reference values are parsed directly from method/ref/multisize_ED_table_kuku.md; nothing is
typed in by hand, and nothing is rescaled or fitted to them.

Usage:  python run_check.py [--report-only] [--only 4x4]
"""
import argparse, csv, json, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)  # package root
REF = os.path.join(ROOT, "third_party", "multisize_ED_table_kuku.md")  # fetched by fetch_official.sh (sha256-checked)
REF2 = os.path.join(ROOT, "third_party", "comparison_4x4_ED_vs_VMC.md")
RUNS = os.path.join(HERE, "runs")
RES = os.path.join(HERE, "results")
PY = sys.executable

POINTS = (
    [("2x4", n, "plain") for n in range(1, 5)]
    + [("sqrt8", n, "plain") for n in range(1, 5)]
    + [("sqrt10", n, "plain") for n in range(1, 6)]
    + [("3x4", n, "plain") for n in range(1, 7)]
    + [("4x4", n, "plain") for n in range(1, 5)]
    + [("4x4", n, "k") for n in range(3, 9)]
)

E_TOL = 2e-8          # Ku-ku E0 printed to 8 decimals (truncated)
P_REL = 1e-5          # Ku-ku P values printed to ~6 significant digits (truncated)
P_ABS = 1e-8


def parse_ref():
    ref = {}
    for line in open(REF):
        if not line.startswith("| ") or line.startswith("| cluster") or line.startswith("|---"):
            continue
        cells = [x.strip() for x in line.strip().strip("|").split("|")]
        if len(cells) < 15 or cells[0] not in ("2x4", "sqrt8", "sqrt10", "3x4", "4x4"):
            continue
        def num(x):
            try:
                return float(x)
            except ValueError:
                return None
        ref[(cells[0], int(cells[3]))] = dict(
            n=num(cells[2]), dim=int(cells[4]), E0=num(cells[5]), deg=int(cells[6]),
            deg_U0=int(cells[7]), R_max=num(cells[8]), P_d=num(cells[9]), P_org=num(cells[10]),
            P_org_U0=num(cells[12]), raw_P_org=cells[10], raw_E0=cells[5])
    return ref


def run_point(cluster, nup, method):
    tag = f"{cluster}_nup{nup}_U8_{method}"
    out = os.path.join(RUNS, tag + ".json")
    if os.path.isfile(out):
        try:
            if json.load(open(out)).get("status") == "ok":
                return
        except Exception:
            pass
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               PYTHONHASHSEED="0")
    logf = open(os.path.join(RUNS, tag + ".log"), "w")
    print(time.strftime("%Y-%m-%d %H:%M:%S %Z"), "START", tag, flush=True)
    rc = subprocess.call(["nice", "-n", "19", PY, os.path.join(HERE, "run_point.py"),
                          "--cluster", cluster, "--nup", str(nup), "--method", method],
                         cwd=HERE, env=env, stdout=logf, stderr=subprocess.STDOUT)
    print(time.strftime("%Y-%m-%d %H:%M:%S %Z"), "END", tag, "rc", rc, flush=True)
    if not os.path.isfile(out):
        json.dump(dict(tag=tag, status=f"failed rc={rc}"), open(out, "w"), indent=1)


def fmt(x, nd=10):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)


def build_report():
    ref = parse_ref()
    rows = []
    for cluster, nup, method in POINTS:
        tag = f"{cluster}_nup{nup}_U8_{method}"
        p = os.path.join(RUNS, tag + ".json")
        r = json.load(open(p)) if os.path.isfile(p) else dict(status="not run yet")
        k = ref.get((cluster, nup), {})
        import math as _m
        true_dim = _m.comb(r["N"], nup) ** 2 if r.get("status") == "ok" else None
        row = dict(cluster=cluster, N_up=nup, method=method, status=r.get("status"),
                   n=r.get("n", k.get("n")), dim=true_dim, ref_dim=k.get("dim"),
                   E0=r.get("E0"), ref_E0=k.get("E0"), deg=r.get("deg"), ref_deg=k.get("deg"),
                   R_max=r.get("R_max"), ref_R_max=k.get("R_max"),
                   P_org=r.get("P_org"), ref_P_org=k.get("P_org"),
                   residual_max=r.get("residual_max"), solver=r.get("solver"),
                   runtime_s=r.get("runtime_s"), peak_rss_gb=r.get("peak_rss_gb"))
        if r.get("status") == "ok" and k:
            dE = r["E0"] - k["E0"]
            dP = r["P_org"] - k["P_org"]
            row.update(dE0=dE, dP_org=dP)
            okE = abs(dE) <= E_TOL
            okP = abs(dP) <= P_REL * abs(k["P_org"]) + P_ABS
            okD = r["deg"] == k["deg"]
            okR = abs(r["R_max"] - k["R_max"]) < 1e-3
            row["match"] = "yes" if (okE and okP and okD and okR) else (
                "NO (" + ",".join(x for x, ok in (("E0", okE), ("P_org", okP), ("deg", okD), ("R_max", okR)) if not ok) + ")")
        rows.append(row)
    os.makedirs(RES, exist_ok=True)
    keys = ["cluster", "n", "N_up", "method", "status", "dim", "ref_dim", "E0", "ref_E0", "dE0", "deg", "ref_deg",
            "R_max", "ref_R_max", "P_org", "ref_P_org", "dP_org", "match", "residual_max", "solver",
            "runtime_s", "peak_rss_gb"]
    with open(os.path.join(RES, "ed_check.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, keys, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    done = [r for r in rows if r.get("match")]
    nmatch = sum(1 for r in done if r["match"] == "yes")
    L = ["# ED solver correctness check vs Ku-ku multi-size ED (U/t = 8)", "",
         f"Generated {time.strftime('%Y-%m-%d %H:%M:%S %Z')}. Reference: `method/ref/multisize_ED_table_kuku.md` "
         "(values truncated there). This is a correctness check only; nothing was tuned or rescaled.", "",
         "Operator and normalization: Ku-ku pair operator, P_org = 2 P_d at R_max, degenerate ground states averaged "
         "with equal weight (|E - E0| < 1e-7). Seed 0 for every start vector.", "",
         f"Match rule (set before comparing): |ΔE0| ≤ {E_TOL:g}; |ΔP_org| ≤ {P_REL:g}·|ref| + {P_ABS:g}; same degeneracy and R_max.", "",
         f"Points compared so far: {len(done)}; matching: {nmatch}.", "",
         "| cluster | n | N_up | method | dim | E0 ours | E0 Ku-ku | ΔE0 | deg ours/ref | P_org ours | P_org Ku-ku | ΔP_org | match | residual | runtime s | peak GB | status |",
         "|---|---:|---:|---|---:|---:|---:|---:|---|---:|---:|---:|---|---:|---:|---:|---|"]
    for r in rows:
        L.append(f"| {r['cluster']} | {fmt(r['n'],4)} | {r['N_up']} | {r['method']} | {fmt(r['dim'])} | "
                 f"{fmt(r['E0'],11)} | {fmt(r['ref_E0'],11)} | {fmt(r.get('dE0'),3)} | {fmt(r['deg'])}/{fmt(r['ref_deg'])} | "
                 f"{fmt(r['P_org'],8)} | {fmt(r['ref_P_org'],8)} | {fmt(r.get('dP_org'),3)} | {r.get('match','')} | "
                 f"{fmt(r['residual_max'],2)} | {fmt(r['runtime_s'])} | {fmt(r['peak_rss_gb'])} | {r['status']} |")
    L += ["", "Methods: `plain` = real-space basis, dense diagonalization for dim ≤ 6000, otherwise scipy eigsh "
          "(matrix-free LinearOperator) plus Lanczos deflation; `k` = crystal-momentum sectors (C4v stars on 4x4) "
          "with a low-memory two-pass Lanczos plus deflation. Peak GB is the process max RSS. Every process ran "
          "with nice 19, one thread, oom_score_adj 1000, and a MemAvailable watchdog."]
    open(os.path.join(RES, "ed_check.md"), "w").write("\n".join(L) + "\n")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report-only", action="store_true")
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    os.makedirs(RUNS, exist_ok=True)
    if not a.report_only:
        for cluster, nup, method in POINTS:
            if a.only and cluster != a.only:
                continue
            run_point(cluster, nup, method)
            build_report()
    build_report()
    print("wrote", os.path.join(RES, "ed_check.md"))


if __name__ == "__main__":
    main()
