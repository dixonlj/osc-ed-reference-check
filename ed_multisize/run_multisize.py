#!/usr/bin/env python3
"""Resumable orchestrator for the pre-registered 18t/20t points (PREREGISTRATION.md, sha256 in LOG.md).
Runs points one at a time (cheapest first) via run_point_ms.py; a point is redone only if its json is
missing, aborted, or a retryable memory-wait skip (partial sector checkpoints are reused).
--report-only rebuilds results/multisize_runs.{csv,md} from runs/*.json (and the ed_solver 8-16 site runs)."""
import argparse, csv, glob, json, os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import geometry  # noqa
from lattice import Cluster  # noqa
from run_point import estimate_gb  # noqa
RUNS = os.path.join(HERE, "runs"); RES = os.path.join(HERE, "results")
EDRUNS = os.path.join(HERE, "..", "ed_solver", "runs")
POINTS = [("18t", n) for n in range(1, 10)] + [("20t", n) for n in range(1, 11)]


def status(cl, nup):
    p = os.path.join(RUNS, f"{cl}_nup{nup}_U8_k.json")
    return json.load(open(p)) if os.path.exists(p) else None


def done(r):
    return r is not None and (r.get("status") == "ok" or str(r.get("status", "")).startswith("skipped: estimated"))


def report():
    rows = []
    for cl, nup in POINTS:
        c = Cluster(cl); r = status(cl, nup) or {}
        est = estimate_gb(c, nup, "k")
        rows.append(dict(cluster=cl, N=c.N, nup=nup, n=round(2 * nup / c.N, 4), sector_dim=(__import__("math").comb(c.N, nup) ** 2) // c.N,
                         est_gb=round(est, 2), plan="run" if est <= 5 else "skip (memory)", status=r.get("status", "pending"),
                         E0=r.get("E0", ""), deg=r.get("deg", ""), P_org=r.get("P_org", ""), runtime_s=r.get("runtime_s", ""),
                         peak_rss_gb=r.get("peak_rss_gb", "")))
    # existing 8-16 site results from the Ku-ku-checked solver (read only)
    for p in sorted(glob.glob(os.path.join(EDRUNS, "*_U8_*.json"))):
        r = json.load(open(p))
        if r.get("status") != "ok" or r.get("cluster") == "2x2":
            continue
        rows.append(dict(cluster=r["cluster"], N=r["N"], nup=r["nup"], n=round(r["n"], 4), sector_dim="", est_gb=r.get("est_gb", ""),
                         plan=f"ed_solver ({r['method']})", status="ok", E0=r["E0"], deg=r["deg"], P_org=r["P_org"],
                         runtime_s=r.get("runtime_s", ""), peak_rss_gb=r.get("peak_rss_gb", "")))
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "multisize_runs.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(os.path.join(RES, "multisize_runs.md"), "w") as f:
        f.write(f"# Multi-size ED runs (U/t=8), rebuilt {time.strftime('%Y-%m-%d %H:%M %Z')}\n\n")
        f.write("Pre-registered plan: PREREGISTRATION.md. Raw P_org only; extrapolation/A_d mapping is a separate step.\n\n")
        f.write("| " + " | ".join(rows[0]) + " |\n|" + "---|" * len(rows[0]) + "\n")
        for r in rows:
            f.write("| " + " | ".join(str(v) for v in r.values()) + " |\n")
    n_ok = sum(1 for r in rows[:len(POINTS)] if r["status"] == "ok")
    n_skip = sum(1 for r in rows[:len(POINTS)] if str(r["status"]).startswith("skipped"))
    return n_ok, n_skip


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--report-only", action="store_true"); a = ap.parse_args()
    if a.report_only:
        print(report()); return
    order = sorted(POINTS, key=lambda p: estimate_gb(Cluster(p[0]), p[1], "k"))
    for pas in (1, 2):
        for cl, nup in order:
            if done(status(cl, nup)):
                continue
            tag = f"{cl}_nup{nup}_U8_k"
            print(time.strftime("%Y-%m-%d %H:%M:%S %Z"), f"START {tag} (pass {pas})", flush=True)
            with open(os.path.join(RUNS, tag + ".log"), "a") as lg:
                rc = subprocess.call([sys.executable, os.path.join(HERE, "run_point_ms.py"), "--cluster", cl, "--nup", str(nup)],
                                     stdout=lg, stderr=subprocess.STDOUT)
            print(time.strftime("%Y-%m-%d %H:%M:%S %Z"), f"END {tag} rc {rc} status {(status(cl, nup) or {}).get('status')}", flush=True)
            report()
    print(time.strftime("%Y-%m-%d %H:%M:%S %Z"), "ALL DONE", report(), flush=True)


if __name__ == "__main__":
    main()
