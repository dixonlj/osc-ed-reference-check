#!/usr/bin/env python3
"""One multi-size ED point -> runs/<tag>.json. Same physics/safety as ed_solver/run_point.py, plus:
  * general K-star reduction from the automatically detected point group (geometry.stars),
  * per-sector checkpoint runs/<tag>.partial.json, so a killed/restarted point resumes,
  * memory wait uses 1.15 x estimate + 2.5 GB (stricter than ed_solver; budget rule unchanged: est <= 5 GB).
Safety: caller sets nice 19 and 1 thread; oom_score_adj 1000; watchdog aborts below 1.2 GB free.
"""
import argparse, json, os, resource, sys, threading, time, math, platform

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import geometry  # noqa: E402  (registers 18t/20t, adds ed_solver to path)
from run_point import mem_available_gb, estimate_gb  # noqa: E402  (ed_solver helpers, read-only import)

RUNS = os.path.join(HERE, "runs")
WATCH_MIN_GB, RESERVE_GB, MARGIN = 1.2, 2.5, 1.15


def watchdog(tag, out_path):
    while True:
        m = mem_available_gb()
        if m < WATCH_MIN_GB:
            json.dump(dict(tag=tag, status=f"aborted: MemAvailable {m:.2f} GB < {WATCH_MIN_GB} GB (watchdog)",
                           time=time.strftime("%Y-%m-%d %H:%M:%S %Z")), open(out_path, "w"), indent=1)
            os._exit(3)
        time.sleep(1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cluster", required=True)
    ap.add_argument("--nup", type=int, required=True)
    ap.add_argument("--U", type=float, default=8.0)
    ap.add_argument("--max-gb", type=float, default=5.0)
    ap.add_argument("--wait-min", type=float, default=720.0)
    ap.add_argument("--no-stars", action="store_true", help="all K individually (used by tests)")
    ap.add_argument("--runs", default=RUNS)
    a = ap.parse_args()

    from lattice import Cluster
    c = Cluster(a.cluster)
    tag = f"{a.cluster}_nup{a.nup}_U{a.U:g}_k"
    os.makedirs(a.runs, exist_ok=True)
    out_path = os.path.join(a.runs, tag + ".json")
    part_path = os.path.join(a.runs, tag + ".partial.json")
    try:
        open("/proc/self/oom_score_adj", "w").write("1000")
    except OSError:
        pass
    est = estimate_gb(c, a.nup, "k")
    need = MARGIN * est + RESERVE_GB
    base = dict(tag=tag, cluster=a.cluster, N=c.N, nup=a.nup, ndn=a.nup, n=2 * a.nup / c.N, U=a.U, method="k",
                est_gb=round(est, 2), need_free_gb=round(need, 2), seed=0, python=platform.python_version(),
                sector_dim=int(math.comb(c.N, a.nup) ** 2 // c.N))
    if est > a.max_gb:
        base.update(status=f"skipped: estimated {est:.1f} GB > budget {a.max_gb} GB (pre-registered rule)")
        json.dump(base, open(out_path, "w"), indent=1); print(base["status"]); return 0
    t_wait = time.time()
    while mem_available_gb() < need:
        if time.time() - t_wait > a.wait_min * 60:
            base.update(status=f"skipped: MemAvailable {mem_available_gb():.1f} GB < need {need:.1f} GB after {a.wait_min} min (retryable)")
            json.dump(base, open(out_path, "w"), indent=1); print(base["status"]); return 0
        time.sleep(60)
    threading.Thread(target=watchdog, args=(tag, out_path), daemon=True).start()

    import numpy as np
    from hubbard import configs
    from hubbard_k import KSector, translation_tables

    def log(msg):
        print(time.strftime("%H:%M:%S"), msg, flush=True)
    part = json.load(open(part_path)) if os.path.exists(part_path) else {"energy": {}, "ground": {}, "elapsed_s": 0.0}
    save = lambda: json.dump(part, open(part_path, "w"), indent=1)
    t0 = time.time() - part.get("elapsed_s", 0.0)
    R = c.r_max()
    rec = dict(base, R_max=R, shell0=c.shell(0, R), mem_available_at_start_gb=round(mem_available_gb(), 2),
               point_group=geometry.point_group(c))
    tabs = translation_tables(c, configs(c.N, a.nup))
    if a.no_stars:
        ks = c.k_points(); sector_list = [(i, [i]) for i in range(len(ks))]
        rec["k_reduction"] = "none (all K)"
    else:
        ks, sector_list = geometry.stars(c)
        rec["k_reduction"] = f"stars under detected point group {rec['point_group']}"
    total_mv = sum(v["matvecs"] for v in part["energy"].values())
    for i, members in sector_list:
        key = str(i)
        if key in part["energy"]:
            log(f"resume: sector {key} already done E={part['energy'][key]['E']:.10f}"); continue
        s = KSector(c, a.nup, a.nup, a.U, ks[i], tables=tabs)
        ts = time.time()
        E, nmv = s.energy_lanczos(log=log)
        total_mv += nmv
        part["energy"][key] = dict(K_over_pi=[round(float(x) / math.pi, 6) for x in ks[i]], star_size=len(members),
                                   E=E, sector_dim=s.dim, matvecs=nmv, seconds=round(time.time() - ts, 1))
        part["elapsed_s"] = time.time() - t0; save()
        log(f"sector K/pi={part['energy'][key]['K_over_pi']} star={len(members)} dim={s.dim} E={E:.10f}")
        del s
    Emin = min(v["E"] for v in part["energy"].values())
    for i, members in sector_list:
        key = str(i)
        if part["energy"][key]["E"] - Emin > 1e-5 or key in part["ground"]:
            continue
        s = KSector(c, a.nup, a.nup, a.U, ks[i], tables=tabs)
        E0s, Xs, ress, nmv = s.lowest_lanczos(log=log)
        total_mv += nmv
        g = dict(E_refined=E0s, matvecs=nmv)
        if not (abs(E0s - Emin) >= 1e-7 and E0s > Emin):
            g.update(deg_in_sector=len(Xs), P_d_each=[s.pd(X, R)[1] for X in Xs], residuals=ress)
        part["ground"][key] = g
        part["elapsed_s"] = time.time() - t0; save()
        log(f"refined sector {key}: {g}")
        del s, Xs
    E0 = min(g["E_refined"] for g in part["ground"].values())
    P_sum, deg, res_all, sectors = 0.0, 0, [], []
    for i, members in sector_list:
        sec = dict(part["energy"][str(i)])
        g = part["ground"].get(str(i))
        if g:
            sec["E_refined"] = g["E_refined"]
            inman = "P_d_each" in g and abs(g["E_refined"] - E0) < 1e-7
            sec["in_ground_manifold"] = inman
            if inman:
                sec.update(deg_in_sector=g["deg_in_sector"], P_d_each=g["P_d_each"], residuals=g["residuals"])
                P_sum += len(members) * sum(g["P_d_each"]); deg += len(members) * g["deg_in_sector"]
                res_all += g["residuals"]
        sectors.append(sec)
    rec.update(E0=E0, deg=deg, solver="K-sector Lanczos + deflation (ed_solver/hubbard_k.py)", residual_max=max(res_all),
               P_d=P_sum / deg, P_org=2 * P_sum / deg, sectors=sectors, matvecs=total_mv,
               dim_total=int(math.comb(c.N, a.nup) ** 2))
    rec["runtime_s"] = round(time.time() - t0, 1)
    rec["peak_rss_gb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 3)
    rec["status"] = "ok"
    rec["time"] = time.strftime("%Y-%m-%d %H:%M:%S %Z")
    json.dump(rec, open(out_path, "w"), indent=1)
    log(f"DONE {tag} E0={E0:.10f} deg={deg} P_org={rec['P_org']:.8g} t={rec['runtime_s']}s peak={rec['peak_rss_gb']}GB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
