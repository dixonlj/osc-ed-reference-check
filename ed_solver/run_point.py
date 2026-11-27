#!/usr/bin/env python3
"""Run one ED point and write runs/<tag>.json (with runtime and peak memory).

Safety for the shared box: single thread (set by caller), oom_score_adj=1000 so the kernel
kills this process before the DFT queue, a start-time memory check, and a watchdog that
aborts if MemAvailable falls below WATCH_MIN_GB.
"""
import argparse, json, os, resource, sys, threading, time, math, platform

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RUNS = os.path.join(HERE, "runs")
WATCH_MIN_GB = 1.2
RESERVE_GB = 2.5


def mem_available_gb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) / 1e6
    return 0.0


def watchdog(tag, out_path):
    while True:
        m = mem_available_gb()
        if m < WATCH_MIN_GB:
            json.dump(dict(tag=tag, status=f"aborted: MemAvailable {m:.2f} GB < {WATCH_MIN_GB} GB (protecting DFT queue)",
                           time=time.strftime("%Y-%m-%d %H:%M:%S %Z")), open(out_path, "w"), indent=1)
            os._exit(3)
        time.sleep(1.0)


def estimate_gb(cluster, nup, method):
    from math import comb
    N = cluster.N
    du = comb(N, nup)
    if method == "plain":
        dim = du * du
        return dim * 8 * 26 / 1e9 + 0.2  # eigsh ncv<=20 + deflation + temps
    dim = du * du / N
    return dim * 16 * 9 / 1e9 + 0.3      # Lanczos vectors + deflation + matvec temps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cluster", required=True)
    ap.add_argument("--nup", type=int, required=True)
    ap.add_argument("--U", type=float, default=8.0)
    ap.add_argument("--method", choices=["plain", "k"], default="plain")
    ap.add_argument("--max-gb", type=float, default=5.0)
    ap.add_argument("--wait-min", type=float, default=30.0, help="minutes to wait for free memory before skipping")
    a = ap.parse_args()

    from lattice import Cluster
    c = Cluster(a.cluster)
    tag = f"{a.cluster}_nup{a.nup}_U{a.U:g}_{a.method}"
    os.makedirs(RUNS, exist_ok=True)
    out_path = os.path.join(RUNS, tag + ".json")
    try:
        open("/proc/self/oom_score_adj", "w").write("1000")
    except OSError:
        pass
    est = estimate_gb(c, a.nup, a.method)
    base = dict(tag=tag, cluster=a.cluster, N=c.N, nup=a.nup, ndn=a.nup, n=2 * a.nup / c.N, U=a.U,
                method=a.method, est_gb=round(est, 2), seed=0, python=platform.python_version())
    if est > a.max_gb:
        base.update(status=f"skipped: estimated {est:.1f} GB > budget {a.max_gb} GB")
        json.dump(base, open(out_path, "w"), indent=1); print(base["status"]); return 0
    t_wait = time.time()
    while mem_available_gb() < est + RESERVE_GB:
        if time.time() - t_wait > a.wait_min * 60:
            base.update(status=f"skipped: MemAvailable {mem_available_gb():.1f} GB < need {est + RESERVE_GB:.1f} GB after {a.wait_min} min")
            json.dump(base, open(out_path, "w"), indent=1); print(base["status"]); return 0
        time.sleep(60)
    threading.Thread(target=watchdog, args=(tag, out_path), daemon=True).start()

    import numpy as np
    def log(msg):
        print(time.strftime("%H:%M:%S"), msg, flush=True)
    t0 = time.time()
    R = c.r_max()
    rec = dict(base, R_max=R, shell0=c.shell(0, R), mem_available_at_start_gb=round(mem_available_gb(), 2))
    if a.method == "plain":
        from hubbard import Sector, pd_correlation
        s = Sector(c, a.nup, a.nup, a.U)
        E0, w, vecs, how = s.ground_states()
        res = [s.residual(X, E0) for X in vecs]
        Ps = [pd_correlation(s, X, R)[1] for X in vecs]
        rec.update(dim=s.dim, E0=E0, deg=len(vecs), solver=how, residual_max=max(res),
                   P_d=float(np.mean(Ps)), P_org=2 * float(np.mean(Ps)), P_d_each=Ps,
                   low_energies=[float(x) for x in np.sort(w)[:8]])
    else:
        from hubbard import configs
        from hubbard_k import KSector, translation_tables, c4v_stars
        tabs = translation_tables(c, configs(c.N, a.nup))
        if a.cluster in ("4x4", "2x2"):
            ks, stars = c4v_stars(c)
            sector_list = [(i, members) for i, members in stars]
            rec["k_reduction"] = "C4v stars (energy and P_d(R_max) equal within a star; checked against plain solver at low filling)"
        else:
            ks = c.k_points()
            sector_list = [(i, [i]) for i in range(len(ks))]
            rec["k_reduction"] = "none (all K)"
        sectors = []
        total_mv = 0
        for i, members in sector_list:
            s = KSector(c, a.nup, a.nup, a.U, ks[i], tables=tabs)
            ts = time.time()
            E, nmv = s.energy_lanczos(log=log)
            total_mv += nmv
            sectors.append(dict(K_over_pi=[round(float(x) / math.pi, 6) for x in ks[i]], star_size=len(members),
                                E=E, sector_dim=s.dim, matvecs=nmv, seconds=round(time.time() - ts, 1)))
            log(f"sector K/pi={sectors[-1]['K_over_pi']} star={len(members)} dim={s.dim} E={E:.10f}")
            del s
        Emin = min(x["E"] for x in sectors)
        P_sum, deg, res_all = 0.0, 0, []
        for sec, (i, members) in zip(sectors, sector_list):
            if sec["E"] - Emin > 1e-5:
                continue
            s = KSector(c, a.nup, a.nup, a.U, ks[i], tables=tabs)
            E0s, Xs, ress, nmv = s.lowest_lanczos(log=log)
            total_mv += nmv
            sec["E_refined"] = E0s
            if abs(E0s - Emin) >= 1e-7 and E0s > Emin:
                sec["in_ground_manifold"] = False
                continue
            Ps = [s.pd(X, R)[1] for X in Xs]
            sec.update(in_ground_manifold=True, deg_in_sector=len(Xs), P_d_each=Ps, residuals=ress)
            P_sum += len(members) * sum(Ps)
            deg += len(members) * len(Xs)
            res_all += ress
            log(f"ground sector K/pi={sec['K_over_pi']} deg_in_sector={len(Xs)} P_d={Ps}")
            del s
        # recompute E0 from refined ground sectors
        E0 = min(x.get("E_refined", x["E"]) for x in sectors)
        rec.update(E0=E0, deg=deg, solver="K-sector Lanczos + deflation", residual_max=max(res_all),
                   P_d=P_sum / deg, P_org=2 * P_sum / deg, sectors=sectors, matvecs=total_mv,
                   dim=int(round(sum(x["sector_dim"] * x["star_size"] for x in sectors))))
    rec["runtime_s"] = round(time.time() - t0, 1)
    rec["peak_rss_gb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 3)
    rec["status"] = "ok"
    rec["time"] = time.strftime("%Y-%m-%d %H:%M:%S %Z")
    json.dump(rec, open(out_path, "w"), indent=1)
    log(f"DONE {tag} E0={rec['E0']:.10f} deg={rec['deg']} P_org={rec['P_org']:.8g} t={rec['runtime_s']}s peak={rec['peak_rss_gb']}GB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
