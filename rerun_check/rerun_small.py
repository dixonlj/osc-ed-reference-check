#!/usr/bin/env python3
"""Independent-seed re-run of the small stored points; compares E0, degeneracy and P_org with the stored
production results (seed 0). The seed comes from OSC_ED_SEED (read by ed_solver/hubbard.py at import).
Nothing in ed_solver/runs or ed_multisize/runs is written; K-sector re-runs go to rerun_check/runs_seed<S>/."""
import csv, glob, json, os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "ed_multisize"))
import geometry  # noqa: registers 18t/20t, puts ed_solver on the path
import hubbard
from lattice import Cluster
from hubbard import Sector, pd_correlation
import numpy as np
SEED = hubbard.SEED
OUT = os.path.join(HERE, f"runs_seed{SEED}"); os.makedirs(OUT, exist_ok=True)
E_TOL, P_ABS, P_REL = 1e-9, 1e-8, 1e-6
rows = []


def add(stored, E0, deg, P, how, secs):
    dE, dP = E0 - stored["E0"], P - stored["P_org"]
    ok = abs(dE) <= E_TOL and deg == stored["deg"] and abs(dP) <= P_ABS + P_REL * abs(stored["P_org"])
    rows.append(dict(cluster=stored["cluster"], nup=stored["nup"], n=round(stored["n"], 4), how=how, seed=SEED,
                     E0_stored=stored["E0"], E0_rerun=E0, dE=dE, deg_stored=stored["deg"], deg_rerun=deg,
                     P_org_stored=stored["P_org"], P_org_rerun=P, dP=dP, match="yes" if ok else "NO", seconds=round(secs, 1)))
    print(rows[-1]["cluster"], rows[-1]["n"], how, rows[-1]["match"], f"dE={dE:.1e} dP={dP:.1e}", flush=True)


for p in sorted(glob.glob(os.path.join(ROOT, "ed_solver", "runs", "*_plain.json"))):
    r = json.load(open(p))
    if r.get("status") != "ok" or r["N"] > 12 and not (r["cluster"] == "4x4" and r["nup"] <= 3):
        continue
    t0 = time.time(); c = Cluster(r["cluster"]); R = c.r_max()
    s = Sector(c, r["nup"], r["nup"], r["U"]); E0, w, vecs, how = s.ground_states()
    P = 2 * float(np.mean([pd_correlation(s, X, R)[1] for X in vecs]))
    add(r, E0, len(vecs), P, "plain", time.time() - t0)

k_points = [("4x4", 3), ("4x4", 4)] + [(cl, n) for cl in ("18t", "20t") for n in (1, 2, 3, 4)]
for cl, nup in k_points:
    d = "ed_solver" if cl == "4x4" else "ed_multisize"
    r = json.load(open(os.path.join(ROOT, d, "runs", f"{cl}_nup{nup}_U8_k.json")))
    t0 = time.time()
    subprocess.run([sys.executable, os.path.join(ROOT, "ed_multisize", "run_point_ms.py"), "--cluster", cl, "--nup", str(nup),
                    "--runs", OUT, "--wait-min", "10"], check=True, capture_output=True)
    q = json.load(open(os.path.join(OUT, f"{cl}_nup{nup}_U8_k.json")))
    add(r, q["E0"], q["deg"], q["P_org"], "k (stars)", time.time() - t0)

with open(os.path.join(HERE, "rerun_small.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
n_ok = sum(r["match"] == "yes" for r in rows)
with open(os.path.join(HERE, "rerun_small.md"), "w") as f:
    f.write(f"# Independent-seed re-run of small points ({time.strftime('%Y-%m-%d %H:%M %Z')})\n\n")
    f.write(f"Seed {SEED} (production used seed 0). Rule: |dE0| <= {E_TOL}, same degeneracy, |dP_org| <= {P_ABS} + {P_REL}|P|.\n\n")
    f.write(f"**{n_ok}/{len(rows)} points match.** Max |dE0| = {max(abs(r['dE']) for r in rows):.2e}, max |dP_org| = {max(abs(r['dP']) for r in rows):.2e}.\n\n")
    f.write("| cluster | n | method | E0 stored | E0 re-run | dE0 | deg | P_org stored | P_org re-run | dP | match | s |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        f.write(f"| {r['cluster']} | {r['n']} | {r['how']} | {r['E0_stored']:.10f} | {r['E0_rerun']:.10f} | {r['dE']:.1e} | {r['deg_stored']}/{r['deg_rerun']} | {r['P_org_stored']:.10g} | {r['P_org_rerun']:.10g} | {r['dP']:.1e} | {r['match']} | {r['seconds']} |\n")
print(f"{n_ok}/{len(rows)} match")
sys.exit(0 if n_ok == len(rows) else 1)
