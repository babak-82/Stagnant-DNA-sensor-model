#!/usr/bin/env python3
"""Extra checks added for the ABEC revision. Every number they print is written to results/claims.jsonl and claims.md.

Usage:  python extra_checks.py                      # all stages (about 2 minutes on one CPU core)
        python extra_checks.py --stages tolerance,kon

Run it in the same folder as reproduce_all.py, after reproduce_all.py has run once (it reuses its helper functions).

Stages
  tolerance   baseline detection limit against solver tolerances and grid (Sec. 2.3)
  laplace     Laplace-domain benchmark in the semi-infinite and the finite domain, two linearisations (Sec. 2.3)
  beta        slowdown ratio at equal Da* for beta = 0.01, 0.1, 1, 10 at K = 1.78, f = 0.5 (Sec. 3.2)
  kon         detection limit and stagnant / well-mixed ratio for k_on,0 scaled by 0.01 to 1 (Sec. 3.5 and 3.7)
  groups      Da, K, beta at the baseline detection limits (Sec. 3.4)
  optimum     optimum probe density by direct minimisation, and the gain over the baseline probe density (Sec. 3.6)
  thermal     critical temperature difference of the column against its height (Sec. 3.7)
  ratio       stagnant / well-mixed ratio range with four decimals (needs results/ratio_range.json)
"""
import argparse, json, os, time
import numpy as np
import mpmath as mp
from scipy.optimize import minimize_scalar
import stagnant_cell_model as SM
from stagnant_cell_model import solve
import reproduce_all as R
from lod_model import BASE, lod, theta_stagnant, theta_wellmixed, groups

SREL = 0.02


def stage_tolerance():
    def th(c, p, n, rt, at):
        Da, K, b = groups(c, p); tt = p["D"] * 2700.0 / p["h"] ** 2
        return solve(Da, K, b, t_end=tt, t_eval=np.array([tt]), n=n, rtol=rt, atol=at)["theta"][-1]
    out = []
    for n, rt, at in [(110, 1e-6, 1e-10), (110, 1e-9, 1e-12), (300, 1e-9, 1e-12)]:
        out.append(f"n={n}, rtol={rt:g}: {lod(lambda c: th(c, BASE, n, rt, at), BASE):.5f}")
    R.claim("baseline LOD(45 min) against grid and solver tolerance (nM)", "; ".join(out))


def stage_laplace():
    orig = SM.surface_state
    def make(mode):
        def f(u0, th, x0, Da, K):
            lam = K if mode == "K" else 1 + K
            us = (u0 / x0 + Da * lam * th) / (1.0 / x0 + Da); return us, (u0 - us) / x0
        return f
    cases = [(20, 0.5, 5.0), (5, 2.0, 5.0), (200, 1.0, 10.0)]
    try:
        for dom, ts in (("semi-infinite", np.array([1e-4, 1e-3, 3e-3, 1e-2, 2e-2])), ("finite", np.array([0.01, 0.1, 0.5, 1.0, 2.0]))):
            for mode in ("K", "1+K"):
                SM.surface_state = make(mode); worst = 0.0
                for Da, K, b in cases:
                    lam = K if mode == "K" else 1 + K
                    r = solve(Da, K, b, t_end=ts[-1], t_eval=ts, n=300, stretch=10.0, rtol=1e-8, atol=1e-11)
                    if dom == "finite":
                        F = lambda s: Da / (s * (b * s + Da * b * mp.sqrt(s) / mp.tanh(mp.sqrt(s)) + Da * lam))
                    else:
                        F = lambda s: Da / (s * (b * s + Da * b * mp.sqrt(s) + Da * lam))
                    ref = np.array([float(mp.invertlaplace(F, t, method="talbot")) for t in ts])
                    worst = max(worst, float(np.max(np.abs(r["theta"] - ref) / ref)))
                R.claim(f"Laplace benchmark, {dom} domain, linearisation lambda = {mode}: max relative error", f"{worst:.1e}")
    finally:
        SM.surface_state = orig


def stage_beta():
    K = 1.78; Dastar = [1, 3.16, 10, 31.6, 100, 316, 1000]; tab = {}
    for b in (0.01, 0.1, 1.0, 10.0):
        tab[b] = [R.slowdown(d / b, K, b, fracs=(0.5,), n=150)[0] for d in Dastar]
        R.claim(f"slowdown ratio f=0.5, K=1.78, beta={b}, at Da* = {Dastar}", [round(float(x), 2) for x in tab[b]])
    for b in (0.01, 1.0, 10.0):
        R.claim(f"ratio beta={b} / beta=0.1 at equal Da*, K=1.78", [round(float(x / y), 2) for x, y in zip(tab[b], tab[0.1])])
    x = np.log10(Dastar); y1 = np.log10(tab[0.1]); y10 = np.log10(tab[10.0])
    sh = []
    for target in (0.3, 0.5, 0.8):
        if y10.max() >= target and y1.max() >= target:
            sh.append(10 ** (np.interp(target, y10, x) - np.interp(target, y1, x)))
    R.claim("shift of the beta = 10 curve along Da*, relative to beta = 0.1 (range over ratios 2 to 6.3)", f"{min(sh):.1f} to {max(sh):.1f}")


def stage_kon():
    for fac in (0.01, 0.03, 0.1, 0.3, 1.0):
        q = dict(BASE); q["kon0"] = BASE["kon0"] * fac
        s = lod(lambda c: theta_stagnant(c, q, 2700.0), q, lo=-4, hi=6); w = lod(lambda c: theta_wellmixed(c, q, 2700.0), q, lo=-4, hi=6)
        Da, _, _ = groups(s, q)
        R.claim(f"k_on,0 x {fac}: Da, stagnant LOD, well-mixed LOD (nM), ratio at 45 min", f"{Da:.1f}, {s:.3f}, {w:.3f}, {s / w:.2f}")


def stage_groups():
    for name, fn in (("stagnant", theta_stagnant), ("well-mixed", theta_wellmixed)):
        c = lod(lambda x: fn(x, BASE, 2700.0), BASE)
        Da, K, b = groups(c, BASE)
        R.claim(f"baseline {name} LOD (nM), then Da, K, beta, Da* at that concentration", f"{c:.4f}, {Da:.1f}, {K:.3f}, {b:.3f}, {Da * b:.1f}")


def stage_optimum():
    def f(model, rho, lsp):
        q = dict(BASE); q["SP"] = 10 ** lsp; q["sigma"] = float(np.sqrt((rho * SREL * 2e10 / q["SP"]) ** 2 + SREL ** 2))
        fn = (lambda c: theta_stagnant(c, q, 2700.0)) if model == "stag" else (lambda c: theta_wellmixed(c, q, 2700.0))
        return lod(fn, q, lo=-4, hi=6, it=22)
    res = {}
    for rho in (0.25, 1.0, 4.0):
        for model in ("stag", "wm"):
            m = minimize_scalar(lambda x: f(model, rho, x), bounds=(9.3, 11.0), method="bounded", options={"xatol": 2e-3})
            base = f(model, rho, np.log10(2e10)); res[(rho, model)] = 10 ** m.x * 1e6
            R.claim(f"direct minimisation, sigma_a/sigma_r = {rho}, {model}: optimum S_P (m^-2), min LOD (nM), LOD at baseline S_P, gain",
                    f"{10 ** m.x * 1e6:.2e}, {m.fun:.3f}, {base:.3f}, {base / m.fun:.2f}")
    R.claim("direct minimisation: optimum S_P ratio well-mixed / stagnant for sigma_a/sigma_r = 0.25, 1, 4",
            ", ".join(f"{res[(r, 'wm')] / res[(r, 'stag')]:.2f}" for r in (0.25, 1.0, 4.0)))


def stage_thermal():
    g, nu, kappa, beta_T, Ra_c = 9.81, 0.95e-6, 1.43e-7, 2.3e-4, 1708.0
    for h in (9.9e-3, 7.5e-3, 5.0e-3, 3.3e-3):
        R.claim(f"critical vertical temperature difference (K) of a water column of height {h * 1e3:.1f} mm (Ra_c = 1708)", f"{Ra_c * nu * kappa / (g * beta_T * h ** 3):.2f}")


def stage_ratio():
    path = os.path.join(R.OUT, "ratio_range.json")
    if not os.path.exists(path):
        print("ratio_range.json not found: run reproduce_all.py --stages ratio_range first"); return
    r = [x[4] for x in json.load(open(path))]
    R.claim("stagnant / well-mixed LOD ratio range at 45 min, four decimals", f"{min(r):.4f}, {max(r):.4f}")


STAGES = {"tolerance": stage_tolerance, "laplace": stage_laplace, "beta": stage_beta, "kon": stage_kon, "groups": stage_groups,
          "optimum": stage_optimum, "thermal": stage_thermal, "ratio": stage_ratio}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--stages", default="all"); a = ap.parse_args()
    names = list(STAGES) if a.stages == "all" else a.stages.split(",")
    for n in names:
        s = time.time(); print(f"=== stage {n}", flush=True); STAGES[n](); print(f"=== {n} done in {time.time() - s:.0f} s", flush=True)
    rows = list({r["claim"]: r for r in (json.loads(l) for l in open(os.path.join(R.OUT, "claims.jsonl")))}.values())
    with open(os.path.join(R.OUT, "claims.md"), "w") as f:
        f.write("| claim | value |\n|---|---|\n")
        for r in rows: f.write(f"| {r['claim']} | {r['value']} |\n")
    print("claims.md updated")
