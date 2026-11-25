#!/usr/bin/env python3
"""Reproduce every number, table and figure of the paper
'A reaction-diffusion model shows that diffusion limits stagnant electrochemical DNA sensors at nanomolar concentrations'.

Usage:  python reproduce_all.py                 # everything, from scratch (about 15 min on one CPU core)
        python reproduce_all.py --fast          # coarse grids, pipeline test (about 3 min)
        python reproduce_all.py --stages validate,lod,fig1

Stages (paper location):
  validate        Sec. 2.3  five numerical checks
  convergence     Sec. 2.3  grid refinement (slowdown ratio and baseline detection limit)
  maps            Fig. 2, Sec. 3.1  slowdown ratio over Da and K
  scaling         Fig. 3, Sec. 3.2  collapse onto Da*, large-Da slope
  scaling2        Eq. (4), Sec. 3.2  computed prefactors against the Cottrell estimate
  kpeak           Sec. 3.1  why the f = 0.9 slowdown is not monotonic in K (needs maps)
  closure         Fig. 4e-f, Sec. 3.3  36-combination test of the growing layer
  twocomp         Table 2, Fig. 4a-d, Sec. 3.3  fixed and growing layers (uses closure results)
  lod             Fig. 5, Table 3, Sec. 3.4  calibration curves and detection limits
  tornado         Fig. 6, Sec. 3.5  one-at-a-time sensitivity
  ratio_range     Sec. 3.5  stagnant / well-mixed ratio over the noise and probe-density ranges
  convection      Sec. 4  Rayleigh-number estimate that the column stays stagnant
  sobol           Fig. 7, Sec. 3.5  Sobol indices (resumes from results/sobol_evals.npy if interrupted)
  variance_split  Sec. 3.5  analytic variance split of ln LOD from the Cottrell estimate
  noise           Fig. 8, Sec. 3.6  relative, absolute and mixed noise
  optimum         Fig. 9, Sec. 3.6  optimum probe density
  edge            Sec. 4  radial (edge) correction, Shoup and Szabo
  pereiro         Sec. 4  static microtiter-well example of Pereiro et al.
  thesis          Sec. 4  saturation time and 1D validity ratio
  fig1            Fig. 1  schematic
Outputs: results/*.json, results/claims.md (every number quoted in the paper), results/figures/*.png|pdf
Requires: see requirements.txt
"""
import argparse, json, os, platform, sys, time
import numpy as np, scipy, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
import stagnant_cell_model as SM
from stagnant_cell_model import solve, wellmixed, theta_eq, t_fraction
import lod_model as LM
from lod_model import BASE, lod, theta_stagnant, theta_wellmixed, theta_equilibrium, lod_depletion_bound, derived, groups

OUT = "results"; FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)
from matplotlib import font_manager as _fm
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator
INK, MUTED, GRID = "#1f2933", "#66727f", "#e4e9ef"
BLUE, ORANGE, TEAL, PURPLE, GOLD = "#2f6db3", "#e0702f", "#1f9d8a", "#8d5fb0", "#d2a82e"
BETA_COL = {0.01: TEAL, 0.1: BLUE, 1.0: ORANGE, 10.0: PURPLE}
CMAP_SEQ = LinearSegmentedColormap.from_list("abec_seq", ["#f4f8fc", "#d3e3f3", "#9cc3e4", "#5f97cf", "#2f6db3", "#1e4784", "#12284f"])
FONT = "Inter" if "Inter" in {f.name for f in _fm.fontManager.ttflist} else "DejaVu Sans"
_math = ({"mathtext.fontset": "custom", "mathtext.rm": FONT, "mathtext.it": FONT + ":italic", "mathtext.bf": FONT + ":bold", "mathtext.sf": FONT, "mathtext.default": "it"}
         if FONT == "Inter" else {"mathtext.fontset": "dejavusans"})
plt.rcParams.update({**_math, "font.family": "sans-serif", "font.sans-serif": [FONT, "DejaVu Sans"], "font.size": 8.5, "axes.labelsize": 8.5,
                     "axes.titlesize": 8.5, "axes.titlelocation": "left", "legend.fontsize": 7.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                     "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": "#4a5560", "axes.linewidth": 0.8,
                     "xtick.color": "#4a5560", "ytick.color": "#4a5560", "xtick.direction": "out", "ytick.direction": "out",
                     "xtick.major.size": 3.0, "ytick.major.size": 3.0, "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
                     "axes.axisbelow": True, "legend.frameon": False, "lines.linewidth": 1.9, "lines.markersize": 4.2,
                     "savefig.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42})
def _enc(o):
    return o.tolist() if hasattr(o, "tolist") else float(o)
CL = {}
def claim(k, v):
    CL[k] = v
    print(f"[claim] {k}: {v}", flush=True)
    with open(os.path.join(OUT, "claims.jsonl"), "a") as fh:
        fh.write(json.dumps({"claim": k, "value": v}, default=_enc) + "\n")
def save(name, obj):
    json.dump(obj, open(os.path.join(OUT, name), "w"), default=_enc)
def load(name):
    return json.load(open(os.path.join(OUT, name)))
def savefig(fig, name, rect=None):
    fig.tight_layout(rect=rect) if rect else fig.tight_layout(); fig.savefig(os.path.join(FIG, name + ".png"), dpi=300); fig.savefig(os.path.join(FIG, name + ".pdf")); plt.close(fig)
def panel(ax, s):
    ax.annotate(s, xy=(0, 1), xycoords="axes fraction", xytext=(-36, 7), textcoords="offset points", fontsize=10.5, fontweight="bold", color=INK, ha="left", va="bottom")

# ------------------------------------------------------------------ helpers
def slowdown(Da, K, beta, fracs=(0.5, 0.9), n=120, tmin=1e-10, stretch=10.0):
    """Time to reach fraction f of equilibrium coverage: stagnant PDE / well-mixed ODE (same sample depletion)."""
    phi = Da / beta; th = theta_eq(K, beta)
    t_end = min(60 * max(1, 1 / (phi * (1 + K))), 5e4)
    t = np.logspace(np.log10(tmin), np.log10(t_end), 700); t[-1] = t_end
    r = solve(Da, K, beta, t_end=t_end, t_eval=t, n=n, stretch=stretch, rtol=1e-7, atol=1e-11)
    wm = wellmixed(K, beta, phi, t)
    out = []
    for f in fracs:
        a = t_fraction(r["t"], r["theta"], f * th); b = t_fraction(t, wm, f * th)
        out.append(a / b if np.isfinite(a) and np.isfinite(b) and b > 0 else np.nan)
    return out

def two_comp(Da, K, beta, t, delta):
    def f(tt, y):
        th = y[0]; ub = 1 - beta * th
        d = min(np.sqrt(np.pi * max(tt, 1e-14)), 1.0) if delta == "grow" else delta
        us = (ub / d + Da * K * th) / (1 / d + Da * (1 - th))
        return [(ub - us) / d / beta]
    sol = solve_ivp(f, (0, t[-1]), [0.0], t_eval=t, method="LSODA", rtol=1e-9, atol=1e-12)
    if not sol.success or sol.y.shape[1] != len(t):
        raise RuntimeError(f"Two-compartment solver failed: {sol.message}")
    return sol.y[0]

# ------------------------------------------------------------------ stages
def stage_validate(fast):
    # V1: reaction-limited limit vs well-mixed ODE
    Da, K, beta = 1e-3, 0.5, 0.1; t = np.logspace(-1, 3, 60)
    r = solve(Da, K, beta, t_end=1e3, t_eval=t); ref = wellmixed(K, beta, Da / beta, t)
    claim("V1 max |dtheta| vs well-mixed ODE", f"{np.max(np.abs(r['theta'] - ref)):.2e}")
    # V2: Cottrell limit
    Da, K, beta = 1e5, 0.0, 1.0; t = np.array([1e-5, 1e-4, 1e-3, 1e-2])
    r = solve(Da, K, beta, t_end=1e-2, t_eval=t, n=400, stretch=12.0)
    ana = (2 / beta) * np.sqrt(t / np.pi); rel = np.abs(r["theta"] - ana) / ana
    claim("V2 Cottrell relative error (min, max)", f"{rel.min():.2e}, {rel.max():.2e}")
    # V3, V4
    r = solve(5.0, 0.3, 0.5, t_end=5.0); claim("V3 mass conservation max |M-1|", f"{np.max(np.abs(r['mass'] - 1)):.1e}")
    r = solve(5.0, 0.3, 0.5, t_end=50.0); claim("V4 |theta_end - theta_eq|", f"{abs(r['theta'][-1] - theta_eq(0.3, 0.5)):.1e}")
    # V5: Laplace-domain solution for reversible linearised kinetics
    import mpmath as mp
    def surf_lin(u0, th, x0, Da, K):
        us = (u0 / x0 + Da * K * th) / (1.0 / x0 + Da); return us, (u0 - us) / x0
    orig = SM.surface_state; SM.surface_state = surf_lin
    worst = 0.0
    for (Da, K, beta) in [(20, 0.5, 5.0), (5, 2.0, 5.0), (200, 1.0, 10.0)]:
        ts = np.array([1e-4, 1e-3, 3e-3, 1e-2, 2e-2])
        r = solve(Da, K, beta, t_end=ts[-1], t_eval=ts, n=400, stretch=12.0, rtol=1e-9, atol=1e-12)
        F = lambda s: Da / (s * (beta * s + Da * beta * mp.sqrt(s) + Da * K))
        ref = np.array([float(mp.invertlaplace(F, t, method="talbot")) for t in ts])
        worst = max(worst, float(np.max(np.abs(r["theta"] - ref) / ref)))
    SM.surface_state = orig
    claim("V5 Laplace benchmark max relative error", f"{worst:.1e}")

def stage_convergence(fast):
    cases = [(1000, 1.78, 0.1), (100, 0.1, 0.1), (10, 31.6, 0.1)]
    worst = {0.5: 0.0, 0.9: 0.0}
    for (Da, K, b) in cases:
        vals = [slowdown(Da, K, b, n=n, stretch=s) for n, s in [(100, 8.0), (150, 10.0), (300, 10.0), (600, 12.0)]]
        for k, f in enumerate((0.5, 0.9)):
            v = np.array([x[k] for x in vals]); worst[f] = max(worst[f], 100 * (v.max() - v.min()) / v[-1])
    claim("grid refinement 100->600 cells, max change of slowdown ratio f=0.5 (%)", f"{worst[0.5]:.2f}")
    claim("grid refinement 100->600 cells, max change of slowdown ratio f=0.9 (%)", f"{worst[0.9]:.2f}")
    p = dict(BASE); ts = 2700.0
    def th_n(c, n):
        Da, K, beta = groups(c, p); tt = p["D"] * ts / p["h"] ** 2
        return solve(Da, K, beta, t_end=tt, t_eval=np.array([tt]), n=n, rtol=1e-8, atol=1e-11)["theta"][-1]
    claim("baseline LOD(45 min), n=110 / 300 / 600 cells (nM)", " / ".join(f"{lod(lambda c: th_n(c, n_), p):.4f}" for n_ in (110, 300, 600)))

def grids(fast):
    return (np.logspace(-2, 3, 6 if fast else 11), np.logspace(-1, 1.5, 3 if fast else 7))

def stage_maps(fast):
    Das, Ks = grids(fast); res = {}
    for b in (0.01, 0.1, 1.0, 10.0):
        P = np.full((2, len(Ks), len(Das)), np.nan)
        for i, K in enumerate(Ks):
            for j, Da in enumerate(Das):
                P[:, i, j] = slowdown(Da, K, b)
        res[str(b)] = P; print("maps beta", b, "done", flush=True)
    save("maps.json", {"Das": Das, "Ks": Ks, "P": res})
    fig_maps()

def fig_maps():
    M = load("maps.json"); Das = np.array(M["Das"]); Ks = np.array(M["Ks"]); P = {k: np.array(v) for k, v in M["P"].items()}
    # claims about the f=0.5 ratio at beta=0.1, K~1.8
    iK = int(np.argmin(abs(Ks - 1.78)))
    row = P["0.1"][0][iK]
    claim("slowdown f=0.5, beta=0.1, K=%.2f: Da grid" % Ks[iK], np.round(Das, 3).tolist())
    claim("slowdown f=0.5, beta=0.1, K=%.2f: ratio" % Ks[iK], np.round(row, 2).tolist())
    claim("slowdown f=0.9, beta=0.1, K=%.2f: ratio" % Ks[iK], np.round(P["0.1"][1][iK], 1).tolist())
    # collapse on Da*beta, beta=0.01 vs 0.1, all K
    step = np.log10(Das[1] / Das[0]); s = int(round(1.0 / step))
    for i, K in enumerate(Ks):
        r = [P["0.01"][0][i, j] / P["0.1"][0][i, j - s] for j in range(s, len(Das))
             if np.isfinite(P["0.01"][0][i, j]) and np.isfinite(P["0.1"][0][i, j - s])]
        claim("collapse beta=0.01 vs 0.1 at equal Da*beta, K=%.2f (min,max ratio)" % K, f"{min(r):.3f}, {max(r):.3f}")
    # beta=1 relative to beta=0.1 at equal Da*beta (K~1.8)
    r1 = []
    for Dab in (1, 3.16, 10, 31.6, 100):
        j1 = int(np.argmin(abs(Das * 0.1 - Dab))); j2 = int(np.argmin(abs(Das * 1.0 - Dab)))
        r1.append(P["1.0"][0][iK, j2] / P["0.1"][0][iK, j1])
    claim("beta=1 / beta=0.1 at equal Da*beta, K~1.8 (range)", f"{min(r1):.2f} to {max(r1):.2f}")
    # Figure 2
    fig, ax = plt.subplots(figsize=(3.9, 3.2))
    Z = P["0.1"][0]; zl = np.log10(Z); lv = np.linspace(0, float(np.nanmax(zl)), 15)
    c = ax.contourf(Das, Ks, zl, levels=lv, cmap=CMAP_SEQ)
    cs = ax.contour(Das, Ks, Z, levels=[1.5, 3, 10], colors=INK, linewidths=1.0)
    for tx in ax.clabel(cs, fmt=lambda v: f"{v:g}\u00d7", fontsize=7.5): tx.set_path_effects([pe.withStroke(linewidth=2.4, foreground="white")])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(False)
    ax.set_xlabel(r"$Da = k_{on}\Gamma_{max}h/D$"); ax.set_ylabel(r"$K = k_{off}/(k_{on}C_0)$")
    cb = fig.colorbar(c, ax=ax, pad=0.03, ticks=np.log10([1, 2, 5, 10, 20])); cb.ax.set_yticklabels(["1", "2", "5", "10", "20"])
    cb.set_label("slowdown ratio"); cb.outline.set_visible(False); cb.ax.grid(False); cb.ax.minorticks_off(); savefig(fig, "fig2_slowdown_map_f0.5_beta0.1")
    # Figure 3
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)
    for b, col in BETA_COL.items():
        y = P[str(b)][0][iK]
        for a, xx in zip(axs, (Das, Das * b)):
            a.loglog(xx, y, "-o", ms=4.2, color=col, mfc="white", mec=col, mew=1.2, label=fr"$\beta$ = {b:g}")
    axs[0].set_ylabel("slowdown ratio, f = 0.5")
    for a, xl in zip(axs, (r"$Da$", r"$Da\cdot\beta$")):
        a.set_xlabel(xl); a.axhline(1.5, color=MUTED, lw=0.8, ls=":")
    axs[0].legend(loc="upper left")
    panel(axs[0], "(a)"); panel(axs[1], "(b)"); savefig(fig, "fig3_beta_collapse")

def stage_scaling(fast):
    Das = (1e2, 1e3, 1e4) if fast else (1e2, 1e3, 1e4, 1e5, 1e6)
    Es, Ps, rows = [], [], []
    for K in (0.3, 1.0, 3.0):
        for Da in Das:
            p = slowdown(Da, K, 0.1, fracs=(0.5,), n=150)[0]; th = theta_eq(K, 0.1)
            E = th ** 2 * (1 + K) * Da * 0.1; Es.append(E); Ps.append(p); rows.append((K, Da, E, p, p / E))
    for r in rows: print("scaling row K=%.1f Da=%.0e E=%.3g ratio=%.3g penalty/E=%.3f" % r, flush=True)
    Es, Ps = np.array(Es), np.array(Ps); m = Es >= 30
    claim("scaling slope (points with E>=30)", f"{np.polyfit(np.log(Es[m]), np.log(Ps[m]), 1)[0]:.2f}")
    for K in (0.3, 1.0, 3.0):
        big = [r[4] for r in rows if r[0] == K and r[1] >= 1e4]
        if big: claim(f"penalty/E at Da>=1e4, K={K}", f"{np.mean(big):.2f}")
    claim("perfect-sink estimate (pi/4) f^2 / (-ln(1-f)) at f=0.5", f"{(np.pi/4)*0.25/(-np.log(0.5)):.3f}")
    save("scaling.json", rows)

def stage_twocomp(fast):
    cases = [(100, 1.0, 0.1), (10, 1.0, 1.0), (1000, 1.0, 0.1), (30, 1.0, 0.01), (100, 0.3, 0.1), (100, 3.0, 0.1)]
    rows = []; fig = plt.figure(figsize=(10.5, 5.8)); gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.05]); axs = [fig.add_subplot(gs[0, k]) for k in range(4)]
    for k, (Da, K, beta) in enumerate(cases):
        phi = Da / beta; th_eq = theta_eq(K, beta); t_end = min(60 * max(1, 1 / (phi * (1 + K))), 5e4)
        t = np.logspace(-8, np.log10(t_end), 300); t[-1] = t_end
        pde = solve(Da, K, beta, t_end=t_end, t_eval=t, n=200, rtol=1e-8, atol=1e-11)["theta"]
        grow = two_comp(Da, K, beta, t, "grow"); eg = float(np.max(np.abs(grow - pde)) / th_eq)
        if k < 4:
            e1 = float(np.max(np.abs(two_comp(Da, K, beta, t, 1.0) - pde)) / th_eq)
            obj = lambda ld: np.max(np.abs(two_comp(Da, K, beta, t, np.exp(ld)) - pde)) / th_eq
            m = minimize_scalar(obj, bounds=(np.log(0.003), np.log(3.0)), method="bounded", options={"xatol": 1e-3})
            db, eb = float(np.exp(m.x)), float(m.fun)
            rows.append((Da, K, beta, 100 * e1, db, 100 * eb, 100 * eg))
            ax = axs[k]
            ax.semilogx(t, pde / th_eq, color=INK, lw=2.3, label="PDE (full solution)")
            ax.semilogx(t, two_comp(Da, K, beta, t, 1.0) / th_eq, "--", color=PURPLE, lw=1.8, label=r"fixed layer, $\delta^*=1$")
            ax.semilogx(t, two_comp(Da, K, beta, t, db) / th_eq, ":", color=ORANGE, lw=2.3, label="best fixed layer")
            ax.semilogx(t, grow / th_eq, "-", color=TEAL, lw=2.0, label=r"growing layer, $\delta^*=\sqrt{\pi t^*}$")
            ax.text(0.97, 0.07, fr"best $\delta^*$ = {db:.2g}", transform=ax.transAxes, ha="right", fontsize=7, color=ORANGE); ax.set_ylim(-0.03, 1.06)
            fr = pde / th_eq; lo = t[np.argmax(fr > 0.01)] / 10; hi = t[np.argmax(fr > 0.999)] * 10 if np.any(fr > 0.999) else t[-1]
            ax.set_xlim(lo, hi)
            ax.set_title(fr"$Da$ = {Da}, $K$ = {K:g}, $\beta$ = {beta:g}"); ax.set_xlabel(r"$t^*$"); panel(ax, "(%s)" % "abcd"[k])
        else:
            rows.append((Da, K, beta, None, None, None, 100 * eg))
    if not os.path.exists(os.path.join(OUT, "closure.json")): stage_closure(fast)
    C = np.array(load("closure.json"), dtype=float)
    axe = fig.add_subplot(gs[1, 0:2]); axr = fig.add_subplot(gs[1, 2:4])
    for b, col in BETA_COL.items():
        e = np.sort(C[C[:, 0] == b, 3]); axe.step(e, np.arange(1, len(e) + 1) / len(e), where="post", color=col, lw=1.8, label=fr"$\beta$ = {b:g}")
    e = np.sort(C[:, 3]); axe.step(e, np.arange(1, len(e) + 1) / len(e), where="post", color=INK, lw=2.4, label="all 36 cases")
    axe.set_xlabel(r"maximum error of the growing-layer model (% of $\theta_{eq}$)"); axe.set_ylabel("ECDF"); axe.legend(loc="lower right"); panel(axe, "(e)")
    q = np.sort(C[:, 4]); axr.step(q, np.arange(1, len(q) + 1) / len(q), where="post", color=INK, lw=2.4); axr.axvline(1.0, color=MUTED, ls=":")
    axr.set_xscale("log"); axr.set_xticks([0.8, 1, 1.5, 2]); axr.set_xticklabels(["0.8", "1", "1.5", "2"]); axr.minorticks_off(); axr.set_xlabel("time to half coverage, growing layer / PDE"); axr.set_ylabel("ECDF"); panel(axr, "(f)")
    axs[0].set_ylabel(r"$\theta/\theta_{eq}$"); fig.legend(*axs[0].get_legend_handles_labels(), loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.0)); savefig(fig, "fig4_two_compartment", rect=(0, 0, 1, 0.94))
    for r in rows: print("two-compartment row (Da,K,beta, err d*=1 %, best d*, err best %, err growing %):", r, flush=True)
    claim("growing-layer error range over six cases (%)", f"{min(r[6] for r in rows):.1f} to {max(r[6] for r in rows):.1f}")
    claim("fixed delta*=1 error range, four cases (%)", f"{min(r[3] for r in rows[:4]):.0f} to {max(r[3] for r in rows[:4]):.0f}")
    claim("best fixed delta* error range / delta* range", f"{min(r[5] for r in rows[:4]):.0f} to {max(r[5] for r in rows[:4]):.0f} %, delta* {min(r[4] for r in rows[:4]):.3f} to {max(r[4] for r in rows[:4]):.3f}")
    save("twocomp.json", rows)


TIMES = [5, 15, 45, 120]
def stage_lod(fast):
    p = dict(BASE); kon, koff, Gm = derived(p)
    claim("baseline k_on (M-1 s-1), Kd (nM), Gamma_max (cm-2), h (cm)", f"{kon/1e3:.3g}, {koff/kon*1e12:.3g}, {Gm*LM.NA:.3g}, {p['h']}")
    stag = [lod(lambda c: theta_stagnant(c, p, t * 60.0), p) for t in TIMES]
    wm = [lod(lambda c: theta_wellmixed(c, p, t * 60.0), p) for t in TIMES]
    eq = lod(lambda c: theta_equilibrium(c, p), p); dep = lod_depletion_bound(p)
    claim("LOD stagnant 5/15/45/120 min (nM)", [round(x, 3) for x in stag]); claim("LOD well-mixed (nM)", [round(x, 3) for x in wm])
    claim("LOD ratio stagnant/well-mixed", [round(a / b, 2) for a, b in zip(stag, wm)])
    claim("LOD equilibrium, depletion bound (nM)", f"{eq:.3f}, {dep:.3f}")
    claim("log-log slope of stagnant LOD, 5 to 120 min", f"{np.log(stag[-1]/stag[0])/np.log(TIMES[-1]/TIMES[0]):.2f}")
    cott = 3 * p["sigma"] * Gm * np.sqrt(np.pi) / (2 * np.sqrt(p["D"] * 2700.0)) * 1e12
    claim("Cottrell LOD estimate at 45 min (nM)", f"{cott:.2f}")
    claim("LOD from K_d alone (nM), 3 sigma threshold", f"{0.06/0.94*koff/kon*1e12:.3f}")
    save("lod.json", {"times": TIMES, "stag": stag, "wm": wm, "eq": eq, "dep": dep})
    C = np.logspace(-2, 3, 26); cols = {5: PURPLE, 15: BLUE, 45: ORANGE, 120: TEAL}
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.1)); ax = axs[0]
    for t, c in cols.items():
        ax.loglog(C, [theta_stagnant(x, p, t * 60.0) for x in C], "-", color=c)
        ax.loglog(C, [theta_wellmixed(x, p, t * 60.0) for x in C], "--", color=c, alpha=0.85)
    ax.axhline(3 * p["sigma"], color=MUTED, lw=0.9, ls=":"); ax.text(C[0] * 1.3, 3 * p["sigma"] * 1.25, r"$3\sigma_\theta$", fontsize=7.5, color=MUTED)
    ax.set_ylim(1e-4, 1.5); ax.set_xlabel("target concentration (nM)"); ax.set_ylabel(r"$\theta(t_{meas})$")
    l1 = ax.legend(handles=[Line2D([0], [0], color=c, lw=2.2, label=f"{t} min") for t, c in cols.items()], title=r"$t_{meas}$", loc="lower right", ncol=2)
    ax.add_artist(l1)
    ax.legend(handles=[Line2D([0], [0], color=INK, lw=1.8, ls="-", label="stagnant"), Line2D([0], [0], color=INK, lw=1.8, ls="--", label="well-mixed")], loc="upper left")
    panel(ax, "(a)"); ax = axs[1]
    ax.loglog(TIMES, stag, "-o", color=ORANGE, mec="white", mew=1.0, ms=6.5, label="stagnant (PDE)")
    ax.loglog(TIMES, wm, "--s", color=TEAL, mec="white", mew=1.0, ms=6, label="well-mixed")
    ax.axhline(eq, color=INK, lw=1.0, ls="-.", label="equilibrium limit"); ax.axhline(dep, color=MUTED, lw=1.1, ls=":", label="sample-depletion bound")
    ax.set_xticks(TIMES); ax.set_xticklabels([str(x) for x in TIMES]); ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel(r"$t_{meas}$ (min)"); ax.set_ylabel("LOD (nM)"); ax.legend(loc="upper right"); panel(ax, "(b)")
    savefig(fig, "fig5_calibration_lod")

def stage_tornado(fast):
    ts = 2700.0; base = dict(BASE)
    ranges = {"D (1.0 to 2.0e-6 cm2/s)": ("D", 1.0e-6, 2.0e-6), "k_on,0 (/3 to x3)": ("kon0", base["kon0"] / 3, base["kon0"] * 3),
              "k_off (/3 to x3)": ("koff", base["koff"] / 3, base["koff"] * 3), "gamma (7.5 to 11.7e-12 mm2)": ("gamma", 7.5e-12, 11.7e-12),
              "S_P (1 to 5e10 mm-2)": ("SP", 1e10, 5e10), "yield (0.5 to 0.9)": ("y", 0.5, 0.9),
              "h (0.25 to 0.99 cm)": ("h", 0.25, 0.99), "sigma (0.01 to 0.05)": ("sigma", 0.01, 0.05)}
    a0 = lod(lambda c: theta_stagnant(c, base, ts), base); b0 = lod(lambda c: theta_wellmixed(c, base, ts), base)
    claim("baseline LOD 45 min stagnant, well-mixed, ratio", f"{a0:.3f}, {b0:.3f}, {a0/b0:.2f}")
    rows = []
    for name, (k, lo, hi) in ranges.items():
        v = []
        for x in (lo, hi):
            q = dict(base); q[k] = x; v.append(lod(lambda c: theta_stagnant(c, q, ts), q))
        rows.append((name, v)); claim("OAT " + name, f"{v[0]:.3f} -> {v[1]:.3f} nM")
    save("tornado.json", rows)
    rs = sorted(rows, key=lambda r: -abs(np.log(r[1][1] / r[1][0]))); fig, ax = plt.subplots(figsize=(5.0, 2.9))
    for i, (n, v) in enumerate(rs[::-1]):
        ax.barh(i, np.log10(v[1] / a0), color=ORANGE, height=0.6, zorder=3); ax.barh(i, np.log10(v[0] / a0), color=TEAL, height=0.6, zorder=3)
    DISP = {"D (1.0 to 2.0e-6 cm2/s)": r"$D$ (1.0 to 2.0 $\times$ 10$^{-10}$ m$^2$ s$^{-1}$)", "k_on,0 (/3 to x3)": r"$k_{on,0}$ ($\div$3 to $\times$3)",
            "k_off (/3 to x3)": r"$k_{off}$ ($\div$3 to $\times$3)", "gamma (7.5 to 11.7e-12 mm2)": r"$\gamma$ (7.5 to 11.7 $\times$ 10$^{-18}$ m$^2$)",
            "S_P (1 to 5e10 mm-2)": r"$S_P$ (1 to 5 $\times$ 10$^{16}$ m$^{-2}$)", "yield (0.5 to 0.9)": r"$y$ (0.5 to 0.9)",
            "h (0.25 to 0.99 cm)": r"$h$ (2.5 to 9.9 mm)", "sigma (0.01 to 0.05)": r"$\sigma_\theta$ (0.01 to 0.05)"}
    ax.set_yticks(range(len(rs))); ax.set_yticklabels([DISP[r[0]] for r in rs[::-1]]); ax.axvline(0, color=INK, lw=0.9, zorder=4); ax.grid(axis="y", visible=False)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=TEAL, label="low end of range"), Patch(color=ORANGE, label="high end of range")], loc="lower right")
    ax.set_xlabel(r"$\log_{10}$(LOD / LOD$_{base}$), 45 min"); savefig(fig, "fig6_tornado")

def stage_noise(fast):
    SREL, SABS = 0.02, 0.02 * 2e10
    def sig(SP, mode):
        return SREL if mode == "relative" else (SABS / SP if mode == "absolute" else np.sqrt((SABS / SP) ** 2 + SREL ** 2))
    SPs = np.array([3e9, 1e10, 2e10, 5e10, 1e11]); tm = [5, 45, 120]; out = {}
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7), sharey=True)
    for ax, mode, lab in zip(axs, ("relative", "absolute", "mixed"), "abc"):
        out[mode] = {}
        for t, c in zip(tm, (PURPLE, ORANGE, TEAL)):
            s, w = [], []
            for SP in SPs:
                q = dict(BASE); q["SP"] = SP; q["sigma"] = sig(SP, mode)
                s.append(lod(lambda x: theta_stagnant(x, q, t * 60.0), q, lo=-4, hi=6)); w.append(lod(lambda x: theta_wellmixed(x, q, t * 60.0), q, lo=-4, hi=6))
            out[mode][t] = (s, w)
            ax.loglog(SPs * 1e6, s, "o-", color=c, ms=5, mec="white", mew=0.8, label=f"stagnant {t} min"); ax.loglog(SPs * 1e6, w, "s--", color=c, ms=4.6, mec="white", mew=0.8, alpha=0.85, label=f"well-mixed {t} min")
            claim(f"noise {mode} {t} min stagnant / well-mixed LOD (nM)", f"{[round(x,3) for x in s]} / {[round(x,3) for x in w]}")
        ax.set_title(f"{mode} noise"); ax.set_xlabel(r"$S_P$ (m$^{-2}$)"); panel(ax, f"({lab})")
    axs[0].set_ylabel("LOD (nM)"); fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=6, fontsize=6, frameon=False); savefig(fig, "fig8_noise_models", rect=(0, 0.08, 1, 1))
    save("noise.json", {m: {str(t): v for t, v in d.items()} for m, d in out.items()})
    ib = int(np.where(SPs == 2e10)[0][0])
    claim("mixed noise: sigma_theta at the baseline S_P of 2e16 m^-2 (sqrt(2) x 0.02)", f"{sig(2e10, 'mixed'):.4f}")
    claim("baseline S_P, stagnant 45 min LOD: sigma = 0.02 (Table 1) / mixed noise (nM)", f"{lod(lambda x: theta_stagnant(x, BASE, 2700.0), BASE, lo=-4, hi=6):.3f} / {out['mixed'][45][0][ib]:.3f}")

def stage_thesis(fast):
    NA = LM.NA; D = 1.5e-6; V = 0.150; A = np.pi * (0.44 / 2) ** 2; h = V / A
    claim("opening area (cm2), h (cm), h^2/D (days)", f"{A:.3f}, {h:.2f}, {h*h/D/86400:.1f}")
    claim("thesis check: fixed parameters", "K_d = 0.5 nM, C0 = 2 uM, D = 150 um2/s, k_off = K_d k_on, h from 150 uL over the 4.4 mm opening")
    rough = 3.76e-4 / 1.96e-5   # electrochemically active area of nanoporous gold / planar gold in the thesis (cyclic voltammetry)
    for tag, rf in (("planar electrode", 1.0), (f"roughness factor {rough:.1f} (nanoporous / planar active area), pore diffusion ignored", rough)):
        t90, combos = [], []
        for kon in (9.34e4, 9.34e3):
            for G in (1e12, 5e12):
                C0 = 2000e-12; Gm = rf * G / NA; k = kon * 1e3; ko = k * 0.5e-9 * 1e-3
                Da = k * Gm * h / D; K = ko / (k * C0); beta = Gm / (C0 * h); th = theta_eq(K, beta)
                tau = h * h / D; phi = Da / beta; t_end = 60 * max(1, 1 / (phi * (1 + K)))
                t = np.logspace(-9, np.log10(t_end), 600); t[-1] = t_end
                r = solve(Da, K, beta, t_end=t_end, t_eval=t, n=150, rtol=1e-7, atol=1e-11)
                t90.append(t_fraction(r["t"], r["theta"], 0.9 * th) * tau); combos.append((kon, G / 100.0, Da * beta))
        claim(f"thesis check, {tag}: time to 90% of equilibrium coverage; (k_on in 1/(M s), Gamma_max per projected area in m^-2, Da*) -> seconds",
              "; ".join(f"({k:.3g}, {g*1e6:.3g}, {d:.3g}) -> {x:.1f}" for (k, g, d), x in zip(combos, t90)))
    for tm in (5, 15, 45, 120, 480, 1440):
        claim(f"1D validity ratio sqrt(Dt)/a at {tm} min", f"{np.sqrt(D*tm*60)/0.22:.2f}")

def stage_fig1(fast):
    from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
    from scipy.special import erf
    fig, ax = plt.subplots(figsize=(6.2, 4.6)); ax.set_xlim(-0.6, 12.2); ax.set_ylim(-2.7, 8.9); ax.set_aspect("equal"); ax.axis("off"); ax.grid(False)
    x0, x1, H, dl = 1.3, 6.7, 7.2, 2.0
    yy = np.linspace(0, H, 240); conc = erf(yy / (0.75 * dl))
    ax.imshow(conc[:, None], extent=(x0, x1, 0, H), origin="lower", aspect="auto", zorder=0, vmin=0, vmax=1,
              cmap=LinearSegmentedColormap.from_list("sol", ["#ffffff", "#dcecf7"]))
    ax.add_patch(Rectangle((x0, 0), x1 - x0, H, fill=False, ec=INK, lw=1.3, zorder=3))
    ax.plot([x0, x1], [H, H], color=INK, lw=3.4, solid_capstyle="butt", zorder=4)
    ax.add_patch(FancyBboxPatch((x0 - 0.35, -0.75), x1 - x0 + 0.7, 0.75, boxstyle="round,pad=0,rounding_size=0.12", fc=GOLD, ec="#a8861c", lw=1.0, zorder=4))
    ax.text((x0 + x1) / 2, -0.38, "planar gold electrode  (x = 0)", ha="center", va="center", fontsize=8, color=INK, zorder=5)
    rng = np.random.default_rng(7); n = 0
    while n < 70:
        x = rng.uniform(x0 + 0.2, x1 - 0.2); y = rng.uniform(1.2, H - 0.3)
        if rng.random() < conc[int(y / H * 239)] ** 1.3:
            ax.plot(x, y, "o", ms=4.4, color=TEAL, mec="white", mew=0.5, zorder=2); n += 1
    xs = np.linspace(x0 + 0.45, x1 - 0.45, 8)
    for k, xp in enumerate(xs):
        ax.plot([xp, xp], [0, 0.72], color="#4a5560", lw=1.3, zorder=3, solid_capstyle="round")
        ax.plot(xp, 0.82, "o", ms=5.2, color=ORANGE, mec="white", mew=0.6, zorder=4)
        if k in (1, 4, 5, 7): ax.plot(xp, 1.22, "o", ms=4.4, color=TEAL, mec="white", mew=0.5, zorder=4)
    ax.annotate("", xy=(0.45, H), xytext=(0.45, 0), arrowprops=dict(arrowstyle="<->", lw=1.1, color=INK, shrinkA=0, shrinkB=0))
    ax.text(0.2, H / 2, "h", fontsize=11, ha="right", va="center", style="italic")
    ax.annotate("no-flux wall, x = h", xy=((x0 + x1) / 2, H + 0.05), xytext=((x0 + x1) / 2, H + 0.85), ha="center", fontsize=8, color=INK)
    ax.text((x0 + x1) / 2, H - 0.75, r"target in solution: $C_0$, $D$", ha="center", fontsize=8, zorder=5, bbox=dict(fc="white", ec="none", pad=1.6, alpha=0.9))
    ax.annotate("", xy=(x1 + 0.35, dl), xytext=(x1 + 0.35, 0), arrowprops=dict(arrowstyle="<->", lw=1.1, color=BLUE, shrinkA=0, shrinkB=0))
    ax.text(x1 + 0.6, dl / 2, "depletion\nlayer\n" + r"$\delta\approx\sqrt{\pi D t}$", fontsize=7.8, color=BLUE, va="center", ha="left", linespacing=1.25)
    px0, px1 = 9.6, 11.9
    ax.plot([px0, px0], [0, H], color="#4a5560", lw=0.9); ax.plot([px0, px1 + 0.3], [0, 0], color="#4a5560", lw=0.9)
    ax.plot(px0 + (px1 - px0) * conc, yy, color=BLUE, lw=2.2, zorder=3)
    ax.plot([px1, px1], [0, H], color=MUTED, lw=0.9, ls=(0, (3, 2)))
    ax.text(px1, H + 0.35, r"$C_0$", fontsize=8.5, ha="center", color=MUTED); ax.text(px0 + 0.1, H + 0.35, r"$C(x,t)$", fontsize=8.5, ha="left", color=BLUE)
    ax.text(px0 + (px1 - px0) / 2, -0.55, r"concentration", fontsize=7.5, ha="center", va="center", color=MUTED)
    ax.annotate("x", xy=(8.9, 2.6), xytext=(8.9, 0.0), arrowprops=dict(arrowstyle="->", lw=1.1, color=INK, shrinkA=0, shrinkB=0), ha="center", va="bottom", fontsize=9, style="italic")
    ax.plot([x1 + 0.35, 8.9], [0, 0], color=MUTED, lw=0.6, ls=":")
    ax.plot(1.05, -1.6, "o", ms=5.2, color=ORANGE, mec="white"); ax.text(1.35, -1.6, "probe", fontsize=7.6, va="center")
    ax.plot(2.75, -1.6, "o", ms=4.6, color=TEAL, mec="white"); ax.text(3.05, -1.6, "target", fontsize=7.6, va="center")
    ax.text(4.4, -1.6, r"$\Gamma_{max}$, $k_{on}$, $k_{off}$", fontsize=7.6, va="center")
    ax.text((x0 + x1) / 2 + 2.6, -2.35, r"$D\,\partial C/\partial x\,|_{x=0} = k_{on}C_s(\Gamma_{max}-\Gamma) - k_{off}\Gamma$", ha="center", va="center", fontsize=8.2, color=INK,
            bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.5", fc="#f3f6f9", ec="#cdd5de", lw=0.8))
    savefig(fig, "fig1_cell_schematic")


def r_pred(K, f):
    """Series-driving-force estimate of penalty/(theta_eq^2 (1+K) Da*) for beta -> 0 and Da -> infinity:
    Cottrell balance with reversible driving force C0 - K_d*theta/(1-theta)."""
    a = 1 + K
    return (np.pi / 4) * (f - K * np.log(1 - f)) ** 2 / (a ** 2 * (-np.log(1 - f)))

def stage_scaling2(fast):
    rows = []
    for beta in (0.01, 0.1, 1.0):
        for K in (0.3, 1.0, 3.0):
            for f in (0.3, 0.5, 0.7):
                Da = 1e6; th = theta_eq(K, beta)
                pen = slowdown(Da, K, beta, fracs=(f,), n=150)[0]
                E = th ** 2 * (1 + K) * Da * beta
                rows.append((beta, K, f, pen / E, r_pred(K, f)))
                print("scaling2 beta=%.2f K=%.1f f=%.1f  obs=%.3f  pred=%.3f  obs/pred=%.2f" % (beta, K, f, pen / E, r_pred(K, f), pen / E / r_pred(K, f)), flush=True)
    save("scaling2.json", rows)
    for beta in (0.01, 0.1, 1.0):
        sub = [r for r in rows if r[0] == beta]
        claim(f"scaling2 obs/pred range, beta={beta}", f"{min(r[3]/r[4] for r in sub):.2f} to {max(r[3]/r[4] for r in sub):.2f}")
    sub = [r for r in rows if r[0] == 0.01]
    claim("scaling2 prefactor obs range over K,f at beta=0.01", f"{min(r[3] for r in sub):.2f} to {max(r[3] for r in sub):.2f}")

def stage_kpeak(fast):
    """Eq. (4) implies, at fixed Da, a slowdown proportional to r(K,f)/(1+K). For f = 0.5 this decreases with K,
    for f = 0.9 the prefactor r grows faster than 1/(1+K) falls, so the slowdown has a maximum at intermediate K."""
    Kg = np.logspace(-1.5, 1.5, 601)
    for f in (0.5, 0.9):
        v = r_pred(Kg, f) / (1 + Kg); i = int(np.argmax(v))
        claim(f"Eq. (4) slowdown ~ r/(1+K), f={f}: maximum at K", f"{Kg[i]:.2f}" + (" (lower grid edge: decreasing in K)" if i == 0 else " (interior maximum)"))
    if os.path.exists(os.path.join(OUT, "maps.json")):
        M = load("maps.json"); Ks = np.array(M["Ks"]); P = np.array(M["P"]["0.1"])
        for k, f in enumerate((0.5, 0.9)):
            claim(f"PDE slowdown at Da={M['Das'][-1]:.0f}, beta=0.1, f={f}: K of the maximum on the map grid", f"{Ks[int(np.argmax(P[k][:, -1]))]:.2f}")

def stage_variance_split(fast):
    """Variance of ln(LOD) for the Cottrell estimate LOD ~ sigma * y * S_P / sqrt(D), with the log-uniform ranges of the Sobol analysis."""
    b = BASE
    R = {"D": (1.0e-6, 2.0e-6, -0.5), "k_on,0": (b["kon0"] / 3, b["kon0"] * 3, 0.0), "k_off": (b["koff"] / 3, b["koff"] * 3, 0.0),
         "gamma": (7.5e-12, 11.7e-12, 0.0), "S_P": (1e10, 5e10, 1.0), "y": (0.5, 0.9, 1.0), "h": (0.25, 0.99, 0.0), "sigma": (0.01, 0.05, 1.0)}
    v = {k: e ** 2 * np.log(hi / lo) ** 2 / 12 for k, (lo, hi, e) in R.items()}; tot = sum(v.values())
    claim("Cottrell estimate: share of var(ln LOD) from sigma / S_P / y / D (%)", " / ".join(f"{100 * v[k] / tot:.1f}" for k in ("sigma", "S_P", "y", "D")))

def stage_closure(fast):
    rows = []
    for beta in (0.01, 0.1, 1.0, 10.0):
        for K in (0.1, 1.0, 10.0):
            for Da in (10, 100, 1000):
                phi = Da / beta; th_eq = theta_eq(K, beta); t_end = min(60 * max(1, 1 / (phi * (1 + K))), 5e4)
                t = np.logspace(-9, np.log10(t_end), 400); t[-1] = t_end
                pde = solve(Da, K, beta, t_end=t_end, t_eval=t, n=120, rtol=1e-7, atol=1e-11)["theta"]
                cl = two_comp(Da, K, beta, t, "grow")
                err = float(np.max(np.abs(cl - pde)) / th_eq)
                t5p = t_fraction(t, pde, 0.5 * th_eq); t5c = t_fraction(t, cl, 0.5 * th_eq)
                rows.append((beta, K, Da, 100 * err, t5c / t5p))
    save("closure.json", rows); e = np.array([r[3] for r in rows]); q = np.array([r[4] for r in rows])
    claim("closure delta=sqrt(pi t): cases", len(rows))
    claim("closure max error % of theta_eq: median, 90th pct, max", f"{np.median(e):.1f}, {np.percentile(e,90):.1f}, {e.max():.1f}")
    claim("closure t50 ratio closure/PDE: min, median, max", f"{q.min():.2f}, {np.median(q):.2f}, {q.max():.2f}")
    worst = sorted(rows, key=lambda r: -r[3])[:4]
    claim("closure worst four cases (beta,K,Da,err%)", [tuple(round(x, 2) for x in r[:4]) for r in worst])

def _sobol_qoi(x):
    p = dict(BASE); p.update(dict(zip(("D", "kon0", "koff", "gamma", "SP", "y", "h", "sigma"), x)))
    return float(np.log10(lod(lambda c: theta_stagnant(c, p, 2700.0), p, lo=-4, hi=5, it=14)))

def stage_ratio_range(fast):
    """Stagnant / well-mixed detection-limit ratio at 45 min over the noise and probe-density ranges of the sensitivity analysis."""
    rows = []
    for sg in (0.01, 0.02, 0.05):
        for SP in (1e10, 2e10, 5e10):
            p = dict(BASE); p["sigma"] = sg; p["SP"] = SP
            a = lod(lambda c: theta_stagnant(c, p, 2700.0), p, lo=-4, hi=6); b = lod(lambda c: theta_wellmixed(c, p, 2700.0), p, lo=-4, hi=6)
            rows.append((sg, SP, a, b, a / b))
            print(f"ratio_range sigma={sg} S_P={SP:.0e}: stagnant {a:.3f} nM, well-mixed {b:.3f} nM, ratio {a/b:.2f}", flush=True)
    save("ratio_range.json", rows)
    r = [x[4] for x in rows]
    claim("stagnant/well-mixed LOD ratio at 45 min over sigma = 0.01-0.05 and S_P = 1e16-5e16 m^-2: min, max", f"{min(r):.2f}, {max(r):.2f}")
    claim("ratio at the baseline S_P = 2e16 m^-2 over sigma = 0.01, 0.02, 0.05", ", ".join(f"{x[4]:.2f}" for x in rows if x[1] == 2e10))

def stage_convection(fast):
    """Order-of-magnitude check that the column stays stagnant (water at 23 C, 150 mM NaCl taken as water)."""
    g, rho, nu, kappa, beta_T, D = 9.81, 1000.0, 0.95e-6, 1.43e-7, 2.3e-4, 1.5e-10
    M = 7.6; h = 9.9e-3                       # kg/mol for a 23-mer (about 330 g/mol per nucleotide); column height in m
    for C0_nM in (2.4, 10.0):
        drho = C0_nM * 1e-6 * M               # kg/m3, upper bound: the depleted layer is free of target
        for tmin in (15, 45, 120):
            dl = np.sqrt(np.pi * D * tmin * 60.0)
            Ra = g * drho * dl ** 3 / (rho * nu * D)
            print(f"convection C0={C0_nM} nM t={tmin} min: depletion layer {dl*1e3:.2f} mm, solutal Rayleigh number {Ra:.2f}", flush=True)
        dl = np.sqrt(np.pi * D * 45 * 60.0); claim(f"solutal Rayleigh number of the depletion layer, C0 = {C0_nM} nM, 45 min", f"{g * drho * dl ** 3 / (rho * nu * D):.1f}")
    for dT in (0.01, 0.1, 1.0):
        claim(f"thermal Rayleigh number of the {h*1e3:.1f} mm column for a vertical temperature difference of {dT} K", f"{g * beta_T * dT * h ** 3 / (nu * kappa):.0f}")
    claim("critical Rayleigh numbers for comparison (textbook values, rigid walls)", "about 1.7e3 (fluid heated from below)")

def stage_sobol(fast, budget=float("inf")):
    from scipy.stats import qmc
    names = ["D", "k_on,0", "k_off", "gamma", "S_P", "y", "h", "sigma"]
    b = BASE; lo = np.array([1.0e-6, b["kon0"] / 3, b["koff"] / 3, 7.5e-12, 1e10, 0.5, 0.25, 0.01])
    hi = np.array([2.0e-6, b["kon0"] * 3, b["koff"] * 3, 11.7e-12, 5e10, 0.9, 0.99, 0.05])
    N = 16 if fast else 64; d = len(names)
    U = qmc.Sobol(d=2 * d, scramble=True, seed=7).random(N); A, B = U[:, :d], U[:, d:]
    rows = [A, B] + [np.where(np.arange(d) == i, B, A) for i in range(d)]
    X = np.vstack(rows); X = np.exp(np.log(lo) + X * (np.log(hi) - np.log(lo)))
    path = os.path.join(OUT, "sobol_evals.npy")
    y = np.load(path) if os.path.exists(path) and not fast else np.full(len(X), np.nan)
    t0 = time.time()
    for i in range(len(X)):
        if np.isnan(y[i]):
            y[i] = _sobol_qoi(X[i])
            if time.time() - t0 > budget: break
        if i % 50 == 0:
            print("sobol", i, "of", len(X), flush=True)
            if not fast: np.save(path, y)
    if not fast: np.save(path, y)
    done = int(np.sum(~np.isnan(y))); print(f"sobol evaluations done: {done} of {len(X)}", flush=True)
    if done < len(X): print("INCOMPLETE: run the stage again to continue", flush=True); return
    fA, fB = y[:N], y[N:2 * N]; V = np.var(np.concatenate([fA, fB]))
    S1, ST = [], []
    for i in range(d):
        fAB = y[(2 + i) * N:(3 + i) * N]
        S1.append(float(np.mean(fB * (fAB - fA)) / V)); ST.append(float(0.5 * np.mean((fA - fAB) ** 2) / V))
    claim("Sobol N, evaluations, variance of log10 LOD", f"{N}, {len(X)}, {V:.4f}")
    rng = np.random.default_rng(11); NB = 2000; S1b = np.zeros((NB, d)); STb = np.zeros((NB, d))
    for k in range(NB):
        idx = rng.integers(0, N, N); a_, b_ = fA[idx], fB[idx]; Vb = np.var(np.concatenate([a_, b_]))
        for i in range(d):
            ab = y[(2 + i) * N:(3 + i) * N][idx]
            S1b[k, i] = np.mean(b_ * (ab - a_)) / Vb; STb[k, i] = 0.5 * np.mean((a_ - ab) ** 2) / Vb
    l1, h1 = np.percentile(S1b, [2.5, 97.5], axis=0); lT, hT = np.percentile(STb, [2.5, 97.5], axis=0)
    for n_, s1, st, a1, b1, aT, bT in zip(names, S1, ST, l1, h1, lT, hT):
        claim(f"Sobol {n_}: S1 [95% CI], ST [95% CI]", f"{s1:.2f} [{a1:.2f}, {b1:.2f}], {st:.2f} [{aT:.2f}, {bT:.2f}]")
    claim("Sobol sum of first-order indices", f"{sum(S1):.2f}")
    save("sobol.json", {"names": names, "S1": S1, "ST": ST, "S1_ci": [l1, h1], "ST_ci": [lT, hT], "N": N})
    fig, ax = plt.subplots(figsize=(5.0, 3.2)); yy = np.arange(d)[::-1]
    ax.errorbar(S1, yy + 0.12, xerr=[np.array(S1) - l1, h1 - np.array(S1)], fmt="o", color=BLUE, capsize=2.5, elinewidth=1.3, ms=5.5, mec="white", mew=0.8, label="first order")
    ax.errorbar(ST, yy - 0.12, xerr=[np.array(ST) - lT, hT - np.array(ST)], fmt="s", color=ORANGE, capsize=2.5, elinewidth=1.3, ms=5.5, mec="white", mew=0.8, label="total")
    ax.axvline(0, color=INK, lw=0.9); ax.grid(axis="y", visible=False); ax.set_yticks(yy); ax.set_yticklabels([r"$D$", r"$k_{on,0}$", r"$k_{off}$", r"$\gamma$", r"$S_P$", r"$y$", r"$h$", r"$\sigma_\theta$"]); ax.set_xlabel("Sobol index of log$_{10}$ LOD (95% bootstrap interval)"); ax.legend(loc="lower right"); savefig(fig, "fig7_sobol")

def stage_optimum(fast):
    SREL = 0.02; SPs = np.logspace(9.3, 11.0, 6 if fast else 14); rows = {}
    fig, ax = plt.subplots(figsize=(4.4, 2.9))
    h_s, h_w = [], []
    for rho, col in ((0.25, TEAL), (1.0, ORANGE), (4.0, PURPLE)):
        SABS = rho * 0.02 * 2e10; ls, lw = [], []
        for SP in SPs:
            q = dict(BASE); q["SP"] = SP; q["sigma"] = float(np.sqrt((SABS / SP) ** 2 + SREL ** 2))
            ls.append(lod(lambda c: theta_stagnant(c, q, 2700.0), q, lo=-4, hi=6)); lw.append(lod(lambda c: theta_wellmixed(c, q, 2700.0), q, lo=-4, hi=6))
        ls, lw = np.array(ls), np.array(lw)
        def opt(y):
            ok = np.isfinite(y); xs, ys = SPs[ok], y[ok]
            i = int(np.argmin(ys)); j0, j1 = max(i - 2, 0), min(i + 3, len(ys))
            c = np.polyfit(np.log10(xs[j0:j1]), np.log10(ys[j0:j1]), 2)
            x = -c[1] / (2 * c[0]) if (c[0] > 0 and j1 - j0 >= 3) else np.log10(xs[i])
            return 10 ** x, float(ys[i])
        os_, ws_ = opt(ls), opt(lw); rows[str(rho)] = (os_, ws_)
        claim(f"optimum S_P (m-2) and min LOD (nM) at 45 min, sigma_a/sigma_r = {rho}: stagnant | well-mixed", f"{os_[0]*1e6:.2g}, {os_[1]:.2f} | {ws_[0]*1e6:.2g}, {ws_[1]:.2f}")
        h1, = ax.loglog(SPs * 1e6, ls, "o-", ms=4.6, color=col, mec="white", mew=0.8, label=fr"stagnant, $\sigma_a/\sigma_r$ = {rho:g}")
        h2, = ax.loglog(SPs * 1e6, lw, "s--", ms=4.2, color=col, mec="white", mew=0.8, alpha=0.85, label=fr"well-mixed, $\sigma_a/\sigma_r$ = {rho:g}")
        h_s.append(h1); h_w.append(h2)
        for (xo, yo), arr in ((os_, ls), (ws_, lw)):
            ok = np.isfinite(arr); yi = 10 ** np.interp(np.log10(xo), np.log10(SPs[ok]), np.log10(arr[ok]))
            ax.plot(xo * 1e6, yi, marker="*", ms=11, color=col, mec="white", mew=0.9, zorder=6, ls="none")
    ax.set_xlabel(r"$S_P$ (m$^{-2}$)"); ax.set_ylabel("LOD at 45 min (nM)"); ax.set_ylim(top=250); ax.legend(h_s + h_w, [h.get_label() for h in h_s + h_w], ncol=2, loc="upper center", fontsize=6.8); savefig(fig, "fig9_optimum_probe_density")
    save("optimum.json", rows)

def stage_edge(fast):
    import edge_effect as EE
    from scipy.integrate import quad
    out = []
    for tm in (5, 15, 45, 120, 480, 1440):
        t = tm * 60.0; Nd = quad(EE.i_disk, 0, t, limit=200)[0]; Np = quad(EE.i_planar, 0, t, limit=200)[0]
        out.append((tm, Nd / Np)); claim(f"edge enhancement N_disk/N_planar at {tm} min (laterally unbounded sample)", f"{Nd/Np:.3f}")
    save("edge.json", out)

def stage_pereiro(fast):
    """Literature illustration: static 96-well IgG assay, Box 3 of Pereiro et al. (2020); parameters from that paper (no measured curves used)."""
    NA = LM.NA; D = 6.5e-11; kon = 1000.0; koff = 1e-3; c0 = 13e-6
    V = 100e-9; A = np.pi * (3.3e-3) ** 2; h = V / A; G = 18000 * 1e12 / NA
    Da = kon * G * h / D; K = koff / (kon * c0); beta = G / (c0 * h); th = theta_eq(K, beta); tau = h * h / D
    claim("Pereiro Box 3 static case: h (mm), Da, K, beta, Da*", f"{h*1e3:.2f}, {Da:.0f}, {K:.3f}, {beta:.2f}, {Da*beta:.0f}")
    tR = 1 / (kon * c0 + koff)
    claim("Pereiro Box 3: Da * delta/h with delta = sqrt(2 D t_R) (their quoted Da = 44)", f"{Da*np.sqrt(2*D*tR)/h:.1f}")
    phi = Da / beta; t_end = 60 * max(1, 1 / (phi * (1 + K))); t = np.logspace(-9, np.log10(t_end), 800); t[-1] = t_end
    r = solve(Da, K, beta, t_end=t_end, t_eval=t, n=200, rtol=1e-8, atol=1e-11); wm = wellmixed(K, beta, phi, t)
    a = t_fraction(r["t"], r["theta"], 0.5 * th) * tau; b = t_fraction(t, wm, 0.5 * th) * tau
    E = th ** 2 * (1 + K) * Da * beta
    claim("Pereiro Box 3: t50 stagnant (h), well-mixed (s), ratio", f"{a/3600:.2f}, {b:.0f}, {a/b:.0f}")
    claim("Pereiro Box 3: observed prefactor ratio/E vs Eq. (4) r(K, 0.5)", f"{a/b/E:.3f} vs {r_pred(K, 0.5):.3f}")

STAGES = {"validate": stage_validate, "convergence": stage_convergence, "maps": stage_maps, "scaling": stage_scaling, "scaling2": stage_scaling2,
          "kpeak": stage_kpeak, "closure": stage_closure, "twocomp": stage_twocomp, "lod": stage_lod, "tornado": stage_tornado, "ratio_range": stage_ratio_range, "convection": stage_convection,
          "sobol": stage_sobol, "variance_split": stage_variance_split, "noise": stage_noise, "optimum": stage_optimum,
          "edge": stage_edge, "pereiro": stage_pereiro, "thesis": stage_thesis, "fig1": stage_fig1}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--fast", action="store_true"); ap.add_argument("--stages", default="all")
    a = ap.parse_args()
    print("python", platform.python_version(), "numpy", np.__version__, "scipy", scipy.__version__, "matplotlib", matplotlib.__version__, flush=True)
    names = list(STAGES) if a.stages == "all" else a.stages.split(",")
    if a.stages == "all" and os.path.exists(os.path.join(OUT, "claims.jsonl")): os.remove(os.path.join(OUT, "claims.jsonl"))
    t0 = time.time()
    for n in names:
        s = time.time(); print(f"=== stage {n}", flush=True); STAGES[n](a.fast); print(f"=== {n} done in {time.time()-s:.0f} s", flush=True)
    rows = list({r["claim"]: r for r in (json.loads(l) for l in open(os.path.join(OUT, "claims.jsonl")))}.values())
    with open(os.path.join(OUT, "claims.md"), "w") as f:
        f.write("| claim | value |\n|---|---|\n")
        for r in rows: f.write(f"| {r['claim']} | {r['value']} |\n")
    print(f"all done in {time.time()-t0:.0f} s", flush=True)
