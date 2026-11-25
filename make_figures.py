"""Build the paper figures from results/*.json. Usage: python make_figures.py 1 2 ..."""
import os, sys, json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Circle, Polygon
from matplotlib.colors import to_rgb
from scipy.interpolate import RectBivariateSpline
import figstyle as S
from figstyle import CM, INK, MUTED, GRID

HERE = os.path.dirname(os.path.abspath(__file__))
def load(name): return json.load(open(os.path.join(HERE, "results", name)))


# ----------------------------------------------------------------------------------------------
def fig1():
    from stagnant_cell_model import solve, surface_state
    W, H = 16 * CM, 7.2 * CM
    fig = plt.figure(figsize=(W, H))
    A = fig.add_axes([0, 0, 1, 1]); A.set_xlim(0, W); A.set_ylim(0, H); A.axis("off")
    cx0, cx1, y0, y1 = 1.60, 2.36, 0.95, 2.60            # column and electrode surface (inches)
    dep = 0.40                                             # drawn thickness of the depletion layer
    BLUE, TEAL, VERM = "#1D6FA3", "#009E73", "#D55E00"

    # solution
    A.add_patch(Rectangle((cx0, y0), cx1 - cx0, y1 - y0, fc="#EEF6FB", ec="none", zorder=1))
    grad = np.zeros((64, 1, 4)); grad[..., :3] = to_rgb("#BFDFF2"); grad[:, 0, 3] = np.linspace(0.95, 0.0, 64)
    A.imshow(grad, extent=(cx0, cx1, y0, y0 + dep), origin="lower", aspect="auto", zorder=2, interpolation="bilinear")
    A.plot([cx0, cx1], [y0 + dep] * 2, color="#4A9CCF", lw=0.8, ls=(0, (3, 2)), zorder=3)
    for x in (cx0, cx1):
        A.plot([x, x], [y0, y1], color="#94A3B8", lw=1.3, zorder=4, solid_capstyle="butt")
    # no-flux wall with hatching
    A.plot([cx0 - 0.03, cx1 + 0.03], [y1, y1], color="#475569", lw=2.2, zorder=5, solid_capstyle="butt")
    for xs in np.arange(cx0 - 0.03, cx1 + 0.03, 0.065):
        A.plot([xs, xs + 0.045], [y1, y1 + 0.055], color="#94A3B8", lw=0.7, zorder=5)
    A.text((cx0 + cx1) / 2, y1 + 0.12, r"no-flux wall, $x=h$", ha="center", va="bottom", fontsize=7)
    # electrode (gold gradient)
    gold = np.zeros((64, 1, 3)); top, bot = np.array(to_rgb("#EFCB72")), np.array(to_rgb("#C4921A"))
    for i, t in enumerate(np.linspace(0, 1, 64)): gold[i, 0] = top * (1 - t) + bot * t
    A.imshow(gold, extent=(1.50, 2.46, 0.78, y0), origin="upper", aspect="auto", zorder=3)
    A.add_patch(Rectangle((1.50, 0.78), 0.96, y0 - 0.78, fc="none", ec="#8A6A12", lw=0.6, zorder=4))
    # probes and bound targets
    xs = np.linspace(cx0 + 0.07, cx1 - 0.07, 8); bound = {1, 3, 4, 6}
    for i, x in enumerate(xs):
        A.plot([x, x], [y0, y0 + 0.19], color="#374151", lw=1.1, zorder=6, solid_capstyle="round")
        A.add_patch(Circle((x, y0 + 0.215), 0.030, fc=VERM, ec="white", lw=0.5, zorder=7))
        if i in bound:
            A.add_patch(Circle((x, y0 + 0.285), 0.027, fc=TEAL, ec="white", lw=0.5, zorder=8))
    # free targets: uniform in the bulk, thinning toward the electrode inside the depletion layer
    rng = np.random.default_rng(1); pts = []
    while len(pts) < 46:
        x = rng.uniform(cx0 + 0.05, cx1 - 0.05); y = rng.uniform(y0 + 0.36, y1 - 0.05)
        p = 1.0 if y > y0 + 0.62 else ((y - (y0 + 0.36)) / 0.26) ** 1.2
        if rng.uniform() > p: continue
        if all((x - a) ** 2 + (y - b) ** 2 > 0.075 ** 2 for a, b in pts): pts.append((x, y))
    for x, y in pts: A.add_patch(Circle((x, y), 0.027, fc=TEAL, ec="white", lw=0.4, zorder=6))

    # dimension arrows on the left
    def dim(x, ya, yb, color, lw=0.9):
        A.add_patch(FancyArrowPatch((x, ya), (x, yb), arrowstyle="<|-|>", mutation_scale=6, lw=lw, color=color, zorder=9, shrinkA=0, shrinkB=0))
    dim(0.30, y0, y1, MUTED)
    A.text(0.20, (y0 + y1) / 2, "h", ha="right", va="center", fontsize=9, style="italic")
    for yy in (y0, y1): A.plot([0.30, cx0 - 0.02], [yy, yy], color="#CBD5E1", lw=0.6, ls=":", zorder=1)
    dim(1.40, y0, y0 + dep, BLUE)
    A.plot([1.40, cx0], [y0 + dep] * 2, color="#9CC4DF", lw=0.6, ls=":", zorder=1)
    A.text(1.33, y0 + dep / 2, r"$\delta$", ha="right", va="center", fontsize=9, color=BLUE)
    # x axis arrow
    A.add_patch(FancyArrowPatch((0.68, y0), (0.68, y0 + 0.62), arrowstyle="-|>", mutation_scale=7, lw=1.0, color=INK, zorder=9, shrinkA=0, shrinkB=0))
    A.text(0.68, y0 + 0.69, r"$x$", ha="center", va="bottom", fontsize=9)
    A.plot([0.60, 0.76], [y0, y0], color=INK, lw=0.8)
    A.text(0.76, y0 - 0.025, r"$x=0$", ha="left", va="top", fontsize=6.5, color=MUTED)

    # right-hand labels with leaders
    LX = 2.74
    def lab(txt, y, tx, ty, color=INK, **kw):
        A.text(LX, y, txt, ha="left", va="center", fontsize=7, color=color, **kw)
        A.plot([LX - 0.03, tx], [y, ty], color="#9CA3AF", lw=0.6, zorder=10)
        A.add_patch(Circle((tx, ty), 0.012, fc="#6B7280", ec="none", zorder=10))
    lab(r"bulk: $C_0$, $D$", 2.30, pts[0][0] if pts[0][0] > 2.0 else 2.25, 2.30 if False else pts[0][1])
    lab(r"depletion layer, $\delta\approx\sqrt{\pi D t}$", 1.52, cx1 - 0.02, y0 + 0.30, color=BLUE)
    lab(r"probes: $\Gamma_{max}$, $k_{on}$, $k_{off}$", 1.33, xs[-1], y0 + 0.14)
    lab(r"$C_s$ (surface)", 1.14, cx1 - 0.02, y0 + 0.02)
    lab("gold electrode", 0.87, 2.40, 0.87)
    # reaction inset
    A.add_patch(FancyBboxPatch((LX, 1.68), 1.04, 0.50, boxstyle="round,pad=0.0,rounding_size=0.05", fc="white", ec="#CBD5E1", lw=0.8, zorder=2))
    A.text(LX + 0.52, 2.07, r"$k_{on}$", ha="center", va="center", fontsize=7, color=VERM)
    A.text(LX + 0.52, 1.93, r"$\mathrm{T}+\mathrm{P}\;\rightleftharpoons\;\mathrm{TP}$", ha="center", va="center", fontsize=7.5)
    A.text(LX + 0.52, 1.77, r"$k_{off}$", ha="center", va="center", fontsize=7, color=VERM)
    A.plot([LX, xs[4] + 0.04], [1.93, y0 + 0.28], color="#9CA3AF", lw=0.6, zorder=10)
    A.add_patch(Circle((xs[4] + 0.04, y0 + 0.28), 0.012, fc="#6B7280", ec="none", zorder=10))
    # boundary condition
    A.add_patch(FancyBboxPatch((0.45, 0.20), 3.40, 0.36, boxstyle="round,pad=0.0,rounding_size=0.06", fc="#F3F4F6", ec="#E5E7EB", lw=0.8))
    A.text(2.15, 0.38, r"$D\,\partial C/\partial x\,|_{x=0}\;=\;k_{on}\,C_s\,(\Gamma_{max}-\Gamma)\;-\;k_{off}\,\Gamma$", ha="center", va="center", fontsize=7.5)
    A.text(0.45, 0.08, r"baseline: $h$ = 9.9 mm, opening diameter 4.4 mm (150 µL)", ha="left", va="center", fontsize=6.3, color=MUTED)
    A.text(0.05, H - 0.06, "(a)", fontsize=9.5, fontweight="bold", va="top")

    # ------------------------------ panel (b): real concentration profiles
    Da, K, beta = 100, 1.0, 0.1; ts = np.array([1e-4, 1e-3, 1e-2])
    r = solve(Da, K, beta, t_end=1e-2, t_eval=ts, n=300)
    axb = fig.add_axes([4.50 / W, 0.62 / H, 1.62 / W, 1.98 / H])
    cols = [S.STAG[5], S.STAG[45], S.STAG[120]]; names = [r"$t^{*}=10^{-4}$", r"$t^{*}=10^{-3}$", r"$t^{*}=10^{-2}$"]
    for i, (t, c) in enumerate(zip(ts, cols)):
        u = r["u"][:, i]; us, _ = surface_state(u[0], r["theta"][i], r["xc"][0], Da, K)
        xx = np.concatenate([[0.0], r["xc"]]); uu = np.concatenate([[us], u])
        d = np.sqrt(np.pi * t)
        axb.axvline(d, color=c, lw=0.8, ls=(0, (1.5, 2.0)), alpha=0.9, zorder=2)
        lbl = rf"$\delta^{{*}}=\sqrt{{\pi t^{{*}}}}={d:.3f}$" if i == 2 else f"{d:.3f}"
        axb.text(d + 0.007, 0.03, lbl, fontsize=6.3, color=c, rotation=90, ha="left", va="bottom")
        axb.plot(xx, uu, color=c, lw=2.0, solid_capstyle="round", zorder=4)
        axb.plot([0], [us], "o", color=c, mec="white", mew=0.8, ms=5, zorder=6, clip_on=False)
        yk = 0.50 - 0.085 * i
        axb.plot([0.27, 0.33], [yk, yk], color=c, lw=2.0, solid_capstyle="round")
        axb.text(0.345, yk, names[i], fontsize=6.8, va="center", ha="left")
    axb.axhline(1.0, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=2)
    axb.text(0.495, 1.03, r"bulk, $u=1$", fontsize=6.5, color=MUTED, ha="right", va="bottom")
    axb.text(0.495, 0.03, "Da = 100, K = 1,\n" + r"$\beta=0.1$", fontsize=6.3, color=MUTED, ha="right", va="bottom", linespacing=1.15)
    axb.text(0.495, 0.80, "column top at\n" + r"$x^{*}=1$ (not shown)", fontsize=6.3, color=MUTED, ha="right", va="center", linespacing=1.15)
    axb.set_xlim(0, 0.5); axb.set_ylim(0, 1.12); axb.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5]); axb.set_yticks([0, 0.5, 1])
    axb.set_xlabel(r"$x^{*}=x/h$"); axb.set_ylabel(r"$u=C/C_0$")
    fig.text(4.02 / W, (H - 0.06) / H, "(b)", fontsize=9.5, fontweight="bold", va="top")
    S.save(fig, "Figure_1")


# ----------------------------------------------------------------------------------------------
def fig2():
    M = load("maps.json"); Da = np.array(M["Das"]); K = np.array(M["Ks"]); P = np.array(M["P"]["0.1"][0])
    lx, ly = np.log10(Da), np.log10(K)
    spl = RectBivariateSpline(ly, lx, np.log10(P), kx=3, ky=3)
    gx = np.linspace(lx[0], lx[-1], 600); gy = np.linspace(ly[0], ly[-1], 400)
    Z = np.maximum(spl(gy, gx), 0.0)
    X, Y = np.meshgrid(10 ** gx, 10 ** gy)
    W, H = 14 * CM, 8.2 * CM
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0.66 / W, 0.55 / H, 3.45 / W, 2.52 / H])
    cax = fig.add_axes([4.40 / W, 0.55 / H, 0.14 / W, 2.52 / H])
    vmax = np.log10(30)
    pm = ax.pcolormesh(X, Y, Z, cmap=S.DEEPSEA, vmin=0, vmax=vmax, shading="gouraud", rasterized=True, zorder=1)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-2, 1e3); ax.set_ylim(0.1, 31.623); ax.grid(False)
    cs = ax.contour(X, Y, 10 ** Z, levels=[1.5, 3, 10], colors="white", linewidths=1.2, zorder=3)
    cs.set_path_effects([pe.withStroke(linewidth=2.6, foreground="#0B0405", alpha=0.45)])
    txt = ax.clabel(cs, fmt={1.5: "1.5×", 3: "3×", 10: "10×"}, inline=True, fontsize=7.5, colors="white",
                    manual=[(15, 0.9), (70, 0.35), (420, 0.35)], inline_spacing=6)
    for t in txt: t.set_fontweight("bold"); t.set_path_effects([pe.withStroke(linewidth=2, foreground="#0B0405", alpha=0.6)])
    ax.axvline(0.3, color=INK, lw=0.9, ls=(0, (4, 3)), alpha=0.75, zorder=4)
    ax.text(0.26, 1.15, "Da < 0.3:\nstagnant ≈\nwell-mixed", fontsize=6.8, ha="right", va="center", color=INK, linespacing=1.15)
    ax.text(0.013, 19, "reaction-limited", fontsize=7.3, style="italic", ha="left", va="center", color=INK)
    ax.text(900, 0.135, "transport-limited", fontsize=7.3, style="italic", ha="right", va="center", color="white")
    ax.text(0.013, 0.135, r"$f=0.5,\ \beta=0.1$", fontsize=7, ha="left", va="center", color=INK)
    ax.set_xlabel(r"$Da=k_{on}\,\Gamma_{max}\,h\,/\,D$"); ax.set_ylabel(r"$K=k_{off}\,/\,(k_{on}\,C_0)$")
    cb = fig.colorbar(pm, cax=cax); cb.set_ticks([0, 0.25, 0.5, 0.75, 1.0, 1.25]); cb.outline.set_linewidth(0.5)
    cax.set_ylabel(r"$\log_{10}$ slowdown ratio", fontsize=8, labelpad=3); cax.yaxis.set_label_position("left")
    cax.tick_params(length=2.5, labelsize=7); cax.grid(False)
    cax.text(4.3, 1.50, "ratio", transform=cax.get_yaxis_transform(), fontsize=6.8, color=MUTED, va="bottom", ha="left")
    for rr in (1, 1.8, 3.2, 5.6, 10, 18):
        cax.text(4.3, np.log10(rr), f"{rr:g}×", transform=cax.get_yaxis_transform(), fontsize=6.8, color=MUTED, va="center", ha="left")
    S.save(fig, "Figure_2")


# ----------------------------------------------------------------------------------------------
def fig3():
    from matplotlib.lines import Line2D
    from stagnant_cell_model import theta_eq
    M = load("maps.json"); Da = np.array(M["Das"]); Ks = np.array(M["Ks"]); kidx = int(np.argmin(abs(Ks - 1.778))); K = Ks[kidx]
    W, H = 16 * CM, 6.6 * CM
    fig = plt.figure(figsize=(W, H))
    axa = fig.add_axes([0.62 / W, 0.55 / H, 2.55 / W, 1.78 / H]); axb = fig.add_axes([3.62 / W, 0.55 / H, 2.55 / W, 1.78 / H], sharey=axa)
    betas = [0.01, 0.1, 1.0, 10.0]; Y = {b: np.array(M["P"][str(b)][0][kidx]) for b in betas}
    f = 0.5; r4 = (np.pi / 4) * (f - K * np.log(1 - f)) ** 2 / ((1 + K) ** 2 * (-np.log(1 - f)))
    for b in betas:
        c = S.BETA[b]
        axa.plot(Da, Y[b], "-o", color=c, lw=1.8, ms=4, mec="white", mew=0.6, zorder=4)
        axb.plot(Da * b, Y[b], "-o", color=c, lw=1.8, ms=4, mec="white", mew=0.6, zorder=4)
        axa.text(Da[-1] * 1.35, Y[b][-1], rf"$\beta={b:g}$", color=c, fontsize=7.2, va="center", ha="left", fontweight="bold")
    for b in (0.1, 1.0):                         # Eq. (4) asymptote, drawn where theta_eq^2 (1+K) Da* >= 10
        te = theta_eq(K, b); pref = r4 * te ** 2 * (1 + K)
        xs = np.logspace(np.log10(10 / (te ** 2 * (1 + K))), np.log10(Da[-1] * b), 40)
        axb.plot(xs, pref * xs, color=S.BETA[b], lw=1.1, ls=(0, (1.5, 1.8)), zorder=3)
    xg = np.array([3.0, 40.0]); axb.plot(xg, 12 * xg / 3.0, color=MUTED, lw=0.9, ls=(0, (5, 3)), zorder=2)
    for ax in (axa, axb):
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.9, 450)
        ax.axhline(1.5, color=MUTED, lw=0.8, ls=(0, (1.5, 2.0)), zorder=1)
        ax.set_yticks([1, 10, 100]); ax.set_yticklabels(["1", "10", "100"])
    axa.set_xlim(7e-3, 5e3); axb.set_xlim(7e-5, 3e4)
    axa.text(8.5e-3, 1.58, "1.5×", fontsize=6.6, color=MUTED, ha="left", va="bottom")
    axa.text(9e-3, 330, r"$K=1.8$", fontsize=7.2, color=MUTED, ha="left", va="top")
    axa.set_xticks([1e-2, 1e-1, 1, 1e1, 1e2, 1e3]); axb.set_xticks([1e-4, 1e-2, 1, 1e2, 1e4])
    axa.set_xlabel(r"$Da$"); axb.set_xlabel(r"$Da\,\beta = k_{on}\,\Gamma_{max}^{2}/(C_0\,D)$")
    axa.set_ylabel(r"slowdown ratio, $f=0.5$"); plt.setp(axb.get_yticklabels(), visible=False); axb.tick_params(axis="y", length=0)
    axb.text(1.2e-4, 400, "$\\beta=0.01$ and $0.1$ coincide within 4%\n$\\beta=1$ lies 14 to 32% lower\n$\\beta=10$ is shifted about 10×".replace("\\\\", "\\"), fontsize=6.6, color=INK, ha="left", va="top", linespacing=1.3)
    axb.text(42, 12 * 42 / 3.0 * 1.0, "slope 1", fontsize=6.6, color=MUTED, ha="left", va="center")
    handles = [Line2D([0], [0], color=S.BETA[b], lw=1.8, marker="o", ms=3.5, mec="white", mew=0.5, label=rf"$\beta={b:g}$") for b in betas]
    handles += [Line2D([0], [0], color=MUTED, lw=1.1, ls=(0, (1.5, 1.8)), label=r"Eq. (4), $\beta\leq1$")]
    axb.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.60), fontsize=6.6, handlelength=2.2, labelspacing=0.35, borderaxespad=0)
    S.panel_letter(axa, "(a)", x=-0.17); S.panel_letter(axb, "(b)", x=-0.06)
    S.save(fig, "Figure_3")

# ----------------------------------------------------------------------------------------------
def _two_comp(Da, K, beta, t, delta):
    from scipy.integrate import solve_ivp
    def f(tt, y):
        th = y[0]; ub = 1 - beta * th
        d = min(np.sqrt(np.pi * max(tt, 1e-14)), 1.0) if delta == "grow" else delta
        us = (ub / d + Da * K * th) / (1 / d + Da * (1 - th))
        return [(ub - us) / d / beta]
    sol = solve_ivp(f, (0, t[-1]), [0.0], t_eval=t, method="LSODA", rtol=1e-9, atol=1e-12)
    if not sol.success or sol.y.shape[1] != len(t):
        raise RuntimeError(f"Two-compartment solver failed: {sol.message}")
    return sol.y[0]


def fig4():
    from scipy.optimize import minimize_scalar
    from matplotlib.lines import Line2D
    from stagnant_cell_model import solve, theta_eq
    W, H = 16 * CM, 12.3 * CM
    fig = plt.figure(figsize=(W, H))
    cases = [(100, 1.0, 0.1), (10, 1.0, 1.0), (1000, 1.0, 0.1), (30, 1.0, 0.01)]
    top_w, gap = 1.18, 0.32
    for k, (Da, K, beta) in enumerate(cases):
        ax = fig.add_axes([(0.52 + k * (top_w + gap)) / W, 3.14 / H, top_w / W, 1.50 / H])
        phi = Da / beta; te = theta_eq(K, beta); t_end = min(60 * max(1, 1 / (phi * (1 + K))), 5e4)
        t = np.logspace(-8, np.log10(t_end), 300); t[-1] = t_end
        pde = solve(Da, K, beta, t_end=t_end, t_eval=t, n=200, rtol=1e-8, atol=1e-11)["theta"]
        g1 = _two_comp(Da, K, beta, t, 1.0); gg = _two_comp(Da, K, beta, t, "grow")
        obj = lambda ld: np.max(np.abs(_two_comp(Da, K, beta, t, np.exp(ld)) - pde)) / te
        mres = minimize_scalar(obj, bounds=(np.log(0.003), np.log(3.0)), method="bounded", options={"xatol": 1e-3})
        db = float(np.exp(mres.x)); gb = _two_comp(Da, K, beta, t, db)
        e1, eb, eg = (100 * np.max(np.abs(c - pde)) / te for c in (g1, gb, gg))
        ax.semilogx(t, pde / te, color=S.PDE, lw=1.9, zorder=5)
        ax.semilogx(t, g1 / te, color=S.FIXED, lw=1.4, ls=(0, (4, 2.2)), zorder=4)
        ax.semilogx(t, gb / te, color=S.BESTFIXED, lw=1.7, ls=(0, (1.2, 1.6)), zorder=4)
        ax.semilogx(t, gg / te, color=S.GROW, lw=1.7, zorder=6)
        fr = pde / te; lo = t[np.argmax(fr > 0.01)] / 8; hi = t[np.argmax(fr > 0.999)] * 8 if np.any(fr > 0.999) else t[-1]
        ax.set_xlim(lo, hi); ax.set_ylim(0, 1.45); ax.set_yticks([0, 0.5, 1.0]); ax.set_yticklabels(["0", "0.5", "1"])
        ax.text(0.04, 0.985, rf"$Da={Da},\ K={K:g},\ \beta={beta:g}$", transform=ax.transAxes, fontsize=6.4, color=INK, ha="left", va="top")
        ax.text(0.04, 0.895, rf"best fixed $\delta^{{*}}={db:.2g}$", transform=ax.transAxes, fontsize=6.2, color=S.BESTFIXED, ha="left", va="top", fontweight="bold")
        ax.spines["left"].set_bounds(0, 1); ax.set_xlabel(r"$t^{*}$")
        if k == 0: ax.set_ylabel(r"$\theta/\theta_{eq}$")
        S.panel_letter(ax, "(" + "abcd"[k] + ")", x=-0.30 if k == 0 else -0.16)
    hs = [Line2D([0], [0], color=S.PDE, lw=1.9, label="PDE (full solution)"),
          Line2D([0], [0], color=S.FIXED, lw=1.4, ls=(0, (4, 2.2)), label=r"fixed layer, $\delta^{*}=1$"),
          Line2D([0], [0], color=S.BESTFIXED, lw=1.7, ls=(0, (1.2, 1.6)), label=r"best fixed layer"),
          Line2D([0], [0], color=S.GROW, lw=1.7, label=r"growing layer, $\delta^{*}=\sqrt{\pi t^{*}}$")]
    fig.legend(handles=hs, loc="center", bbox_to_anchor=(0.5, 2.58 / H), ncol=4, fontsize=7, handlelength=2.6, columnspacing=1.6, frameon=False)

    # strip plots from the 36 combinations
    C = np.array(load("closure.json"), dtype=float)          # beta, K, Da, err %, t50 ratio
    rows = [("all 36", None, INK), (r"$\beta=0.01$", 0.01, S.BETA[0.01]), (r"$\beta=0.1$", 0.1, S.BETA[0.1]), (r"$\beta=1$", 1.0, S.BETA[1.0]), (r"$\beta=10$", 10.0, S.BETA[10.0])]
    mk = {10: "o", 100: "s", 1000: "^"}
    axe = fig.add_axes([0.98 / W, 0.55 / H, 2.15 / W, 1.62 / H]); axf = fig.add_axes([3.90 / W, 0.55 / H, 2.15 / W, 1.62 / H], sharey=axe)
    rng = np.random.default_rng(1)
    for ax, col, xlim in ((axe, 3, (0, 25)), (axf, 4, (0.75, 2.05))):
        for ri, (lab, b, c) in enumerate(rows):
            y = 4 - ri; sel = C if b is None else C[C[:, 0] == b]; v = sel[:, col]
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.add_patch(Rectangle((q1, y - 0.27), q3 - q1, 0.54, fc=c, ec="none", alpha=0.13, zorder=2))
            ax.plot([med, med], [y - 0.27, y + 0.27], color=c, lw=1.9, zorder=3, solid_capstyle="butt")
            for da, mm in mk.items():
                s2 = sel[sel[:, 2] == da]
                ax.plot(s2[:, col], y + rng.uniform(-0.15, 0.15, len(s2)), mm, color=c, ms=3.6, mec="white", mew=0.4, ls="none", zorder=4, alpha=0.95)
        ax.set_xlim(*xlim); ax.set_ylim(-0.75, 6.45); ax.grid(axis="y", visible=False)
    axe.set_yticks([4, 3, 2, 1, 0]); axe.set_yticklabels([r[0] for r in rows]); axe.spines["left"].set_bounds(-0.5, 4.5)
    for tl, r in zip(axe.get_yticklabels(), rows): tl.set_color(r[2] if r[1] is not None else INK)
    plt.setp(axf.get_yticklabels(), visible=False); axf.tick_params(axis="y", length=0)
    e = C[:, 3]; q = C[:, 4]
    axe.axvline(10, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1); axe.text(10.5, -0.55, "10%", fontsize=6.5, color=MUTED, ha="left", va="center")
    axe.text(24.7, 6.12, f"median {np.median(e):.1f}, 90th percentile {np.percentile(e, 90):.0f}, worst {e.max():.0f}", fontsize=6.0, color=MUTED, ha="right", va="center")
    axf.axvline(1, color=MUTED, lw=0.9, zorder=1); axf.text(1.02, -0.55, "equal", fontsize=6.5, color=MUTED, ha="left", va="center")
    axf.text(2.03, 6.12, f"range {q.min():.2f} to {q.max():.2f}, median {np.median(q):.2f}", fontsize=6.0, color=MUTED, ha="right", va="center")
    kh = [Line2D([0], [0], marker=mm, color=MUTED, ls="none", ms=3.6, label=f"$Da={da}$") for da, mm in mk.items()]
    axe.legend(handles=kh, loc="upper left", bbox_to_anchor=(-0.02, 0.905), ncol=3, fontsize=6.0, handletextpad=0.1, columnspacing=0.8, borderaxespad=0)
    axe.set_xlabel(r"maximum error of the growing layer (% of $\theta_{eq}$)", fontsize=7.8); axf.set_xlabel(r"$t_{50}$, growing layer / PDE", fontsize=7.8)
    S.panel_letter(axe, "(e)", x=-0.36, y=1.01); S.panel_letter(axf, "(f)", x=-0.10, y=1.01)
    S.save(fig, "Figure_4")

# ----------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------
def fig5():
    from lod_model import BASE, theta_stagnant, theta_wellmixed
    L = load("lod.json"); times = L["times"]; stag, wm, eq, dep = L["stag"], L["wm"], L["eq"], L["dep"]
    W, H = 16 * CM, 7.0 * CM
    fig = plt.figure(figsize=(W, H))
    axa = fig.add_axes([0.62 / W, 0.62 / H, 2.62 / W, 1.90 / H]); axb = fig.add_axes([3.80 / W, 0.62 / H, 2.33 / W, 1.90 / H])
    p = dict(BASE); C = np.logspace(-2, 3, 26); thr = 3 * p["sigma"]
    for t in times:
        ts = np.array([theta_stagnant(c, p, t * 60.0) for c in C]); tw = np.array([theta_wellmixed(c, p, t * 60.0) for c in C])
        axa.plot(C, ts, color=S.STAG[t], lw=1.9, zorder=5, solid_capstyle="round")
        axa.plot(C, tw, color=S.WM[t], lw=1.7, ls=(0, (4, 2.2)), zorder=4)
    axa.set_xscale("log"); axa.set_yscale("log"); axa.set_xlim(1e-2, 1e3); axa.set_ylim(1e-4, 1.5)
    axa.axhspan(1e-4, thr, color="#F3F4F6", zorder=0); axa.axhline(thr, color=MUTED, lw=0.9, ls=(0, (1.5, 2.0)), zorder=2)
    axa.text(950, thr * 0.80, r"below $3\sigma_\theta$", fontsize=6.8, color=MUTED, ha="right", va="top")
    axa.text(950, thr * 1.18, r"$3\sigma_\theta=0.06$", fontsize=6.8, color=MUTED, ha="right", va="bottom")
    # 45 min detection limits
    ls, lw_ = stag[2], wm[2]
    for x, c, mk in ((ls, S.STAG[45], "o"), (lw_, S.WM[45], "s")):
        axa.plot([x, x], [1e-4, thr], color=c, lw=0.9, ls=(0, (3, 2)), zorder=3)
        axa.plot([x], [thr], mk, color=c, ms=5.2, mec="white", mew=0.8, zorder=8)
    ya = 1.4e-3
    axa.add_patch(FancyArrowPatch((lw_, ya), (ls, ya), arrowstyle="<|-|>", mutation_scale=6, lw=0.9, color=INK, zorder=7, shrinkA=0, shrinkB=0))
    axa.text(np.sqrt(ls * lw_), ya * 0.62, f"{ls / lw_:.1f}×", fontsize=7.2, color=INK, ha="center", va="top", fontweight="bold")
    axa.text(ls * 1.12, 1.6e-4, f"{ls:.1f} nM", fontsize=6.6, color=S.STAG[120], ha="left", va="bottom", fontweight="bold")
    axa.text(lw_ * 0.90, 1.6e-4, f"{lw_:.2f} nM", fontsize=6.6, color=S.WM[120], ha="right", va="bottom", fontweight="bold")
    # key: time shades (stagnant solid, well-mixed dashed) and cell types
    kx0, kx1 = 0.665, 0.760; ky = [0.30, 0.225, 0.15, 0.075]
    for k, t in enumerate(times):
        y = ky[k]
        axa.plot([kx0, kx0 + 0.075], [y, y], transform=axa.transAxes, color=S.STAG[t], lw=1.9, solid_capstyle="round", clip_on=False)
        axa.plot([kx1, kx1 + 0.075], [y, y], transform=axa.transAxes, color=S.WM[t], lw=1.7, ls=(0, (4, 2.2)), clip_on=False)
        axa.text(kx1 + 0.100, y, f"{t} min", transform=axa.transAxes, fontsize=6.4, color=INK, ha="left", va="center")
    from matplotlib.lines import Line2D
    axa.legend(handles=[Line2D([0], [0], color=S.STAG[45], lw=1.9, label="stagnant"), Line2D([0], [0], color=S.WM[45], lw=1.7, ls=(0, (4, 2.2)), label="well-mixed")],
               loc="upper left", bbox_to_anchor=(0.02, 0.985), fontsize=6.8, handlelength=2.4, labelspacing=0.3, borderaxespad=0)
    axa.set_xlabel(r"target concentration, $C_0$ (nM)"); axa.set_ylabel(r"$\theta(t_{meas})$")
    axa.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
    # panel (b)
    axb.plot(times, stag, "-o", color=S.STAG[45], lw=1.9, ms=4.6, mec="white", mew=0.7, zorder=5)
    axb.plot(times, wm, "--s", color=S.WM[45], lw=1.7, ms=4.4, mec="white", mew=0.7, zorder=5, dashes=(4, 2.2))
    for t, s_, w_ in zip(times, stag, wm):
        axb.plot([t, t], [w_, s_], color="#CBD5E1", lw=1.0, zorder=1)
        left = t > 5
        axb.text(t * (0.93 if left else 1.07), np.sqrt(s_ * w_), f"{s_ / w_:.1f}×", fontsize=6.6, color=MUTED, ha="right" if left else "left", va="center", fontweight="bold")
    axb.axhline(eq, color=INK, lw=1.0, ls=(0, (6, 2, 1.2, 2)), zorder=2); axb.axhline(dep, color=MUTED, lw=1.1, ls=(0, (1.2, 1.8)), zorder=2)
    axb.text(3.9, eq * 1.07, f"equilibrium limit, {eq:.2f} nM", fontsize=6.4, color=INK, ha="left", va="bottom")
    axb.text(3.9, dep * 0.93, f"sample-depletion bound, {dep:.2f} nM", fontsize=6.4, color=MUTED, ha="left", va="top")
    xg = np.array([12.0, 120.0]); axb.plot(xg, stag[2] * 1.9 * (xg / 45.0) ** -0.5, color=MUTED, lw=0.9, ls=(0, (5, 3)), zorder=2)
    axb.text(12.5, stag[2] * 1.9 * (12.5 / 45.0) ** -0.5 * 1.10, r"$\propto t^{-1/2}$ (Cottrell)", fontsize=6.4, color=MUTED, ha="left", va="bottom")
    axb.text(122, 21, r"stagnant: log-log slope $-0.6$", fontsize=6.4, color=S.STAG[120], ha="right", va="top")
    axb.set_xscale("log"); axb.set_yscale("log"); axb.set_xlim(3.6, 130); axb.set_ylim(0.1, 22)
    axb.set_xticks(times); axb.set_xticklabels([str(t) for t in times]); axb.xaxis.set_minor_locator(plt.NullLocator())
    axb.set_yticks([0.1, 1, 10]); axb.set_yticklabels(["0.1", "1", "10"])
    axb.set_xlabel(r"$t_{meas}$ (min)"); axb.set_ylabel("LOD (nM)")
    S.panel_letter(axa, "(a)", x=-0.17); S.panel_letter(axb, "(b)", x=-0.17)
    S.save(fig, "Figure_5")

# ----------------------------------------------------------------------------------------------
def fig6():
    from matplotlib.patches import Patch
    rows = load("tornado.json"); L = load("lod.json"); base = L["stag"][2]
    DISP = {"D (1.0 to 2.0e-6 cm2/s)": r"$D$ (1.0 to 2.0 $\times$ 10$^{-10}$ m$^2$ s$^{-1}$)", "k_on,0 (/3 to x3)": r"$k_{on,0}$ ($\div$3 to $\times$3)",
            "k_off (/3 to x3)": r"$k_{off}$ ($\div$3 to $\times$3)", "gamma (7.5 to 11.7e-12 mm2)": r"$\gamma$ (7.5 to 11.7 $\times$ 10$^{-18}$ m$^2$)",
            "S_P (1 to 5e10 mm-2)": r"$S_P$ (1 to 5 $\times$ 10$^{16}$ m$^{-2}$)", "yield (0.5 to 0.9)": r"$y$ (0.5 to 0.9)",
            "h (0.25 to 0.99 cm)": r"$h$ (2.5 to 9.9 mm)", "sigma (0.01 to 0.05)": r"$\sigma_\theta$ (0.01 to 0.05)"}
    rows = sorted(rows, key=lambda r: -abs(np.log(r[1][1] / r[1][0])))
    LOW, HIGH = "#2A9D8F", "#E9A23B"
    W, H = 12.5 * CM, 7.3 * CM
    fig = plt.figure(figsize=(W, H)); ax = fig.add_axes([2.12 / W, 0.55 / H, 2.62 / W, 1.78 / H])
    n = len(rows)
    for i, (name, (lo, hi)) in enumerate(rows):
        y = n - 1 - i
        if i % 2 == 0: ax.axhspan(y - 0.5, y + 0.5, color="#F7F8FA", zorder=0, lw=0)
        for val, col in ((lo, LOW), (hi, HIGH)):
            x = np.log10(val / base)
            if abs(x) < 1e-6: continue
            ax.barh(y, x, height=0.56, color=col, ec="white", lw=0.4, zorder=3)
            ax.text(x + (0.012 if x > 0 else -0.012), y, f"{val:.2f}", fontsize=6.4, color=INK, ha="left" if x > 0 else "right", va="center", zorder=5)
        if abs(np.log10(lo / base)) < 1e-6 and abs(np.log10(hi / base)) < 1e-6:
            ax.plot([0, 0], [y - 0.28, y + 0.28], color=INK, lw=1.4, zorder=4); ax.text(0.015, y, "no effect", fontsize=6.4, color=MUTED, ha="left", va="center")
    ax.axvline(0, color=INK, lw=0.9, zorder=4)
    ax.set_yticks(range(n)); ax.set_yticklabels([DISP[r[0]] for r in rows[::-1]], fontsize=7.2)
    ax.set_xlim(-0.40, 0.50); ax.set_ylim(-1.45, n - 0.45); ax.grid(axis="y", visible=False)
    ax.set_xticks([-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5]); ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    ax.set_xlabel(r"$\log_{10}(\mathrm{LOD}/\mathrm{LOD}_{base})$ at 45 min")
    sec = ax.secondary_xaxis("top", functions=(lambda x: base * 10 ** x, lambda v: np.log10(np.maximum(v / base, 1e-12))))
    sec.set_xticks([1, 2, base, 4, 6]); sec.set_xticklabels(["1", "2", f"{base:.1f}", "4", "6"])
    sec.get_xticklabels()[2].set_fontweight("bold"); sec.set_xlabel("LOD (nM)", fontsize=8.2, labelpad=3)
    sec.tick_params(length=3, labelsize=7.5); sec.spines["top"].set_linewidth(0.7)
    ax.legend(handles=[Patch(color=LOW, label="low end of range"), Patch(color=HIGH, label="high end of range")], loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=6.8, handlelength=1.2, columnspacing=1.6, borderaxespad=0.1)
    S.save(fig, "Figure_6")

# ----------------------------------------------------------------------------------------------
def fig7():
    from matplotlib.lines import Line2D
    D = load("sobol.json"); names = D["names"]
    DISP = {"D": r"$D$", "k_on,0": r"$k_{on,0}$", "k_off": r"$k_{off}$", "gamma": r"$\gamma$", "S_P": r"$S_P$", "y": r"$y$", "h": r"$h$", "sigma": r"$\sigma_\theta$"}
    S1, ST = np.array(D["S1"]), np.array(D["ST"]); S1c, STc = np.array(D["S1_ci"]), np.array(D["ST_ci"])
    order = np.argsort(-ST); n = len(order)
    CS1, CST = "#D08A00", "#0B3C5D"
    W, H = 12 * CM, 6.5 * CM
    fig = plt.figure(figsize=(W, H)); ax = fig.add_axes([0.70 / W, 0.58 / H, 2.55 / W, 1.78 / H])
    for r, idx in enumerate(order):
        y = n - 1 - r
        if r % 2 == 0: ax.axhspan(y - 0.5, y + 0.5, color="#F7F8FA", zorder=0, lw=0)
        ax.plot([S1c[0][idx], S1c[1][idx]], [y + 0.15, y + 0.15], color=CS1, lw=1.0, zorder=3, solid_capstyle="butt")
        ax.plot([STc[0][idx], STc[1][idx]], [y - 0.15, y - 0.15], color=CST, lw=1.5, zorder=3, solid_capstyle="butt")
        for xx in (STc[0][idx], STc[1][idx]): ax.plot([xx, xx], [y - 0.25, y - 0.05], color=CST, lw=1.0, zorder=3)
        ax.plot([S1[idx]], [y + 0.15], "o", mfc="white", mec=CS1, mew=1.3, ms=5.2, zorder=5)
        ax.plot([ST[idx]], [y - 0.15], "s", color=CST, mec="white", mew=0.5, ms=5.0, zorder=5)
        txt = f"{ST[idx]:.2f}  [{STc[0][idx]:.2f}, {STc[1][idx]:.2f}]"
        ax.text(1.17, y - 0.15, txt, fontsize=6.4, color=CST, ha="left", va="center", clip_on=False, fontweight="bold" if r < 2 else "normal")
    ax.text(1.17, n - 0.35, r"$S_T$ [95% interval]", fontsize=6.6, color=CST, ha="left", va="bottom", clip_on=False, fontweight="bold")
    ax.axvline(0, color=INK, lw=0.8, zorder=2)
    ax.set_yticks([n - 1 - r for r in range(n)]); ax.set_yticklabels([DISP[names[i]] for i in order], fontsize=8)
    for k, tl in enumerate(ax.get_yticklabels()):
        if k < 2: tl.set_fontweight("bold")
    ax.set_xlim(-0.2, 1.1); ax.set_ylim(-0.6, n - 0.4); ax.grid(axis="y", visible=False); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    ax.set_xticks([-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax.set_xlabel(r"Sobol index of $\log_{10}$ LOD at 45 min")
    ax.legend(handles=[Line2D([0], [0], marker="o", mfc="white", mec=CS1, mew=1.3, color=CS1, lw=1.0, ms=5, label=r"first order, $S_1$"),
                       Line2D([0], [0], marker="s", color=CST, lw=1.5, ms=5, label=r"total, $S_T$")],
              loc="upper right", bbox_to_anchor=(1.0, 0.585), fontsize=6.8, handlelength=2.2, labelspacing=0.35, borderaxespad=0)
    ax.text(1.09, 2.35, "64 base samples (640 model runs),\n2000 bootstrap resamples.\nFirst-order intervals are about\n" + r"$\pm0.5$: the split between $\sigma_\theta$" + "\nand $S_P$ is not resolved.", fontsize=6.0, color=MUTED, ha="right", va="top", linespacing=1.3)
    S.save(fig, "Figure_7")

# ----------------------------------------------------------------------------------------------
def fig8():
    from matplotlib.lines import Line2D
    N = load("noise.json"); SP = np.array([3e9, 1e10, 2e10, 5e10, 1e11]) * 1e6        # mm^-2 -> m^-2
    W, H = 16 * CM, 6.7 * CM
    fig = plt.figure(figsize=(W, H)); axs = []
    for k in range(3):
        axs.append(fig.add_axes([(0.58 + k * 1.90) / W, 0.95 / H, 1.78 / W, 1.50 / H], sharey=axs[0] if k else None))
    modes = [("relative", "relative noise"), ("absolute", "absolute noise"), ("mixed", "mixed noise")]; times = [5, 45, 120]
    for ax, (key, title), letter in zip(axs, modes, "abc"):
        for t in times:
            s_, w_ = N[key][str(t)]
            ax.plot(SP, w_, "--s", color=S.WM[t], lw=1.7, ms=3.6, mec="white", mew=0.5, zorder=4, dashes=(4, 2.2))
            ax.plot(SP, s_, "-o", color=S.STAG[t], lw=1.8, ms=3.8, mec="white", mew=0.5, zorder=5)
        ax.axvline(2e16, color=MUTED, lw=0.8, ls=(0, (1.5, 2.0)), zorder=1)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(2e15, 1.5e17); ax.set_ylim(0.1, 80)
        ax.set_xticks([1e16, 1e17]); ax.set_yticks([0.1, 1, 10]); ax.set_yticklabels(["0.1", "1", "10"])
        ax.text(0.04, 0.97, title, transform=ax.transAxes, fontsize=7.6, color=INK, ha="left", va="top", fontweight="bold")
        ax.set_xlabel(r"$S_P$ (m$^{-2}$)")
        S.panel_letter(ax, f"({letter})", x=-0.20 if letter == "a" else -0.08)
        if letter != "a": plt.setp(ax.get_yticklabels(), visible=False); ax.tick_params(axis="y", length=0)
    axs[0].set_ylabel("LOD (nM)")
    axs[0].text(2.2e16, 62, "baseline", fontsize=6.2, color=MUTED, ha="left", va="top")
    s45 = N["absolute"]["45"][0]
    axs[1].annotate(f"plateau {min(s45[2:]):.1f} to {max(s45[2:]):.1f} nM", xy=(6e16, 2.55), xytext=(6e16, 4.6), fontsize=6.2, color=S.STAG[120], ha="center", va="center", fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=S.STAG[120], lw=0.7, shrinkA=0, shrinkB=1))
    s, w = N["mixed"]["45"]; i_s, i_w = int(np.argmin(s)), int(np.argmin(w))
    axs[2].plot([SP[i_s]], [s[i_s]], "*", color=S.STAG[45], ms=10, mec="white", mew=0.7, zorder=8)
    axs[2].plot([SP[i_w]], [w[i_w]], "*", color=S.WM[45], ms=10, mec="white", mew=0.7, zorder=8)
    axs[2].text(1.3e17, 0.108, "baseline $S_P$: $\\sigma_\\theta=\\sqrt{2}\\times0.02$,\nstagnant LOD 3.4 nM at 45 min", fontsize=5.9, color=MUTED, ha="right", va="bottom", linespacing=1.25)
    hs = []
    for t in times:
        hs += [Line2D([0], [0], color=S.STAG[t], lw=1.8, marker="o", ms=3.6, mec="white", mew=0.4, label=f"stagnant, {t} min"),
               Line2D([0], [0], color=S.WM[t], lw=1.7, ls=(0, (4, 2.2)), marker="s", ms=3.4, mec="white", mew=0.4, label=f"well-mixed, {t} min")]
    hs += [Line2D([0], [0], marker="*", color=MUTED, ls="none", ms=8, mec="white", label="minimum at 45 min"), Line2D([0], [0], alpha=0, label="")]
    fig.legend(handles=hs, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=4, fontsize=6.6, handlelength=2.4, columnspacing=1.3, labelspacing=0.3, frameon=False)
    for model_name in ("absolute", "mixed"):
        if "45" not in N[model_name]:
            raise ValueError(f"Figure 8 data are missing the 45-min entry for {model_name}")
        series = N[model_name]["45"]
        if len(series) != 2 or any(not np.all(np.isfinite(np.asarray(v, dtype=float))) for v in series):
            raise ValueError(f"Figure 8 contains invalid 45-min data for {model_name}")
    S.save(fig, "Figure_8")

# ----------------------------------------------------------------------------------------------
def _optimum_curves():
    """LOD at 45 min against S_P for three sigma_a/sigma_r ratios (as in stage_optimum of reproduce_all.py); cached."""
    path = os.path.join(HERE, "results", "optimum_curves.json")
    if os.path.exists(path): return json.load(open(path))
    from lod_model import BASE, lod, theta_stagnant, theta_wellmixed
    SREL = 0.02; SPs = np.logspace(9.3, 11.0, 14); out = {"SP_mm2": SPs.tolist()}
    for rho in (0.25, 1.0, 4.0):
        ls, lw = [], []
        for sp in SPs:
            q = dict(BASE); q["SP"] = sp; q["sigma"] = float(np.sqrt((rho * SREL * 2e10 / sp) ** 2 + SREL ** 2))
            ls.append(lod(lambda c: theta_stagnant(c, q, 2700.0), q, lo=-4, hi=6)); lw.append(lod(lambda c: theta_wellmixed(c, q, 2700.0), q, lo=-4, hi=6))
        out[str(rho)] = {"stag": [None if not np.isfinite(v) else float(v) for v in ls], "wm": [None if not np.isfinite(v) else float(v) for v in lw]}
    json.dump(out, open(path, "w")); return out


def fig9():
    from matplotlib.lines import Line2D
    C = _optimum_curves(); O = load("optimum.json"); SP = np.array(C["SP_mm2"]) * 1e6
    W, H = 11 * CM, 7.6 * CM
    fig = plt.figure(figsize=(W, H)); ax = fig.add_axes([0.66 / W, 0.62 / H, 3.45 / W, 2.18 / H])
    rhos = [0.25, 1.0, 4.0]; ya = 0.215
    for rho in rhos:
        c = S.RATIO[rho]; ls = np.array([np.nan if v is None else v for v in C[str(rho)]["stag"]]); lw = np.array([np.nan if v is None else v for v in C[str(rho)]["wm"]])
        ax.plot(SP, lw, "--s", color=c, lw=1.6, ms=3.8, mec="white", mew=0.5, zorder=4, dashes=(4, 2.2), alpha=0.92)
        ax.plot(SP, ls, "-o", color=c, lw=1.9, ms=4.0, mec="white", mew=0.5, zorder=5)
        (xs_, ys_), (xw_, yw_) = O[str(rho)]
        for xo, arr in ((xs_ * 1e6, ls), (xw_ * 1e6, lw)):
            ok = np.isfinite(arr); yi = 10 ** np.interp(np.log10(xo), np.log10(SP[ok]), np.log10(arr[ok]))
            ax.plot([xo], [yi], "*", color=c, ms=10.5, mec="white", mew=0.8, zorder=8)
        ax.add_patch(FancyArrowPatch((xw_ * 1e6, ya), (xs_ * 1e6, ya), arrowstyle="-|>", mutation_scale=6, lw=1.0, color=c, zorder=7, shrinkA=0, shrinkB=0))
        ax.text(np.sqrt(xs_ * xw_) * 1e6, ya * 1.16, f"{xw_ / xs_:.1f}×", fontsize=6.8, color=c, ha="center", va="bottom", fontweight="bold")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1.5e15, 1.3e17); ax.set_ylim(0.15, 80)
    ax.set_yticks([0.3, 1, 3, 10, 30]); ax.set_yticklabels(["0.3", "1", "3", "10", "30"]); ax.set_xticks([1e16, 1e17])
    ax.set_xlabel(r"$S_P$ (m$^{-2}$)"); ax.set_ylabel("LOD at 45 min (nM)")
    hs = [Line2D([0], [0], color=S.RATIO[r], lw=1.9, label=rf"$\sigma_a/\sigma_r={r:g}$") for r in rhos]
    hs += [Line2D([0], [0], color=MUTED, lw=1.9, marker="o", ms=3.8, mec="white", mew=0.5, label="stagnant"),
           Line2D([0], [0], color=MUTED, lw=1.6, ls=(0, (4, 2.2)), marker="s", ms=3.6, mec="white", mew=0.5, label="well-mixed"),
           Line2D([0], [0], color=MUTED, marker="*", ls="none", ms=9, mec="white", label="optimum")]
    ax.legend(handles=hs, loc="upper right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=6.6, handlelength=2.3, columnspacing=1.3, labelspacing=0.35, borderaxespad=0.1)
    ax.text(1.7e15, 0.215 * 0.80, "", fontsize=6)
    S.save(fig, "Figure_9")

FIGS = {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5, "6": fig6, "7": fig7, "8": fig8, "9": fig9}

if __name__ == "__main__":
    S.setup()
    args = sys.argv[1:] or ["all"]
    for a in (list(FIGS) if args == ["all"] else args):
        FIGS[a](); print("built Figure", a, "->", S.OUT)
