"""Calibration curves, LOD and one-at-a-time sensitivity for a stagnant electrochemical DNA sensor.
Parameters are tied to published values (see PARAMS notes); sigma and h are assumptions to be varied."""
import numpy as np
from scipy.integrate import solve_ivp
from stagnant_cell_model import solve, theta_eq

NA = 6.022e23

# Baseline parameters. Sources / status:
#  D      : 150 um2/s for ~20-nt ssDNA (Sheehan & Whitman 2005, citing Stellwagen 2002)
#  kon0   : 9.34e4 M-1 s-1 at low probe density (Vanjur 2021, RPI optical, hydrogel, stirred)
#  gamma  : 9.6e-12 mm2 ; kon = kon0*exp(-2*gamma*S_P)  (Vanjur 2021)
#  koff   : 0.4e-4 s-1 (MCP4, Vanjur 2021), independent of S_P within their uncertainty
#  yield  : 0.75 (accessible fraction of probes, Vanjur 2021)
#  h      : ~0.99 cm from 150 uL over a 4.4 mm O-ring (thesis geometry, own calculation)
#  sigma  : ASSUMED noise (fraction of full-scale signal); LOD threshold = 3*sigma
BASE = dict(D=1.5e-6, kon0=9.34e4, gamma=9.6e-12, koff=0.4e-4, y=0.75,
            SP=2e10, h=0.99, sigma=0.02)

def derived(p):
    kon = p["kon0"] * np.exp(-2 * p["gamma"] * p["SP"])          # M-1 s-1
    Gm = p["y"] * p["SP"] * 100.0 / NA                           # mol/cm2 (mm-2 -> cm-2)
    return kon * 1e3, p["koff"], Gm                               # kon in cm3/mol/s

def groups(C0_nM, p):
    kon, koff, Gm = derived(p)
    C0 = C0_nM * 1e-12                                           # mol/cm3
    Da = kon * Gm * p["h"] / p["D"]
    K = koff / (kon * C0)
    beta = Gm / (C0 * p["h"])
    return Da, K, beta

def theta_stagnant(C0_nM, p, t_s):
    Da, K, beta = groups(C0_nM, p)
    ts = p["D"] * t_s / p["h"] ** 2
    r = solve(Da, K, beta, t_end=ts, t_eval=np.array([ts]), n=110, rtol=1e-6, atol=1e-10)
    return r["theta"][-1]

def theta_wellmixed(C0_nM, p, t_s):
    Da, K, beta = groups(C0_nM, p)
    phi = Da / beta
    ts = p["D"] * t_s / p["h"] ** 2
    s = solve_ivp(lambda t, y: phi * ((1 - beta * y[0]) * (1 - y[0]) - K * y[0]),
                  (0, ts), [0.0], t_eval=[ts], method="LSODA", rtol=1e-9, atol=1e-12)
    return s.y[0][-1]

def theta_equilibrium(C0_nM, p):
    Da, K, beta = groups(C0_nM, p)
    return theta_eq(K, beta)

def lod(fn, p, lo=-4.0, hi=5.0, it=18):
    """Smallest C0 (nM) with theta >= 3*sigma; bisection in log10 C0. nan if unreachable."""
    target = 3 * p["sigma"]
    if fn(10 ** hi) < target:
        return np.nan
    if fn(10 ** lo) >= target:
        return 10 ** lo
    a, b = lo, hi
    for _ in range(it):
        m = 0.5 * (a + b)
        if fn(10 ** m) >= target:
            b = m
        else:
            a = m
    return 10 ** (0.5 * (a + b))

def lod_depletion_bound(p):
    """Even perfect capture cannot exceed theta = C0*h/Gmax; C0_min = 3*sigma*Gmax/h (nM)."""
    _, _, Gm = derived(p)
    return 3 * p["sigma"] * Gm / p["h"] / 1e-12

if __name__ == "__main__":
    kon, koff, Gm = derived(BASE)
    print(f"baseline: kon={kon/1e3:.3g} M-1s-1, Kd={koff/kon*1e9*1e3:.3g} nM, Gmax={Gm*NA:.3g} cm-2")
    print("depletion bound on LOD (nM):", lod_depletion_bound(BASE))
