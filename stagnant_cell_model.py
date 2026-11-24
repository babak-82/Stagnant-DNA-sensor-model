"""
Stagnant (no-flow) electrochemical cell: transport-limited reversible binding.

Physical picture
----------------
A liquid column of height h sits on a functionalised electrode (x = 0).
The analyte (bulk concentration C0 at t = 0) diffuses (coefficient D) toward
the electrode and binds reversibly to immobilised receptors (Langmuir):

    dGamma/dt = k_on * C_s * (Gamma_max - Gamma) - k_off * Gamma

The top boundary (x = h) is impermeable, so the sample is finite and can be
depleted by capture.

Dimensionless variables
-----------------------
    x* = x/h,  t* = D t / h^2,  u = C/C0,  theta = Gamma/Gamma_max

Dimensionless groups
--------------------
    Da   = k_on * Gamma_max * h / D     (surface reaction vs diffusion)
    K    = k_off / (k_on * C0) = Kd/C0  (relative affinity)
    beta = Gamma_max / (C0 * h)         (receptor capacity vs sample amount)

Equations
---------
    u_t = u_xx                              0 < x* < 1
    u_x(1) = 0
    u_x(0) = beta * theta_t = Da * [u_s (1 - theta) - K theta]
    beta * theta_t = Da * [u_s (1 - theta) - K theta]

Conservation: integral(u dx) + beta * theta = 1.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix


def make_grid(n=300, stretch=10.0):
    """Cell faces clustered toward the electrode (x = 0)."""
    i = np.arange(n + 1)
    xf = (np.exp(stretch * i / n) - 1.0) / (np.exp(stretch) - 1.0)
    xc = 0.5 * (xf[:-1] + xf[1:])
    dx = np.diff(xf)
    return xf, xc, dx


def surface_state(u0, theta, x0, Da, K):
    """Surface concentration from flux continuity between cell 0 and the surface."""
    us = (u0 / x0 + Da * K * theta) / (1.0 / x0 + Da * (1.0 - theta))
    flux = (u0 - us) / x0  # flux into the surface (dimensionless, D = 1)
    return us, flux


def solve(Da, K, beta, t_end=10.0, n=300, stretch=10.0, t_eval=None,
          rtol=1e-8, atol=1e-11):
    xf, xc, dx = make_grid(n, stretch)
    dxc = np.diff(xc)  # distances between neighbouring cell centres

    def rhs(t, y):
        u = y[:n]
        theta = y[n]
        us, flux = surface_state(u[0], theta, xc[0], Da, K)
        du = np.zeros(n)
        f = (u[1:] - u[:-1]) / dxc  # diffusive flux between cells (+ means toward +x)
        du[:-1] += f / dx[:-1]
        du[1:] -= f / dx[1:]
        du[0] -= flux / dx[0]
        dtheta = flux / beta
        return np.concatenate([du, [dtheta]])

    sp = lil_matrix((n + 1, n + 1))
    for i in range(n):
        for j in (i - 1, i, i + 1):
            if 0 <= j < n:
                sp[i, j] = 1
    sp[0, n] = 1
    sp[n, 0] = 1
    sp[n, n] = 1

    y0 = np.concatenate([np.ones(n), [0.0]])
    if t_eval is None:
        t_eval = np.logspace(-8, np.log10(t_end), 200)
        t_eval[-1] = t_end
    sol = solve_ivp(rhs, (0.0, t_end), y0, method="BDF", t_eval=t_eval,
                    jac_sparsity=sp, rtol=rtol, atol=atol)
    u = sol.y[:n]
    theta = sol.y[n]
    mass = (u * dx[:, None]).sum(axis=0) + beta * theta
    return dict(t=sol.t, theta=theta, u=u, xc=xc, dx=dx, mass=mass)


def theta_eq(K, beta):
    """Equilibrium coverage with a finite sample: (1 - beta*theta)(1 - theta) = K*theta."""
    a = beta
    b = -(1.0 + beta + K)
    c = 1.0
    disc = b * b - 4 * a * c
    return (-b - np.sqrt(disc)) / (2 * a)


def wellmixed(K, beta, phi, t):
    """Reaction-limited reference (uniform bulk): dtheta/dt = phi*[(1-beta*theta)(1-theta) - K*theta]."""
    sol = solve_ivp(lambda tt, th: phi * ((1 - beta * th[0]) * (1 - th[0]) - K * th[0]),
                    (0, t[-1]), [0.0], t_eval=t, method="LSODA", rtol=1e-9, atol=1e-12)
    return sol.y[0]


def t_fraction(t, theta, target):
    """First time theta reaches `target` (linear interpolation); nan if never."""
    idx = np.where(theta >= target)[0]
    if len(idx) == 0:
        return np.nan
    j = idx[0]
    if j == 0:
        return t[0]
    return t[j - 1] + (target - theta[j - 1]) * (t[j] - t[j - 1]) / (theta[j] - theta[j - 1])


# ----------------------------------------------------------------------------
# Validation against analytical limits
# ----------------------------------------------------------------------------
def validate():
    print("=== Validation ===")

    # V1: reaction-limited limit (Da -> 0) should match well-mixed ODE
    Da, K, beta = 1e-3, 0.5, 0.1
    phi = Da / beta
    t = np.logspace(-1, 3, 60)
    r = solve(Da, K, beta, t_end=1e3, t_eval=t)
    ref = wellmixed(K, beta, phi, t)
    err = np.max(np.abs(r["theta"] - ref))
    print(f"V1 reaction-limited vs well-mixed ODE: max |dtheta| = {err:.2e}")

    # V2: Cottrell limit (irreversible, perfect sink, semi-infinite): theta = (2/beta)*sqrt(t/pi)
    Da, K, beta = 1e5, 0.0, 1.0
    t = np.array([1e-5, 1e-4, 1e-3, 1e-2])
    r = solve(Da, K, beta, t_end=1e-2, t_eval=t, n=400, stretch=12.0)
    ana = (2.0 / beta) * np.sqrt(t / np.pi)
    rel = np.abs(r["theta"] - ana) / ana
    for ti, th, an, re in zip(t, r["theta"], ana, rel):
        print(f"V2 Cottrell  t*={ti:.0e}: numeric={th:.5f} analytic={an:.5f} rel.err={re:.2e}")

    # V3: mass conservation
    Da, K, beta = 5.0, 0.3, 0.5
    r = solve(Da, K, beta, t_end=5.0)
    print(f"V3 mass conservation: max |M-1| = {np.max(np.abs(r['mass'] - 1)):.2e}")

    # V4: equilibrium coverage
    r = solve(Da, K, beta, t_end=50.0)
    print(f"V4 equilibrium: numeric={r['theta'][-1]:.6f} analytic={theta_eq(K, beta):.6f}")


if __name__ == "__main__":
    validate()
