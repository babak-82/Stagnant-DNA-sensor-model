"""Edge (radial diffusion) enhancement of the diffusion-limited accumulation at a disk electrode of radius a
in a laterally unbounded sample, relative to planar (1D, Cottrell) accumulation.
Uses the Shoup-Szabo (1982) disk-current expression i = 4 n F D C a f(tau), tau = 4 D t / a^2,
f = 0.7854 + 0.8862 tau^-1/2 + 0.2146 exp(-0.7823 tau^-1/2).  A perfect sink is assumed (upper bound on the enhancement).
For a closed column with the electrode's own cross-section the 1D model is exact and the enhancement is 1."""
import numpy as np
from scipy.integrate import quad

# Units: D in cm^2/s; a in cm; t in s.
D = 1.5e-6  # cm^2/s
a = 0.22    # cm (electrode radius)

def f(tau):
    tau = np.asarray(tau, dtype=float)
    if np.any(tau <= 0):
        raise ValueError("tau must be strictly positive")
    return 0.7854 + 0.8862 / np.sqrt(tau) + 0.2146 * np.exp(-0.7823 / np.sqrt(tau))

def i_disk(t):
    if t <= 0:
        raise ValueError("t must be strictly positive")
    return 4 * D * a * f(4 * D * t / a ** 2)  # proportional to disk current per (n F C)

def i_planar(t):
    if t <= 0:
        raise ValueError("t must be strictly positive")
    return np.pi * a ** 2 * np.sqrt(D / (np.pi * t))  # proportional to planar current per (n F C)
if __name__ == "__main__":
    # checks of the expression itself
    print("limit checks: f(tau->inf) =", round(f(1e12), 4), "| ratio i_disk/i_planar at tau=1e-6:", round(i_disk(1e-6 * a**2 / (4*D)) / i_planar(1e-6 * a**2 / (4*D)), 4))
    print("t(min)  sqrt(Dt)/a  N_disk/N_planar (accumulation enhancement)  LOD correction factor")
    for tm in (5, 15, 45, 120, 480, 1440):
        t = tm * 60.0
        Nd = quad(i_disk, 0, t, limit=200)[0]; Np = quad(i_planar, 0, t, limit=200)[0]
        print(f"{tm:5d}   {np.sqrt(D*t)/a:6.2f}     {Nd/Np:6.3f}                                   {Np/Nd:6.3f}")
