| claim | value |
|---|---|
| V1 max |dtheta| vs well-mixed ODE | 6.13e-05 |
| V2 Cottrell relative error (min, max) | 1.52e-04, 2.86e-03 |
| V3 mass conservation max |M-1| | 8.9e-16 |
| V4 |theta_end - theta_eq| | 7.9e-11 |
| V5 Laplace benchmark max relative error | 5.4e-05 |
| grid refinement 100->600 cells, max change of slowdown ratio f=0.5 (%) | 0.07 |
| grid refinement 100->600 cells, max change of slowdown ratio f=0.9 (%) | 0.07 |
| baseline LOD(45 min), n=110 / 300 / 600 cells (nM) | 2.4090 / 2.4081 / 2.4081 |
| slowdown f=0.5, beta=0.1, K=1.78: Da grid | [0.01, 0.032, 0.1, 0.316, 1.0, 3.162, 10.0, 31.623, 100.0, 316.228, 1000.0] |
| slowdown f=0.5, beta=0.1, K=1.78: ratio | [1.0, 1.01, 1.02, 1.06, 1.11, 1.23, 1.46, 2.01, 3.5, 8.01, 22.2] |
| slowdown f=0.9, beta=0.1, K=1.78: ratio | [1.0, 1.0, 1.0, 1.1, 1.2, 1.6, 2.6, 6.0, 16.8, 51.2, 160.0] |
| collapse beta=0.01 vs 0.1 at equal Da*beta, K=0.10 (min,max ratio) | 1.010, 1.029 |
| collapse beta=0.01 vs 0.1 at equal Da*beta, K=0.26 (min,max ratio) | 1.010, 1.033 |
| collapse beta=0.01 vs 0.1 at equal Da*beta, K=0.68 (min,max ratio) | 1.009, 1.039 |
| collapse beta=0.01 vs 0.1 at equal Da*beta, K=1.78 (min,max ratio) | 1.008, 1.034 |
| collapse beta=0.01 vs 0.1 at equal Da*beta, K=4.64 (min,max ratio) | 1.005, 1.019 |
| collapse beta=0.01 vs 0.1 at equal Da*beta, K=12.12 (min,max ratio) | 1.003, 1.007 |
| collapse beta=0.01 vs 0.1 at equal Da*beta, K=31.62 (min,max ratio) | 1.001, 1.002 |
| beta=1 / beta=0.1 at equal Da*beta, K~1.8 (range) | 0.68 to 0.86 |
| scaling slope (points with E>=30) | 0.99 |
| penalty/E at Da>=1e4, K=0.3 | 0.38 |
| penalty/E at Da>=1e4, K=1.0 | 0.52 |
| penalty/E at Da>=1e4, K=3.0 | 0.68 |
| perfect-sink estimate (pi/4) f^2 / (-ln(1-f)) at f=0.5 | 0.283 |
| scaling2 obs/pred range, beta=0.01 | 1.06 to 2.07 |
| scaling2 obs/pred range, beta=0.1 | 1.07 to 1.98 |
| scaling2 obs/pred range, beta=1.0 | 1.16 to 1.53 |
| scaling2 prefactor obs range over K,f at beta=0.01 | 0.23 to 1.57 |
| Eq. (4) slowdown ~ r/(1+K), f=0.5: maximum at K | 0.03 (lower grid edge: decreasing in K) |
| Eq. (4) slowdown ~ r/(1+K), f=0.9: maximum at K | 0.83 (interior maximum) |
| PDE slowdown at Da=1000, beta=0.1, f=0.5: K of the maximum on the map grid | 0.26 |
| PDE slowdown at Da=1000, beta=0.1, f=0.9: K of the maximum on the map grid | 1.78 |
| closure delta=sqrt(pi t): cases | 36 |
| closure max error % of theta_eq: median, 90th pct, max | 9.5, 21.2, 22.4 |
| closure t50 ratio closure/PDE: min, median, max | 0.80, 1.14, 1.98 |
| closure worst four cases (beta,K,Da,err%) | [[10.0, 0.1, 10, 22.4], [10.0, 0.1, 100, 22.1], [10.0, 0.1, 1000, 21.65], [10.0, 1.0, 10, 21.58]] |
| growing-layer error range over six cases (%) | 3.8 to 11.0 |
| fixed delta*=1 error range, four cases (%) | 36 to 77 |
| best fixed delta* error range / delta* range | 5 to 19 %, delta* 0.017 to 0.232 |
| baseline k_on (M-1 s-1), Kd (nM), Gamma_max (cm-2), h (cm) | 6.36e+04, 0.629, 1.5e+12, 0.99 |
| LOD stagnant 5/15/45/120 min (nM) | [9.103, 4.549, 2.409, 1.415] |
| LOD well-mixed (nM) | [3.339, 1.18, 0.465, 0.253] |
| LOD ratio stagnant/well-mixed | [2.73, 3.86, 5.18, 5.6] |
| LOD equilibrium, depletion bound (nM) | 0.191, 0.151 |
| log-log slope of stagnant LOD, 5 to 120 min | -0.59 |
| Cottrell LOD estimate at 45 min (nM) | 2.08 |
| LOD from K_d alone (nM), 3 sigma threshold | 0.040 |
| baseline LOD 45 min stagnant, well-mixed, ratio | 2.409, 0.465, 5.18 |
| OAT D (1.0 to 2.0e-6 cm2/s) | 2.875 -> 2.132 nM |
| OAT k_on,0 (/3 to x3) | 3.105 -> 2.189 nM |
| OAT k_off (/3 to x3) | 2.390 -> 2.468 nM |
| OAT gamma (7.5 to 11.7e-12 mm2) | 2.382 -> 2.439 nM |
| OAT S_P (1 to 5e10 mm-2) | 1.315 -> 5.782 nM |
| OAT yield (0.5 to 0.9) | 1.719 -> 2.824 nM |
| OAT h (0.25 to 0.99 cm) | 2.409 -> 2.409 nM |
| OAT sigma (0.01 to 0.05) | 1.201 -> 6.077 nM |
| stagnant/well-mixed LOD ratio at 45 min over sigma = 0.01-0.05 and S_P = 1e16-5e16 m^-2: min, max | 3.58, 6.55 |
| ratio at the baseline S_P = 2e16 m^-2 over sigma = 0.01, 0.02, 0.05 | 5.24, 5.18, 5.00 |
| solutal Rayleigh number of the depletion layer, C0 = 2.4 nM, 45 min | 1.8 |
| solutal Rayleigh number of the depletion layer, C0 = 10.0 nM, 45 min | 7.5 |
| thermal Rayleigh number of the 9.9 mm column for a vertical temperature difference of 0.01 K | 161 |
| thermal Rayleigh number of the 9.9 mm column for a vertical temperature difference of 0.1 K | 1612 |
| thermal Rayleigh number of the 9.9 mm column for a vertical temperature difference of 1.0 K | 16115 |
| critical Rayleigh numbers for comparison (textbook values, rigid walls) | about 1.7e3 (fluid heated from below) |
| Sobol N, evaluations, variance of log10 LOD | 64, 640, 0.0751 |
| Sobol D: S1 [95% CI], ST [95% CI] | 0.01 [-0.08, 0.09], 0.02 [0.01, 0.03] |
| Sobol k_on,0: S1 [95% CI], ST [95% CI] | -0.02 [-0.14, 0.10], 0.03 [0.02, 0.05] |
| Sobol k_off: S1 [95% CI], ST [95% CI] | 0.00 [-0.01, 0.01], 0.00 [0.00, 0.00] |
| Sobol gamma: S1 [95% CI], ST [95% CI] | 0.00 [-0.01, 0.01], 0.00 [0.00, 0.00] |
| Sobol S_P: S1 [95% CI], ST [95% CI] | 0.53 [0.04, 1.03], 0.52 [0.38, 0.71] |
| Sobol y: S1 [95% CI], ST [95% CI] | 0.04 [-0.10, 0.17], 0.05 [0.03, 0.07] |
| Sobol h: S1 [95% CI], ST [95% CI] | 0.00 [0.00, 0.00], 0.00 [0.00, 0.00] |
| Sobol sigma: S1 [95% CI], ST [95% CI] | 0.50 [0.02, 1.02], 0.55 [0.40, 0.77] |
| Sobol sum of first-order indices | 1.07 |
| Cottrell estimate: share of var(ln LOD) from sigma / S_P / y / D (%) | 45.9 / 45.9 / 6.1 / 2.1 |
| noise relative 5 min stagnant / well-mixed LOD (nM) | [3.181, 5.553, 9.103, 20.601, 44.446] / [2.365, 2.731, 3.339, 5.996, 15.541] |
| noise relative 45 min stagnant / well-mixed LOD (nM) | [0.564, 1.315, 2.409, 5.782, 11.927] / [0.287, 0.355, 0.465, 0.893, 2.191] |
| noise relative 120 min stagnant / well-mixed LOD (nM) | [0.296, 0.754, 1.415, 3.436, 7.025] / [0.126, 0.175, 0.253, 0.532, 1.214] |
| noise absolute 5 min stagnant / well-mixed LOD (nM) | [25.145, 11.299, 9.103, 8.194, 8.809] / [19.533, 5.64, 3.339, 2.355, 3.033] |
| noise absolute 45 min stagnant / well-mixed LOD (nM) | [4.258, 2.654, 2.409, 2.307, 2.376] / [2.374, 0.733, 0.465, 0.352, 0.429] |
| noise absolute 120 min stagnant / well-mixed LOD (nM) | [2.194, 1.518, 1.415, 1.372, 1.401] / [1.046, 0.359, 0.253, 0.21, 0.239] |
| noise mixed 5 min stagnant / well-mixed LOD (nM) | [25.496, 12.687, 12.94, 22.204, 45.337] / [19.818, 6.356, 4.784, 6.474, 15.858] |
| noise mixed 45 min stagnant / well-mixed LOD (nM) | [4.315, 2.974, 3.415, 6.23, 12.164] / [2.409, 0.825, 0.665, 0.963, 2.236] |
| noise mixed 120 min stagnant / well-mixed LOD (nM) | [2.223, 1.7, 2.005, 3.702, 7.165] / [1.062, 0.405, 0.361, 0.574, 1.239] |
| mixed noise: sigma_theta at the baseline S_P of 2e16 m^-2 (sqrt(2) x 0.02) | 0.0283 |
| baseline S_P, stagnant 45 min LOD: sigma = 0.02 (Table 1) / mixed noise (nM) | 2.409 / 3.415 |
| optimum S_P (m-2) and min LOD (nM) at 45 min, sigma_a/sigma_r = 0.25: stagnant | well-mixed | 4e+15, 1.09 | 8.9e+15, 0.40 |
| optimum S_P (m-2) and min LOD (nM) at 45 min, sigma_a/sigma_r = 1.0: stagnant | well-mixed | 9.8e+15, 2.98 | 2e+16, 0.67 |
| optimum S_P (m-2) and min LOD (nM) at 45 min, sigma_a/sigma_r = 4.0: stagnant | well-mixed | 2.5e+16, 10.08 | 4e+16, 1.69 |
| edge enhancement N_disk/N_planar at 5 min (laterally unbounded sample) | 1.086 |
| edge enhancement N_disk/N_planar at 15 min (laterally unbounded sample) | 1.150 |
| edge enhancement N_disk/N_planar at 45 min (laterally unbounded sample) | 1.266 |
| edge enhancement N_disk/N_planar at 120 min (laterally unbounded sample) | 1.450 |
| edge enhancement N_disk/N_planar at 480 min (laterally unbounded sample) | 1.952 |
| edge enhancement N_disk/N_planar at 1440 min (laterally unbounded sample) | 2.712 |
| Pereiro Box 3 static case: h (mm), Da, K, beta, Da* | 2.92, 1344, 0.077, 0.79, 1057 |
| Pereiro Box 3: Da * delta/h with delta = sqrt(2 D t_R) (their quoted Da = 44) | 44.3 |
| Pereiro Box 3: t50 stagnant (h), well-mixed (s), ratio | 3.26, 52, 227 |
| Pereiro Box 3: observed prefactor ratio/E vs Eq. (4) r(K, 0.5) | 0.295 vs 0.299 |
| opening area (cm2), h (cm), h^2/D (days) | 0.152, 0.99, 7.5 |
| thesis check: fixed parameters | K_d = 0.5 nM, C0 = 2 uM, D = 150 um2/s, k_off = K_d k_on, h from 150 uL over the 4.4 mm opening |
| thesis check, planar electrode: time to 90% of equilibrium coverage; (k_on in 1/(M s), Gamma_max per projected area in m^-2, Da*) -> seconds | (9.34e+04, 1e+16, 0.0859) -> 14.4; (9.34e+04, 5e+16, 2.15) -> 26.6; (9.34e+03, 1e+16, 0.00859) -> 129.6; (9.34e+03, 5e+16, 0.215) -> 158.1 |
| thesis check, roughness factor 19.2 (nanoporous / planar active area), pore diffusion ignored: time to 90% of equilibrium coverage; (k_on in 1/(M s), Gamma_max per projected area in m^-2, Da*) -> seconds | (9.34e+04, 1e+16, 31.6) -> 134.0; (9.34e+04, 5e+16, 790) -> 2725.1; (9.34e+03, 1e+16, 3.16) -> 310.7; (9.34e+03, 5e+16, 79) -> 2974.1 |
| 1D validity ratio sqrt(Dt)/a at 5 min | 0.10 |
| 1D validity ratio sqrt(Dt)/a at 15 min | 0.17 |
| 1D validity ratio sqrt(Dt)/a at 45 min | 0.29 |
| 1D validity ratio sqrt(Dt)/a at 120 min | 0.47 |
| 1D validity ratio sqrt(Dt)/a at 480 min | 0.94 |
| 1D validity ratio sqrt(Dt)/a at 1440 min | 1.64 |
