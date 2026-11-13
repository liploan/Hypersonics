# Independent Physics, Topology & Signature Validation Report

Every mathematical module in `+waverider/` (MATLAB) and `docs/physics.js` (Web Studio) is continuously verified against closed-form analytical theorems, historical supersonic wind-tunnel charts, and independent numerical solvers.

---

## 1. How to Run All Test Suites

```bash
# 1. Node.js Web Studio validation suite (36 checks, ~6 seconds)
node tests/web/validate_physics.mjs
```

```matlab
% 2. MATLAB unit, manifold & export suite (16 checks)
setup_waverider();
run('tests/test_waverider_suite.m');

% 3. MATLAB physics, generalized winding number, canonical CEM, polarimetry & RV suite (33 checks)
run('tests/test_physics_validation.m');
```

---

## 2. Quantitative Benchmark Results

### A. Supersonic & Hypersonic Shock-Layer Solvers (`A1–A5`)
| Check | Reference / Theorem | Tolerance | MATLAB Result | Web JS Result | Status |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **A1** Oblique shock $\delta, M_2$ | 7-point SciPy Rankine–Hugoniot table ($M \in [2, 10]$) | $< 10^{-4}$ | $4.9\times 10^{-6\circ}$ | $4.93\times 10^{-6\circ}$ | **PASS** |
| **A2** Taylor–Maccoll $\theta_c$ | Independent SciPy `solve_ivp` RK45 (`rtol=1e-12`) | $< 2\times 10^{-3\circ}$ | $1.9\times 10^{-6\circ}$ | $4.29\times 10^{-6\circ}$ | **PASS** |
| **A3** Cone surface $C_{p,c}$ | Independent SciPy `solve_ivp` conical pressure | $< 2\times 10^{-4}$ | $4.7\times 10^{-6}$ | $4.69\times 10^{-6}$ | **PASS** |
| **A4** NACA Report 1135 | $M_\infty = 2.0, \theta_c = 20.0^\circ \implies \beta = 37.80^\circ$ [1] | $\pm 0.05^\circ$ | $37.796^\circ$ | $37.796^\circ$ | **PASS** |
| **A5** Sub-Mach $\beta \le \mu$ guard | Rejects $M_\infty = 3.5, \beta = 15^\circ < \mu = 16.60^\circ$ | Throws error | Verified | `RangeError` | **PASS** |

---

### B. 42-Configuration Manifold, Winding-Number & Volume Sweep (`B1–B5`)
Evaluates **42 configurations**: all 5 presets (`conical`, `osculating_cone`, `caret`, `cst`, `booster_mounted`) + 24 deterministic randomized configurations across all families/blends + 13 extreme validator-boundary edge cases ($P = 0.25$, $P = 2.5$, $\Gamma_{\text{anh}} = \pm 35^\circ$, near-Mach $\beta$, extreme $R_b/L_f$, negative $\theta_r$, high canopy $h/L = 0.35$).

| Check | Invariant Tested | Tolerance | MATLAB Result | Web JS Result | Status |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **B1** Manifold closure | 0 boundary edges, 0 non-manifold edges, 0 conflicting half-edges, $\chi = 2$ | 0 failures | 42 / 42 clean | 42 / 42 clean | **PASS** |
| **B2** Self-intersection free | Generalized winding number $W(\mathbf{q}) \in \{0, 1\}$ [22, 23] | $< 0.05$ | $0.000$ | $0.000$ | **PASS** |
| **B3** No surface crossover | Upper/lower separation $Z_{\text{low}} - Z_{\text{up}} \le 0$ | $< 10^{-6}\text{ m}$ | $0.00\text{ m}$ | $0.00\text{ m}$ | **PASS** |
| **B4** Non-degenerate mesh | Strictly positive triangle areas across all patches | 0 degenerate | 0 faces | 0 faces | **PASS** |
| **B5** Enclosed volume | Monte-Carlo winding volume vs. Divergence Theorem [24] | $< 4.0\,\sigma$ | $< 2.1\,\sigma$ | $1.44\,\sigma$ | **PASS** |

---

### C. Canonical Electromagnetic Targets: PO & GO-PO (`C1–C5`)
| Check | Canonical Target & Exact Closed-Form Solution [17] | Tolerance | MATLAB Error | Web JS Error | Status |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **C1** Flat PEC plate ($\theta = 0^\circ$) | $1\text{ m}\times 1\text{ m}$ plate at $10\text{ GHz}$: $\sigma_0 = 4\pi A^2/\lambda^2 = +41.46\text{ dBsm}$ | $\pm 0.10\text{ dB}$ | $0.000\text{ dB}$ | $0.000\text{ dB}$ | **PASS** |
| **C2** Coarse plate $\text{sinc}^2$ | **Only 2 triangles** ($33\lambda$ wide) at $\theta \in \{0.5^\circ, 1^\circ, 2^\circ, 5^\circ, 10^\circ, 20^\circ\}$ | $\pm 0.50\text{ dB}$ | $0.00\text{ dB}$ max | $0.00\text{ dB}$ max | **PASS** |
| **C3** PEC sphere ($a = 1\text{ m}$) | Optical limit $\sigma = \pi a^2 = +4.9715\text{ dBsm}$ at $1\text{ GHz}$ & $3\text{ GHz}$ | $\pm 1.00\text{ dB}$ | $0.12 / 0.04\text{ dB}$ | $0.12 / 0.04\text{ dB}$ | **PASS** |
| **C4** $90^\circ$ dihedral ($\phi = 45^\circ$) | Symmetric corner reflector ($a=0.5\text{ m}, b=1\text{ m}, 10\text{ GHz}$): $+38.45\text{ dBsm}$ | $\pm 1.00\text{ dB}$ | $+0.03\text{ dB}$ | $+0.03\text{ dB}$ | **PASS** |
| **C5** $90^\circ$ dihedral ($\phi = 30^\circ$) | $15^\circ$ off-symmetric incidence ($\sin^2(2\phi) = 0.75$): $+35.44\text{ dBsm}$ | $\pm 1.00\text{ dB}$ | $+0.08\text{ dB}$ | $+0.08\text{ dB}$ | **PASS** |

---

### D. Wideband HRRP & Spectral Window Point-Spread Functions (`D1–D4`)
Two small point-like facets separated by $2.000\text{ m}$ along $X$ at $f_c = 10\text{ GHz}$, $B = 1000\text{ MHz}$ ($\Delta R = 14.99\text{ cm}$):

| Check | Theoretical Harris (1978) Benchmark [19] | Tolerance | MATLAB Result | Web JS Result | Status |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **D1** Peak positions (`rect` / `hann`) | Exact scatterer down-range coordinates | $< 1\text{ cm}$ | $0.04\text{ cm}$ | $0.03\text{--}0.04\text{ cm}$ | **PASS** |
| **D2** $-3\text{ dB}$ width (`rect` / `hann`) | $0.886\,\Delta R = 13.28\text{ cm}$ / $1.44\,\Delta R = 21.59\text{ cm}$ | $\pm 8\%$ | $13.41 / 21.74\text{ cm}$ | $13.41 / 21.74\text{ cm}$ | **PASS** |
| **D3** Peak sidelobe (`rect` / `hann`) | $-13.26\text{ dB}$ / $-31.5\text{ dB}$ | $\pm 2.0\text{ dB}$ | $-13.3 / -31.4\text{ dB}$ | $-13.3 / -31.4\text{ dB}$ | **PASS** |
| **D4** Unphysical $B > f_c$ guard | Rejects $B = 1000\text{ MHz}$ at $f_c = 0.5\text{ GHz}$ | Throws error | Verified | `RangeError` | **PASS** |

---

### E. Exact Aerodynamic Theorems, Bilateral Symmetry & Latency (`E1–E4`)
| Check | Exact Identity / Requirement | Tolerance | MATLAB Result | Web JS Result | Status |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **E1** Nonweiler caret theorem | $C_L = C_{p,2}$, $C_{D,w} = C_{p,2}\tan\delta$, $L/D = \cot\delta$ [6] | $< 10^{-11}$ | $1.2\times 10^{-14}$ | $3.82\times 10^{-14}$ | **PASS** |
| **E2** Off-design calibration | `computeAeroPolar(0) == AeroMetrics` across all presets | $< 10^{-13}$ | $0.00$ | $0.00$ | **PASS** |
| **E3** Bilateral symmetry | $\lvert C_Y\rvert = 0$ and $\sigma(+\psi) = \sigma(-\psi)$ | $< 10^{-11}\text{ dB}$ | $9.3\times 10^{-13}\text{ dB}$ | $1.76\times 10^{-12}\text{ dB}$ | **PASS** |
| **E4** Full Web Studio latency | Mesh + 5-band Az + 3-band El + 61-angle RTI (`booster_mounted`, all polarizations) | $< 250\text{ ms}$ | N/A | **$144.2\text{ ms}$** | **PASS** |

---

### F. Full Polarimetric Sinclair Matrix (`HH`, `VV`, `HV`, `VH`) & Circular Basis (`PP`, `OP`) (`F1–F5`)
| Check | Canonical Polarimetric Theorem / Invariant [17] | Tolerance | MATLAB Result | Web JS Result | Status |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **F1** PEC plate polarimetry | $1\text{ m}\times 1\text{ m}$ plate: $HH = VV = PP = \sigma_0$, $HV = VH = OP = 0$ ($-50\text{ dBsm}$ floor) | $< 0.05\text{ dB}$ | $3.11\times 10^{-9}\text{ dB}$ | $3.11\times 10^{-9}\text{ dB}$ | **PASS** |
| **F2** $0^\circ$ aligned $90^\circ$ dihedral | $S_{HH}^{\text{MB}} = -S_{VV}^{\text{MB}}$, $OP = +38.45\text{ dBsm}$, $\lvert S_{PP}^{\text{MB}}\rvert / \lvert S_{OP}^{\text{MB}}\rvert = 0$, $HV = 0$ | $< 0.2\text{ dB}$, $< 10^{-12}$ | $6.80\times 10^{-16}$, $OP - PP = 37.7\text{ dB}$ | $2.55\times 10^{-16}$, $OP - PP = 37.7\text{ dB}$ | **PASS** |
| **F3** $45^\circ$-canted $90^\circ$ dihedral | $100\%$ cross-pol conversion ($HV = VH = +38.45\text{ dBsm}$, $S_{HH}^{\text{MB}} = S_{VV}^{\text{MB}} = 0$) & roll-invariant $OP$ | $< 0.2\text{ dB}$, $< 10^{-10}\text{ dB}$ | $\lvert OP(45^\circ)-OP(0^\circ)\rvert = 2.13\times 10^{-14}\text{ dB}$ | $\lvert OP(45^\circ)-OP(0^\circ)\rvert = 1.42\times 10^{-14}\text{ dB}$ | **PASS** |
| **F4** 3D waverider polarimetry | Lorentz reciprocity $\lvert S_{HV}-S_{VH}\rvert = 0$, Parseval span conservation, & mirror handedness $OP(+\psi) = LL(-\psi)$ | $< 10^{-13}$, $< 10^{-9}\text{ dB}$ | $0.00$, span $2.70\times 10^{-16}$, mirror $6.54\times 10^{-13}\text{ dB}$ | $0.00$, span $2.89\times 10^{-16}$, mirror $4.40\times 10^{-11}\text{ dB}$ | **PASS** |
| **F5** Polarimetric HRRP / RTI | `OP` & `HV` reject $100\%$ of 1-bounce PO clutter ($-50\text{ dBsm}$) while isolating delayed 2-bounce cavity returns | $< 10^{-12}\text{ dB}$ | 1-bounce $-50.0\text{ dBsm}$, 2-bounce $-19.3\text{ dBsm}$ | 1-bounce $-50.0\text{ dBsm}$, 2-bounce $-9.4\text{ dBsm}$ | **PASS** |

---

### G. Conventional Axisymmetric Reentry Vehicles (`RV` Topology, Volumes, Aero & Concave Aft RCS) (`G1–G5`)
| Check | Analytical / Topological Benchmark | Tolerance | MATLAB Result | Web JS Result | Status |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **G1** 19-config RV manifold & cavity winding | $3\times 3\times 2$ (`mono`/`bi`/`tri` $\times$ `conical`/`spherical`/`ogive` $\times$ `flat`/`concave`) + booster RV: $\chi = 2$, $W_{\text{in}}=1$, $W_{\text{out}}=W_{\text{cav}}=0$ | 0 failures, $< 10^{-6}$ | 19 / 19 clean, $\max\lvert\Delta W\rvert = 1.80\times 10^{-14}$ | 19 / 19 clean, $\max\lvert\Delta W\rvert = 2.53\times 10^{-14}$ | **PASS** |
| **G2** Cone & paraboloid cavity volumes | Exact cone $V = \frac{1}{3}\pi R_b^2 L$ & recessed $p=2$ dome $\Delta V = \frac{1}{2}\pi R_b^2 d_{\text{aft}}$ | $< 0.2\%$, $< 1.0\%$ | Cone $0.008\%$, cap $0.156\%$ | Cone $0.008\%$, cap $0.156\%$ | **PASS** |
| **G3** Sharp-cone TM drag & RV axisymmetry | $C_{D,\text{wave}} S_{\text{plan}}/(\pi R_b^2) = C_{p,\text{cone}}$, $\lvert C_L(0^\circ)\rvert = \lvert C_Y(0^\circ)\rvert = \lvert C_m(0^\circ)\rvert = 0$, $C_L(+\alpha) = -C_L(-\alpha)$ | $< 0.2\%$, $< 10^{-13}$ | $0.008\%$, $\max = 8.64\times 10^{-17}$ | $0.008\%$, $\max = 1.33\times 10^{-16}$ | **PASS** |
| **G4** RV radar symmetry & circular disk tail-on | $\sigma(+\psi)=\sigma(-\psi)$, $\sigma(+\theta_{\text{el}})=\sigma(-\theta_{\text{el}})$, flat aft tail-on ($\psi=180^\circ$) $\sigma = 4\pi(\pi R_b^2)^2/\lambda^2$, $OP=0$ | $< 0.05\text{ dB}$ | $+40.34\text{ dBsm}$ ($0.001\text{ dB}$ err) | $+40.34\text{ dBsm}$ ($0.001\text{ dB}$ err) | **PASS** |
| **G5** Concave vs. flat RV aft plane (`OP` & HRRP) | Off-tail ($\psi=168^\circ$, $10\text{ GHz}$, $B=2000\text{ MHz}$): flat `OP` $= -50\text{ dBsm}$ vs. concave `OP` $> +10\text{ dBsm}$ & delayed cavity HRRP peak | $> +10\text{ dBsm}$, $\Delta R_{\text{cav}} > 20\text{ cm}$ | Flat $-50.0 \to$ Concave $+13.29\text{ dBsm}$, $+56.0\text{ cm}$ | Flat $-50.0 \to$ Concave $+13.29\text{ dBsm}$, $+55.7\text{ cm}$ | **PASS** |
