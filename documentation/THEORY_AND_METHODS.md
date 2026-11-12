# Mathematical Theory, Numerical Methods & Modeling Scope

This manual documents the exact governing equations, coordinate systems, geometric lofting theorems, topological invariants, aerodynamic models, high-frequency electromagnetic scattering formulations, and wideband radar signal-processing algorithms implemented identically across `+waverider/` (MATLAB) and `docs/physics.js` (Web Studio).

---

## 1. Coordinate Frames & Angular Conventions

### 1.1 Vehicle Body Frame
All surface grids $(X, Y, Z)$ and watertight triangulations $(\mathbf{V}, \mathbf{F})$ are defined in a right-handed Cartesian body frame with origin at the nose apex $(0, 0, 0)$:
- **$+X$**: Axial aft direction ($X = 0$ at nose apex, $X = L$ at waverider trailing-edge base, $X = L + L_f$ at conformal booster adapter ring, $X = L + L_f + L_r$ at rocket stage base).
- **$+Y$**: Starboard spanwise direction ($Y \in [-b/2, +b/2]$ on the forebody; $Y = 0$ is the bilateral symmetry plane).
- **$+Z$**: Dorsal upward direction ($Z \le 0$ on shock-capturing compression surfaces below the apex plane).

### 1.2 Aerodynamic Freestream Vector
At angle of attack $\alpha$ (positive nose-up) and zero sideslip, the unit freestream velocity vector in body coordinates is:

$$
\hat{\mathbf{V}}_\infty(\alpha) = \begin{bmatrix} \cos\alpha \\ 0 \\ \sin\alpha \end{bmatrix}, \qquad \hat{\mathbf{L}}(\alpha) = \begin{bmatrix} -\sin\alpha \\ 0 \\ \cos\alpha \end{bmatrix}.
$$

### 1.3 Monostatic Radar Look-Direction Convention
For radar azimuth $\psi$ and look elevation $\theta_{\text{el}}$, the unit vector pointing **from the vehicle origin toward the monostatic radar** is:

$$
\hat{\mathbf{u}}_{\text{radar}}(\psi, \theta_{\text{el}}) = \begin{bmatrix} -\cos\theta_{\text{el}}\cos\psi \\ \cos\theta_{\text{el}}\sin\psi \\ \sin\theta_{\text{el}} \end{bmatrix}, \qquad \hat{\mathbf{k}}_i = -\hat{\mathbf{u}}_{\text{radar}}.
$$
- $\psi = 0^\circ, \theta_{\text{el}} = 0^\circ$: Nose-on incidence ($\hat{\mathbf{u}}_{\text{radar}} = [-1, 0, 0]^T$).
- $\psi = 90^\circ, \theta_{\text{el}} = 0^\circ$: Starboard broadside ($\hat{\mathbf{u}}_{\text{radar}} = [0, +1, 0]^T$).
- $\psi = 180^\circ, \theta_{\text{el}} = 0^\circ$: Tail-on incidence ($\hat{\mathbf{u}}_{\text{radar}} = [+1, 0, 0]^T$).
- $\theta_{\text{el}} < 0^\circ$: **Ventral** illumination (looking upward into the concave compression underside).
- $\theta_{\text{el}} > 0^\circ$: **Dorsal** illumination (looking downward onto the upper canopy).

---

## 2. Supersonic & Hypersonic Shock-Layer Solvers

### 2.1 Planar Oblique-Shock Relations (`+flow/obliqueShock.m`)
For calorically perfect gas ($\gamma = 1.4$ default), freestream Mach number $M_1 > 1$, and shock wave angle $\beta \in (\mu, 90^\circ)$ where $\mu = \arcsin(1/M_1)$ is the Mach wave angle, the normal Mach number is $M_{1n} = M_1 \sin\beta > 1$ and the flow deflection angle $\delta$ follows from the $\theta$-$\beta$-$M$ relation [1, 3]:

$$
\tan\delta = 2\cot\beta \,\frac{M_1^2 \sin^2\beta - 1}{M_1^2 (\gamma + \cos 2\beta) + 2}.
$$

Post-shock Mach number $M_2$, static pressure ratio $p_2/p_1$, density ratio $\rho_2/\rho_1$, pressure coefficient $C_{p,2}$, and nondimensional adiabatic limit speed $\bar{V}_2 = V_2 / V_{\max}$ are:

$$
M_{2n}^2 = \frac{1 + \frac{\gamma-1}{2}M_{1n}^2}{\gamma M_{1n}^2 - \frac{\gamma-1}{2}}, \qquad M_2 = \frac{M_{2n}}{\sin(\beta - \delta)},
$$

$$
\frac{p_2}{p_1} = 1 + \frac{2\gamma}{\gamma+1}(M_{1n}^2 - 1), \qquad C_{p,2} = \frac{p_2/p_1 - 1}{\frac{1}{2}\gamma M_1^2}, \qquad \bar{V}_2 = \left(1 + \frac{2}{(\gamma-1)M_2^2}\right)^{-1/2}.
$$

### 2.2 Axisymmetric Taylor–Maccoll Conical Flow (`+flow/taylorMaccoll.m`)
In spherical polar coordinates $(r, \theta)$ about a semi-infinite circular cone axis, conical self-similarity ($\partial/\partial r = 0$) reduces the steady Euler equations to the second-order **Taylor–Maccoll ODE** [2, 3] for $(u_r(\theta), u_\theta(\theta)) = (V_r, V_\theta)/V_{\max}$:

$$
\frac{d u_r}{d\theta} = u_\theta, \qquad \frac{d u_\theta}{d\theta} = \frac{\bar{a}^2(2 u_r + u_\theta \cot\theta) - u_r u_\theta^2}{u_\theta^2 - \bar{a}^2}, \qquad \bar{a}^2 = \frac{\gamma-1}{2}(1 - u_r^2 - u_\theta^2).
$$

Boundary conditions immediately behind the conical shock at $\theta = \beta$ are:

$$
u_r(\beta) = \bar{V}_2 \cos(\beta - \delta), \qquad u_\theta(\beta) = -\bar{V}_2 \sin(\beta - \delta) < 0.
$$

The system is integrated inward ($\theta$ decreasing from $\beta$) until the zero-normal-velocity cone boundary $u_\theta(\theta_c) = 0$:
- **MATLAB (`taylorMaccoll.m`)**: Dormand–Prince `ode45` [4, 5] with `RelTol = 1e-10`, `AbsTol = 1e-12`, and exact event detection `Events = @coneSurfaceEvent`.
- **Web Studio (`docs/physics.js`)**: Fixed-step RK4 ($\Delta\theta = -2\times 10^{-4}\text{ rad}$) with 4-iteration cubic Hermite Newton root refinement on the terminal bracket, matching SciPy/MATLAB $\theta_c$ to $4.3\times 10^{-6\circ}$.

Since the post-shock conical shock layer is homentropic, local Mach number $M(\theta)$ and pressure coefficient $C_p(\theta)$ follow isentropically from $M_2$ and $p_2/p_1$:

$$
M^2(\theta) = \frac{2}{\gamma-1}\frac{u_r^2 + u_\theta^2}{1 - u_r^2 - u_\theta^2}, \qquad \frac{p(\theta)}{p_1} = \frac{p_2}{p_1}\left(\frac{1 + \frac{\gamma-1}{2}M_2^2}{1 + \frac{\gamma-1}{2}M^2(\theta)}\right)^{\gamma/(\gamma-1)}.
$$

### 2.3 Conical Streamline Integration (`+flow/traceStreamline.m`)
In any meridian plane of an axisymmetric conical flowfield, cylindrical velocity components are $u_x = u_r\cos\theta - u_\theta\sin\theta$ and $u_{\text{cyl}} = u_r\sin\theta + u_\theta\cos\theta$. Each streamline originating on the shock cone at $(x_{\text{shk}}, r_{\text{shk}} = x_{\text{shk}}\tan\beta)$ satisfies:

$$
\frac{dr}{dx} = \frac{u_{\text{cyl}}(\theta)}{u_x(\theta)}, \qquad \theta(x, r) = \arctan\!\left(\frac{r}{x}\right),
$$

integrated downstream to $x_{\text{end}}$ via 4-substep classical RK4 with binary-search/pchip shock-layer interpolation.

---

## 3. Waverider Forebody Lofting & Non-Folding Geometry Guards

All four families share the nondimensional spanwise coordinate $\eta \in [-1, +1]$ and power-law leading-edge planform projection [9]:

$$
Y_{\text{LE}}(\eta) = \frac{b}{2}\,\eta, \qquad X_{\text{LE}}(\eta) = L\,|\eta|^{1/P}, \qquad P \in [0.25, 2.5].
$$

### 3.1 Axisymmetric Cone-Derived (`'conical'`, Jones et al. [7] & Rasmussen [8])
A single generating shock cone of half-angle $\beta$ is placed with virtual apex at $(-x_0, 0, z_v)$ where $z_v = x_0\tan\beta$ puts the vehicle nose $(0,0,0)$ on the top generator of the shock cone. To ensure real shock intersection $Y_{\text{LE}}^2 < R_{\text{shk}}^2$ across all $|\eta| \le 1$, $x_0$ is chosen to satisfy both $R_{\text{shk}}(X_{\text{LE}}) \ge 1.20\,|Y_{\text{LE}}|$ and target anhedral drop at the wingtip, so every leading-edge point lies **exactly** on the shock cone:

$$
Z_{\text{LE}}(\eta) = z_v - \sqrt{(X_{\text{LE}}(\eta) + x_0)^2\tan^2\beta - Y_{\text{LE}}(\eta)^2}, \qquad \phi(\eta) = \mathrm{atan2}\!\big(Y_{\text{LE}}(\eta),\, z_v - Z_{\text{LE}}(\eta)\big),
$$

and streamlines are integrated strictly along radial unit vector $\hat{\mathbf{e}}_r(\phi) = [0, \sin\phi, -\cos\phi]^T$ without artificial spanwise projection.

### 3.2 Inverse Osculating Cones (`'osculating_cone'`, Sobieczky et al. [10] & Ding et al. [12])
At the base plane $X = L$, the Inlet Capture Curve (ICC) shock trace $Z_s(Y_s)$ has a flat centerbody ($|Y_s| \le f_{\text{flat}} b/2$) blended via a $C^2$ regularization kernel $u_{\text{eff}} = \sqrt{u^2 + \varepsilon_u^2} - \varepsilon_u$ into a concave power-law outer wing ($m = 2.2$). With unit tangent $\mathbf{t}_s$ and inward normal $\mathbf{n}_s = [-Z_s', 1]^T / \sqrt{1 + (Z_s')^2}$:
- Local shock curvature is $\kappa = -Z_s'' / (1 + (Z_s')^2)^{3/2} \ge 0$, radius of curvature $R_c = 1/\kappa$, and osculating center $\mathbf{C}(Y_s) = \mathbf{r}_s(Y_s) + R_c\,\mathbf{n}_s(Y_s)$.
- **Evolute Caustic Guard**: By planar differential geometry, adjacent osculating planes $\mathbf{r}(s, d) = \mathbf{r}_s(s) + d\,\mathbf{n}_s(s)$ have Jacobian $(1 - \kappa d)$ and cross at the evolute $d = R_c$. Enforcing

  $$
  d_{\text{LE}}(Y_s) = \big(L - X_{\text{LE}}(Y_s)\big)\tan\beta \le 0.82\,R_c(Y_s)
  $$

  (clamping $X_{\text{LE}}$ aft only where extreme curvature/span combinations would violate $0.82\,R_c$) guarantees $1 - \kappa d \ge 0.18 > 0$ and strictly monotonic spanwise coordinates $Y_{\text{low}}(i, k+1) > Y_{\text{low}}(i, k)$ across the entire parameter space.

### 3.3 Nonweiler Caret (`'caret'`, Nonweiler [6]) & Kulfan CST (`'cst'`, Kulfan [11])
- **Caret**: Leading-edge height $Z_{\text{LE}}(\eta) = -X_{\text{LE}}(\eta)\tan\beta - |Y_{\text{LE}}(\eta)|\tan\Gamma_{\text{anh}}$, lower compression surface $Z_{\text{low}}(X, \eta) = Z_{\text{LE}}(\eta) - (X - X_{\text{LE}}(\eta))\tan\delta$, and uniform exact 2D oblique-shock pressure $C_p = C_{p,2}$. When $P = 1$ and $\text{UpperSurfaceType} = \text{'freestream'}$, each starboard and port wing panel is a single exact plane containing the attached planar shock along its leading edge.
- **CST Boost-Glide**: Longitudinal Bernstein shape polynomial $S_x(\xi) = \sum_{r=0}^n A_r \binom{n}{r}\xi^r(1-\xi)^{n-r}$ paired with class function $C_x(\xi) = \xi^{N_1}(1-\xi)^{N_2}$ and spanwise arch bell $(1-\eta^2)^2$, evaluated with modified Newtonian windward + Prandtl–Meyer leeward pressures.

### 3.4 Conformal Aft Booster Frustum (`+geometry/appendConformalFrustum.m`)
For axial parameter $\tau = (X - L)/L_f \in [0, 1]$, smooth blending function $\sigma(\tau) \in [0, 1]$ (`'linear'`, $C^1$ `'cubic'` $3\tau^2 - 2\tau^3$, $C^2$ `'quintic'` $10\tau^3 - 15\tau^4 + 6\tau^5$, or `'conic'`) interpolates upper and lower trailing-edge stations into semicircular arcs of radius $R_b$ centered at the exact area-weighted base centroid:

$$
z_c = \frac{\int_{-b/2}^{b/2} \frac{1}{2}\big(Z_{\text{up}}^2(L,Y) - Z_{\text{low}}^2(L,Y)\big)\,dY}{\int_{-b/2}^{b/2} \big(Z_{\text{up}}(L,Y) - Z_{\text{low}}(L,Y)\big)\,dY}.
$$

---

## 4. Mirror-Symmetric 2-Manifold Stitching & Topological Verification

1. **Mirror-Symmetric Stitcher (`+geometry/stitchWatertightMesh.m`)**: Collapses coincident leading-edge and wingtip vertices to unique shared indices and splits port ($c < N_{\text{half}}$) vs. starboard ($c \ge N_{\text{half}}$) surface quads along mirror-image diagonals so triangulation is strictly invariant under $Y \mapsto -Y$.
2. **Topological & Geometric Invariants (`+mesh/verifyWatertight.m`, [22–24])**:
   - Zero boundary edges ($e_{ij}$ count $= 2$ for every undirected edge) and zero conflicting half-edges ($\sum \mathrm{sgn}(e) = 0$).
   - Euler–Poincaré characteristic $\chi = V - E + F = 2$ (genus-0 topological sphere).
   - Enclosed volume via Divergence Theorem [24]: $V = \frac{1}{6}\sum_{f=1}^F \mathbf{v}_{f,1} \cdot (\mathbf{v}_{f,2} \times \mathbf{v}_{f,3}) > 0$.
   - **Generalized Winding Number** ([22, 23]): For any query point $\mathbf{q} \notin \partial\mathcal{M}$, summing exact triangle solid angles $\Omega_f(\mathbf{q})$ satisfies $W(\mathbf{q}) = \frac{1}{4\pi}\sum_f \Omega_f(\mathbf{q}) \in \{0, 1\}$, certifying absence of self-intersections or inverted pockets.

---

## 5. Calibrated Engineering Aerodynamics (`WaveriderBody.m`)

- **Design & Off-Design Inviscid Pressure**: At $\alpha = 0^\circ$, forebody compression facets integrate the exact flowfield pressure $C_{p,\text{design}}$ (Taylor–Maccoll or oblique-shock $C_{p,2}$). For general $\alpha$, each windward compression facet is calibrated via $k_f = \mathrm{clamp}(C_{p,\text{design},f} / (C_{p,\max}\sin^2\theta_{f,0}),\, 0.25,\, 3.0)$ so $C_{p,f}(\alpha) = k_f C_{p,\max}\sin^2\theta_f(\alpha)$ reproduces the exact shock-layer pressure at $\alpha = 0^\circ$ while varying smoothly with $\alpha$ [13]; leeward facets ($\sin\theta_f < 0$) use centered Prandtl–Meyer expansion $C_p = \frac{2}{\gamma M_\infty^2}\big[\max(1 - \frac{\gamma-1}{2}M_\infty \sin\theta_{\text{lee}}, 0)^{2\gamma/(\gamma-1)} - 1\big]$.
- **Viscous & Base Drag (`LD_viscous_est`)**: Uses Eckert's compressible reference temperature [14] at $30\text{ km}$ standard atmosphere ($T_\infty = 226.5\text{ K}, \rho_\infty = 0.01841\text{ kg/m}^3, T_w = 1000\text{ K}, r = 0.89$):

  $$
  \frac{T^*}{T_\infty} = 0.5\!\left(1 + \frac{T_w}{T_\infty}\right) + 0.22\,r\,\frac{\gamma-1}{2}M_\infty^2, \qquad C_{f,\text{turb}}^* = \frac{0.074}{(Re_L^*)^{1/5}}\frac{\rho^*}{\rho_\infty}\left(\frac{\mu^*}{\mu_\infty}\right)^{1/5},
  $$

  plus hypersonic base drag $C_{D,\text{base}} = \frac{1}{\gamma M_\infty^2}\frac{S_{\text{base}}}{S_{\text{plan}}}$ [3].

---

## 6. High-Frequency Electromagnetics: Gordon PO & GO-PO Double-Bounce

### 6.1 Exact Gordon Contour Physical Optics (`+rcs/facetPhaseIntegral.m` [15])
For wavenumber $k_0 = 2\pi/\lambda$ and monostatic phase vector $\mathbf{w} = 2 k_0 \hat{\mathbf{u}}_{\text{radar}}$, each illuminated PEC triangle $T_f$ ($\hat{\mathbf{n}}_f \cdot \hat{\mathbf{u}}_{\text{radar}} > 0$) contributes complex square-root RCS amplitude (`sigma = |S|^2`):

$$
S_f^{\text{PO}} = \frac{i k_0}{\sqrt{\pi}} (\hat{\mathbf{n}}_f \cdot \hat{\mathbf{u}}_{\text{radar}}) \iint_{T_f} e^{i \mathbf{w}\cdot\mathbf{r}}\,dA.
$$

Let $\mathbf{w}_t = \mathbf{w} - (\mathbf{w}\cdot\hat{\mathbf{n}}_f)\hat{\mathbf{n}}_f$ be the in-plane projection. When $|\mathbf{w}_t| L_{\max} \ge 10^{-2}$, Stokes' theorem reduces the 2D surface integral to an exact closed-form contour sum over the three directed CCW edges $\Delta\mathbf{r}_e = \mathbf{v}_{e,b} - \mathbf{v}_{e,a}$ with outward in-plane normal $\hat{\mathbf{m}}_e = (\Delta\mathbf{r}_e \times \hat{\mathbf{n}}_f)/|\Delta\mathbf{r}_e|$:

$$
\iint_{T_f} e^{i\mathbf{w}\cdot\mathbf{r}}\,dA = \frac{1}{i|\mathbf{w}_t|^2} \sum_{e=1}^3 (\mathbf{w}\cdot\Delta\mathbf{r}_e\times\hat{\mathbf{n}}_f)\, \mathrm{sinc}_{\text{rad}}\!\left(\frac{\mathbf{w}\cdot\Delta\mathbf{r}_e}{2}\right) e^{i\mathbf{w}\cdot(\mathbf{v}_{e,a}+\mathbf{v}_{e,b})/2}.
$$

In `docs/physics.js`, using $\mathrm{sinc}_{\text{rad}}((\phi_b-\phi_a)/2)e^{i(\phi_a+\phi_b)/2} = (e^{i\phi_b}-e^{i\phi_a})/[i(\phi_b-\phi_a)]$ evaluates $\cos\phi_v, \sin\phi_v$ **once per unique lit vertex**, cutting transcendental calls by $>8\times$ with zero loss of accuracy.

### 6.2 GO-PO Concave Ventral Double-Bounce (`+rcs/doubleBounceGOPO.m` [16, 21])
Lower compression quads are grouped into mirror-symmetric planar CCW super-cells. For each mutually visible concave pair $(P, Q)$:
1. Primary cell $P$ reflects incident wavevector $\hat{\mathbf{k}}_i = -\hat{\mathbf{u}}_{\text{radar}}$ into specular direction $\hat{\mathbf{k}}_r = \hat{\mathbf{k}}_i - 2(\hat{\mathbf{n}}_P\cdot\hat{\mathbf{k}}_i)\hat{\mathbf{n}}_P$.
2. If $\hat{\mathbf{n}}_Q \cdot \hat{\mathbf{k}}_r < 0$, the four vertices of $P$ are projected along $+\hat{\mathbf{k}}_r$ onto the plane of $Q$ and clipped against $Q$'s four edges via **Sutherland–Hodgman convex clipping** [21] to form exact illuminated polygon $\Omega_{PQ}$.
3. For any $\mathbf{q} \in \Omega_{PQ}$ traced from $\mathbf{p}(\mathbf{q}) = \mathbf{q} - t(\mathbf{q})\hat{\mathbf{k}}_r \in P$, the unfolded geometric path identity gives linear phase $\Phi(\mathbf{q}) = k_0(\hat{\mathbf{u}}_{\text{radar}} - \hat{\mathbf{k}}_r)\cdot\mathbf{q} - k_0 c_{\text{path}}$ with constant $c_{\text{path}} = (\hat{\mathbf{k}}_i - \hat{\mathbf{k}}_r)\cdot\mathbf{r}_P$, integrated analytically over $\Omega_{PQ}$ via Gordon's polygon contour formula:

   $$
   S_{P\to Q}^{\text{GO-PO}} = \frac{i k_0}{2\sqrt{\pi}} \Gamma\,\big[\hat{\mathbf{n}}_Q\cdot(\hat{\mathbf{u}}_{\text{radar}} - \hat{\mathbf{k}}_r)\big] e^{-i k_0 c_{\text{path}}} \iint_{\Omega_{PQ}} e^{i k_0(\hat{\mathbf{u}}_{\text{radar}} - \hat{\mathbf{k}}_r)\cdot\mathbf{q}}\,dA_Q, \qquad \Gamma = -1.
   $$

---

## 7. Wideband Stepped-Frequency HRRP & RTI (`+rcs/computeRangeProfile.m`)

For waveform center frequency $f_c$ and bandwidth $B \in (0, f_c]$, unwindowed slant-range resolution is $\Delta R = c_0 / (2B)$ [20]. Applying a two-term cosine spectral window $W(f) = a_0 + 2a_1\cos(2\pi(f - f_c)/B)$ over $f \in [f_c - B/2, f_c + B/2]$ produces the exact analytical Harris (1978) point-spread function [19] in normalized range $x = R/\Delta R$:

$$
h_{\text{PSF}}(R) = \frac{a_0\,\mathrm{sinc}(x) + a_1\big[\mathrm{sinc}(x - 1) + \mathrm{sinc}(x + 1)\big]}{a_0}, \qquad \mathrm{sinc}(x) = \frac{\sin(\pi x)}{\pi x},
$$
with $(a_0, a_1) = (1, 0)$ for Rectangular ($-3\text{ dB}$ width $0.886\,\Delta R$, sidelobe $-13.26\text{ dB}$), $(0.5, 0.25)$ for Hann ($1.44\,\Delta R$, $-31.5\text{ dB}$), and $(0.54, 0.23)$ for Hamming ($1.30\,\Delta R$, $-42.7\text{ dB}$). Double-bounce cavity returns appear at true delayed down-range $R_D = \frac{1}{2}\big[-(\hat{\mathbf{u}}_{\text{radar}} - \hat{\mathbf{k}}_r)\cdot\mathbf{r}_O + c_{\text{path}}\big] + x_{\text{ref}}\hat{u}_x$, shifting behind direct surface speculars by $+\frac{1}{2}|\mathbf{r}_O - \mathbf{p}_O|$.

---

## 8. Explicit Physical Scope & Modeling Limitations

1. **Inviscid Perfect-Gas Design Flow**: Taylor–Maccoll and oblique-shock solvers assume steady, attached shock waves in calorically perfect air ($\gamma = 1.4$). High-enthalpy vibrational excitation/dissociation ($M_\infty \gtrsim 10$), finite leading-edge radius entropy swallowing, and shock–boundary-layer interactions require Navier–Stokes / DSMC solvers.
2. **High-Frequency Asymptotic PEC CEM**: Gordon PO and two-bounce GO-PO are optical-regime approximations ($L/\lambda \gg 1$) on a PEC surface. At lower UHF frequencies where vehicle features approach $\lambda$, Rayleigh/resonance effects, creeping waves, traveling surface waves, and higher-order ($N \ge 3$) bounces or full multi-body shadow z-buffering require full-wave Method of Moments (MoM) / MLFMM solvers.
3. **Quasi-Narrowband Scatterer Amplitudes Across $B$**: `computeRangeProfile` evaluates complex scatterer amplitudes at $f_c$ and convolves with the exact spectral window PSF over $B \le f_c$; for ultra-wideband waveforms ($B / f_c \gtrsim 0.5$), frequency-dependent amplitude dispersion $\propto f$ across the band can slightly tilt point-spread sidelobes unless multi-frequency DFT synthesis is used.
