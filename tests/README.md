# Automated Verification & Physics Validation Suites (`tests/`)

This directory contains three independent test suites verifying topological watertightness, gas-dynamic ODE accuracy, Nonweiler & Newtonian/Prandtl–Meyer aerodynamics, canonical PEC electromagnetic targets, and wideband HRRP/RTI signal processing across both MATLAB and Node.js.

---

## Test Suites

### 1. [`test_waverider_suite.m`](test_waverider_suite.m) — MATLAB Unit & Manifold Suite (**12 / 12 PASS**)
Tests core package workflows:
- `WaveriderConfig` preset loading and validator rejection of invalid parameters
- `obliqueShock` and `taylorMaccoll` sanity checks
- 2-manifold watertightness (`NumBoundaryEdges == 0`, `NumNonManifoldEdges == 0`, `EulerCharacteristic == 2`, `EnclosedVolume > 0`) across all 4 forebody families + conformal booster frustum
- Binary/ASCII `.stl` and grouped `.obj` export round-trip file checks

```matlab
setup_waverider();
run('tests/test_waverider_suite.m');
```

### 2. [`test_physics_validation.m`](test_physics_validation.m) — MATLAB Physics & Canonical CEM Suite (**23 / 23 PASS**)
Rigorous independent physics audit:
- **Flow (`A1–A4`)**: 7-case comparison against independent SciPy RK45 reference table + NACA Report 1135 conical shock chart ($M_\infty = 2.0, \theta_c = 20^\circ \implies \beta = 37.80^\circ$).
- **Topology & Geometry (`B1–B5`)**: 42-configuration sweep (presets + 24 randomized + 13 validator-boundary extremes) checking directed half-edges, generalized winding numbers $W(\mathbf{q}) \in \{0, 1\}$, zero upper/lower surface crossover, zero degenerate faces, and Monte-Carlo vs. Divergence-Theorem enclosed volume.
- **Electromagnetics (`C1–C5`)**: $1\text{ m}\times 1\text{ m}$ PEC plate normal incidence & coarse 2-triangle off-normal $\text{sinc}^2$ lobes ($33\lambda$), $a=1\text{ m}$ PEC sphere optical limit $\pi a^2$, and $90^\circ$ dihedral corner reflector at symmetric ($\phi = 45^\circ$) and off-symmetric ($\phi = 30^\circ$) incidence.
- **Wideband HRRP (`D1–D6`)**: Two-scatterer range accuracy, $-3\text{ dB}$ mainlobe widths, and Harris sidelobe levels (`'rect'` $-13.26\text{ dB}$, `'hann'` $-31.5\text{ dB}$).
- **Aerodynamics & Symmetry (`E1–E3`)**: Nonweiler caret theorem $L/D = \cot\delta$, `computeAeroPolar(0) == AeroMetrics`, and bilateral symmetry ($|C_Y| = 0$, $\sigma(+\psi) = \sigma(-\psi)$).

```matlab
setup_waverider();
run('tests/test_physics_validation.m');
```

### 3. [`web/validate_physics.mjs`](web/validate_physics.mjs) — Node.js Web Studio Physics Suite (**26 / 26 PASS**)
Executes all 23 physics benchmarks directly against [`docs/physics.js`](../docs/physics.js) plus Web-specific input guardrails (`A5` sub-Mach shock rejection, `D4` unphysical $B > f_c$ bandwidth rejection) and `E4` real-time interactive latency benchmark ($< 250\text{ ms}$ full model + 5-band Az + 3-band El + 61-angle RTI budget; achieved $\approx 159\text{ ms}$).

```bash
node tests/web/validate_physics.mjs
```

See [`documentation/VALIDATION_REPORT.md`](../documentation/VALIDATION_REPORT.md) for the complete benchmark table.
