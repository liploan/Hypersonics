# Parametric Watertight Hypersonic Waveriding Glide Body Generator (MATLAB)

[![Live Interactive Web Studio](https://img.shields.io/badge/GitHub%20Pages-Live%20Interactive%20Waverider%20%26%20RCS%20Studio-0284c7?style=for-the-badge&logo=github)](https://liploan.github.io/Hypersonics/)

A modular MATLAB codebase (`+waverider`) and interactive client-side web studio (`docs/`) for generating, verifying, analyzing, visualizing, and exporting **2-manifold watertight 3D hypersonic waverider glide bodies**, **conformal aft rocket booster transition frustums**, and **UHF-to-X-band Radar Cross Section (RCS) & concave cavity multi-bounce signatures**.

* **Live Interactive Web Studio (GitHub Pages)**: [https://liploan.github.io/Hypersonics/](https://liploan.github.io/Hypersonics/)
* **GitHub Repository**: [https://github.com/liploan/Hypersonics](https://github.com/liploan/Hypersonics)

![Booster Integration Overview](docs/assets/waverider_booster_integration.png)
![UHF to X-Band RCS & Concave Underside Analysis](docs/assets/waverider_rcs_uhf_to_xband_caret.png)

---

## Key Capabilities

1. **Four Waverider Design Methodologies**:
   - **Axisymmetric Cone-Derived (`'conical'`)**: Solves the exact axisymmetric supersonic/hypersonic **Taylor-Maccoll ODE** (`+waverider/+flow/taylorMaccoll.m`) and integrates 3D shock-capturing streamlines from the conical shock wave downstream to the base plane.
   - **Sobieczky Osculating Cones (`'osculating_cone'`)**: Parameterizes the Inlet Capture Curve (ICC) with a configurable flat central shock segment and power-law curved outer wing segments, tracing local Taylor-Maccoll streamlines in osculating planes.
   - **Nonweiler Caret / 2D Oblique-Shock Wedge (`'caret'`)**: Uses exact $\theta\text{-}\beta\text{-}M$ oblique shock relations (`+waverider/+flow/obliqueShock.m`) to generate planar/anhedral compression ramps with attached oblique shock sheets.
   - **Kulfan CST Boost-Glide Waverider (`'cst'`)**: Uses Class-Shape Transformation (CST) Bernstein polynomials longitudinally and transversely (`+waverider/+geometry/generateCSTBody.m`) with a wave-trapping concave lower arch.

2. **Configurable Upper Surface, Conformal Aft Frustum & Watertight Closure**:
   - **`'freestream'`**: Zero-pressure-gradient freestream streamline tracing aft from the sharp shock-capturing leading edge.
   - **`'canopy'`**: Volumetric cambered dorsal canopy smoothly blended into the sharp leading edge and wingtips to maximize internal packaging volume and volumetric efficiency $\eta_V = V^{2/3} / S_{\text{plan}}$.
   - **Conformal Aft Transition Frustum & Conical Rocket Booster Adapter (`+waverider/+geometry/appendConformalFrustum.m`)**: Smoothly morphs the non-axisymmetric waverider aft base at $x = L$ into a circular rocket booster ring (radius `BoosterRadius` over axial length `FrustumLength`) using selectable blending laws (`'cubic'`, `'quintic'`, `'linear'`, `'conic'`), optionally followed by a conical/cylindrical rocket stage (`EnableRocketStage`).
   - **Guaranteed 2-Manifold Watertight Stitching (`+waverider/+geometry/stitchWatertightMesh.m`)**: Shares exact vertex indices along the leading edge, singular wingtip corners, lateral frustum/booster seams, and aft base cap, guaranteeing:
     - `size(freeBoundary(TR), 1) == 0` (zero boundary edges)
     - Every edge shared by **exactly 2 triangles** with opposite half-edge directions
     - Consistent outward-pointing normals and positive enclosed volume $V > 0$
     - Euler characteristic $\chi = V - E + F = 2$

3. **Mesh Export & Engineering Visualization**:
   - Binary & ASCII `.stl` export (`body.exportSTL(...)`) and grouped Wavefront `.obj` export (`body.exportOBJ(...)`).
   - **4-Panel Engineering Visualization & Attribute Dashboards (`+waverider/+viz`)**:
     - `body.plot()` (`plotWaverider.m`) — 3D watertight mesh patches, compression $C_p$ + attached shock sheet, aft base cross-section ($x = L$), and symmetry plane profile ($y = 0$).
     - `body.plotGeometryLines()` (`plotGeometryLines.m`) — 3D parametric rib/streamline lofting frame, multi-station cross-sections ($x/L = 0.25, 0.50, 0.75, 1.00$), planform thickness contour map, and longitudinal area-rule $A(x/L)$ & cumulative volume $V(x/L)$ curves.
     - `body.plotKeyAttributes()` (`plotKeyAttributes.m`) — Taylor-Maccoll shock-layer flowfield profiles ($C_p(\theta), M(\theta), |u_\theta(\theta)|$), angle-of-attack aerodynamic polar ($C_L, C_{D,\text{wave}}, L/D$ vs. $\alpha$), spanwise load distribution, and 3D face $C_p$ with outward unit normal quiver vectors.
     - `waverider.viz.plotFamilyComparison(bodies)` & `waverider.viz.plotTradeSpace(preset)` — Side-by-side comparison of all 4 waverider families and 2D contour maps of the $(\beta, h_{\text{canopy}})$ trade space & Pareto front.

---

## Repository Structure

```text
Hypersonics/
├── setup_waverider.m                 % Adds project paths
├── README.md                         % Documentation and quickstart
├── +waverider/
│   ├── WaveriderConfig.m             % Parametric configuration class & presets
│   ├── WaveriderBody.m               % Main generator, verifier, aero & export class
│   ├── +flow/
│   │   ├── obliqueShock.m            % Exact 2D oblique shock relations (theta-beta-M)
│   │   ├── taylorMaccoll.m           % Exact axisymmetric Taylor-Maccoll ODE solver
│   │   └── traceStreamline.m         % Vectorized RK4 conical streamline integrator
│   ├── +geometry/
│   │   ├── generateOsculatingCone.m  % Conical & Osculating Cones surface generator
│   │   ├── generateCaretWedge.m      % Nonweiler Caret / 2D wedge surface generator
│   │   ├── generateCSTBody.m         % Kulfan CST boost-glide surface generator
│   │   ├── buildUpperSurface.m       % Freestream & volumetric dorsal canopy builder
│   │   └── stitchWatertightMesh.m    % Manifold watertight mesh stitcher + base cap
│   ├── +mesh/
│   │   ├── verifyWatertight.m        % Topology, manifold, Euler & volume verification
│   │   ├── exportSTL.m               % Binary/ASCII STL writer
│   │   └── exportOBJ.m               % Grouped Wavefront OBJ writer
│   └── +viz/
│       ├── plotWaverider.m           % 4-panel 3D & cross-sectional visualizer
│       ├── plotGeometryLines.m       % 3D lofting ribs, station slices & area-rule plot
│       ├── plotKeyAttributes.m       % Flowfield, AoA polar, spanwise load & normals plot
│       ├── plotFamilyComparison.m    % Side-by-side 4-family comparison dashboard
│       └── plotTradeSpace.m          % 2D parametric contour maps & Pareto front
├── examples/
│   ├── demo_all_waveriders.m         % Generates, plots, and exports all 4 presets
│   ├── demo_parametric_sweep.m       % Parametric trade study (shock angle vs canopy)
│   ├── plot_geometry_blueprint.m     % Standalone script for geometry lofting & lines plan
│   ├── plot_key_attributes.m         % Standalone script for flowfield & aero attributes
│   └── plot_family_comparison.m      % Standalone script for family & trade-space plots
└── tests/
    └── test_waverider_suite.m        % Automated unit & manifold verification suite
```

---

## Quickstart Usage in MATLAB

```matlab
% 1. Initialize paths
setup_waverider();

% 2. Generate a preset Osculating Cones waverider (or 'conical', 'caret', 'cst')
body = waverider.WaveriderBody('osculating_cone');
body.printSummary();

% 3. Plot 4-panel 3D overview
body.plot();

% 4. Customize any parameter via WaveriderConfig
cfg = waverider.WaveriderConfig( ...
    'Method', 'osculating_cone', ...
    'Mach', 8.0, ...
    'Length', 6.0, ...
    'HalfSpan', 1.8, ...
    'ShockAngleDeg', 13.0, ...
    'PlanformPower', 0.60, ...
    'FlatShockFraction', 0.30, ...
    'UpperSurfaceType', 'canopy', ...
    'CanopyHeightRatio', 0.065, ...
    'NumStreamwise', 80, ...
    'NumHalfSpan', 50);

customBody = waverider.WaveriderBody(cfg);
customBody.printSummary();
customBody.exportSTL('output/custom_waverider.stl');
customBody.exportOBJ('output/custom_waverider.obj');
```
