# MATLAB Package API Reference (`+waverider`)

Initialize the package path from the repository root before calling any class or function:
```matlab
setup_waverider();
```

---

## 1. Core Classes

### 1.1 `waverider.WaveriderConfig`
Value-validated parametric configuration object defining flow conditions, forebody lofting, upper canopy, conformal booster frustum, and mesh resolution.

#### Constructor & Factory Methods
```matlab
cfg = waverider.WaveriderConfig('Name', Value, ...);
cfg = waverider.WaveriderConfig.fromPreset(presetName);
% presetName: 'conical' | 'osculating_cone' | 'caret' | 'cst' | 'booster_mounted'
cfg.validate(); % Checksbeta > asin(1/Mach) + 0.15 deg and geometric consistency
```

#### Key Properties
| Property | Default | Domain / Options | Description |
| :--- | :---: | :--- | :--- |
| `Method` | `'osculating_cone'` | `'conical'`, `'osculating_cone'`, `'caret'`, `'cst'` | Forebody shock/surface generation method |
| `Mach` | `6.0` | `(1.05, 30]` | Design freestream Mach number $M_\infty$ |
| `Gamma` | `1.4` | `(1.05, 1.67]` | Specific heat ratio $\gamma$ |
| `ShockAngleDeg` | `14.0` | `> asind(1/Mach) + 0.15` | Design shock wave angle $\beta$ [deg] |
| `Length` | `5.0` | `> 0` | Forebody axial length $L$ [m] |
| `HalfSpan` | `1.6` | `> 0` | Forebody half-span $b/2$ [m] |
| `PlanformPower` | `0.65` | `[0.25, 2.5]` | Leading-edge power exponent $P$ in $Y_{\text{LE}} \propto X_{\text{LE}}^P$ |
| `LeadingEdgeAnhedralDeg` | `4.0` | `[-35, 35]` | Ventral anhedral / concavity angle $\Gamma_{\text{anh}}$ [deg] |
| `FlatShockFraction` | `0.25` | `[0, 0.85]` | Flat centerbody shock fraction (`'osculating_cone'`) |
| `UpperSurfaceType` | `'canopy'` | `'freestream'`, `'canopy'` | Upper surface lofting mode |
| `CanopyHeightRatio` | `0.045` | `[0, 0.35]` | Peak dorsal canopy height $h/L$ |
| `EnableAftFrustum` | `false` | `logical` | Append conformal transition frustum aft of $X = L$ |
| `FrustumLength` | `1.60` | `> 0` | Axial length $L_f$ of conformal frustum [m] |
| `BoosterRadius` | `0.70` | `> 0` | Circular booster adapter radius $R_b$ [m] |
| `FrustumBlendLaw` | `'cubic'` | `'linear'`, `'cubic'`, `'quintic'`, `'conic'` | Axial cross-section morphing law $\sigma(\tau)$ |
| `FrustumZOffset` | `NaN` | `double` (`NaN` = auto centroid) | Booster axis vertical coordinate $z_c$ [m] |
| `EnableRocketStage` | `false` | `logical` | Append conical rocket stage aft of frustum |
| `RocketStageLength` | `2.20` | `> 0` | Rocket stage length $L_r$ [m] |
| `RocketConeAngleDeg` | `4.5` | `[-15, 25]` | Rocket half-angle $\theta_r$ [deg] (`0` = cylinder) |
| `NumStreamwise` | `65` | `>= 15` | Forebody streamwise stations $N_x$ |
| `NumHalfSpan` | `45` | `>= 11` | Half-span stations $N_{\text{half}}$ ($N_y = 2N_{\text{half}} - 1$) |

---

### 1.2 `waverider.WaveriderBody`
Constructs structured surfaces, stitches the 2-manifold watertight triangulation, verifies topological closure, and evaluates aerodynamic and electromagnetic signatures.

#### Constructor & Key Methods
```matlab
body = waverider.WaveriderBody(cfgOrPresetName);
body.printSummary();
polar = body.computeAeroPolar(alphaDegVec); % Off-design AoA polar struct
rcs   = body.computeRCS('FrequenciesGHz', [0.5 1.5 3 6 10], 'ElevationDeg', -22);
fig   = body.plot();                        % 4-panel 3D + Cp + base + symmetry view
fig   = body.plotGeometryLines();           % 3-view lofting ribs & area-rule curve
fig   = body.plotKeyAttributes();           % TM shock layer, AoA polar & face normals
fig   = body.plotRCS();                     % Multi-band UHF-to-X-band RCS dashboard
body.exportSTL('output/model.stl', 'binary');
body.exportOBJ('output/model.obj');
```

#### Key Read-Only Properties
- **Structured Grids**: `X_up, Y_up, Z_up`, `X_low, Y_low, Z_low`, `X_shk, Y_shk, Z_shk`, `Cp_low`.
- **Watertight Mesh**: `Vertices` ($N_V \times 3$), `Faces` ($N_F \times 3$), `FacePatchID` (`1`=Upper, `2`=Lower, `3`=Base, `4`=Frustum, `5`=Rocket), `FaceCpDesign` ($N_F \times 1$), `TR` (`triangulation`).
- **Metrics Structs**:
  - `MeshReport`: `IsWatertight`, `NumBoundaryEdges`, `NumNonManifoldEdges`, `NumConflictingHalfEdges`, `NumDegenerateFaces`, `HasSurfaceFold`, `EulerCharacteristic`, `EnclosedVolume`, `VolumeCentroid`, `PlanformArea`, `WettedArea`, `VolumetricEfficiency`.
  - `AeroMetrics`: `CL`, `CD_wave`, `CY`, `LD_inviscid`, `CD_friction_est`, `CD_base_est`, `LD_viscous_est`, `Cm_nose`.

---

## 2. Subpackage Functions

### 2.1 `+waverider/+flow`
- `shock = waverider.flow.obliqueShock(M1, betaDeg, gamma)` — Exact planar Rankine–Hugoniot oblique-shock properties (`deltaDeg`, `M2`, `p2_p1`, `rho2_rho1`, `T2_T1`, `Cp2`, `V2_bar`).
- `sol = waverider.flow.taylorMaccoll(M1, betaDeg, gamma)` — High-accuracy `ode45` conical-flow solver returning `thetaConeDeg`, `CpCone`, `MachCone`, radial arrays `thetaRad`, `ur`, `ut`, `Cp`, `Mach`, and `slopeRatio`.
- `[x_s, r_s, Cp_s, M_s] = waverider.flow.traceStreamline(sol, x_shk, r_shk, x_end, nPts)` — Vectorized 4-substep RK4 conical streamline integrator.

### 2.2 `+waverider/+geometry`
- `[X_low, Y_low, Z_low, X_shk, Y_shk, Z_shk, Cp_low, flowSol, meta] = waverider.geometry.generateOsculatingCone(cfg)` — Generates `'conical'` and `'osculating_cone'` forebodies (`meta.NumClampedStations` reports evolute caustic clamping).
- `[X_low, Y_low, Z_low, X_shk, Y_shk, Z_shk, Cp_low, flowSol] = waverider.geometry.generateCaretWedge(cfg)` — Generates exact Nonweiler `'caret'` compression surfaces.
- `[X_low, Y_low, Z_low, X_shk, Y_shk, Z_shk, Cp_low, flowSol] = waverider.geometry.generateCSTBody(cfg)` — Generates Kulfan `'cst'` boost-glide surfaces.
- `[X_up, Y_up, Z_up] = waverider.geometry.buildUpperSurface(X_low, Y_low, Z_low, cfg)` — Builds freestream or dorsal-canopy upper surfaces sharing exact leading-edge and wingtip boundaries.
- `[X_up, Y_up, Z_up, X_low, Y_low, Z_low, Cp_low, info] = waverider.geometry.appendConformalFrustum(...)` — Appends conformal transition frustum and optional conical rocket stage.
- `[V, F, patchID, faceCpDesign] = waverider.geometry.stitchWatertightMesh(X_up, Y_up, Z_up, X_low, Y_low, Z_low, numBaseLayers, splitIndices, Cp_low)` — Mirror-symmetric 2-manifold triangulation stitcher.

### 2.3 `+waverider/+mesh`
- `report = waverider.mesh.verifyWatertight(V, F, patchID, surfStruct)` — Full directed half-edge, Euler, degenerate-triangle, structured-fold, and Divergence-Theorem volume audit.
- `waverider.mesh.exportSTL(V, F, filename, mode)` — Binary (`'binary'`) or ASCII (`'ascii'`) STL exporter.
- `waverider.mesh.exportOBJ(V, F, patchID, filename)` — Patch-grouped Wavefront `.obj` exporter.

### 2.4 `+waverider/+rcs`
- `I = waverider.rcs.facetPhaseIntegral(V, F, w)` — Exact vectorized Gordon (1975) triangle contour phase integral $\iint_{T_f} e^{i\mathbf{w}\cdot\mathbf{r}}\,dA$ for an $M \times 3$ wavevector array `w`.
- `pairs = waverider.rcs.buildConcavePairList(body)` — Builds mirror-symmetric planar CCW super-cells on the lower compression surface and enumerates mutually visible concave pairs.
- `[E_mb, contrib] = waverider.rcs.doubleBounceGOPO(pairs, u_hat, k0)` — Exact Sutherland–Hodgman clipped GO-PO two-bounce cavity solver returning total complex scattering amplitude `E_mb` and per-patch `[amp, L_half, r_O]` rows.
- `rcs = waverider.rcs.computeRCS(body, 'FrequenciesGHz', ..., 'AzimuthDeg', ..., 'ElevationDeg', ..., 'IncludeMultiBounce', true, 'IncludeEdgeDiffraction', true, 'Compute2DSphere', true)` — Full multi-band monostatic RCS sweep.
- `rti = waverider.rcs.computeRangeProfile(body, 'CenterFreqGHz', 10.0, 'BandwidthMHz', 1000, 'ElevationDeg', -22, 'AzimuthDeg', 0:2:180, 'Window', 'hann')` — Wideband stepped-frequency / matched-filter RTI sinogram and 1D HRRP solver (enforces $B \le f_c$).

### 2.5 `+waverider/+viz`
- `plotWaverider`, `plotGeometryLines`, `plotKeyAttributes`, `plotFamilyComparison`, `plotTradeSpace`, `plotBoosterIntegration`, `plotRCSAnalysis`, `plotRCSConcavityComparison`, `plotRangeAngleIntensity`.
