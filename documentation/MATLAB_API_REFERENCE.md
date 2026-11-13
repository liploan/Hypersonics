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
% Waverider presets: 'conical' | 'osculating_cone' | 'caret' | 'cst' | 'booster_mounted'
% Reentry Vehicle (RV) presets: 'rv_monoconic' | 'rv_biconic' | 'rv_triconic_flare' | 'rv_sharp_conic'
cfg.validate(); % Checks beta > asind(1/Mach) + 0.15 deg (waveriders) or conic/stage bounds (RVs)
```

#### Key Properties
| Property | Default | Domain / Options | Description |
| :--- | :---: | :--- | :--- |
| `Method` | `'osculating_cone'` | `'conical'`, `'osculating_cone'`, `'caret'`, `'cst'`, `'rv'` | Forebody shock/surface or axisymmetric RV generation method |
| `Mach` | `6.0` | `(1.05, 30]` | Design freestream Mach number $M_\infty$ |
| `Gamma` | `1.4` | `(1.05, 1.67]` | Specific heat ratio $\gamma$ |
| `ShockAngleDeg` | `14.0` | `> asind(1/Mach) + 0.15` | Design shock wave angle $\beta$ [deg] (waverider methods) |
| `Length` | `5.0` | `> 0` | Forebody axial length $L$ [m] |
| `HalfSpan` | `1.6` | `> 0` | Forebody half-span $b/2$ [m] (waverider methods) |
| `PlanformPower` | `0.65` | `[0.25, 2.5]` | Leading-edge power exponent $P$ in $Y_{\text{LE}} \propto X_{\text{LE}}^P$ |
| `LeadingEdgeAnhedralDeg` | `4.0` | `[-35, 35]` | Ventral anhedral / concavity angle $\Gamma_{\text{anh}}$ [deg] |
| `FlatShockFraction` | `0.25` | `[0, 0.85]` | Flat centerbody shock fraction (`'osculating_cone'`) |
| `UpperSurfaceType` | `'canopy'` | `'freestream'`, `'canopy'` | Upper surface lofting mode (waverider methods) |
| `CanopyHeightRatio` | `0.045` | `[0, 0.35]` | Peak dorsal canopy height $h/L$ |
| `RVBodyType` | `'bi_conic'` | `'mono_conic'`, `'bi_conic'`, `'tri_conic'` | Axisymmetric RV conic stage architecture (`Method = 'rv'`) |
| `RVNoseType` | `'spherical'` | `'conical'`, `'spherical'`, `'ogive'` | Axisymmetric RV nose tip geometry (`Method = 'rv'`) |
| `RVNoseBluntnessRatio` | `0.18` | `[0.0, 0.80]` | Spherical nose radius $R_N / R_{\text{base}}$ or ogive fullness (`'rv'`) |
| `RVConeAnglesDeg` | `[11.5, 6.0, 12.0]` | `theta_1 in [1, 55]`, `theta_2,3 in [-12, 55]` | Stage half-angles $[\theta_1, \theta_2, \theta_3]$ [deg] (`'rv'`) |
| `RVStageFractions` | `[0.42, 0.58, 0.0]` | Nonnegative, positive active sum | Relative axial frustum stage weights (`'rv'`) |
| `RVAftType` | `'flat'` | `'flat'`, `'concave'` | Aft base bulkhead geometry (`'rv'`) |
| `RVAftCavityDepthRatio` | `0.35` | `[0.0, 0.85]` | Recessed concave aft depth $d_{\text{aft}} / R_{\text{base}}$ (`'rv'`) |
| `RVAftCavityPower` | `2.0` | `[1.0, 4.0]` | Recessed aft dome exponent $p$ (`2` = paraboloid, `'rv'`) |
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
rcs   = body.computeRCS([0.5 1.5 3 6 10], linspace(0,180,181), -22, 'Polarization', 'HH');
rti   = body.computeRangeProfile(10.0, 1000, linspace(-90,90,121), -22, 'Polarization', 'OP');
fig   = body.plot();                        % 4-panel 3D + Cp + base + symmetry view
fig   = body.plotGeometryLines();           % 3-view lofting ribs & area-rule curve
fig   = body.plotKeyAttributes();           % TM shock layer, AoA polar & face normals
fig   = body.plotRCS('Polarization', 'HH'); % Multi-band UHF-to-X-band RCS dashboard
fig   = body.plotRangeAngleIntensity('Polarization', 'OP', 'WindowType', 'hann');
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
- `surfData = waverider.geometry.generateOsculatingCone(cfg)` — Generates `'conical'` and `'osculating_cone'` forebodies (`surfData.NumClampedStations` reports evolute caustic clamping).
- `surfData = waverider.geometry.generateCaretWedge(cfg)` — Generates exact Nonweiler `'caret'` compression surfaces.
- `surfData = waverider.geometry.generateCSTBody(cfg)` — Generates Kulfan `'cst'` boost-glide surfaces.
- `surfData = waverider.geometry.generateReentryVehicle(cfg)` — Generates axisymmetric `'rv'` mono-conic, bi-conic, and tri-conic reentry vehicles with `'conical'`, $C^1$ `'spherical'`, or $C^1$ `'ogive'` noses, inverse Taylor–Maccoll sharp-cone shock/pressure matching, and `'flat'` or positive-wall-thickness `'concave'` aft bulkheads.
- `surfData = waverider.geometry.buildUpperSurface(surfData, cfg)` — Builds freestream or dorsal-canopy upper surfaces sharing exact leading-edge and wingtip boundaries.
- `surfData = waverider.geometry.appendConformalFrustum(surfData, cfg)` — Appends conformal transition frustum and optional conical rocket stage to either waveriders or RVs.
- `meshData = waverider.geometry.stitchWatertightMesh(surfData, cfg)` — Mirror-symmetric waverider and polar-to-polar axisymmetric RV 2-manifold triangulation stitcher.

### 2.3 `+waverider/+mesh`
- `report = waverider.mesh.verifyWatertight(V, F, patchID, surfStruct)` — Full directed half-edge, Euler, degenerate-triangle, structured-fold, and Divergence-Theorem volume audit.
- `waverider.mesh.exportSTL(V, F, filename, mode)` — Binary (`'binary'`) or ASCII (`'ascii'`) STL exporter.
- `waverider.mesh.exportOBJ(V, F, patchID, filename)` — Patch-grouped Wavefront `.obj` exporter.

### 2.4 `+waverider/+rcs`
- `I = waverider.rcs.facetPhaseIntegral(v1, v2, v3, nHat, w)` — Exact vectorized Gordon (1975) triangle contour phase integral $\iint_{T_f} e^{i\mathbf{w}\cdot\mathbf{r}}\,dA$ for an $N_D \times 3$ wavevector array `w`.
- `pairs = waverider.rcs.buildConcavePairList(body)` — Builds mirror-symmetric planar CCW super-cells on the lower compression surface and enumerates mutually visible concave pairs.
- `[S_mb, contrib, S_mb_pol] = waverider.rcs.doubleBounceGOPO(pairs, uRadar, k0, wantContrib, eHOpt)` — Exact Sutherland–Hodgman clipped GO-PO two-bounce cavity solver returning scalar complex amplitude `S_mb`, full vector polarimetric struct `S_mb_pol` (`HH`, `VV`, `HV`, `VH`), and per-patch `[ampMean, L_half, r_O(1:3), ampHH, ampVV, ampHV, ampVH]` (`[N x 9]`) rows.
- `rcs = waverider.rcs.computeRCS(body, freqGHz, azDeg, elDeg, 'IncludeMultiBounce', true, 'IncludeEdgeDiffraction', false, 'Polarization', 'Mean')` — Full multi-band monostatic RCS sweep supporting `'Polarization'` in `{'Mean', 'HH', 'VV', 'HV', 'VH', 'PP', 'RL', 'LR', 'OP', 'RR', 'LL'}`. Always populates `rcs.Polarization` and `rcs.Polarimetric` (`S_Mean`, `S_HH`, `S_VV`, `S_HV`, `S_VH`, `S_PP`, `S_OP`, `S_LL`, `*_dBsm`, `Span_m2`, `Span_dBsm`) in a single pass.
- `rti = waverider.rcs.computeRangeProfile(body, freqGHz, bandwidthMHz, azVecDeg, elDeg, 'IncludeMultiBounce', true, 'WindowType', 'hann', 'Polarization', 'Mean')` — Wideband stepped-frequency / matched-filter RTI sinogram and 1D HRRP solver (enforces $B \le f_c$; `'HV'`, `'VH'`, `'OP'`, `'LL'` suppress single-bounce PO clutter and isolate delayed cavity double-bounce returns).

### 2.5 `+waverider/+viz`
- `plotWaverider`, `plotGeometryLines`, `plotKeyAttributes`, `plotFamilyComparison`, `plotTradeSpace`, `plotBoosterIntegration`, `plotRCSAnalysis` (accepts `'Polarization'`), `plotRCSConcavityComparison`, `plotRangeAngleIntensity` (accepts `'Polarization'` and `'WindowType'`).
