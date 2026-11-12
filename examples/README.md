# MATLAB Example & Figure-Generation Scripts (`examples/`)

All scripts in `examples/` are self-contained and automatically call `setup_waverider()` to configure the MATLAB path and write publication-resolution PNG figures and STL/OBJ exports into `output/`.

---

## Script Catalog

| Script | Description & Generated Outputs |
| :--- | :--- |
| [`demo_all_waveriders.m`](demo_all_waveriders.m) | Instantiates all 4 baseline presets (`'conical'`, `'osculating_cone'`, `'caret'`, `'cst'`), prints topological & aerodynamic summaries, saves 4-panel overview plots (`output/waverider_*.png`), and exports watertight `.stl` and `.obj` meshes. |
| [`demo_parametric_sweep.m`](demo_parametric_sweep.m) | Sweeps shock wave angle $\beta$ and dorsal canopy height $h/L$ across the Bowcutt viscous/inviscid trade space, identifying the Pareto front between enclosed volume efficiency $\eta_V$ and $L/D$. |
| [`plot_geometry_blueprint.m`](plot_geometry_blueprint.m) | Generates 3-view parametric rib/streamline lofting blueprints, multi-station cross-sections ($X/L = 0.25, 0.50, 0.75, 1.00$), planform thickness contours, and longitudinal area-rule progression (`output/geometry_blueprint_*.png`). |
| [`plot_key_attributes.m`](plot_key_attributes.m) | Plots Taylor–Maccoll conical shock-layer profiles ($C_p(\theta), M(\theta), u_\theta(\theta)$), calibrated angle-of-attack aerodynamic polars, spanwise lift distribution, and 3D outward unit normals (`output/key_attributes_*.png`). |
| [`plot_family_comparison.m`](plot_family_comparison.m) | Renders a side-by-side 4-family geometric, volumetric, and aerodynamic comparison dashboard (`output/waverider_family_comparison.png`). |
| [`plot_booster_integration.m`](plot_booster_integration.m) | Demonstrates conformal aft transition frustums (`'cubic'`, `'quintic'`, `'linear'`, `'conic'`) mating an osculating-cone waverider to a circular booster adapter and conical rocket stage (`output/waverider_booster_integration.png`). |
| [`explore_waverider_rcs.m`](explore_waverider_rcs.m) | Performs full UHF-to-X-band (0.5–10 GHz) Gordon PO + GO-PO concave underside multi-bounce RCS sweeps, 2D Az–El RCS spheres, and wideband Range–Angle–Intensity (RTI) / HRRP bandwidth comparisons (`output/waverider_rcs_*.png`, `output/waverider_range_angle_intensity.png`). |

---

## Running Any Example

From the repository root in MATLAB:
```matlab
setup_waverider();
run('examples/explore_waverider_rcs.m');
```
Or headless from a terminal:
```bash
matlab -batch "setup_waverider(); run('examples/demo_all_waveriders.m');"
```
