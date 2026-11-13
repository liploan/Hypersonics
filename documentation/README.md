# Technical Documentation Hub

Welcome to the technical documentation suite for the **Parametric Watertight Hypersonic Waverider & UHF–X Band RCS Studio**. Every mathematical formulation, coordinate frame, numerical integrator, API signature, validation test, and bibliographic citation used in `+waverider/` (MATLAB) and `docs/` (Web Studio) is documented in the manuals below.

---

## Documentation Index

| Document | Contents & Intended Audience |
| :--- | :--- |
| [**`THEORY_AND_METHODS.md`**](THEORY_AND_METHODS.md) | Coordinate conventions, exact 2D oblique-shock & 3D axisymmetric Taylor–Maccoll ODEs, streamline integration, all 4 waverider lofting families, Sobieczky evolute caustic guard ($d_{\text{LE}} \le 0.82\,R_c$), conventional axisymmetric Reentry Vehicles (`mono_conic`, `bi_conic`, `tri_conic`; `conical`, $C^1$ `spherical`, $C^1$ `ogive` noses; `flat` & `concave` aft bulkheads), conformal booster frustum blending, 2-manifold stitching & generalized winding numbers, calibrated Newtonian + Prandtl–Meyer + Eckert $T^*$ aerodynamics, Gordon (1975) contour PO, vector dyadic Sutherland–Hodgman GO-PO double-bounce across linear (`HH`, `VV`, `HV`, `VH`) and circular (`PP`, `OP`) polarizations, Harris (1978) wideband HRRP/RTI signal processing, and explicit modeling limitations. |
| [**`MATLAB_API_REFERENCE.md`**](MATLAB_API_REFERENCE.md) | Complete property, method, input/output struct, and name-value option reference for `waverider.WaveriderConfig` (including RV presets & properties), `waverider.WaveriderBody`, `+flow`, `+geometry`, `+mesh`, `+rcs`, and `+viz`. |
| [**`WEB_STUDIO_GUIDE.md`**](WEB_STUDIO_GUIDE.md) | Architecture of the zero-dependency GitHub Pages Web Studio (`docs/index.html`, `docs/physics.js`, `docs/app.js`, `docs/references.html`, `docs/styles.css`), Waverider / Reentry Vehicle tab switcher & meridional blueprint view, UMD Node/Browser API, vertex-phasor PO acceleration, polarimetric channel switching, dynamic physical slider bounds, tiered `requestAnimationFrame` scheduling, and WCAG AA accessibility design. |
| [**`VALIDATION_REPORT.md`**](VALIDATION_REPORT.md) | Methodology, exact pass/fail tolerances, reproduction commands, and full quantitative results for `tests/test_waverider_suite.m` (16/16), `tests/test_physics_validation.m` (33/33), and `tests/web/validate_physics.mjs` (36/36). |
| [**`REFERENCES.md`**](REFERENCES.md) | Complete Markdown bibliography of all **24 primary-source verified references** (NACA Report 1135, Taylor & Maccoll, Anderson, Dormand & Prince, Shampine & Reichelt, Nonweiler, Jones et al., Rasmussen, Bowcutt et al., Sobieczky et al., Kulfan, Ding et al., Lees, Eckert, Gordon, Ling et al., Knott et al., Ufimtsev, Harris, Wehner, Sutherland & Hodgman, van Oosterom & Strackee, Jacobson et al., Zhang & Chen) with verified DOIs/NTRS links, code cross-references, and BibTeX entries. |

---

## Subdirectory Guides

* [**`examples/README.md`**](../examples/README.md) — Walkthrough of the 7 standalone MATLAB demonstration and publication-figure scripts.
* [**`tests/README.md`**](../tests/README.md) — Instructions for running the MATLAB and Node.js automated test suites and interpreting output diagnostics.
* [**`docs/README.md`**](../docs/README.md) — Guide to the live GitHub Pages deployment, local static previewing, and headless Node.js usage of `docs/physics.js`.
