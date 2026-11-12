# Interactive Web Studio & JavaScript Physics Engine Guide

The interactive Web Studio in [`docs/`](../docs/) runs 100% client-side in the browser (deployed via GitHub Pages at **https://liploan.github.io/Hypersonics/**) with zero runtime dependencies or build steps, while exposing its entire physics core (`docs/physics.js`) as a UMD module for deterministic headless verification in Node.js.

---

## 1. Module Architecture (`docs/`)

| File | Role |
| :--- | :--- |
| [`docs/physics.js`](../docs/physics.js) | Standalone UMD physics, topology, aerodynamic, CEM, and wideband DSP module (`window.WaveriderPhysics` in browsers, `module.exports` in CommonJS/Node.js). |
| [`docs/app.js`](../docs/app.js) | Tiered `requestAnimationFrame` UI controller, physical slider guardrail enforcer, 6-panel HTML5 Canvas renderer, STL/OBJ exporter, and accessible clipboard handler. |
| [`docs/index.html`](../docs/index.html) | Main 6-panel parametric studio with native `<dialog closedby="any">` MATLAB code generator and R2026b verification gallery. |
| [`docs/references.html`](../docs/references.html) | Searchable, discipline-filterable 24-citation bibliography with verified DOI/NTRS links, source-file cross-reference chips, and 1-click BibTeX copy buttons. |
| [`docs/styles.css`](../docs/styles.css) | Dual-theme (`prefers-color-scheme: dark / light`) stylesheet with WCAG AA contrast ($\ge 4.5:1$), brand-consistent native range sliders, custom scrollbars, and $\ge 24\times 24\text{ px}$ interactive targets. |

---

## 2. JavaScript Physics API (`docs/physics.js`)

`docs/physics.js` can be imported directly in Node.js scripts or tests:

```javascript
import { createRequire } from "node:module";
const require = createRequire(import.meta.url);
const phys = require("./docs/physics.js");

// 1. Flow solvers
const shock = phys.solveObliqueShock(7.0, 13.5, 1.4);
const tm    = phys.solveTaylorMaccoll(7.0, 13.5, 1.4);

// 2. Build 2-manifold watertight waverider + conformal booster model
const model = phys.buildWaveriderModel({
  method: "osculating_cone",
  mach: 7.0,
  beta: 13.5,
  length: 5.5,
  span: 1.75,
  power: 0.62,
  anhedral: 8.5,
  canopy: 0.045,
  frustum: true,
  flength: 1.65,
  bradius: 0.72,
  blend: "cubic",
  rocket: true,
  rcone: 4.5
});
console.log(model.isWatertight, model.euler, model.enclosedVol, model.LD_inviscid);

// 3. Multi-frequency monostatic RCS (shares culling & polygon clipping across bands)
const rcsBands = phys.evaluateRCSMultiFreq(model, [0.5, 1.5, 3.0, 6.0, 10.0], 45.0, -22.0, true);

// 4. Wideband Range-Angle-Intensity (RTI) & HRRP with Harris spectral window
const rti = phys.computeRangeAngleIntensity(model, 10.0, 1000.0, -22.0, true, "hann");
```

---

## 3. Real-Time Browser Optimization Techniques

To evaluate **exact** Gordon (1975) contour Physical Optics + Sutherland–Hodgman clipped GO-PO double-bounce across 5 radar bands (305 azimuthal + 183 elevation + 61 wideband RTI evaluations = **549 full-mesh CEM solves**) inside an interactive ~160 ms frame budget:

1. **Shared Vertex-Phasor Gordon Integration (`prepareLitPOGeometry` + `evaluateLitPOAtK`)**:
   Using the exact analytical edge identity
   $$\operatorname{sinc}_{\text{rad}}\!\left(\frac{\phi_b - \phi_a}{2}\right) e^{i(\phi_a + \phi_b)/2} = \frac{e^{i\phi_b} - e^{i\phi_a}}{i(\phi_b - \phi_a)}$$
   (with a 3-term machine-precision Taylor expansion when $|\phi_b - \phi_a| < 10^{-3}$), `Math.cos` and `Math.sin` are evaluated **once per unique illuminated vertex** rather than 9 times per triangle—eliminating $>85\%$ of transcendental operations with zero approximation error.
2. **Zero-Allocation Ping-Pong Polygon Clipping (`clipQuadIntoFlatBuffer`)**:
   Pre-allocated static `Float64Array` ping-pong buffers, exact ray-cylinder culling ($|(\mathbf{r}_Q - \mathbf{r}_P)\times\hat{\mathbf{k}}_r| \le r_P + r_Q$), and 3D AABB pre-rejection eliminate heap garbage collection during multi-angle GO-PO cavity sweeps.
3. **Tiered `requestAnimationFrame` Coalescing (`scheduleUpdate` in `docs/app.js`)**:
   - `'view3d'` (~1 ms): Camera orbit/zoom or 3D wireframe/shock/ray checkboxes only repaint `#canvas-3d`.
   - `'cut'` (~2 ms): Dragging the RTI heatmap or `#rng-rtiaz` only repaints `#canvas-rti` and `#canvas-hrrp`.
   - `'rcs'`: Changing radar frequency, bandwidth, elevation, or window updates RCS/RTI without rebuilding the 3D mesh.
   - `'geometry'`: Rebuilds the watertight mesh and updates all 6 panels.

---

## 4. Dynamic Physical Guardrails & Accessibility

- **Mach-Cone Shock Floor**: Whenever `M∞` changes, `#rng-beta` dynamically updates its minimum to $\lceil 2(\arcsin(1/M_\infty) + 1^\circ)\rceil / 2$, preventing unphysical detached/sub-Mach shock angles.
- **Positive-Frequency Bandwidth Ceiling**: Whenever `f_c` changes, `#rng-bw` dynamically caps $B \le 0.95\,f_c$ (e.g., $450\text{ MHz}$ at $0.5\text{ GHz}$ UHF), guaranteeing strictly positive spectral support $[f_c - B/2, f_c + B/2]$.
- **WCAG AA & Assistive Technology**: Semantic landmarks (`<header>`, `<nav aria-label="Primary">`, `<aside>`, `<main>`, `<footer>`), descriptive `role="img"` + `aria-label` on all 6 canvases, non-color-only topology status badges (`✓ PASS` vs `⚠ FAIL`), and external sibling `role="status"` copy tooltips with VoiceOver repeat-announcement deduplication.
