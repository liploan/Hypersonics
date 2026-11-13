/**
 * Independent Node.js Validation Suite for the Web Studio Physics Engine (docs/physics.js)
 *
 * Mirrors all 23 checks in tests/test_physics_validation.m plus Web-specific
 * guardrails and real-time latency benchmarks:
 *   A1–A4: Oblique shock & Taylor-Maccoll ODE vs SciPy & NACA Report 1135
 *   B1–B5: 42 configurations (5 presets + 24 randomized + 13 validator-boundary
 *          edge cases) verified for 2-manifold directed half-edges, Euler = 2,
 *          zero degenerate faces, zero surface folds/crossovers, generalized
 *          winding number W in {0, 1}, and Monte-Carlo vs divergence volume
 *   C1–C5: Gordon triangle PO & GO-PO Sutherland-Hodgman double-bounce vs
 *          exact canonical targets (flat plate normal + sinc^2, PEC sphere,
 *          90-deg dihedral symmetric & off-symmetric)
 *   D1–D6: Wideband HRRP / RTI positions, -3 dB widths, and Harris sidelobes
 *          for rectangular and Hann windows + B <= f_c guard
 *   E1–E4: Nonweiler caret L/D = cot(delta), polar(0) == design aero, bilateral
 *          symmetry (|CY| = 0, RCS(+az) == RCS(-az)), and < 250 ms frame budget
 *   F1–F5: Full polarimetric Sinclair scattering matrix (HH, VV, HV, VH) and
 *          circular (PP, OP) benchmarks (PEC plate, 0/45 deg dihedral, reciprocity,
 *          Parseval span invariance, and OP/HV 1-bounce clutter rejection)
 *   G1–G5: Conventional axisymmetric Reentry Vehicles (RVs): 19-config 2-manifold
 *          topology & recessed-cavity winding audit, exact cone/paraboloid volumes,
 *          Taylor-Maccoll cone drag & zero-lift axisymmetry, circular-disk tail-on
 *          PO RCS, and concave aft cavity multi-bounce OP + delayed HRRP peak
 */

import { createRequire } from "node:module";
import { performance } from "node:perf_hooks";

const require = createRequire(import.meta.url);
const phys = require("../../docs/physics.js");

const results = [];
function record(name, pass, detail) {
  results.push({ name, pass: Boolean(pass), detail });
  console.log(`${pass ? "[PASS]" : "[FAIL]"} ${name.padEnd(58)} ${detail}`);
  if (!pass) {
    process.exitCode = 1;
  }
}

// Deterministic PRNG (Mulberry32) for reproducible randomized sweep
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0;
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

console.log("\n=====================================================================");
console.log("  WEB STUDIO PHYSICS / TOPOLOGY / SIGNATURE VALIDATION (Node.js)");
console.log("=====================================================================");

// ---------------------------------------------------------------------------
// A. FLOW SOLVERS vs. INDEPENDENT SCIPY & NACA-1135
// ---------------------------------------------------------------------------
const refTable = [
  [6.0, 15.0, 7.18494, 5.02407, 2.64680, 10.73711, 0.08164],
  [7.0, 13.5, 6.94351, 5.76651, 2.94873,  9.95551, 0.06930],
  [8.0, 13.0, 7.40949, 6.32931, 3.61169, 10.01540, 0.06867],
  [5.0, 18.0, 8.48714, 4.14905, 2.61850, 12.79755, 0.11583],
  [10.0, 10.0, 5.53325, 8.10685, 3.35126, 7.60744, 0.04003],
  [2.0, 40.0, 10.62291, 1.61732, 1.76149, 22.43032, 0.39437],
  [3.5, 24.0, 9.57164, 2.92968, 2.19767, 15.92470, 0.18411]
];

let maxErrDelta = 0, maxErrM2 = 0, maxErrCone = 0, maxErrCp = 0;
for (const [M, beta, dRef, m2Ref, , tcRef, cpRef] of refTable) {
  const sh = phys.solveObliqueShock(M, beta);
  const tm = phys.solveTaylorMaccoll(M, beta);
  maxErrDelta = Math.max(maxErrDelta, Math.abs(sh.deltaDeg - dRef));
  maxErrM2 = Math.max(maxErrM2, Math.abs(sh.M2 - m2Ref));
  maxErrCone = Math.max(maxErrCone, Math.abs(tm.thetaConeDeg - tcRef));
  maxErrCp = Math.max(maxErrCp, Math.abs(tm.CpCone - cpRef));
}
record(
  "A1 Oblique shock delta & M2 vs SciPy",
  maxErrDelta < 1e-4 && maxErrM2 < 1e-4,
  `max|d delta| = ${maxErrDelta.toExponential(2)} deg, max|d M2| = ${maxErrM2.toExponential(2)}`
);
record(
  "A2 Taylor-Maccoll cone angle vs SciPy",
  maxErrCone < 2e-3,
  `max|d theta_c| = ${maxErrCone.toExponential(2)} deg`
);
record(
  "A3 Taylor-Maccoll cone-surface Cp vs SciPy",
  maxErrCp < 2e-4,
  `max|d Cp_c| = ${maxErrCp.toExponential(2)}`
);

// Bisection for NACA 1135: M=2.0, theta_c = 20.0 deg -> beta = 37.80 deg
let bLo = 32.0, bHi = 50.0;
for (let it = 0; it < 40; it++) {
  const bMid = 0.5 * (bLo + bHi);
  if (phys.solveTaylorMaccoll(2.0, bMid).thetaConeDeg < 20.0) bLo = bMid;
  else bHi = bMid;
}
const betaNACA = 0.5 * (bLo + bHi);
record(
  "A4 NACA-1135 chart: M=2, 20deg cone -> beta=37.80",
  Math.abs(betaNACA - 37.80) < 0.05,
  `beta = ${betaNACA.toFixed(3)} deg`
);

// Sub-Mach shock guardrail check
let subMachThrew = false;
try {
  phys.buildWaveriderModel({ method: "osculating_cone", mach: 3.5, beta: 15.0 });
} catch {
  subMachThrew = true;
}
record("A5 Sub-Mach shock angle guard (M=3.5, beta=15 deg)", subMachThrew, "RangeError thrown as expected");

// ---------------------------------------------------------------------------
// B. 42-CONFIG GEOMETRY SWEEP + WINDING NUMBER + MONTE-CARLO VOLUME
// ---------------------------------------------------------------------------
const PRESETS = {
  conical: {
    method: "conical", mach: 6.0, beta: 14.5, length: 5.0, span: 1.55,
    power: 0.68, anhedral: 6.0, canopy: 0.0, frustum: false, flength: 1.5,
    bradius: 0.65, blend: "cubic", rocket: false, rcone: 4.0
  },
  osculating_cone: {
    method: "osculating_cone", mach: 7.0, beta: 13.5, length: 5.5, span: 1.75,
    power: 0.62, anhedral: 8.0, canopy: 0.045, frustum: false, flength: 1.65,
    bradius: 0.72, blend: "cubic", rocket: false, rcone: 4.5
  },
  caret: {
    method: "caret", mach: 6.5, beta: 15.0, length: 5.0, span: 1.5,
    power: 0.85, anhedral: 11.0, canopy: 0.035, frustum: false, flength: 1.5,
    bradius: 0.65, blend: "cubic", rocket: false, rcone: 4.5
  },
  cst: {
    method: "cst", mach: 8.0, beta: 12.5, length: 6.0, span: 1.8,
    power: 0.6, anhedral: 7.0, canopy: 0.042, frustum: false, flength: 1.7,
    bradius: 0.75, blend: "cubic", rocket: false, rcone: 4.5
  },
  booster_mounted: {
    method: "osculating_cone", mach: 7.0, beta: 13.5, length: 5.5, span: 1.75,
    power: 0.62, anhedral: 8.5, canopy: 0.045, frustum: true, flength: 1.65,
    bradius: 0.72, blend: "cubic", rocket: true, rcone: 4.5
  }
};

const rng = mulberry32(7);
const methods = ["conical", "osculating_cone", "caret", "cst"];
const blends = ["cubic", "quintic", "linear", "conic"];
const cfgCases = [];

for (const [k, p] of Object.entries(PRESETS)) {
  cfgCases.push({ name: `preset:${k}`, cfg: { ...p } });
}
for (let n = 1; n <= 24; n++) {
  const method = methods[Math.floor(rng() * 4)];
  const mach = 3.5 + rng() * 11.5;
  const muDeg = (Math.asin(1 / mach) * 180) / Math.PI;
  const beta = Math.max(10 + rng() * 14, muDeg + 2.0);
  const frustum = rng() > 0.4;
  const canopy = rng() * 0.12;
  cfgCases.push({
    name: `random#${String(n).padStart(2, "0")}:${method}`,
    cfg: {
      method, mach, beta,
      length: 2.5 + rng() * 7.5,
      span: 0.8 + rng() * 2.7,
      power: 0.35 + rng() * 0.9,
      anhedral: rng() * 20,
      canopy: canopy < 0.01 ? 0 : canopy,
      frustum,
      flength: 0.5 + rng() * 3.0,
      bradius: 0.3 + rng() * 1.1,
      blend: blends[Math.floor(rng() * 4)],
      rocket: frustum && rng() > 0.4,
      rcone: rng() * 12
    }
  });
}

const edgeDefs = [
  ["conical", { power: 0.25 }],
  ["osculating_cone", { power: 0.25 }],
  ["caret", { power: 0.25 }],
  ["cst", { power: 0.25 }],
  ["conical", { power: 2.5 }],
  ["osculating_cone", { power: 2.5, anhedral: 35 }],
  ["caret", { anhedral: 35 }],
  ["caret", { anhedral: -35 }],
  ["osculating_cone", { mach: 3.5, beta: 18.5 }],
  ["booster_mounted", { bradius: 0.3, flength: 0.5 }],
  ["booster_mounted", { bradius: 1.4, flength: 0.5 }],
  ["booster_mounted", { rcone: -15 }],
  ["cst", { canopy: 0.35 }]
];
for (const [baseName, overrides] of edgeDefs) {
  cfgCases.push({
    name: `edge:${baseName}${JSON.stringify(overrides)}`,
    cfg: { ...PRESETS[baseName], ...overrides }
  });
}

const geoFails = [];
let worstWinding = 0;
let worstCross = 0;
let totalDegen = 0;
let worstVolSigma = 0;

for (const item of cfgCases) {
  let m;
  try {
    m = phys.buildWaveriderModel(item.cfg);
  } catch (err) {
    geoFails.push(`${item.name}: threw ${err.message}`);
    continue;
  }

  let cross = -Infinity;
  for (let i = 0; i < m.NxWave; i++) {
    for (let k = 0; k < m.Ny; k++) {
      cross = Math.max(cross, m.Z_low[i][k] - m.Z_up[i][k]);
    }
  }
  worstCross = Math.max(worstCross, cross);
  totalDegen += m.degenerateFaces;

  const isPreset = item.name.startsWith("preset:");
  const nPts = isPreset ? 1200 : 250;
  let insideCount = 0;
  let maxDevLocal = 0;
  for (let p = 0; p < nPts; p++) {
    const q = [
      m.xMin + rng() * (m.xMax - m.xMin),
      m.yMin + rng() * (m.yMax - m.yMin),
      m.zMin + rng() * (m.zMax - m.zMin)
    ];
    const W = phys.generalizedWindingNumber(m.vertices, m.faces, q);
    const dev = Math.min(Math.abs(W), Math.abs(W - 1.0));
    if (dev > maxDevLocal) maxDevLocal = dev;
    if (Math.abs(W - 1.0) < 0.25) insideCount++;
  }
  worstWinding = Math.max(worstWinding, maxDevLocal);

  const bboxVol = (m.xMax - m.xMin) * (m.yMax - m.yMin) * (m.zMax - m.zMin);
  const pIn = insideCount / nPts;
  const volMC = pIn * bboxVol;
  const volSig = Math.sqrt(Math.max(pIn * (1 - pIn), 1e-6) / nPts) * bboxVol;
  const zScore = Math.abs(volMC - m.enclosedVol) / Math.max(volSig, 1e-9);
  if (isPreset) worstVolSigma = Math.max(worstVolSigma, zScore);

  const probs = [];
  if (!m.isWatertight) probs.push("isWatertight=false");
  if (m.boundaryEdges !== 0) probs.push(`${m.boundaryEdges} boundary edges`);
  if (m.nonManifoldEdges !== 0) probs.push(`${m.nonManifoldEdges} non-manifold edges`);
  if (m.conflictingHalfEdges !== 0) probs.push(`${m.conflictingHalfEdges} conflicting half-edges`);
  if (m.euler !== 2) probs.push(`euler=${m.euler}`);
  if (m.hasSurfaceFold) probs.push("hasSurfaceFold=true");
  if (maxDevLocal > 0.05) probs.push(`winding dev=${maxDevLocal.toFixed(3)}`);
  if (cross > 1e-6) probs.push(`crossover=${cross.toExponential(2)}`);
  if (probs.length > 0) geoFails.push(`${item.name}: ${probs.join(", ")}`);
}

record("B1 Geometry sweep: all 42 configs generate & close", geoFails.length === 0, `${cfgCases.length} configs, ${geoFails.length} issues`);
for (const f of geoFails) console.log(`     -> ${f}`);
record("B2 No self-intersection (winding number in {0,1})", worstWinding < 0.05, `worst |W - {0,1}| = ${worstWinding.toFixed(3)}`);
record("B3 No upper/lower surface crossover", worstCross < 1e-6, `worst Z_low - Z_up = ${worstCross.toExponential(2)} m`);
record("B4 No degenerate (zero-area) triangles", totalDegen === 0, `${totalDegen} degenerate faces across 42 configs`);
record("B5 Monte-Carlo volume == divergence volume (presets)", worstVolSigma < 4.0, `worst deviation = ${worstVolSigma.toFixed(2)} sigma`);

// ---------------------------------------------------------------------------
// C. CANONICAL ELECTROMAGNETIC TARGETS (PO & GO-PO)
// ---------------------------------------------------------------------------
// C1/C2: 1m x 1m Square PEC Plate in X-Y plane (2 coarse triangles, 33 lambda @ 10 GHz)
{
  const a = 1.0;
  const Vp = [[-0.5, -0.5, 0], [0.5, -0.5, 0], [0.5, 0.5, 0], [-0.5, 0.5, 0]];
  const Fp = [[0, 1, 2], [0, 2, 3]];
  const plate = phys.buildCanonicalTargetModel(Vp, Fp, null, 1.0);
  const fGHz = 10.0;
  const lam = (phys.C0 * 1e-9) / fGHz;
  const k = (2 * Math.PI) / lam;
  const sig0 = (4 * Math.PI * Math.pow(a, 4)) / (lam * lam);

  const r0 = phys.evaluateRCSDirection(plate, fGHz, 0, 90, false);
  const err0 = r0.total_dBsm - 10 * Math.log10(sig0);
  record("C1 Flat plate normal incidence: 4 pi A^2 / lambda^2", Math.abs(err0) < 0.1, `err = ${err0.toFixed(3)} dB`);

  const thList = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0];
  const errs = thList.map((thDeg) => {
    const th = (thDeg * Math.PI) / 180;
    const rr = phys.evaluateRCSDirection(plate, fGHz, 0, 90 - thDeg, false);
    const x = k * a * Math.sin(th);
    const exact = sig0 * Math.pow(Math.cos(th), 2) * Math.pow(Math.sin(x) / x, 2);
    return rr.total_dBsm - 10 * Math.log10(Math.max(exact, 1e-5));
  });
  const maxSincErr = Math.max(...errs.map(Math.abs));
  record(
    "C2 Flat plate off-normal sinc^2 (2 coarse facets, 33 lambda)",
    maxSincErr < 0.5,
    `errors [dB] = [${errs.map((e) => e.toFixed(2)).join(", ")}]`
  );
}

// C3: PEC Sphere radius 1 m (icosphere subdivision 4 -> 5120 faces) -> pi a^2 = 4.9715 dBsm
{
  function buildIcosphere(nSub) {
    const t = (1 + Math.sqrt(5)) / 2;
    let V = [
      [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
      [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
      [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]
    ].map((v) => {
      const l = Math.hypot(v[0], v[1], v[2]);
      return [v[0] / l, v[1] / l, v[2] / l];
    });
    let F = [
      [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
      [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
      [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
      [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]
    ];
    for (let s = 0; s < nSub; s++) {
      const midCache = new Map();
      function midIdx(i, j) {
        const key = i < j ? `${i}_${j}` : `${j}_${i}`;
        if (midCache.has(key)) return midCache.get(key);
        const mx = 0.5 * (V[i][0] + V[j][0]);
        const my = 0.5 * (V[i][1] + V[j][1]);
        const mz = 0.5 * (V[i][2] + V[j][2]);
        const l = Math.hypot(mx, my, mz);
        V.push([mx / l, my / l, mz / l]);
        const idx = V.length - 1;
        midCache.set(key, idx);
        return idx;
      }
      const nextF = [];
      for (const [a, b, c] of F) {
        const ab = midIdx(a, b), bc = midIdx(b, c), ca = midIdx(c, a);
        nextF.push([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]);
      }
      F = nextF;
    }
    return { V, F };
  }

  const { V, F } = buildIcosphere(4); // 5120 faces
  const sphere = phys.buildCanonicalTargetModel(V, F, null, 2.0);
  const exactSphereDB = 10 * Math.log10(Math.PI);
  const err1 = Math.abs(phys.evaluateRCSDirection(sphere, 1.0, 37, 21, false).total_dBsm - exactSphereDB);
  const err3 = Math.abs(phys.evaluateRCSDirection(sphere, 3.0, 90, -50, false).total_dBsm - exactSphereDB);
  record(
    "C3 PEC sphere a=1m (1, 3 GHz) -> pi a^2 = 4.97 dBsm",
    err1 < 1.0 && err3 < 1.0,
    `err = ${err1.toFixed(2)} dB (1 GHz), ${err3.toFixed(2)} dB (3 GHz)`
  );
}

// C4/C5: 90-deg Dihedral Trough opening downward (a=0.5m, b=1.0m, 10 GHz)
{
  const aD = 0.5, bD = 1.0, fD = 10.0;
  const lamD = (phys.C0 * 1e-9) / fD;
  const yMax = aD / Math.SQRT2;
  const nX = 24, nY = 25; // Odd nY so fold apex Y=0 is an exact grid column!
  const X = Array.from({ length: nX }, (_, i) => Array.from({ length: nY }, () => (i * bD) / (nX - 1)));
  const Y = Array.from({ length: nX }, () => Array.from({ length: nY }, (_, j) => -yMax + (2 * yMax * j) / (nY - 1)));
  const Z = Array.from({ length: nX }, (_, i) => Array.from({ length: nY }, (_, j) => -Math.abs(Y[i][j])));
  const Vd = [];
  const idGrid = Array.from({ length: nX }, () => new Int32Array(nY));
  for (let i = 0; i < nX; i++) {
    for (let j = 0; j < nY; j++) {
      idGrid[i][j] = Vd.length;
      Vd.push([X[i][j], Y[i][j], Z[i][j]]);
    }
  }
  const Fd = [];
  for (let i = 0; i < nX - 1; i++) {
    for (let j = 0; j < nY - 1; j++) {
      const p1 = idGrid[i][j], p2 = idGrid[i + 1][j], p3 = idGrid[i + 1][j + 1], p4 = idGrid[i][j + 1];
      Fd.push([p1, p3, p2], [p1, p4, p3]);
    }
  }
  const dih = phys.buildCanonicalTargetModel(Vd, Fd, { X, Y, Z }, bD);
  const rSym = phys.evaluateRCSDirection(dih, fD, 90, -90, true);
  const sigSym = (16 * Math.PI * aD * aD * bD * bD * Math.pow(Math.sin(Math.PI / 4), 2)) / (lamD * lamD);
  const errSym = rSym.total_dBsm - 10 * Math.log10(sigSym);
  record(
    "C4 90deg dihedral double-bounce, symmetric (phi=45 deg)",
    Math.abs(errSym) < 1.0,
    `model ${rSym.total_dBsm.toFixed(2)} vs exact ${(10 * Math.log10(sigSym)).toFixed(2)} dBsm (err ${errSym >= 0 ? "+" : ""}${errSym.toFixed(2)} dB)`
  );

  const rOff = phys.evaluateRCSDirection(dih, fD, 90, -75, true);
  const sigOff = (16 * Math.PI * aD * aD * bD * bD * Math.pow(Math.sin(Math.PI / 6), 2)) / (lamD * lamD);
  const errOff = rOff.total_dBsm - 10 * Math.log10(sigOff);
  record(
    "C5 90deg dihedral double-bounce, 15deg off-symmetric (phi=30)",
    Math.abs(errOff) < 1.0,
    `model ${rOff.total_dBsm.toFixed(2)} vs exact ${(10 * Math.log10(sigOff)).toFixed(2)} dBsm (err ${errOff >= 0 ? "+" : ""}${errOff.toFixed(2)} dB)`
  );
}

// ---------------------------------------------------------------------------
// D. RANGE PROFILE vs. ANALYTIC HARRIS WINDOW PSFs
// ---------------------------------------------------------------------------
{
  const s = 0.04;
  const sq = [[0, -s, -s], [0, s, -s], [0, s, s], [0, -s, s]];
  const Vr = [
    ...sq.map((p) => [p[0] + 2.0, p[1], p[2]]),
    ...sq.map((p) => [p[0] + 4.0, p[1], p[2]])
  ];
  const Fr = [[0, 2, 1], [0, 3, 2], [4, 6, 5], [4, 7, 6]];
  const pair = phys.buildCanonicalTargetModel(Vr, Fr, null, 4.0);
  const B = 1000.0;
  const dR = phys.C0 / (2 * B * 1e6);

  for (const win of ["rect", "hann"]) {
    const rp = phys.computeRangeAngleIntensity(pair, 10.0, B, 0.0, false, win, {
      azVec: [0],
      minRangeBins: 4001,
      maxRangeBins: 8001
    });
    const R = rp.rangeVec;
    const profDB = rp.RTI_total_dBsm[0];
    const profLin = profDB.map((db) => Math.pow(10, db / 10));
    const Rexp = [2.0 - rp.xRef, 4.0 - rp.xRef];

    let pk1 = -1, i1 = 0, pk2 = -1, i2 = 0;
    for (let i = 0; i < R.length; i++) {
      if (Math.abs(R[i] - Rexp[0]) < 0.5 && profLin[i] > pk1) { pk1 = profLin[i]; i1 = i; }
      if (Math.abs(R[i] - Rexp[1]) < 0.5 && profLin[i] > pk2) { pk2 = profLin[i]; i2 = i; }
    }
    const posErr = Math.max(Math.abs(R[i1] - Rexp[0]), Math.abs(R[i2] - Rexp[1]));
    const half = pk1 / 2;
    let lo = i1, hi = i1;
    while (lo > 0 && profLin[lo] >= half) lo--;
    while (hi < R.length - 1 && profLin[hi] >= half) hi++;
    const w3 = R[hi] - R[lo];
    const wExp = (win === "rect" ? 0.886 : 1.44) * dR;
    const slExp = win === "rect" ? -13.26 : -31.5;
    const mainNull = (win === "rect" ? 1.0 : 2.0) * dR;

    let maxSL = 0;
    for (let i = 0; i < R.length; i++) {
      const dist = Math.abs(R[i] - R[i1]);
      if (dist > mainNull && dist < 6 * dR && profLin[i] > maxSL) maxSL = profLin[i];
    }
    const sl = 10 * Math.log10(maxSL / pk1);

    record(`D1 Range profile (${win}): scatterer positions`, posErr < 0.01, `max pos err = ${posErr.toFixed(4)} m`);
    record(`D2 Range profile (${win}): 3-dB width`, Math.abs(w3 - wExp) / wExp < 0.08, `${(100 * w3).toFixed(2)} cm vs ${(100 * wExp).toFixed(2)} cm`);
    record(`D3 Range profile (${win}): peak sidelobe`, Math.abs(sl - slExp) < 2.0, `${sl.toFixed(1)} dB vs ${slExp.toFixed(1)} dB`);
  }

  let bwThrew = false;
  try {
    phys.computeRangeAngleIntensity(pair, 0.5, 1000.0, 0.0, false, "hann", { strictBandwidth: true });
  } catch {
    bwThrew = true;
  }
  record("D4 Bandwidth B > f_c guard (1000 MHz @ 0.5 GHz)", bwThrew, "RangeError thrown in strict mode");
}

// ---------------------------------------------------------------------------
// E. EXACT AERODYNAMICS, SYMMETRY & REAL-TIME FRAME BUDGET
// ---------------------------------------------------------------------------
{
  const mCaret = phys.buildWaveriderModel({
    method: "caret", mach: 6.5, beta: 15.0, power: 1.0, canopy: 0.0, upperSurfaceType: "freestream"
  });
  const sh = phys.solveObliqueShock(6.5, 15.0, 1.4);
  const clEx = sh.Cp2;
  const cdEx = sh.Cp2 * Math.tan(sh.deltaRad);
  const ldEx = 1.0 / Math.tan(sh.deltaRad);
  const errCaret = Math.max(Math.abs(mCaret.CL - clEx), Math.abs(mCaret.CD_wave - cdEx), Math.abs(mCaret.LD_inviscid - ldEx));
  record(
    "E1 Nonweiler caret theorem: CL=Cp2, L/D=cot(delta)",
    errCaret < 1e-11,
    `L/D = ${mCaret.LD_inviscid.toFixed(6)} vs cot(delta) = ${ldEx.toFixed(6)} (err ${errCaret.toExponential(2)})`
  );

  let maxPolDiff = 0;
  let maxCY = 0;
  for (const p of Object.values(PRESETS)) {
    const mp = phys.buildWaveriderModel(p);
    const pol0 = phys.computeAeroPolar(mp, [0.0])[0];
    maxPolDiff = Math.max(maxPolDiff, Math.abs(pol0.CL - mp.CL), Math.abs(pol0.CD_wave - mp.CD_wave));
    maxCY = Math.max(maxCY, Math.abs(mp.CY));
  }
  record("E2 AeroMetrics == computeAeroPolar(0) identity", maxPolDiff < 1e-13, `max diff = ${maxPolDiff.toExponential(2)}`);

  const mSym = phys.buildWaveriderModel(PRESETS.osculating_cone);
  const rNeg55 = phys.evaluateRCSDirection(mSym, 10.0, -55, -20, true).total_dBsm;
  const rPos55 = phys.evaluateRCSDirection(mSym, 10.0, +55, -20, true).total_dBsm;
  const rNeg25 = phys.evaluateRCSDirection(mSym, 10.0, -25, -20, true).total_dBsm;
  const rPos25 = phys.evaluateRCSDirection(mSym, 10.0, +25, -20, true).total_dBsm;
  const asymDB = Math.max(Math.abs(rNeg55 - rPos55), Math.abs(rNeg25 - rPos25));
  record(
    "E3 Bilateral mesh symmetry: |CY|=0 & RCS(+az)==RCS(-az)",
    maxCY < 1e-13 && asymDB < 1e-11,
    `max|CY| = ${maxCY.toExponential(2)}, max|dRCS(+/-az)| = ${asymDB.toExponential(2)} dB`
  );

  // Real-time interactive frame budget benchmark (full model build + 5-band az + 3-band el + 61-az RTI)
  const t0 = performance.now();
  const mBench = phys.buildWaveriderModel(PRESETS.booster_mounted);
  const t1 = performance.now();
  const azBands = [0.5, 1.5, 3.0, 6.0, 10.0];
  for (let az = 0; az <= 180; az += 3) {
    phys.evaluateRCSMultiFreq(mBench, azBands, az, -22, true, "OP");
  }
  const t2 = performance.now();
  const elBands = [0.5, 3.0, 10.0];
  for (let el = -90; el <= 90; el += 3) {
    phys.evaluateRCSMultiFreq(mBench, elBands, 45, el, true, "OP");
  }
  const t3 = performance.now();
  phys.computeRangeAngleIntensity(mBench, 10.0, 1000.0, -22, true, "hann", { polarization: "OP" });
  const t4 = performance.now();
  const dtMs = t4 - t0;
  record(
    "E4 Full Web Studio rebuild + multi-band RCS + RTI latency",
    dtMs < 250,
    `total ${dtMs.toFixed(1)} ms (mesh ${(t1 - t0).toFixed(1)}, az ${(t2 - t1).toFixed(1)}, el ${(t3 - t2).toFixed(1)}, rti ${(t4 - t3).toFixed(1)})`
  );
}

// ---------------------------------------------------------------------------
// F. POLARIMETRIC SCATTERING MATRIX (HH, VV, HV, VH) & CIRCULAR (PP, OP)
// ---------------------------------------------------------------------------
{
  // F1: Flat PEC plate polarimetry: HH == VV == PP == exact, HV == VH == OP == 0 (-50 dBsm floor)
  const a = 1.0, fGHz = 10.0;
  const lam = (phys.C0 * 1e-9) / fGHz;
  const sig0DB = 10 * Math.log10((4 * Math.PI * Math.pow(a, 4)) / (lam * lam));
  const Vp = [[-0.5, -0.5, 0], [0.5, -0.5, 0], [0.5, 0.5, 0], [-0.5, 0.5, 0]];
  const Fp = [[0, 1, 2], [0, 2, 3]];
  const plate = phys.buildCanonicalTargetModel(Vp, Fp, null, 1.0);

  const pHH = phys.evaluateRCSDirection(plate, fGHz, 0, 90, false, "HH");
  const pVV = phys.evaluateRCSDirection(plate, fGHz, 0, 90, false, "VV");
  const pHV = phys.evaluateRCSDirection(plate, fGHz, 0, 90, false, "HV");
  const pPP = phys.evaluateRCSDirection(plate, fGHz, 0, 90, false, "PP");
  const pOP = phys.evaluateRCSDirection(plate, fGHz, 0, 90, false, "OP");

  const coErr = Math.max(
    Math.abs(pHH.total_dBsm - sig0DB),
    Math.abs(pVV.total_dBsm - sig0DB),
    Math.abs(pPP.total_dBsm - sig0DB)
  );
  const nullFloorErr = Math.max(Math.abs(pHV.total_dBsm - -50.0), Math.abs(pOP.total_dBsm - -50.0));
  record(
    "F1 PEC plate polarimetry: HH=VV=PP=exact, HV=OP=0",
    coErr < 0.05 && nullFloorErr < 1e-12,
    `co-pol err = ${coErr.toExponential(2)} dB, HV/OP = ${pOP.total_dBsm.toFixed(2)} dBsm`
  );

  // Build helper for 90-deg PEC dihedral rotated in X-Y plane by roll angle alphaDeg
  function buildRotatedDihedral(alphaDeg) {
    const aD = 0.5, bD = 1.0, nX = 24, nY = 25;
    const yMax = aD / Math.SQRT2;
    const ca = Math.cos((alphaDeg * Math.PI) / 180);
    const sa = Math.sin((alphaDeg * Math.PI) / 180);
    const X = Array.from({ length: nX }, () => new Float64Array(nY));
    const Y = Array.from({ length: nX }, () => new Float64Array(nY));
    const Z = Array.from({ length: nX }, () => new Float64Array(nY));
    const Vd = [];
    const idGrid = Array.from({ length: nX }, () => new Int32Array(nY));
    for (let i = 0; i < nX; i++) {
      const x0 = (i * bD) / (nX - 1) - 0.5 * bD;
      for (let j = 0; j < nY; j++) {
        const y0 = -yMax + (2 * yMax * j) / (nY - 1);
        const z0 = -Math.abs(y0);
        const xr = x0 * ca - y0 * sa;
        const yr = x0 * sa + y0 * ca;
        X[i][j] = xr; Y[i][j] = yr; Z[i][j] = z0;
        idGrid[i][j] = Vd.length;
        Vd.push([xr, yr, z0]);
      }
    }
    const Fd = [];
    for (let i = 0; i < nX - 1; i++) {
      for (let j = 0; j < nY - 1; j++) {
        const p1 = idGrid[i][j], p2 = idGrid[i + 1][j], p3 = idGrid[i + 1][j + 1], p4 = idGrid[i][j + 1];
        Fd.push([p1, p3, p2], [p1, p4, p3]);
      }
    }
    return phys.buildCanonicalTargetModel(Vd, Fd, { X, Y, Z }, bD);
  }

  // F2: 0-deg aligned 90-deg dihedral: S_HH^MB == -S_VV^MB, OP matches exact, HV=0, PP rejects 2-bounce
  const fD = 10.0, aD = 0.5, bD = 1.0;
  const lamD = (phys.C0 * 1e-9) / fD;
  const exactDihDB = 10 * Math.log10((8 * Math.PI * aD * aD * bD * bD) / (lamD * lamD));
  const dih0 = buildRotatedDihedral(0.0);
  const d0HH = phys.evaluateRCSDirection(dih0, fD, 0, -90, true, "HH");
  const d0VV = phys.evaluateRCSDirection(dih0, fD, 0, -90, true, "VV");
  const d0HV = phys.evaluateRCSDirection(dih0, fD, 0, -90, true, "HV");
  const d0PP = phys.evaluateRCSDirection(dih0, fD, 0, -90, true, "PP");
  const d0OP = phys.evaluateRCSDirection(dih0, fD, 0, -90, true, "OP");

  const oppPhaseResidual =
    Math.hypot(d0HH.reMB + d0VV.reMB, d0HH.imMB + d0VV.imMB) / Math.hypot(d0HH.reMB, d0HH.imMB);
  const ppMBResidual = Math.hypot(d0PP.reMB, d0PP.imMB) / Math.hypot(d0OP.reMB, d0OP.imMB);
  const ppRejectionDB = d0OP.total_dBsm - d0PP.total_dBsm;
  record(
    "F2 0deg dihedral: S_HH=-S_VV, OP=exact, PP cancels 2-bounce",
    Math.abs(d0OP.total_dBsm - exactDihDB) < 0.2 &&
      oppPhaseResidual < 1e-12 &&
      ppMBResidual < 1e-12 &&
      d0HV.channels.HV.m2 < 1e-20 &&
      ppRejectionDB > 35.0,
    `OP=${d0OP.total_dBsm.toFixed(2)} dBsm (exact ${exactDihDB.toFixed(2)}), |S_PP^MB|/|S_OP^MB|=${ppMBResidual.toExponential(2)}, OP-PP=${ppRejectionDB.toFixed(1)} dB`
  );

  // F3: 45-deg canted 90-deg dihedral: 100% linear depolarization into HV/VH & roll-invariant OP
  const dih45 = buildRotatedDihedral(45.0);
  const d45HH = phys.evaluateRCSDirection(dih45, fD, 0, -90, true, "HH");
  const d45VV = phys.evaluateRCSDirection(dih45, fD, 0, -90, true, "VV");
  const d45HV = phys.evaluateRCSDirection(dih45, fD, 0, -90, true, "HV");
  const d45VH = phys.evaluateRCSDirection(dih45, fD, 0, -90, true, "VH");
  const d45OP = phys.evaluateRCSDirection(dih45, fD, 0, -90, true, "OP");

  const hvErrDB = Math.abs(d45HV.total_dBsm - exactDihDB);
  const coMBResidual =
    Math.max(Math.hypot(d45HH.reMB, d45HH.imMB), Math.hypot(d45VV.reMB, d45VV.imMB)) /
    Math.hypot(d45HV.reMB, d45HV.imMB);
  const coSuppressionDB = d45HV.total_dBsm - d45HH.total_dBsm;
  const opRollDiffDB = Math.abs(d45OP.total_dBsm - d0OP.total_dBsm);
  record(
    "F3 45deg canted dihedral: pure HV/VH & roll-invariant OP",
    hvErrDB < 0.2 &&
      Math.abs(d45HV.total_dBsm - d45VH.total_dBsm) < 1e-13 &&
      coMBResidual < 1e-12 &&
      coSuppressionDB > 35.0 &&
      opRollDiffDB < 1e-10,
    `HV=${d45HV.total_dBsm.toFixed(2)} dBsm, |S_HH^MB|/|S_HV^MB|=${coMBResidual.toExponential(2)}, |OP(45)-OP(0)|=${opRollDiffDB.toExponential(2)} dB`
  );

  // F4: 3D waverider Lorentz reciprocity (S_HV == S_VH), Parseval Span identity & mirror symmetry
  const mPol = phys.buildWaveriderModel(PRESETS.booster_mounted);
  let maxRecipErr = 0.0;
  let maxSpanRelErr = 0.0;
  for (const az of [-55, -25, 0, 25, 45, 55]) {
    for (const f of [1.5, 10.0]) {
      const res = phys.evaluateRCSDirection(mPol, f, az, -22.0, true, "HH");
      const ch = res.channels;
      maxRecipErr = Math.max(maxRecipErr, Math.hypot(ch.HV.re - ch.VH.re, ch.HV.im - ch.VH.im));
      const circPower = ch.PP.m2 + 0.5 * (ch.OP.m2 + ch.LL.m2);
      const linHalfSpan = 0.5 * ch.span_m2;
      maxSpanRelErr = Math.max(maxSpanRelErr, Math.abs(circPower - linHalfSpan) / Math.max(linHalfSpan, 1e-15));
    }
  }
  const rPos = phys.evaluateRCSDirection(mPol, 10.0, +45.0, -22.0, true, "OP");
  const rNeg = phys.evaluateRCSDirection(mPol, 10.0, -45.0, -22.0, true, "LL");
  const mirrorHandErr = Math.max(
    Math.abs(rPos.channels.HH.dBsm - rNeg.channels.HH.dBsm),
    Math.abs(rPos.channels.HV.dBsm - rNeg.channels.HV.dBsm),
    Math.abs(rPos.channels.OP.dBsm - rNeg.channels.LL.dBsm)
  );
  record(
    "F4 3D waverider reciprocity, Parseval span & mirror symmetry",
    maxRecipErr < 1e-14 && maxSpanRelErr < 1e-13 && mirrorHandErr < 1e-9,
    `|S_HV-S_VH|=${maxRecipErr.toExponential(2)}, span err=${maxSpanRelErr.toExponential(2)}, mirror=${mirrorHandErr.toExponential(2)} dB`
  );

  // F5: Polarimetric HRRP/RTI: OP & HV reject 100% of 1-bounce clutter while isolating delayed 2-bounce
  const mCav = phys.buildWaveriderModel({ ...PRESETS.caret, anhedral: 20.0 });
  const rtiHH = phys.computeRangeAngleIntensity(mCav, 10.0, 1000.0, -45.0, true, "hann", {
    azVec: [75, 84],
    polarization: "HH"
  });
  const rtiOP = phys.computeRangeAngleIntensity(mCav, 10.0, 1000.0, -45.0, true, "hann", {
    azVec: [75, 84],
    polarization: "OP"
  });
  const rtiHV = phys.computeRangeAngleIntensity(mCav, 10.0, 1000.0, -45.0, true, "hann", {
    azVec: [75, 84],
    polarization: "HV"
  });
  const peakSngHH = Math.max(...rtiHH.RTI_single_dBsm.flatMap((r) => Array.from(r)));
  const maxSngOPFloorDev = Math.max(...rtiOP.RTI_single_dBsm.flatMap((r) => Array.from(r)).map((v) => Math.abs(v - -50.0)));
  const maxSngHVFloorDev = Math.max(...rtiHV.RTI_single_dBsm.flatMap((r) => Array.from(r)).map((v) => Math.abs(v - -50.0)));
  const peakDblOP = Math.max(...rtiOP.RTI_double_dBsm.flatMap((r) => Array.from(r)));
  const opTotVsDblDiff = Math.max(
    ...rtiOP.RTI_total_dBsm.flatMap((r, ia) => Array.from(r).map((v, ir) => Math.abs(v - rtiOP.RTI_double_dBsm[ia][ir])))
  );
  record(
    "F5 Polarimetric HRRP/RTI: OP & HV isolate delayed 2-bounce",
    peakSngHH > -25.0 &&
      maxSngOPFloorDev < 1e-12 &&
      maxSngHVFloorDev < 1e-12 &&
      peakDblOP > -20.0 &&
      opTotVsDblDiff < 1e-12,
    `HH 1-bounce=${peakSngHH.toFixed(1)} dBsm -> OP/HV 1-bounce=-50.0 dBsm, OP 2-bounce=${peakDblOP.toFixed(1)} dBsm`
  );
}

// ---------------------------------------------------------------------------
// G. CONVENTIONAL AXISYMMETRIC REENTRY VEHICLES (RVs): TOPOLOGY, AERO & RCS
// ---------------------------------------------------------------------------
{
  // G1: 19-configuration RV topology & winding audit (3 body x 3 nose x 2 aft + booster-mounted RV)
  const rvBodies = ["mono_conic", "bi_conic", "tri_conic"];
  const rvNoses = ["conical", "spherical", "ogive"];
  const rvAfts = ["flat", "concave"];
  const rvConfigs = [];
  for (const rvBodyType of rvBodies) {
    for (const rvNoseType of rvNoses) {
      for (const rvAftType of rvAfts) {
        rvConfigs.push({
          label: `${rvBodyType}/${rvNoseType}/${rvAftType}`,
          cfg: {
            method: "rv",
            mach: 7.0,
            length: 3.2,
            rvBodyType,
            rvNoseType,
            rvBluntness: rvNoseType === "conical" ? 0.0 : 0.2,
            rvConeAngles: rvBodyType === "tri_conic" ? [12.0, 6.5, 13.0] : [11.5, 6.5, 10.0],
            rvStageFractions: [0.38, 0.34, 0.28],
            rvAftType,
            rvAftDepth: 0.35,
            rvAftPower: 2.2
          }
        });
      }
    }
  }
  rvConfigs.push({
    label: "bi_conic/spherical/booster_mounted",
    cfg: {
      method: "rv",
      mach: 7.0,
      length: 3.2,
      rvBodyType: "bi_conic",
      rvNoseType: "spherical",
      rvBluntness: 0.18,
      rvConeAngles: [11.5, 6.0, 10.0],
      rvStageFractions: [0.42, 0.58, 0.0],
      rvAftType: "flat",
      frustum: true,
      flength: 1.1,
      bradius: 0.55,
      blend: "cubic",
      rocket: true,
      rocketLength: 1.4,
      rcone: 3.5
    }
  });

  const rvFails = [];
  let maxWindingErr = 0.0;
  for (const item of rvConfigs) {
    const m = phys.buildWaveriderModel(item.cfg);
    const wIn = phys.generalizedWindingNumber(m.vertices, m.faces, [0.52 * item.cfg.length, 0.0, 0.0]);
    const wOut = phys.generalizedWindingNumber(m.vertices, m.faces, [-0.2 * item.cfg.length, 0.0, 0.0]);
    let wCav = 0.0;
    if (m.rvMeta && m.rvMeta.aftCavityDepth_m > 1e-5 && !item.cfg.frustum) {
      wCav = phys.generalizedWindingNumber(
        m.vertices,
        m.faces,
        [item.cfg.length - 0.5 * m.rvMeta.aftCavityDepth_m, 0.0, 0.0]
      );
    }
    const wErr = Math.max(Math.abs(wIn - 1.0), Math.abs(wOut), Math.abs(wCav));
    maxWindingErr = Math.max(maxWindingErr, wErr);
    if (
      !m.isWatertight ||
      m.euler !== 2 ||
      m.boundaryEdges !== 0 ||
      m.nonManifoldEdges !== 0 ||
      m.conflictingHalfEdges !== 0 ||
      m.degenerateFaces !== 0 ||
      wErr > 1e-6
    ) {
      rvFails.push(`${item.label} (euler=${m.euler}, wErr=${wErr.toExponential(2)})`);
    }
  }
  record(
    "G1 RV sweep: all 19 configs watertight & cavity exterior",
    rvFails.length === 0,
    `${rvConfigs.length} RV configs, ${rvFails.length} issues, max|dW|=${maxWindingErr.toExponential(2)}`
  );

  // G2: Analytical right-circular cone & recessed paraboloid cap volumetric convergence
  const Lcone = 3.0;
  const thConeDeg = 10.0;
  const mConeFlat = phys.buildWaveriderModel({
    method: "rv",
    mach: 7.0,
    length: Lcone,
    rvBodyType: "mono_conic",
    rvNoseType: "conical",
    rvConeAngles: [thConeDeg, thConeDeg, thConeDeg],
    rvAftType: "flat",
    numStreamwise: 48,
    numHalfSpan: 73,
    numBaseLayers: 16
  });
  const mConeConc = phys.buildWaveriderModel({
    method: "rv",
    mach: 7.0,
    length: Lcone,
    rvBodyType: "mono_conic",
    rvNoseType: "conical",
    rvConeAngles: [thConeDeg, thConeDeg, thConeDeg],
    rvAftType: "concave",
    rvAftDepth: 0.4,
    rvAftPower: 2.0,
    numStreamwise: 48,
    numHalfSpan: 73,
    numBaseLayers: 16
  });
  const RbCone = Lcone * Math.tan((thConeDeg * Math.PI) / 180.0);
  const VconeExact = (Math.PI * RbCone * RbCone * Lcone) / 3.0;
  const dAftExact = mConeConc.rvMeta.aftCavityDepth_m;
  const dVparabExact = 0.5 * Math.PI * RbCone * RbCone * dAftExact;
  const coneVolRelErr = Math.abs(mConeFlat.enclosedVol - VconeExact) / VconeExact;
  const cavVolRelErr = Math.abs((mConeFlat.enclosedVol - mConeConc.enclosedVol) - dVparabExact) / dVparabExact;
  record(
    "G2 RV exact cone & paraboloid aft-cavity volumes",
    coneVolRelErr < 2e-3 && cavVolRelErr < 1e-2,
    `cone err=${(100 * coneVolRelErr).toFixed(3)}%, paraboloid cap err=${(100 * cavVolRelErr).toFixed(3)}%`
  );

  // G3: Exact Taylor-Maccoll sharp-cone wave drag & RV zero-lift axisymmetry
  const cpConeExact = mConeFlat.flowSolution.CpCone;
  const cdWaveBaseRef = (mConeFlat.CD_wave * mConeFlat.planformArea) / (Math.PI * RbCone * RbCone);
  const cdConeRelErr = Math.abs(cdWaveBaseRef - cpConeExact) / cpConeExact;

  const rvPresetConfigs = [
    { method: "rv", mach: 7.0, length: 3.1, rvBodyType: "mono_conic", rvNoseType: "spherical", rvBluntness: 0.16, rvConeAngles: [9.0, 9.0, 9.0], rvAftType: "flat" },
    { method: "rv", mach: 7.0, length: 3.4, rvBodyType: "bi_conic", rvNoseType: "spherical", rvBluntness: 0.18, rvConeAngles: [11.5, 6.0, 6.0], rvStageFractions: [0.42, 0.58, 0.0], rvAftType: "flat" },
    { method: "rv", mach: 7.0, length: 3.5, rvBodyType: "tri_conic", rvNoseType: "ogive", rvBluntness: 0.25, rvConeAngles: [12.0, 6.5, 12.5], rvStageFractions: [0.35, 0.4, 0.25], rvAftType: "concave", rvAftDepth: 0.35, rvAftPower: 2.0 },
    { method: "rv", mach: 7.0, length: 3.0, rvBodyType: "mono_conic", rvNoseType: "conical", rvConeAngles: [8.0, 8.0, 8.0], rvAftType: "flat" }
  ];
  let maxZeroAero = 0.0;
  let maxAntisymErr = 0.0;
  let minPosLift = Infinity;
  for (const cfgRV of rvPresetConfigs) {
    const mRV = phys.buildWaveriderModel(cfgRV);
    maxZeroAero = Math.max(maxZeroAero, Math.abs(mRV.CL), Math.abs(mRV.CY), Math.abs(mRV.Cm));
    const pPos5 = phys.evaluateAeroState(mRV, +5.0);
    const pNeg5 = phys.evaluateAeroState(mRV, -5.0);
    minPosLift = Math.min(minPosLift, pPos5.CL);
    maxAntisymErr = Math.max(maxAntisymErr, Math.abs(pPos5.CL + pNeg5.CL), Math.abs(pPos5.CD_wave - pNeg5.CD_wave));
  }
  record(
    "G3 RV Taylor-Maccoll cone drag & zero-lift axisymmetry",
    cdConeRelErr < 2e-3 && maxZeroAero < 1e-13 && minPosLift > 0.02 && maxAntisymErr < 1e-13,
    `TM CD err=${(100 * cdConeRelErr).toFixed(3)}%, max|CL0,CY0,Cm0|=${maxZeroAero.toExponential(2)}, antisym=${maxAntisymErr.toExponential(2)}`
  );

  // G4: Axisymmetric RV radar symmetries & flat aft circular-disk tail-on PO RCS
  const rAzPos = phys.evaluateRCSDirection(mConeFlat, 10.0, +42.0, -18.0, true, "HH").total_dBsm;
  const rAzNeg = phys.evaluateRCSDirection(mConeFlat, 10.0, -42.0, -18.0, true, "HH").total_dBsm;
  const rElPos = phys.evaluateRCSDirection(mConeFlat, 10.0, 35.0, +24.0, true, "VV").total_dBsm;
  const rElNeg = phys.evaluateRCSDirection(mConeFlat, 10.0, 35.0, -24.0, true, "VV").total_dBsm;
  const symRadErr = Math.max(Math.abs(rAzPos - rAzNeg), Math.abs(rElPos - rElNeg));

  const lam10 = (phys.C0 * 1e-9) / 10.0;
  const Abase = Math.PI * RbCone * RbCone;
  const sigDiskExactDB = 10 * Math.log10((4 * Math.PI * Abase * Abase) / (lam10 * lam10));
  const tailHH = phys.evaluateRCSDirection(mConeFlat, 10.0, 180.0, 0.0, true, "HH");
  const tailVV = phys.evaluateRCSDirection(mConeFlat, 10.0, 180.0, 0.0, true, "VV");
  const tailPP = phys.evaluateRCSDirection(mConeFlat, 10.0, 180.0, 0.0, true, "PP");
  const tailOP = phys.evaluateRCSDirection(mConeFlat, 10.0, 180.0, 0.0, true, "OP");
  const diskErrDB = Math.max(
    Math.abs(tailHH.total_dBsm - sigDiskExactDB),
    Math.abs(tailVV.total_dBsm - sigDiskExactDB),
    Math.abs(tailPP.total_dBsm - sigDiskExactDB)
  );
  record(
    "G4 RV radar symmetry & flat-aft circular disk tail-on RCS",
    symRadErr < 1e-9 && diskErrDB < 0.05 && Math.abs(tailOP.total_dBsm - -50.0) < 1e-12,
    `tail HH=${tailHH.total_dBsm.toFixed(2)} dBsm vs disk=${sigDiskExactDB.toFixed(2)} dBsm (err ${diskErrDB.toFixed(3)} dB), sym=${symRadErr.toExponential(2)} dB`
  );

  // G5: Flat vs. concave RV aft plane off-tail multi-bounce OP & delayed HRRP cavity peak
  const rvBaseCfg = {
    method: "rv",
    mach: 7.0,
    length: 3.2,
    rvBodyType: "bi_conic",
    rvNoseType: "spherical",
    rvBluntness: 0.18,
    rvConeAngles: [11.5, 6.0, 6.0],
    rvStageFractions: [0.42, 0.58, 0.0]
  };
  const mAftFlat = phys.buildWaveriderModel({ ...rvBaseCfg, rvAftType: "flat" });
  const mAftConc = phys.buildWaveriderModel({ ...rvBaseCfg, rvAftType: "concave", rvAftDepth: 0.45, rvAftPower: 2.6 });
  const flatOffOP = phys.evaluateRCSDirection(mAftFlat, 10.0, 168.0, 0.0, true, "OP");
  const concOffOP = phys.evaluateRCSDirection(mAftConc, 10.0, 168.0, 0.0, true, "OP");
  const rtiFlatAft = phys.computeRangeAngleIntensity(mAftFlat, 10.0, 2000.0, 0.0, true, "hann", {
    azVec: [168.0],
    polarization: "HH"
  });
  const rtiConcAft = phys.computeRangeAngleIntensity(mAftConc, 10.0, 2000.0, 0.0, true, "hann", {
    azVec: [168.0],
    polarization: "HH"
  });

  let idxSng = 0, valSng = -Infinity, idxDbl = 0, valDbl = -Infinity;
  for (let i = 0; i < rtiConcAft.rangeVec.length; i++) {
    if (rtiConcAft.RTI_single_dBsm[0][i] > valSng) {
      valSng = rtiConcAft.RTI_single_dBsm[0][i];
      idxSng = i;
    }
    if (rtiConcAft.RTI_double_dBsm[0][i] > valDbl) {
      valDbl = rtiConcAft.RTI_double_dBsm[0][i];
      idxDbl = i;
    }
  }
  const rSngPk = rtiConcAft.rangeVec[idxSng];
  const rDblPk = rtiConcAft.rangeVec[idxDbl];
  const rBaseRim = rtiConcAft.R_base[0];
  const flatMaxDbl = Math.max(...rtiFlatAft.RTI_double_dBsm[0]);
  const cavDelayM = rDblPk - rBaseRim;
  record(
    "G5 Concave RV aft cavity: OP retroreflection & delayed HRRP",
    Math.abs(flatOffOP.total_dBsm - -50.0) < 1e-12 &&
      Math.abs(flatMaxDbl - -50.0) < 1e-12 &&
      concOffOP.total_dBsm > 10.0 &&
      valDbl > 10.0 &&
      cavDelayM > 0.2 &&
      rDblPk > rSngPk + 0.15,
    `flat OP=${flatOffOP.total_dBsm.toFixed(1)} dBsm -> concave OP=${concOffOP.total_dBsm.toFixed(2)} dBsm, cavity delay +${(100 * cavDelayM).toFixed(1)} cm`
  );
}

const nPass = results.filter((r) => r.pass).length;
console.log("=====================================================================");
console.log(`  ${nPass} / ${results.length} WEB STUDIO VALIDATION CHECKS PASSED`);
console.log("=====================================================================\n");
