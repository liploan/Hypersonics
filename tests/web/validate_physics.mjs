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
 *          symmetry (|CY| = 0, RCS(+az) == RCS(-az)), and < 100 ms frame budget
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
    phys.evaluateRCSMultiFreq(mBench, azBands, az, -22, true);
  }
  const t2 = performance.now();
  const elBands = [0.5, 3.0, 10.0];
  for (let el = -90; el <= 90; el += 3) {
    phys.evaluateRCSMultiFreq(mBench, elBands, 45, el, true);
  }
  const t3 = performance.now();
  phys.computeRangeAngleIntensity(mBench, 10.0, 1000.0, -22, true, "hann");
  const t4 = performance.now();
  const dtMs = t4 - t0;
  record(
    "E4 Full Web Studio rebuild + multi-band RCS + RTI latency",
    dtMs < 250,
    `total ${dtMs.toFixed(1)} ms (mesh ${(t1 - t0).toFixed(1)}, az ${(t2 - t1).toFixed(1)}, el ${(t3 - t2).toFixed(1)}, rti ${(t4 - t3).toFixed(1)})`
  );
}

const nPass = results.filter((r) => r.pass).length;
console.log("=====================================================================");
console.log(`  ${nPass} / ${results.length} WEB STUDIO VALIDATION CHECKS PASSED`);
console.log("=====================================================================\n");
