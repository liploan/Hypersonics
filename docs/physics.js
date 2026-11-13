/* ============================================================================
 * Hypersonic Waverider & UHF-to-X-Band RCS Physics Engine (docs/physics.js)
 *
 * Exact mathematical parity with the verified MATLAB package (+waverider/):
 *   - Flow Solvers: 2D oblique shock & 3D axisymmetric Taylor-Maccoll ODE
 *   - Waverider Lofting: Cone-derived (Rasmussen/Jones), Sobieczky inverse
 *     osculating-cones with non-folding evolute guard, Nonweiler planar caret,
 *     and Kulfan CST boost-glide bodies + conformal booster frustum
 *   - Topology: Mirror-symmetric 2-manifold watertight stitcher & half-edge,
 *     Euler, degenerate-area, surface-fold, and generalized-winding verifier
 *   - Aerodynamics: Calibrated embedded-Newtonian + Prandtl-Meyer expansion
 *     and Eckert (1955) 30 km reference-temperature turbulent skin friction
 *   - Electromagnetics: Exact Gordon (1975) triangle contour Physical Optics
 *     integral + Sutherland-Hodgman GO-PO concave double-bounce cavity solver
 *   - Radar Signal Processing: Stepped-frequency / matched-filter RTI & HRRP
 *     with analytical Harris spectral window PSFs (rect, hann, hamming)
 * ============================================================================ */

(function (root, factory) {
  if (typeof module === "object" && module.exports) {
    module.exports = factory();
  } else {
    root.WaveriderPhysics = factory();
  }
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  const C0 = 299792458.0; // Speed of light in vacuum [m/s]
  const DEG2RAD = Math.PI / 180.0;
  const RAD2DEG = 180.0 / Math.PI;

  // --------------------------------------------------------------------------
  // 1. Compressible & Conical Flow Solvers
  // --------------------------------------------------------------------------

  function solveObliqueShock(M1, betaDeg, gamma = 1.4) {
    if (!(M1 > 1.0)) {
      throw new RangeError(`Freestream Mach M1 = ${M1} must exceed 1.0.`);
    }
    const muRad = Math.asin(1.0 / M1);
    const muDeg = muRad * RAD2DEG;
    if (!(betaDeg > muDeg && betaDeg < 90.0)) {
      throw new RangeError(
        `Shock angle beta = ${betaDeg.toFixed(2)} deg must exceed Mach angle mu = ${muDeg.toFixed(2)} deg (M = ${M1.toFixed(2)}).`
      );
    }
    const betaRad = betaDeg * DEG2RAD;
    const M1n = M1 * Math.sin(betaRad);
    const cotBeta = 1.0 / Math.tan(betaRad);
    const num = 2.0 * cotBeta * (M1n * M1n - 1.0);
    const den = M1 * M1 * (gamma + Math.cos(2.0 * betaRad)) + 2.0;
    const deltaRad = Math.atan(num / den);
    const deltaDeg = deltaRad * RAD2DEG;

    const M2n_sq = (1.0 + 0.5 * (gamma - 1.0) * M1n * M1n) / (gamma * M1n * M1n - 0.5 * (gamma - 1.0));
    const M2n = Math.sqrt(Math.max(M2n_sq, 1e-10));
    const M2 = M2n / Math.sin(Math.max(betaRad - deltaRad, 1e-8));

    const p2_p1 = 1.0 + ((2.0 * gamma) / (gamma + 1.0)) * (M1n * M1n - 1.0);
    const rho2_rho1 = ((gamma + 1.0) * M1n * M1n) / (2.0 + (gamma - 1.0) * M1n * M1n);
    const T2_T1 = p2_p1 / rho2_rho1;
    const V1_bar = Math.pow(1.0 + 2.0 / ((gamma - 1.0) * M1 * M1), -0.5);
    const V2_bar = Math.pow(1.0 + 2.0 / ((gamma - 1.0) * M2 * M2), -0.5);
    const Cp2 = (p2_p1 - 1.0) / (0.5 * gamma * M1 * M1);

    return {
      M1, gamma, betaDeg, betaRad, deltaDeg, deltaRad, muDeg,
      M1n, M2, p2_p1, rho2_rho1, T2_T1, V1_bar, V2_bar, Cp2
    };
  }

  function solveTaylorMaccoll(M1, betaDeg, gamma = 1.4) {
    const shock = solveObliqueShock(M1, betaDeg, gamma);
    const betaRad = shock.betaRad;
    let ur = shock.V2_bar * Math.cos(betaRad - shock.deltaRad);
    let ut = -shock.V2_bar * Math.sin(betaRad - shock.deltaRad);
    let th = betaRad;
    const dth = -0.0002; // ~0.0115 deg step for high-accuracy RK4 integration

    const thetaList = [];
    const urList = [];
    const utList = [];
    const slopeList = [];
    const cpList = [];
    const machList = [];

    function pushState(theta, u_r, u_t) {
      const ux = u_r * Math.cos(theta) - u_t * Math.sin(theta);
      const ucyl = u_r * Math.sin(theta) + u_t * Math.cos(theta);
      const Vsq = u_r * u_r + u_t * u_t;
      const MlocSq = (2.0 / (gamma - 1.0)) * (Vsq / Math.max(1.0 - Vsq, 1e-12));
      const T0_T2 = 1.0 + 0.5 * (gamma - 1.0) * shock.M2 * shock.M2;
      const T0_Tloc = 1.0 + 0.5 * (gamma - 1.0) * MlocSq;
      const p_p1 = shock.p2_p1 * Math.pow(T0_T2 / T0_Tloc, gamma / (gamma - 1.0));

      thetaList.push(theta);
      urList.push(u_r);
      utList.push(u_t);
      slopeList.push(ucyl / Math.max(ux, 1e-12));
      cpList.push((p_p1 - 1.0) / (0.5 * gamma * M1 * M1));
      machList.push(Math.sqrt(Math.max(MlocSq, 0.0)));
    }

    function tmDeriv(theta, u_r, u_t) {
      const a_sq = 0.5 * (gamma - 1.0) * (1.0 - u_r * u_r - u_t * u_t);
      const num = a_sq * (2.0 * u_r + u_t / Math.tan(theta)) - u_r * u_t * u_t;
      let den = u_t * u_t - a_sq;
      if (Math.abs(den) < 1e-12) {
        den = (den >= 0 ? 1.0 : -1.0) * 1e-12;
      }
      return [u_t, num / den];
    }

    pushState(th, ur, ut);

    const maxSteps = Math.ceil(betaRad / Math.abs(dth)) + 5;
    for (let step = 0; step < maxSteps; step++) {
      if (th + dth <= 1e-4) break;
      const [k1r, k1t] = tmDeriv(th, ur, ut);
      const [k2r, k2t] = tmDeriv(th + 0.5 * dth, ur + 0.5 * dth * k1r, ut + 0.5 * dth * k1t);
      const [k3r, k3t] = tmDeriv(th + 0.5 * dth, ur + 0.5 * dth * k2r, ut + 0.5 * dth * k2t);
      const [k4r, k4t] = tmDeriv(th + dth, ur + dth * k3r, ut + dth * k3t);

      const urNext = ur + (dth / 6.0) * (k1r + 2.0 * k2r + 2.0 * k3r + k4r);
      const utNext = ut + (dth / 6.0) * (k1t + 2.0 * k2t + 2.0 * k3t + k4t);

      if (utNext >= 0.0) {
        // Refine u_theta(theta) = 0 root using cubic Hermite interpolation on [th, th + dth]
        let s = -ut / Math.max(utNext - ut, 1e-15);
        for (let nr = 0; nr < 4; nr++) {
          const h00 = (1 + 2 * s) * (1 - s) * (1 - s);
          const h10 = s * (1 - s) * (1 - s);
          const h01 = s * s * (3 - 2 * s);
          const h11 = s * s * (s - 1);
          const fVal = h00 * ut + h10 * dth * k1t + h01 * utNext + h11 * dth * k4t;
          const dfds =
            6 * s * (s - 1) * (ut - utNext) +
            (1 - 4 * s + 3 * s * s) * dth * k1t +
            s * (3 * s - 2) * dth * k4t;
          s = Math.max(0.0, Math.min(1.0, s - fVal / (Math.abs(dfds) > 1e-15 ? dfds : (utNext - ut))));
        }
        const h00 = (1 + 2 * s) * (1 - s) * (1 - s);
        const h10 = s * (1 - s) * (1 - s);
        const h01 = s * s * (3 - 2 * s);
        const h11 = s * s * (s - 1);
        th = th + s * dth;
        ur = h00 * ur + h10 * dth * k1r + h01 * urNext + h11 * dth * k4r;
        ut = 0.0;
        pushState(th, ur, ut);
        break;
      }

      th += dth;
      ur = urNext;
      ut = utNext;
      pushState(th, ur, ut);
    }

    thetaList.reverse();
    urList.reverse();
    utList.reverse();
    slopeList.reverse();
    cpList.reverse();
    machList.reverse();

    const nPts = thetaList.length;
    const thetaRad = Float64Array.from(thetaList);
    const slopeRatio = Float64Array.from(slopeList);
    const Cp = Float64Array.from(cpList);
    const Mach = Float64Array.from(machList);

    function interpTM(thetaQuery) {
      if (thetaQuery <= thetaRad[0]) {
        return { slope: slopeRatio[0], cp: Cp[0] };
      }
      if (thetaQuery >= thetaRad[nPts - 1]) {
        return { slope: slopeRatio[nPts - 1], cp: Cp[nPts - 1] };
      }
      let lo = 0;
      let hi = nPts - 1;
      while (hi - lo > 1) {
        const mid = (lo + hi) >> 1;
        if (thetaRad[mid] <= thetaQuery) lo = mid;
        else hi = mid;
      }
      const w = (thetaQuery - thetaRad[lo]) / Math.max(thetaRad[hi] - thetaRad[lo], 1e-15);
      return {
        slope: slopeRatio[lo] * (1.0 - w) + slopeRatio[hi] * w,
        cp: Cp[lo] * (1.0 - w) + Cp[hi] * w
      };
    }

    return {
      M1, gamma, betaDeg, betaRad,
      thetaConeRad: thetaRad[0],
      thetaConeDeg: thetaRad[0] * RAD2DEG,
      CpCone: Cp[0],
      MachCone: Mach[0],
      shock,
      thetaRad,
      slopeRatio,
      Cp,
      Mach,
      interpTM
    };
  }

  function traceStreamline(tm, xStartLocal, xEndLocal, Nx) {
    const nCols = xStartLocal.length;
    const xGrid = Array.from({ length: Nx }, () => new Float64Array(nCols));
    const rGrid = Array.from({ length: Nx }, () => new Float64Array(nCols));
    const thetaGrid = Array.from({ length: Nx }, () => new Float64Array(nCols));
    const cpGrid = Array.from({ length: Nx }, () => new Float64Array(nCols));

    const tanBeta = Math.tan(tm.betaRad);
    const thMin = tm.thetaConeRad;
    const thMax = tm.betaRad;
    const numSubSteps = 4;

    function clampTh(v) {
      return v < thMin ? thMin : v > thMax ? thMax : v;
    }

    for (let c = 0; c < nCols; c++) {
      const xs = xStartLocal[c];
      const xe = xEndLocal[c];
      const dxTot = Math.max(xe - xs, 0.0);
      let rCurr = xs * tanBeta;

      xGrid[0][c] = xs;
      rGrid[0][c] = rCurr;
      thetaGrid[0][c] = thMax;
      cpGrid[0][c] = tm.interpTM(thMax).cp;

      const dxStep = Nx > 1 ? dxTot / ((Nx - 1) * numSubSteps) : 0.0;
      for (let i = 0; i < Nx - 1; i++) {
        let xSub = xs + (i / (Nx - 1)) * dxTot;
        let rSub = rCurr;
        if (dxTot > 1e-12) {
          for (let s = 0; s < numSubSteps; s++) {
            const th1 = clampTh(Math.atan2(rSub, Math.max(xSub, 1e-9)));
            const k1 = tm.interpTM(th1).slope;

            const x2 = xSub + 0.5 * dxStep;
            const r2 = rSub + 0.5 * dxStep * k1;
            const th2 = clampTh(Math.atan2(r2, Math.max(x2, 1e-9)));
            const k2 = tm.interpTM(th2).slope;

            const r3 = rSub + 0.5 * dxStep * k2;
            const th3 = clampTh(Math.atan2(r3, Math.max(x2, 1e-9)));
            const k3 = tm.interpTM(th3).slope;

            const x4 = xSub + dxStep;
            const r4 = rSub + dxStep * k3;
            const th4 = clampTh(Math.atan2(r4, Math.max(x4, 1e-9)));
            const k4 = tm.interpTM(th4).slope;

            rSub += (dxStep / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            xSub = x4;
          }
        }
        rCurr = rSub;
        const xNext = xs + ((i + 1) / (Nx - 1)) * dxTot;
        const thNext = clampTh(Math.atan2(rCurr, Math.max(xNext, 1e-9)));
        xGrid[i + 1][c] = xNext;
        rGrid[i + 1][c] = rCurr;
        thetaGrid[i + 1][c] = thNext;
        cpGrid[i + 1][c] = tm.interpTM(thNext).cp;
      }
    }

    return { xGrid, rGrid, thetaGrid, cpGrid };
  }

  // --------------------------------------------------------------------------
  // 2. Configuration Normalization & Validation
  // --------------------------------------------------------------------------

  function normalizeConfig(raw = {}) {
    const method = String(raw.method ?? raw.Method ?? "osculating_cone").toLowerCase();
    const mach = Number(raw.mach ?? raw.Mach ?? 7.0);
    const beta = Number(raw.beta ?? raw.ShockAngleDeg ?? 13.5);
    const gamma = Number(raw.gamma ?? raw.Gamma ?? 1.4);
    const length = Number(raw.length ?? raw.Length ?? 5.5);
    const span = Number(raw.span ?? raw.HalfSpan ?? 1.75);
    const power = Number(raw.power ?? raw.PlanformPower ?? 0.62);
    const anhedral = Number(raw.anhedral ?? raw.LeadingEdgeAnhedralDeg ?? 8.0);
    const canopy = Number(raw.canopy ?? raw.CanopyHeightRatio ?? 0.045);
    const upperSurfaceType = String(
      raw.upperSurfaceType ?? raw.UpperSurfaceType ?? (canopy > 0 ? "canopy" : "freestream")
    ).toLowerCase();
    const canopySpanPower = Number(raw.canopySpanPower ?? raw.CanopySpanPower ?? 2.2);

    const flatShockFraction = Number(raw.flatShockFraction ?? raw.FlatShockFraction ?? 0.25);
    const shockCurvaturePower = Number(raw.shockCurvaturePower ?? raw.ShockCurvaturePower ?? 2.2);
    const vertexOffsetRatio = Number(raw.vertexOffsetRatio ?? raw.VertexOffsetRatio ?? 0.35);

    const cstUpperWeights = Array.isArray(raw.cstUpperWeights ?? raw.CSTUpperWeights)
      ? (raw.cstUpperWeights ?? raw.CSTUpperWeights).map(Number)
      : [0.16, 0.22, 0.18, 0.12];
    const cstLowerArchRatio = Number(raw.cstLowerArchRatio ?? raw.CSTLowerArchRatio ?? 0.025);
    const cstNoseExponent = Number(raw.cstNoseExponent ?? raw.CSTNoseExponent ?? 0.65);
    const cstAftExponent = Number(raw.cstAftExponent ?? raw.CSTAftExponent ?? 1.0);

    const frustum = Boolean(raw.frustum ?? raw.EnableAftFrustum ?? false);
    const flength = Number(raw.flength ?? raw.FrustumLength ?? 1.65);
    const bradius = Number(raw.bradius ?? raw.BoosterRadius ?? 0.72);
    const blend = String(raw.blend ?? raw.FrustumBlendLaw ?? "cubic").toLowerCase();
    const frustumZOffset = raw.frustumZOffset ?? raw.FrustumZOffset ?? NaN;
    const rocket = Boolean(raw.rocket ?? raw.EnableRocketStage ?? false);
    const rocketLength = Number(raw.rocketLength ?? raw.RocketStageLength ?? 1.85);
    const rcone = Number(raw.rcone ?? raw.RocketConeAngleDeg ?? 4.5);

    const numStreamwise = Math.max(8, Math.round(Number(raw.numStreamwise ?? raw.NumStreamwise ?? 24)));
    const numHalfSpan = Math.max(5, Math.round(Number(raw.numHalfSpan ?? raw.NumHalfSpan ?? 19)));
    const numBaseLayers = Math.max(2, Math.round(Number(raw.numBaseLayers ?? raw.NumBaseLayers ?? 6)));
    const numFrustumStations = Math.max(4, Math.round(Number(raw.numFrustumStations ?? raw.NumFrustumStations ?? 12)));
    const numRocketStations = Math.max(4, Math.round(Number(raw.numRocketStations ?? raw.NumRocketStations ?? 9)));

    // Conventional Reentry Vehicle (RV) parameters (used when method === 'rv')
    const rvBodyType = String(raw.rvBodyType ?? raw.RVBodyType ?? "bi_conic").toLowerCase().trim();
    const rvNoseType = String(raw.rvNoseType ?? raw.RVNoseType ?? "spherical").toLowerCase().trim();
    const rvBluntness = Number(raw.rvBluntness ?? raw.RVNoseBluntnessRatio ?? 0.18);
    const rvConeAngles = Array.isArray(raw.rvConeAngles ?? raw.RVConeAnglesDeg)
      ? (raw.rvConeAngles ?? raw.RVConeAnglesDeg).map(Number)
      : [
          Number(raw.rvCone1 ?? 11.0),
          Number(raw.rvCone2 ?? 6.5),
          Number(raw.rvCone3 ?? 10.0)
        ];
    const rvStageFractions = Array.isArray(raw.rvStageFractions ?? raw.RVStageFractions)
      ? (raw.rvStageFractions ?? raw.RVStageFractions).map(Number)
      : [
          Number(raw.rvFrac1 ?? 0.42),
          Number(raw.rvFrac2 ?? 0.33),
          Number(raw.rvFrac3 ?? Math.max(0.05, 1.0 - Number(raw.rvFrac1 ?? 0.42) - Number(raw.rvFrac2 ?? 0.33)))
        ];
    const rvAftType = String(raw.rvAftType ?? raw.RVAftType ?? "flat").toLowerCase().trim();
    const rvAftDepth = Number(raw.rvAftDepth ?? raw.RVAftCavityDepthRatio ?? 0.28);
    const rvAftPower = Number(raw.rvAftPower ?? raw.RVAftCavityPower ?? 2.0);

    return {
      method, mach, beta, gamma, length, span, power, anhedral,
      canopy: upperSurfaceType === "freestream" ? 0.0 : canopy,
      upperSurfaceType, canopySpanPower,
      flatShockFraction, shockCurvaturePower, vertexOffsetRatio,
      cstUpperWeights, cstLowerArchRatio, cstNoseExponent, cstAftExponent,
      rvBodyType, rvNoseType, rvBluntness, rvConeAngles, rvStageFractions,
      rvAftType, rvAftDepth, rvAftPower,
      frustum, flength, bradius, blend, frustumZOffset,
      rocket, rocketLength, rcone,
      numStreamwise, numHalfSpan, numBaseLayers, numFrustumStations, numRocketStations
    };
  }

  function validateConfig(cfg) {
    if (!(cfg.mach > 1.2 && Number.isFinite(cfg.mach))) {
      throw new RangeError(`Mach (${cfg.mach}) must exceed 1.2.`);
    }
    if (!(cfg.length > 0 && Number.isFinite(cfg.length))) {
      throw new RangeError("Length must be strictly positive.");
    }
    if (cfg.method === "rv") {
      if (!["mono_conic", "bi_conic", "tri_conic"].includes(cfg.rvBodyType)) {
        throw new RangeError(`Unsupported rvBodyType "${cfg.rvBodyType}".`);
      }
      if (!["conical", "spherical", "ogive"].includes(cfg.rvNoseType)) {
        throw new RangeError(`Unsupported rvNoseType "${cfg.rvNoseType}".`);
      }
      if (!["flat", "concave"].includes(cfg.rvAftType)) {
        throw new RangeError(`Unsupported rvAftType "${cfg.rvAftType}".`);
      }
      if (!(cfg.rvConeAngles[0] >= 1.0 && cfg.rvConeAngles[0] <= 55.0)) {
        throw new RangeError("Stage 1 RV cone half-angle must lie in [1, 55] deg.");
      }
      for (let s = 1; s < cfg.rvConeAngles.length; s++) {
        if (!(cfg.rvConeAngles[s] >= -15.0 && cfg.rvConeAngles[s] <= 55.0)) {
          throw new RangeError("Upper RV conic stage half-angles must lie in [-15, 55] deg.");
        }
      }
      return true;
    }
    const muDeg = Math.asin(1.0 / cfg.mach) * RAD2DEG;
    if (!(cfg.beta > muDeg + 0.15 && cfg.beta < 65.0)) {
      throw new RangeError(
        `ShockAngleDeg (${cfg.beta.toFixed(2)} deg) must exceed Mach angle (${muDeg.toFixed(2)} deg) by at least 0.15 deg.`
      );
    }
    if (!(cfg.span > 0)) {
      throw new RangeError("Length and HalfSpan must be strictly positive.");
    }
    if (!(cfg.power >= 0.25 && cfg.power <= 2.5)) {
      throw new RangeError("PlanformPower must be in [0.25, 2.5].");
    }
    return true;
  }

  // --------------------------------------------------------------------------
  // 3. Waverider Surface Generation (4 Families + Conformal Booster Frustum)
  // --------------------------------------------------------------------------

  function buildUpperSurface(X_low, Y_low, Z_low, cfg) {
    const Nx = X_low.length;
    const Ny = X_low[0].length;
    const xLE = X_low[0];
    const yLE = Y_low[0];
    const zLE = Z_low[0];

    const X_up = Array.from({ length: Nx }, (_, i) => Array.from(X_low[i]));
    const Y_up = Array.from({ length: Nx }, () => Array.from(yLE));
    const Z_up = Array.from({ length: Nx }, () => Array.from(zLE));

    if (cfg.upperSurfaceType !== "freestream" && cfg.canopy > 0) {
      const L = cfg.length;
      let bHalf = 0.0;
      for (let k = 0; k < Ny; k++) bHalf = Math.max(bHalf, Math.abs(yLE[k]));
      if (bHalf < 1e-9) bHalf = cfg.span;

      for (let k = 0; k < Ny; k++) {
        const dxMax = L - xLE[k];
        if (dxMax < 1e-7) continue;
        const etaSpan = Math.min(Math.abs(yLE[k]) / bHalf, 1.0);
        const spanProf = Math.pow(1.0 - Math.pow(etaSpan, cfg.canopySpanPower), 2.0);
        const chordScale = Math.pow(Math.max(1.0 - xLE[k] / L, 0.0), 0.75);

        for (let i = 0; i < Nx; i++) {
          const sNorm = Math.max(X_up[i][k] - xLE[k], 0.0) / dxMax;
          const longProf = Math.pow(sNorm, 0.85) * (1.0 - 0.25 * sNorm);
          Z_up[i][k] = zLE[k] + cfg.canopy * L * chordScale * spanProf * longProf;
        }
      }
    }

    // Exact coincidence along LE (row 0) and wingtips (col 0, Ny-1)
    for (let k = 0; k < Ny; k++) {
      X_up[0][k] = X_low[0][k];
      Y_up[0][k] = Y_low[0][k];
      Z_up[0][k] = Z_low[0][k];
    }
    for (let i = 0; i < Nx; i++) {
      X_up[i][0] = X_low[i][0]; Y_up[i][0] = Y_low[i][0]; Z_up[i][0] = Z_low[i][0];
      X_up[i][Ny - 1] = X_low[i][Ny - 1]; Y_up[i][Ny - 1] = Y_low[i][Ny - 1]; Z_up[i][Ny - 1] = Z_low[i][Ny - 1];
    }
    return { X_up, Y_up, Z_up };
  }

  function binomCoeff(n, r) {
    if (r < 0 || r > n) return 0;
    let res = 1;
    for (let i = 1; i <= r; i++) res = (res * (n - i + 1)) / i;
    return res;
  }

  function computeStructuredLowerNormals(X, Y, Z) {
    const Nx = X.length;
    const Ny = X[0].length;
    const Nx_l = Array.from({ length: Nx }, () => new Float64Array(Ny));
    const Ny_l = Array.from({ length: Nx }, () => new Float64Array(Ny));
    const Nz_l = Array.from({ length: Nx }, () => new Float64Array(Ny));

    for (let i = 0; i < Nx; i++) {
      const i0 = Math.max(0, i - 1);
      const i1 = Math.min(Nx - 1, i + 1);
      for (let k = 0; k < Ny; k++) {
        const k0 = Math.max(0, k - 1);
        const k1 = Math.min(Ny - 1, k + 1);
        // u along +X (row), v along +Y (col); v x u points downward (-Z)
        let ux = X[i1][k] - X[i0][k];
        let uy = Y[i1][k] - Y[i0][k];
        let uz = Z[i1][k] - Z[i0][k];
        if (ux * ux + uy * uy + uz * uz < 1e-18) {
          // Wingtip singularity: borrow neighbor column streamwise tangent
          const kn = k === 0 ? 1 : Ny - 2;
          ux = X[i1][kn] - X[i0][kn];
          uy = Y[i1][kn] - Y[i0][kn];
          uz = Z[i1][kn] - Z[i0][kn];
        }
        const vx = X[i][k1] - X[i][k0];
        const vy = Y[i][k1] - Y[i][k0];
        const vz = Z[i][k1] - Z[i][k0];

        let nx = vy * uz - vz * uy;
        let ny = vz * ux - vx * uz;
        let nz = vx * uy - vy * ux;
        if (nz > 0) {
          nx = -nx; ny = -ny; nz = -nz;
        }
        const mag = Math.sqrt(nx * nx + ny * ny + nz * nz) + 1e-12;
        Nx_l[i][k] = nx / mag;
        Ny_l[i][k] = ny / mag;
        Nz_l[i][k] = nz / mag;
      }
    }
    return { Nx_l, Ny_l, Nz_l };
  }

  function generateForebodySurfaces(cfg) {
    const L = cfg.length;
    const bHalf = cfg.span;
    const Nx = cfg.numStreamwise;
    const Ny = 2 * cfg.numHalfSpan - 1;
    const iC = cfg.numHalfSpan - 1;
    const betaRad = cfg.beta * DEG2RAD;
    const tanBeta = Math.tan(betaRad);
    const invP = 1.0 / cfg.power;

    const eta = new Float64Array(Ny);
    const xLE = new Float64Array(Ny);
    for (let k = 0; k < Ny; k++) {
      eta[k] = -1.0 + (2.0 * k) / (Ny - 1);
      xLE[k] = L * Math.pow(Math.abs(eta[k]), invP);
    }
    xLE[0] = L;
    xLE[Ny - 1] = L;
    xLE[iC] = 0.0;

    const X_low = Array.from({ length: Nx }, () => new Array(Ny).fill(0));
    const Y_low = Array.from({ length: Nx }, () => new Array(Ny).fill(0));
    const Z_low = Array.from({ length: Nx }, () => new Array(Ny).fill(0));
    const X_shk = Array.from({ length: Nx }, () => new Array(Ny).fill(0));
    const Y_shk = Array.from({ length: Nx }, () => new Array(Ny).fill(0));
    const Z_shk = Array.from({ length: Nx }, () => new Array(Ny).fill(0));
    const Cp_low = Array.from({ length: Nx }, () => new Array(Ny).fill(0));

    let flowSolution = null;
    let planformClampCount = 0;

    if (cfg.method === "conical" || cfg.method === "osculating_cone") {
      const tm = solveTaylorMaccoll(cfg.mach, cfg.beta, cfg.gamma);
      flowSolution = tm;
      const tanDelta = Math.tan(tm.shock.deltaRad);
      const CpWedge = tm.shock.Cp2;

      if (cfg.method === "conical") {
        // Axisymmetric cone-derived waverider (Jones et al. 1968; Rasmussen 1980)
        const minR_L = 1.28 * bHalf;
        const x0_monotone = Math.max((invP - 1.0 + 0.08) * L, 0.0);
        const x0_min = Math.max(minR_L / tanBeta - L, 0.15 * L, x0_monotone);
        const vOff = Math.min(Math.max(cfg.vertexOffsetRatio, 0.05), 0.92);
        const x0_cfg = (vOff / (1.0 - vOff)) * L;
        const x0 = Math.max(x0_cfg, x0_min);
        const zv = x0 * tanBeta;

        const xStartLocal = new Float64Array(Ny);
        const xEndLocal = new Float64Array(Ny);
        const phiLE = new Float64Array(Ny);

        for (let k = 0; k < Ny; k++) {
          const yk = bHalf * eta[k];
          xStartLocal[k] = xLE[k] + x0;
          xEndLocal[k] = L + x0;
          const R_shk = xStartLocal[k] * tanBeta;
          const zDepth = Math.sqrt(Math.max(R_shk * R_shk - yk * yk, 1e-12));
          phiLE[k] = Math.atan2(yk, zDepth);
        }
        phiLE[iC] = 0.0;

        const traced = traceStreamline(tm, xStartLocal, xEndLocal, Nx);
        for (let i = 0; i < Nx; i++) {
          for (let k = 0; k < Ny; k++) {
            const xLoc = traced.xGrid[i][k];
            const rVal = traced.rGrid[i][k];
            const rShk = xLoc * tanBeta;
            const sPhi = k === iC ? 0.0 : Math.sin(phiLE[k]);
            const cPhi = k === iC ? 1.0 : Math.cos(phiLE[k]);

            X_low[i][k] = xLoc - x0;
            Y_low[i][k] = rVal * sPhi;
            Z_low[i][k] = zv - rVal * cPhi;
            Cp_low[i][k] = traced.cpGrid[i][k];

            X_shk[i][k] = X_low[i][k];
            Y_shk[i][k] = rShk * sPhi;
            Z_shk[i][k] = zv - rShk * cPhi;
          }
          Y_low[i][iC] = 0.0;
          Y_shk[i][iC] = 0.0;
        }
      } else {
        // Sobieczky, Dougherty & Jones (1990) inverse osculating-cones method
        const gammaS = Math.max(cfg.anhedral, 0.0) * DEG2RAD;
        const yF = cfg.flatShockFraction * bHalf;
        const nExp = cfg.shockCurvaturePower;
        const wOut = Math.max(bHalf - yF, 1e-9);
        const hRise = wOut * Math.tan(gammaS);
        const zC = -L * tanBeta;

        const yS = new Float64Array(Ny);
        const zS = new Float64Array(Ny);
        const dzdy = new Float64Array(Ny);
        const d2z = new Float64Array(Ny);
        const flat = new Array(Ny).fill(false);
        const Rc = new Float64Array(Ny);
        const ncY = new Float64Array(Ny);
        const ncZ = new Float64Array(Ny);

        const s0 = 0.18;
        const denomG = Math.pow(1.0 + s0 * s0, 0.5 * nExp) - Math.pow(s0, nExp);

        for (let k = 0; k < Ny; k++) {
          yS[k] = bHalf * eta[k];
          const inFlatSeg = Math.abs(yS[k]) <= yF + 1e-12;
          const s = Math.max(Math.abs(yS[k]) - yF, 0.0) / wOut;
          let gVal = 0.0, dgds = 0.0, d2gds2 = 0.0;
          if (!inFlatSeg) {
            if (nExp >= 2.0) {
              gVal = Math.pow(s, nExp);
              dgds = nExp * Math.pow(s, nExp - 1.0);
              d2gds2 = nExp * (nExp - 1.0) * Math.pow(s, nExp - 2.0);
            } else {
              const u2 = s * s + s0 * s0;
              gVal = (Math.pow(u2, 0.5 * nExp) - Math.pow(s0, nExp)) / denomG;
              dgds = (nExp * s * Math.pow(u2, 0.5 * nExp - 1.0)) / denomG;
              d2gds2 = (nExp * Math.pow(u2, 0.5 * nExp - 2.0) * ((nExp - 1.0) * s * s + s0 * s0)) / denomG;
            }
          }
          zS[k] = zC + hRise * gVal;
          const sgn = yS[k] > 0 ? 1.0 : yS[k] < 0 ? -1.0 : 0.0;
          dzdy[k] = sgn * (hRise / wOut) * dgds;
          d2z[k] = (hRise / (wOut * wOut)) * d2gds2;

          let isFlat = inFlatSeg || d2z[k] <= 1e-12 || hRise <= 0.0;
          let rcVal = Infinity;
          if (!isFlat) {
            rcVal = Math.pow(1.0 + dzdy[k] * dzdy[k], 1.5) / d2z[k];
            if (rcVal > 1e4 * L) {
              isFlat = true;
              rcVal = Infinity;
            }
          }
          flat[k] = isFlat;
          Rc[k] = rcVal;

          const nrm = Math.sqrt(1.0 + dzdy[k] * dzdy[k]);
          ncY[k] = -dzdy[k] / nrm;
          ncZ[k] = 1.0 / nrm;
        }

        // Evolute non-folding guard: d_LE <= 0.82 * Rc and outboard-monotone x_LE
        for (let k = iC + 1; k < Ny; k++) {
          let dLE = (L - xLE[k]) * tanBeta;
          if (!flat[k] && dLE > 0.82 * Rc[k]) {
            dLE = 0.82 * Rc[k];
            xLE[k] = L - dLE / tanBeta;
            planformClampCount++;
          }
          if (xLE[k] < xLE[k - 1]) {
            xLE[k] = xLE[k - 1];
            planformClampCount++;
          }
        }
        for (let k = 0; k < iC; k++) {
          xLE[k] = xLE[Ny - 1 - k];
        }
        xLE[0] = L;
        xLE[Ny - 1] = L;
        xLE[iC] = 0.0;

        const curvedCols = [];
        for (let k = 0; k < Ny; k++) {
          if (flat[k]) {
            const zLEf = -xLE[k] * tanBeta;
            for (let i = 0; i < Nx; i++) {
              const xi = i / (Nx - 1);
              const xv = xLE[k] + xi * (L - xLE[k]);
              X_low[i][k] = xv;
              Y_low[i][k] = yS[k];
              Z_low[i][k] = zLEf - (xv - xLE[k]) * tanDelta;
              Cp_low[i][k] = CpWedge;
              X_shk[i][k] = xv;
              Y_shk[i][k] = yS[k];
              Z_shk[i][k] = -xv * tanBeta;
            }
          } else {
            curvedCols.push(k);
          }
        }

        if (curvedCols.length > 0) {
          const nCurv = curvedCols.length;
          const xApexArr = new Float64Array(nCurv);
          const xStartArr = new Float64Array(nCurv);
          const xEndArr = new Float64Array(nCurv);
          for (let m = 0; m < nCurv; m++) {
            const k = curvedCols[m];
            xApexArr[m] = L - Rc[k] / tanBeta;
            xStartArr[m] = xLE[k] - xApexArr[m];
            xEndArr[m] = L - xApexArr[m];
          }
          const traced = traceStreamline(tm, xStartArr, xEndArr, Nx);
          for (let m = 0; m < nCurv; m++) {
            const k = curvedCols[m];
            const cY = yS[k] + Rc[k] * ncY[k];
            const cZ = zS[k] + Rc[k] * ncZ[k];
            const xAp = xApexArr[m];
            for (let i = 0; i < Nx; i++) {
              const xv = traced.xGrid[i][m] + xAp;
              const rv = traced.rGrid[i][m];
              const rS = (xv - xAp) * tanBeta;
              X_low[i][k] = xv;
              Y_low[i][k] = cY - rv * ncY[k];
              Z_low[i][k] = cZ - rv * ncZ[k];
              Cp_low[i][k] = traced.cpGrid[i][m];
              X_shk[i][k] = xv;
              Y_shk[i][k] = cY - rS * ncY[k];
              Z_shk[i][k] = cZ - rS * ncZ[k];
            }
          }
        }

        for (let i = 0; i < Nx; i++) {
          Y_low[i][iC] = 0.0;
          Y_shk[i][iC] = 0.0;
          Y_low[i][0] = Y_low[0][0];
          Z_low[i][0] = Z_low[0][0];
          Y_low[i][Ny - 1] = Y_low[0][Ny - 1];
          Z_low[i][Ny - 1] = Z_low[0][Ny - 1];
        }
      }

      const { X_up, Y_up, Z_up } = buildUpperSurface(X_low, Y_low, Z_low, cfg);
      return { X_up, Y_up, Z_up, X_low, Y_low, Z_low, X_shk, Y_shk, Z_shk, Cp_low, flowSolution, planformClampCount };
    }

    if (cfg.method === "caret") {
      // Exact Nonweiler (1959, 1963) planar oblique-shock Caret waverider
      const shock = solveObliqueShock(cfg.mach, cfg.beta, cfg.gamma);
      flowSolution = shock;
      const tanDelta = Math.tan(shock.deltaRad);
      for (let k = 0; k < Ny; k++) {
        const yk = k === iC ? 0.0 : bHalf * eta[k];
        const zLEk = k === iC ? 0.0 : -xLE[k] * tanBeta;
        for (let i = 0; i < Nx; i++) {
          const xi = i / (Nx - 1);
          const xv = xLE[k] + xi * (L - xLE[k]);
          const dx = xv - xLE[k];
          X_low[i][k] = xv;
          Y_low[i][k] = yk;
          Z_low[i][k] = zLEk - dx * tanDelta;
          Cp_low[i][k] = shock.Cp2;
          X_shk[i][k] = xv;
          Y_shk[i][k] = yk;
          Z_shk[i][k] = zLEk - dx * tanBeta;
        }
      }
      const { X_up, Y_up, Z_up } = buildUpperSurface(X_low, Y_low, Z_low, cfg);
      return { X_up, Y_up, Z_up, X_low, Y_low, Z_low, X_shk, Y_shk, Z_shk, Cp_low, flowSolution, planformClampCount };
    }

    if (cfg.method === "cst") {
      // Kulfan (2008) CST boost-glide waverider body
      const shock = solveObliqueShock(cfg.mach, cfg.beta, cfg.gamma);
      flowSolution = shock;
      const tanDelta = Math.tan(shock.deltaRad);
      const anhedralRad = cfg.anhedral * DEG2RAD;

      const weights = cfg.cstUpperWeights;
      const nDeg = weights.length - 1;
      const cstLongProfile = new Float64Array(Nx);
      let maxCST = 1e-12;
      for (let i = 0; i < Nx; i++) {
        const s = i / (Nx - 1);
        const C_long = Math.pow(s, cfg.cstNoseExponent) * Math.pow(1.0 - 0.30 * s, cfg.cstAftExponent);
        let S_b = 0.0;
        for (let r = 0; r <= nDeg; r++) {
          S_b += weights[r] * binomCoeff(nDeg, r) * Math.pow(s, r) * Math.pow(1.0 - s, nDeg - r);
        }
        cstLongProfile[i] = C_long * S_b;
        if (cstLongProfile[i] > maxCST) maxCST = cstLongProfile[i];
      }

      const X_up = Array.from({ length: Nx }, () => new Array(Ny).fill(0));
      const Y_up = Array.from({ length: Nx }, () => new Array(Ny).fill(0));
      const Z_up = Array.from({ length: Nx }, () => new Array(Ny).fill(0));

      for (let k = 0; k < Ny; k++) {
        const yk = k === iC ? 0.0 : bHalf * eta[k];
        const zLEk = k === iC
          ? 0.0
          : -0.65 * xLE[k] * tanBeta - Math.abs(yk) * Math.tan(anhedralRad) * Math.pow(Math.abs(eta[k]), 0.8);
        const C_span = Math.pow(1.0 - Math.pow(Math.abs(eta[k]), cfg.canopySpanPower), 1.5);
        const chordScale = Math.pow(Math.max(1.0 - xLE[k] / L, 0.0), 0.65);
        const archBell = Math.pow(1.0 - eta[k] * eta[k], 2.0);
        const canopyAmpl = (cfg.canopy * L) / maxCST;

        for (let i = 0; i < Nx; i++) {
          const s = i / (Nx - 1);
          const xv = xLE[k] + s * (L - xLE[k]);
          const dx = xv - xLE[k];
          const dzArch = cfg.cstLowerArchRatio * L * Math.pow(s, 0.9) * archBell * chordScale;
          let zLowVal = zLEk - dx * tanDelta + dzArch;
          zLowVal = Math.min(zLowVal, zLEk - 0.25 * dx * tanDelta);
          if (i === 0) zLowVal = zLEk;

          X_low[i][k] = xv;
          Y_low[i][k] = yk;
          Z_low[i][k] = zLowVal;

          X_up[i][k] = xv;
          Y_up[i][k] = yk;
          if (cfg.upperSurfaceType === "freestream" || cfg.canopy <= 0) {
            Z_up[i][k] = zLEk;
          } else {
            Z_up[i][k] = zLEk + canopyAmpl * cstLongProfile[i] * C_span * chordScale;
          }

          X_shk[i][k] = xv;
          Y_shk[i][k] = yk;
          Z_shk[i][k] = zLEk - dx * tanBeta;
        }
        Z_up[0][k] = Z_low[0][k];
      }
      for (let i = 0; i < Nx; i++) {
        X_up[i][0] = X_low[i][0]; Y_up[i][0] = Y_low[i][0]; Z_up[i][0] = Z_low[i][0];
        X_up[i][Ny - 1] = X_low[i][Ny - 1]; Y_up[i][Ny - 1] = Y_low[i][Ny - 1]; Z_up[i][Ny - 1] = Z_low[i][Ny - 1];
      }

      // Windward Modified Newtonian + Leeward Prandtl-Meyer Cp_low
      const g = cfg.gamma, M = cfg.mach;
      const CpMax = (2.0 / (g * M * M)) * (
        Math.pow(((g + 1.0) * (g + 1.0) * M * M) / (4.0 * g * M * M - 2.0 * (g - 1.0)), g / (g - 1.0)) *
        ((1.0 - g + 2.0 * g * M * M) / (g + 1.0)) - 1.0
      );
      const { Nx_l } = computeStructuredLowerNormals(X_low, Y_low, Z_low);
      const pmExp = (2.0 * g) / (g - 1.0);
      for (let i = 0; i < Nx; i++) {
        for (let k = 0; k < Ny; k++) {
          const sinImp = Math.max(-Nx_l[i][k], 0.0);
          if (sinImp > 0.0) {
            Cp_low[i][k] = CpMax * sinImp * sinImp;
          } else {
            const sinLee = Math.max(Nx_l[i][k], 0.0);
            Cp_low[i][k] = (2.0 / (g * M * M)) * (Math.pow(Math.max(1.0 - 0.5 * (g - 1.0) * M * sinLee, 0.0), pmExp) - 1.0);
          }
        }
      }

      return { X_up, Y_up, Z_up, X_low, Y_low, Z_low, X_shk, Y_shk, Z_shk, Cp_low, flowSolution, planformClampCount };
    }

    if (cfg.method === "rv") {
      return generateReentryVehicleSurfaces(cfg);
    }

    throw new Error(`Unsupported waverider method "${cfg.method}".`);
  }

  function solveInverseTaylorMaccoll(M, targetConeDeg, gamma = 1.4) {
    const muDeg = Math.asin(1.0 / M) * RAD2DEG;
    let bLo = Math.max(muDeg + 0.12, targetConeDeg + 0.05);
    let bHi = Math.min(64.0, Math.max(targetConeDeg * 2.4 + muDeg + 4.0, bLo + 8.0));
    if (bLo >= bHi) return { ok: false, betaDeg: NaN, tmSol: null };
    let sLo, sHi;
    try {
      sLo = solveTaylorMaccoll(M, bLo, gamma);
      sHi = solveTaylorMaccoll(M, bHi, gamma);
    } catch (_) {
      return { ok: false, betaDeg: NaN, tmSol: null };
    }
    let fLo = sLo.thetaConeDeg - targetConeDeg;
    let fHi = sHi.thetaConeDeg - targetConeDeg;
    if (fLo > 0 || fHi < 0) return { ok: false, betaDeg: NaN, tmSol: null };

    for (let it = 1; it <= 18; it++) {
      let bMid;
      if (it % 3 === 0) {
        bMid = 0.5 * (bLo + bHi);
      } else {
        const w = Math.max(0.15, Math.min(0.85, -fLo / (fHi - fLo)));
        bMid = bLo + w * (bHi - bLo);
      }
      const sMid = solveTaylorMaccoll(M, bMid, gamma);
      const fMid = sMid.thetaConeDeg - targetConeDeg;
      if (Math.abs(fMid) < 1e-6 || bHi - bLo < 1e-6) {
        return { ok: true, betaDeg: bMid, tmSol: sMid };
      }
      if (fMid < 0) {
        bLo = bMid; fLo = fMid;
      } else {
        bHi = bMid; fHi = fMid;
      }
    }
    const betaDeg = 0.5 * (bLo + bHi);
    return { ok: true, betaDeg, tmSol: solveTaylorMaccoll(M, betaDeg, gamma) };
  }

  function interp1Linear(xArr, yArr, xQuery) {
    const n = xArr.length;
    if (n === 0) return 0.0;
    if (n === 1 || xQuery <= xArr[0]) return yArr[0];
    if (xQuery >= xArr[n - 1]) return yArr[n - 1];
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (xArr[mid] <= xQuery) lo = mid;
      else hi = mid;
    }
    const dx = xArr[hi] - xArr[lo];
    if (Math.abs(dx) < 1e-15) return yArr[lo];
    return yArr[lo] + ((xQuery - xArr[lo]) / dx) * (yArr[hi] - yArr[lo]);
  }

  function generateReentryVehicleSurfaces(cfg) {
    const L = cfg.length;
    const M = cfg.mach;
    const g = cfg.gamma;
    const bodyType = cfg.rvBodyType;
    let noseType = cfg.rvNoseType;
    const aftType = cfg.rvAftType;

    const nStages = bodyType === "mono_conic" ? 1 : (bodyType === "bi_conic" ? 2 : 3);
    const coneAnglesDeg = cfg.rvConeAngles.slice(0, nStages);
    const rawWeights = cfg.rvStageFractions.slice(0, nStages).map((w) => Math.max(w, 0.0));
    const sumW = rawWeights.reduce((acc, v) => acc + v, 0.0);
    const stageWeights = sumW > 0 ? rawWeights.map((w) => w / sumW) : Array(nStages).fill(1.0 / nStages);

    const th1Deg = Math.max(1.0, Math.min(52.0, coneAnglesDeg[0]));
    const th1Rad = th1Deg * DEG2RAD;
    const tanStages = coneAnglesDeg.map((a) => Math.tan(a * DEG2RAD));
    let meanTanPostNose = 0.0;
    for (let s = 0; s < nStages; s++) meanTanPostNose += stageWeights[s] * tanStages[s];
    meanTanPostNose = Math.max(meanTanPostNose, Math.tan(2.0 * DEG2RAD));

    let etaNose = Math.max(0.0, Math.min(0.80, cfg.rvBluntness));
    if (noseType === "conical" || etaNose < 1e-4) {
      noseType = "conical";
      etaNose = 0.0;
    }

    let RN = 0.0;
    let rhoOgive = 0.0;
    let th0OgiveRad = th1Rad;
    let x_t = 0.0;
    let r_t = 0.0;

    if (noseType === "spherical") {
      const denom = 1.0 - etaNose * (Math.cos(th1Rad) - (1.0 - Math.sin(th1Rad)) * meanTanPostNose);
      const R_base_nom = (L * meanTanPostNose) / Math.max(denom, 0.15);
      const RN_max = (0.50 * L) / Math.max(1.0 - Math.sin(th1Rad), 0.05);
      RN = Math.min(etaNose * R_base_nom, RN_max);
      x_t = RN * (1.0 - Math.sin(th1Rad));
      r_t = RN * Math.cos(th1Rad);
    } else if (noseType === "ogive") {
      const lamNose = Math.min(0.56, Math.max(0.14, 0.14 + 0.50 * etaNose));
      x_t = lamNose * L;
      const dThetaDeg = 4.0 + 18.0 * etaNose;
      const thSecDeg = Math.min(th1Deg + dThetaDeg, 62.0);
      r_t = x_t * Math.tan(thSecDeg * DEG2RAD);
      const denomOg = 2.0 * (r_t * Math.cos(th1Rad) - x_t * Math.sin(th1Rad));
      rhoOgive = (x_t * x_t + r_t * r_t) / Math.max(denomOg, 1e-9);
      const sinTh0 = Math.min(Math.sin(th1Rad) + x_t / rhoOgive, 0.995);
      th0OgiveRad = Math.asin(sinTh0);
    }

    const NxTarget = Math.max(cfg.numStreamwise, 24);
    const xVec = [];
    const rVec = [];
    const stageIDPerRow = [];

    if (noseType === "conical") {
      xVec.push(0.0);
      rVec.push(0.0);
      stageIDPerRow.push(0);
    } else {
      const nNosePts = Math.max(10, Math.round(0.32 * NxTarget));
      if (noseType === "spherical") {
        const alphaMax = 0.5 * Math.PI - th1Rad;
        for (let i = 0; i < nNosePts; i++) {
          const al = (alphaMax * i) / (nNosePts - 1);
          xVec.push(RN * (1.0 - Math.cos(al)));
          rVec.push(RN * Math.sin(al));
          stageIDPerRow.push(0);
        }
      } else {
        for (let i = 0; i < nNosePts; i++) {
          const psi = th0OgiveRad + ((th1Rad - th0OgiveRad) * i) / (nNosePts - 1);
          xVec.push(rhoOgive * (Math.sin(th0OgiveRad) - Math.sin(psi)));
          rVec.push(rhoOgive * (Math.cos(psi) - Math.cos(th0OgiveRad)));
          stageIDPerRow.push(0);
        }
      }
      xVec[0] = 0.0; rVec[0] = 0.0;
      xVec[xVec.length - 1] = x_t; rVec[rVec.length - 1] = r_t;
    }

    const L_rem = Math.max(L - x_t, 1e-4);
    const nRemPts = Math.max(NxTarget - xVec.length + 1, 3 * nStages + 1);
    const stageBreaksX = new Array(nStages + 1).fill(0.0);
    const stageBreaksR = new Array(nStages + 1).fill(0.0);
    stageBreaksX[0] = x_t;
    stageBreaksR[0] = r_t;

    const rMinFloor = Math.max(0.025 * L, 0.20 * Math.max(r_t, 0.02 * L));
    for (let s = 0; s < nStages; s++) {
      const Ls = stageWeights[s] * L_rem;
      const xStart = stageBreaksX[s];
      const xEnd = s === nStages - 1 ? L : xStart + Ls;
      const rStart = stageBreaksR[s];
      const rEnd = Math.max(rStart + (xEnd - xStart) * tanStages[s], rMinFloor);
      const nSeg = Math.max(3, Math.round(stageWeights[s] * (nRemPts - 1)));

      for (let j = 1; j <= nSeg; j++) {
        const tau = j / nSeg;
        xVec.push(xStart + tau * (xEnd - xStart));
        rVec.push(rStart + tau * (rEnd - rStart));
        stageIDPerRow.push(s + 1);
      }
      stageBreaksX[s + 1] = xEnd;
      stageBreaksR[s + 1] = rEnd;
    }

    const Nx = xVec.length;
    const Ny = 2 * cfg.numHalfSpan - 1;
    const R_base = rVec[Nx - 1];

    // Safe concave aft cavity depth clamping (positive wall thickness from center to rim)
    let dAft_m = 0.0;
    const aftPower = cfg.rvAftPower;
    if (aftType === "concave" && cfg.rvAftDepth > 0) {
      let dSafe = Math.min(cfg.rvAftDepth * R_base, 0.55 * L);
      for (let iter = 0; iter < 14; iter++) {
        let ok = L - dSafe > 0.04 * L;
        for (let m = 0; m < 40 && ok; m++) {
          const u = (0.98 * m) / 39.0;
          const xC = L - dSafe * (1.0 - Math.pow(u, aftPower));
          const rSkin = xC >= xVec[0] ? interp1Linear(xVec, rVec, xC) : 0.0;
          const wallMargin = 0.04 * R_base * (1.0 - u);
          if (u * R_base > rSkin - wallMargin + 1e-12) {
            ok = false;
          }
        }
        if (ok) break;
        dSafe *= 0.82;
      }
      dAft_m = Math.max(dSafe, 0.0);
    }

    // Exact inverse Taylor-Maccoll on sharp conical noses or bow-shock envelope
    const CpMax = (2.0 / (g * M * M)) * (
      Math.pow(((g + 1.0) * (g + 1.0) * M * M) / (4.0 * g * M * M - 2.0 * (g - 1.0)), g / (g - 1.0)) *
      ((1.0 - g + 2.0 * g * M * M) / (g + 1.0)) - 1.0
    );
    const muRad = Math.asin(1.0 / M);
    let isSharpMonoTM = false;
    let flowSolution = { M1: M, gamma: g, CpMax };
    let betaAsymRad = Math.max(muRad + 2.0 * DEG2RAD, Math.min(1.18 * th1Rad + 0.65 * muRad, 64.0 * DEG2RAD));

    if (noseType === "conical") {
      const invTM = solveInverseTaylorMaccoll(M, th1Deg, g);
      if (invTM.ok && invTM.tmSol) {
        betaAsymRad = invTM.betaDeg * DEG2RAD;
        flowSolution = invTM.tmSol;
        isSharpMonoTM = nStages === 1;
      }
    }

    const sinG = new Float64Array(Ny);
    const cosG = new Float64Array(Ny);
    for (let k = 0; k < Ny; k++) {
      const eta = -1.0 + (2.0 * k) / (Ny - 1);
      const gam = 0.5 * Math.PI * eta;
      sinG[k] = Math.sin(gam);
      cosG[k] = Math.cos(gam);
    }
    sinG[cfg.numHalfSpan - 1] = 0.0;
    cosG[0] = 0.0;
    cosG[Ny - 1] = 0.0;

    const standoff = RN > 0 ? 0.143 * RN * Math.exp(3.24 / (M * M)) : 0.0;
    const rShkVec = new Array(Nx);

    const X_up = Array.from({ length: Nx }, () => new Array(Ny));
    const Y_up = Array.from({ length: Nx }, () => new Array(Ny));
    const Z_up = Array.from({ length: Nx }, () => new Array(Ny));
    const X_low = Array.from({ length: Nx }, () => new Array(Ny));
    const Y_low = Array.from({ length: Nx }, () => new Array(Ny));
    const Z_low = Array.from({ length: Nx }, () => new Array(Ny));
    const X_shk = Array.from({ length: Nx }, () => new Array(Ny));
    const Y_shk = Array.from({ length: Nx }, () => new Array(Ny));
    const Z_shk = Array.from({ length: Nx }, () => new Array(Ny));
    const Cp_low = Array.from({ length: Nx }, () => new Array(Ny));
    const Cp_up = Array.from({ length: Nx }, () => new Array(Ny));

    const pmExp = (2.0 * g) / (g - 1.0);
    for (let i = 0; i < Nx; i++) {
      const xi = xVec[i];
      const ri = rVec[i];
      const rShk = Math.max((xi + standoff) * Math.tan(betaAsymRad), 1.06 * ri + (xi > 0 ? 0.015 * R_base : 0.0));
      rShkVec[i] = rShk;
      const xShk = xi - standoff * Math.exp((-2.5 * xi) / Math.max(L, 1e-6));

      let cpRow = 0.0;
      if (isSharpMonoTM) {
        cpRow = flowSolution.CpCone;
      } else {
        const i0 = Math.max(0, i - 1);
        const i1 = Math.min(Nx - 1, i + 1);
        const drdx = (rVec[i1] - rVec[i0]) / Math.max(xVec[i1] - xVec[i0], 1e-9);
        const sinLoc = drdx / Math.sqrt(1.0 + drdx * drdx);
        if (noseType === "conical" && stageIDPerRow[i] <= 1 && Number.isFinite(flowSolution.CpCone)) {
          cpRow = flowSolution.CpCone;
        } else if (sinLoc >= 0) {
          cpRow = CpMax * sinLoc * sinLoc;
        } else {
          const sinLee = -sinLoc;
          cpRow = (2.0 / (g * M * M)) * (Math.pow(Math.max(1.0 - 0.5 * (g - 1.0) * M * sinLee, 0.0), pmExp) - 1.0);
        }
      }

      for (let k = 0; k < Ny; k++) {
        const yv = ri * sinG[k];
        const zv = ri * cosG[k];
        X_up[i][k] = xi; Y_up[i][k] = yv; Z_up[i][k] = zv;
        X_low[i][k] = xi; Y_low[i][k] = yv; Z_low[i][k] = k === 0 || k === Ny - 1 ? 0.0 : -zv;
        X_shk[i][k] = xShk; Y_shk[i][k] = rShk * sinG[k]; Z_shk[i][k] = -rShk * cosG[k];
        Cp_low[i][k] = cpRow;
        Cp_up[i][k] = cpRow;
      }
    }

    return {
      X_up, Y_up, Z_up, X_low, Y_low, Z_low, X_shk, Y_shk, Z_shk,
      Cp_low, Cp_up,
      flowSolution,
      planformClampCount: 0,
      isAxisymmetric: true,
      useExactDesignCp: isSharpMonoTM,
      rvMeta: {
        bodyType,
        noseType,
        aftType,
        aftCavityDepth_m: dAft_m,
        aftCavityPower: aftPower,
        rNose_m: RN,
        rBase_m: R_base,
        xTangency_m: x_t,
        rTangency_m: r_t,
        xProfile: xVec,
        rProfile: rVec,
        rShockProfile: rShkVec,
        stageIDPerRow,
        stageBreaksX,
        stageBreaksR,
        coneAnglesDeg
      }
    };
  }

  // --------------------------------------------------------------------------
  // 4. Conformal Booster Frustum, Mirror-Symmetric Stitcher & 2-Manifold Audit
  // --------------------------------------------------------------------------

  function buildWaveriderModel(rawCfg = {}) {
    const cfg = normalizeConfig(rawCfg);
    validateConfig(cfg);

    const surf = generateForebodySurfaces(cfg);
    const { X_up, Y_up, Z_up, X_low, Y_low, Z_low, X_shk, Y_shk, Z_shk, Cp_low } = surf;
    const NxWave = X_low.length;
    const Ny = X_low[0].length;
    const L = cfg.length;

    let NxFrustEnd = NxWave;
    let NxTotal = NxWave;

    const yTE_up = Y_up[NxWave - 1].slice();
    const zTE_up = Z_up[NxWave - 1].slice();
    const yTE_low = Y_low[NxWave - 1].slice();
    const zTE_low = Z_low[NxWave - 1].slice();

    // Area-weighted vertical centroid of waverider base at X = L
    let zCenterBooster = -0.35;
    if (Number.isFinite(cfg.frustumZOffset)) {
      zCenterBooster = cfg.frustumZOffset;
    } else {
      let baseAreaTrap = 0.0;
      let zMomTrap = 0.0;
      for (let k = 0; k < Ny - 1; k++) {
        const dy = yTE_up[k + 1] - yTE_up[k];
        const dz0 = Math.max(zTE_up[k] - zTE_low[k], 0.0);
        const dz1 = Math.max(zTE_up[k + 1] - zTE_low[k + 1], 0.0);
        const m0 = 0.5 * (zTE_up[k] * zTE_up[k] - zTE_low[k] * zTE_low[k]);
        const m1 = 0.5 * (zTE_up[k + 1] * zTE_up[k + 1] - zTE_low[k + 1] * zTE_low[k + 1]);
        baseAreaTrap += 0.5 * dy * (dz0 + dz1);
        zMomTrap += 0.5 * dy * (m0 + m1);
      }
      if (baseAreaTrap > 1e-9) {
        zCenterBooster = zMomTrap / baseAreaTrap;
      } else {
        zCenterBooster = 0.5 * (zTE_up[Math.floor(Ny / 2)] + zTE_low[Math.floor(Ny / 2)]);
      }
    }

    const eta = new Float64Array(Ny);
    for (let k = 0; k < Ny; k++) eta[k] = -1.0 + (2.0 * k) / (Ny - 1);

    if (cfg.frustum && cfg.flength > 0) {
      const Nf = cfg.numFrustumStations;
      const Rb = cfg.bradius;
      const cpLast = Cp_low[NxWave - 1].slice();

      for (let m = 1; m < Nf; m++) {
        const tau = m / (Nf - 1);
        let sigma = tau;
        if (cfg.blend === "cubic") {
          sigma = 3.0 * tau * tau - 2.0 * tau * tau * tau;
        } else if (cfg.blend === "quintic") {
          sigma = 10.0 * Math.pow(tau, 3) - 15.0 * Math.pow(tau, 4) + 6.0 * Math.pow(tau, 5);
        } else if (cfg.blend === "conic") {
          const pE = 1.45;
          sigma = Math.pow(tau, pE) / (Math.pow(tau, pE) + Math.pow(1.0 - tau, pE));
        }
        const xVal = L + tau * cfg.flength;
        const rowX_u = new Array(Ny), rowY_u = new Array(Ny), rowZ_u = new Array(Ny);
        const rowX_l = new Array(Ny), rowY_l = new Array(Ny), rowZ_l = new Array(Ny);
        const rowCp = new Array(Ny);

        for (let k = 0; k < Ny; k++) {
          const thUp = 0.5 * Math.PI * (1.0 - eta[k]);
          const thLow = -0.5 * Math.PI * (1.0 - eta[k]);
          const yRingU = Rb * Math.cos(thUp);
          const zRingU = zCenterBooster + Rb * Math.sin(thUp);
          const yRingL = Rb * Math.cos(thLow);
          const zRingL = zCenterBooster + Rb * Math.sin(thLow);

          rowX_u[k] = xVal;
          rowY_u[k] = (1.0 - sigma) * yTE_up[k] + sigma * yRingU;
          rowZ_u[k] = (1.0 - sigma) * zTE_up[k] + sigma * zRingU;

          rowX_l[k] = xVal;
          rowY_l[k] = (1.0 - sigma) * yTE_low[k] + sigma * yRingL;
          rowZ_l[k] = (1.0 - sigma) * zTE_low[k] + sigma * zRingL;
          rowCp[k] = (1.0 - 0.65 * sigma) * cpLast[k];
        }
        rowY_l[0] = rowY_u[0]; rowZ_l[0] = rowZ_u[0];
        rowY_l[Ny - 1] = rowY_u[Ny - 1]; rowZ_l[Ny - 1] = rowZ_u[Ny - 1];

        X_up.push(rowX_u); Y_up.push(rowY_u); Z_up.push(rowZ_u);
        X_low.push(rowX_l); Y_low.push(rowY_l); Z_low.push(rowZ_l);
        Cp_low.push(rowCp);
      }
      NxFrustEnd = X_low.length;
      NxTotal = NxFrustEnd;

      if (cfg.rocket && cfg.rocketLength > 0) {
        const Nr = cfg.numRocketStations;
        const Lr = cfg.rocketLength;
        const tanCone = Math.tan(cfg.rcone * DEG2RAD);
        const sinCone = Math.max(Math.sin(cfg.rcone * DEG2RAD), 0.0);
        const cpCone = 2.0 * sinCone * sinCone;

        for (let m = 1; m < Nr; m++) {
          const tauR = m / (Nr - 1);
          const dxR = tauR * Lr;
          const Rloc = Math.max(Rb + dxR * tanCone, 0.05 * Rb);
          const xVal = L + cfg.flength + dxR;

          const rowX_u = new Array(Ny), rowY_u = new Array(Ny), rowZ_u = new Array(Ny);
          const rowX_l = new Array(Ny), rowY_l = new Array(Ny), rowZ_l = new Array(Ny);
          const rowCp = new Array(Ny);

          for (let k = 0; k < Ny; k++) {
            const thUp = 0.5 * Math.PI * (1.0 - eta[k]);
            const thLow = -0.5 * Math.PI * (1.0 - eta[k]);
            rowX_u[k] = xVal;
            rowY_u[k] = Rloc * Math.cos(thUp);
            rowZ_u[k] = zCenterBooster + Rloc * Math.sin(thUp);
            rowX_l[k] = xVal;
            rowY_l[k] = Rloc * Math.cos(thLow);
            rowZ_l[k] = zCenterBooster + Rloc * Math.sin(thLow);
            rowCp[k] = cpCone;
          }
          rowY_l[0] = rowY_u[0]; rowZ_l[0] = rowZ_u[0];
          rowY_l[Ny - 1] = rowY_u[Ny - 1]; rowZ_l[Ny - 1] = rowZ_u[Ny - 1];

          X_up.push(rowX_u); Y_up.push(rowY_u); Z_up.push(rowZ_u);
          X_low.push(rowX_l); Y_low.push(rowY_l); Z_low.push(rowZ_l);
          Cp_low.push(rowCp);
        }
        NxTotal = X_low.length;
      }
    }

    // Stitch 2-Manifold Watertight Triangulation (mirror-symmetric diagonals)
    const vertices = [];
    function addVertex(x, y, z) {
      vertices.push([x, y, z]);
      return vertices.length - 1;
    }

    const vid_up = Array.from({ length: NxTotal }, () => new Int32Array(Ny));
    const vid_low = Array.from({ length: NxTotal }, () => new Int32Array(Ny));
    const faces = [];
    const facePatch = [];    // 1=Upper, 2=Lower, 3=Base, 4=Frustum, 5=Rocket
    const faceCpDesign = []; // Exact nodal average on calibrated patches, NaN elsewhere
    const halfCols = Math.floor((Ny - 1) / 2);

    function triangulateSubGrid(grid, rStart, rEnd, reverse, patchCode, cpGrid = null) {
      for (let r = rStart; r < rEnd - 1; r++) {
        for (let c = 0; c < Ny - 1; c++) {
          const A = grid[r][c];
          const B = grid[r + 1][c];
          const C = grid[r + 1][c + 1];
          const D = grid[r][c + 1];

          if (c < halfCols) {
            // Port half: diagonal A-C
            if (A !== B && B !== C && C !== A) {
              faces.push(reverse ? [A, C, B] : [A, B, C]);
              facePatch.push(patchCode);
              faceCpDesign.push(cpGrid ? (cpGrid[r][c] + cpGrid[r + 1][c] + cpGrid[r + 1][c + 1]) / 3.0 : NaN);
            }
            if (A !== C && C !== D && D !== A) {
              faces.push(reverse ? [A, D, C] : [A, C, D]);
              facePatch.push(patchCode);
              faceCpDesign.push(cpGrid ? (cpGrid[r][c] + cpGrid[r + 1][c + 1] + cpGrid[r][c + 1]) / 3.0 : NaN);
            }
          } else {
            // Starboard half: mirror diagonal B-D
            if (A !== B && B !== D && D !== A) {
              faces.push(reverse ? [A, D, B] : [A, B, D]);
              facePatch.push(patchCode);
              faceCpDesign.push(cpGrid ? (cpGrid[r][c] + cpGrid[r + 1][c] + cpGrid[r][c + 1]) / 3.0 : NaN);
            }
            if (B !== C && C !== D && D !== B) {
              faces.push(reverse ? [B, D, C] : [B, C, D]);
              facePatch.push(patchCode);
              faceCpDesign.push(cpGrid ? (cpGrid[r + 1][c] + cpGrid[r + 1][c + 1] + cpGrid[r][c + 1]) / 3.0 : NaN);
            }
          }
        }
      }
    }

    let aftAperture = null;
    if (surf.isAxisymmetric) {
      const hasAftExt = cfg.frustum && cfg.flength > 0;
      const dAft = !hasAftExt && surf.rvMeta ? surf.rvMeta.aftCavityDepth_m : 0.0;
      const pAft = surf.rvMeta ? surf.rvMeta.aftCavityPower : 2.0;
      const Nb = dAft > 0 ? Math.max(cfg.numBaseLayers, 24) : Math.max(cfg.numBaseLayers, 3);
      const pGrid = dAft > 0 ? Math.max(1.0, Math.min(pAft, 2.5)) : 1.0;
      const midCol = Math.floor(Ny / 2);

      // 1. Single shared nose apex pole at row 0
      const vApex = addVertex(X_up[0][midCol], 0.0, 0.0);
      for (let k = 0; k < Ny; k++) {
        vid_up[0][k] = vApex;
        vid_low[0][k] = vApex;
      }

      // 2. Shared port/starboard lateral seam along rows 1..NxTotal-1
      for (let i = 1; i < NxTotal; i++) {
        const vL = addVertex(X_up[i][0], Y_up[i][0], Z_up[i][0]);
        const vR = addVertex(X_up[i][Ny - 1], Y_up[i][Ny - 1], Z_up[i][Ny - 1]);
        vid_up[i][0] = vL; vid_low[i][0] = vL;
        vid_up[i][Ny - 1] = vR; vid_low[i][Ny - 1] = vR;
      }

      // 3. Interior circumferential stations along skin rows 1..NxTotal-1
      for (let k = 1; k < Ny - 1; k++) {
        for (let i = 1; i < NxTotal; i++) {
          vid_up[i][k] = addVertex(X_up[i][k], Y_up[i][k], Z_up[i][k]);
        }
        for (let i = 1; i < NxTotal; i++) {
          vid_low[i][k] = addVertex(X_low[i][k], Y_low[i][k], Z_low[i][k]);
        }
      }

      // 4. Polar aft base cap (flat bulkhead or recessed concave dome)
      const vid_base_up = Array.from({ length: Nb }, () => new Int32Array(Ny));
      const vid_base_low = Array.from({ length: Nb }, () => new Int32Array(Ny));
      const xRim = X_up[NxTotal - 1][midCol];
      const zCenterAft = 0.5 * (Z_up[NxTotal - 1][midCol] + Z_low[NxTotal - 1][midCol]);
      if (dAft > 1e-6 && surf.rvMeta) {
        aftAperture = {
          active: true,
          xRim,
          zCenter: zCenterAft,
          rBaseSq: surf.rvMeta.rBase_m * surf.rvMeta.rBase_m
        };
      }

      for (let k = 0; k < Ny; k++) {
        vid_base_up[0][k] = vid_up[NxTotal - 1][k];
        vid_base_low[0][k] = vid_low[NxTotal - 1][k];
      }
      const xBaseCenter = dAft > 0 ? xRim - dAft : xRim;
      const vBaseCenter = addVertex(xBaseCenter, 0.0, zCenterAft);
      for (let k = 0; k < Ny; k++) {
        vid_base_up[Nb - 1][k] = vBaseCenter;
        vid_base_low[Nb - 1][k] = vBaseCenter;
      }

      for (let m = 1; m < Nb - 1; m++) {
        const u = Math.pow(1.0 - m / (Nb - 1), 1.0 / pGrid);
        const xRing = dAft > 0 ? xRim - dAft * (1.0 - Math.pow(u, pAft)) : xRim;
        const vBL = addVertex(xRing, u * Y_up[NxTotal - 1][0], zCenterAft + u * (Z_up[NxTotal - 1][0] - zCenterAft));
        const vBR = addVertex(xRing, u * Y_up[NxTotal - 1][Ny - 1], zCenterAft + u * (Z_up[NxTotal - 1][Ny - 1] - zCenterAft));
        vid_base_up[m][0] = vBL; vid_base_low[m][0] = vBL;
        vid_base_up[m][Ny - 1] = vBR; vid_base_low[m][Ny - 1] = vBR;
        for (let k = 1; k < Ny - 1; k++) {
          vid_base_up[m][k] = addVertex(
            xRing,
            u * Y_up[NxTotal - 1][k],
            zCenterAft + u * (Z_up[NxTotal - 1][k] - zCenterAft)
          );
          vid_base_low[m][k] = addVertex(
            xRing,
            u * Y_low[NxTotal - 1][k],
            zCenterAft + u * (Z_low[NxTotal - 1][k] - zCenterAft)
          );
        }
      }

      const cpUpGrid = surf.useExactDesignCp ? surf.Cp_up : null;
      const cpLowGrid = surf.useExactDesignCp ? Cp_low : null;
      triangulateSubGrid(vid_up, 0, NxWave, false, 1, cpUpGrid);
      triangulateSubGrid(vid_low, 0, NxWave, true, 2, cpLowGrid);
      if (NxFrustEnd > NxWave) {
        triangulateSubGrid(vid_up, NxWave - 1, NxFrustEnd, false, 4, null);
        triangulateSubGrid(vid_low, NxWave - 1, NxFrustEnd, true, 4, null);
      }
      if (NxTotal > NxFrustEnd) {
        triangulateSubGrid(vid_up, NxFrustEnd - 1, NxTotal, false, 5, null);
        triangulateSubGrid(vid_low, NxFrustEnd - 1, NxTotal, true, 5, null);
      }
      triangulateSubGrid(vid_base_up, 0, Nb, false, 3, null);
      triangulateSubGrid(vid_base_low, 0, Nb, true, 3, null);
    } else {
      const Nb = cfg.numBaseLayers;
      const vid_base = Array.from({ length: Nb }, () => new Int32Array(Ny));

      const vLeftTip = addVertex(X_up[0][0], Y_up[0][0], Z_up[0][0]);
      const vRightTip = addVertex(X_up[0][Ny - 1], Y_up[0][Ny - 1], Z_up[0][Ny - 1]);

      for (let i = 0; i < NxWave; i++) {
        vid_up[i][0] = vLeftTip; vid_low[i][0] = vLeftTip;
        vid_up[i][Ny - 1] = vRightTip; vid_low[i][Ny - 1] = vRightTip;
      }
      for (let i = NxWave; i < NxTotal; i++) {
        const vL = addVertex(X_up[i][0], Y_up[i][0], Z_up[i][0]);
        const vR = addVertex(X_up[i][Ny - 1], Y_up[i][Ny - 1], Z_up[i][Ny - 1]);
        vid_up[i][0] = vL; vid_low[i][0] = vL;
        vid_up[i][Ny - 1] = vR; vid_low[i][Ny - 1] = vR;
      }
      for (let m = 0; m < Nb; m++) {
        vid_base[m][0] = vid_up[NxTotal - 1][0];
        vid_base[m][Ny - 1] = vid_up[NxTotal - 1][Ny - 1];
      }

      for (let k = 1; k < Ny - 1; k++) {
        const vLE = addVertex(X_up[0][k], Y_up[0][k], Z_up[0][k]);
        vid_up[0][k] = vLE; vid_low[0][k] = vLE;

        for (let i = 1; i < NxTotal; i++) {
          vid_up[i][k] = addVertex(X_up[i][k], Y_up[i][k], Z_up[i][k]);
        }
        for (let i = 1; i < NxTotal; i++) {
          vid_low[i][k] = addVertex(X_low[i][k], Y_low[i][k], Z_low[i][k]);
        }

        vid_base[0][k] = vid_up[NxTotal - 1][k];
        for (let m = 1; m < Nb - 1; m++) {
          const tau = m / (Nb - 1);
          const xb = X_up[NxTotal - 1][k] * (1.0 - tau) + X_low[NxTotal - 1][k] * tau;
          const yb = Y_up[NxTotal - 1][k] * (1.0 - tau) + Y_low[NxTotal - 1][k] * tau;
          const zb = Z_up[NxTotal - 1][k] * (1.0 - tau) + Z_low[NxTotal - 1][k] * tau;
          vid_base[m][k] = addVertex(xb, yb, zb);
        }
        vid_base[Nb - 1][k] = vid_low[NxTotal - 1][k];
      }

      triangulateSubGrid(vid_up, 0, NxWave, false, 1, null);
      triangulateSubGrid(vid_low, 0, NxWave, true, 2, Cp_low);
      if (NxFrustEnd > NxWave) {
        triangulateSubGrid(vid_up, NxWave - 1, NxFrustEnd, false, 4, null);
        triangulateSubGrid(vid_low, NxWave - 1, NxFrustEnd, true, 4, null);
      }
      if (NxTotal > NxFrustEnd) {
        triangulateSubGrid(vid_up, NxFrustEnd - 1, NxTotal, false, 5, null);
        triangulateSubGrid(vid_low, NxFrustEnd - 1, NxTotal, true, 5, null);
      }
      triangulateSubGrid(vid_base, 0, Nb, false, 3, null);
    }

    // Precompute typed-array triangle geometry + full 2-manifold topology audit
    const nF = faces.length;
    const faceNormals = new Array(nF);
    const faceCenters = new Array(nF);
    const faceAreas = new Float64Array(nF);
    const crossVecs = new Array(nF);

    // Flat typed arrays for fast vertex-phasor Gordon PO facet phase integration
    const nV = vertices.length;
    const poGeom = {
      nF,
      nV,
      facePatch: Int32Array.from(facePatch),
      aftAperture,
      vx: new Float64Array(nV), vy: new Float64Array(nV), vz: new Float64Array(nV),
      f1: new Int32Array(nF), f2: new Int32Array(nF), f3: new Int32Array(nF),
      uDotV: new Float64Array(nV), cosV: new Float64Array(nV), sinV: new Float64Array(nV),
      vStamp: new Int32Array(nV), stampCounter: 0,
      litVerts: new Int32Array(nV), litFaces: new Int32Array(nF),
      cosLit: new Float64Array(nF), ut2Lit: new Float64Array(nF),
      um1Lit: new Float64Array(nF), um2Lit: new Float64Array(nF), um3Lit: new Float64Array(nF),
      nx: new Float64Array(nF), ny: new Float64Array(nF), nz: new Float64Array(nF),
      rcx: new Float64Array(nF), rcy: new Float64Array(nF), rcz: new Float64Array(nF),
      area: new Float64Array(nF), Lmax2: new Float64Array(nF),
      m1x: new Float64Array(nF), m1y: new Float64Array(nF), m1z: new Float64Array(nF),
      m2x: new Float64Array(nF), m2y: new Float64Array(nF), m2z: new Float64Array(nF),
      m3x: new Float64Array(nF), m3y: new Float64Array(nF), m3z: new Float64Array(nF)
    };

    let xMin = Infinity, xMax = -Infinity;
    let yMin = Infinity, yMax = -Infinity;
    let zMin = Infinity, zMax = -Infinity;
    for (let v = 0; v < nV; v++) {
      const pt = vertices[v];
      poGeom.vx[v] = pt[0]; poGeom.vy[v] = pt[1]; poGeom.vz[v] = pt[2];
      if (pt[0] < xMin) xMin = pt[0]; if (pt[0] > xMax) xMax = pt[0];
      if (pt[1] < yMin) yMin = pt[1]; if (pt[1] > yMax) yMax = pt[1];
      if (pt[2] < zMin) zMin = pt[2]; if (pt[2] > zMax) zMax = pt[2];
    }
    const bboxSpan = Math.max(xMax - xMin, yMax - yMin, zMax - zMin);
    const degAreaTol = Math.max(1e-13, 1e-12 * bboxSpan * bboxSpan);

    const edgeCountMap = new Map();
    const edgeSignMap = new Map();
    let enclosedVol = 0.0;
    let momX = 0.0, momY = 0.0, momZ = 0.0;
    let planformArea = 0.0;
    let upperArea = 0.0, lowerArea = 0.0, frustumArea = 0.0, rocketArea = 0.0, baseArea = 0.0, projBaseArea = 0.0;
    let degenerateFaces = 0;

    for (let f = 0; f < nF; f++) {
      const [i1, i2, i3] = faces[f];
      const v1 = vertices[i1], v2 = vertices[i2], v3 = vertices[i3];

      const he = [[i1, i2], [i2, i3], [i3, i1]];
      for (let e = 0; e < 3; e++) {
        const u = he[e][0], v = he[e][1];
        const key = u < v ? `${u}_${v}` : `${v}_${u}`;
        const sgn = u < v ? 1 : -1;
        edgeCountMap.set(key, (edgeCountMap.get(key) || 0) + 1);
        edgeSignMap.set(key, (edgeSignMap.get(key) || 0) + sgn);
      }

      const ux = v2[0] - v1[0], uy = v2[1] - v1[1], uz = v2[2] - v1[2];
      const wx = v3[0] - v1[0], wy = v3[1] - v1[1], wz = v3[2] - v1[2];
      const cx = 0.5 * (uy * wz - uz * wy);
      const cy = 0.5 * (uz * wx - ux * wz);
      const cz = 0.5 * (ux * wy - uy * wx);
      const trueArea = Math.sqrt(cx * cx + cy * cy + cz * cz);
      if (trueArea < degAreaTol) degenerateFaces++;
      const area = trueArea + 1e-15;

      const nx = cx / area, ny = cy / area, nz = cz / area;
      const rcx = (v1[0] + v2[0] + v3[0]) / 3.0;
      const rcy = (v1[1] + v2[1] + v3[1]) / 3.0;
      const rcz = (v1[2] + v2[2] + v3[2]) / 3.0;

      faceAreas[f] = trueArea;
      faceNormals[f] = [nx, ny, nz];
      faceCenters[f] = [rcx, rcy, rcz];
      crossVecs[f] = [cx, cy, cz];

      // Populate Gordon PO precomputed vertex indices & in-plane edge normals m_e = dr_e x n
      const e1x = v2[0] - v1[0], e1y = v2[1] - v1[1], e1z = v2[2] - v1[2];
      const e2x = v3[0] - v2[0], e2y = v3[1] - v2[1], e2z = v3[2] - v2[2];
      const e3x = v1[0] - v3[0], e3y = v1[1] - v3[1], e3z = v1[2] - v3[2];
      const l1Sq = e1x * e1x + e1y * e1y + e1z * e1z;
      const l2Sq = e2x * e2x + e2y * e2y + e2z * e2z;
      const l3Sq = e3x * e3x + e3y * e3y + e3z * e3z;

      poGeom.f1[f] = i1; poGeom.f2[f] = i2; poGeom.f3[f] = i3;
      poGeom.nx[f] = nx; poGeom.ny[f] = ny; poGeom.nz[f] = nz;
      poGeom.rcx[f] = rcx; poGeom.rcy[f] = rcy; poGeom.rcz[f] = rcz;
      poGeom.area[f] = trueArea; poGeom.Lmax2[f] = Math.max(l1Sq, l2Sq, l3Sq);

      poGeom.m1x[f] = e1y * nz - e1z * ny; poGeom.m1y[f] = e1z * nx - e1x * nz; poGeom.m1z[f] = e1x * ny - e1y * nx;
      poGeom.m2x[f] = e2y * nz - e2z * ny; poGeom.m2y[f] = e2z * nx - e2x * nz; poGeom.m2z[f] = e2x * ny - e2y * nx;
      poGeom.m3x[f] = e3y * nz - e3z * ny; poGeom.m3y[f] = e3z * nx - e3x * nz; poGeom.m3z[f] = e3x * ny - e3y * nx;

      // Signed tetrahedral volume: (1/6) * v1 . (v2 x v3)
      const v2xv3_x = v2[1] * v3[2] - v2[2] * v3[1];
      const v2xv3_y = v2[2] * v3[0] - v2[0] * v3[2];
      const v2xv3_z = v2[0] * v3[1] - v2[1] * v3[0];
      const tetVol = (v1[0] * v2xv3_x + v1[1] * v2xv3_y + v1[2] * v2xv3_z) / 6.0;
      enclosedVol += tetVol;
      momX += 0.25 * (v1[0] + v2[0] + v3[0]) * tetVol;
      momY += 0.25 * (v1[1] + v2[1] + v3[1]) * tetVol;
      momZ += 0.25 * (v1[2] + v2[2] + v3[2]) * tetVol;

      const pCode = facePatch[f];
      if (pCode === 1) { upperArea += trueArea; planformArea += Math.abs(cz); }
      else if (pCode === 2) lowerArea += trueArea;
      else if (pCode === 3) { baseArea += trueArea; projBaseArea += Math.max(cx, 0.0); }
      else if (pCode === 4) frustumArea += trueArea;
      else if (pCode === 5) rocketArea += trueArea;
    }

    if (planformArea < 1e-9) {
      for (let f = 0; f < nF; f++) {
        if (facePatch[f] === 2) planformArea += Math.abs(crossVecs[f][2]);
      }
    }

    let boundaryEdges = 0;
    let nonManifoldEdges = 0;
    let conflictingHalfEdges = 0;
    for (const [key, cnt] of edgeCountMap.entries()) {
      if (cnt === 1) boundaryEdges++;
      if (cnt !== 2) nonManifoldEdges++;
      if (edgeSignMap.get(key) !== 0) conflictingHalfEdges++;
    }

    // Structured-surface fold check on forebody (rows 0 .. NxWave-1)
    const foldTol = 1e-9 * Math.max(bboxSpan, 1.0);
    let hasSurfaceFold = false;
    for (let i = 0; i < NxWave && !hasSurfaceFold; i++) {
      for (let k = 0; k < Ny - 1; k++) {
        if (Y_low[i][k + 1] - Y_low[i][k] < -foldTol || Y_up[i][k + 1] - Y_up[i][k] < -foldTol) {
          hasSurfaceFold = true;
          break;
        }
      }
      if (i >= 1) {
        for (let k = 1; k < Ny - 1; k++) {
          if (Z_up[i][k] - Z_low[i][k] < -foldTol) {
            hasSurfaceFold = true;
            break;
          }
        }
      }
    }

    const numEdges = edgeCountMap.size;
    const euler = vertices.length - numEdges + nF;
    const volumeCentroid = Math.abs(enclosedVol) > 1e-12
      ? [momX / enclosedVol, momY / enclosedVol, momZ / enclosedVol]
      : [0.5 * (xMin + xMax), 0.0, 0.5 * (zMin + zMax)];
    const wettedArea = upperArea + lowerArea + frustumArea + rocketArea;
    const etaV = Math.pow(Math.max(enclosedVol, 0.0), 2.0 / 3.0) / Math.max(planformArea, 1e-12);
    const isWatertight =
      boundaryEdges === 0 &&
      nonManifoldEdges === 0 &&
      conflictingHalfEdges === 0 &&
      degenerateFaces === 0 &&
      !hasSurfaceFold &&
      enclosedVol > 0;

    // Build Concave Super-Cells for GO-PO Double-Bounce (Waverider Underside or RV Aft Cavity / Flared Cones)
    const concaveData = surf.isAxisymmetric
      ? buildAxisymmetricConcaveSuperCells(surf.rvMeta, cfg.frustum && cfg.flength > 0)
      : buildConcaveSuperCells(X_low, Y_low, Z_low, NxWave, Ny);

    // Precompute Base-Rim Wedge Edges for Polarimetric Ufimtsev/Mitzner PTD Edge Diffraction
    const rimEdges = buildRimPTDEdges(vertices, faces, facePatch, faceNormals, Boolean(aftAperture && aftAperture.active));

    // Precompute Aero Cache & Design-Point (alpha = 0) Metrics
    const aeroCache = buildAeroCache({
      cfg, faces, facePatch, faceCpDesign, crossVecs, faceNormals, faceCenters,
      volumeCentroid, planformArea, wettedArea, baseArea, projBaseArea
    });
    const aero0 = evaluateAeroFromCache(aeroCache, 0.0);

    return {
      cfg,
      isAxisymmetric: Boolean(surf.isAxisymmetric),
      rvMeta: surf.rvMeta || null,
      X_up, Y_up, Z_up,
      X_low, Y_low, Z_low,
      X_shk, Y_shk, Z_shk,
      Cp_low,
      flowSolution: surf.flowSolution,
      planformClampCount: surf.planformClampCount,
      NxWave, NxFrustEnd, NxTotal, Ny, zCenterBooster,
      xMin, xMax, yMin, yMax, zMin, zMax,
      vertices, faces, facePatch, faceCpDesign,
      faceNormals, faceCenters, faceAreas, crossVecs, poGeom,
      boundaryEdges, nonManifoldEdges, conflictingHalfEdges, degenerateFaces,
      hasSurfaceFold, isWatertight,
      numEdges, euler,
      enclosedVol, volumeCentroid, planformArea, wettedArea, baseArea, projBaseArea, etaV,
      aeroCache,
      CL: aero0.CL,
      CD: aero0.CD_wave,
      CD_wave: aero0.CD_wave,
      CY: aero0.CY,
      LD: aero0.LD_inviscid,
      LD_inviscid: aero0.LD_inviscid,
      CD_friction_est: aero0.CD_friction_est,
      CD_base_est: aero0.CD_base_est,
      LD_viscous_est: aero0.LD_viscous_est,
      Cm: aero0.Cm,
      rimEdges,
      concaveData,
      concavePairs: concaveData.rayPairs,
      totalConcavePairs: concaveData.numPairs
    };
  }

  // --------------------------------------------------------------------------
  // 5. Calibrated Inviscid + Prandtl-Meyer + Eckert Reference-Temp Aerodynamics
  // --------------------------------------------------------------------------

  function buildAeroCache(m) {
    const { cfg, faces, facePatch, faceCpDesign, crossVecs, faceNormals, faceCenters, volumeCentroid, planformArea, wettedArea, baseArea, projBaseArea } = m;
    const M = cfg.mach;
    const g = cfg.gamma;
    const CpMax = (2.0 / (g * M * M)) * (
      Math.pow(((g + 1.0) * (g + 1.0) * M * M) / (4.0 * g * M * M - 2.0 * (g - 1.0)), g / (g - 1.0)) *
      ((1.0 - g + 2.0 * g * M * M) / (g + 1.0)) - 1.0
    );

    const nF = faces.length;
    const kCal = new Float64Array(nF);
    const validWindward0 = new Uint8Array(nF);
    const hasCalCp = new Uint8Array(nF);

    for (let f = 0; f < nF; f++) {
      if ((facePatch[f] === 1 || facePatch[f] === 2) && Number.isFinite(faceCpDesign[f])) {
        hasCalCp[f] = 1;
        const sinImp0 = Math.max(-faceNormals[f][0], 0.0);
        if (sinImp0 > 1e-4) {
          validWindward0[f] = 1;
          kCal[f] = faceCpDesign[f] / (CpMax * sinImp0 * sinImp0);
        } else {
          kCal[f] = 1.0;
        }
      } else {
        kCal[f] = 1.0;
      }
    }

    // Eckert (1955) reference-temperature turbulent skin-friction at h = 30 km (US Std Atm)
    const T_inf = 226.51;
    const rho_inf = 0.01841;
    const a_inf = Math.sqrt(g * 287.05 * T_inf);
    const U_inf = M * a_inf;
    const T_wall = 1000.0;
    const rRec = Math.cbrt(0.72);
    const T_aw = T_inf * (1.0 + 0.5 * rRec * (g - 1.0) * M * M);
    const T_star = T_inf + 0.5 * (T_wall - T_inf) + 0.22 * (T_aw - T_inf);
    const mu_star = (1.458e-6 * Math.pow(T_star, 1.5)) / (T_star + 110.4);
    const rho_star = rho_inf * (T_inf / T_star);
    const L_char = Math.max(0.67 * cfg.length, 0.1);
    const Re_star = Math.max((rho_star * U_inf * L_char) / mu_star, 1e4);
    const Cf_star = 0.074 / Math.pow(Re_star, 0.2);
    const Cf_inf = Cf_star * (rho_star / rho_inf);

    const Sref = Math.max(planformArea, 1e-12);
    const CD_friction_est = Cf_inf * (wettedArea / Sref);
    const Cp_base = -1.0 / (g * M * M);
    const baseEff = Number.isFinite(projBaseArea) ? projBaseArea : baseArea;
    const CD_base_est = -Cp_base * (baseEff / Sref);

    return {
      M, g, CpMax, Sref, Lref: cfg.length,
      facePatch, faceCpDesign, crossVecs, faceNormals, faceCenters, volumeCentroid,
      hasCalCp, kCal, validWindward0, CD_friction_est, CD_base_est
    };
  }

  function evaluateAeroFromCache(cache, alphaDeg) {
    const aRad = alphaDeg * DEG2RAD;
    const ux = Math.cos(aRad), uz = Math.sin(aRad);
    const lx = -Math.sin(aRad), lz = Math.cos(aRad);
    const { M, g, CpMax, Sref, Lref, facePatch, faceCpDesign, crossVecs, faceNormals, faceCenters, volumeCentroid, hasCalCp, kCal, validWindward0 } = cache;
    const pmFactor = 2.0 / (g * M * M);
    const pmExp = (2.0 * g) / (g - 1.0);
    const atDesign = Math.abs(alphaDeg) < 1e-13;

    let fx = 0.0, fy = 0.0, fz = 0.0;
    let my = 0.0;

    for (let f = 0; f < facePatch.length; f++) {
      const pCode = facePatch[f];
      if (pCode === 3) continue; // Base excluded from forebody/wetted wave pressure integration
      const n = faceNormals[f];
      const nDotU = n[0] * ux + n[2] * uz;
      const sinImp = Math.max(-nDotU, 0.0);
      const sinLee = Math.max(nDotU, 0.0);
      const cpPM = pmFactor * (Math.pow(Math.max(1.0 - 0.5 * (g - 1.0) * M * sinLee, 0.0), pmExp) - 1.0);

      let cp = 0.0;
      if (hasCalCp[f]) {
        if (atDesign) {
          cp = validWindward0[f] ? faceCpDesign[f] : (sinLee > 0 ? cpPM : faceCpDesign[f]);
        } else {
          cp = sinImp > 0 ? kCal[f] * CpMax * sinImp * sinImp : cpPM;
        }
      } else {
        cp = sinImp > 0 ? CpMax * sinImp * sinImp : cpPM;
      }

      const cv = crossVecs[f];
      const dfx = -cp * cv[0];
      const dfy = -cp * cv[1];
      const dfz = -cp * cv[2];
      fx += dfx; fy += dfy; fz += dfz;

      const rx = faceCenters[f][0] - volumeCentroid[0];
      const rz = faceCenters[f][2] - volumeCentroid[2];
      // Pitching moment (r x dF)_y = rz * dfx - rx * dfz
      my += rz * dfx - rx * dfz;
    }

    const Fx = fx / Sref, Fy = fy / Sref, Fz = fz / Sref;
    const CD_wave = Fx * ux + Fz * uz;
    const CL = Fx * lx + Fz * lz;
    const CY = Fy;
    const LD_inviscid = CL / Math.max(CD_wave, 1e-9);
    const CD_total_viscous_est = CD_wave + cache.CD_friction_est;
    const LD_viscous_est = CL / Math.max(CD_total_viscous_est, 1e-9);
    const Cm = my / (Sref * Lref);

    return {
      alphaDeg, CL, CD_wave, CY, LD_inviscid,
      CD_friction_est: cache.CD_friction_est,
      CD_base_est: cache.CD_base_est,
      CD_total_viscous_est,
      LD_viscous_est,
      Cm
    };
  }

  function evaluateAeroState(model, alphaDeg = 0.0) {
    return evaluateAeroFromCache(model.aeroCache, alphaDeg);
  }

  function computeAeroPolar(model, alphaDegVec) {
    const alphas = alphaDegVec || Array.from({ length: 41 }, (_, k) => -4.0 + 0.5 * k);
    return alphas.map((a) => evaluateAeroFromCache(model.aeroCache, a));
  }

  // --------------------------------------------------------------------------
  // 6. Exact Gordon (1975) Triangle PO & GO-PO Sutherland-Hodgman Multi-Bounce
  // --------------------------------------------------------------------------

  /**
   * Standalone exact Gordon (1975) triangle contour integral I = \int_T exp(i w.r) dA
   * over an array of triangles, used by both GO-PO overlap polygons and direct benchmarks.
   */
  function facetPhaseIntegralTriangles(v1List, v2List, v3List, nList, wx, wy, wz) {
    const nT = v1List.length;
    const outRe = new Float64Array(nT);
    const outIm = new Float64Array(nT);
    const w2 = wx * wx + wy * wy + wz * wz;

    for (let t = 0; t < nT; t++) {
      const v1 = v1List[t], v2 = v2List[t], v3 = v3List[t], n = nList[t];
      const e1x = v2[0] - v1[0], e1y = v2[1] - v1[1], e1z = v2[2] - v1[2];
      const e2x = v3[0] - v2[0], e2y = v3[1] - v2[1], e2z = v3[2] - v2[2];
      const e3x = v1[0] - v3[0], e3y = v1[1] - v3[1], e3z = v1[2] - v3[2];

      const cx = e1y * (-e3z) - e1z * (-e3y);
      const cy = e1z * (-e3x) - e1x * (-e3z);
      const cz = e1x * (-e3y) - e1y * (-e3x);
      const area = 0.5 * Math.hypot(cx, cy, cz);
      const Lmax = Math.max(Math.hypot(e1x, e1y, e1z), Math.hypot(e2x, e2y, e2z), Math.hypot(e3x, e3y, e3z));

      const nDotW = n[0] * wx + n[1] * wy + n[2] * wz;
      const wt2 = Math.max(w2 - nDotW * nDotW, 0.0);

      if (Math.sqrt(wt2) * Lmax < 1e-3) {
        const rcx = (v1[0] + v2[0] + v3[0]) / 3.0;
        const rcy = (v1[1] + v2[1] + v3[1]) / 3.0;
        const rcz = (v1[2] + v2[2] + v3[2]) / 3.0;
        const ph = rcx * wx + rcy * wy + rcz * wz;
        outRe[t] = area * Math.cos(ph);
        outIm[t] = area * Math.sin(ph);
      } else {
        const ex = [e1x, e2x, e3x], ey = [e1y, e2y, e3y], ez = [e1z, e2z, e3z];
        const mx = [0.5 * (v1[0] + v2[0]), 0.5 * (v2[0] + v3[0]), 0.5 * (v3[0] + v1[0])];
        const my = [0.5 * (v1[1] + v2[1]), 0.5 * (v2[1] + v3[1]), 0.5 * (v3[1] + v1[1])];
        const mz = [0.5 * (v1[2] + v2[2]), 0.5 * (v2[2] + v3[2]), 0.5 * (v3[2] + v1[2])];
        let sumRe = 0.0, sumIm = 0.0;
        for (let e = 0; e < 3; e++) {
          const inNx = ey[e] * n[2] - ez[e] * n[1];
          const inNy = ez[e] * n[0] - ex[e] * n[2];
          const inNz = ex[e] * n[1] - ey[e] * n[0];
          const wm = inNx * wx + inNy * wy + inNz * wz;
          const x = 0.5 * (ex[e] * wx + ey[e] * wy + ez[e] * wz);
          const wr = mx[e] * wx + my[e] * wy + mz[e] * wz;
          const sinc = Math.abs(x) > 1e-12 ? Math.sin(x) / x : 1.0;
          const coef = wm * sinc;
          sumRe += coef * Math.cos(wr);
          sumIm += coef * Math.sin(wr);
        }
        // Divide by (i * wt2): (sumRe + i sumIm) / (i wt2) = (sumIm - i sumRe) / wt2
        outRe[t] = sumIm / wt2;
        outIm[t] = -sumRe / wt2;
      }
    }
    return { re: outRe, im: outIm };
  }

  /**
   * Prepare frequency-independent lit PO faces, edge projections, and unique lit vertex projections
   * for look direction u = [ux, uy, uz]. Returns { nLit, nLV }.
   */
  function prepareLitPOGeometry(g, ux, uy, uz) {
    const nF = g.nF;
    const stamp = ++g.stampCounter;
    const { nx, ny, nz, f1, f2, f3, vx, vy, vz, vStamp, uDotV, litVerts, litFaces, cosLit, ut2Lit, m1x, m1y, m1z, m2x, m2y, m2z, m3x, m3y, m3z, um1Lit, um2Lit, um3Lit, facePatch, rcx, rcy, rcz, aftAperture } = g;
    const hasAp = Boolean(aftAperture && aftAperture.active && facePatch);
    const xRim = hasAp ? aftAperture.xRim : 0.0;
    const zApCtr = hasAp ? aftAperture.zCenter : 0.0;
    const rBaseSq = hasAp ? aftAperture.rBaseSq : 0.0;
    let nLit = 0;
    let nLV = 0;

    for (let f = 0; f < nF; f++) {
      const cTh = nx[f] * ux + ny[f] * uy + nz[f] * uz;
      if (cTh <= 0.0) continue;
      if (hasAp && facePatch[f] === 3) {
        if (ux <= 1e-7) continue;
        const dxAp = xRim - rcx[f];
        if (dxAp > 1e-7) {
          const tAp = dxAp / ux;
          const yAp = rcy[f] + tAp * uy;
          const zAp = rcz[f] - zApCtr + tAp * uz;
          if (yAp * yAp + zAp * zAp > rBaseSq) continue;
        }
      }
      litFaces[nLit] = f;
      cosLit[nLit] = cTh;
      ut2Lit[nLit] = Math.max(1.0 - cTh * cTh, 0.0);
      um1Lit[nLit] = m1x[f] * ux + m1y[f] * uy + m1z[f] * uz;
      um2Lit[nLit] = m2x[f] * ux + m2y[f] * uy + m2z[f] * uz;
      um3Lit[nLit] = m3x[f] * ux + m3y[f] * uy + m3z[f] * uz;
      nLit++;

      const i1 = f1[f], i2 = f2[f], i3 = f3[f];
      if (vStamp[i1] !== stamp) { vStamp[i1] = stamp; litVerts[nLV++] = i1; uDotV[i1] = vx[i1] * ux + vy[i1] * uy + vz[i1] * uz; }
      if (vStamp[i2] !== stamp) { vStamp[i2] = stamp; litVerts[nLV++] = i2; uDotV[i2] = vx[i2] * ux + vy[i2] * uy + vz[i2] * uz; }
      if (vStamp[i3] !== stamp) { vStamp[i3] = stamp; litVerts[nLV++] = i3; uDotV[i3] = vx[i3] * ux + vy[i3] * uy + vz[i3] * uz; }
    }
    return { nLit, nLV };
  }

  const _poSumOut = new Float64Array(2);

  /**
   * Exact Gordon (1975) triangle contour PO integration via shared vertex phasors:
   * evaluates cos(2 k0 u.v), sin(2 k0 u.v) ONCE per unique lit vertex and integrates
   * all lit triangles with zero additional transcendental calls!
   * Optionally splats single-bounce scattering centers into HRRP range bins when accRe != null.
   */
  function evaluateLitPOAtK(g, nLit, nLV, k0, accRe = null, accIm = null, ux = 0, xRef = 0, rLim = 0, invDBin = 0, nR = 0) {
    const twoK = 2.0 * k0;
    const fourK2 = twoK * twoK;
    const pref = k0 / Math.sqrt(Math.PI);
    const { litVerts, litFaces, uDotV, cosV, sinV, f1, f2, f3, area, Lmax2, cosLit, ut2Lit, um1Lit, um2Lit, um3Lit } = g;

    for (let j = 0; j < nLV; j++) {
      const v = litVerts[j];
      const ph = twoK * uDotV[v];
      cosV[v] = Math.cos(ph);
      sinV[v] = Math.sin(ph);
    }

    let reSingle = 0.0, imSingle = 0.0;
    const splat = accRe !== null;
    const xRefProj = ux * xRef;

    for (let idx = 0; idx < nLit; idx++) {
      const f = litFaces[idx];
      const i1 = f1[f], i2 = f2[f], i3 = f3[f];
      const c1 = cosV[i1], s1 = sinV[i1];
      const c2 = cosV[i2], s2 = sinV[i2];
      const c3 = cosV[i3], s3 = sinV[i3];

      const ut2 = ut2Lit[idx];
      const wt2 = fourK2 * ut2;
      const scale = pref * cosLit[idx];
      let aRe, aIm;

      if (wt2 * Lmax2[f] < 1e-4) {
        // |w_t| L_max < 0.01 -> 3-vertex average is exact to < 1.7e-5 (0.00014 dB) with zero trig calls
        const A3 = area[f] / 3.0;
        const IRe = A3 * (c1 + c2 + c3);
        const IIm = A3 * (s1 + s2 + s3);
        aRe = -scale * IIm;
        aIm = scale * IRe;
        if (splat) {
          const uDotRc = (uDotV[i1] + uDotV[i2] + uDotV[i3]) * (1.0 / 3.0);
          const RS = -(uDotRc - xRefProj);
          const fBin = (RS + rLim) * invDBin;
          const i0 = Math.floor(fBin);
          if (i0 >= 0 && i0 < nR - 1) {
            const t = fBin - i0;
            accRe[i0] += aRe * (1.0 - t);
            accIm[i0] += aIm * (1.0 - t);
            accRe[i0 + 1] += aRe * t;
            accIm[i0 + 1] += aIm * t;
          }
        }
      } else {
        // Exact Stokes/Gordon 3-edge contour sum; symmetric midpoint expansion for |d| < 1e-3
        // guarantees bitwise interior-edge cancellation on coplanar plates in wideband HRRP!
        const d12 = twoK * (uDotV[i2] - uDotV[i1]);
        let j1Re, j1Im;
        if (d12 >= 1e-3 || d12 <= -1e-3) {
          const inv = 1.0 / d12;
          j1Re = (s2 - s1) * inv;
          j1Im = (c1 - c2) * inv;
        } else {
          const fac = 0.5 * (1.0 + d12 * d12 * (1.0 / 12.0));
          j1Re = fac * (c1 + c2);
          j1Im = fac * (s1 + s2);
        }

        const d23 = twoK * (uDotV[i3] - uDotV[i2]);
        let j2Re, j2Im;
        if (d23 >= 1e-3 || d23 <= -1e-3) {
          const inv = 1.0 / d23;
          j2Re = (s3 - s2) * inv;
          j2Im = (c2 - c3) * inv;
        } else {
          const fac = 0.5 * (1.0 + d23 * d23 * (1.0 / 12.0));
          j2Re = fac * (c2 + c3);
          j2Im = fac * (s2 + s3);
        }

        const d31 = twoK * (uDotV[i1] - uDotV[i3]);
        let j3Re, j3Im;
        if (d31 >= 1e-3 || d31 <= -1e-3) {
          const inv = 1.0 / d31;
          j3Re = (s1 - s3) * inv;
          j3Im = (c3 - c1) * inv;
        } else {
          const fac = 0.5 * (1.0 + d31 * d31 * (1.0 / 12.0));
          j3Re = fac * (c3 + c1);
          j3Im = fac * (s3 + s1);
        }

        const um1 = um1Lit[idx], um2 = um2Lit[idx], um3 = um3Lit[idx];
        const edgeScale = scale / (twoK * ut2);
        const e1Re = edgeScale * um1 * j1Re, e1Im = edgeScale * um1 * j1Im;
        const e2Re = edgeScale * um2 * j2Re, e2Im = edgeScale * um2 * j2Im;
        const e3Re = edgeScale * um3 * j3Re, e3Im = edgeScale * um3 * j3Im;
        aRe = e1Re + e2Re + e3Re;
        aIm = e1Im + e2Im + e3Im;

        if (splat) {
          const r1 = -(0.5 * (uDotV[i1] + uDotV[i2]) - xRefProj);
          const fB1 = (r1 + rLim) * invDBin;
          const k1 = Math.floor(fB1);
          if (k1 >= 0 && k1 < nR - 1) {
            const t = fB1 - k1;
            accRe[k1] += e1Re * (1.0 - t); accIm[k1] += e1Im * (1.0 - t);
            accRe[k1 + 1] += e1Re * t;     accIm[k1 + 1] += e1Im * t;
          }
          const r2 = -(0.5 * (uDotV[i2] + uDotV[i3]) - xRefProj);
          const fB2 = (r2 + rLim) * invDBin;
          const k2 = Math.floor(fB2);
          if (k2 >= 0 && k2 < nR - 1) {
            const t = fB2 - k2;
            accRe[k2] += e2Re * (1.0 - t); accIm[k2] += e2Im * (1.0 - t);
            accRe[k2 + 1] += e2Re * t;     accIm[k2 + 1] += e2Im * t;
          }
          const r3 = -(0.5 * (uDotV[i3] + uDotV[i1]) - xRefProj);
          const fB3 = (r3 + rLim) * invDBin;
          const k3 = Math.floor(fB3);
          if (k3 >= 0 && k3 < nR - 1) {
            const t = fB3 - k3;
            accRe[k3] += e3Re * (1.0 - t); accIm[k3] += e3Im * (1.0 - t);
            accRe[k3 + 1] += e3Re * t;     accIm[k3 + 1] += e3Im * t;
          }
        }
      }

      reSingle += aRe;
      imSingle += aIm;
    }

    _poSumOut[0] = reSingle;
    _poSumOut[1] = imSingle;
  }

  /**
   * Build mirror-symmetric planar CCW super-cells on the lower compression surface
   * and find all mutually visible concave pairs (P < Q).
   */
  function buildConcaveSuperCells(X_low, Y_low, Z_low, NxWave, Ny) {
    const nI = Math.min(12, NxWave);
    const nJHalf = Math.min(9, Math.floor((Ny + 1) / 2)); // Odd total nJ = 2*nJHalf - 1 for exact symmetry
    const iIdx = [];
    for (let m = 0; m < nI; m++) {
      const v = Math.round((m * (NxWave - 1)) / (nI - 1));
      if (iIdx.length === 0 || iIdx[iIdx.length - 1] !== v) iIdx.push(v);
    }
    const midCol = Math.floor(Ny / 2);
    const rightHalf = [];
    for (let m = 0; m < nJHalf; m++) {
      const v = midCol + Math.round((m * (Ny - 1 - midCol)) / (nJHalf - 1));
      if (rightHalf.length === 0 || rightHalf[rightHalf.length - 1] !== v) rightHalf.push(v);
    }
    const leftHalf = [];
    for (let m = rightHalf.length - 1; m >= 1; m--) {
      leftHalf.push(Ny - 1 - rightHalf[m]);
    }
    const jIdx = leftHalf.concat(rightHalf);

    const Nr = iIdx.length;
    const Nc = jIdx.length;
    const cellCenter = [];
    const cellNormal = [];
    const cellArea = [];
    const cellVerts = [];
    const cellRadius = [];

    for (let r = 0; r < Nr - 1; r++) {
      const r0 = iIdx[r], r1 = iIdx[r + 1];
      for (let c = 0; c < Nc - 1; c++) {
        const c0 = jIdx[c], c1 = jIdx[c + 1];
        const p1 = [X_low[r0][c0], Y_low[r0][c0], Z_low[r0][c0]];
        const p2 = [X_low[r1][c0], Y_low[r1][c0], Z_low[r1][c0]];
        const p3 = [X_low[r1][c1], Y_low[r1][c1], Z_low[r1][c1]];
        const p4 = [X_low[r0][c1], Y_low[r0][c1], Z_low[r0][c1]];

        const a1x = p3[0] - p1[0], a1y = p3[1] - p1[1], a1z = p3[2] - p1[2];
        const b1x = p2[0] - p1[0], b1y = p2[1] - p1[1], b1z = p2[2] - p1[2];
        const cr1x = a1y * b1z - a1z * b1y;
        const cr1y = a1z * b1x - a1x * b1z;
        const cr1z = a1x * b1y - a1y * b1x;

        const a2x = p4[0] - p1[0], a2y = p4[1] - p1[1], a2z = p4[2] - p1[2];
        const cr2x = a2y * a1z - a2z * a1y;
        const cr2y = a2z * a1x - a2x * a1z;
        const cr2z = a2x * a1y - a2y * a1x;

        const nvx = 0.5 * (cr1x + cr2x);
        const nvy = 0.5 * (cr1y + cr2y);
        const nvz = 0.5 * (cr1z + cr2z);
        const areaVal = Math.hypot(nvx, nvy, nvz);
        if (areaVal > 1e-6) {
          const ctr = [
            0.25 * (p1[0] + p2[0] + p3[0] + p4[0]),
            0.25 * (p1[1] + p2[1] + p3[1] + p4[1]),
            0.25 * (p1[2] + p2[2] + p3[2] + p4[2])
          ];
          const nHat = [nvx / areaVal, nvy / areaVal, nvz / areaVal];
          // Planar CCW quad about outward normal: p1 -> p4 -> p3 -> p2
          const rawQuad = [p1, p4, p3, p2];
          let maxRad = 0.0;
          const quad = rawQuad.map((q) => {
            const dotN = (q[0] - ctr[0]) * nHat[0] + (q[1] - ctr[1]) * nHat[1] + (q[2] - ctr[2]) * nHat[2];
            const qp = [q[0] - dotN * nHat[0], q[1] - dotN * nHat[1], q[2] - dotN * nHat[2]];
            const rad = Math.hypot(qp[0] - ctr[0], qp[1] - ctr[1], qp[2] - ctr[2]);
            if (rad > maxRad) maxRad = rad;
            return qp;
          });
          cellCenter.push(ctr);
          cellNormal.push(nHat);
          cellArea.push(areaVal);
          cellVerts.push(quad);
          cellRadius.push(maxRad);
        }
      }
    }

    const nCells = cellCenter.length;
    const pIdx = [];
    const qIdx = [];
    const rayPairs = [];

    for (let p = 0; p < nCells - 1; p++) {
      const rP = cellCenter[p], nP = cellNormal[p];
      for (let q = p + 1; q < nCells; q++) {
        const rQ = cellCenter[q], nQ = cellNormal[q];
        const dx = rQ[0] - rP[0], dy = rQ[1] - rP[1], dz = rQ[2] - rP[2];
        const dist = Math.hypot(dx, dy, dz);
        if (dist < 1e-4) continue;
        const dHat = [dx / dist, dy / dist, dz / dist];
        const cP = nP[0] * dHat[0] + nP[1] * dHat[1] + nP[2] * dHat[2];
        const cQ = -(nQ[0] * dHat[0] + nQ[1] * dHat[1] + nQ[2] * dHat[2]);
        if (cP > 0.015 && cQ > 0.015) {
          pIdx.push(p);
          qIdx.push(q);
          rayPairs.push({
            rP, rQ, nP, nQ, dist, dHat, cP, cQ,
            effA: Math.sqrt(cellArea[p] * cellArea[q] * cP * cQ)
          });
        }
      }
    }

    rayPairs.sort((a, b) => b.effA - a.effA);

    // Directed pairs (P -> Q) and (Q -> P) stored in flat typed arrays for zero-chasing GO-PO culling
    const nUndir = pIdx.length;
    const nDir = 2 * nUndir;
    const dirP = new Int32Array(nDir);
    const dirQ = new Int32Array(nDir);
    const pRx = new Float64Array(nDir), pRy = new Float64Array(nDir), pRz = new Float64Array(nDir);
    const pNx = new Float64Array(nDir), pNy = new Float64Array(nDir), pNz = new Float64Array(nDir);
    const qRx = new Float64Array(nDir), qRy = new Float64Array(nDir), qRz = new Float64Array(nDir);
    const qNx = new Float64Array(nDir), qNy = new Float64Array(nDir), qNz = new Float64Array(nDir);
    const pRad = new Float64Array(nDir), qRad = new Float64Array(nDir);

    for (let m = 0; m < nUndir; m++) {
      dirP[m] = pIdx[m]; dirQ[m] = qIdx[m];
      dirP[m + nUndir] = qIdx[m]; dirQ[m + nUndir] = pIdx[m];
    }
    for (let m = 0; m < nDir; m++) {
      const P = dirP[m], Q = dirQ[m];
      pRx[m] = cellCenter[P][0]; pRy[m] = cellCenter[P][1]; pRz[m] = cellCenter[P][2];
      pNx[m] = cellNormal[P][0]; pNy[m] = cellNormal[P][1]; pNz[m] = cellNormal[P][2];
      qRx[m] = cellCenter[Q][0]; qRy[m] = cellCenter[Q][1]; qRz[m] = cellCenter[Q][2];
      qNx[m] = cellNormal[Q][0]; qNy[m] = cellNormal[Q][1]; qNz[m] = cellNormal[Q][2];
      pRad[m] = cellRadius[P]; qRad[m] = cellRadius[Q];
    }

    return {
      numPairs: nUndir,
      nCells,
      cellCenter,
      cellNormal,
      cellArea,
      cellVerts,
      cellRadius,
      pIdx,
      qIdx,
      dirP, dirQ, dirGroup: null, aftAperture: null,
      pRx, pRy, pRz, pNx, pNy, pNz,
      qRx, qRy, qRz, qNx, qNy, qNz,
      pRad, qRad,
      rayPairs: rayPairs.slice(0, 60)
    };
  }

  /**
   * Build GO-PO super-cells for axisymmetric Reentry Vehicles (RVs):
   *   1) Recessed concave aft cavity (`aftCavityDepth_m > 0 && !hasAftFrustum`)
   *   2) Flared multi-conic forebody stage corners (`coneAnglesDeg[s+1] > coneAnglesDeg[s]`)
   */
  function buildAxisymmetricConcaveSuperCells(rvMeta, hasAftFrustum) {
    const cellCenter = [];
    const cellNormal = [];
    const cellArea = [];
    const cellVerts = [];
    const cellRadius = [];
    const cellGroup = []; // 1 = concave aft cavity, 2 = flared exterior skin
    const cellEr = [];

    const hasConcaveAft = Boolean(rvMeta && !hasAftFrustum && rvMeta.aftCavityDepth_m > 1e-5);
    let hasFlaredCones = false;
    if (rvMeta && Array.isArray(rvMeta.coneAnglesDeg) && rvMeta.coneAnglesDeg.length >= 2) {
      for (let s = 0; s < rvMeta.coneAnglesDeg.length - 1; s++) {
        if (rvMeta.coneAnglesDeg[s + 1] - rvMeta.coneAnglesDeg[s] > 0.25) {
          hasFlaredCones = true;
          break;
        }
      }
    }

    const nPhi = 16; // Multiple of 4 -> exact Y <-> -Y and Z <-> -Z symmetry
    const phiEdges = new Float64Array(nPhi + 1);
    for (let c = 0; c <= nPhi; c++) phiEdges[c] = -Math.PI + (2.0 * Math.PI * c) / nPhi;

    function pushQuadCell(p1, p2, p3, p4, grp, phMid) {
      const a1x = p2[0] - p1[0], a1y = p2[1] - p1[1], a1z = p2[2] - p1[2];
      const b1x = p3[0] - p1[0], b1y = p3[1] - p1[1], b1z = p3[2] - p1[2];
      const c1x = p4[0] - p1[0], c1y = p4[1] - p1[1], c1z = p4[2] - p1[2];

      const cr1x = a1y * b1z - a1z * b1y;
      const cr1y = a1z * b1x - a1x * b1z;
      const cr1z = a1x * b1y - a1y * b1x;

      const cr2x = b1y * c1z - b1z * c1y;
      const cr2y = b1z * c1x - b1x * c1z;
      const cr2z = b1x * c1y - b1y * c1x;

      const nvx = 0.5 * (cr1x + cr2x);
      const nvy = 0.5 * (cr1y + cr2y);
      const nvz = 0.5 * (cr1z + cr2z);
      const areaVal = Math.hypot(nvx, nvy, nvz);
      if (areaVal <= 1e-7) return;

      const ctr = [
        0.25 * (p1[0] + p2[0] + p3[0] + p4[0]),
        0.25 * (p1[1] + p2[1] + p3[1] + p4[1]),
        0.25 * (p1[2] + p2[2] + p3[2] + p4[2])
      ];
      const nHat = [nvx / areaVal, nvy / areaVal, nvz / areaVal];
      const rawQuad = [p1, p2, p3, p4];
      let maxRad = 0.0;
      const quad = rawQuad.map((q) => {
        const dotN = (q[0] - ctr[0]) * nHat[0] + (q[1] - ctr[1]) * nHat[1] + (q[2] - ctr[2]) * nHat[2];
        const qp = [q[0] - dotN * nHat[0], q[1] - dotN * nHat[1], q[2] - dotN * nHat[2]];
        const rad = Math.hypot(qp[0] - ctr[0], qp[1] - ctr[1], qp[2] - ctr[2]);
        if (rad > maxRad) maxRad = rad;
        return qp;
      });
      cellCenter.push(ctr);
      cellNormal.push(nHat);
      cellArea.push(areaVal);
      cellVerts.push(quad);
      cellRadius.push(maxRad);
      cellGroup.push(grp);
      cellEr.push([0.0, Math.sin(phMid), Math.cos(phMid)]);
    }

    let aftAperture = null;
    if (hasConcaveAft) {
      const R_b = rvMeta.rBase_m;
      const L_aft = rvMeta.xProfile[rvMeta.xProfile.length - 1];
      const dAft = rvMeta.aftCavityDepth_m;
      const pAft = rvMeta.aftCavityPower;
      const pGrid = Math.max(1.0, Math.min(pAft, 2.5));
      aftAperture = { active: true, xRim: L_aft, zCenter: 0.0, rBaseSq: R_b * R_b };
      const nRad = 6;
      for (let r = 0; r < nRad; r++) {
        const u1 = Math.pow(1.0 - ((1.0 - 0.04) * r) / nRad, 1.0 / pGrid);
        const u2 = Math.pow(1.0 - ((1.0 - 0.04) * (r + 1)) / nRad, 1.0 / pGrid);
        const x1 = L_aft - dAft * (1.0 - Math.pow(u1, pAft));
        const x2 = L_aft - dAft * (1.0 - Math.pow(u2, pAft));
        for (let c = 0; c < nPhi; c++) {
          const ph1 = phiEdges[c], ph2 = phiEdges[c + 1];
          const p1 = [x1, R_b * u1 * Math.sin(ph1), R_b * u1 * Math.cos(ph1)];
          const p2 = [x2, R_b * u2 * Math.sin(ph1), R_b * u2 * Math.cos(ph1)];
          const p3 = [x2, R_b * u2 * Math.sin(ph2), R_b * u2 * Math.cos(ph2)];
          const p4 = [x1, R_b * u1 * Math.sin(ph2), R_b * u1 * Math.cos(ph2)];
          pushQuadCell(p1, p2, p3, p4, 1, 0.5 * (ph1 + ph2));
        }
      }
    }

    if (hasFlaredCones) {
      const xProf = rvMeta.xProfile;
      const rProf = rvMeta.rProfile;
      let startIdx = 1;
      while (startIdx < rProf.length - 1 && rProf[startIdx] <= 1e-5) startIdx++;
      const idxSet = new Set();
      for (let k = 0; k < 12; k++) {
        idxSet.add(Math.round(startIdx + ((rProf.length - 1 - startIdx) * k) / 11));
      }
      if (Array.isArray(rvMeta.stageBreaksX)) {
        for (const xb of rvMeta.stageBreaksX) {
          let bestI = startIdx, bestD = Infinity;
          for (let i = startIdx; i < xProf.length; i++) {
            const d = Math.abs(xProf[i] - xb);
            if (d < bestD) { bestD = d; bestI = i; }
          }
          if (rProf[bestI] > 1e-5) idxSet.add(bestI);
        }
      }
      const iIdx = Array.from(idxSet).sort((a, b) => a - b);
      for (let r = 0; r < iIdx.length - 1; r++) {
        const x1 = xProf[iIdx[r]], rad1 = rProf[iIdx[r]];
        const x2 = xProf[iIdx[r + 1]], rad2 = rProf[iIdx[r + 1]];
        for (let c = 0; c < nPhi; c++) {
          const ph1 = phiEdges[c], ph2 = phiEdges[c + 1];
          const p1 = [x1, rad1 * Math.sin(ph1), rad1 * Math.cos(ph1)];
          const p2 = [x2, rad2 * Math.sin(ph1), rad2 * Math.cos(ph1)];
          const p3 = [x2, rad2 * Math.sin(ph2), rad2 * Math.cos(ph2)];
          const p4 = [x1, rad1 * Math.sin(ph2), rad1 * Math.cos(ph2)];
          pushQuadCell(p1, p2, p3, p4, 2, 0.5 * (ph1 + ph2));
        }
      }
    }

    const nCells = cellCenter.length;
    const pIdx = [];
    const qIdx = [];
    const rayPairs = [];

    for (let p = 0; p < nCells - 1; p++) {
      const rP = cellCenter[p], nP = cellNormal[p], gP = cellGroup[p];
      for (let q = p + 1; q < nCells; q++) {
        if (cellGroup[q] !== gP) continue;
        if (gP === 2) {
          const dotEr = cellEr[p][1] * cellEr[q][1] + cellEr[p][2] * cellEr[q][2];
          if (dotEr <= 0.25) continue;
        }
        const rQ = cellCenter[q], nQ = cellNormal[q];
        const dx = rQ[0] - rP[0], dy = rQ[1] - rP[1], dz = rQ[2] - rP[2];
        const dist = Math.hypot(dx, dy, dz);
        if (dist < 1e-4) continue;
        const dHat = [dx / dist, dy / dist, dz / dist];
        const cP = nP[0] * dHat[0] + nP[1] * dHat[1] + nP[2] * dHat[2];
        const cQ = -(nQ[0] * dHat[0] + nQ[1] * dHat[1] + nQ[2] * dHat[2]);
        if (cP > 0.015 && cQ > 0.015) {
          pIdx.push(p);
          qIdx.push(q);
          rayPairs.push({
            rP, rQ, nP, nQ, dist, dHat, cP, cQ,
            effA: Math.sqrt(cellArea[p] * cellArea[q] * cP * cQ)
          });
        }
      }
    }

    rayPairs.sort((a, b) => b.effA - a.effA);

    const nUndir = pIdx.length;
    const nDir = 2 * nUndir;
    const dirP = new Int32Array(nDir);
    const dirQ = new Int32Array(nDir);
    const dirGroup = new Int32Array(nDir);
    const pRx = new Float64Array(nDir), pRy = new Float64Array(nDir), pRz = new Float64Array(nDir);
    const pNx = new Float64Array(nDir), pNy = new Float64Array(nDir), pNz = new Float64Array(nDir);
    const qRx = new Float64Array(nDir), qRy = new Float64Array(nDir), qRz = new Float64Array(nDir);
    const qNx = new Float64Array(nDir), qNy = new Float64Array(nDir), qNz = new Float64Array(nDir);
    const pRad = new Float64Array(nDir), qRad = new Float64Array(nDir);

    for (let m = 0; m < nUndir; m++) {
      dirP[m] = pIdx[m]; dirQ[m] = qIdx[m];
      dirP[m + nUndir] = qIdx[m]; dirQ[m + nUndir] = pIdx[m];
    }
    for (let m = 0; m < nDir; m++) {
      const P = dirP[m], Q = dirQ[m];
      dirGroup[m] = cellGroup[P];
      pRx[m] = cellCenter[P][0]; pRy[m] = cellCenter[P][1]; pRz[m] = cellCenter[P][2];
      pNx[m] = cellNormal[P][0]; pNy[m] = cellNormal[P][1]; pNz[m] = cellNormal[P][2];
      qRx[m] = cellCenter[Q][0]; qRy[m] = cellCenter[Q][1]; qRz[m] = cellCenter[Q][2];
      qNx[m] = cellNormal[Q][0]; qNy[m] = cellNormal[Q][1]; qNz[m] = cellNormal[Q][2];
      pRad[m] = cellRadius[P]; qRad[m] = cellRadius[Q];
    }

    return {
      numPairs: nUndir,
      nCells,
      cellCenter,
      cellNormal,
      cellArea,
      cellVerts,
      cellRadius,
      pIdx,
      qIdx,
      dirP, dirQ, dirGroup, aftAperture,
      pRx, pRy, pRz, pNx, pNy, pNz,
      qRx, qRy, qRz, qNx, qNy, qNz,
      pRad, qRad,
      rayPairs: rayPairs.slice(0, 60)
    };
  }

  /**
   * Sutherland-Hodgman clipping of a planar polygon by a convex polygon (CCW about nHat).
   */
  function clipConvexPolygon(subject, clip, nHat) {
    let out = subject;
    const nC = clip.length;
    for (let e = 0; e < nC; e++) {
      if (out.length === 0) return out;
      const a = clip[e];
      const b = clip[(e + 1) % nC];
      const ebx = b[0] - a[0], eby = b[1] - a[1], ebz = b[2] - a[2];
      // inward = nHat x (b - a)
      const inx = nHat[1] * ebz - nHat[2] * eby;
      const iny = nHat[2] * ebx - nHat[0] * ebz;
      const inz = nHat[0] * eby - nHat[1] * ebx;

      const inp = out;
      const nI = inp.length;
      const dist = new Float64Array(nI);
      for (let i = 0; i < nI; i++) {
        dist[i] = (inp[i][0] - a[0]) * inx + (inp[i][1] - a[1]) * iny + (inp[i][2] - a[2]) * inz;
      }
      out = [];
      for (let i = 0; i < nI; i++) {
        const j = (i + 1) % nI;
        const di = dist[i], dj = dist[j];
        if (di >= 0.0) {
          out.push(inp[i]);
        }
        if ((di >= 0.0) !== (dj >= 0.0)) {
          const t = di / (di - dj);
          out.push([
            inp[i][0] + t * (inp[j][0] - inp[i][0]),
            inp[i][1] + t * (inp[j][1] - inp[i][1]),
            inp[i][2] + t * (inp[j][2] - inp[i][2])
          ]);
        }
      }
    }
    return out;
  }

  // Zero-allocation scratch buffers for Sutherland-Hodgman quad-against-quad clipping (max 8 vertices = 24 floats)
  const _projBuf = new Float64Array(12);
  const _clipBufA = new Float64Array(48);
  const _clipBufB = new Float64Array(48);
  const _clipDist = new Float64Array(16);

  function clipQuadIntoFlatBuffer(proj12, qv, nQx, nQy, nQz) {
    // Fast 3D AABB reject before Sutherland-Hodgman edge clipping
    let pMinX = proj12[0], pMaxX = proj12[0];
    let pMinY = proj12[1], pMaxY = proj12[1];
    let pMinZ = proj12[2], pMaxZ = proj12[2];
    for (let i = 3; i < 12; i += 3) {
      const x = proj12[i], y = proj12[i + 1], z = proj12[i + 2];
      if (x < pMinX) pMinX = x; if (x > pMaxX) pMaxX = x;
      if (y < pMinY) pMinY = y; if (y > pMaxY) pMaxY = y;
      if (z < pMinZ) pMinZ = z; if (z > pMaxZ) pMaxZ = z;
    }
    const q0 = qv[0], q1 = qv[1], q2 = qv[2], q3 = qv[3];
    const qMinX = Math.min(q0[0], q1[0], q2[0], q3[0]) - 1e-9;
    if (pMaxX < qMinX) return 0;
    const qMaxX = Math.max(q0[0], q1[0], q2[0], q3[0]) + 1e-9;
    if (pMinX > qMaxX) return 0;
    const qMinY = Math.min(q0[1], q1[1], q2[1], q3[1]) - 1e-9;
    if (pMaxY < qMinY) return 0;
    const qMaxY = Math.max(q0[1], q1[1], q2[1], q3[1]) + 1e-9;
    if (pMinY > qMaxY) return 0;
    const qMinZ = Math.min(q0[2], q1[2], q2[2], q3[2]) - 1e-9;
    if (pMaxZ < qMinZ) return 0;
    const qMaxZ = Math.max(q0[2], q1[2], q2[2], q3[2]) + 1e-9;
    if (pMinZ > qMaxZ) return 0;

    let src = _clipBufA;
    let dst = _clipBufB;
    for (let i = 0; i < 12; i++) src[i] = proj12[i];
    let nSrc = 4;

    for (let e = 0; e < 4; e++) {
      if (nSrc === 0) return 0;
      const a = qv[e];
      const b = qv[(e + 1) & 3];
      const ebx = b[0] - a[0], eby = b[1] - a[1], ebz = b[2] - a[2];
      const inx = nQy * ebz - nQz * eby;
      const iny = nQz * ebx - nQx * ebz;
      const inz = nQx * eby - nQy * ebx;

      for (let i = 0; i < nSrc; i++) {
        const idx = i * 3;
        _clipDist[i] = (src[idx] - a[0]) * inx + (src[idx + 1] - a[1]) * iny + (src[idx + 2] - a[2]) * inz;
      }

      let nDst = 0;
      for (let i = 0; i < nSrc; i++) {
        const j = i + 1 === nSrc ? 0 : i + 1;
        const di = _clipDist[i], dj = _clipDist[j];
        const i3 = i * 3;
        if (di >= 0.0) {
          const o3 = nDst * 3;
          dst[o3] = src[i3]; dst[o3 + 1] = src[i3 + 1]; dst[o3 + 2] = src[i3 + 2];
          nDst++;
        }
        if ((di >= 0.0) !== (dj >= 0.0)) {
          const t = di / (di - dj);
          const j3 = j * 3;
          const o3 = nDst * 3;
          dst[o3] = src[i3] + t * (src[j3] - src[i3]);
          dst[o3 + 1] = src[i3 + 1] + t * (src[j3 + 1] - src[i3 + 1]);
          dst[o3 + 2] = src[i3 + 2] + t * (src[j3 + 2] - src[i3 + 2]);
          nDst++;
        }
      }
      const tmp = src; src = dst; dst = tmp;
      nSrc = nDst;
    }

    if (src !== _clipBufA) {
      for (let i = 0; i < nSrc * 3; i++) _clipBufA[i] = src[i];
    }
    return nSrc;
  }

  /**
   * Normalize polarization channel identifier (IEEE Backscatter Alignment convention).
   * Supported: Mean, HH, VV, HV, VH, PP (aliases RL, LR), OP (alias RR), LL.
   */
  function normalizePolarization(pol) {
    const s = String(pol == null ? "Mean" : pol).trim().toUpperCase();
    switch (s) {
      case "":
      case "MEAN":
      case "SCALAR":
      case "UNPOL":
        return "Mean";
      case "HH":
      case "VV":
      case "HV":
      case "VH":
      case "LL":
        return s;
      case "PP":
      case "RL":
      case "LR":
        return "PP";
      case "OP":
      case "RR":
        return "OP";
      default:
        throw new RangeError(
          `Unsupported polarization "${pol}". Expected Mean, HH, VV, HV, VH, PP (RL/LR), OP (RR), or LL.`
        );
    }
  }

  /**
   * Singularity-free IEEE BSA polarization basis (e_H, e_V) for look vector u = [ux, uy, uz]
   * pointing FROM target TO radar (propagation k_i = -u, so e_H x e_V = k_i).
   */
  function computePolarizationBasis(ux, uy, uz, azRadOpt) {
    let eHx, eHy;
    if (typeof azRadOpt === "number" && Number.isFinite(azRadOpt)) {
      eHx = Math.sin(azRadOpt);
      eHy = Math.cos(azRadOpt);
    } else {
      const rhoXY = Math.hypot(ux, uy);
      if (rhoXY > 1e-15) {
        eHx = uy / rhoXY;
        eHy = -ux / rhoXY;
      } else {
        eHx = 0.0;
        eHy = 1.0;
      }
    }
    const eHz = 0.0;
    const eVx = eHy * uz;
    const eVy = -eHx * uz;
    const eVz = eHx * uy - eHy * ux;
    return [eHx, eHy, eHz, eVx, eVy, eVz];
  }

  /**
   * Collect clipped GO-PO overlap polygons and 3D vector dyadic polarization weights
   * for a look direction u = [ux, uy, uz] (frequency-independent!).
   */
  function collectDoubleBouncePatches(cp, ux, uy, uz, azRadOpt) {
    if (!cp || cp.numPairs === 0) return [];
    const { cellNormal, cellVerts, dirP, dirQ, dirGroup, aftAperture, pRx, pRy, pRz, pNx, pNy, pNz, qRx, qRy, qRz, qNx, qNy, qNz, pRad, qRad } = cp;
    const hasAp = Boolean(aftAperture && aftAperture.active && dirGroup);
    const xRim = hasAp ? aftAperture.xRim : 0.0;
    const zApCtr = hasAp ? aftAperture.zCenter : 0.0;
    const rBaseSq = hasAp ? aftAperture.rBaseSq : 0.0;
    const nDir = dirP.length;
    const kix = -ux, kiy = -uy, kiz = -uz;
    const [eHx, eHy, eHz, eVx, eVy, eVz] = computePolarizationBasis(ux, uy, uz, azRadOpt);
    const patches = [];

    for (let m = 0; m < nDir; m++) {
      const isAftCav = hasAp && dirGroup[m] === 1;
      if (isAftCav) {
        if (ux <= 1e-7) continue;
        const dxP = xRim - pRx[m];
        if (dxP > 1e-7) {
          const tP = dxP / ux;
          const yP = pRy[m] + tP * uy;
          const zP = pRz[m] - zApCtr + tP * uz;
          if (yP * yP + zP * zP > rBaseSq) continue;
        }
      }
      const cosP = pNx[m] * ux + pNy[m] * uy + pNz[m] * uz;
      if (cosP <= 1e-6) continue;

      const krx = kix + 2.0 * cosP * pNx[m];
      const kry = kiy + 2.0 * cosP * pNy[m];
      const krz = kiz + 2.0 * cosP * pNz[m];

      const nQkr = qNx[m] * krx + qNy[m] * kry + qNz[m] * krz;
      if (nQkr >= -1e-6) continue;

      const dPQx = qRx[m] - pRx[m], dPQy = qRy[m] - pRy[m], dPQz = qRz[m] - pRz[m];
      // Exact ray-cylinder necessary condition: |(rQ - rP) x k_r| <= pRad + qRad
      const crx = dPQy * krz - dPQz * kry;
      const cry = dPQz * krx - dPQx * krz;
      const crz = dPQx * kry - dPQy * krx;
      const sumRad = pRad[m] + qRad[m];
      if (crx * crx + cry * cry + crz * crz > sumRad * sumRad) continue;

      const tHit = (dPQx * qNx[m] + dPQy * qNy[m] + dPQz * qNz[m]) / nQkr;
      if (tHit <= 0.0) continue;

      const pv = cellVerts[dirP[m]];
      const nQx = qNx[m], nQy = qNy[m], nQz = qNz[m];
      const qrx = qRx[m], qry = qRy[m], qrz = qRz[m];
      let behind = false;
      for (let v = 0; v < 4; v++) {
        const pvv = pv[v];
        const tv = ((qrx - pvv[0]) * nQx + (qry - pvv[1]) * nQy + (qrz - pvv[2]) * nQz) / nQkr;
        if (tv < -1e-9) { behind = true; break; }
        const v3 = v * 3;
        _projBuf[v3] = pvv[0] + tv * krx;
        _projBuf[v3 + 1] = pvv[1] + tv * kry;
        _projBuf[v3 + 2] = pvv[2] + tv * krz;
      }
      if (behind) continue;

      const nV = clipQuadIntoFlatBuffer(_projBuf, cellVerts[dirQ[m]], nQx, nQy, nQz);
      if (nV < 3) continue;

      // Signed area & area-weighted centroid of convex polygon in _clipBufA
      let avx = 0.0, avy = 0.0, avz = 0.0;
      let sumA = 0.0, ox = 0.0, oy = 0.0, oz = 0.0;
      const p0x = _clipBufA[0], p0y = _clipBufA[1], p0z = _clipBufA[2];
      for (let i = 1; i < nV - 1; i++) {
        const i3 = i * 3, j3 = (i + 1) * 3;
        const pix = _clipBufA[i3], piy = _clipBufA[i3 + 1], piz = _clipBufA[i3 + 2];
        const pjx = _clipBufA[j3], pjy = _clipBufA[j3 + 1], pjz = _clipBufA[j3 + 2];
        const ax = pix - p0x, ay = piy - p0y, az = piz - p0z;
        const bx = pjx - p0x, by = pjy - p0y, bz = pjz - p0z;
        const cx = ay * bz - az * by;
        const cy = az * bx - ax * bz;
        const cz = ax * by - ay * bx;
        avx += cx; avy += cy; avz += cz;
        const tA = 0.5 * Math.hypot(cx, cy, cz);
        sumA += tA;
        ox += ((p0x + pix + pjx) / 3.0) * tA;
        oy += ((p0y + piy + pjy) / 3.0) * tA;
        oz += ((p0z + piz + pjz) / 3.0) * tA;
      }
      const signed2A = avx * nQx + avy * nQy + avz * nQz;
      if (Math.abs(signed2A) < 1e-12) continue;

      const polyFlat = new Float64Array(nV * 3);
      if (signed2A >= 0.0) {
        for (let k = 0; k < nV * 3; k++) polyFlat[k] = _clipBufA[k];
      } else {
        for (let v = 0; v < nV; v++) {
          const src3 = (nV - 1 - v) * 3, dst3 = v * 3;
          polyFlat[dst3] = _clipBufA[src3];
          polyFlat[dst3 + 1] = _clipBufA[src3 + 1];
          polyFlat[dst3 + 2] = _clipBufA[src3 + 2];
        }
      }

      let Lmax = 0.0;
      for (let i = 0; i < nV; i++) {
        const i3 = i * 3, j3 = ((i + 1) % nV) * 3;
        const el = Math.hypot(polyFlat[j3] - polyFlat[i3], polyFlat[j3 + 1] - polyFlat[i3 + 1], polyFlat[j3 + 2] - polyFlat[i3 + 2]);
        if (el > Lmax) Lmax = el;
      }

      const invA = 1.0 / Math.max(sumA, 1e-15);
      const rOx = ox * invA, rOy = oy * invA, rOz = oz * invA;
      if (isAftCav) {
        const dxO = xRim - rOx;
        if (dxO > 1e-7) {
          const tO = dxO / ux;
          const yO = rOy + tO * uy;
          const zO = rOz - zApCtr + tO * uz;
          if (yO * yO + zO * zO > rBaseSq) continue;
        }
      }
      const cPath = (kix - krx) * pRx[m] + (kiy - kry) * pRy[m] + (kiz - krz) * pRz[m];
      const dkx = ux - krx, dky = uy - kry, dkz = uz - krz;
      const nQdk = nQx * dkx + nQy * dky + nQz * dkz;
      const dk2 = dkx * dkx + dky * dky + dkz * dkz;
      const uwt2 = Math.max(dk2 - nQdk * nQdk, 0.0);
      const Lhalf = 0.5 * (-dkx * rOx - dky * rOy - dkz * rOz + cPath);

      // Full 3D vector PEC reflection off P and PO radiation off Q (IEEE BSA)
      const nPx = pNx[m], nPy = pNy[m], nPz = pNz[m];
      const npEh = 2.0 * (nPx * eHx + nPy * eHy + nPz * eHz);
      const erHx = npEh * nPx - eHx, erHy = npEh * nPy - eHy, erHz = npEh * nPz - eHz;
      const npEv = 2.0 * (nPx * eVx + nPy * eVy + nPz * eVz);
      const erVx = npEv * nPx - eVx, erVy = npEv * nPy - eVy, erVz = npEv * nPz - eVz;

      const eH_kr = eHx * krx + eHy * kry + eHz * krz;
      const eV_kr = eVx * krx + eVy * kry + eVz * krz;
      const nq_erH = nQx * erHx + nQy * erHy + nQz * erHz;
      const nq_erV = nQx * erVx + nQy * erVy + nQz * erVz;

      const wHH = 2.0 * (nq_erH * eH_kr - nQkr * (eHx * erHx + eHy * erHy + eHz * erHz));
      const wVV = 2.0 * (nq_erV * eV_kr - nQkr * (eVx * erVx + eVy * erVy + eVz * erVz));
      const wHV = 2.0 * (nq_erV * eH_kr - nQkr * (eHx * erVx + eHy * erVy + eHz * erVz));
      const wVH = 2.0 * (nq_erH * eV_kr - nQkr * (eVx * erHx + eVy * erHy + eVz * erHz));
      const wCross = 0.5 * (wHV + wVH);

      patches.push({
        polyFlat, nV, nQ: cellNormal[dirQ[m]], area: sumA, rOx, rOy, rOz, Lmax,
        dkx, dky, dkz, nQdk, uwt2, cPath, Lhalf,
        wHH, wVV, wCross
      });
    }
    return patches;
  }

  /**
   * Evaluate exact Gordon contour integral over pre-clipped GO-PO patches at wavenumber k0
   * for scalar Mean and all 4 Sinclair channels simultaneously.
   */
  function evaluatePatchesAtK(patches, k0, wantContrib = false) {
    const Gamma = -1.0;
    const prefScale = k0 / (2.0 * Math.sqrt(Math.PI));
    let totRe = 0.0, totIm = 0.0;
    let hhRe = 0.0, hhIm = 0.0;
    let vvRe = 0.0, vvIm = 0.0;
    let hvRe = 0.0, hvIm = 0.0;
    const contrib = wantContrib ? [] : null;

    for (let p = 0; p < patches.length; p++) {
      const pt = patches[p];
      const wx = k0 * pt.dkx, wy = k0 * pt.dky, wz = k0 * pt.dkz;
      const wt2 = k0 * k0 * pt.uwt2;
      let IRe = 0.0, IIm = 0.0;

      if (wt2 * pt.Lmax * pt.Lmax < 0.04) {
        const ph = pt.rOx * wx + pt.rOy * wy + pt.rOz * wz;
        IRe = pt.area * Math.cos(ph);
        IIm = pt.area * Math.sin(ph);
      } else {
        const pf = pt.polyFlat;
        const nQ = pt.nQ;
        const nV = pt.nV;
        let sumRe = 0.0, sumIm = 0.0;
        for (let e = 0; e < nV; e++) {
          const a3 = e * 3, b3 = (e + 1 === nV ? 0 : e + 1) * 3;
          const ax = pf[a3], ay = pf[a3 + 1], az = pf[a3 + 2];
          const bx = pf[b3], by = pf[b3 + 1], bz = pf[b3 + 2];
          const ex = bx - ax, ey = by - ay, ez = bz - az;
          const inNx = ey * nQ[2] - ez * nQ[1];
          const inNy = ez * nQ[0] - ex * nQ[2];
          const inNz = ex * nQ[1] - ey * nQ[0];
          const wm = inNx * wx + inNy * wy + inNz * wz;
          const x = 0.5 * (ex * wx + ey * wy + ez * wz);
          const wr = 0.5 * ((ax + bx) * wx + (ay + by) * wy + (az + bz) * wz);
          const coef = wm * (Math.abs(x) > 1e-12 ? Math.sin(x) / x : 1.0);
          sumRe += coef * Math.cos(wr);
          sumIm += coef * Math.sin(wr);
        }
        IRe = sumIm / wt2;
        IIm = -sumRe / wt2;
      }

      const cosC = Math.cos(-k0 * pt.cPath);
      const sinC = Math.sin(-k0 * pt.cPath);
      const pRe = IRe * cosC - IIm * sinC;
      const pIm = IRe * sinC + IIm * cosC;
      const baseRe = -prefScale * pIm;
      const baseIm = prefScale * pRe;

      const wMean = pt.nQdk * Gamma;
      const ampRe = baseRe * wMean;
      const ampIm = baseIm * wMean;
      const aHHRe = baseRe * pt.wHH, aHHIm = baseIm * pt.wHH;
      const aVVRe = baseRe * pt.wVV, aVVIm = baseIm * pt.wVV;
      const aHVRe = baseRe * pt.wCross, aHVIm = baseIm * pt.wCross;

      totRe += ampRe; totIm += ampIm;
      hhRe += aHHRe;  hhIm += aHHIm;
      vvRe += aVVRe;  vvIm += aVVIm;
      hvRe += aHVRe;  hvIm += aHVIm;
      if (wantContrib) {
        contrib.push([ampRe, ampIm, pt.Lhalf, aHHRe, aHHIm, aVVRe, aVVIm, aHVRe, aHVIm]);
      }
    }
    return {
      re: totRe, im: totIm,
      hhRe, hhIm, vvRe, vvIm, hvRe, hvIm, vhRe: hvRe, vhIm: hvIm,
      contrib
    };
  }

  const _ZERO_MB = Object.freeze({
    re: 0.0, im: 0.0,
    hhRe: 0.0, hhIm: 0.0,
    vvRe: 0.0, vvIm: 0.0,
    hvRe: 0.0, hvIm: 0.0,
    vhRe: 0.0, vhIm: 0.0,
    contrib: []
  });

  /**
   * Exact GO-PO (SBR-style) monostatic two-bounce evaluator for one look direction u = [ux, uy, uz].
   * Matches +waverider/+rcs/doubleBounceGOPO.m.
   */
  function doubleBounceGOPO(cp, ux, uy, uz, k0, wantContrib = false, azRadOpt) {
    const patches = collectDoubleBouncePatches(cp, ux, uy, uz, azRadOpt);
    if (patches.length === 0) {
      return {
        re: 0.0, im: 0.0,
        hhRe: 0.0, hhIm: 0.0, vvRe: 0.0, vvIm: 0.0, hvRe: 0.0, hvIm: 0.0, vhRe: 0.0, vhIm: 0.0,
        contrib: wantContrib ? [] : null
      };
    }
    return evaluatePatchesAtK(patches, k0, wantContrib);
  }

  const _ZERO_PTD = Object.freeze({
    re: 0.0, im: 0.0,
    hhRe: 0.0, hhIm: 0.0,
    vvRe: 0.0, vvIm: 0.0,
    hvRe: 0.0, hvIm: 0.0
  });

  /**
   * Extract exterior sharp base-perimeter wedge edges shared between skin (facePatch != 3)
   * and aft base (facePatch == 3) for polarimetric Ufimtsev/Mitzner PTD fringe diffraction.
   */
  function buildRimPTDEdges(vertices, faces, facePatch, faceNormals, hasConcaveAft = false) {
    if (!facePatch || faces.length === 0) return { count: 0, hasConcaveAft: false };
    const edgeMap = new Map();
    for (let f = 0; f < faces.length; f++) {
      const [i0, i1, i2] = faces[f];
      const isBase = facePatch[f] === 3;
      const pairs = [[i0, i1, i2], [i1, i2, i0], [i2, i0, i1]];
      for (let e = 0; e < 3; e++) {
        const u = pairs[e][0], v = pairs[e][1], ic = pairs[e][2];
        const key = u < v ? `${u}_${v}` : `${v}_${u}`;
        let rec = edgeMap.get(key);
        if (!rec) {
          rec = { ia: u, ib: v, fSkin: -1, icSkin: -1, fBase: -1, icBase: -1 };
          edgeMap.set(key, rec);
        }
        if (isBase) {
          rec.fBase = f;
          rec.icBase = ic;
        } else {
          rec.fSkin = f;
          rec.icSkin = ic;
        }
      }
    }

    const edges = [];
    for (const rec of edgeMap.values()) {
      if (rec.fSkin < 0 || rec.fBase < 0) continue;
      const va = vertices[rec.ia], vb = vertices[rec.ib];
      const vc1 = vertices[rec.icSkin], vc2 = vertices[rec.icBase];
      const rex = 0.5 * (va[0] + vb[0]), rey = 0.5 * (va[1] + vb[1]), rez = 0.5 * (va[2] + vb[2]);
      const dex = vb[0] - va[0], dey = vb[1] - va[1], dez = vb[2] - va[2];
      const len = Math.hypot(dex, dey, dez);
      if (len < 1e-9) continue;

      const n1 = faceNormals[rec.fSkin];
      const n2 = faceNormals[rec.fBase];
      // Strict exterior convex wedge test: Face 2's interior vertex must lie behind Face 1's plane
      const convCheck = (vc2[0] - rex) * n1[0] + (vc2[1] - rey) * n1[1] + (vc2[2] - rez) * n1[2];
      if (convCheck >= -1e-9) continue;

      const cos12 = Math.max(-1.0, Math.min(1.0, n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]));
      const nWedge = 1.0 + Math.acos(cos12) / Math.PI;
      if (nWedge <= 1.02) continue;

      const dx = dex / len, dy = dey / len, dz = dez / len;
      // q1 = n1 x dHat lies in Face 1 perpendicular to the edge
      const q1x = n1[1] * dz - n1[2] * dy;
      const q1y = n1[2] * dx - n1[0] * dz;
      const q1z = n1[0] * dy - n1[1] * dx;
      const sgn = ((vc1[0] - rex) * q1x + (vc1[1] - rey) * q1y + (vc1[2] - rez) * q1z) >= 0.0 ? 1.0 : -1.0;

      // Right-handed wedge basis: tau1 x n1 = tHat, Face 1 at phi = 0, exterior in (0, n*pi)
      const tau1x = sgn * q1x, tau1y = sgn * q1y, tau1z = sgn * q1z;
      const tx = sgn * dx, ty = sgn * dy, tz = sgn * dz;

      const Xwedge = -1.0 / (nWedge * Math.tan(Math.PI / (2.0 * nWedge)));
      const sinPiOverN = Math.sin(Math.PI / nWedge) / nWedge;
      const cosPiOverN = Math.cos(Math.PI / nWedge);
      const poleRes = -0.5 / (nWedge * Math.tan(Math.PI / nWedge));

      edges.push({
        rex, rey, rez, len,
        tx, ty, tz,
        tau1x, tau1y, tau1z,
        n1x: n1[0], n1y: n1[1], n1z: n1[2],
        n2x: n2[0], n2y: n2[1], n2z: n2[2],
        nWedge, nPi: nWedge * Math.PI,
        Xwedge, sinPiOverN, cosPiOverN, poleRes
      });
    }

    return { count: edges.length, hasConcaveAft, edges };
  }

  /**
   * Singularity-free polarimetric Ufimtsev/Mitzner ILDC monostatic PTD fringe evaluator
   * across all sharp aft base-rim wedge segments (and optional HRRP range bin deposition).
   */
  function evaluateRimPTD(rimData, ux, uy, uz, k0, azRadOpt, splatOpts = null) {
    if (!rimData || rimData.count === 0) return _ZERO_PTD;
    const edges = rimData.edges;
    const hasConcaveAft = rimData.hasConcaveAft;
    const [eHx, eHy, eHz, eVx, eVy, eVz] = computePolarizationBasis(ux, uy, uz, azRadOpt);
    const SQRT_PI = Math.sqrt(Math.PI);
    const HALF_PI = 0.5 * Math.PI;

    let re = 0.0, im = 0.0;
    let hhRe = 0.0, hhIm = 0.0;
    let vvRe = 0.0, vvIm = 0.0;
    let hvRe = 0.0, hvIm = 0.0;

    const doSplat = splatOpts !== null;
    const polTag = doSplat ? splatOpts.polTag : "";
    const accRe = doSplat ? splatOpts.accRe : null;
    const accIm = doSplat ? splatOpts.accIm : null;
    const xRefProj = doSplat ? ux * splatOpts.xRef : 0.0;
    const rLim = doSplat ? splatOpts.rLim : 0.0;
    const invDBin = doSplat ? splatOpts.invDBin : 0.0;
    const nR = doSplat ? splatOpts.nR : 0;

    for (let i = 0; i < edges.length; i++) {
      const ed = edges[i];
      const lit1 = (ed.n1x * ux + ed.n1y * uy + ed.n1z * uz) > 1e-9;
      const lit2 = ((ed.n2x * ux + ed.n2y * uy + ed.n2z * uz) > 1e-9) && (!hasConcaveAft || ux > 1e-7);
      if (!lit1 && !lit2) continue;

      const tu = ed.tx * ux + ed.ty * uy + ed.tz * uz;
      const sinB2 = 1.0 - tu * tu;
      if (sinB2 < 1e-6) continue;
      const sinB = Math.sqrt(sinB2);
      const invSinB = 1.0 / sinB;

      const upx = (ux - tu * ed.tx) * invSinB;
      const upy = (uy - tu * ed.ty) * invSinB;
      const upz = (uz - tu * ed.tz) * invSinB;

      const xw = upx * ed.tau1x + upy * ed.tau1y + upz * ed.tau1z;
      const yw = upx * ed.n1x + upy * ed.n1y + upz * ed.n1z;
      let phi = Math.atan2(yw, xw);
      if (phi < 0.0) phi += 2.0 * Math.PI;
      if (phi <= 1e-5 || phi >= ed.nPi - 1e-5) continue;

      const phi2 = ed.nPi - phi;
      let Yf;
      if (Math.abs(phi - HALF_PI) < 1e-5) {
        Yf = ed.poleRes - (lit2 && phi2 > 0.0 && phi2 < Math.PI ? -0.5 * Math.tan(phi2) : 0.0);
      } else if (lit2 && Math.abs(phi2 - HALF_PI) < 1e-5) {
        Yf = ed.poleRes - (phi > 0.0 && phi < Math.PI ? -0.5 * Math.tan(phi) : 0.0);
      } else {
        const symSum = ed.sinPiOverN / (ed.cosPiOverN - Math.cos((2.0 * phi) / ed.nWedge));
        const po1 = phi > 0.0 && phi < Math.PI ? -0.5 * Math.tan(phi) : 0.0;
        const po2 = lit2 && phi2 > 0.0 && phi2 < Math.PI ? -0.5 * Math.tan(phi2) : 0.0;
        Yf = symSum - po1 - po2;
      }

      const fs = Yf - ed.Xwedge;
      const fh = Yf + ed.Xwedge;

      const hPar = (eHx * ed.tx + eHy * ed.ty + eHz * ed.tz) * invSinB;
      const vPar = (eVx * ed.tx + eVy * ed.ty + eVz * ed.tz) * invSinB;
      const hPar2 = hPar * hPar;
      const vPar2 = vPar * vPar;

      const wHH = fs * hPar2 + fh * vPar2;
      const wVV = fs * vPar2 + fh * hPar2;
      const wHV = (fs - fh) * hPar * vPar;

      const xSinc = k0 * ed.len * tu;
      const sincVal = Math.abs(xSinc) > 1e-8 ? Math.sin(xSinc) / xSinc : 1.0 - (xSinc * xSinc) / 6.0;
      const ampE = ((ed.len * sinB) / SQRT_PI) * sincVal;
      const uDotRe = ed.rex * ux + ed.rey * uy + ed.rez * uz;
      const ph = 2.0 * k0 * uDotRe;
      const eRe = ampE * Math.cos(ph);
      const eIm = ampE * Math.sin(ph);

      const mRe = eRe * Yf,  mIm = eIm * Yf;
      const hRe = eRe * wHH, hIm = eIm * wHH;
      const vRe = eRe * wVV, vIm = eIm * wVV;
      const cRe = eRe * wHV, cIm = eIm * wHV;

      re += mRe;   im += mIm;
      hhRe += hRe; hhIm += hIm;
      vvRe += vRe; vvIm += vIm;
      hvRe += cRe; hvIm += cIm;

      if (doSplat) {
        let sRe = 0.0, sIm = 0.0;
        switch (polTag) {
          case "Mean":
          case "PP":
            sRe = mRe; sIm = mIm;
            break;
          case "HH":
            sRe = hRe; sIm = hIm;
            break;
          case "VV":
            sRe = vRe; sIm = vIm;
            break;
          case "HV":
          case "VH":
            sRe = cRe; sIm = cIm;
            break;
          case "OP":
            sRe = 0.5 * (hRe - vRe) + cIm;
            sIm = 0.5 * (hIm - vIm) - cRe;
            break;
          case "LL":
            sRe = 0.5 * (hRe - vRe) - cIm;
            sIm = 0.5 * (hIm - vIm) + cRe;
            break;
        }
        const RS = -(uDotRe - xRefProj);
        const fBin = (RS + rLim) * invDBin;
        const i0 = Math.floor(fBin);
        if (i0 >= 0 && i0 < nR - 1) {
          const t = fBin - i0;
          accRe[i0] += sRe * (1.0 - t);
          accIm[i0] += sIm * (1.0 - t);
          accRe[i0 + 1] += sRe * t;
          accIm[i0 + 1] += sIm * t;
        }
      }
    }

    return { re, im, hhRe, hhIm, vvRe, vvIm, hvRe, hvIm };
  }

  /**
   * Assemble single-bounce PO + optional PTD rim fringe + GO-PO double-bounce into the
   * requested polarization channel plus the complete linear/circular polarimetric dictionary.
   */
  function assemblePolarimetricPoint(rePO, imPO, mb, polTag, noiseFloor = 1e-5, ptd = _ZERO_PTD) {
    const sMeanRe = rePO + ptd.re + mb.re,     sMeanIm = imPO + ptd.im + mb.im;
    const sHHRe   = rePO + ptd.hhRe + mb.hhRe, sHHIm   = imPO + ptd.hhIm + mb.hhIm;
    const sVVRe   = rePO + ptd.vvRe + mb.vvRe, sVVIm   = imPO + ptd.vvIm + mb.vvIm;
    const sHVRe   = ptd.hvRe + mb.hvRe,        sHVIm   = ptd.hvIm + mb.hvIm;
    const sPPRe   = 0.5 * (sHHRe + sVVRe), sPPIm = 0.5 * (sHHIm + sVVIm);
    const sOPRe   = 0.5 * (sHHRe - sVVRe) + sHVIm, sOPIm = 0.5 * (sHHIm - sVVIm) - sHVRe;
    const sLLRe   = 0.5 * (sHHRe - sVVRe) - sHVIm, sLLIm = 0.5 * (sHHIm - sVVIm) + sHVRe;

    let reSingle = 0.0, imSingle = 0.0, reMB = 0.0, imMB = 0.0;
    switch (polTag) {
      case "Mean":
        reSingle = rePO + ptd.re; imSingle = imPO + ptd.im; reMB = mb.re; imMB = mb.im;
        break;
      case "HH":
        reSingle = rePO + ptd.hhRe; imSingle = imPO + ptd.hhIm; reMB = mb.hhRe; imMB = mb.hhIm;
        break;
      case "VV":
        reSingle = rePO + ptd.vvRe; imSingle = imPO + ptd.vvIm; reMB = mb.vvRe; imMB = mb.vvIm;
        break;
      case "HV":
      case "VH":
        reSingle = ptd.hvRe; imSingle = ptd.hvIm; reMB = mb.hvRe; imMB = mb.hvIm;
        break;
      case "PP":
        reSingle = rePO + 0.5 * (ptd.hhRe + ptd.vvRe);
        imSingle = imPO + 0.5 * (ptd.hhIm + ptd.vvIm);
        reMB = 0.5 * (mb.hhRe + mb.vvRe); imMB = 0.5 * (mb.hhIm + mb.vvIm);
        break;
      case "OP":
        reSingle = 0.5 * (ptd.hhRe - ptd.vvRe) + ptd.hvIm;
        imSingle = 0.5 * (ptd.hhIm - ptd.vvIm) - ptd.hvRe;
        reMB = 0.5 * (mb.hhRe - mb.vvRe) + mb.hvIm;
        imMB = 0.5 * (mb.hhIm - mb.vvIm) - mb.hvRe;
        break;
      case "LL":
        reSingle = 0.5 * (ptd.hhRe - ptd.vvRe) - ptd.hvIm;
        imSingle = 0.5 * (ptd.hhIm - ptd.vvIm) + ptd.hvRe;
        reMB = 0.5 * (mb.hhRe - mb.vvRe) - mb.hvIm;
        imMB = 0.5 * (mb.hhIm - mb.vvIm) + mb.hvRe;
        break;
    }

    const reTotal = reSingle + reMB;
    const imTotal = imSingle + imMB;
    const sigSingle = reSingle * reSingle + imSingle * imSingle + noiseFloor;
    const sigDouble = reMB * reMB + imMB * imMB + noiseFloor;
    const sigTotal = reTotal * reTotal + imTotal * imTotal + noiseFloor;

    const m2HH = sHHRe * sHHRe + sHHIm * sHHIm;
    const m2VV = sVVRe * sVVRe + sVVIm * sVVIm;
    const m2HV = sHVRe * sHVRe + sHVIm * sHVIm;
    const m2PP = sPPRe * sPPRe + sPPIm * sPPIm;
    const m2OP = sOPRe * sOPRe + sOPIm * sOPIm;
    const m2LL = sLLRe * sLLRe + sLLIm * sLLIm;
    const m2Mean = sMeanRe * sMeanRe + sMeanIm * sMeanIm;
    const span_m2 = m2HH + m2VV + 2.0 * m2HV;

    return {
      polarization: polTag,
      reSingle, imSingle, reMB, imMB, reTotal, imTotal,
      single_m2: sigSingle,
      double_m2: sigDouble,
      total_m2: sigTotal,
      single_dBsm: 10.0 * Math.log10(sigSingle),
      double_dBsm: 10.0 * Math.log10(sigDouble),
      total_dBsm: 10.0 * Math.log10(sigTotal),
      channels: {
        Mean: { re: sMeanRe, im: sMeanIm, m2: m2Mean, dBsm: 10.0 * Math.log10(m2Mean + noiseFloor) },
        HH:   { re: sHHRe,   im: sHHIm,   m2: m2HH,   dBsm: 10.0 * Math.log10(m2HH + noiseFloor) },
        VV:   { re: sVVRe,   im: sVVIm,   m2: m2VV,   dBsm: 10.0 * Math.log10(m2VV + noiseFloor) },
        HV:   { re: sHVRe,   im: sHVIm,   m2: m2HV,   dBsm: 10.0 * Math.log10(m2HV + noiseFloor) },
        VH:   { re: sHVRe,   im: sHVIm,   m2: m2HV,   dBsm: 10.0 * Math.log10(m2HV + noiseFloor) },
        PP:   { re: sPPRe,   im: sPPIm,   m2: m2PP,   dBsm: 10.0 * Math.log10(m2PP + noiseFloor) },
        OP:   { re: sOPRe,   im: sOPIm,   m2: m2OP,   dBsm: 10.0 * Math.log10(m2OP + noiseFloor) },
        LL:   { re: sLLRe,   im: sLLIm,   m2: m2LL,   dBsm: 10.0 * Math.log10(m2LL + noiseFloor) },
        span_m2
      }
    };
  }

  /**
   * Evaluate coherent monostatic RCS across multiple frequencies at one (azDeg, elDeg)
   * while sharing vertex projections, polygon clipping, and polarimetric weights across bands!
   */
  function evaluateRCSMultiFreq(model, freqGHzArray, azDeg, elDeg, includeMB = true, pol = "Mean", opts = null) {
    const polTag = normalizePolarization(pol);
    const includePTD = Boolean(opts && opts.includeEdgeDiffraction);
    const az = azDeg * DEG2RAD;
    const el = elDeg * DEG2RAD;
    const ux = -Math.cos(el) * Math.cos(az);
    const uy = Math.cos(el) * Math.sin(az);
    const uz = Math.sin(el);

    const g = model.poGeom;
    const { nLit, nLV } = prepareLitPOGeometry(g, ux, uy, uz);
    const patches = (includeMB && model.concaveData && model.concaveData.numPairs > 0)
      ? collectDoubleBouncePatches(model.concaveData, ux, uy, uz, az)
      : [];

    const noiseFloor = 1e-5;
    return freqGHzArray.map((freqGHz) => {
      const lam = (C0 * 1e-9) / freqGHz;
      const k0 = (2.0 * Math.PI) / lam;

      evaluateLitPOAtK(g, nLit, nLV, k0);
      const rePO = _poSumOut[0];
      const imPO = _poSumOut[1];
      const ptd = (includePTD && model.rimEdges && model.rimEdges.count > 0)
        ? evaluateRimPTD(model.rimEdges, ux, uy, uz, k0, az, null)
        : _ZERO_PTD;
      const mb = patches.length > 0 ? evaluatePatchesAtK(patches, k0, false) : _ZERO_MB;
      return assemblePolarimetricPoint(rePO, imPO, mb, polTag, noiseFloor, ptd);
    });
  }

  /**
   * Evaluate coherent monostatic RCS (PO + optional PTD rim fringe + optional GO-PO double-bounce).
   */
  function evaluateRCSDirection(model, freqGHz, azDeg, elDeg, includeMB = true, pol = "Mean", opts = null) {
    const polTag = normalizePolarization(pol);
    const includePTD = Boolean(opts && opts.includeEdgeDiffraction);
    const lam = (C0 * 1e-9) / freqGHz;
    const k0 = (2.0 * Math.PI) / lam;
    const az = azDeg * DEG2RAD;
    const el = elDeg * DEG2RAD;
    const ux = -Math.cos(el) * Math.cos(az);
    const uy = Math.cos(el) * Math.sin(az);
    const uz = Math.sin(el);

    const g = model.poGeom;
    const { nLit, nLV } = prepareLitPOGeometry(g, ux, uy, uz);
    evaluateLitPOAtK(g, nLit, nLV, k0);
    const rePO = _poSumOut[0];
    const imPO = _poSumOut[1];

    const ptd = (includePTD && model.rimEdges && model.rimEdges.count > 0)
      ? evaluateRimPTD(model.rimEdges, ux, uy, uz, k0, az, null)
      : _ZERO_PTD;

    const mb = (includeMB && model.concaveData && model.concaveData.numPairs > 0)
      ? doubleBounceGOPO(model.concaveData, ux, uy, uz, k0, false, az)
      : _ZERO_MB;

    return assemblePolarimetricPoint(rePO, imPO, mb, polTag, 1e-5, ptd);
  }

  // --------------------------------------------------------------------------
  // 7. Wideband Range-Angle-Intensity (RTI / Sinogram) & HRRP Solver
  // --------------------------------------------------------------------------

  function nsinc(x) {
    return Math.abs(x) > 1e-12 ? Math.sin(Math.PI * x) / (Math.PI * x) : 1.0;
  }

  function computeRangeAngleIntensity(model, freqGHz, bwMHz, elDeg, includeMB = true, winType = "hann", opts = {}) {
    if (!(bwMHz > 0)) {
      throw new RangeError("Bandwidth must be strictly positive.");
    }
    const polTag = normalizePolarization(opts.polarization || "Mean");
    const hasSingleBouncePol =
      polTag === "Mean" || polTag === "HH" || polTag === "VV" || polTag === "PP";
    const maxAllowedBW = 1000.0 * freqGHz;
    if (opts.strictBandwidth && bwMHz > maxAllowedBW) {
      throw new RangeError(
        `Bandwidth ${bwMHz} MHz exceeds center frequency ${freqGHz} GHz (band would include f <= 0).`
      );
    }
    const effectiveBwMHz = Math.min(bwMHz, maxAllowedBW);
    const B_Hz = effectiveBwMHz * 1e6;
    const deltaR_m = C0 / (2.0 * B_Hz);

    const wName = String(winType || "hann").toLowerCase();
    let a0 = 0.50, a1 = 0.25, w3 = 1.44;
    if (wName === "hamming") {
      a0 = 0.54; a1 = 0.23; w3 = 1.30;
    } else if (wName === "rect") {
      a0 = 1.00; a1 = 0.00; w3 = 0.886;
    } else if (wName !== "hann") {
      throw new RangeError(`Unsupported window "${winType}". Expected rect, hann, or hamming.`);
    }
    const effDeltaR_m = w3 * deltaR_m;

    const lam = (C0 * 1e-9) / freqGHz;
    const k0 = (2.0 * Math.PI) / lam;
    const pref = k0 / Math.sqrt(Math.PI);

    // True axial extent from mesh vertices
    const xMin = model.xMin;
    const xMax = model.xMax;
    const totalLength = xMax - xMin;
    const xRef = 0.5 * (xMin + xMax);
    const maxSpan = Math.max(Math.abs(model.yMin), Math.abs(model.yMax));
    const rLim = opts.rLim || (Math.max(0.6 * totalLength, maxSpan * 1.1) + 4.0 * effDeltaR_m + 0.3);

    // Fine internal range grid (<= deltaR / 8) for faithful PSF mainlobe & sidelobes
    const maxBins = opts.maxRangeBins || 2401;
    const minBins = opts.minRangeBins || 241;
    let nR = Math.max(minBins, Math.ceil((2.0 * rLim) / (deltaR_m / 8.0)) + 1);
    if (nR > maxBins) nR = maxBins;
    if (nR % 2 === 0) nR += 1;

    const rangeVec = new Float64Array(nR);
    const dBin = (2.0 * rLim) / (nR - 1);
    for (let ir = 0; ir < nR; ir++) rangeVec[ir] = -rLim + ir * dBin;

    const halfK = Math.min(Math.ceil((6.0 * deltaR_m) / dBin), Math.floor(nR / 2));
    const kern = new Float64Array(2 * halfK + 1);
    for (let dk = -halfK; dk <= halfK; dk++) {
      const x = (dk * dBin) / deltaR_m;
      kern[dk + halfK] = (a0 * nsinc(x) + a1 * (nsinc(x - 1.0) + nsinc(x + 1.0))) / a0;
    }

    let azVec;
    if (opts.azVec) {
      azVec = Array.from(opts.azVec);
    } else {
      azVec = [];
      for (let a = 0; a <= 171; a += 3) azVec.push(a);
      for (let a = 173; a <= 180; a += 1) azVec.push(a);
    }
    const includePTD = Boolean(opts.includeEdgeDiffraction);
    const nAz = azVec.length;
    const elRad = elDeg * DEG2RAD;
    const g = model.poGeom;
    const nF = g.nF;

    const RTI_total_dBsm = new Array(nAz);
    const RTI_single_dBsm = new Array(nAz);
    const RTI_double_dBsm = new Array(nAz);
    const R_nose = new Float64Array(nAz);
    const R_base = new Float64Array(nAz);
    const R_aft = new Float64Array(nAz);

    const accSingleRe = new Float64Array(nR);
    const accSingleIm = new Float64Array(nR);
    const accDoubleRe = new Float64Array(nR);
    const accDoubleIm = new Float64Array(nR);
    const noiseFloor = 1e-5; // -50 dBsm floor

    for (let ia = 0; ia < nAz; ia++) {
      const azRad = azVec[ia] * DEG2RAD;
      const ux = -Math.cos(elRad) * Math.cos(azRad);
      const uy = Math.cos(elRad) * Math.sin(azRad);
      const uz = Math.sin(elRad);

      R_nose[ia] = -((xMin - xRef) * ux);
      R_base[ia] = -((model.cfg.length - xRef) * ux);
      R_aft[ia] = -((xMax - xRef) * ux);

      accSingleRe.fill(0.0);
      accSingleIm.fill(0.0);
      accDoubleRe.fill(0.0);
      accDoubleIm.fill(0.0);

      // 1. Single-bounce PO scattering centers (zero for HV, VH, OP, and LL channels)
      if (hasSingleBouncePol) {
        const { nLit, nLV } = prepareLitPOGeometry(g, ux, uy, uz);
        evaluateLitPOAtK(g, nLit, nLV, k0, accSingleRe, accSingleIm, ux, xRef, rLim, 1.0 / dBin, nR);
      }

      // 1b. Optional polarimetric Ufimtsev/Mitzner PTD base-rim edge diffraction centers
      if (includePTD && model.rimEdges && model.rimEdges.count > 0) {
        evaluateRimPTD(model.rimEdges, ux, uy, uz, k0, azRad, {
          polTag, accRe: accSingleRe, accIm: accSingleIm, xRef, rLim, invDBin: 1.0 / dBin, nR
        });
      }

      // 2. GO-PO concave double-bounce scattering centers with true cavity path delay
      if (includeMB && model.concaveData && model.concaveData.numPairs > 0) {
        const mb = doubleBounceGOPO(model.concaveData, ux, uy, uz, k0, true, azRad);
        const rows = mb.contrib;
        for (let r = 0; r < rows.length; r++) {
          const row = rows[r];
          const Lhalf = row[2];
          let aRe = 0.0, aIm = 0.0;
          switch (polTag) {
            case "Mean":
              aRe = row[0]; aIm = row[1];
              break;
            case "HH":
              aRe = row[3]; aIm = row[4];
              break;
            case "VV":
              aRe = row[5]; aIm = row[6];
              break;
            case "HV":
            case "VH":
              aRe = row[7]; aIm = row[8];
              break;
            case "PP":
              aRe = 0.5 * (row[3] + row[5]);
              aIm = 0.5 * (row[4] + row[6]);
              break;
            case "OP":
              aRe = 0.5 * (row[3] - row[5]) + row[8];
              aIm = 0.5 * (row[4] - row[6]) - row[7];
              break;
            case "LL":
              aRe = 0.5 * (row[3] - row[5]) - row[8];
              aIm = 0.5 * (row[4] - row[6]) + row[7];
              break;
          }
          const RD = Lhalf + ux * xRef;
          const fBin = (RD + rLim) / dBin;
          const i0 = Math.floor(fBin);
          if (i0 >= 0 && i0 < nR - 1) {
            const t = fBin - i0;
            accDoubleRe[i0] += aRe * (1.0 - t);
            accDoubleIm[i0] += aIm * (1.0 - t);
            accDoubleRe[i0 + 1] += aRe * t;
            accDoubleIm[i0 + 1] += aIm * t;
          }
        }
      }

      // 3. Convolve deposited scatterers with exact Harris spectral window PSF
      const rowTot = new Float64Array(nR);
      const rowSng = new Float64Array(nR);
      const rowDbl = new Float64Array(nR);

      // Sparse-aware convolution: gather non-zero source indices once per azimuth
      const nzSingle = [];
      const nzDouble = [];
      for (let s = 0; s < nR; s++) {
        if (accSingleRe[s] !== 0.0 || accSingleIm[s] !== 0.0) nzSingle.push(s);
        if (accDoubleRe[s] !== 0.0 || accDoubleIm[s] !== 0.0) nzDouble.push(s);
      }

      const convSRe = new Float64Array(nR);
      const convSIm = new Float64Array(nR);
      const convDRe = new Float64Array(nR);
      const convDIm = new Float64Array(nR);

      for (let idx = 0; idx < nzSingle.length; idx++) {
        const s = nzSingle[idx];
        const sr = accSingleRe[s], si = accSingleIm[s];
        const irMin = Math.max(0, s - halfK);
        const irMax = Math.min(nR - 1, s + halfK);
        for (let ir = irMin; ir <= irMax; ir++) {
          const kv = kern[ir - s + halfK];
          convSRe[ir] += sr * kv;
          convSIm[ir] += si * kv;
        }
      }

      for (let idx = 0; idx < nzDouble.length; idx++) {
        const s = nzDouble[idx];
        const dr = accDoubleRe[s], di = accDoubleIm[s];
        const irMin = Math.max(0, s - halfK);
        const irMax = Math.min(nR - 1, s + halfK);
        for (let ir = irMin; ir <= irMax; ir++) {
          const kv = kern[ir - s + halfK];
          convDRe[ir] += dr * kv;
          convDIm[ir] += di * kv;
        }
      }

      for (let ir = 0; ir < nR; ir++) {
        const sRe = convSRe[ir], sIm = convSIm[ir];
        const dRe = convDRe[ir], dIm = convDIm[ir];
        const sigSng = sRe * sRe + sIm * sIm + noiseFloor;
        const sigDbl = dRe * dRe + dIm * dIm + noiseFloor;
        const sigTot = (sRe + dRe) * (sRe + dRe) + (sIm + dIm) * (sIm + dIm) + noiseFloor;
        rowSng[ir] = 10.0 * Math.log10(sigSng);
        rowDbl[ir] = 10.0 * Math.log10(sigDbl);
        rowTot[ir] = 10.0 * Math.log10(sigTot);
      }

      RTI_single_dBsm[ia] = rowSng;
      RTI_double_dBsm[ia] = rowDbl;
      RTI_total_dBsm[ia] = rowTot;
    }

    return {
      freqGHz,
      bwMHz: effectiveBwMHz,
      requestedBwMHz: bwMHz,
      clampedBW: effectiveBwMHz < bwMHz,
      fractionalBW: effectiveBwMHz / (1000.0 * freqGHz),
      deltaR_m,
      effDeltaR_m,
      polarization: polTag,
      winType: wName,
      elDeg,
      azVec,
      rangeVec: Array.from(rangeVec),
      rLim,
      xRef,
      RTI_total_dBsm,
      RTI_single_dBsm,
      RTI_double_dBsm,
      R_nose: Array.from(R_nose),
      R_base: Array.from(R_base),
      R_aft: Array.from(R_aft)
    };
  }

  // --------------------------------------------------------------------------
  // 8. Generalized Winding Number Solid-Angle Verifier (van Oosterom & Strackee 1983)
  // --------------------------------------------------------------------------

  function generalizedWindingNumber(vertices, faces, queryPt) {
    const qx = queryPt[0], qy = queryPt[1], qz = queryPt[2];
    let omegaSum = 0.0;
    for (let f = 0; f < faces.length; f++) {
      const [i1, i2, i3] = faces[f];
      const v1 = vertices[i1], v2 = vertices[i2], v3 = vertices[i3];
      const ax = v1[0] - qx, ay = v1[1] - qy, az = v1[2] - qz;
      const bx = v2[0] - qx, by = v2[1] - qy, bz = v2[2] - qz;
      const cx = v3[0] - qx, cy = v3[1] - qy, cz = v3[2] - qz;

      const la = Math.hypot(ax, ay, az);
      const lb = Math.hypot(bx, by, bz);
      const lc = Math.hypot(cx, cy, cz);
      if (la < 1e-12 || lb < 1e-12 || lc < 1e-12) continue;

      const det = ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx);
      const den =
        la * lb * lc +
        (ax * bx + ay * by + az * bz) * lc +
        (bx * cx + by * cy + bz * cz) * la +
        (cx * ax + cy * ay + cz * az) * lb;
      omegaSum += 2.0 * Math.atan2(det, den);
    }
    return omegaSum / (4.0 * Math.PI);
  }

  /**
   * Build a minimal model object from arbitrary triangle mesh (and optional structured lower surface)
   * for canonical electromagnetic & HRRP analytical benchmark validation.
   */
  function buildCanonicalTargetModel(vertices, faces, structuredLow = null, lengthVal = 1.0) {
    const nF = faces.length;
    const nV = vertices.length;
    const poGeom = {
      nF,
      nV,
      facePatch: null,
      aftAperture: null,
      vx: new Float64Array(nV), vy: new Float64Array(nV), vz: new Float64Array(nV),
      f1: new Int32Array(nF), f2: new Int32Array(nF), f3: new Int32Array(nF),
      uDotV: new Float64Array(nV), cosV: new Float64Array(nV), sinV: new Float64Array(nV),
      vStamp: new Int32Array(nV), stampCounter: 0,
      litVerts: new Int32Array(nV), litFaces: new Int32Array(nF),
      cosLit: new Float64Array(nF), ut2Lit: new Float64Array(nF),
      um1Lit: new Float64Array(nF), um2Lit: new Float64Array(nF), um3Lit: new Float64Array(nF),
      nx: new Float64Array(nF), ny: new Float64Array(nF), nz: new Float64Array(nF),
      rcx: new Float64Array(nF), rcy: new Float64Array(nF), rcz: new Float64Array(nF),
      area: new Float64Array(nF), Lmax2: new Float64Array(nF),
      m1x: new Float64Array(nF), m1y: new Float64Array(nF), m1z: new Float64Array(nF),
      m2x: new Float64Array(nF), m2y: new Float64Array(nF), m2z: new Float64Array(nF),
      m3x: new Float64Array(nF), m3y: new Float64Array(nF), m3z: new Float64Array(nF)
    };

    let xMin = Infinity, xMax = -Infinity;
    let yMin = Infinity, yMax = -Infinity;
    let zMin = Infinity, zMax = -Infinity;
    for (let v = 0; v < nV; v++) {
      const pt = vertices[v];
      poGeom.vx[v] = pt[0]; poGeom.vy[v] = pt[1]; poGeom.vz[v] = pt[2];
      if (pt[0] < xMin) xMin = pt[0]; if (pt[0] > xMax) xMax = pt[0];
      if (pt[1] < yMin) yMin = pt[1]; if (pt[1] > yMax) yMax = pt[1];
      if (pt[2] < zMin) zMin = pt[2]; if (pt[2] > zMax) zMax = pt[2];
    }

    for (let f = 0; f < nF; f++) {
      const [i1, i2, i3] = faces[f];
      const v1 = vertices[i1], v2 = vertices[i2], v3 = vertices[i3];
      const e1x = v2[0] - v1[0], e1y = v2[1] - v1[1], e1z = v2[2] - v1[2];
      const e2x = v3[0] - v2[0], e2y = v3[1] - v2[1], e2z = v3[2] - v2[2];
      const e3x = v1[0] - v3[0], e3y = v1[1] - v3[1], e3z = v1[2] - v3[2];
      const wx = -e3x, wy = -e3y, wz = -e3z;
      const cx = 0.5 * (e1y * wz - e1z * wy);
      const cy = 0.5 * (e1z * wx - e1x * wz);
      const cz = 0.5 * (e1x * wy - e1y * wx);
      const trueArea = Math.hypot(cx, cy, cz);
      const area = trueArea + 1e-15;
      const nx = cx / area, ny = cy / area, nz = cz / area;
      const l1Sq = e1x * e1x + e1y * e1y + e1z * e1z;
      const l2Sq = e2x * e2x + e2y * e2y + e2z * e2z;
      const l3Sq = e3x * e3x + e3y * e3y + e3z * e3z;

      poGeom.f1[f] = i1; poGeom.f2[f] = i2; poGeom.f3[f] = i3;
      poGeom.nx[f] = nx; poGeom.ny[f] = ny; poGeom.nz[f] = nz;
      poGeom.rcx[f] = (v1[0] + v2[0] + v3[0]) / 3.0;
      poGeom.rcy[f] = (v1[1] + v2[1] + v3[1]) / 3.0;
      poGeom.rcz[f] = (v1[2] + v2[2] + v3[2]) / 3.0;
      poGeom.area[f] = trueArea;
      poGeom.Lmax2[f] = Math.max(l1Sq, l2Sq, l3Sq);

      poGeom.m1x[f] = e1y * nz - e1z * ny; poGeom.m1y[f] = e1z * nx - e1x * nz; poGeom.m1z[f] = e1x * ny - e1y * nx;
      poGeom.m2x[f] = e2y * nz - e2z * ny; poGeom.m2y[f] = e2z * nx - e2x * nz; poGeom.m2z[f] = e2x * ny - e2y * nx;
      poGeom.m3x[f] = e3y * nz - e3z * ny; poGeom.m3y[f] = e3z * nx - e3x * nz; poGeom.m3z[f] = e3x * ny - e3y * nx;
    }

    const concaveData = structuredLow
      ? buildConcaveSuperCells(structuredLow.X, structuredLow.Y, structuredLow.Z, structuredLow.X.length, structuredLow.X[0].length)
      : { numPairs: 0 };

    return {
      cfg: { length: lengthVal },
      vertices,
      faces,
      poGeom,
      concaveData,
      xMin, xMax, yMin, yMax, zMin, zMax
    };
  }

  return {
    VERSION: "20261113-base",
    C0,
    solveObliqueShock,
    solveTaylorMaccoll,
    solveInverseTaylorMaccoll,
    traceStreamline,
    normalizeConfig,
    validateConfig,
    buildWaveriderModel,
    buildCanonicalTargetModel,
    buildConcaveSuperCells,
    buildAxisymmetricConcaveSuperCells,
    evaluateAeroState,
    computeAeroPolar,
    facetPhaseIntegralTriangles,
    clipConvexPolygon,
    normalizePolarization,
    computePolarizationBasis,
    doubleBounceGOPO,
    evaluateRCSDirection,
    evaluateRCSMultiFreq,
    computeRangeAngleIntensity,
    generalizedWindingNumber
  };
});
