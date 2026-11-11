/* ============================================================================
 * Hypersonic Parametric Watertight Waverider & UHF-to-X-Band RCS Web Studio
 * Real-time client-side Taylor-Maccoll, Conformal Frustum, Watertight Mesh
 * Stitcher, Topological Verifier, and Coherent PO + Multi-Bounce RCS Engine.
 * ============================================================================ */

(function () {
  "use strict";

  // 1. Modern Web Guidance Fallbacks: Range Slider CSS Variables & Dialog Light-Dismiss
  const rangeInputs = document.querySelectorAll('input[type="range"]');
  function syncSliderVars(slider) {
    slider.style.setProperty("--attr-min", slider.min || 0);
    slider.style.setProperty("--attr-max", slider.max || 100);
    slider.style.setProperty("--control-value", slider.value);
    slider.setAttribute("aria-valuetext", slider.value);
  }
  rangeInputs.forEach((s) => {
    syncSliderVars(s);
    s.addEventListener("input", () => syncSliderVars(s));
  });

  // Fallback for browsers without native <dialog closedby="any"> support
  if (typeof HTMLDialogElement !== "undefined" && !("closedBy" in HTMLDialogElement.prototype)) {
    const outside = (e, dialog) => {
      const r = dialog.getBoundingClientRect();
      return e.clientX < r.left || e.clientX > r.right || e.clientY < r.top || e.clientY > r.bottom;
    };
    for (const dialog of document.querySelectorAll('dialog[closedby="any"]')) {
      let pressedOutside = false;
      dialog.addEventListener("pointerdown", (e) => {
        pressedOutside = outside(e, dialog);
      });
      dialog.addEventListener("pointerup", (e) => {
        if (pressedOutside && outside(e, dialog)) dialog.close();
      });
    }
  }

  // 2. Application State & Presets
  const PRESETS = {
    osculating_cone: {
      method: "osculating_cone", mach: 7.0, beta: 13.5, length: 5.5, span: 1.75,
      power: 0.62, anhedral: 8.0, canopy: 0.045, frustum: false, flength: 1.65,
      bradius: 0.72, blend: "cubic", rocket: false, rcone: 4.5
    },
    conical: {
      method: "conical", mach: 6.0, beta: 14.5, length: 5.0, span: 1.55,
      power: 0.68, anhedral: 6.0, canopy: 0.000, frustum: false, flength: 1.50,
      bradius: 0.65, blend: "cubic", rocket: false, rcone: 4.0
    },
    caret: {
      method: "caret", mach: 6.5, beta: 15.0, length: 5.0, span: 1.50,
      power: 0.85, anhedral: 11.0, canopy: 0.035, frustum: false, flength: 1.50,
      bradius: 0.65, blend: "cubic", rocket: false, rcone: 4.5
    },
    cst: {
      method: "cst", mach: 8.0, beta: 12.5, length: 6.0, span: 1.80,
      power: 0.60, anhedral: 7.0, canopy: 0.042, frustum: false, flength: 1.70,
      bradius: 0.75, blend: "cubic", rocket: false, rcone: 4.5
    },
    booster_mounted: {
      method: "osculating_cone", mach: 7.0, beta: 13.5, length: 5.5, span: 1.75,
      power: 0.62, anhedral: 8.5, canopy: 0.045, frustum: true, flength: 1.65,
      bradius: 0.72, blend: "cubic", rocket: true, rcone: 4.5
    }
  };

  let currentModel = null;
  let camera = { yaw: -2.32, pitch: 0.38, zoom: 1.0, dragging: false, lastX: 0, lastY: 0 };

  // 3. Exact Oblique Shock & Taylor-Maccoll Conical Flow Solvers in JS
  function solveObliqueShock(M1, betaDeg, gamma = 1.4) {
    const muRad = Math.asin(Math.min(1.0 / M1, 0.99));
    const betaRad = Math.max((betaDeg * Math.PI) / 180.0, muRad + 0.02);
    const M1n = M1 * Math.sin(betaRad);
    const cotBeta = 1.0 / Math.tan(betaRad);
    const num = 2.0 * cotBeta * (M1n * M1n - 1.0);
    const den = M1 * M1 * (gamma + Math.cos(2.0 * betaRad)) + 2.0;
    const deltaRad = Math.atan(num / den);
    const M2n_sq = (1.0 + 0.5 * (gamma - 1.0) * M1n * M1n) / (gamma * M1n * M1n - 0.5 * (gamma - 1.0));
    const M2 = Math.sqrt(Math.max(M2n_sq, 1e-6)) / Math.sin(Math.max(betaRad - deltaRad, 1e-4));
    const p2_p1 = 1.0 + ((2.0 * gamma) / (gamma + 1.0)) * (M1n * M1n - 1.0);
    const V2_bar = Math.pow(1.0 + 2.0 / ((gamma - 1.0) * M2 * M2), -0.5);
    const Cp2 = (p2_p1 - 1.0) / (0.5 * gamma * M1 * M1);
    return { M1, gamma, betaDeg, betaRad, deltaRad, deltaDeg: (deltaRad * 180) / Math.PI, M2, p2_p1, V2_bar, Cp2 };
  }

  function solveTaylorMaccoll(M1, betaDeg, gamma = 1.4) {
    const shock = solveObliqueShock(M1, betaDeg, gamma);
    const betaRad = shock.betaRad;
    let ur = shock.V2_bar * Math.cos(betaRad - shock.deltaRad);
    let ut = -shock.V2_bar * Math.sin(betaRad - shock.deltaRad);
    let th = betaRad;
    const dth = -0.0005;

    const thetaArr = [th];
    const slopeArr = [];
    const cpArr = [];

    function pushState(theta, u_r, u_t) {
      const ux = u_r * Math.cos(theta) - u_t * Math.sin(theta);
      const ucyl = u_r * Math.sin(theta) + u_t * Math.cos(theta);
      slopeArr.push(ucyl / Math.max(ux, 1e-6));
      const Vsq = u_r * u_r + u_t * u_t;
      const MlocSq = (2.0 / (gamma - 1.0)) * (Vsq / Math.max(1.0 - Vsq, 1e-8));
      const T0_T2 = 1.0 + 0.5 * (gamma - 1.0) * shock.M2 * shock.M2;
      const T0_Tloc = 1.0 + 0.5 * (gamma - 1.0) * MlocSq;
      const p_p1 = shock.p2_p1 * Math.pow(T0_T2 / T0_Tloc, gamma / (gamma - 1.0));
      cpArr.push((p_p1 - 1.0) / (0.5 * gamma * M1 * M1));
    }
    pushState(th, ur, ut);

    function tmDeriv(theta, u_r, u_t) {
      const a_sq = 0.5 * (gamma - 1.0) * (1.0 - u_r * u_r - u_t * u_t);
      const num = a_sq * (2.0 * u_r + u_t / Math.tan(theta)) - u_r * u_t * u_t;
      let den = u_t * u_t - a_sq;
      if (Math.abs(den) < 1e-9) den = -1e-9;
      return [u_t, num / den];
    }

    for (let step = 0; step < 2500; step++) {
      if (th + dth <= 0.01) break;
      const [k1r, k1t] = tmDeriv(th, ur, ut);
      const [k2r, k2t] = tmDeriv(th + 0.5 * dth, ur + 0.5 * dth * k1r, ut + 0.5 * dth * k1t);
      const [k3r, k3t] = tmDeriv(th + 0.5 * dth, ur + 0.5 * dth * k2r, ut + 0.5 * dth * k2t);
      const [k4r, k4t] = tmDeriv(th + dth, ur + dth * k3r, ut + dth * k3t);
      const urNext = ur + (dth / 6.0) * (k1r + 2 * k2r + 2 * k3r + k4r);
      const utNext = ut + (dth / 6.0) * (k1t + 2 * k2t + 2 * k3t + k4t);
      if (utNext >= 0.0) {
        const frac = -ut / Math.max(utNext - ut, 1e-9);
        th = th + frac * dth;
        ur = ur + frac * (urNext - ur);
        ut = 0.0;
        thetaArr.push(th);
        pushState(th, ur, ut);
        break;
      }
      th += dth;
      ur = urNext;
      ut = utNext;
      thetaArr.push(th);
      pushState(th, ur, ut);
    }

    thetaArr.reverse();
    slopeArr.reverse();
    cpArr.reverse();

    function interpTM(thetaQuery) {
      const tClamped = Math.max(thetaArr[0], Math.min(thetaQuery, thetaArr[thetaArr.length - 1]));
      let lo = 0, hi = thetaArr.length - 1;
      while (hi - lo > 1) {
        const mid = (lo + hi) >> 1;
        if (thetaArr[mid] <= tClamped) lo = mid;
        else hi = mid;
      }
      const w = (tClamped - thetaArr[lo]) / Math.max(thetaArr[hi] - thetaArr[lo], 1e-9);
      return {
        slope: slopeArr[lo] * (1 - w) + slopeArr[hi] * w,
        cp: cpArr[lo] * (1 - w) + cpArr[hi] * w
      };
    }

    return {
      shock, betaRad, thetaConeRad: thetaArr[0],
      thetaConeDeg: (thetaArr[0] * 180) / Math.PI,
      interpTM
    };
  }

  // 4. Generate 3D Waverider + Conformal Frustum + Conical Booster + Watertight Mesh
  function buildWaveriderModel(cfg) {
    const L = cfg.length;
    const bHalf = cfg.span;
    const tm = solveTaylorMaccoll(cfg.mach, cfg.beta, 1.4);
    const tanBeta = Math.tan(tm.betaRad);
    const tanDelta = Math.tan(tm.shock.deltaRad);
    const anhedralRad = (cfg.anhedral * Math.PI) / 180.0;

    const NxWave = 24;
    const NyHalf = 19;
    const Ny = 2 * NyHalf - 1;
    const midCol = NyHalf - 1;

    const eta = new Array(Ny);
    const yLE = new Array(Ny);
    const xLE = new Array(Ny);
    const zLE = new Array(Ny);
    const invP = 1.0 / Math.max(cfg.power, 0.25);

    const x0_ref = Math.max(0.58 * L, (1.28 * bHalf) / tanBeta - L);
    const zv = x0_ref * tanBeta;

    for (let k = 0; k < Ny; k++) {
      eta[k] = -1.0 + (2.0 * k) / (Ny - 1);
      yLE[k] = bHalf * eta[k];
      xLE[k] = L * Math.pow(Math.abs(eta[k]), invP);
    }
    xLE[0] = L;
    xLE[Ny - 1] = L;
    xLE[midCol] = 0.0;
    yLE[midCol] = 0.0;

    for (let k = 0; k < Ny; k++) {
      if (cfg.method === "conical") {
        const R_shk = (xLE[k] + x0_ref) * tanBeta;
        const zDepth = Math.sqrt(Math.max(R_shk * R_shk - yLE[k] * yLE[k], 1e-8));
        zLE[k] = zv - zDepth - Math.abs(yLE[k]) * Math.tan(0.35 * anhedralRad);
      } else if (cfg.method === "osculating_cone") {
        const fFlat = 0.28;
        const outerEta = Math.max(Math.abs(eta[k]) - fFlat, 0.0) / (1.0 - fFlat);
        zLE[k] = -0.65 * xLE[k] * tanBeta - bHalf * Math.tan(anhedralRad) * Math.pow(outerEta, 2.1);
      } else if (cfg.method === "caret") {
        zLE[k] = -xLE[k] * tanBeta - Math.abs(yLE[k]) * Math.tan(anhedralRad);
      } else {
        zLE[k] = -0.65 * xLE[k] * tanBeta - Math.abs(yLE[k]) * Math.tan(anhedralRad) * Math.pow(Math.abs(eta[k]), 0.85);
      }
    }
    zLE[midCol] = 0.0;

    // Allocate forebody grids [NxWave][Ny]
    const X_up = [], Y_up = [], Z_up = [];
    const X_low = [], Y_low = [], Z_low = [];
    const X_shk = [], Y_shk = [], Z_shk = [];
    const Cp_low = [];

    for (let i = 0; i < NxWave; i++) {
      X_up.push(new Array(Ny)); Y_up.push(new Array(Ny)); Z_up.push(new Array(Ny));
      X_low.push(new Array(Ny)); Y_low.push(new Array(Ny)); Z_low.push(new Array(Ny));
      X_shk.push(new Array(Ny)); Y_shk.push(new Array(Ny)); Z_shk.push(new Array(Ny));
      Cp_low.push(new Array(Ny));
    }

    for (let k = 0; k < Ny; k++) {
      const dxMax = Math.max(L - xLE[k], 0.0);
      let rCurr = (xLE[k] + x0_ref) * tanBeta;
      const r0 = rCurr;
      const phiLE = Math.atan2(yLE[k], Math.max(zv - zLE[k], 1e-5));
      const cosPhi = Math.cos(phiLE);

      for (let i = 0; i < NxWave; i++) {
        const s = i / (NxWave - 1);
        const xVal = xLE[k] + s * dxMax;
        const dx = xVal - xLE[k];

        let zLowVal = zLE[k];
        let cpVal = tm.shock.Cp2;

        if (cfg.method === "conical" || cfg.method === "osculating_cone") {
          if (i > 0 && dxMax > 1e-7) {
            const dxStep = dxMax / (NxWave - 1);
            const xLoc = xLE[k] + x0_ref + (i - 0.5) * dxStep;
            const thLoc = Math.atan2(rCurr, xLoc);
            const tmEval = tm.interpTM(thLoc);
            rCurr += dxStep * tmEval.slope;
            cpVal = tmEval.cp;
          } else {
            cpVal = tm.interpTM(tm.betaRad).cp;
          }
          const dr = rCurr - r0;
          zLowVal = zLE[k] - dr * cosPhi;
        } else if (cfg.method === "caret") {
          zLowVal = zLE[k] - dx * tanDelta;
          cpVal = tm.shock.Cp2;
        } else {
          // CST Boost-Glide
          const archBell = Math.pow(1.0 - eta[k] * eta[k], 2.0);
          const chordScale = Math.pow(Math.max(1.0 - xLE[k] / L, 0.0), 0.65);
          const dzArch = 0.025 * L * Math.pow(s, 0.9) * archBell * chordScale;
          zLowVal = Math.min(zLE[k] - dx * tanDelta + dzArch, zLE[k] - 0.25 * dx * tanDelta);
          cpVal = tm.shock.Cp2 * (1.0 + 0.25 * archBell);
        }

        // Upper Surface (Freestream + Dorsal Canopy)
        const spanBell = Math.pow(1.0 - Math.pow(Math.abs(eta[k]), 2.2), 2.0);
        const chordScale = Math.pow(Math.max(1.0 - xLE[k] / L, 0.0), 0.75);
        const longProf = Math.pow(s, 0.85) * (1.0 - 0.25 * s);
        const hCanopy = (dxMax > 1e-6) ? (cfg.canopy * L * chordScale * spanBell * longProf) : 0.0;
        const zUpVal = zLE[k] + hCanopy;

        X_low[i][k] = xVal; Y_low[i][k] = yLE[k]; Z_low[i][k] = zLowVal;
        X_up[i][k]  = xVal; Y_up[i][k]  = yLE[k]; Z_up[i][k]  = zUpVal;
        X_shk[i][k] = xVal; Y_shk[i][k] = yLE[k]; Z_shk[i][k] = zLE[k] - dx * tanBeta;
        Cp_low[i][k] = cpVal;
      }
      // Enforce exact coincidence at LE and wingtips
      Z_up[0][k] = Z_low[0][k];
    }
    for (let i = 0; i < NxWave; i++) {
      Z_up[i][0] = Z_low[i][0];
      Z_up[i][Ny - 1] = Z_low[i][Ny - 1];
    }

    // 5. Append Conformal Transition Frustum & Conical Rocket Stage if enabled
    let NxFrustEnd = NxWave;
    let NxTotal = NxWave;
    let zCenterBooster = -0.35;

    const yTE_up = Y_up[NxWave - 1].slice();
    const zTE_up = Z_up[NxWave - 1].slice();
    const yTE_low = Y_low[NxWave - 1].slice();
    const zTE_low = Z_low[NxWave - 1].slice();

    // Area-weighted vertical centroid of waverider base at X = L
    let areaSum = 0.0, zMomSum = 0.0;
    for (let k = 0; k < Ny - 1; k++) {
      const dy = yTE_up[k + 1] - yTE_up[k];
      const dzAvg = 0.5 * ((zTE_up[k] - zTE_low[k]) + (zTE_up[k + 1] - zTE_low[k + 1]));
      const zMid = 0.25 * (zTE_up[k] + zTE_low[k] + zTE_up[k + 1] + zTE_low[k + 1]);
      areaSum += dy * dzAvg;
      zMomSum += dy * dzAvg * zMid;
    }
    if (areaSum > 1e-6) zCenterBooster = zMomSum / areaSum;

    if (cfg.frustum && cfg.flength > 0) {
      const Nf = 14;
      const Rb = cfg.bradius;
      for (let m = 1; m < Nf; m++) {
        const tau = m / (Nf - 1);
        let sigma = tau;
        if (cfg.blend === "cubic") sigma = 3 * tau * tau - 2 * tau * tau * tau;
        else if (cfg.blend === "quintic") sigma = 10 * Math.pow(tau, 3) - 15 * Math.pow(tau, 4) + 6 * Math.pow(tau, 5);
        else if (cfg.blend === "conic") {
          const pE = 1.45;
          sigma = Math.pow(tau, pE) / (Math.pow(tau, pE) + Math.pow(1 - tau, pE));
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
          rowY_u[k] = (1 - sigma) * yTE_up[k] + sigma * yRingU;
          rowZ_u[k] = (1 - sigma) * zTE_up[k] + sigma * zRingU;

          rowX_l[k] = xVal;
          rowY_l[k] = (1 - sigma) * yTE_low[k] + sigma * yRingL;
          rowZ_l[k] = (1 - sigma) * zTE_low[k] + sigma * zRingL;
          rowCp[k] = (1 - 0.65 * sigma) * Cp_low[NxWave - 1][k];
        }
        rowY_l[0] = rowY_u[0]; rowZ_l[0] = rowZ_u[0];
        rowY_l[Ny - 1] = rowY_u[Ny - 1]; rowZ_l[Ny - 1] = rowZ_u[Ny - 1];

        X_up.push(rowX_u); Y_up.push(rowY_u); Z_up.push(rowZ_u);
        X_low.push(rowX_l); Y_low.push(rowY_l); Z_low.push(rowZ_l);
        Cp_low.push(rowCp);
      }
      NxFrustEnd = X_low.length;
      NxTotal = NxFrustEnd;

      if (cfg.rocket) {
        const Nr = 10;
        const Lr = 1.85;
        const tanCone = Math.tan((cfg.rcone * Math.PI) / 180.0);
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
            rowCp[k] = 2.0 * Math.pow(Math.sin((cfg.rcone * Math.PI) / 180.0), 2);
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

    // 6. Stitch 2-Manifold Watertight Triangulation
    const Nb = 6;
    const vertices = [];
    function addVertex(x, y, z) {
      vertices.push([x, y, z]);
      return vertices.length - 1;
    }

    const vid_up = Array.from({ length: NxTotal }, () => new Array(Ny));
    const vid_low = Array.from({ length: NxTotal }, () => new Array(Ny));
    const vid_base = Array.from({ length: Nb }, () => new Array(Ny));

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
        const xb = X_up[NxTotal - 1][k] * (1 - tau) + X_low[NxTotal - 1][k] * tau;
        const yb = Y_up[NxTotal - 1][k] * (1 - tau) + Y_low[NxTotal - 1][k] * tau;
        const zb = Z_up[NxTotal - 1][k] * (1 - tau) + Z_low[NxTotal - 1][k] * tau;
        vid_base[m][k] = addVertex(xb, yb, zb);
      }
      vid_base[Nb - 1][k] = vid_low[NxTotal - 1][k];
    }

    const faces = [];
    const facePatch = []; // 1=Upper, 2=Lower, 3=Base, 4=Frustum, 5=Rocket

    function triangulateSubGrid(grid, rStart, rEnd, reverse, patchCode) {
      for (let r = rStart; r < rEnd - 1; r++) {
        for (let c = 0; c < Ny - 1; c++) {
          const A = grid[r][c], B = grid[r + 1][c], C = grid[r + 1][c + 1], D = grid[r][c + 1];
          if (A !== B && B !== C && C !== A) {
            faces.push(reverse ? [A, C, B] : [A, B, C]);
            facePatch.push(patchCode);
          }
          if (A !== C && C !== D && D !== A) {
            faces.push(reverse ? [A, D, C] : [A, C, D]);
            facePatch.push(patchCode);
          }
        }
      }
    }

    triangulateSubGrid(vid_up, 0, NxWave, false, 1);
    triangulateSubGrid(vid_low, 0, NxWave, true, 2);
    if (NxFrustEnd > NxWave) {
      triangulateSubGrid(vid_up, NxWave - 1, NxFrustEnd, false, 4);
      triangulateSubGrid(vid_low, NxWave - 1, NxFrustEnd, true, 4);
    }
    if (NxTotal > NxFrustEnd) {
      triangulateSubGrid(vid_up, NxFrustEnd - 1, NxTotal, false, 5);
      triangulateSubGrid(vid_low, NxFrustEnd - 1, NxTotal, true, 5);
    }
    triangulateSubGrid(vid_base, 0, Nb, false, 3);

    // 7. Verify Manifold Watertightness, Volume, Planform Area & Inviscid Aero
    const edgeMap = new Map();
    let enclosedVol = 0.0;
    let planformArea = 0.0;

    const M = cfg.mach, g = 1.4;
    const CpMax = (2.0 / (g * M * M)) * (
      Math.pow(((g + 1) * (g + 1) * M * M) / (4 * g * M * M - 2 * (g - 1)), g / (g - 1)) *
      ((1 - g + 2 * g * M * M) / (g + 1)) - 1.0
    );
    const meanFlowCp = tm.shock.Cp2;

    let forceX = 0.0, forceZ = 0.0;
    const faceNormals = new Array(faces.length);
    const faceCenters = new Array(faces.length);
    const faceAreas = new Array(faces.length);

    for (let f = 0; f < faces.length; f++) {
      const [i1, i2, i3] = faces[f];
      const v1 = vertices[i1], v2 = vertices[i2], v3 = vertices[i3];

      // Half-edges for manifold check
      const pairs = [[i1, i2], [i2, i3], [i3, i1]];
      for (const [u, v] of pairs) {
        const key = u < v ? `${u}_${v}` : `${v}_${u}`;
        edgeMap.set(key, (edgeMap.get(key) || 0) + 1);
      }

      const ux = v2[0] - v1[0], uy = v2[1] - v1[1], uz = v2[2] - v1[2];
      const wx = v3[0] - v1[0], wy = v3[1] - v1[1], wz = v3[2] - v1[2];
      const cx = 0.5 * (uy * wz - uz * wy);
      const cy = 0.5 * (uz * wx - ux * wz);
      const cz = 0.5 * (ux * wy - uy * wx);
      const area = Math.sqrt(cx * cx + cy * cy + cz * cz) + 1e-15;

      faceAreas[f] = area;
      faceNormals[f] = [cx / area, cy / area, cz / area];
      faceCenters[f] = [(v1[0] + v2[0] + v3[0]) / 3, (v1[1] + v2[1] + v3[1]) / 3, (v1[2] + v2[2] + v3[2]) / 3];

      // Signed tetrahedral volume: (1/6) * v1 . (v2 x v3)
      const v2xv3_x = v2[1] * v3[2] - v2[2] * v3[1];
      const v2xv3_y = v2[2] * v3[0] - v2[0] * v3[2];
      const v2xv3_z = v2[0] * v3[1] - v2[1] * v3[0];
      enclosedVol += (v1[0] * v2xv3_x + v1[1] * v2xv3_y + v1[2] * v2xv3_z) / 6.0;

      if (facePatch[f] === 1) {
        planformArea += Math.abs(cz);
      }

      if (facePatch[f] === 1 || facePatch[f] === 2) {
        const nx = cx / area;
        const sinImpact = Math.max(-nx, 0.0);
        let cp = 0.0;
        if (facePatch[f] === 2) {
          cp = 0.65 * meanFlowCp + 0.35 * CpMax * sinImpact * sinImpact;
        } else {
          const sinLee = Math.max(nx, 0.0);
          cp = sinImpact > 0 ? CpMax * sinImpact * sinImpact : -Math.min(1.2 / (M * M), (2.0 / M) * sinLee);
        }
        forceX += -cp * cx;
        forceZ += -cp * cz;
      }
    }

    let boundaryEdges = 0;
    for (const count of edgeMap.values()) {
      if (count !== 2) boundaryEdges++;
    }
    const euler = vertices.length - edgeMap.size + faces.length;
    const etaV = Math.pow(Math.max(enclosedVol, 0), 2 / 3) / Math.max(planformArea, 1e-6);
    const CL = forceZ / Math.max(planformArea, 1e-6);
    const CD = forceX / Math.max(planformArea, 1e-6);
    const LD = CL / Math.max(CD, 1e-6);

    // 8. Build Mutually Visible Concave Lower-Surface Pairs for Double-Bounce RCS
    const lowerIndices = [];
    for (let f = 0; f < faces.length; f++) {
      if (facePatch[f] === 2 && (f % 3 === 0)) lowerIndices.push(f);
    }
    const concavePairs = [];
    for (let a = 0; a < lowerIndices.length; a++) {
      const fP = lowerIndices[a];
      const rP = faceCenters[fP], nP = faceNormals[fP];
      for (let b = a + 1; b < lowerIndices.length; b++) {
        const fQ = lowerIndices[b];
        const rQ = faceCenters[fQ], nQ = faceNormals[fQ];
        const dx = rQ[0] - rP[0], dy = rQ[1] - rP[1], dz = rQ[2] - rP[2];
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (dist < 0.08) continue;
        const cP = (nP[0] * dx + nP[1] * dy + nP[2] * dz) / dist;
        const cQ = -(nQ[0] * dx + nQ[1] * dy + nQ[2] * dz) / dist;
        if (cP > 0.02 && cQ > 0.02) {
          const effA = 3.0 * Math.sqrt(faceAreas[fP] * faceAreas[fQ] * cP * cQ);
          concavePairs.push({
            rP, rQ, nP, nQ, dist,
            dHat: [dx / dist, dy / dist, dz / dist],
            cP, cQ, effA
          });
        }
      }
    }
    concavePairs.sort((p1, p2) => p2.effA - p1.effA);
    const activePairs = concavePairs.slice(0, 180);

    return {
      cfg, X_up, Y_up, Z_up, X_low, Y_low, Z_low, X_shk, Y_shk, Z_shk,
      NxWave, NxFrustEnd, NxTotal, Ny, zCenterBooster,
      vertices, faces, facePatch, faceNormals, faceCenters, faceAreas,
      boundaryEdges, numEdges: edgeMap.size, euler,
      enclosedVol, planformArea, etaV, CL, CD, LD,
      concavePairs: activePairs, totalConcavePairs: concavePairs.length
    };
  }

  // 5. Fast Multi-Band Coherent PO + Concave Double-Bounce RCS Evaluator
  function evaluateRCSDirection(model, freqGHz, azDeg, elDeg, includeMB) {
    const c0 = 0.299792458; // m * GHz
    const lam = c0 / freqGHz;
    const k0 = (2.0 * Math.PI) / lam;

    const az = (azDeg * Math.PI) / 180.0;
    const el = (elDeg * Math.PI) / 180.0;
    const ux = -Math.cos(el) * Math.cos(az);
    const uy = Math.cos(el) * Math.sin(az);
    const uz = Math.sin(el);

    const { faces, faceNormals, faceCenters, faceAreas, concavePairs } = model;
    const prefactor = k0 / Math.sqrt(Math.PI);

    let reSingle = 0.0, imSingle = 0.0;
    for (let f = 0; f < faces.length; f += 2) {
      const n = faceNormals[f];
      const cosTh = n[0] * ux + n[1] * uy + n[2] * uz;
      if (cosTh <= 1e-4) continue;
      const rc = faceCenters[f];
      const rDotU = rc[0] * ux + rc[1] * uy + rc[2] * uz;
      const damp = 1.0 / (1.0 + 0.08 * k0 * Math.sqrt(faceAreas[f]) * (1.0 - cosTh));
      const amp = 2.0 * prefactor * faceAreas[f] * cosTh * damp;
      const phase = 2.0 * k0 * rDotU;
      reSingle += -amp * Math.sin(phase);
      imSingle += amp * Math.cos(phase);
    }

    let reMB = 0.0, imMB = 0.0;
    if (includeMB && concavePairs.length > 0) {
      const lobeSharp = Math.min(Math.max(0.35 * k0, 3.5), 28.0);
      const elecScale = Math.min(Math.pow(1.0 / (0.4 * lam), 0.55), 2.6);
      for (let p = 0; p < concavePairs.length; p++) {
        const cp = concavePairs[p];
        const cosP = cp.nP[0] * ux + cp.nP[1] * uy + cp.nP[2] * uz;
        const cosQ = cp.nQ[0] * ux + cp.nQ[1] * uy + cp.nQ[2] * uz;
        if (cosP <= 0.02 || cosQ <= 0.02) continue;

        const dDotU = cp.dHat[0] * ux + cp.dHat[1] * uy + cp.dHat[2] * uz;
        const alignPQ = -dDotU + 2.0 * cosP * cp.cP;
        const alignQP = dDotU + 2.0 * cosQ * cp.cQ;
        const alignMax = Math.max(alignPQ, alignQP);
        if (alignMax <= 0.62) continue;

        const wLobe = Math.exp(-lobeSharp * (1.0 - alignMax));
        const amp = prefactor * cp.effA * elecScale * Math.sqrt(cosP * cosQ) * wLobe;
        const phase = k0 * (
          (cp.rP[0] + cp.rQ[0]) * ux +
          (cp.rP[1] + cp.rQ[1]) * uy +
          (cp.rP[2] + cp.rQ[2]) * uz - cp.dist
        );
        reMB += -amp * Math.sin(phase);
        imMB += amp * Math.cos(phase);
      }
    }

    const noiseFloor = 1e-4;
    const sigSingle = reSingle * reSingle + imSingle * imSingle + noiseFloor;
    const sigTotal = (reSingle + reMB) * (reSingle + reMB) + (imSingle + imMB) * (imSingle + imMB) + noiseFloor;
    return {
      single_dBsm: 10 * Math.log10(sigSingle),
      total_dBsm: 10 * Math.log10(sigTotal),
      total_m2: sigTotal
    };
  }

  // 6. Canvas 1: Interactive 3D Viewport
  const canvas3D = document.getElementById("canvas-3d");
  function render3DViewport(model) {
    const rect = canvas3D.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    canvas3D.width = rect.width * dpr;
    canvas3D.height = rect.height * dpr;
    const ctx = canvas3D.getContext("2d");
    ctx.scale(dpr, dpr);
    const W = rect.width, H = rect.height;

    ctx.fillStyle = "#090e1a";
    ctx.fillRect(0, 0, W, H);

    const { vertices, faces, facePatch, faceNormals, faceCenters, cfg } = model;
    const xMid = 0.5 * cfg.length + (cfg.frustum ? 0.5 * cfg.flength : 0);
    const scale = Math.min(W, H) * 0.13 * camera.zoom * (5.5 / cfg.length);

    const cy = Math.cos(camera.yaw), sy = Math.sin(camera.yaw);
    const cp = Math.cos(camera.pitch), sp = Math.sin(camera.pitch);

    function projectPt(pt) {
      const x = pt[0] - xMid, y = pt[1], z = pt[2] + 0.35;
      const rx = x * cy - y * sy;
      const ry = x * sy + y * cy;
      const rz = z;
      const px = rx;
      const py = ry * sp - rz * cp;
      const depth = ry * cp + rz * sp;
      return [W * 0.5 + px * scale, H * 0.52 + py * scale, depth];
    }

    const projVerts = vertices.map(projectPt);
    const showWire = document.getElementById("chk-wireframe").checked;
    const showShock = document.getElementById("chk-shock").checked;
    const showRays = document.getElementById("chk-rays").checked;

    // Sort faces back-to-front (Painter's algorithm)
    const order = faces.map((_, idx) => idx);
    order.sort((a, b) => {
      const da = projVerts[faces[a][0]][2] + projVerts[faces[a][1]][2] + projVerts[faces[a][2]][2];
      const db = projVerts[faces[b][0]][2] + projVerts[faces[b][1]][2] + projVerts[faces[b][2]][2];
      return da - db;
    });

    const baseColors = {
      1: [52, 148, 226],  // Waverider Upper (Blue)
      2: [232, 112, 54],  // Waverider Lower Compression (Orange)
      3: [74, 196, 115],  // Watertight Base Cap (Green)
      4: [148, 163, 184], // Conformal Transition Frustum (Slate)
      5: [215, 222, 232]  // Conical Rocket Booster (Silver)
    };

    for (let idx = 0; idx < order.length; idx++) {
      const f = order[idx];
      const [i1, i2, i3] = faces[f];
      const p1 = projVerts[i1], p2 = projVerts[i2], p3 = projVerts[i3];
      const n = faceNormals[f];
      const shade = Math.max(0.32, Math.min(1.0, 0.58 + 0.35 * n[2] - 0.22 * n[0] + 0.15 * n[1]));
      const rgb = baseColors[facePatch[f]] || [180, 180, 180];

      ctx.beginPath();
      ctx.moveTo(p1[0], p1[1]);
      ctx.lineTo(p2[0], p2[1]);
      ctx.lineTo(p3[0], p3[1]);
      ctx.closePath();
      ctx.fillStyle = `rgb(${Math.round(rgb[0] * shade)}, ${Math.round(rgb[1] * shade)}, ${Math.round(rgb[2] * shade)})`;
      ctx.fill();
      if (showWire) {
        ctx.strokeStyle = "rgba(15, 23, 42, 0.45)";
        ctx.lineWidth = 0.55;
        ctx.stroke();
      }
    }

    // Optional Attached Shock Wave Sheet
    if (showShock) {
      ctx.strokeStyle = "rgba(244, 63, 94, 0.45)";
      ctx.lineWidth = 1.1;
      for (let i = 0; i < model.NxWave; i += 4) {
        ctx.beginPath();
        for (let k = 0; k < model.Ny; k++) {
          const p = projectPt([model.X_shk[i][k], model.Y_shk[i][k], model.Z_shk[i][k]]);
          if (k === 0) ctx.moveTo(p[0], p[1]);
          else ctx.lineTo(p[0], p[1]);
        }
        ctx.stroke();
      }
    }

    // Highlight Leading Edge Curve
    ctx.strokeStyle = "#f8fafc";
    ctx.lineWidth = 2.0;
    ctx.beginPath();
    for (let k = 0; k < model.Ny; k++) {
      const p = projectPt([model.X_low[0][k], model.Y_low[0][k], model.Z_low[0][k]]);
      if (k === 0) ctx.moveTo(p[0], p[1]);
      else ctx.lineTo(p[0], p[1]);
    }
    ctx.stroke();

    // Optional Ventral Concave Double-Bounce Coupling Rays
    if (showRays && model.concavePairs.length > 0) {
      ctx.strokeStyle = "rgba(250, 204, 21, 0.78)";
      ctx.lineWidth = 1.25;
      const nShow = Math.min(35, model.concavePairs.length);
      for (let m = 0; m < nShow; m++) {
        const cp = model.concavePairs[m];
        const pA = projectPt(cp.rP);
        const pB = projectPt(cp.rQ);
        ctx.beginPath();
        ctx.moveTo(pA[0], pA[1]);
        ctx.lineTo(pB[0], pB[1]);
        ctx.stroke();
      }
    }
  }

  // 7. Canvas 2: Aft Cross-Section & Conformal Frustum Morphing View
  const canvasSec = document.getElementById("canvas-sections");
  function renderCrossSections(model) {
    const rect = canvasSec.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    canvasSec.width = rect.width * dpr;
    canvasSec.height = rect.height * dpr;
    const ctx = canvasSec.getContext("2d");
    ctx.scale(dpr, dpr);
    const W = rect.width, H = rect.height;

    ctx.fillStyle = "#090e1a";
    ctx.fillRect(0, 0, W, H);

    const bMax = Math.max(model.cfg.span * 1.15, model.cfg.bradius * 1.3);
    const scale = Math.min(W / (2.2 * bMax), H / (1.6 * bMax));
    const cx = W * 0.5, cy = H * 0.36;

    function toXY(y, z) {
      return [cx + y * scale, cy - z * scale];
    }

    // Draw grid lines
    ctx.strokeStyle = "rgba(148, 163, 184, 0.14)";
    ctx.lineWidth = 1;
    for (let y = -3.0; y <= 3.0; y += 0.5) {
      const [px] = toXY(y, 0);
      ctx.beginPath(); ctx.moveTo(px, 0); ctx.lineTo(px, H); ctx.stroke();
    }
    for (let z = -2.0; z <= 1.0; z += 0.25) {
      const [, py] = toXY(0, z);
      ctx.beginPath(); ctx.moveTo(0, py); ctx.lineTo(W, py); ctx.stroke();
    }

    const rWave = model.NxWave - 1;
    // Fill Waverider Base Cross-Section at X = L
    ctx.beginPath();
    for (let k = 0; k < model.Ny; k++) {
      const [px, py] = toXY(model.Y_up[rWave][k], model.Z_up[rWave][k]);
      if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    for (let k = model.Ny - 1; k >= 0; k--) {
      const [px, py] = toXY(model.Y_low[rWave][k], model.Z_low[rWave][k]);
      ctx.lineTo(px, py);
    }
    ctx.closePath();
    ctx.fillStyle = "rgba(56, 189, 248, 0.16)";
    ctx.fill();

    // Draw Upper TE (Blue) & Lower Compression TE (Orange)
    ctx.strokeStyle = "#38bdf8";
    ctx.lineWidth = 2.4;
    ctx.beginPath();
    for (let k = 0; k < model.Ny; k++) {
      const [px, py] = toXY(model.Y_up[rWave][k], model.Z_up[rWave][k]);
      if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();

    ctx.strokeStyle = "#fb923c";
    ctx.lineWidth = 2.4;
    ctx.beginPath();
    for (let k = 0; k < model.Ny; k++) {
      const [px, py] = toXY(model.Y_low[rWave][k], model.Z_low[rWave][k]);
      if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();

    // Draw Shock Wave Trace at X = L (Dashed Red)
    ctx.save();
    ctx.setLineDash([6, 4]);
    ctx.strokeStyle = "#f43f5e";
    ctx.lineWidth = 1.8;
    ctx.beginPath();
    for (let k = 0; k < model.Ny; k++) {
      const [px, py] = toXY(model.Y_shk[rWave][k], model.Z_shk[rWave][k]);
      if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();
    ctx.restore();

    // If Conformal Frustum is enabled, draw intermediate morphing contours + Booster Ring
    if (model.NxFrustEnd > model.NxWave) {
      const rMid = Math.round(0.5 * (rWave + model.NxFrustEnd - 1));
      const rEnd = model.NxFrustEnd - 1;
      const morphRows = [
        { r: rMid, color: "rgba(192, 132, 252, 0.75)", width: 1.6 },
        { r: rEnd, color: "#4ade80", width: 2.4 }
      ];
      for (const m of morphRows) {
        ctx.strokeStyle = m.color;
        ctx.lineWidth = m.width;
        ctx.beginPath();
        for (let k = 0; k < model.Ny; k++) {
          const [px, py] = toXY(model.Y_up[m.r][k], model.Z_up[m.r][k]);
          if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
        }
        for (let k = model.Ny - 1; k >= 0; k--) {
          const [px, py] = toXY(model.Y_low[m.r][k], model.Z_low[m.r][k]);
          ctx.lineTo(px, py);
        }
        ctx.closePath();
        ctx.stroke();
      }
    }
  }

  // 8. Canvas 3 & 4: Multi-Band UHF-to-X-Band RCS Charts
  function drawLineChart(canvas, xVec, seriesList, xLabel, yLabel, yMin, yMax) {
    const rect = canvas.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    canvas.width = rect.width * dpr;
    canvas.height = rect.height * dpr;
    const ctx = canvas.getContext("2d");
    ctx.scale(dpr, dpr);
    const W = rect.width, H = rect.height;

    ctx.fillStyle = "#090e1a";
    ctx.fillRect(0, 0, W, H);

    const padL = 48, padR = 18, padT = 28, padB = 36;
    const plotW = W - padL - padR, plotH = H - padT - padB;
    const xMin = xVec[0], xMax = xVec[xVec.length - 1];

    function toPx(x, y) {
      const clampedY = Math.max(yMin, Math.min(yMax, y));
      const px = padL + ((x - xMin) / (xMax - xMin)) * plotW;
      const py = padT + (1.0 - (clampedY - yMin) / (yMax - yMin)) * plotH;
      return [px, py];
    }

    // Grid & Axes
    ctx.strokeStyle = "rgba(148, 163, 184, 0.16)";
    ctx.fillStyle = "#94a3b8";
    ctx.font = "11px ui-monospace, monospace";
    for (let y = yMin; y <= yMax; y += 10) {
      const [, py] = toPx(xMin, y);
      ctx.beginPath(); ctx.moveTo(padL, py); ctx.lineTo(W - padR, py); ctx.stroke();
      ctx.fillText(`${y}`, 10, py + 4);
    }

    const xStep = (xMax - xMin) <= 180 ? 30 : 30;
    for (let x = xMin; x <= xMax; x += xStep) {
      const [px] = toPx(x, yMin);
      ctx.beginPath(); ctx.moveTo(px, padT); ctx.lineTo(px, H - padB); ctx.stroke();
      ctx.fillText(`${x}°`, px - 10, H - 14);
    }

    // Plot Series
    for (const s of seriesList) {
      ctx.save();
      if (s.dashed) ctx.setLineDash([5, 4]);
      ctx.strokeStyle = s.color;
      ctx.lineWidth = s.width || 2.0;
      ctx.beginPath();
      for (let i = 0; i < xVec.length; i++) {
        const [px, py] = toPx(xVec[i], s.data[i]);
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.stroke();
      ctx.restore();
    }

    // Legend bar at top
    let lx = padL + 6;
    ctx.font = "11px -apple-system, sans-serif";
    for (const s of seriesList) {
      ctx.fillStyle = s.color;
      ctx.fillRect(lx, 10, 14, 3);
      ctx.fillStyle = "#e2e8f0";
      ctx.fillText(s.label, lx + 18, 15);
      lx += ctx.measureText(s.label).width + 32;
    }
  }

  function updateRCSAndPlots() {
    const elLook = parseFloat(document.getElementById("rng-rcsel").value);
    const includeMB = document.getElementById("chk-multibounce").checked;

    document.getElementById("hint-az-el").innerHTML =
      `Evaluated at Look Elevation &theta;<sub>el</sub> = ${elLook}&deg; (${elLook < 0 ? "Ventral Cavity" : "Dorsal Canopy"})`;

    // 1. Azimuthal Sweep (0 to 180 deg) across UHF (0.5), L (1.5), S (3.0), C (6.0), X (10.0 GHz)
    const azVec = [];
    for (let a = 0; a <= 180; a += 3) azVec.push(a);

    const bands = [
      { f: 0.5, label: "UHF 0.5GHz", color: "#38bdf8" },
      { f: 1.5, label: "L 1.5GHz",   color: "#4ade80" },
      { f: 3.0, label: "S 3.0GHz",   color: "#fbbf24" },
      { f: 6.0, label: "C 6.0GHz",   color: "#fb923c" },
      { f: 10.0, label: "X 10GHz",   color: "#c084fc" }
    ];

    const azSeries = bands.map((b) => {
      const raw = azVec.map((az) => evaluateRCSDirection(currentModel, b.f, az, elLook, includeMB).total_dBsm);
      // 3-point smoothing
      const smooth = raw.map((v, i) => {
        const prev = raw[Math.max(0, i - 1)];
        const next = raw[Math.min(raw.length - 1, i + 1)];
        return 0.25 * prev + 0.5 * v + 0.25 * next;
      });
      return { label: b.label, color: b.color, data: smooth, width: 2.0 };
    });

    drawLineChart(
      document.getElementById("canvas-rcs-az"),
      azVec, azSeries, "Azimuth [deg]", "RCS [dBsm]", -35, 25
    );

    // 2. Elevation Cut (-90 deg Ventral to +90 deg Dorsal) at Oblique Azimuth = 45 deg
    const elVec = [];
    for (let e = -90; e <= 90; e += 3) elVec.push(e);

    const elBands = [
      { f: 0.5, name: "UHF 0.5GHz", color: "#38bdf8" },
      { f: 3.0, name: "S 3.0GHz",   color: "#fbbf24" },
      { f: 10.0, name: "X 10GHz",   color: "#c084fc" }
    ];

    const elSeries = [];
    let maxMBGain = 0.0;
    let ventralXSum = 0.0, dorsalXSum = 0.0, ventralUHFSum = 0.0, nSec = 0;

    elBands.forEach((b) => {
      const totArr = [];
      const sngArr = [];
      elVec.forEach((el) => {
        const res = evaluateRCSDirection(currentModel, b.f, 45.0, el, includeMB);
        totArr.push(res.total_dBsm);
        sngArr.push(res.single_dBsm);
        const gain = res.total_dBsm - res.single_dBsm;
        if (gain > maxMBGain) maxMBGain = gain;
      });
      elSeries.push({ label: `${b.name} Total`, color: b.color, data: totArr, width: 2.2, dashed: false });
      elSeries.push({ label: `${b.name} 1-Bounce`, color: b.color, data: sngArr, width: 1.3, dashed: true });
    });

    for (let az = 25; az <= 75; az += 10) {
      ventralXSum += evaluateRCSDirection(currentModel, 10.0, az, -22.0, includeMB).total_m2;
      dorsalXSum  += evaluateRCSDirection(currentModel, 10.0, az, +22.0, includeMB).total_m2;
      ventralUHFSum += evaluateRCSDirection(currentModel, 0.5, az, -22.0, includeMB).total_m2;
      nSec++;
    }

    drawLineChart(
      document.getElementById("canvas-rcs-el"),
      elVec, elSeries, "Elevation [deg]", "RCS [dBsm]", -35, 35
    );

    const vX_dB = 10 * Math.log10(ventralXSum / nSec);
    const dX_dB = 10 * Math.log10(dorsalXSum / nSec);
    const vUHF_dB = 10 * Math.log10(ventralUHFSum / nSec);

    document.getElementById("m-pairs").textContent = `${currentModel.totalConcavePairs} Pairs`;
    document.getElementById("m-mbgain").textContent = `Max Double-Bounce +${maxMBGain.toFixed(1)} dB`;
    document.getElementById("m-rcs-x").textContent = `${vX_dB >= 0 ? "+" : ""}${vX_dB.toFixed(1)} dBsm`;
    document.getElementById("m-rcs-diff").textContent =
      `Dorsal: ${dX_dB >= 0 ? "+" : ""}${dX_dB.toFixed(1)} | UHF: ${vUHF_dB.toFixed(1)} dBsm`;
  }

  // 9. Read UI Controls & Rebuild Model
  function readConfigFromUI() {
    return {
      method: document.getElementById("sel-method").value,
      mach: parseFloat(document.getElementById("rng-mach").value),
      beta: parseFloat(document.getElementById("rng-beta").value),
      length: parseFloat(document.getElementById("rng-length").value),
      span: parseFloat(document.getElementById("rng-span").value),
      power: parseFloat(document.getElementById("rng-power").value),
      anhedral: parseFloat(document.getElementById("rng-anhedral").value),
      canopy: parseFloat(document.getElementById("rng-canopy").value),
      frustum: document.getElementById("chk-frustum").checked,
      flength: parseFloat(document.getElementById("rng-flength").value),
      bradius: parseFloat(document.getElementById("rng-bradius").value),
      blend: document.getElementById("sel-blend").value,
      rocket: document.getElementById("chk-rocket").checked,
      rcone: parseFloat(document.getElementById("rng-rcone").value)
    };
  }

  function updateAll() {
    const cfg = readConfigFromUI();
    document.getElementById("val-mach").textContent = cfg.mach.toFixed(1);
    document.getElementById("val-beta").innerHTML = `${cfg.beta.toFixed(1)}&deg;`;
    document.getElementById("val-length").textContent = `${cfg.length.toFixed(2)} m`;
    document.getElementById("val-span").textContent = `${cfg.span.toFixed(2)} m`;
    document.getElementById("val-power").textContent = cfg.power.toFixed(2);
    document.getElementById("val-anhedral").innerHTML = `${cfg.anhedral.toFixed(1)}&deg;`;
    document.getElementById("val-canopy").textContent = cfg.canopy.toFixed(3);
    document.getElementById("val-flength").textContent = `${cfg.flength.toFixed(2)} m`;
    document.getElementById("val-bradius").textContent = `${cfg.bradius.toFixed(2)} m`;
    document.getElementById("val-rcone").innerHTML = `${cfg.rcone.toFixed(1)}&deg;`;

    const elVal = parseFloat(document.getElementById("rng-rcsel").value);
    document.getElementById("val-rcsel").innerHTML =
      `${elVal >= 0 ? "+" : ""}${elVal}&deg; (${elVal < 0 ? "Ventral" : "Dorsal"})`;

    currentModel = buildWaveriderModel(cfg);

    // Update Telemetry Ribbon
    document.getElementById("m-topology").textContent = `${currentModel.boundaryEdges} Open Edges`;
    document.getElementById("m-faces").innerHTML =
      `${currentModel.faces.length.toLocaleString()} Faces &middot; &#967; = ${currentModel.euler}`;
    document.getElementById("m-volume").innerHTML = `${currentModel.enclosedVol.toFixed(2)} m&sup3;`;
    document.getElementById("m-etav").innerHTML =
      `&eta;<sub>V</sub> = ${currentModel.etaV.toFixed(3)} &middot; S<sub>plan</sub> = ${currentModel.planformArea.toFixed(1)} m&sup2;`;
    document.getElementById("m-ld").textContent = currentModel.LD.toFixed(2);
    document.getElementById("m-clcd").innerHTML =
      `C<sub>L</sub> = ${currentModel.CL.toFixed(3)} &middot; C<sub>D</sub> = ${currentModel.CD.toFixed(4)}`;

    render3DViewport(currentModel);
    renderCrossSections(currentModel);
    updateRCSAndPlots();
  }

  // 10. STL & OBJ Exporters + MATLAB Snippet Generator
  function downloadBlob(content, filename, mimeType) {
    const blob = new Blob([content], { type: mimeType });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }

  document.getElementById("btn-export-stl").addEventListener("click", () => {
    if (!currentModel) return;
    const { vertices, faces, faceNormals, cfg } = currentModel;
    const lines = [`solid waverider_${cfg.method}`];
    for (let f = 0; f < faces.length; f++) {
      const n = faceNormals[f];
      const [i1, i2, i3] = faces[f];
      const v1 = vertices[i1], v2 = vertices[i2], v3 = vertices[i3];
      lines.push(`  facet normal ${n[0].toExponential(6)} ${n[1].toExponential(6)} ${n[2].toExponential(6)}`);
      lines.push("    outer loop");
      lines.push(`      vertex ${v1[0].toExponential(6)} ${v1[1].toExponential(6)} ${v1[2].toExponential(6)}`);
      lines.push(`      vertex ${v2[0].toExponential(6)} ${v2[1].toExponential(6)} ${v2[2].toExponential(6)}`);
      lines.push(`      vertex ${v3[0].toExponential(6)} ${v3[1].toExponential(6)} ${v3[2].toExponential(6)}`);
      lines.push("    endloop");
      lines.push("  endfacet");
    }
    lines.push(`endsolid waverider_${cfg.method}`);
    downloadBlob(lines.join("\n"), `waverider_${cfg.method}_watertight.stl`, "model/stl");
  });

  document.getElementById("btn-export-obj").addEventListener("click", () => {
    if (!currentModel) return;
    const { vertices, faces, facePatch, cfg } = currentModel;
    const lines = [
      "# Watertight Hypersonic Waverider Glide Body",
      `# Vertices: ${vertices.length} | Faces: ${faces.length}`,
      ""
    ];
    for (const v of vertices) {
      lines.push(`v ${v[0].toFixed(6)} ${v[1].toFixed(6)} ${v[2].toFixed(6)}`);
    }
    const groups = {
      1: "UpperSurface",
      2: "LowerCompressionSurface",
      4: "ConformalTransitionFrustum",
      5: "ConicalRocketBooster",
      3: "AftBaseCap"
    };
    for (const [idStr, name] of Object.entries(groups)) {
      const id = parseInt(idStr, 10);
      lines.push(`\ng ${name}`);
      for (let f = 0; f < faces.length; f++) {
        if (facePatch[f] === id) {
          const [i1, i2, i3] = faces[f];
          lines.push(`f ${i1 + 1} ${i2 + 1} ${i3 + 1}`);
        }
      }
    }
    downloadBlob(lines.join("\n"), `waverider_${cfg.method}_watertight.obj`, "model/obj");
  });

  // Modals: MATLAB Code & Figure Gallery
  const dlgMatlab = document.getElementById("dlg-matlab");
  const dlgGallery = document.getElementById("dlg-gallery");

  document.getElementById("btn-open-matlab").addEventListener("click", () => {
    const cfg = readConfigFromUI();
    const code = [
      "setup_waverider();",
      "",
      "cfg = waverider.WaveriderConfig( ...",
      `    'Method', '${cfg.method}', ...`,
      `    'Mach', ${cfg.mach.toFixed(2)}, ...`,
      `    'ShockAngleDeg', ${cfg.beta.toFixed(2)}, ...`,
      `    'Length', ${cfg.length.toFixed(2)}, ...`,
      `    'HalfSpan', ${cfg.span.toFixed(2)}, ...`,
      `    'PlanformPower', ${cfg.power.toFixed(2)}, ...`,
      `    'LeadingEdgeAnhedralDeg', ${cfg.anhedral.toFixed(2)}, ...`,
      `    'UpperSurfaceType', '${cfg.canopy > 0 ? "canopy" : "freestream"}', ...`,
      `    'CanopyHeightRatio', ${cfg.canopy.toFixed(3)}, ...`,
      `    'EnableAftFrustum', ${cfg.frustum}, ...`,
      `    'FrustumLength', ${cfg.flength.toFixed(2)}, ...`,
      `    'BoosterRadius', ${cfg.bradius.toFixed(2)}, ...`,
      `    'FrustumBlendLaw', '${cfg.blend}', ...`,
      `    'EnableRocketStage', ${cfg.rocket}, ...`,
      `    'RocketConeAngleDeg', ${cfg.rcone.toFixed(2)});`,
      "",
      "body = waverider.WaveriderBody(cfg);",
      "body.printSummary();",
      "body.plot();",
      "body.plotRCS();",
      "body.exportSTL('output/custom_waverider.stl');"
    ].join("\n");
    document.getElementById("matlab-code-output").textContent = code;
    dlgMatlab.showModal();
  });

  document.getElementById("btn-copy-matlab").addEventListener("click", () => {
    const text = document.getElementById("matlab-code-output").textContent;
    navigator.clipboard.writeText(text);
  });

  document.getElementById("btn-open-gallery").addEventListener("click", () => {
    dlgGallery.showModal();
  });

  // Preset Buttons
  document.querySelectorAll("[data-preset]").forEach((btn) => {
    btn.addEventListener("click", () => {
      const p = PRESETS[btn.getAttribute("data-preset")];
      if (!p) return;
      document.getElementById("sel-method").value = p.method;
      document.getElementById("rng-mach").value = p.mach;
      document.getElementById("rng-beta").value = p.beta;
      document.getElementById("rng-length").value = p.length;
      document.getElementById("rng-span").value = p.span;
      document.getElementById("rng-power").value = p.power;
      document.getElementById("rng-anhedral").value = p.anhedral;
      document.getElementById("rng-canopy").value = p.canopy;
      document.getElementById("chk-frustum").checked = p.frustum;
      document.getElementById("rng-flength").value = p.flength;
      document.getElementById("rng-bradius").value = p.bradius;
      document.getElementById("sel-blend").value = p.blend;
      document.getElementById("chk-rocket").checked = p.rocket;
      document.getElementById("rng-rcone").value = p.rcone;
      rangeInputs.forEach(syncSliderVars);
      updateAll();
    });
  });

  // Bind all form controls
  document.querySelectorAll("input, select").forEach((el) => {
    el.addEventListener("input", updateAll);
    el.addEventListener("change", updateAll);
  });

  // 3D Orbit Camera Interaction
  canvas3D.addEventListener("pointerdown", (e) => {
    camera.dragging = true;
    camera.lastX = e.clientX;
    camera.lastY = e.clientY;
    canvas3D.setPointerCapture(e.pointerId);
  });
  canvas3D.addEventListener("pointermove", (e) => {
    if (!camera.dragging) return;
    const dx = e.clientX - camera.lastX;
    const dy = e.clientY - camera.lastY;
    camera.lastX = e.clientX;
    camera.lastY = e.clientY;
    camera.yaw += dx * 0.0085;
    camera.pitch = Math.max(-1.35, Math.min(1.35, camera.pitch + dy * 0.0085));
    if (currentModel) render3DViewport(currentModel);
  });
  canvas3D.addEventListener("pointerup", () => {
    camera.dragging = false;
  });
  canvas3D.addEventListener("wheel", (e) => {
    e.preventDefault();
    camera.zoom = Math.max(0.45, Math.min(2.6, camera.zoom * (e.deltaY < 0 ? 1.08 : 0.92)));
    if (currentModel) render3DViewport(currentModel);
  }, { passive: false });
  canvas3D.addEventListener("dblclick", () => {
    camera.yaw = -2.32;
    camera.pitch = 0.38;
    camera.zoom = 1.0;
    if (currentModel) render3DViewport(currentModel);
  });

  window.addEventListener("resize", () => {
    if (!currentModel) return;
    render3DViewport(currentModel);
    renderCrossSections(currentModel);
    updateRCSAndPlots();
  });

  // Initialize
  updateAll();
})();
