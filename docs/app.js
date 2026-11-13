/* ============================================================================
 * Hypersonic Parametric Watertight Waverider & UHF-to-X-Band RCS Web Studio
 * UI Controller & 6-Panel Canvas Renderer powered by window.WaveriderPhysics
 * (docs/physics.js).
 * ============================================================================ */

(function () {
  "use strict";

  const phys = window.WaveriderPhysics;
  if (!phys) {
    throw new Error("WaveriderPhysics (physics.js) must be loaded before app.js.");
  }

  // 1. Modern Web Guidance Fallbacks: Range Slider CSS Variables, Dialog Light-Dismiss & Copy Button
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

  // Accessible Copy Button handler (modern-web-guidance 'copy-button'):
  // Register only ONE path (command event OR click fallback) so a click never copies twice,
  // announce in external sibling [role="status"], and append \u00A0 on repeated messages for VoiceOver.
  const copyTimers = new WeakMap();
  let lastAnnouncedCopyMsg = "";

  async function performCopy(btn, targetEl) {
    if (!btn || !targetEl) return;
    const tip = btn.nextElementSibling;
    const text = "value" in targetEl ? targetEl.value : targetEl.textContent;
    let ok = false;
    try {
      await navigator.clipboard.writeText(text);
      ok = true;
    } catch {
      ok = false;
    }

    clearTimeout(copyTimers.get(btn));
    delete btn.dataset.copied;
    delete btn.dataset.failed;
    btn.dataset[ok ? "copied" : "failed"] = "";

    if (tip) {
      let msg = ok ? "Copied!" : "Copy failed";
      if (msg === lastAnnouncedCopyMsg) msg += "\u00A0";
      lastAnnouncedCopyMsg = msg;
      tip.textContent = msg;
    }

    copyTimers.set(
      btn,
      setTimeout(() => {
        delete btn.dataset.copied;
        delete btn.dataset.failed;
        if (tip) tip.textContent = "";
      }, 2000)
    );
  }

  if ("CommandEvent" in window) {
    document.addEventListener(
      "command",
      (e) => {
        if (e.command === "--copy") {
          performCopy(e.source, e.target);
        }
      },
      true
    );
  } else {
    document.addEventListener("click", (e) => {
      const btn = e.target.closest('button[command="--copy"]');
      if (!btn) return;
      const targetId = btn.getAttribute("commandfor");
      const targetEl = targetId ? document.getElementById(targetId) : null;
      performCopy(btn, targetEl);
    });
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
      arch: "waverider",
      method: "osculating_cone", mach: 7.0, beta: 13.5, length: 5.5, span: 1.75,
      power: 0.62, anhedral: 8.5, canopy: 0.045, frustum: true, flength: 1.65,
      bradius: 0.72, blend: "cubic", rocket: true, rcone: 4.5
    },
    rv_monoconic: {
      arch: "rv",
      rvBodyType: "mono_conic", rvNoseType: "spherical", rvBluntness: 0.16,
      rvMach: 12.0, rvLength: 2.60, rvCone1: 8.0, rvCone2: 8.0, rvCone3: 8.0, rvFrac1: 0.45,
      rvAftType: "flat", rvAftDepth: 0.35, rvAftPower: 2.0,
      frustum: false, flength: 1.20, bradius: 0.55, blend: "cubic", rocket: false, rcone: 4.0
    },
    rv_biconic: {
      arch: "rv",
      rvBodyType: "bi_conic", rvNoseType: "spherical", rvBluntness: 0.22,
      rvMach: 14.0, rvLength: 3.10, rvCone1: 11.5, rvCone2: 6.0, rvCone3: 6.0, rvFrac1: 0.42,
      rvAftType: "concave", rvAftDepth: 0.38, rvAftPower: 2.0,
      frustum: false, flength: 1.20, bradius: 0.55, blend: "cubic", rocket: false, rcone: 4.0
    },
    rv_triconic_flare: {
      arch: "rv",
      rvBodyType: "tri_conic", rvNoseType: "ogive", rvBluntness: 0.25,
      rvMach: 10.0, rvLength: 3.40, rvCone1: 9.5, rvCone2: 3.5, rvCone3: 13.5, rvFrac1: 0.38,
      rvAftType: "concave", rvAftDepth: 0.42, rvAftPower: 2.6,
      frustum: false, flength: 1.20, bradius: 0.55, blend: "cubic", rocket: false, rcone: 4.0
    },
    rv_sharp_conic: {
      arch: "rv",
      rvBodyType: "mono_conic", rvNoseType: "conical", rvBluntness: 0.16,
      rvMach: 15.0, rvLength: 2.80, rvCone1: 7.5, rvCone2: 7.5, rvCone3: 7.5, rvFrac1: 0.50,
      rvAftType: "flat", rvAftDepth: 0.35, rvAftPower: 2.0,
      frustum: false, flength: 1.20, bradius: 0.55, blend: "cubic", rocket: false, rcone: 4.0
    }
  };

  let activeArch = "waverider";
  let hasInitializedRVTab = false;

  function getCheckedRadio(name, fallback) {
    const checked = document.querySelector(`input[name="${name}"]:checked`);
    return checked ? checked.value : fallback;
  }

  function setCheckedRadio(name, val) {
    const target = document.querySelector(`input[name="${name}"][value="${val}"]`);
    if (target) target.checked = true;
  }

  function syncRVConditionalInertStates() {
    const bodyType = getCheckedRadio("rv-body-type", "bi_conic");
    const noseType = getCheckedRadio("rv-nose-type", "spherical");
    const aftType = getCheckedRadio("rv-aft-type", "concave");
    const frustumOn = document.getElementById("chk-frustum").checked;

    const bluntGroup = document.getElementById("rv-bluntness-group");
    const stage2Group = document.getElementById("rv-stage2-group");
    const stage3Group = document.getElementById("rv-stage3-group");
    const aftCavGroup = document.getElementById("rv-aft-cavity-group");
    const bluntLbl = document.getElementById("lbl-rv-bluntness");

    if (bluntGroup) bluntGroup.inert = noseType === "conical";
    if (bluntLbl) {
      bluntLbl.innerHTML = noseType === "ogive"
        ? "Ogive Curvature / Length Ratio"
        : "Nose Bluntness (R<sub>N</sub> / R<sub>base</sub>)";
    }
    if (stage2Group) stage2Group.inert = bodyType === "mono_conic";
    if (stage3Group) stage3Group.inert = bodyType !== "tri_conic";
    if (aftCavGroup) aftCavGroup.inert = aftType !== "concave" || frustumOn;
  }

  const archTablist = document.querySelector(".arch-tabs");
  const archTabs = archTablist ? Array.from(archTablist.querySelectorAll('[role="tab"]')) : [];

  function selectVehicleTab(nextTab, moveFocus = false) {
    if (!nextTab) return;
    const nextArch = nextTab.id === "tab-rv" ? "rv" : "waverider";
    archTabs.forEach((tab) => {
      const selected = tab === nextTab;
      tab.setAttribute("aria-selected", String(selected));
      tab.tabIndex = selected ? 0 : -1;
      const panel = document.getElementById(tab.getAttribute("aria-controls"));
      if (panel) {
        panel.hidden = selected ? false : "until-found";
      }
    });
    if (nextArch === "rv" && !hasInitializedRVTab) {
      hasInitializedRVTab = true;
      // Default standalone RV exposes its flat/concave aft plane unless user re-enables booster mount
      document.getElementById("chk-frustum").checked = false;
      document.getElementById("chk-rocket").checked = false;
    }
    activeArch = nextArch;
    syncRVConditionalInertStates();
    if (moveFocus) nextTab.focus();
    scheduleUpdate("geometry");
  }

  if (archTablist) {
    archTabs.forEach((tab) => {
      tab.addEventListener("click", () => selectVehicleTab(tab, false));
      const panel = document.getElementById(tab.getAttribute("aria-controls"));
      if (panel) {
        panel.addEventListener("beforematch", () => selectVehicleTab(tab, false));
      }
    });
    archTablist.addEventListener("keydown", (e) => {
      const idx = archTabs.indexOf(e.target);
      if (idx === -1) return;
      const dir = getComputedStyle(archTablist).direction === "rtl" ? -1 : 1;
      let nextIdx = -1;
      if (e.key === "ArrowRight") nextIdx = (idx + dir + archTabs.length) % archTabs.length;
      else if (e.key === "ArrowLeft") nextIdx = (idx - dir + archTabs.length) % archTabs.length;
      else if (e.key === "Home") nextIdx = 0;
      else if (e.key === "End") nextIdx = archTabs.length - 1;
      if (nextIdx !== -1) {
        e.preventDefault();
        selectVehicleTab(archTabs[nextIdx], true);
      }
    });
  }

  let currentModel = null;
  let currentRTIData = null;
  let camera = { yaw: -2.32, pitch: 0.38, zoom: 1.0, dragging: false, lastX: 0, lastY: 0 };

  // Dynamic physical guardrails on Shock Angle (beta > Mach angle mu) and Bandwidth (B <= 0.95 * f_c)
  function enforcePhysicalSliderBounds() {
    const rngMach = document.getElementById("rng-mach");
    const rngBeta = document.getElementById("rng-beta");
    const mach = parseFloat(rngMach.value);
    const muDeg = (Math.asin(Math.min(1.0 / mach, 0.999)) * 180.0) / Math.PI;
    const minBeta = Math.max(10.0, Math.ceil((muDeg + 1.0) * 2.0) / 2.0);
    const maxBeta = Math.max(24.0, minBeta + 6.0);
    rngBeta.min = minBeta.toFixed(1);
    rngBeta.max = maxBeta.toFixed(1);
    if (parseFloat(rngBeta.value) < minBeta) {
      rngBeta.value = minBeta.toFixed(1);
    } else if (parseFloat(rngBeta.value) > maxBeta) {
      rngBeta.value = maxBeta.toFixed(1);
    }
    syncSliderVars(rngBeta);

    const rngFc = document.getElementById("rng-fc");
    const rngBw = document.getElementById("rng-bw");
    const fcGHz = parseFloat(rngFc.value);
    const maxBwMHz = Math.max(50, Math.min(4000, Math.floor((fcGHz * 950.0) / 50.0) * 50));
    rngBw.max = String(maxBwMHz);
    if (parseFloat(rngBw.value) > maxBwMHz) {
      rngBw.value = String(maxBwMHz);
    }
    syncSliderVars(rngBw);
  }

  // 3. Canvas 1: Interactive 3D Viewport
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

    const { vertices, faces, facePatch, faceNormals, cfg } = model;
    const xMaxTotal = model.xMax;
    const xMid = 0.5 * (model.xMin + xMaxTotal);
    const zMid = 0.5 * (model.zMin + model.zMax);
    const scale = Math.min(W, H) * 0.13 * camera.zoom * (5.5 / Math.max(xMaxTotal, 2.5));

    const cy = Math.cos(camera.yaw), sy = Math.sin(camera.yaw);
    const cp = Math.cos(camera.pitch), sp = Math.sin(camera.pitch);

    function projectPt(pt) {
      const x = pt[0] - xMid, y = pt[1], z = pt[2] - zMid;
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
      ctx.strokeStyle = "rgba(244, 63, 94, 0.48)";
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

    // Highlight Leading Edge Curve (Waverider) or Stage Break / Aft Rim Rings (RV)
    ctx.strokeStyle = "#f8fafc";
    ctx.lineWidth = 2.0;
    if (model.isAxisymmetric) {
      const rAft = model.NxWave - 1;
      ctx.beginPath();
      for (let k = 0; k < model.Ny; k++) {
        const p = projectPt([model.X_up[rAft][k], model.Y_up[rAft][k], model.Z_up[rAft][k]]);
        if (k === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
      }
      for (let k = model.Ny - 1; k >= 0; k--) {
        const p = projectPt([model.X_low[rAft][k], model.Y_low[rAft][k], model.Z_low[rAft][k]]);
        ctx.lineTo(p[0], p[1]);
      }
      ctx.closePath();
      ctx.stroke();
    } else {
      ctx.beginPath();
      for (let k = 0; k < model.Ny; k++) {
        const p = projectPt([model.X_low[0][k], model.Y_low[0][k], model.Z_low[0][k]]);
        if (k === 0) ctx.moveTo(p[0], p[1]);
        else ctx.lineTo(p[0], p[1]);
      }
      ctx.stroke();
    }

    // Optional Ventral / Cavity Concave Double-Bounce Coupling Rays
    if (showRays && model.concavePairs.length > 0) {
      ctx.strokeStyle = "rgba(250, 204, 21, 0.78)";
      ctx.lineWidth = 1.25;
      const nShow = Math.min(35, model.concavePairs.length);
      for (let m = 0; m < nShow; m++) {
        const cpPair = model.concavePairs[m];
        const pA = projectPt(cpPair.rP);
        const pB = projectPt(cpPair.rQ);
        ctx.beginPath();
        ctx.moveTo(pA[0], pA[1]);
        ctx.lineTo(pB[0], pB[1]);
        ctx.stroke();
      }
    }

    // 3D Reference Coordinate Bounding Frame & Labeled Tick Marks in Meters
    const bHalfRef = model.isAxisymmetric && model.rvMeta ? model.rvMeta.rBase_m : cfg.span;
    const yAxisRef = -Math.max(bHalfRef, Math.abs(model.yMin)) * 1.18;
    const zAxisRef = model.zMin - 0.18;

    ctx.strokeStyle = "rgba(148, 163, 184, 0.45)";
    ctx.fillStyle = "#cbd5e1";
    ctx.font = "11px ui-monospace, SFMono-Regular, Menlo, monospace";
    ctx.lineWidth = 1.1;

    // Axial X-axis ruler [0 .. xMaxTotal]
    const pX0 = projectPt([0, yAxisRef, zAxisRef]);
    const pX1 = projectPt([xMaxTotal, yAxisRef, zAxisRef]);
    ctx.beginPath(); ctx.moveTo(pX0[0], pX0[1]); ctx.lineTo(pX1[0], pX1[1]); ctx.stroke();
    const xTickStep = xMaxTotal > 6 ? 2.0 : 1.0;
    for (let xv = 0; xv <= xMaxTotal + 0.01; xv += xTickStep) {
      const pt = projectPt([xv, yAxisRef, zAxisRef]);
      ctx.fillRect(pt[0] - 1.5, pt[1] - 1.5, 3, 3);
      ctx.fillText(`${xv.toFixed(0)}m`, pt[0] - 8, pt[1] + 14);
    }
    const pXMid = projectPt([0.5 * xMaxTotal, yAxisRef, zAxisRef]);
    ctx.font = "bold 11.5px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.fillStyle = "#38bdf8";
    ctx.fillText("Axial Length X [m]", pXMid[0] - 42, pXMid[1] + 28);

    // Spanwise / Radial Y-axis ruler at aft plane [-bHalfRef .. +bHalfRef]
    ctx.font = "11px ui-monospace, SFMono-Regular, Menlo, monospace";
    ctx.fillStyle = "#cbd5e1";
    const pY0 = projectPt([xMaxTotal, -bHalfRef, zAxisRef]);
    const pY1 = projectPt([xMaxTotal, bHalfRef, zAxisRef]);
    ctx.beginPath(); ctx.moveTo(pY0[0], pY0[1]); ctx.lineTo(pY1[0], pY1[1]); ctx.stroke();
    [-bHalfRef, 0, bHalfRef].forEach((yv) => {
      const pt = projectPt([xMaxTotal, yv, zAxisRef]);
      ctx.fillRect(pt[0] - 1.5, pt[1] - 1.5, 3, 3);
      ctx.fillText(`${yv >= 0 ? "+" : ""}${yv.toFixed(2)}m`, pt[0] + 5, pt[1] + 4);
    });
    ctx.font = "bold 11.5px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.fillStyle = "#4ade80";
    ctx.fillText(model.isAxisymmetric ? "Radius Y [m]" : "Span Y [m]", pY1[0] + 6, pY1[1] - 8);

    // 3D Orientation Triad in Bottom-Left Corner
    const ox = 54, oy = H - 48, triLen = 34;
    const axesTriad = [
      { dir: [1, 0, 0], label: "+X (Aft)", color: "#f87171" },
      { dir: [0, 1, 0], label: "+Y (Span)", color: "#4ade80" },
      { dir: [0, 0, 1], label: "+Z (Up)", color: "#38bdf8" }
    ];
    axesTriad.forEach((ax) => {
      const rx = ax.dir[0] * cy - ax.dir[1] * sy;
      const ry = ax.dir[0] * sy + ax.dir[1] * cy;
      const rz = ax.dir[2];
      const dx = rx * triLen;
      const dy = (ry * sp - rz * cp) * triLen;
      ctx.strokeStyle = ax.color;
      ctx.fillStyle = ax.color;
      ctx.lineWidth = 2.0;
      ctx.beginPath();
      ctx.moveTo(ox, oy);
      ctx.lineTo(ox + dx, oy + dy);
      ctx.stroke();
      ctx.font = "bold 10.5px ui-monospace, monospace";
      ctx.fillText(ax.label, ox + dx * 1.18 - 14, oy + dy * 1.18 + 3);
    });

    // Top-Left Surface Patch Legend
    const legItems = [
      { label: "Upper Surface", color: "rgb(52, 148, 226)" },
      { label: "Lower Compression", color: "rgb(232, 112, 54)" },
      ...(cfg.frustum ? [{ label: "Conformal Frustum", color: "rgb(148, 163, 184)" }] : []),
      ...(cfg.frustum && cfg.rocket ? [{ label: "Rocket Booster", color: "rgb(215, 222, 232)" }] : []),
      { label: "Watertight Base", color: "rgb(74, 196, 115)" }
    ];
    let ly = 16;
    ctx.font = "11px -apple-system, BlinkMacSystemFont, sans-serif";
    legItems.forEach((item) => {
      ctx.fillStyle = item.color;
      ctx.fillRect(14, ly - 8, 12, 10);
      ctx.fillStyle = "#e2e8f0";
      ctx.fillText(item.label, 32, ly);
      ly += 17;
    });
  }

  // 4. Canvas 2: Aft Cross-Section (Waverider) OR Meridional X-R Blueprint Profile (RV)
  const canvasSec = document.getElementById("canvas-sections");
  function renderRVMeridionalBlueprint(ctx, W, H, model) {
    const padL = 64, padR = 24, padT = 34, padB = 48;
    const plotW = Math.max(W - padL - padR, 100);
    const plotH = Math.max(H - padT - padB, 100);

    const rv = model.rvMeta;
    const midCol = Math.floor(model.Ny / 2);
    const xMaxPlot = Math.max(model.xMax * 1.06, 1.0);
    const xMinPlot = -0.05 * xMaxPlot;
    let rMaxAbs = rv.rBase_m * 1.35;
    for (let i = 0; i < rv.rShockProfile.length; i++) {
      rMaxAbs = Math.max(rMaxAbs, rv.rShockProfile[i] * 1.05);
    }
    for (let i = 0; i < model.NxTotal; i++) {
      rMaxAbs = Math.max(rMaxAbs, Math.abs(model.Z_up[i][midCol]) * 1.15, Math.abs(model.Z_low[i][midCol]) * 1.15);
    }

    const scale = Math.min(plotW / (xMaxPlot - xMinPlot), plotH / (2.0 * rMaxAbs));
    const xSpanActual = plotW / scale;
    const xOrigin = 0.5 * (xMinPlot + xMaxPlot) - 0.5 * xSpanActual;
    const cy = padT + plotH * 0.5;

    function toPx(x, z) {
      return [padL + (x - xOrigin) * scale, cy - z * scale];
    }

    ctx.save();
    ctx.beginPath();
    ctx.rect(padL, padT, plotW, plotH);
    ctx.clip();

    const xStep = xMaxPlot > 5.0 ? 1.0 : 0.5;
    const zStep = rMaxAbs > 1.2 ? 0.5 : 0.25;
    ctx.strokeStyle = "rgba(148, 163, 184, 0.16)";
    ctx.lineWidth = 1;
    for (let xv = Math.ceil(xOrigin / xStep) * xStep; xv <= xOrigin + xSpanActual; xv += xStep) {
      const [px] = toPx(xv, 0);
      ctx.beginPath(); ctx.moveTo(px, padT); ctx.lineTo(px, padT + plotH); ctx.stroke();
    }
    const zLim = (plotH * 0.5) / scale;
    for (let zv = -Math.floor(zLim / zStep) * zStep; zv <= zLim + 1e-6; zv += zStep) {
      const [, py] = toPx(0, zv);
      ctx.beginPath(); ctx.moveTo(padL, py); ctx.lineTo(padL + plotW, py); ctx.stroke();
    }

    // Centerline axis Z = 0
    ctx.strokeStyle = "rgba(148, 163, 184, 0.35)";
    ctx.beginPath(); ctx.moveTo(padL, cy); ctx.lineTo(padL + plotW, cy); ctx.stroke();

    // Build aft bulkhead path (flat or recessed concave cavity)
    const hasAftExt = model.NxTotal > model.NxWave;
    const aftPts = [];
    const rTopAft = model.Z_up[model.NxTotal - 1][midCol];
    const rBotAft = model.Z_low[model.NxTotal - 1][midCol];
    const xRimAft = model.X_up[model.NxTotal - 1][midCol];
    if (!hasAftExt && rv.aftCavityDepth_m > 0) {
      const nDome = 29;
      for (let m = 0; m < nDome; m++) {
        const s = 1.0 - (2.0 * m) / (nDome - 1); // +1 (top) -> -1 (bottom)
        const u = Math.abs(s);
        const xD = xRimAft - rv.aftCavityDepth_m * (1.0 - Math.pow(u, rv.aftCavityPower));
        const zD = s >= 0 ? u * rTopAft : u * Math.abs(rBotAft) * -1.0;
        aftPts.push([xD, zD]);
      }
    } else {
      aftPts.push([xRimAft, rTopAft], [xRimAft, rBotAft]);
    }

    // Fill RV Meridional Solid Cross-Section
    ctx.beginPath();
    for (let i = 0; i < model.NxTotal; i++) {
      const [px, py] = toPx(model.X_up[i][midCol], model.Z_up[i][midCol]);
      if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    for (let m = 0; m < aftPts.length; m++) {
      const [px, py] = toPx(aftPts[m][0], aftPts[m][1]);
      ctx.lineTo(px, py);
    }
    for (let i = model.NxTotal - 1; i >= 0; i--) {
      const [px, py] = toPx(model.X_low[i][midCol], model.Z_low[i][midCol]);
      ctx.lineTo(px, py);
    }
    ctx.closePath();
    ctx.fillStyle = "rgba(56, 189, 248, 0.16)";
    ctx.fill();

    // Stage break & nose tangency markers
    ctx.save();
    ctx.setLineDash([4, 4]);
    ctx.lineWidth = 1.2;
    for (let b = 0; b < rv.stageBreaksX.length - 1; b++) {
      const xb = rv.stageBreaksX[b];
      const rb = rv.stageBreaksR[b];
      if (xb <= 1e-4) continue;
      const [px, pyTop] = toPx(xb, rb);
      const [, pyBot] = toPx(xb, -rb);
      ctx.strokeStyle = b === 0 ? "rgba(250, 204, 21, 0.75)" : "rgba(192, 132, 252, 0.75)";
      ctx.beginPath(); ctx.moveTo(px, pyTop); ctx.lineTo(px, pyBot); ctx.stroke();
      ctx.fillStyle = b === 0 ? "#fde047" : "#d8b4fe";
      ctx.font = "10px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.textAlign = "center";
      ctx.fillText(b === 0 ? "Tangency" : `Stage ${b + 1}`, px, Math.max(padT + 12, pyTop - 5));
    }
    ctx.restore();

    // Upper (+R) & Lower (-R) Outer Mold Line
    ctx.strokeStyle = "#38bdf8";
    ctx.lineWidth = 2.4;
    ctx.beginPath();
    for (let i = 0; i < model.NxTotal; i++) {
      const [px, py] = toPx(model.X_up[i][midCol], model.Z_up[i][midCol]);
      if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();

    ctx.strokeStyle = "#fb923c";
    ctx.lineWidth = 2.4;
    ctx.beginPath();
    for (let i = 0; i < model.NxTotal; i++) {
      const [px, py] = toPx(model.X_low[i][midCol], model.Z_low[i][midCol]);
      if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();

    // Aft Bulkhead / Concave Cavity Contour (Green)
    ctx.strokeStyle = "#4ade80";
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    for (let m = 0; m < aftPts.length; m++) {
      const [px, py] = toPx(aftPts[m][0], aftPts[m][1]);
      if (m === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();

    // Bow / Conical Shock Envelope (+R_shk and -R_shk)
    ctx.save();
    ctx.setLineDash([6, 4]);
    ctx.strokeStyle = "#f43f5e";
    ctx.lineWidth = 1.7;
    for (const sgn of [+1, -1]) {
      ctx.beginPath();
      for (let i = 0; i < model.NxWave; i++) {
        const [px, py] = toPx(model.X_shk[i][midCol], sgn * Math.abs(model.Z_shk[i][midCol]));
        if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }
    ctx.restore();
    ctx.restore();

    ctx.strokeStyle = "#475569";
    ctx.lineWidth = 1.2;
    ctx.strokeRect(padL, padT, plotW, plotH);

    ctx.fillStyle = "#cbd5e1";
    ctx.font = "11px ui-monospace, SFMono-Regular, Menlo, monospace";
    ctx.textAlign = "center";
    for (let xv = Math.ceil(xOrigin / xStep) * xStep; xv <= xOrigin + xSpanActual; xv += xStep) {
      const [px] = toPx(xv, 0);
      if (px >= padL + 10 && px <= padL + plotW - 10) {
        ctx.fillText(`${xv.toFixed(1)}`, px, padT + plotH + 16);
      }
    }
    ctx.textAlign = "right";
    for (let zv = -Math.floor(zLim / zStep) * zStep; zv <= zLim + 1e-6; zv += zStep) {
      const [, py] = toPx(0, zv);
      if (py >= padT + 8 && py <= padT + plotH - 6) {
        ctx.fillText(`${zv.toFixed(2)}`, padL - 8, py + 4);
      }
    }

    ctx.fillStyle = "#f0f4fc";
    ctx.font = "bold 12px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText(`Axial Station X [m]  (D_base = ${(2.0 * rv.rBase_m).toFixed(2)} m, ${rv.aftCavityDepth_m > 0 && !hasAftExt ? `Concave d = ${(rv.aftCavityDepth_m * 100).toFixed(1)} cm` : "Flat Aft Plane"})`, padL + plotW * 0.5, H - 10);

    ctx.save();
    ctx.translate(16, padT + plotH * 0.5);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.fillText("Radial Contour \u00B1R(x) [m]", 0, 0);
    ctx.restore();

    const rvLegend = [
      { label: "Dorsal +R(x)", color: "#38bdf8" },
      { label: "Ventral -R(x)", color: "#fb923c" },
      { label: rv.aftCavityDepth_m > 0 && !hasAftExt ? "Concave Aft Cavity" : "Flat Aft Bulkhead", color: "#4ade80" },
      { label: "Shock Envelope", color: "#f43f5e" }
    ];
    let lx = padL + 4;
    ctx.textAlign = "left";
    ctx.font = "11px -apple-system, BlinkMacSystemFont, sans-serif";
    rvLegend.forEach((s) => {
      ctx.fillStyle = s.color;
      ctx.fillRect(lx, 12, 14, 3);
      ctx.fillStyle = "#e2e8f0";
      ctx.fillText(s.label, lx + 18, 17);
      lx += ctx.measureText(s.label).width + 26;
    });
  }

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

    const titleEl = document.getElementById("title-panel-2");
    const hintEl = document.getElementById("hint-panel-2");
    if (model.isAxisymmetric && model.rvMeta) {
      if (titleEl) {
        titleEl.innerHTML = "2. Meridional RV Blueprint (Nose Tangency, Conic Stages &amp; Aft Bulkhead)";
      }
      if (hintEl) {
        hintEl.innerHTML = "Blue/Orange = Outer Skin &middot; Green = Flat/Concave Aft Plane &middot; Dashed = Shock";
      }
      renderRVMeridionalBlueprint(ctx, W, H, model);
      return;
    }
    if (titleEl) {
      titleEl.innerHTML = "2. Cross-Sectional Morphing (Waverider Base X=L &rarr; Booster Ring &rarr; Rocket)";
    }
    if (hintEl) {
      hintEl.innerHTML = "Blue/Orange = Waverider &middot; Purple/Green = Frustum Morph &middot; Dashed = Shock";
    }

    const padL = 64, padR = 24, padT = 34, padB = 48;
    const plotW = Math.max(W - padL - padR, 100);
    const plotH = Math.max(H - padT - padB, 100);

    const rWave = model.NxWave - 1;
    let zMinData = -0.8, zMaxData = 0.25;
    for (let r = 0; r < model.NxTotal; r++) {
      for (let k = 0; k < model.Ny; k++) {
        zMinData = Math.min(zMinData, model.Z_low[r][k]);
        zMaxData = Math.max(zMaxData, model.Z_up[r][k]);
      }
    }
    for (let k = 0; k < model.Ny; k++) {
      zMinData = Math.min(zMinData, model.Z_shk[rWave][k]);
    }

    const yMaxAbs = Math.max(model.cfg.span * 1.12, model.cfg.bradius * 1.25, 1.0);
    const zSpan = Math.max(zMaxData - zMinData, 0.8);
    const zMid = 0.5 * (zMaxData + zMinData);

    const scale = Math.min(plotW / (2.0 * yMaxAbs), plotH / (zSpan * 1.22));
    const cx = padL + plotW * 0.5;
    const cy = padT + plotH * 0.5 + zMid * scale;

    const yLeftVal = -(plotW * 0.5) / scale;
    const yRightVal = (plotW * 0.5) / scale;
    const zTopVal = (cy - padT) / scale;
    const zBotVal = (cy - (padT + plotH)) / scale;

    function toXY(y, z) {
      return [cx + y * scale, cy - z * scale];
    }

    ctx.save();
    ctx.beginPath();
    ctx.rect(padL, padT, plotW, plotH);
    ctx.clip();

    const yStep = yMaxAbs > 2.2 ? 1.0 : 0.5;
    const zStep = (zTopVal - zBotVal) > 2.0 ? 0.5 : 0.25;

    ctx.strokeStyle = "rgba(148, 163, 184, 0.16)";
    ctx.lineWidth = 1;
    for (let y = Math.ceil(yLeftVal / yStep) * yStep; y <= yRightVal; y += yStep) {
      const [px] = toXY(y, 0);
      ctx.beginPath(); ctx.moveTo(px, padT); ctx.lineTo(px, padT + plotH); ctx.stroke();
    }
    for (let z = Math.ceil(zBotVal / zStep) * zStep; z <= zTopVal; z += zStep) {
      const [, py] = toXY(0, z);
      ctx.beginPath(); ctx.moveTo(padL, py); ctx.lineTo(padL + plotW, py); ctx.stroke();
    }

    ctx.strokeStyle = "rgba(148, 163, 184, 0.35)";
    ctx.beginPath(); ctx.moveTo(cx, padT); ctx.lineTo(cx, padT + plotH); ctx.stroke();

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

    // Upper TE (Blue) & Lower Compression TE (Orange)
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

    // Shock Wave Trace at X = L (Dashed Red)
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

    if (model.NxFrustEnd > model.NxWave) {
      const rMid = Math.round(0.5 * (rWave + model.NxFrustEnd - 1));
      const rEnd = model.NxFrustEnd - 1;
      const morphRows = [
        { r: rMid, color: "rgba(192, 132, 252, 0.85)", width: 1.7 },
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
    ctx.restore();

    ctx.strokeStyle = "#475569";
    ctx.lineWidth = 1.2;
    ctx.strokeRect(padL, padT, plotW, plotH);

    ctx.fillStyle = "#cbd5e1";
    ctx.font = "11px ui-monospace, SFMono-Regular, Menlo, monospace";
    ctx.textAlign = "center";
    for (let y = Math.ceil(yLeftVal / yStep) * yStep; y <= yRightVal; y += yStep) {
      const [px] = toXY(y, 0);
      if (px >= padL + 10 && px <= padL + plotW - 10) {
        ctx.fillText(`${y.toFixed(1)}`, px, padT + plotH + 16);
      }
    }

    ctx.textAlign = "right";
    for (let z = Math.ceil(zBotVal / zStep) * zStep; z <= zTopVal; z += zStep) {
      const [, py] = toXY(0, z);
      if (py >= padT + 8 && py <= padT + plotH - 6) {
        ctx.fillText(`${z.toFixed(2)}`, padL - 8, py + 4);
      }
    }

    ctx.fillStyle = "#f0f4fc";
    ctx.font = "bold 12px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Spanwise Coordinate Y [m]  (Port -Y  \u2194  Starboard +Y)", padL + plotW * 0.5, H - 10);

    ctx.save();
    ctx.translate(16, padT + plotH * 0.5);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.fillText("Vertical Coordinate Z [m]", 0, 0);
    ctx.restore();

    const secLegend = [
      { label: "Upper TE (X=L)", color: "#38bdf8" },
      { label: "Lower TE (X=L)", color: "#fb923c" },
      { label: "Shock Wave", color: "#f43f5e" },
      ...(model.NxFrustEnd > model.NxWave ? [
        { label: "Frustum Mid (\u03C4=0.5)", color: "#c084fc" },
        { label: "Booster Ring (\u03C4=1.0)", color: "#4ade80" }
      ] : [])
    ];
    let lx = padL + 4;
    ctx.textAlign = "left";
    ctx.font = "11px -apple-system, BlinkMacSystemFont, sans-serif";
    secLegend.forEach((s) => {
      ctx.fillStyle = s.color;
      ctx.fillRect(lx, 12, 14, 3);
      ctx.fillStyle = "#e2e8f0";
      ctx.fillText(s.label, lx + 18, 17);
      lx += ctx.measureText(s.label).width + 30;
    });
  }

  // 5. Canvas 3, 4 & 6: Labeled Line Charts (Angular RCS & Down-Range HRRP)
  function drawLineChart(canvas, xVec, seriesList, xLabel, yLabel, yMin, yMax, vLines = [], xUnit = "°", xStepOpt = null) {
    const rect = canvas.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    canvas.width = rect.width * dpr;
    canvas.height = rect.height * dpr;
    const ctx = canvas.getContext("2d");
    ctx.scale(dpr, dpr);
    const W = rect.width, H = rect.height;

    ctx.fillStyle = "#090e1a";
    ctx.fillRect(0, 0, W, H);

    const padL = 64, padR = 22, padT = 34, padB = 48;
    const plotW = Math.max(W - padL - padR, 100);
    const plotH = Math.max(H - padT - padB, 100);
    const xMin = xVec[0], xMax = xVec[xVec.length - 1];

    function toPx(x, y) {
      const clampedY = Math.max(yMin, Math.min(yMax, y));
      const px = padL + ((x - xMin) / (xMax - xMin)) * plotW;
      const py = padT + (1.0 - (clampedY - yMin) / (yMax - yMin)) * plotH;
      return [px, py];
    }

    ctx.strokeStyle = "rgba(148, 163, 184, 0.16)";
    ctx.lineWidth = 1;
    ctx.fillStyle = "#cbd5e1";
    ctx.font = "11px ui-monospace, SFMono-Regular, Menlo, monospace";

    ctx.textAlign = "right";
    for (let y = yMin; y <= yMax; y += 10) {
      const [, py] = toPx(xMin, y);
      ctx.beginPath(); ctx.moveTo(padL, py); ctx.lineTo(padL + plotW, py); ctx.stroke();
      ctx.fillText(`${y >= 0 ? "+" : ""}${y}`, padL - 8, py + 4);
    }

    ctx.textAlign = "center";
    const xStep = xStepOpt || 30;
    const xStart = Math.ceil(xMin / xStep) * xStep;
    for (let x = xStart; x <= xMax + 1e-6; x += xStep) {
      const [px] = toPx(x, yMin);
      if (px < padL - 2 || px > padL + plotW + 2) continue;
      ctx.beginPath(); ctx.moveTo(px, padT); ctx.lineTo(px, padT + plotH); ctx.stroke();
      const xTxt = xUnit === " m" ? `${x >= 0 ? "+" : ""}${x.toFixed(1)} m` : `${Math.round(x)}${xUnit}`;
      ctx.fillText(xTxt, px, padT + plotH + 16);
    }

    vLines.forEach((vl) => {
      if (vl.x < xMin || vl.x > xMax) return;
      const [px] = toPx(vl.x, yMin);
      ctx.save();
      ctx.setLineDash([4, 4]);
      ctx.strokeStyle = vl.color || "rgba(248, 250, 252, 0.55)";
      ctx.lineWidth = 1.3;
      ctx.beginPath();
      ctx.moveTo(px, padT);
      ctx.lineTo(px, padT + plotH);
      ctx.stroke();
      if (vl.label) {
        ctx.fillStyle = vl.color || "#e2e8f0";
        ctx.font = "10.5px -apple-system, BlinkMacSystemFont, sans-serif";
        ctx.textAlign = px > padL + plotW * 0.72 ? "right" : "left";
        const dxTxt = px > padL + plotW * 0.72 ? -5 : 5;
        ctx.fillText(vl.label, px + dxTxt, padT + (vl.yOffset || 14));
      }
      ctx.restore();
    });

    ctx.save();
    ctx.beginPath();
    ctx.rect(padL, padT, plotW, plotH);
    ctx.clip();
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
    ctx.restore();

    ctx.strokeStyle = "#475569";
    ctx.lineWidth = 1.2;
    ctx.strokeRect(padL, padT, plotW, plotH);

    ctx.fillStyle = "#f0f4fc";
    ctx.font = "bold 12px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText(xLabel, padL + plotW * 0.5, H - 10);

    ctx.save();
    ctx.translate(16, padT + plotH * 0.5);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.fillText(yLabel, 0, 0);
    ctx.restore();

    let lx = padL + 4;
    ctx.textAlign = "left";
    ctx.font = "11px -apple-system, BlinkMacSystemFont, sans-serif";
    for (const s of seriesList) {
      ctx.fillStyle = s.color;
      ctx.fillRect(lx, 12, 14, 3);
      ctx.fillStyle = "#e2e8f0";
      ctx.fillText(s.label, lx + 18, 17);
      lx += ctx.measureText(s.label).width + 22;
    }
  }

  function turboColormapRGB(t) {
    const x = Math.max(0.0, Math.min(1.0, t));
    const r = Math.round(255 * Math.min(1.0, Math.max(0.0, 1.55 * x - 0.15 + 0.35 * Math.pow(x, 3))));
    const g = Math.round(255 * Math.min(1.0, Math.max(0.0, Math.sin(Math.PI * Math.pow(x, 0.88)) * 0.95 + 0.18 * x)));
    const b = Math.round(255 * Math.min(1.0, Math.max(0.0, 0.22 + 0.85 * Math.sin(Math.PI * (1.0 - x * 0.82)) * (1.0 - 0.65 * x))));
    return [r, g, b];
  }

  // 6. Canvas 5: 2D Range-Angle-Intensity Sinogram (RTI)
  const canvasRTI = document.getElementById("canvas-rti");
  function renderRangeAngleIntensity(rtiData, azCutDeg) {
    const rect = canvasRTI.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    canvasRTI.width = rect.width * dpr;
    canvasRTI.height = rect.height * dpr;
    const ctx = canvasRTI.getContext("2d");
    ctx.scale(dpr, dpr);
    const W = rect.width, H = rect.height;

    ctx.fillStyle = "#090e1a";
    ctx.fillRect(0, 0, W, H);

    const padL = 64, padR = 76, padT = 34, padB = 48;
    const plotW = Math.max(W - padL - padR, 100);
    const plotH = Math.max(H - padT - padB, 100);

    const { azVec, rangeVec, RTI_total_dBsm, R_nose, R_base, R_aft } = rtiData;
    const nAz = azVec.length;
    const nR = rangeVec.length;
    const rMin = rangeVec[0], rMax = rangeVec[nR - 1];
    const azMin = azVec[0], azMax = azVec[nAz - 1];
    const cMin = -42.0, cMax = 22.0;

    const offCanvas = document.createElement("canvas");
    offCanvas.width = nR;
    offCanvas.height = nAz;
    const offCtx = offCanvas.getContext("2d");
    const imgData = offCtx.createImageData(nR, nAz);
    const data = imgData.data;

    for (let ia = 0; ia < nAz; ia++) {
      const rowIdx = nAz - 1 - ia;
      const rowVals = RTI_total_dBsm[ia];
      for (let ir = 0; ir < nR; ir++) {
        const norm = (rowVals[ir] - cMin) / (cMax - cMin);
        const [r, g, b] = turboColormapRGB(norm);
        const pIdx = (rowIdx * nR + ir) * 4;
        data[pIdx] = r;
        data[pIdx + 1] = g;
        data[pIdx + 2] = b;
        data[pIdx + 3] = 255;
      }
    }
    offCtx.putImageData(imgData, 0, 0);

    ctx.save();
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(offCanvas, padL, padT, plotW, plotH);
    ctx.restore();

    function toPx(rVal, azVal) {
      const px = padL + ((rVal - rMin) / (rMax - rMin)) * plotW;
      const py = padT + (1.0 - (azVal - azMin) / (azMax - azMin)) * plotH;
      return [px, py];
    }

    ctx.save();
    ctx.beginPath();
    ctx.rect(padL, padT, plotW, plotH);
    ctx.clip();

    ctx.strokeStyle = "rgba(255, 255, 255, 0.14)";
    ctx.lineWidth = 1;
    for (let r = Math.ceil(rMin); r <= Math.floor(rMax); r += 1.0) {
      const [px] = toPx(r, 0);
      ctx.beginPath(); ctx.moveTo(px, padT); ctx.lineTo(px, padT + plotH); ctx.stroke();
    }
    for (let a = 0; a <= 180; a += 30) {
      const [, py] = toPx(0, a);
      ctx.beginPath(); ctx.moveTo(padL, py); ctx.lineTo(padL + plotW, py); ctx.stroke();
    }

    ctx.save();
    ctx.setLineDash([5, 4]);
    ctx.strokeStyle = "#38bdf8";
    ctx.lineWidth = 1.7;
    ctx.beginPath();
    for (let ia = 0; ia < nAz; ia++) {
      const [px, py] = toPx(R_nose[ia], azVec[ia]);
      if (ia === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();

    ctx.strokeStyle = "#fbbf24";
    ctx.beginPath();
    for (let ia = 0; ia < nAz; ia++) {
      const [px, py] = toPx(R_base[ia], azVec[ia]);
      if (ia === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();

    if (Math.abs(R_aft[0] - R_base[0]) > 0.05) {
      ctx.setLineDash([2, 3]);
      ctx.strokeStyle = "#4ade80";
      ctx.beginPath();
      for (let ia = 0; ia < nAz; ia++) {
        const [px, py] = toPx(R_aft[ia], azVec[ia]);
        if (ia === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }
    ctx.restore();

    const [, pyCut] = toPx(0, azCutDeg);
    ctx.save();
    ctx.setLineDash([6, 4]);
    ctx.strokeStyle = "#ffffff";
    ctx.lineWidth = 1.8;
    ctx.beginPath();
    ctx.moveTo(padL, pyCut);
    ctx.lineTo(padL + plotW, pyCut);
    ctx.stroke();
    ctx.fillStyle = "rgba(9, 14, 26, 0.78)";
    ctx.fillRect(padL + 6, Math.max(padT + 2, pyCut - 18), 118, 15);
    ctx.fillStyle = "#ffffff";
    ctx.font = "bold 10.5px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.textAlign = "left";
    ctx.fillText(`HRRP Cut \u03C8 = ${azCutDeg}\u00B0`, padL + 10, Math.max(padT + 13, pyCut - 7));
    ctx.restore();

    ctx.restore();

    ctx.strokeStyle = "#475569";
    ctx.lineWidth = 1.2;
    ctx.strokeRect(padL, padT, plotW, plotH);

    ctx.fillStyle = "#cbd5e1";
    ctx.font = "11px ui-monospace, SFMono-Regular, Menlo, monospace";
    ctx.textAlign = "center";
    const rStep = (rMax - rMin) > 8.5 ? 2.0 : 1.0;
    for (let r = Math.ceil(rMin / rStep) * rStep; r <= rMax; r += rStep) {
      const [px] = toPx(r, 0);
      if (px >= padL + 8 && px <= padL + plotW - 8) {
        ctx.fillText(`${r >= 0 ? "+" : ""}${r.toFixed(1)} m`, px, padT + plotH + 16);
      }
    }

    ctx.textAlign = "right";
    for (let a = 0; a <= 180; a += 30) {
      const [, py] = toPx(0, a);
      ctx.fillText(`${a}\u00B0`, padL - 8, py + 4);
    }

    ctx.fillStyle = "#f0f4fc";
    ctx.font = "bold 12px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Down-Range Coordinate R [m]  (Relative to Vehicle Center:  -R Toward Radar  \u2194  +R Away)", padL + plotW * 0.5, H - 10);

    ctx.save();
    ctx.translate(16, padT + plotH * 0.5);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.fillText("Azimuth Aspect Angle \u03C8 [deg]", 0, 0);
    ctx.restore();

    const rtiLegend = [
      { label: "Nose Locus", color: "#38bdf8" },
      { label: "Waverider Base (X=L)", color: "#fbbf24" },
      ...(Math.abs(R_aft[0] - R_base[0]) > 0.05 ? [{ label: "Booster Base", color: "#4ade80" }] : []),
      { label: "Selected HRRP Cut", color: "#ffffff" }
    ];
    let lx = padL + 4;
    ctx.textAlign = "left";
    ctx.font = "11px -apple-system, BlinkMacSystemFont, sans-serif";
    rtiLegend.forEach((s) => {
      ctx.fillStyle = s.color;
      ctx.fillRect(lx, 12, 14, 3);
      ctx.fillStyle = "#e2e8f0";
      ctx.fillText(s.label, lx + 18, 17);
      lx += ctx.measureText(s.label).width + 24;
    });

    const cbX = padL + plotW + 12;
    const cbW = 12;
    for (let iy = 0; iy < plotH; iy++) {
      const norm = 1.0 - iy / plotH;
      const [r, g, b] = turboColormapRGB(norm);
      ctx.fillStyle = `rgb(${r},${g},${b})`;
      ctx.fillRect(cbX, padT + iy, cbW, 1.5);
    }
    ctx.strokeStyle = "#475569";
    ctx.strokeRect(cbX, padT, cbW, plotH);

    ctx.fillStyle = "#cbd5e1";
    ctx.font = "10px ui-monospace, SFMono-Regular, Menlo, monospace";
    ctx.textAlign = "left";
    for (const dbVal of [-35, -20, -5, 10, 20]) {
      const py = padT + (1.0 - (dbVal - cMin) / (cMax - cMin)) * plotH;
      ctx.fillText(`${dbVal >= 0 ? "+" : ""}${dbVal}`, cbX + cbW + 4, py + 3);
    }
    ctx.save();
    ctx.translate(W - 8, padT + plotH * 0.5);
    ctx.rotate(-Math.PI / 2);
    ctx.fillStyle = "#f0f4fc";
    ctx.font = "bold 11px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Intensity [dBsm]", 0, 0);
    ctx.restore();
  }

  // 7. Canvas 6: 1D High-Resolution Range Profile (HRRP) Cut & Delayed Cavity Bounce
  function renderHRRPChart(rtiData, azCutDeg) {
    const idx = Math.max(0, Math.min(rtiData.azVec.length - 1, Math.round(azCutDeg / 3)));
    const actualAz = rtiData.azVec[idx];
    const deltaR_cm = (rtiData.deltaR_m * 100.0).toFixed(1);

    const hrrpSeries = [
      { label: "Total Profile (PO + 2-Bounce)", color: "#38bdf8", data: rtiData.RTI_total_dBsm[idx], width: 2.3, dashed: false },
      { label: "Single-Bounce PO Only", color: "#fb923c", data: rtiData.RTI_single_dBsm[idx], width: 1.6, dashed: true },
      { label: "Concave 2-Bounce (+0.5 d_PQ Delay)", color: "#c084fc", data: rtiData.RTI_double_dBsm[idx], width: 2.0, dashed: false }
    ];

    const vLines = [
      { x: rtiData.R_nose[idx], label: "Nose Tip", color: "#38bdf8", yOffset: 14 },
      { x: rtiData.R_base[idx], label: "Waverider Base (X=L)", color: "#fbbf24", yOffset: 28 }
    ];
    if (Math.abs(rtiData.R_aft[idx] - rtiData.R_base[idx]) > 0.08) {
      vLines.push({ x: rtiData.R_aft[idx], label: "Booster Base", color: "#4ade80", yOffset: 42 });
    }

    const rStep = (rtiData.rangeVec[rtiData.rangeVec.length - 1] - rtiData.rangeVec[0]) > 8.5 ? 2.0 : 1.0;
    const polLabel = rtiData.polarization || "Mean";
    drawLineChart(
      document.getElementById("canvas-hrrp"),
      rtiData.rangeVec,
      hrrpSeries,
      `Down-Range Coordinate R [m]  (Pol: ${polLabel}, f_c = ${rtiData.freqGHz.toFixed(1)} GHz, B = ${rtiData.bwMHz.toFixed(0)} MHz, \u0394R = ${deltaR_cm} cm, \u03C8 = ${actualAz}\u00B0)`,
      `Range Profile \u03C3_${polLabel}(R) [dBsm]`,
      -45,
      25,
      vLines,
      " m",
      rStep
    );
  }

  // 8. Multi-Band RCS & Wideband RTI / HRRP Evaluation via WaveriderPhysics
  function updateRCSAndPlots() {
    enforcePhysicalSliderBounds();

    const pol = document.getElementById("sel-pol").value;
    const elLook = parseFloat(document.getElementById("rng-rcsel").value);
    const fcGHz = parseFloat(document.getElementById("rng-fc").value);
    const bwMHz = parseFloat(document.getElementById("rng-bw").value);
    const winType = document.getElementById("sel-window").value;
    const rtiAz = parseFloat(document.getElementById("rng-rtiaz").value);
    const includeMB = document.getElementById("chk-multibounce").checked;

    const deltaR_cm = (phys.C0 / (2.0 * bwMHz * 1e6)) * 100.0;
    const bandLabel = fcGHz < 1.0 ? "UHF" : fcGHz < 2.0 ? "L" : fcGHz < 4.0 ? "S" : fcGHz < 8.0 ? "C" : "X";

    document.getElementById("val-rcsel").innerHTML =
      `${elLook >= 0 ? "+" : ""}${elLook}&deg; (${elLook < 0 ? "Ventral" : "Dorsal"})`;
    document.getElementById("val-fc").textContent = `${fcGHz.toFixed(1)} GHz (${bandLabel})`;
    document.getElementById("val-bw").innerHTML = `${bwMHz.toFixed(0)} MHz (&Delta;R=${deltaR_cm.toFixed(1)} cm)`;
    document.getElementById("val-rtiaz").innerHTML = `${rtiAz.toFixed(0)}&deg;`;

    document.getElementById("m-rangeres").textContent = `${deltaR_cm.toFixed(1)} cm`;
    document.getElementById("m-bw-sub").textContent =
      `B = ${bwMHz.toFixed(0)} MHz @ ${fcGHz.toFixed(1)} GHz (${((bwMHz * 0.1) / fcGHz).toFixed(1)}% BW)`;

    document.getElementById("hint-az-el").innerHTML =
      `Pol: <strong>${pol}</strong> &middot; Look Elevation &theta;<sub>el</sub> = ${elLook}&deg; (${elLook < 0 ? "Ventral Cavity" : "Dorsal Canopy"})`;
    document.getElementById("hint-rti-bw").innerHTML =
      `Pol: <strong>${pol}</strong> &middot; f<sub>c</sub> = ${fcGHz.toFixed(1)} GHz &middot; B = ${bwMHz.toFixed(0)} MHz (&Delta;R = ${deltaR_cm.toFixed(1)} cm) &middot; Click Heatmap to Cut`;
    document.getElementById("hint-hrrp-cut").innerHTML =
      `Pol: <strong>${pol}</strong> &middot; Cut at &psi;<sub>cut</sub> = ${rtiAz.toFixed(0)}&deg;, &theta;<sub>el</sub> = ${elLook}&deg; &middot; Window: ${winType.toUpperCase()}`;

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
    const azFreqs = bands.map((b) => b.f);
    const azRaw = bands.map(() => new Float64Array(azVec.length));

    for (let ia = 0; ia < azVec.length; ia++) {
      const multi = phys.evaluateRCSMultiFreq(currentModel, azFreqs, azVec[ia], elLook, includeMB, pol);
      for (let ib = 0; ib < bands.length; ib++) {
        azRaw[ib][ia] = multi[ib].total_dBsm;
      }
    }

    const azSeries = bands.map((b, ib) => ({
      label: b.label,
      color: b.color,
      data: azRaw[ib],
      width: 2.0
    }));

    let azMarkers = [];
    if (currentModel.isAxisymmetric && currentModel.rvMeta) {
      const th1 = currentModel.rvMeta.coneAnglesDeg[0];
      const specAz = Math.max(0.0, Math.min(180.0, 90.0 - th1));
      azMarkers.push({ x: specAz, label: `Cone Specular (~${Math.round(specAz)}\u00B0)` });
      if (currentModel.rvMeta.aftCavityDepth_m > 0 && currentModel.NxTotal === currentModel.NxWave) {
        azMarkers.push({ x: 162, label: "Aft Cavity Sector", color: "#4ade80", yOffset: 28 });
      }
    } else {
      const leSweepNomDeg = (Math.atan2(currentModel.cfg.length, currentModel.cfg.span) * 180.0) / Math.PI;
      azMarkers.push({ x: leSweepNomDeg, label: `LE Specular (~${Math.round(leSweepNomDeg)}\u00B0)` });
    }
    drawLineChart(
      document.getElementById("canvas-rcs-az"),
      azVec,
      azSeries,
      "Azimuth Angle \u03C8 [deg]  (0\u00B0 = Nose-On, 90\u00B0 = Broadside, 180\u00B0 = Tail-On)",
      `Monostatic RCS \u03C3_${pol} [dBsm]`,
      -40,
      25,
      azMarkers
    );

    // 2. Elevation Cut (-90 deg Ventral to +90 deg Dorsal) at Oblique Azimuth = 45 deg
    const elVec = [];
    for (let e = -90; e <= 90; e += 3) elVec.push(e);

    const elBands = [
      { f: 0.5, name: "UHF 0.5GHz", color: "#38bdf8" },
      { f: 3.0, name: "S 3.0GHz",   color: "#fbbf24" },
      { f: 10.0, name: "X 10GHz",   color: "#c084fc" }
    ];
    const elFreqs = elBands.map((b) => b.f);
    const elTot = elBands.map(() => new Float64Array(elVec.length));
    const elSng = elBands.map(() => new Float64Array(elVec.length));
    let maxMBGain = 0.0;

    for (let ie = 0; ie < elVec.length; ie++) {
      const multi = phys.evaluateRCSMultiFreq(currentModel, elFreqs, 45.0, elVec[ie], includeMB, pol);
      for (let ib = 0; ib < elBands.length; ib++) {
        elTot[ib][ie] = multi[ib].total_dBsm;
        elSng[ib][ie] = multi[ib].single_dBsm;
        const gain = multi[ib].total_dBsm - multi[ib].single_dBsm;
        if (gain > maxMBGain) maxMBGain = gain;
      }
    }

    const elSeries = [];
    elBands.forEach((b, ib) => {
      elSeries.push({ label: `${b.name} Total`, color: b.color, data: elTot[ib], width: 2.2, dashed: false });
      elSeries.push({ label: `${b.name} 1-Bounce`, color: b.color, data: elSng[ib], width: 1.3, dashed: true });
    });

    let ventralXSum = 0.0, dorsalXSum = 0.0, ventralUHFSum = 0.0, nSec = 0;
    for (let az = 25; az <= 75; az += 10) {
      const vPair = phys.evaluateRCSMultiFreq(currentModel, [0.5, 10.0], az, -22.0, includeMB, pol);
      ventralUHFSum += vPair[0].total_m2;
      ventralXSum += vPair[1].total_m2;
      dorsalXSum += phys.evaluateRCSDirection(currentModel, 10.0, az, +22.0, includeMB, pol).total_m2;
      nSec++;
    }

    drawLineChart(
      document.getElementById("canvas-rcs-el"),
      elVec,
      elSeries,
      "Elevation Angle \u03B8_el [deg]  (< 0\u00B0 Ventral Concave Underside  |  > 0\u00B0 Dorsal Canopy)",
      `Monostatic RCS \u03C3_${pol} [dBsm]`,
      -40,
      35,
      [{ x: 0, label: "Horizon (\u03B8_el = 0\u00B0)" }]
    );

    const vX_dB = 10.0 * Math.log10(ventralXSum / nSec);
    const dX_dB = 10.0 * Math.log10(dorsalXSum / nSec);
    const vUHF_dB = 10.0 * Math.log10(ventralUHFSum / nSec);

    document.getElementById("m-pairs").textContent = `${currentModel.totalConcavePairs} Pairs`;
    document.getElementById("m-mbgain").textContent =
      pol === "HV" || pol === "VH" || pol === "OP"
        ? `Cavity Gain +${maxMBGain.toFixed(1)} dB (${pol})`
        : `Max Double-Bounce +${maxMBGain.toFixed(1)} dB (${pol})`;
    document.getElementById("m-rcs-x").textContent = `${vX_dB >= 0 ? "+" : ""}${vX_dB.toFixed(1)} dBsm (${pol})`;
    document.getElementById("m-rcs-diff").textContent =
      `Dorsal: ${dX_dB >= 0 ? "+" : ""}${dX_dB.toFixed(1)} | UHF: ${vUHF_dB.toFixed(1)} dBsm`;

    // 3. Compute & Render 2D Range-Angle-Intensity (Sinogram) and 1D High-Resolution Range Profile (HRRP)
    currentRTIData = phys.computeRangeAngleIntensity(currentModel, fcGHz, bwMHz, elLook, includeMB, winType, {
      polarization: pol
    });
    renderRangeAngleIntensity(currentRTIData, rtiAz);
    renderHRRPChart(currentRTIData, rtiAz);
  }

  // 9. Read UI Controls, Rebuild Model & Tiered rAF Update Scheduler
  function readConfigFromUI() {
    const frustum = document.getElementById("chk-frustum").checked;
    const flength = parseFloat(document.getElementById("rng-flength").value);
    const bradius = parseFloat(document.getElementById("rng-bradius").value);
    const blend = document.getElementById("sel-blend").value;
    const rocket = document.getElementById("chk-rocket").checked;
    const rcone = parseFloat(document.getElementById("rng-rcone").value);

    if (activeArch === "rv") {
      const rvBodyType = getCheckedRadio("rv-body-type", "bi_conic");
      const rvNoseType = getCheckedRadio("rv-nose-type", "spherical");
      const rvAftType = getCheckedRadio("rv-aft-type", "concave");
      const rvMach = parseFloat(document.getElementById("rng-rv-mach").value);
      const rvLength = parseFloat(document.getElementById("rng-rv-length").value);
      const rvBluntness = parseFloat(document.getElementById("rng-rv-bluntness").value);
      const c1 = parseFloat(document.getElementById("rng-rv-cone1").value);
      const c2 = parseFloat(document.getElementById("rng-rv-cone2").value);
      const c3 = parseFloat(document.getElementById("rng-rv-cone3").value);
      const f1 = parseFloat(document.getElementById("rng-rv-frac1").value);
      const rvAftDepth = parseFloat(document.getElementById("rng-rv-aftdepth").value);
      const rvAftPower = parseFloat(document.getElementById("rng-rv-aftpower").value);

      let rvStageFractions = [1.0, 0.0, 0.0];
      if (rvBodyType === "bi_conic") {
        rvStageFractions = [f1, Math.max(0.05, 1.0 - f1), 0.0];
      } else if (rvBodyType === "tri_conic") {
        const rem = Math.max(0.10, 1.0 - f1);
        rvStageFractions = [f1, 0.52 * rem, 0.48 * rem];
      }

      return {
        method: "rv",
        mach: rvMach,
        length: rvLength,
        rvBodyType,
        rvNoseType,
        rvBluntness,
        rvConeAngles: [c1, c2, c3],
        rvStageFractions,
        rvAftType,
        rvAftDepth,
        rvAftPower,
        frustum, flength, bradius, blend, rocket, rcone
      };
    }

    return {
      method: document.getElementById("sel-method").value,
      mach: parseFloat(document.getElementById("rng-mach").value),
      beta: parseFloat(document.getElementById("rng-beta").value),
      length: parseFloat(document.getElementById("rng-length").value),
      span: parseFloat(document.getElementById("rng-span").value),
      power: parseFloat(document.getElementById("rng-power").value),
      anhedral: parseFloat(document.getElementById("rng-anhedral").value),
      canopy: parseFloat(document.getElementById("rng-canopy").value),
      frustum, flength, bradius, blend, rocket, rcone
    };
  }

  function updateGeometryAndAll() {
    enforcePhysicalSliderBounds();
    syncRVConditionalInertStates();
    const cfg = readConfigFromUI();

    document.getElementById("val-mach").textContent = parseFloat(document.getElementById("rng-mach").value).toFixed(1);
    document.getElementById("val-beta").innerHTML = `${parseFloat(document.getElementById("rng-beta").value).toFixed(1)}&deg;`;
    document.getElementById("val-length").textContent = `${parseFloat(document.getElementById("rng-length").value).toFixed(2)} m`;
    document.getElementById("val-span").textContent = `${parseFloat(document.getElementById("rng-span").value).toFixed(2)} m`;
    document.getElementById("val-power").textContent = parseFloat(document.getElementById("rng-power").value).toFixed(2);
    document.getElementById("val-anhedral").innerHTML = `${parseFloat(document.getElementById("rng-anhedral").value).toFixed(1)}&deg;`;
    document.getElementById("val-canopy").textContent = parseFloat(document.getElementById("rng-canopy").value).toFixed(3);
    document.getElementById("val-flength").textContent = `${cfg.flength.toFixed(2)} m`;
    document.getElementById("val-bradius").textContent = `${cfg.bradius.toFixed(2)} m`;
    document.getElementById("val-rcone").innerHTML = `${cfg.rcone.toFixed(1)}&deg;`;

    document.getElementById("val-rv-mach").textContent = parseFloat(document.getElementById("rng-rv-mach").value).toFixed(1);
    document.getElementById("val-rv-length").textContent = `${parseFloat(document.getElementById("rng-rv-length").value).toFixed(2)} m`;
    document.getElementById("val-rv-bluntness").textContent = parseFloat(document.getElementById("rng-rv-bluntness").value).toFixed(2);
    document.getElementById("val-rv-cone1").innerHTML = `${parseFloat(document.getElementById("rng-rv-cone1").value).toFixed(1)}&deg;`;
    document.getElementById("val-rv-cone2").innerHTML = `${parseFloat(document.getElementById("rng-rv-cone2").value).toFixed(1)}&deg;`;
    document.getElementById("val-rv-cone3").innerHTML = `${parseFloat(document.getElementById("rng-rv-cone3").value).toFixed(1)}&deg;`;
    document.getElementById("val-rv-frac1").textContent = parseFloat(document.getElementById("rng-rv-frac1").value).toFixed(2);
    document.getElementById("val-rv-aftdepth").textContent = parseFloat(document.getElementById("rng-rv-aftdepth").value).toFixed(2);
    document.getElementById("val-rv-aftpower").textContent = parseFloat(document.getElementById("rng-rv-aftpower").value).toFixed(1);

    currentModel = phys.buildWaveriderModel(cfg);

    // Dynamic Watertight Topology Badge (uses icon + text + border style, never color alone)
    const badge = document.getElementById("watertight-badge");
    if (currentModel.isWatertight) {
      badge.className = "badge-watertight";
      badge.innerHTML = `&#10003; PASS &middot; 2-Manifold Watertight (&#967; = ${currentModel.euler})`;
    } else {
      badge.className = "badge-watertight is-fail";
      badge.innerHTML = `&#9888; FAIL &middot; Topology Alert (${currentModel.boundaryEdges} Open, &#967; = ${currentModel.euler})`;
    }

    // Update Telemetry Ribbon
    document.getElementById("m-topology").textContent = `${currentModel.boundaryEdges} Open Edges`;
    document.getElementById("m-faces").innerHTML =
      `${currentModel.faces.length.toLocaleString()} Faces &middot; &#967; = ${currentModel.euler}`;
    document.getElementById("m-volume").innerHTML = `${currentModel.enclosedVol.toFixed(2)} m&sup3;`;
    document.getElementById("m-etav").innerHTML =
      `&eta;<sub>V</sub> = ${currentModel.etaV.toFixed(3)} &middot; S<sub>plan</sub> = ${currentModel.planformArea.toFixed(1)} m&sup2;`;
    document.getElementById("m-ld").textContent =
      `${currentModel.LD_inviscid.toFixed(2)} / ${currentModel.LD_viscous_est.toFixed(2)}`;
    document.getElementById("m-clcd").innerHTML =
      `C<sub>L</sub> = ${currentModel.CL.toFixed(3)} &middot; C<sub>D,w</sub> = ${currentModel.CD_wave.toFixed(4)}`;

    render3DViewport(currentModel);
    renderCrossSections(currentModel);
    updateRCSAndPlots();
  }

  function updateCutOnly() {
    if (!currentRTIData) {
      updateRCSAndPlots();
      return;
    }
    const pol = document.getElementById("sel-pol").value;
    const rtiAz = parseFloat(document.getElementById("rng-rtiaz").value);
    const elLook = parseFloat(document.getElementById("rng-rcsel").value);
    const winType = document.getElementById("sel-window").value;
    document.getElementById("val-rtiaz").innerHTML = `${rtiAz.toFixed(0)}&deg;`;
    document.getElementById("hint-hrrp-cut").innerHTML =
      `Pol: <strong>${pol}</strong> &middot; Cut at &psi;<sub>cut</sub> = ${rtiAz.toFixed(0)}&deg;, &theta;<sub>el</sub> = ${elLook}&deg; &middot; Window: ${winType.toUpperCase()}`;
    renderRangeAngleIntensity(currentRTIData, rtiAz);
    renderHRRPChart(currentRTIData, rtiAz);
  }

  // Tiered rAF scheduler prevents duplicate work during fast slider drags
  const LEVEL_PRIORITY = { view3d: 1, cut: 2, rcs: 3, geometry: 4 };
  let pendingLevel = null;
  let rafHandle = 0;

  function scheduleUpdate(level) {
    if (!pendingLevel || LEVEL_PRIORITY[level] > LEVEL_PRIORITY[pendingLevel]) {
      pendingLevel = level;
    }
    if (!rafHandle) {
      rafHandle = requestAnimationFrame(() => {
        const task = pendingLevel;
        pendingLevel = null;
        rafHandle = 0;
        if (task === "geometry" || !currentModel) {
          updateGeometryAndAll();
        } else if (task === "rcs") {
          updateRCSAndPlots();
        } else if (task === "cut") {
          updateCutOnly();
        } else if (task === "view3d") {
          render3DViewport(currentModel);
        }
      });
    }
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
    const tag = cfg.method === "rv" ? `rv_${cfg.rvBodyType}_${cfg.rvNoseType}` : `waverider_${cfg.method}`;
    const lines = [`solid ${tag}`];
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
    lines.push(`endsolid ${tag}`);
    downloadBlob(lines.join("\n"), `${tag}_watertight.stl`, "model/stl");
  });

  document.getElementById("btn-export-obj").addEventListener("click", () => {
    if (!currentModel) return;
    const { vertices, faces, facePatch, cfg } = currentModel;
    const tag = cfg.method === "rv" ? `rv_${cfg.rvBodyType}_${cfg.rvNoseType}` : `waverider_${cfg.method}`;
    const lines = [
      cfg.method === "rv" ? "# Watertight Hypersonic Reentry Vehicle (RV)" : "# Watertight Hypersonic Waverider Glide Body",
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
    downloadBlob(lines.join("\n"), `${tag}_watertight.obj`, "model/obj");
  });

  // Modals: MATLAB Code & Figure Gallery
  const dlgMatlab = document.getElementById("dlg-matlab");
  const dlgGallery = document.getElementById("dlg-gallery");

  document.getElementById("btn-open-matlab").addEventListener("click", () => {
    enforcePhysicalSliderBounds();
    const cfg = readConfigFromUI();
    const fcGHz = parseFloat(document.getElementById("rng-fc").value);
    const bwMHz = parseFloat(document.getElementById("rng-bw").value);
    const elLook = parseFloat(document.getElementById("rng-rcsel").value);
    const rtiAz = parseFloat(document.getElementById("rng-rtiaz").value);
    const winType = document.getElementById("sel-window").value;
    const pol = document.getElementById("sel-pol").value;

    const cfgLines = cfg.method === "rv"
      ? [
          "cfg = waverider.WaveriderConfig( ...",
          "    'Method', 'rv', ...",
          `    'RVBodyType', '${cfg.rvBodyType}', ...`,
          `    'RVNoseType', '${cfg.rvNoseType}', ...`,
          `    'RVNoseBluntnessRatio', ${cfg.rvBluntness.toFixed(2)}, ...`,
          `    'Mach', ${cfg.mach.toFixed(2)}, ...`,
          `    'Length', ${cfg.length.toFixed(2)}, ...`,
          `    'RVConeAnglesDeg', [${cfg.rvConeAngles.map((v) => v.toFixed(2)).join(", ")}], ...`,
          `    'RVStageFractions', [${cfg.rvStageFractions.map((v) => v.toFixed(3)).join(", ")}], ...`,
          `    'RVAftType', '${cfg.rvAftType}', ...`,
          `    'RVAftCavityDepthRatio', ${cfg.rvAftDepth.toFixed(2)}, ...`,
          `    'RVAftCavityPower', ${cfg.rvAftPower.toFixed(2)}, ...`,
          `    'EnableAftFrustum', ${cfg.frustum}, ...`,
          `    'FrustumLength', ${cfg.flength.toFixed(2)}, ...`,
          `    'BoosterRadius', ${cfg.bradius.toFixed(2)}, ...`,
          `    'FrustumBlendLaw', '${cfg.blend}', ...`,
          `    'EnableRocketStage', ${cfg.rocket}, ...`,
          `    'RocketConeAngleDeg', ${cfg.rcone.toFixed(2)});`
        ]
      : [
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
          `    'RocketConeAngleDeg', ${cfg.rcone.toFixed(2)});`
        ];

    const stlOut = cfg.method === "rv" ? "output/custom_rv.stl" : "output/custom_waverider.stl";
    const code = [
      "setup_waverider();",
      "",
      ...cfgLines,
      "",
      "body = waverider.WaveriderBody(cfg);",
      "body.printSummary();",
      "body.plot();",
      `body.plotRCS('Polarization', '${pol}', 'VentralElevationDeg', ${elLook.toFixed(1)});`,
      `waverider.viz.plotRangeAngleIntensity(body, 'CenterFreqGHz', ${fcGHz.toFixed(1)}, 'BandwidthMHz', ${bwMHz.toFixed(0)}, 'ElevationDeg', ${elLook.toFixed(1)}, 'CutAzimuthDeg', ${rtiAz.toFixed(0)}, 'WindowType', '${winType}', 'Polarization', '${pol}');`,
      `body.exportSTL('${stlOut}');`
    ].join("\n");
    document.getElementById("matlab-code-output").textContent = code;
    dlgMatlab.showModal();
  });

  document.getElementById("btn-open-gallery").addEventListener("click", () => {
    dlgGallery.showModal();
  });

  // Preset Buttons (both Waverider and Reentry Vehicle presets)
  document.querySelectorAll("[data-preset]").forEach((btn) => {
    btn.addEventListener("click", () => {
      const p = PRESETS[btn.getAttribute("data-preset")];
      if (!p) return;
      if (p.arch === "rv") {
        hasInitializedRVTab = true;
        setCheckedRadio("rv-body-type", p.rvBodyType);
        setCheckedRadio("rv-nose-type", p.rvNoseType);
        setCheckedRadio("rv-aft-type", p.rvAftType);
        document.getElementById("rng-rv-mach").value = p.rvMach;
        document.getElementById("rng-rv-length").value = p.rvLength;
        document.getElementById("rng-rv-bluntness").value = p.rvBluntness;
        document.getElementById("rng-rv-cone1").value = p.rvCone1;
        document.getElementById("rng-rv-cone2").value = p.rvCone2;
        document.getElementById("rng-rv-cone3").value = p.rvCone3;
        document.getElementById("rng-rv-frac1").value = p.rvFrac1;
        document.getElementById("rng-rv-aftdepth").value = p.rvAftDepth;
        document.getElementById("rng-rv-aftpower").value = p.rvAftPower;
        document.getElementById("chk-frustum").checked = p.frustum;
        document.getElementById("chk-rocket").checked = p.rocket;
        const tabRV = document.getElementById("tab-rv");
        if (tabRV && activeArch !== "rv") {
          selectVehicleTab(tabRV, false);
        }
      } else {
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
        const tabWave = document.getElementById("tab-waverider");
        if (tabWave && activeArch !== "waverider") {
          selectVehicleTab(tabWave, false);
        }
      }
      enforcePhysicalSliderBounds();
      syncRVConditionalInertStates();
      rangeInputs.forEach(syncSliderVars);
      scheduleUpdate("geometry");
    });
  });

  // Bind each control to its exact minimum update tier (single event per control)
  const rcsControlIds = new Set(["sel-pol", "rng-rcsel", "rng-fc", "rng-bw", "sel-window", "chk-multibounce"]);
  const view3DControlIds = new Set(["chk-shock", "chk-rays", "chk-wireframe"]);

  document.querySelectorAll(".sidebar-controls input, .sidebar-controls select").forEach((el) => {
    const id = el.id;
    let level = "geometry";
    if (id === "rng-rtiaz") level = "cut";
    else if (rcsControlIds.has(id)) level = "rcs";
    else if (view3DControlIds.has(id)) level = "view3d";

    const evtName = el.type === "range" ? "input" : "change";
    el.addEventListener(evtName, () => scheduleUpdate(level));
  });

  // Interactive Azimuth Cut Selection on Canvas 5 (RTI Sinogram)
  let rtiDragging = false;
  function selectRTIAzimuthFromEvent(e) {
    if (!currentRTIData) return;
    const rect = canvasRTI.getBoundingClientRect();
    const padT = 34, padB = 48;
    const plotH = Math.max(rect.height - padT - padB, 100);
    const yRel = e.clientY - rect.top - padT;
    const normAz = 1.0 - Math.max(0.0, Math.min(1.0, yRel / plotH));
    const azSelected = Math.round((normAz * 180.0) / 3.0) * 3;
    const rngAz = document.getElementById("rng-rtiaz");
    rngAz.value = azSelected;
    syncSliderVars(rngAz);
    scheduleUpdate("cut");
  }
  canvasRTI.addEventListener("pointerdown", (e) => {
    rtiDragging = true;
    canvasRTI.setPointerCapture(e.pointerId);
    selectRTIAzimuthFromEvent(e);
  });
  canvasRTI.addEventListener("pointermove", (e) => {
    if (!rtiDragging) return;
    selectRTIAzimuthFromEvent(e);
  });
  canvasRTI.addEventListener("pointerup", () => {
    rtiDragging = false;
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
    scheduleUpdate("view3d");
  });
  canvas3D.addEventListener("pointerup", () => {
    camera.dragging = false;
  });
  canvas3D.addEventListener("wheel", (e) => {
    e.preventDefault();
    camera.zoom = Math.max(0.45, Math.min(2.6, camera.zoom * (e.deltaY < 0 ? 1.08 : 0.92)));
    scheduleUpdate("view3d");
  }, { passive: false });
  canvas3D.addEventListener("dblclick", () => {
    camera.yaw = -2.32;
    camera.pitch = 0.38;
    camera.zoom = 1.0;
    scheduleUpdate("view3d");
  });

  window.addEventListener("resize", () => {
    if (!currentModel) return;
    render3DViewport(currentModel);
    renderCrossSections(currentModel);
    if (currentRTIData) {
      updateRCSAndPlots();
    }
  });

  // Initialize synchronously on first paint
  updateGeometryAndAll();
})();
