function results = test_physics_validation(varargin)
%TEST_PHYSICS_VALIDATION Independent physics, topology & signature validation suite.
%   results = TEST_PHYSICS_VALIDATION() checks the waverider codebase against
%   references that do NOT come from the code under test:
%
%     A. Flow solvers vs. an independent SciPy oblique-shock / Taylor-Maccoll
%        solution (and the NACA 1135 chart value M=2, 20 deg cone -> 37.80 deg).
%     B. Randomized geometry sweep across the full web-UI parameter ranges with
%        independent topology checks (directed half-edges, degenerate faces),
%        upper/lower surface crossover, and a generalized-winding-number test
%        that detects self-intersections / inverted regions (which the
%        freeBoundary/Euler checks in verifyWatertight cannot see).
%     C. Physical-optics RCS vs. closed-form canonical targets:
%          - square flat plate (normal + off-normal sinc pattern, coarse mesh)
%          - PEC sphere in the optical region (sigma -> pi a^2)
%          - 90 deg dihedral double-bounce (sigma = 16 pi a^2 b^2 sin^2(phi)/lambda^2)
%     D. Range profile vs. analytic point-spread function: scatterer position,
%        3-dB width (0.886 c/2B rectangular, 1.44 c/2B Hann) and peak sidelobe.
%
%   Name-Value options:
%     'NumRandom' (default 24)  number of randomized geometry configurations
%     'Seed'      (default 7)
%     'Strict'    (default true) error on first failure; false = report only

    p = inputParser;
    addParameter(p, 'NumRandom', 24);
    addParameter(p, 'Seed', 7);
    addParameter(p, 'Strict', true);
    parse(p, varargin{:});
    opts = p.Results;

    rootDir = fileparts(fileparts(mfilename('fullpath')));
    addpath(rootDir);

    results = struct('Name', {}, 'Pass', {}, 'Detail', {});
    function record(name, pass, detail)
        results(end + 1) = struct('Name', name, 'Pass', logical(pass), 'Detail', detail);
        tag = '[PASS]';
        if ~pass
            tag = '[FAIL]';
        end
        fprintf('%s %-58s %s\n', tag, name, detail);
        if ~pass && opts.Strict
            error('test_physics_validation:Failure', '%s: %s', name, detail);
        end
    end

    fprintf('\n=====================================================================\n');
    fprintf('  INDEPENDENT PHYSICS / TOPOLOGY / SIGNATURE VALIDATION\n');
    fprintf('=====================================================================\n');

    % =====================================================================
    % A. FLOW SOLVERS vs. INDEPENDENT SCIPY REFERENCE
    % =====================================================================
    %        M     beta   delta_wedge  M2_wedge   p2/p1    theta_cone  Cp_cone
    ref = [ 6.0, 15.00, 7.18494, 5.02407, 2.64680, 10.73711, 0.08164;
            7.0, 13.50, 6.94351, 5.76651, 2.94873,  9.95551, 0.06930;
            8.0, 13.00, 7.40949, 6.32931, 3.61169, 10.01540, 0.06867;
            5.0, 18.00, 8.48714, 4.14905, 2.61850, 12.79755, 0.11583;
           10.0, 10.00, 5.53325, 8.10685, 3.35126,  7.60744, 0.04003;
            2.0, 40.00, 10.62291, 1.61732, 1.76149, 22.43032, 0.39437;
            3.5, 24.00, 9.57164, 2.92968, 2.19767, 15.92470, 0.18411];
    maxErrDelta = 0; maxErrCone = 0; maxErrCp = 0; maxErrM2 = 0;
    for r = 1:size(ref, 1)
        sh = waverider.flow.obliqueShock(ref(r, 1), ref(r, 2));
        tm = waverider.flow.taylorMaccoll(ref(r, 1), ref(r, 2));
        maxErrDelta = max(maxErrDelta, abs(sh.deltaDeg - ref(r, 3)));
        maxErrM2    = max(maxErrM2, abs(sh.M2 - ref(r, 4)));
        maxErrCone  = max(maxErrCone, abs(tm.thetaConeDeg - ref(r, 6)));
        maxErrCp    = max(maxErrCp, abs(tm.Cp(1) - ref(r, 7)));
    end
    record('A1 Oblique shock delta & M2 vs SciPy', maxErrDelta < 1e-4 && maxErrM2 < 1e-4, ...
        sprintf('max|d delta| = %.2e deg, max|d M2| = %.2e', maxErrDelta, maxErrM2));
    record('A2 Taylor-Maccoll cone angle vs SciPy', maxErrCone < 2e-3, ...
        sprintf('max|d theta_c| = %.2e deg', maxErrCone));
    record('A3 Taylor-Maccoll cone-surface Cp vs SciPy', maxErrCp < 2e-4, ...
        sprintf('max|d Cp_c| = %.2e', maxErrCp));
    % NACA 1135 chart: M = 2, 20 deg cone -> beta = 37.80 deg
    betaNACA = fzero(@(b) waverider.flow.taylorMaccoll(2.0, b).thetaConeDeg - 20.0, [32, 50]);
    record('A4 NACA-1135 chart: M=2, 20deg cone -> beta=37.80', abs(betaNACA - 37.80) < 0.05, ...
        sprintf('beta = %.3f deg', betaNACA));

    % =====================================================================
    % B. GEOMETRY: PRESETS + RANDOMIZED SWEEP OVER THE WEB-UI RANGES
    % =====================================================================
    rng(opts.Seed);
    methods = {'conical', 'osculating_cone', 'caret', 'cst'};
    blends = {'cubic', 'quintic', 'linear', 'conic'};
    cfgList = {};
    names = {};
    for pName = {'conical', 'osculating_cone', 'caret', 'cst', 'booster_mounted'}
        cfgList{end + 1} = waverider.WaveriderConfig(pName{1}); %#ok<AGROW>
        names{end + 1} = ['preset:' pName{1}]; %#ok<AGROW>
    end
    for n = 1:opts.NumRandom
        c = waverider.WaveriderConfig();
        c.Method = methods{randi(4)};
        c.Mach = 3.5 + rand * 11.5;
        muDeg = asind(1 / c.Mach);
        c.ShockAngleDeg = max(10 + rand * 14, muDeg + 2.0);
        c.Length = 2.5 + rand * 7.5;
        c.HalfSpan = 0.8 + rand * 2.7;
        c.PlanformPower = 0.35 + rand * 0.90;
        c.LeadingEdgeAnhedralDeg = rand * 20;
        c.CanopyHeightRatio = rand * 0.12;
        c.UpperSurfaceType = 'canopy';
        if c.CanopyHeightRatio < 0.01
            c.UpperSurfaceType = 'freestream';
        end
        c.EnableAftFrustum = rand > 0.4;
        c.FrustumLength = 0.5 + rand * 3.0;
        c.BoosterRadius = 0.30 + rand * 1.10;
        c.FrustumBlendLaw = blends{randi(4)};
        c.EnableRocketStage = c.EnableAftFrustum && (rand > 0.4);
        c.RocketConeAngleDeg = rand * 12;
        cfgList{end + 1} = c; %#ok<AGROW>
        names{end + 1} = sprintf('random#%02d:%s', n, c.Method); %#ok<AGROW>
    end
    % Deterministic edge cases at the limits of the WaveriderConfig validators
    edgeDefs = {
        'conical',         {'PlanformPower', 0.25};
        'osculating_cone', {'PlanformPower', 0.25};
        'caret',           {'PlanformPower', 0.25};
        'cst',             {'PlanformPower', 0.25};
        'conical',         {'PlanformPower', 2.5};
        'osculating_cone', {'PlanformPower', 2.5, 'LeadingEdgeAnhedralDeg', 35};
        'caret',           {'LeadingEdgeAnhedralDeg', 35};
        'caret',           {'LeadingEdgeAnhedralDeg', -35};
        'osculating_cone', {'Mach', 3.5, 'ShockAngleDeg', 18.5};
        'booster_mounted', {'BoosterRadius', 0.3, 'FrustumLength', 0.5};
        'booster_mounted', {'BoosterRadius', 1.4, 'FrustumLength', 0.5};
        'booster_mounted', {'RocketConeAngleDeg', -15};
        'cst',             {'CanopyHeightRatio', 0.35}};
    for e = 1:size(edgeDefs, 1)
        c = waverider.WaveriderConfig(edgeDefs{e, 1});
        nv = edgeDefs{e, 2};
        for j = 1:2:numel(nv)
            c.(nv{j}) = nv{j + 1};
        end
        cfgList{end + 1} = c; %#ok<AGROW>
        names{end + 1} = sprintf('edge:%s{%s}', edgeDefs{e, 1}, strjoin(cellfun(@num2str, nv, 'UniformOutput', false), ',')); %#ok<AGROW>
    end

    nCfg = numel(cfgList);
    geoFail = {};
    worstWinding = 0; worstCross = 0; nDegenTotal = 0; worstVolErr = 0;
    for n = 1:nCfg
        c = cfgList{n};
        try
            body = waverider.WaveriderBody(c);
        catch ME
            geoFail{end + 1} = sprintf('%s: generation error "%s"', names{n}, ME.message); %#ok<AGROW>
            continue;
        end
        V = body.MeshData.TR.Points;
        F = body.MeshData.TR.ConnectivityList;
        chk = independentTopology(V, F);
        sd = body.SurfaceData;
        nxw = sd.NxWaverider;
        cross = max(max(sd.Z_low(1:nxw, :) - sd.Z_up(1:nxw, :)));
        worstCross = max(worstCross, cross);
        nDegenTotal = nDegenTotal + chk.nDegenerate;
        % Spanwise fold check: Y must be monotone across every streamwise row
        dYl = diff(sd.Y_low, 1, 2); dYu = diff(sd.Y_up, 1, 2);
        nFold = sum(any(dYl < -1e-9, 2) & any(dYl > 1e-9, 2)) + sum(any(dYu < -1e-9, 2) & any(dYu > 1e-9, 2));

        % Generalized winding number on random bounding-box points
        nPts = 600;
        if startsWith(names{n}, 'preset')
            nPts = 2500;
        end
        bbMin = min(V, [], 1); bbMax = max(V, [], 1);
        P = bbMin + rand(nPts, 3) .* (bbMax - bbMin);
        W = windingNumber(V, F, P);
        dev = min(abs(W), abs(W - 1));
        worstWinding = max(worstWinding, max(dev));
        insideFrac = mean(abs(W - 1) < 0.25);
        volMC = insideFrac * prod(bbMax - bbMin);
        % binomial standard error of the Monte-Carlo estimate
        volSigma = sqrt(max(insideFrac * (1 - insideFrac), 1e-6) / nPts) * prod(bbMax - bbMin);
        volErr = abs(volMC - chk.volume) / max(volSigma, 1e-9);
        if startsWith(names{n}, 'preset')
            worstVolErr = max(worstVolErr, volErr);
        end

        problems = {};
        if ~chk.manifoldOriented, problems{end + 1} = sprintf('%d bad directed edges', chk.nBadDirected); end %#ok<AGROW>
        if chk.volume <= 0, problems{end + 1} = 'non-positive volume'; end %#ok<AGROW>
        if max(dev) > 0.05, problems{end + 1} = sprintf('winding number %.2f (self-intersection/inversion)', W(find(dev == max(dev), 1))); end %#ok<AGROW>
        if cross > 1e-6, problems{end + 1} = sprintf('lower surface above upper by %.3g m', cross); end %#ok<AGROW>
        if chk.nDegenerate > 0, problems{end + 1} = sprintf('%d degenerate faces (min area %.2e)', chk.nDegenerate, chk.minArea); end %#ok<AGROW>
        if nFold > 0, problems{end + 1} = sprintf('%d streamwise rows fold spanwise (non-monotone Y)', nFold); end %#ok<AGROW>
        if ~body.Report.IsWatertight, problems{end + 1} = 'verifyWatertight reports FAIL'; end %#ok<AGROW>
        if ~isempty(problems)
            geoFail{end + 1} = sprintf('%s: %s', names{n}, strjoin(problems, '; ')); %#ok<AGROW>
        end
    end
    record('B1 Geometry sweep: all configs generate & close', isempty(geoFail), ...
        sprintf('%d configs, %d with problems', nCfg, numel(geoFail)));
    for k = 1:numel(geoFail)
        fprintf('        -> %s\n', geoFail{k});
    end
    record('B2 No self-intersection (winding number in {0,1})', worstWinding < 0.05, ...
        sprintf('worst |W - {0,1}| = %.3f', worstWinding));
    record('B3 No upper/lower surface crossover', worstCross < 1e-6, ...
        sprintf('worst Z_low - Z_up = %.3g m', worstCross));
    record('B4 No degenerate (zero-area) triangles', nDegenTotal == 0, ...
        sprintf('%d degenerate faces across sweep', nDegenTotal));
    record('B5 Monte-Carlo volume == divergence volume (presets)', worstVolErr < 4.0, ...
        sprintf('worst deviation = %.2f sigma (binomial)', worstVolErr));

    % =====================================================================
    % C. PHYSICAL-OPTICS RCS vs. CLOSED-FORM CANONICAL TARGETS
    % =====================================================================
    c0 = 299792458.0;

    % C1/C2: Square PEC plate a x a in the X-Y plane, normal +Z, coarse 2-triangle mesh
    a = 1.0;
    Vp = [0 0 0; a 0 0; a a 0; 0 a 0] - [a/2, a/2, 0];
    Fp = [1 2 3; 1 3 4];
    plate = fakeBody(Vp, Fp, 1);
    fGHz = 10.0; lam = c0 / (fGHz * 1e9); k = 2 * pi / lam;
    r0 = waverider.rcs.computeRCS(plate, fGHz, 0, 90, 'IncludeMultiBounce', false, 'IncludeEdgeDiffraction', false);
    sig0 = 4 * pi * a^4 / lam^2;
    err0 = r0.SigmaTotal_dBsm - 10 * log10(sig0);
    record('C1 Flat plate normal incidence: 4 pi A^2 / lambda^2', abs(err0) < 0.1, ...
        sprintf('model %.2f dBsm vs exact %.2f dBsm (err %+.2f dB)', r0.SigmaTotal_dBsm, 10*log10(sig0), err0));

    % Off-normal: rotate look direction in X-Z plane (az = 0 => u = [-cos(el), 0, sin(el)])
    thList = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0];
    errs = zeros(size(thList));
    for t = 1:numel(thList)
        th = deg2rad(thList(t));
        el = 90 - thList(t);
        rr = waverider.rcs.computeRCS(plate, fGHz, 0, el, 'IncludeMultiBounce', false, 'IncludeEdgeDiffraction', false);
        x = k * a * sin(th);
        exact = sig0 * cos(th)^2 * (sin(x) / x)^2;
        errs(t) = rr.SigmaTotal_dBsm - 10 * log10(max(exact, 1e-5));
    end
    % compare only where exact pattern is above -45 dB relative (avoid deep nulls)
    record('C2 Flat plate off-normal sinc^2 pattern (2 coarse facets, 33 lambda)', max(abs(errs)) < 1.0, ...
        sprintf('errors [dB] = %s', mat2str(round(errs, 2))));

    % C3: PEC sphere radius 1 m (PO optical limit sigma = pi a^2)
    [Vs, Fs] = icosphere(5);  % 20480 faces
    sphere = fakeBody(Vs, Fs, 1);
    sphErr = zeros(1, 2);
    fS = [1.0, 3.0];
    for t = 1:2
        rs = waverider.rcs.computeRCS(sphere, fS(t), [0; 37; 90], [0; 21; -50], 'IncludeMultiBounce', false, 'IncludeEdgeDiffraction', false);
        sphErr(t) = max(abs(rs.SigmaTotal_dBsm - 10 * log10(pi)));
    end
    record('C3 PEC sphere a=1 m (1, 3 GHz) -> pi a^2 = 4.97 dBsm', max(sphErr) < 1.0, ...
        sprintf('max error = %.2f dB (1 GHz), %.2f dB (3 GHz)', sphErr(1), sphErr(2)));

    % C4/C5: 90-degree dihedral trough opening downward (concave "underside")
    aD = 0.5; bD = 1.0; fD = 10.0; lamD = c0 / (fD * 1e9);
    dih = dihedralBody(aD, bD, 40, 31);
    % Symmetric incidence from directly below (el = -90): phi = 45 deg
    rd = waverider.rcs.computeRCS(dih, fD, 90, -90, 'IncludeEdgeDiffraction', false);
    sigD = 16 * pi * aD^2 * bD^2 * sin(pi/4)^2 / lamD^2;
    errD = rd.SigmaTotal_dBsm - 10 * log10(sigD);
    record('C4 90deg dihedral double-bounce, symmetric (8 pi a^2 b^2/lam^2)', abs(errD) < 1.5, ...
        sprintf('model %.2f dBsm vs exact %.2f dBsm (err %+.2f dB)', rd.SigmaTotal_dBsm, 10*log10(sigD), errD));
    % 15 deg off symmetric in the plane normal to the fold: phi = 30 deg
    rd2 = waverider.rcs.computeRCS(dih, fD, 90, -75, 'IncludeEdgeDiffraction', false);
    sigD2 = 16 * pi * aD^2 * bD^2 * sin(deg2rad(30))^2 / lamD^2;
    errD2 = rd2.SigmaTotal_dBsm - 10 * log10(sigD2);
    record('C5 90deg dihedral double-bounce, 15deg off-symmetric (phi=30)', abs(errD2) < 1.5, ...
        sprintf('model %.2f dBsm vs exact %.2f dBsm (err %+.2f dB)', rd2.SigmaTotal_dBsm, 10*log10(sigD2), errD2));

    % =====================================================================
    % D. RANGE PROFILE vs. ANALYTIC POINT-SPREAD FUNCTION
    % =====================================================================
    % Two small plates facing the radar (az = 0 => u = [-1,0,0]) separated by 2.0 m in X
    s = 0.04;
    sq = [0 -s -s; 0 s -s; 0 s s; 0 -s s];
    Vr = [sq + [2.0 0 0]; sq + [4.0 0 0]];
    Fr = [1 3 2; 1 4 3; 5 7 6; 5 8 7];   % normals toward -X
    pair = fakeBody(Vr, Fr, 1);
    pair.Config = struct('Length', 4.0);
    for win = {'rect', 'hann'}
        B = 1000;
        rp = waverider.rcs.computeRangeProfile(pair, 10.0, B, 0, 0, 'WindowType', win{1}, ...
            'IncludeMultiBounce', false, 'NumRangeBins', 4001);
        R = rp.rangeVec_m; prof = rp.RTI_total_m2(1, :);
        xRef = rp.rRef(1);
        % expected down-range of each plate: R = -(x - xRef) * (-1) = x - xRef
        Rexp = [2.0, 4.0] - xRef;
        [pk1, i1] = max(prof .* (abs(R - Rexp(1)) < 0.5));
        [~, i2] = max(prof .* (abs(R - Rexp(2)) < 0.5));
        posErr = max(abs([R(i1), R(i2)] - Rexp));
        % 3-dB width of first peak
        half = pk1 / 2;
        lo = find(prof(1:i1) < half, 1, 'last');
        hi = i1 - 1 + find(prof(i1:end) < half, 1, 'first');
        w3 = R(hi) - R(lo);
        dR = c0 / (2 * B * 1e6);
        if strcmp(win{1}, 'rect')
            wExp = 0.886 * dR; slExp = -13.26;
        else
            wExp = 1.44 * dR; slExp = -31.5;
        end
        % peak sidelobe within +/- 6 dR of peak 1, outside mainlobe (first nulls)
        winMask = abs(R - R(i1)) < 6 * dR & abs(R - R(i1)) > (1.0 + strcmp(win{1}, 'hann')) * dR;
        sl = 10 * log10(max(prof(winMask)) / pk1);
        record(sprintf('D1 Range profile (%s): scatterer positions', win{1}), posErr < 0.01, ...
            sprintf('max position error = %.4f m', posErr));
        record(sprintf('D2 Range profile (%s): 3-dB width', win{1}), abs(w3 - wExp) / wExp < 0.08, ...
            sprintf('%.2f cm vs analytic %.2f cm', 100 * w3, 100 * wExp));
        record(sprintf('D3 Range profile (%s): peak sidelobe', win{1}), abs(sl - slExp) < 2.0, ...
            sprintf('%.1f dB vs analytic %.1f dB', sl, slExp));
    end

    % =====================================================================
    % E. EXACT AERODYNAMICS & BILATERAL SYMMETRY INVARIANCE
    % =====================================================================
    % E1: Nonweiler exact inviscid caret theorem: CL = Cp_wedge, CD = Cp_wedge*tan(delta), L/D = cot(delta)
    cfgCaret = waverider.WaveriderConfig('Method', 'caret', 'Mach', 6.5, 'ShockAngleDeg', 15.0, ...
        'PlanformPower', 1.0, 'UpperSurfaceType', 'freestream', 'CanopyHeightRatio', 0.0);
    bodyCaret = waverider.WaveriderBody(cfgCaret);
    shCaret = waverider.flow.obliqueShock(6.5, 15.0, 1.4);
    cpWedgeExact = (shCaret.p2_p1 - 1.0) / (0.5 * 1.4 * 6.5^2);
    clExact = cpWedgeExact;
    cdExact = cpWedgeExact * tan(shCaret.deltaRad);
    ldExact = 1.0 / tan(shCaret.deltaRad);
    errCaretAero = max([abs(bodyCaret.AeroMetrics.CL - clExact), ...
                        abs(bodyCaret.AeroMetrics.CD_wave - cdExact), ...
                        abs(bodyCaret.AeroMetrics.LD_inviscid - ldExact)]);
    record('E1 Nonweiler caret theorem: CL=Cp2, L/D=cot(delta)', errCaretAero < 1e-11, ...
        sprintf('L/D = %.6f vs cot(delta) = %.6f (max err %.2e)', bodyCaret.AeroMetrics.LD_inviscid, ldExact, errCaretAero));

    % E2: AeroMetrics == computeAeroPolar(0) consistency across all presets
    maxPolarDiff = 0.0;
    maxCY = 0.0;
    for pName = {'conical', 'osculating_cone', 'caret', 'cst', 'booster_mounted'}
        bP = waverider.WaveriderBody(pName{1});
        pol0 = bP.computeAeroPolar(0.0);
        maxPolarDiff = max([maxPolarDiff, abs(pol0.CL - bP.AeroMetrics.CL), abs(pol0.CD_wave - bP.AeroMetrics.CD_wave)]);
        maxCY = max(maxCY, abs(bP.AeroMetrics.CY));
    end
    record('E2 AeroMetrics == computeAeroPolar(0) identity', maxPolarDiff < 1e-13, ...
        sprintf('max|polar(0) - AeroMetrics| = %.2e', maxPolarDiff));

    % E3: Bilateral port/starboard symmetry (CY == 0 and RCS(+az) == RCS(-az))
    bSym = waverider.WaveriderBody('osculating_cone');
    rcsSym = bSym.computeRCS(10.0, [-55; -25; 25; 55], -20.0);
    rcsAsymErr = max(abs(rcsSym.SigmaTotal_dBsm([1, 2]) - rcsSym.SigmaTotal_dBsm([4, 3])));
    record('E3 Bilateral mesh symmetry: |CY|=0 & RCS(+az)==RCS(-az)', maxCY < 1e-13 && rcsAsymErr < 1e-11, ...
        sprintf('max|CY| = %.2e, max|dRCS(+/-az)| = %.2e dB', maxCY, rcsAsymErr));

    nPass = sum([results.Pass]);
    fprintf('=====================================================================\n');
    fprintf('  %d / %d VALIDATION CHECKS PASSED\n', nPass, numel(results));
    fprintf('=====================================================================\n');
end

% -------------------------------------------------------------------------
function chk = independentTopology(V, F)
    nF = size(F, 1);
    he = [F(:, [1 2]); F(:, [2 3]); F(:, [3 1])];
    nV = size(V, 1);
    keyD = (he(:, 1) - 1) * nV + he(:, 2);             % directed
    keyR = (he(:, 2) - 1) * nV + he(:, 1);             % reversed
    [uD, ~, ic] = unique(keyD);
    cntD = accumarray(ic, 1);
    dupDirected = sum(cntD > 1);
    hasTwin = ismember(keyD, keyR);
    chk.nBadDirected = dupDirected + sum(~hasTwin);
    chk.manifoldOriented = (chk.nBadDirected == 0) && numel(uD) == 3 * nF;
    cr = cross(V(F(:, 2), :) - V(F(:, 1), :), V(F(:, 3), :) - V(F(:, 1), :), 2);
    area = 0.5 * sqrt(sum(cr.^2, 2));
    chk.minArea = min(area);
    edgeL = max([vecnorm(V(F(:,2),:)-V(F(:,1),:),2,2), vecnorm(V(F(:,3),:)-V(F(:,2),:),2,2), vecnorm(V(F(:,1),:)-V(F(:,3),:),2,2)], [], 2);
    chk.nDegenerate = sum(area < 1e-9 * max(edgeL.^2, 1e-12) | area < 1e-14);
    chk.volume = sum(dot(V(F(:, 1), :), cross(V(F(:, 2), :), V(F(:, 3), :), 2), 2)) / 6;
end

function W = windingNumber(V, F, P)
    % Generalized winding number via van Oosterom-Strackee solid angles
    nP = size(P, 1);
    W = zeros(nP, 1);
    A = V(F(:, 1), :); B = V(F(:, 2), :); C = V(F(:, 3), :);
    chunk = 40;
    for s = 1:chunk:nP
        idx = s:min(s + chunk - 1, nP);
        for j = idx
            a = A - P(j, :); b = B - P(j, :); c = C - P(j, :);
            la = vecnorm(a, 2, 2); lb = vecnorm(b, 2, 2); lc = vecnorm(c, 2, 2);
            num = dot(a, cross(b, c, 2), 2);
            den = la .* lb .* lc + dot(a, b, 2) .* lc + dot(a, c, 2) .* lb + dot(b, c, 2) .* la;
            W(j) = sum(2 * atan2(num, den)) / (4 * pi);
        end
    end
end

function body = fakeBody(V, F, patchID)
    body = struct();
    body.MeshData.TR = triangulation(F, V);
    body.MeshData.FacePatchID = patchID * ones(size(F, 1), 1);
    body.SurfaceData.X_low = V(1:2, 1).';
    body.SurfaceData.Y_low = V(1:2, 2).';
    body.SurfaceData.Z_low = V(1:2, 3).';
    body.SurfaceData.NxWaverider = 1;
    body.Config = struct('Length', max(V(:, 1)));
end

function body = dihedralBody(a, b, nX, nY)
    % Trough along X, opening downward: Z = -|Y| (faces at +/-45 deg), face width a
    yMax = a / sqrt(2);
    x = linspace(0, b, nX);
    y = linspace(-yMax, yMax, nY);
    [Y, X] = meshgrid(y, x);
    Z = -abs(Y);
    V = [X(:), Y(:), Z(:)];
    id = reshape(1:numel(X), nX, nY);
    F = zeros(2 * (nX - 1) * (nY - 1), 3); m = 0;
    for i = 1:nX - 1
        for j = 1:nY - 1
            p1 = id(i, j); p2 = id(i + 1, j); p3 = id(i + 1, j + 1); p4 = id(i, j + 1);
            m = m + 1; F(m, :) = [p1 p3 p2];   % outward normal downward
            m = m + 1; F(m, :) = [p1 p4 p3];
        end
    end
    body = struct();
    body.MeshData.TR = triangulation(F, V);
    body.MeshData.FacePatchID = 2 * ones(size(F, 1), 1);
    body.SurfaceData.X_low = X; body.SurfaceData.Y_low = Y; body.SurfaceData.Z_low = Z;
    body.SurfaceData.NxWaverider = nX;
    body.Config = struct('Length', b);
    % sanity: normals point downward
    n = faceNormal(body.MeshData.TR);
    assert(all(n(:, 3) < 0), 'dihedral normals must point downward');
end

function [V, F] = icosphere(nSub)
    t = (1 + sqrt(5)) / 2;
    V = [-1 t 0; 1 t 0; -1 -t 0; 1 -t 0; 0 -1 t; 0 1 t; 0 -1 -t; 0 1 -t; t 0 -1; t 0 1; -t 0 -1; -t 0 1];
    V = V ./ vecnorm(V, 2, 2);
    F = [1 12 6; 1 6 2; 1 2 8; 1 8 11; 1 11 12; 2 6 10; 6 12 5; 12 11 3; 11 8 7; 8 2 9;
         4 10 5; 4 5 3; 4 3 7; 4 7 9; 4 9 10; 5 10 6; 3 5 12; 7 3 11; 9 7 8; 10 9 2];
    for s = 1:nSub
        E = sort([F(:, [1 2]); F(:, [2 3]); F(:, [3 1])], 2);
        [Eu, ~, ic] = unique(E, 'rows');
        mid = (V(Eu(:, 1), :) + V(Eu(:, 2), :)) / 2;
        mid = mid ./ vecnorm(mid, 2, 2);
        nV = size(V, 1);
        V = [V; mid]; %#ok<AGROW>
        nF = size(F, 1);
        m12 = nV + ic(1:nF); m23 = nV + ic(nF + 1:2 * nF); m31 = nV + ic(2 * nF + 1:3 * nF);
        F = [F(:, 1) m12 m31; F(:, 2) m23 m12; F(:, 3) m31 m23; m12 m23 m31];
    end
end
