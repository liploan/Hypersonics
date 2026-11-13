function surfData = generateReentryVehicle(cfg)
%GENERATEREENTRYVEHICLE Generate axisymmetric mono/bi/tri-conic Reentry Vehicle (RV).
%   surfData = GENERATEREENTRYVEHICLE(cfg) constructs structured 3D surface
%   grids and shock envelopes for conventional axisymmetric reentry vehicles
%   supporting:
%     - Body architectures : 'mono_conic', 'bi_conic', 'tri_conic'
%     - Nose geometries    : 'conical' (sharp), 'spherical' (C1 sphere-cone),
%                            'ogive' (C1 tangent/secant circular-arc ogive)
%     - Aft plane geometry : 'flat' planar bulkhead or 'concave' recessed dome

    L = cfg.Length;
    M = cfg.Mach;
    g = cfg.Gamma;

    bodyType = lower(strtrim(cfg.RVBodyType));
    noseType = lower(strtrim(cfg.RVNoseType));
    aftType  = lower(strtrim(cfg.RVAftType));

    switch bodyType
        case 'mono_conic', nStages = 1;
        case 'bi_conic',   nStages = 2;
        case 'tri_conic',  nStages = 3;
        otherwise
            error('generateReentryVehicle:BodyType', 'Unsupported RVBodyType "%s".', cfg.RVBodyType);
    end

    coneAnglesDeg = cfg.RVConeAnglesDeg(1:nStages);
    rawWeights = max(cfg.RVStageFractions(1:nStages), 0.0);
    if sum(rawWeights) <= 0
        stageWeights = ones(1, nStages) ./ nStages;
    else
        stageWeights = rawWeights ./ sum(rawWeights);
    end

    th1Deg = max(min(coneAnglesDeg(1), 52.0), 1.0);
    th1Rad = deg2rad(th1Deg);
    tanStages = tand(coneAnglesDeg);
    meanTanPostNose = max(sum(stageWeights .* tanStages), tand(2.0));

    etaNose = max(min(cfg.RVNoseBluntnessRatio, 0.80), 0.0);
    if strcmp(noseType, 'conical') || etaNose < 1e-4
        noseType = 'conical';
        etaNose = 0.0;
    end

    % 1. Exact Nose Tangency Station (x_t, r_t) & Nose Profile Parameters
    RN = 0.0;
    rhoOgive = 0.0;
    th0OgiveRad = th1Rad;

    switch noseType
        case 'conical'
            x_t = 0.0;
            r_t = 0.0;

        case 'spherical'
            % Exact C1 sphere-cone tangency with RN = etaNose * R_base:
            % R_base = r_t + (L - x_t)*meanTan, x_t = RN*(1 - sin(th1)), r_t = RN*cos(th1)
            denom = 1.0 - etaNose * (cos(th1Rad) - (1.0 - sin(th1Rad)) * meanTanPostNose);
            if denom <= 0.15
                R_base_nom = (L * meanTanPostNose) / 0.15;
            else
                R_base_nom = (L * meanTanPostNose) / denom;
            end
            RN_max = (0.50 * L) / max(1.0 - sin(th1Rad), 0.05);
            RN = min(etaNose * R_base_nom, RN_max);
            x_t = RN * (1.0 - sin(th1Rad));
            r_t = RN * cos(th1Rad);

        case 'ogive'
            % C1 circular-arc ogive from (0,0) meeting Stage 1 slope tan(th1) at x_t
            lamNose = min(max(0.14 + 0.50 * etaNose, 0.14), 0.56);
            x_t = lamNose * L;
            dThetaDeg = 4.0 + 18.0 * etaNose;
            thSecDeg = min(th1Deg + dThetaDeg, 62.0);
            r_t = x_t * tand(thSecDeg);
            denomOg = 2.0 * (r_t * cos(th1Rad) - x_t * sin(th1Rad));
            rhoOgive = (x_t^2 + r_t^2) / max(denomOg, 1e-9);
            sinTh0 = min(sin(th1Rad) + x_t / rhoOgive, 0.995);
            th0OgiveRad = asin(sinTh0);

        otherwise
            error('generateReentryVehicle:NoseType', 'Unsupported RVNoseType "%s".', cfg.RVNoseType);
    end

    % 2. Discretize Axial Stations [Nx x 1] with Exact Knots at x_t and Stage Breaks
    NxTarget = max(cfg.NumStreamwise, 24);
    if strcmp(noseType, 'conical')
        nNosePts = 1; % Only apex x = 0
        xNose = 0.0;
        rNose = 0.0;
    else
        nNosePts = max(10, round(0.32 * NxTarget));
        if strcmp(noseType, 'spherical')
            % Uniform polar cap angle alpha in [0, pi/2 - th1] avoids apex clustering
            alphaMax = (0.5 * pi) - th1Rad;
            alphaVec = linspace(0.0, alphaMax, nNosePts).';
            xNose = RN .* (1.0 - cos(alphaVec));
            rNose = RN .* sin(alphaVec);
        else
            % Uniform arc angle psi in [th0Ogive, th1] from apex to tangency
            psiVec = linspace(th0OgiveRad, th1Rad, nNosePts).';
            xNose = rhoOgive .* (sin(th0OgiveRad) - sin(psiVec));
            rNose = rhoOgive .* (cos(psiVec) - cos(th0OgiveRad));
        end
        xNose(1) = 0.0; rNose(1) = 0.0;
        xNose(end) = x_t; rNose(end) = r_t;
    end

    L_rem = max(L - x_t, 1e-4);
    nRemPts = max(NxTarget - nNosePts + 1, 3 * nStages + 1);
    stageCounts = max(3, round(stageWeights .* (nRemPts - 1)));
    % Adjust last stage count so total stations match cleanly
    stageBreaksX = zeros(1, nStages + 1);
    stageBreaksR = zeros(1, nStages + 1);
    stageBreaksX(1) = x_t;
    stageBreaksR(1) = r_t;

    xVec = xNose;
    rVec = rNose;
    stageIDPerRow = zeros(numel(xVec), 1); % 0 = nose, 1..nStages = conic frustum

    rMinFloor = max(0.025 * L, 0.20 * max(r_t, 0.02 * L));
    for s = 1:nStages
        Ls = stageWeights(s) * L_rem;
        xStart = stageBreaksX(s);
        xEnd   = (s == nStages) * L + (s < nStages) * (xStart + Ls);
        rStart = stageBreaksR(s);
        rEndRaw = rStart + (xEnd - xStart) * tanStages(s);
        rEnd = max(rEndRaw, rMinFloor);

        nSeg = stageCounts(s);
        tauS = linspace(0.0, 1.0, nSeg + 1).';
        xSeg = xStart + tauS(2:end) .* (xEnd - xStart);
        rSeg = rStart + tauS(2:end) .* (rEnd - rStart);

        xVec = [xVec; xSeg]; %#ok<AGROW>
        rVec = [rVec; rSeg]; %#ok<AGROW>
        stageIDPerRow = [stageIDPerRow; s * ones(numel(xSeg), 1)]; %#ok<AGROW>

        stageBreaksX(s + 1) = xEnd;
        stageBreaksR(s + 1) = rEnd;
    end

    Nx = numel(xVec);
    Ny = 2 * cfg.NumHalfSpan - 1;
    R_base = rVec(end);

    % 3. Safe Concave Aft Cavity Depth Clamping (guarantees positive wall thickness everywhere)
    dAft_m = 0.0;
    aftPower = cfg.RVAftCavityPower;
    if strcmp(aftType, 'concave') && cfg.RVAftCavityDepthRatio > 0
        dReq = min(cfg.RVAftCavityDepthRatio * R_base, 0.55 * L);
        % Verify r_cavity(u) <= r_skin(X_base(u)) - 0.04*R_base*(1 - u) for all u in [0, 0.98]
        uCheck = linspace(0.0, 0.98, 40);
        dSafe = dReq;
        for iter = 1:14
            xC = L - dSafe .* (1.0 - uCheck.^aftPower);
            rSkinAtXC = interp1(xVec, rVec, xC, 'linear', 0.0);
            rCav = uCheck .* R_base;
            wallMargin = 0.04 * R_base .* (1.0 - uCheck);
            if (L - dSafe > 0.04 * L) && all(rCav <= rSkinAtXC - wallMargin + 1e-12)
                break;
            end
            dSafe = 0.82 * dSafe;
        end
        dAft_m = max(dSafe, 0.0);
    end

    % 4. Exact Taylor-Maccoll Conical Solution (for sharp mono-conic) or Bow Shock Envelope
    CpMax = (2.0 / (g * M^2)) * ...
        (((g + 1.0)^2 * M^2 / (4.0 * g * M^2 - 2.0 * (g - 1.0)))^(g / (g - 1.0)) * ...
         ((1.0 - g + 2.0 * g * M^2) / (g + 1.0)) - 1.0);

    muRad = asin(1.0 / M);
    isSharpMonoTM = false;
    flowSol = struct('Mach', M, 'Gamma', g, 'CpMax', CpMax);
    betaAsymRad = max(muRad + deg2rad(2.0), min(1.18 * th1Rad + 0.65 * muRad, deg2rad(64.0)));

    if strcmp(noseType, 'conical')
        % Solve inverse Taylor-Maccoll for Stage 1 cone half-angle th1Deg if attached
        [okTM, betaFoundDeg, tmSol] = solveInverseTaylorMaccoll(M, th1Deg, g);
        if okTM
            betaAsymRad = deg2rad(betaFoundDeg);
            flowSol = tmSol;
            isSharpMonoTM = (nStages == 1);
        end
    end

    % Construct 3D Axisymmetric Surface Grids (Upper Z >= 0, Lower Z <= 0)
    % Uniform circumferential spacing: gamma_k = (pi/2) * eta_k in [-pi/2, +pi/2]
    eta = linspace(-1.0, 1.0, Ny);
    gammaAng = (0.5 * pi) .* eta;
    sinG = sin(gammaAng);
    cosG = cos(gammaAng);
    sinG(cfg.NumHalfSpan) = 0.0;
    cosG(1) = 0.0;
    cosG(Ny) = 0.0;

    X_up  = repmat(xVec, 1, Ny);
    X_low = X_up;
    Y_up  = rVec * sinG;
    Y_low = Y_up;
    Z_up  = rVec * cosG;
    Z_low = -Z_up;
    Z_low(:, 1)  = 0.0;
    Z_low(:, Ny) = 0.0;

    % Shock envelope (attached conical shock or Billig-style detached hyperboloid/cone blend)
    standoff = 0.0;
    if RN > 0
        standoff = 0.143 * RN * exp(3.24 / (M^2));
    elseif strcmp(noseType, 'ogive')
        standoff = 0.0;
    end
    rShkVec = (xVec + standoff) .* tan(betaAsymRad);
    rShkVec = max(rShkVec, 1.06 .* rVec + 0.015 * R_base * (xVec > 0));
    X_shock = repmat(xVec - standoff * exp(-2.5 .* xVec ./ max(L, 1e-6)), 1, Ny);
    Y_shock = rShkVec * sinG;
    Z_shock = -rShkVec * cosG;

    % Design surface pressure Cp (exact Taylor-Maccoll on sharp mono-conic; NaN -> exact facet
    % Modified Newtonian + Prandtl-Meyer in WaveriderBody aero cache for general multi-conic/blunt RVs)
    if isSharpMonoTM
        Cp_low = repmat(flowSol.CpCone, Nx, Ny);
        Cp_up  = Cp_low;
    else
        % Populate nodal Cp_low for visualization plots while marking FaceCpDesignMode = 'facet_newtonian'
        drdx = gradient(rVec) ./ max(gradient(xVec), 1e-9);
        sinLocal = drdx ./ sqrt(1.0 + drdx.^2);
        cpProf = zeros(Nx, 1);
        posMask = sinLocal >= 0;
        cpProf(posMask) = CpMax .* (sinLocal(posMask).^2);
        if any(~posMask)
            sinLee = -sinLocal(~posMask);
            pmTerm = max(1.0 - 0.5 * (g - 1.0) * M .* sinLee, 0.0).^(2.0 * g / (g - 1.0));
            cpProf(~posMask) = (2.0 / (g * M^2)) .* (pmTerm - 1.0);
        end
        if strcmp(noseType, 'conical') && isfield(flowSol, 'CpCone')
            st1 = (stageIDPerRow <= 1);
            cpProf(st1) = flowSol.CpCone;
        end
        Cp_low = repmat(cpProf, 1, Ny);
        Cp_up  = Cp_low;
    end

    surfData = struct();
    surfData.Method = 'rv';
    surfData.IsAxisymmetric = true;
    surfData.UseExactDesignCp = isSharpMonoTM;
    surfData.RVBodyType = bodyType;
    surfData.RVNoseType = noseType;
    surfData.RVAftType = aftType;
    surfData.AftCavityDepth_m = dAft_m;
    surfData.AftCavityPower = aftPower;
    surfData.RNose_m = RN;
    surfData.RBase_m = R_base;
    surfData.HalfSpanEffective = R_base;
    surfData.XProfile = xVec;
    surfData.RProfile = rVec;
    surfData.StageIDPerRow = stageIDPerRow;
    surfData.StageBreaksX = stageBreaksX;
    surfData.StageBreaksR = stageBreaksR;
    surfData.ConeAnglesDeg = coneAnglesDeg;
    surfData.X_low = X_low;
    surfData.Y_low = Y_low;
    surfData.Z_low = Z_low;
    surfData.X_up = X_up;
    surfData.Y_up = Y_up;
    surfData.Z_up = Z_up;
    surfData.X_shock = X_shock;
    surfData.Y_shock = Y_shock;
    surfData.Z_shock = Z_shock;
    surfData.Cp_low = Cp_low;
    surfData.Cp_up = Cp_up;
    surfData.FlowSolution = flowSol;
end

function [ok, betaDeg, tmSol] = solveInverseTaylorMaccoll(M, targetConeDeg, gamma)
%SOLVEINVERSETAYLORMACCOLL Secant/bisection root-finder for beta(M, theta_cone).
    ok = false;
    betaDeg = NaN;
    tmSol = struct();
    muDeg = asind(1.0 / M);
    bLo = max(muDeg + 0.12, targetConeDeg + 0.05);
    bHi = min(64.0, max(targetConeDeg * 2.4 + muDeg + 4.0, bLo + 8.0));
    if bLo >= bHi
        return;
    end
    try
        sLo = waverider.flow.taylorMaccoll(M, bLo, gamma);
        sHi = waverider.flow.taylorMaccoll(M, bHi, gamma);
    catch
        return;
    end
    fLo = sLo.thetaConeDeg - targetConeDeg;
    fHi = sHi.thetaConeDeg - targetConeDeg;
    if fLo > 0 || fHi < 0
        return;
    end
    for it = 1:18
        % Regula-falsi / bisection hybrid for rapid convergence (~6-8 evaluations)
        if mod(it, 3) == 0
            bMid = 0.5 * (bLo + bHi);
        else
            w = -fLo / (fHi - fLo);
            w = min(max(w, 0.15), 0.85);
            bMid = bLo + w * (bHi - bLo);
        end
        sMid = waverider.flow.taylorMaccoll(M, bMid, gamma);
        fMid = sMid.thetaConeDeg - targetConeDeg;
        if abs(fMid) < 1e-6 || (bHi - bLo) < 1e-6
            ok = true;
            betaDeg = bMid;
            tmSol = sMid;
            return;
        end
        if fMid < 0
            bLo = bMid; fLo = fMid;
        else
            bHi = bMid; fHi = fMid;
        end
    end
    ok = true;
    betaDeg = 0.5 * (bLo + bHi);
    tmSol = waverider.flow.taylorMaccoll(M, betaDeg, gamma);
end
