function surfData = generateOsculatingCone(cfg)
%GENERATEOSCULATINGCONE Cone-derived and osculating-cones waverider surfaces.
%   surfData = GENERATEOSCULATINGCONE(cfg) builds the lower compression
%   surface as a stream surface of an exact supersonic flowfield whose shock
%   passes through a prescribed leading edge.
%
%   'conical' (Jones et al. 1968; Rasmussen 1980):
%     Single axisymmetric Taylor-Maccoll cone flow. The leading edge is the
%     intersection of the planform x_LE(y) with the conical shock; each lower-
%     surface streamline is traced in its meridional plane.
%
%   'osculating_cone' (Sobieczky, Dougherty & Jones 1990 inverse method):
%     1. Inlet capture curve (ICC = shock trace in the base plane x = L):
%          z_s(y) = -L tan(beta)                           |y| <= y_f
%          z_s(y) = -L tan(beta) + h ((|y|-y_f)/(b-y_f))^n  |y| >  y_f
%        y_f = FlatShockFraction*b, n = ShockCurvaturePower,
%        h = (b - y_f) * tan(Gamma_s), Gamma_s = max(LeadingEdgeAnhedralDeg, 0)
%        (outboard shock rise angle; the shock must be convex toward the body).
%     2. At each ICC point: local radius of curvature R = (1+z'^2)^1.5/|z''|,
%        osculating plane spanned by x and the curve normal, local cone with
%        axis through the centre of curvature and apex at x = L - R/tan(beta).
%     3. Leading edge from the planform x_LE(eta), eta = y_s/b, placed on the
%        local shock: radial distance rho_LE = R - (L - x_LE) tan(beta).
%     4. Lower surface: Taylor-Maccoll streamline of the local cone from the
%        leading edge to the base; where the ICC is flat (R -> inf) the exact
%        2-D oblique-shock wedge solution (straight streamlines at delta).
%     5. Upper surface: freestream (or canopy) surface from the leading edge.
%   Osculating planes only intersect at the centres of curvature, so the
%   surfaces cannot fold as long as rho > 0; the planform is locally limited
%   to rho_LE >= 0.05 R where needed (reported in surfData.PlanformClampCount).
%
%   Pressure on the lower surface (Cp_low) is the exact Taylor-Maccoll (or
%   wedge) value at each node; it is the design-point inviscid solution.

    L = cfg.Length;
    bHalf = cfg.HalfSpan;
    betaDeg = cfg.ShockAngleDeg;
    betaRad = deg2rad(betaDeg);
    tanBeta = tan(betaRad);

    tm = waverider.flow.taylorMaccoll(cfg.Mach, betaDeg, cfg.Gamma);
    shock = tm.shock;
    tanDelta = tan(shock.deltaRad);
    CpWedge = (shock.p2_p1 - 1.0) / (0.5 * cfg.Gamma * cfg.Mach^2);

    Ny = 2 * cfg.NumHalfSpan - 1;
    Nx = cfg.NumStreamwise;
    iC = cfg.NumHalfSpan;
    eta = linspace(-1.0, 1.0, Ny);

    invP = 1.0 / cfg.PlanformPower;
    xLE = L .* (abs(eta).^invP);
    xLE([1, Ny]) = L;
    xLE(iC) = 0.0;
    xi = linspace(0.0, 1.0, Nx).';
    clampCount = 0;

    if strcmpi(cfg.Method, 'conical')
        % -------------------------------------------------------------
        % Axisymmetric cone-derived waverider
        % -------------------------------------------------------------
        yLE = bHalf .* eta;
        minR_L = 1.28 * bHalf;
        % Ensure shock cone encloses span AND d(phi_LE)/d|eta| > 0 for all P < 1
        x0_monotone = max((invP - 1.0 + 0.08) * L, 0.0);
        x0_min = max([minR_L / tanBeta - L, 0.15 * L, x0_monotone]);
        x0_cfg = (cfg.VertexOffsetRatio / (1.0 - cfg.VertexOffsetRatio)) * L;
        x0 = max(x0_cfg, x0_min);
        if x0 > x0_cfg + 1e-9
            warning('waverider:VertexOffsetRaised', ...
                'VertexOffsetRatio raised so shock cone encloses span without meridional fold (x0 %.3f -> %.3f m).', x0_cfg, x0);
        end
        zv = x0 * tanBeta;                 % cone axis height: shock passes through the nose

        xStartLocal = xLE + x0;
        xEndLocal = repmat(L + x0, 1, Ny);
        R_shock_LE = xStartLocal .* tanBeta;
        zDepthLE = sqrt(max(R_shock_LE.^2 - yLE.^2, 1e-12));
        phiLE = atan2(yLE, zDepthLE);
        phiLE(iC) = 0.0;

        [xLocalGrid, rGrid, thetaGrid] = waverider.flow.traceStreamline(tm, xStartLocal, xEndLocal, Nx);
        X_low = xLocalGrid - x0;
        Y_low = rGrid .* sin(repmat(phiLE, Nx, 1));
        Z_low = zv - rGrid .* cos(repmat(phiLE, Nx, 1));
        Y_low(:, iC) = 0.0;
        Cp_low = interp1(tm.thetaRad, tm.Cp, thetaGrid, 'pchip', tm.Cp(end));

        R_shock_grid = xLocalGrid .* tanBeta;
        X_shock = X_low;
        Y_shock = R_shock_grid .* sin(repmat(phiLE, Nx, 1));
        Z_shock = zv - R_shock_grid .* cos(repmat(phiLE, Nx, 1));
    else
        % -------------------------------------------------------------
        % Osculating cones (inverse design from a prescribed shock curve)
        % -------------------------------------------------------------
        if cfg.LeadingEdgeAnhedralDeg < 0
            warning('waverider:ConcaveICC', ...
                ['Negative LeadingEdgeAnhedralDeg would make the shock curve concave toward ' ...
                 'the body (no external conical flow exists); using 0 (flat shock).']);
        end
        gammaS = deg2rad(max(cfg.LeadingEdgeAnhedralDeg, 0.0));
        yF = cfg.FlatShockFraction * bHalf;
        nExp = cfg.ShockCurvaturePower;
        wOut = max(bHalf - yF, 1e-9);
        hRise = wOut * tan(gammaS);
        zC = -L * tanBeta;

        yS = bHalf .* eta;                          % ICC parameter (shock-curve y)
        s = max(abs(yS) - yF, 0.0) ./ wOut;
        if nExp >= 2.0
            gVal = s.^nExp;
            dgds = nExp .* (s.^(nExp - 1.0));
            d2gds2 = nExp * (nExp - 1.0) .* (s.^(nExp - 2.0));
        else
            % Smooth C^2 regularization when 1 < nExp < 2 so curvature stays finite at s -> 0+
            s0 = 0.18;
            denomG = (1.0 + s0^2)^(0.5 * nExp) - s0^nExp;
            u2 = s.^2 + s0^2;
            gVal = (u2.^(0.5 * nExp) - s0^nExp) ./ denomG;
            dgds = (nExp .* s .* (u2.^(0.5 * nExp - 1.0))) ./ denomG;
            d2gds2 = (nExp .* (u2.^(0.5 * nExp - 2.0)) .* ((nExp - 1.0) .* (s.^2) + s0^2)) ./ denomG;
        end
        inFlatSeg = (abs(yS) <= yF + 1e-12);
        gVal(inFlatSeg) = 0.0;
        dgds(inFlatSeg) = 0.0;
        d2gds2(inFlatSeg) = 0.0;

        zS = zC + hRise .* gVal;
        dzdy = sign(yS) .* (hRise / wOut) .* dgds;
        d2z = (hRise / (wOut^2)) .* d2gds2;

        flat = inFlatSeg | (d2z <= 1e-12) | (hRise <= 0);
        Rc = inf(1, Ny);
        Rc(~flat) = (1.0 + dzdy(~flat).^2).^1.5 ./ d2z(~flat);
        flat = flat | (Rc > 1e4 * L);               % numerically planar: use the exact wedge solution
        Rc(flat) = inf;

        % Unit normal toward the centre of curvature (body side, upward/inboard)
        nrm = sqrt(1.0 + dzdy.^2);
        ncY = -dzdy ./ nrm;
        ncZ = 1.0 ./ nrm;

        % Evade evolute cusps: require d_LE <= 0.82*R_c AND x_LE monotonically non-decreasing outboard
        % so dy_LE/dy_S = 1 - d_LE/R_c - sin(psi)*(dd_LE/dy_S) >= 0.18 > 0 everywhere.
        dLE = (L - xLE) .* tanBeta;
        for k = iC + 1:Ny
            if ~flat(k) && dLE(k) > 0.82 * Rc(k)
                dLE(k) = 0.82 * Rc(k);
                xLE(k) = L - dLE(k) / tanBeta;
                clampCount = clampCount + 1;
            end
            if xLE(k) < xLE(k - 1)
                xLE(k) = xLE(k - 1);
                dLE(k) = (L - xLE(k)) * tanBeta;
                clampCount = clampCount + 1;
            end
        end
        xLE(1:iC-1) = fliplr(xLE(iC+1:Ny));
        xLE([1, Ny]) = L;
        xLE(iC) = 0.0;

        X_low = repmat(xLE, Nx, 1) + xi .* repmat(L - xLE, Nx, 1);
        Y_low = zeros(Nx, Ny);
        Z_low = zeros(Nx, Ny);
        Cp_low = zeros(Nx, Ny);
        X_shock = X_low; Y_shock = zeros(Nx, Ny); Z_shock = zeros(Nx, Ny);

        % Flat-ICC columns: planar oblique shock z = -x tan(beta), wedge flow
        kf = find(flat);
        zLEf = -xLE(kf) .* tanBeta;
        Y_low(:, kf) = repmat(yS(kf), Nx, 1);
        Z_low(:, kf) = repmat(zLEf, Nx, 1) - (X_low(:, kf) - xLE(kf)) .* tanDelta;
        Cp_low(:, kf) = CpWedge;
        Y_shock(:, kf) = Y_low(:, kf);
        Z_shock(:, kf) = -X_low(:, kf) .* tanBeta;

        % Curved-ICC columns: local Taylor-Maccoll cones in osculating planes
        kc = find(~flat);
        if ~isempty(kc)
            xApex = L - Rc(kc) ./ tanBeta;
            xStartLocal = xLE(kc) - xApex;
            xEndLocal = L - xApex;
            [xLoc, rGrid, thetaGrid] = waverider.flow.traceStreamline(tm, xStartLocal, xEndLocal, Nx);
            cY = yS(kc) + Rc(kc) .* ncY(kc);       % centre of curvature (base plane)
            cZ = zS(kc) + Rc(kc) .* ncZ(kc);
            X_low(:, kc) = xLoc + xApex;
            Y_low(:, kc) = cY - rGrid .* ncY(kc);
            Z_low(:, kc) = cZ - rGrid .* ncZ(kc);
            Cp_low(:, kc) = interp1(tm.thetaRad, tm.Cp, thetaGrid, 'pchip', tm.Cp(end));
            rS = (X_low(:, kc) - xApex) .* tanBeta;
            Y_shock(:, kc) = cY - rS .* ncY(kc);
            Z_shock(:, kc) = cZ - rS .* ncZ(kc);
        end

        % Exact symmetry plane and singular wingtips
        Y_low(:, iC) = 0.0;
        Y_low(:, [1, Ny]) = repmat(Y_low(1, [1, Ny]), Nx, 1);
        Z_low(:, [1, Ny]) = repmat(Z_low(1, [1, Ny]), Nx, 1);
    end

    if clampCount > 0
        warning('waverider:PlanformClamped', ...
            '%d leading-edge stations adjusted to keep stream surfaces from folding (PlanformPower = %.2f).', ...
            clampCount, cfg.PlanformPower);
    end

    [X_up, Y_up, Z_up] = waverider.geometry.buildUpperSurface(X_low, Y_low, Z_low, cfg);

    surfData = struct();
    surfData.Method = cfg.Method;
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
    surfData.FlowSolution = tm;
    surfData.PlanformClampCount = clampCount;
end
