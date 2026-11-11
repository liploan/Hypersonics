function surfData = generateOsculatingCone(cfg)
%GENERATEOSCULATINGCONE Generate conical or osculating-cone waverider surfaces.
%   surfData = GENERATEOSCULATINGCONE(cfg) solves the Taylor-Maccoll conical
%   flowfield and traces shock-capturing streamlines to form the lower
%   compression surface and upper freestream/canopy surface.

    L = cfg.Length;
    bHalf = cfg.HalfSpan;
    betaDeg = cfg.ShockAngleDeg;
    betaRad = deg2rad(betaDeg);
    tanBeta = tan(betaRad);

    % Solve exact axisymmetric Taylor-Maccoll conical flowfield
    tm = waverider.flow.taylorMaccoll(cfg.Mach, betaDeg, cfg.Gamma);

    % Spanwise discretization from -bHalf to +bHalf (2*NumHalfSpan - 1 stations)
    Ny = 2 * cfg.NumHalfSpan - 1;
    Nx = cfg.NumStreamwise;
    eta = linspace(-1.0, 1.0, Ny); % Normalized spanwise coordinate in [-1, 1]
    yLE = bHalf .* eta;

    % Planform leading-edge longitudinal distribution: x_LE / L = |eta|^(1 / P)
    invP = 1.0 / max(cfg.PlanformPower, 0.25);
    xLE = L .* (abs(eta).^invP);
    xLE(1) = L;
    xLE(Ny) = L;
    xLE(cfg.NumHalfSpan) = 0.0; % Exact symmetry apex at (0,0,0)

    if strcmpi(cfg.Method, 'conical')
        % -----------------------------------------------------------------
        % Pure Axisymmetric Cone-Derived Waverider
        % -----------------------------------------------------------------
        % Ensure cone apex offset x0 is large enough that the shock cone at
        % x = L comfortably encloses HalfSpan (R_shock(L) >= 1.25 * bHalf)
        minR_L = 1.28 * bHalf;
        x0_min = max(minR_L / tanBeta - L, 0.15 * L);
        x0_cfg = (cfg.VertexOffsetRatio / (1.0 - cfg.VertexOffsetRatio)) * L;
        x0 = max(x0_cfg, x0_min);
        zv = x0 * tanBeta; % Cone axis height so shock hits (0,0,0) at x=0

        % Shock radius at each leading-edge station xLE(k)
        xStartLocal = xLE + x0;
        xEndLocal = repmat(L + x0, 1, Ny);
        R_shock_LE = xStartLocal .* tanBeta;

        % Leading edge vertical coordinate on the bottom half of the shock cone
        zDepthLE = sqrt(max(R_shock_LE.^2 - yLE.^2, 1e-12));
        zLE = zv - zDepthLE;
        phiLE = atan2(yLE, zDepthLE); % Azimuthal angle from bottom centerline

        % Trace all streamlines in Taylor-Maccoll flowfield
        [xLocalGrid, rGrid, thetaGrid] = waverider.flow.traceStreamline(tm, xStartLocal, xEndLocal, Nx);

        X_low = xLocalGrid - x0;
        Y_low = rGrid .* sin(repmat(phiLE, Nx, 1));
        Z_low = zv - rGrid .* cos(repmat(phiLE, Nx, 1));

        % Shock sheet for visualization
        R_shock_grid = xLocalGrid .* tanBeta;
        X_shock = X_low;
        Y_shock = R_shock_grid .* sin(repmat(phiLE, Nx, 1));
        Z_shock = zv - R_shock_grid .* cos(repmat(phiLE, Nx, 1));
    else
        % -----------------------------------------------------------------
        % Sobieczky Osculating Cones Waverider (Flat-Center + Curved Outer ICC)
        % -----------------------------------------------------------------
        fFlat = cfg.FlatShockFraction;
        mExp = cfg.ShockCurvaturePower;
        anhedralRad = deg2rad(cfg.LeadingEdgeAnhedralDeg);

        % Local effective vertex offset along the span
        x0_base = (cfg.VertexOffsetRatio / (1.0 - cfg.VertexOffsetRatio)) * L;
        x0_min = max(1.25 * bHalf / tanBeta - L, 0.20 * L);
        x0_ref = max(x0_base, x0_min);

        % Osculating shock curve curvature distribution across span
        absEta = abs(eta);
        outerEta = max(absEta - fFlat, 0.0) ./ max(1.0 - fFlat, 1e-6);

        % Leading-edge vertical coordinate combines shock cone depth + anhedral arch
        zShockCenter = -xLE .* tanBeta;
        zAnhedralOffset = (bHalf * tan(anhedralRad)) .* (outerEta.^mExp);
        zLE = zShockCenter + 0.35 * (xLE ./ L) .* tanBeta .* L - zAnhedralOffset;
        zLE(cfg.NumHalfSpan) = 0.0;

        % Local osculating cone axial distance from local virtual apex
        % Larger radius of curvature in flat center -> larger virtual x0;
        % tighter curvature in outer wing -> smaller virtual x0
        curvatureFactor = 1.0 ./ (0.35 + 0.65 .* (outerEta.^(mExp - 1.0)));
        x0_local = x0_ref .* min(curvatureFactor, 3.0);

        xStartLocal = xLE + x0_local;
        xEndLocal = L + x0_local;

        [~, rGrid, thetaGrid] = waverider.flow.traceStreamline(tm, xStartLocal, xEndLocal, Nx);

        % Radial deflection relative to the shock attachment point r(1,:)
        drStream = rGrid - repmat(rGrid(1, :), Nx, 1);

        % Local osculating plane azimuthal normal tilt psi(k)
        dz_dy = zeros(1, Ny);
        if Ny >= 3
            dz_dy = gradient(zLE, yLE);
        end
        psiOsc = atan(0.35 .* dz_dy); % Smooth normal inclination in osculating planes

        xi = linspace(0.0, 1.0, Nx).';
        X_low = repmat(xLE, Nx, 1) + xi .* repmat(L - xLE, Nx, 1);
        Y_low = repmat(yLE, Nx, 1) - drStream .* sin(repmat(psiOsc, Nx, 1));
        Z_low = repmat(zLE, Nx, 1) - (drStream - xi .* repmat(L - xLE, Nx, 1) .* tanBeta) .* cos(repmat(psiOsc, Nx, 1));

        % Wait: let's verify the compression slope on Z_low!
        % In conical flow, rGrid(i) grows slower than x * tan(beta) because
        % dr/dx < tan(beta): specifically, the shock wave drops at x * tan(beta),
        % while the streamline drops at drStream = r(x) - r(xLE) > 0!
        % Therefore, relative to the freestream (z = zLE), the streamline drops
        % downward by drStream * cos(psiOsc)!
        Z_low = repmat(zLE, Nx, 1) - drStream .* cos(repmat(psiOsc, Nx, 1));

        % Enforce exact symmetry plane y = 0 and wingtip coincidence
        Y_low(:, cfg.NumHalfSpan) = 0.0;
        Y_low(:, 1) = yLE(1);
        Y_low(:, Ny) = yLE(Ny);

        % Shock sheet for visualization
        dxGrid = X_low - repmat(xLE, Nx, 1);
        X_shock = X_low;
        Y_shock = repmat(yLE, Nx, 1);
        Z_shock = repmat(zLE, Nx, 1) - dxGrid .* tanBeta;
    end

    % Build Upper Surface ('freestream' or 'canopy')
    [X_up, Y_up, Z_up] = waverider.geometry.buildUpperSurface(X_low, Y_low, Z_low, cfg);

    % Compute local surface pressure coefficient Cp from Taylor-Maccoll solution
    cpInterp = griddedInterpolant(tm.thetaRad, tm.Cp, 'pchip', 'nearest');
    Cp_low = cpInterp(thetaGrid);

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
end
