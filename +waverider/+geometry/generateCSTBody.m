function surfData = generateCSTBody(cfg)
%GENERATECSTBODY Generate parametric CST (Class-Shape Transformation) boost-glide waverider.
%   surfData = GENERATECSTBODY(cfg) constructs a hypersonic boost-glide body
%   using Kulfan Class-Shape Transformation (CST) Bernstein polynomials for
%   the longitudinal and transverse cross-sections, coupled with a wave-
%   trapping concave lower compression surface.

    L = cfg.Length;
    bHalf = cfg.HalfSpan;
    betaDeg = cfg.ShockAngleDeg;
    betaRad = deg2rad(betaDeg);
    tanBeta = tan(betaRad);

    shock = waverider.flow.obliqueShock(cfg.Mach, betaDeg, cfg.Gamma);
    tanDelta = tan(shock.deltaRad);

    Ny = 2 * cfg.NumHalfSpan - 1;
    Nx = cfg.NumStreamwise;
    eta = linspace(-1.0, 1.0, Ny);
    yLE = bHalf .* eta;

    % Power-law / blended planform leading edge
    invP = 1.0 / max(cfg.PlanformPower, 0.25);
    xLE = L .* (abs(eta).^invP);
    xLE(1) = L;
    xLE(Ny) = L;
    xLE(cfg.NumHalfSpan) = 0.0;

    % Leading-edge vertical curve with droop/anhedral for shock containment
    anhedralRad = deg2rad(cfg.LeadingEdgeAnhedralDeg);
    zLE = -0.65 .* xLE .* tanBeta - abs(yLE) .* tan(anhedralRad) .* (abs(eta).^0.8);
    zLE(cfg.NumHalfSpan) = 0.0;

    % Normalized streamwise coordinate s in [0, 1]
    s = linspace(0.0, 1.0, Nx).';
    dxMax = L - xLE; % [1 x Ny]
    X_low = repmat(xLE, Nx, 1) + s .* repmat(dxMax, Nx, 1);
    Y_low = repmat(yLE, Nx, 1);

    % ---------------------------------------------------------------------
    % Kulfan CST Longitudinal Class & Bernstein Shape Functions
    % ---------------------------------------------------------------------
    N1 = cfg.CSTNoseExponent;
    N2 = cfg.CSTAftExponent;
    % Modified longitudinal class function: 0 at LE (s=0), finite at base (s=1)
    C_long = (s.^N1) .* ((1.0 - 0.30 .* s).^N2);

    weights = cfg.CSTUpperWeights(:).';
    nDeg = numel(weights) - 1;
    S_bernstein = zeros(Nx, 1);
    for r = 0:nDeg
        binomCoeff = nchoosek(nDeg, r);
        basis_r = binomCoeff .* (s.^r) .* ((1.0 - s).^(nDeg - r));
        S_bernstein = S_bernstein + weights(r + 1) .* basis_r;
    end
    cstLongProfile = C_long .* S_bernstein; % [Nx x 1], 0 at s=0

    % Spanwise CST class function: 1 at centerline, 0 at wingtips
    C_span = (1.0 - abs(eta).^cfg.CanopySpanPower).^1.5; % [1 x Ny]
    chordScale = max(1.0 - xLE ./ L, 0.0).^0.65;

    % ---------------------------------------------------------------------
    % Lower Compression Surface: Shock deflection + CST concave wave arch
    % ---------------------------------------------------------------------
    dxStream = X_low - repmat(xLE, Nx, 1);
    % Concave transverse arch traps high-pressure shock layer near centerline
    archBell = (1.0 - eta.^2).^2;
    dzArch = (cfg.CSTLowerArchRatio * L) .* (s.^0.9) .* (archBell .* chordScale);
    Z_low = repmat(zLE, Nx, 1) - dxStream .* tanDelta + dzArch;

    % Ensure Z_low remains strictly below zLE for s > 0
    maxAllowedZLow = repmat(zLE, Nx, 1) - 0.25 .* dxStream .* tanDelta;
    Z_low = min(Z_low, maxAllowedZLow);
    Z_low(1, :) = zLE;

    % ---------------------------------------------------------------------
    % Upper Surface: Freestream or CST Dorsal Canopy
    % ---------------------------------------------------------------------
    X_up = X_low;
    Y_up = Y_low;
    Z_fs = repmat(zLE, Nx, 1);

    if strcmpi(cfg.UpperSurfaceType, 'freestream') || cfg.CanopyHeightRatio <= 0
        Z_up = Z_fs;
    else
        % Normalize CST profile so max height equals CanopyHeightRatio * L
        maxCST = max(cstLongProfile) + 1e-12;
        canopyAmpl = (cfg.CanopyHeightRatio * L) / maxCST;
        dzUpperCST = canopyAmpl .* (cstLongProfile * (C_span .* chordScale));
        Z_up = Z_fs + dzUpperCST;
    end

    % Enforce exact coincidence along LE and wingtips
    X_up(1, :) = X_low(1, :);
    Y_up(1, :) = Y_low(1, :);
    Z_up(1, :) = Z_low(1, :);
    for tipCol = [1, Ny]
        X_up(:, tipCol) = X_low(:, tipCol);
        Y_up(:, tipCol) = Y_low(:, tipCol);
        Z_up(:, tipCol) = Z_low(:, tipCol);
    end

    % Reference shock sheet
    X_shock = X_low;
    Y_shock = Y_low;
    Z_shock = repmat(zLE, Nx, 1) - dxStream .* tanBeta;

    % Modified Newtonian pressure coefficient on lower surface (windward/leeward consistent)
    CpMax = (2.0 / (cfg.Gamma * cfg.Mach^2)) * ...
        (((cfg.Gamma + 1.0)^2 * cfg.Mach^2 / (4.0 * cfg.Gamma * cfg.Mach^2 - 2.0 * (cfg.Gamma - 1.0)))^(cfg.Gamma / (cfg.Gamma - 1.0)) * ...
         ((1.0 - cfg.Gamma + 2.0 * cfg.Gamma * cfg.Mach^2) / (cfg.Gamma + 1.0)) - 1.0);
    [Nx_l, Ny_l, Nz_l] = surfnorm(X_low, Y_low, Z_low);
    flipMask = (Nz_l > 0); % Ensure outward lower-surface normal points downward (-Z)
    Nx_l(flipMask) = -Nx_l(flipMask);
    Ny_l(flipMask) = -Ny_l(flipMask);
    Nz_l(flipMask) = -Nz_l(flipMask);
    nMag = sqrt(Nx_l.^2 + Ny_l.^2 + Nz_l.^2) + 1e-12;
    sinImpact = max(-Nx_l ./ nMag, 0.0);
    sinLee = max(Nx_l ./ nMag, 0.0);
    g = cfg.Gamma;
    M = cfg.Mach;
    CpPM = (2.0 / (g * M^2)) .* (max(1.0 - 0.5 * (g - 1.0) * M .* sinLee, 0.0).^(2.0 * g / (g - 1.0)) - 1.0);
    Cp_low = CpMax .* (sinImpact.^2) + CpPM .* (sinImpact == 0);

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
    surfData.FlowSolution = shock;
end
