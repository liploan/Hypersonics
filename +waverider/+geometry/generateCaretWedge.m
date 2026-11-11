function surfData = generateCaretWedge(cfg)
%GENERATECARETWEDGE Generate Nonweiler Caret / 2D oblique-shock wedge waverider.
%   surfData = GENERATECARETWEDGE(cfg) uses exact 2D oblique shock relations
%   (theta-beta-M) to construct a Caret / wedge-derived compression surface
%   with a planar or dihedral-arched shock sheet.

    L = cfg.Length;
    bHalf = cfg.HalfSpan;
    betaDeg = cfg.ShockAngleDeg;
    betaRad = deg2rad(betaDeg);
    tanBeta = tan(betaRad);

    % Exact 2D oblique shock solution
    shock = waverider.flow.obliqueShock(cfg.Mach, betaDeg, cfg.Gamma);
    deltaRad = shock.deltaRad;
    tanDelta = tan(deltaRad);

    Ny = 2 * cfg.NumHalfSpan - 1;
    Nx = cfg.NumStreamwise;
    eta = linspace(-1.0, 1.0, Ny);
    yLE = bHalf .* eta;

    % Leading-edge planform curve: x_LE / L = |eta|^(1 / P)
    invP = 1.0 / max(cfg.PlanformPower, 0.25);
    xLE = L .* (abs(eta).^invP);
    xLE(1) = L;
    xLE(Ny) = L;
    xLE(cfg.NumHalfSpan) = 0.0;

    % Caret leading edge lies on the oblique shock sheet with optional spanwise anhedral
    anhedralRad = deg2rad(cfg.LeadingEdgeAnhedralDeg);
    zLE = -xLE .* tanBeta - abs(yLE) .* tan(anhedralRad);
    zLE(cfg.NumHalfSpan) = 0.0;

    % Streamwise grid from xLE(k) to L
    xi = linspace(0.0, 1.0, Nx).';
    X_low = repmat(xLE, Nx, 1) + xi .* repmat(L - xLE, Nx, 1);
    Y_low = repmat(yLE, Nx, 1);

    % Post-shock streamlines deflect downward by angle delta
    dxStream = X_low - repmat(xLE, Nx, 1);
    Z_low = repmat(zLE, Nx, 1) - dxStream .* tanDelta;

    % Shock sheet trails below at angle beta
    X_shock = X_low;
    Y_shock = Y_low;
    Z_shock = repmat(zLE, Nx, 1) - dxStream .* tanBeta;

    % Build Upper Surface ('freestream' gives classic Caret; 'canopy' adds dorsal loft)
    [X_up, Y_up, Z_up] = waverider.geometry.buildUpperSurface(X_low, Y_low, Z_low, cfg);

    % Uniform post-shock pressure coefficient on planar compression surface
    Cp_val = (shock.p2_p1 - 1.0) / (0.5 * cfg.Gamma * cfg.Mach^2);
    Cp_low = repmat(Cp_val, Nx, Ny);

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
