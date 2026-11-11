function shock = obliqueShock(M1, betaDeg, gamma)
%OBLIQUESHOCK Exact 2D supersonic/hypersonic oblique shock relations.
%   shock = OBLIQUESHOCK(M1, betaDeg, gamma) computes the flow deflection
%   angle delta, downstream Mach number M2, and thermodynamic ratios across
%   an oblique shock wave at angle betaDeg [deg] for freestream Mach M1.

    if nargin < 3
        gamma = 1.4;
    end

    muRad = asin(1.0 / M1);
    muDeg = rad2deg(muRad);
    betaRad = deg2rad(betaDeg);

    if betaDeg <= muDeg
        error('obliqueShock:SubsonicNormalMach', ...
            'Shock angle beta = %.2f deg must exceed Mach angle mu = %.2f deg (M = %.2f).', ...
            betaDeg, muDeg, M1);
    end

    M1n = M1 * sin(betaRad);

    % Exact theta-beta-M relation for flow deflection angle delta
    num = 2.0 * cot(betaRad) * (M1n^2 - 1.0);
    den = M1^2 * (gamma + cos(2.0 * betaRad)) + 2.0;
    deltaRad = atan(num / den);
    deltaDeg = rad2deg(deltaRad);

    % Normal Mach number downstream
    M2n_sq = (1.0 + 0.5 * (gamma - 1.0) * M1n^2) / (gamma * M1n^2 - 0.5 * (gamma - 1.0));
    M2n = sqrt(max(M2n_sq, 1e-8));
    M2 = M2n / sin(betaRad - deltaRad);

    % Thermodynamic ratios across the shock
    p2_p1 = 1.0 + (2.0 * gamma / (gamma + 1.0)) * (M1n^2 - 1.0);
    rho2_rho1 = ((gamma + 1.0) * M1n^2) / (2.0 + (gamma - 1.0) * M1n^2);
    T2_T1 = p2_p1 / rho2_rho1;

    % Nondimensional velocity ratio V / V_max (where V_max = sqrt(2*cp*T0))
    V1_bar = (1.0 + 2.0 / ((gamma - 1.0) * M1^2))^(-0.5);
    V2_bar = (1.0 + 2.0 / ((gamma - 1.0) * M2^2))^(-0.5);

    shock = struct();
    shock.M1 = M1;
    shock.gamma = gamma;
    shock.betaDeg = betaDeg;
    shock.betaRad = betaRad;
    shock.deltaDeg = deltaDeg;
    shock.deltaRad = deltaRad;
    shock.muDeg = muDeg;
    shock.M1n = M1n;
    shock.M2 = M2;
    shock.p2_p1 = p2_p1;
    shock.rho2_rho1 = rho2_rho1;
    shock.T2_T1 = T2_T1;
    shock.V1_bar = V1_bar;
    shock.V2_bar = V2_bar;
end
