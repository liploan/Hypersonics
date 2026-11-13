function tm = taylorMaccoll(M1, betaDeg, gamma)
%TAYLORMACCOLL Solve exact Taylor-Maccoll equations for supersonic conical flow.
%   tm = TAYLORMACCOLL(M1, betaDeg, gamma) integrates the axisymmetric
%   conical flow equations inward from the shock angle betaDeg [deg] to the
%   cone half-angle thetaConeDeg [deg] where normal velocity u_theta = 0.
%
%   Returns a struct containing the flowfield profile (thetaRad, ur, utheta,
%   ux, ucyl, Mach, p_p1) and fast interpolants for streamline tracing.

    if nargin < 3
        gamma = 1.4;
    end

    shock = waverider.flow.obliqueShock(M1, betaDeg, gamma);
    betaRad = shock.betaRad;
    deltaRad = shock.deltaRad;
    V2_bar = shock.V2_bar;

    % Initial conditions immediately behind the conical shock at theta = beta
    ur0 = V2_bar * cos(betaRad - deltaRad);
    utheta0 = -V2_bar * sin(betaRad - deltaRad);
    y0 = [ur0; utheta0];

    % Integrate inward from theta = betaRad toward theta = 1e-4 rad
    odeOpts = odeset('Events', @coneSurfaceEvent, 'RelTol', 1e-9, 'AbsTol', 1e-11, 'MaxStep', 1e-3);
    [thetaSol, ySol] = ode45(@(th, y) tmODE(th, y, gamma), [betaRad, 1e-4], y0, odeOpts);

    % Ensure monotonically increasing theta from cone surface (theta_c) to shock (beta)
    thetaRad = flipud(thetaSol(:));
    ur = flipud(ySol(:, 1));
    utheta = flipud(ySol(:, 2));
    utheta(1) = 0.0; % Enforce exact tangency at cone surface

    % Cylindrical velocity components: axial (x) and radial (r)
    ux = ur .* cos(thetaRad) - utheta .* sin(thetaRad);
    ucyl = ur .* sin(thetaRad) + utheta .* cos(thetaRad);

    % Local Mach number and pressure ratio across the shock layer
    V_sq = ur.^2 + utheta.^2;
    M_local = sqrt((2.0 ./ (gamma - 1.0)) .* (V_sq ./ max(1.0 - V_sq, 1e-12)));

    % Isentropic relation behind the conical shock (constant stagnation pressure behind straight shock)
    T0_T2 = 1.0 + 0.5 * (gamma - 1.0) * shock.M2^2;
    T0_Tlocal = 1.0 + 0.5 * (gamma - 1.0) * M_local.^2;
    p_p2 = (T0_T2 ./ T0_Tlocal).^(gamma / (gamma - 1.0));
    p_p1 = shock.p2_p1 .* p_p2;

    tm = struct();
    tm.M1 = M1;
    tm.gamma = gamma;
    tm.betaDeg = betaDeg;
    tm.betaRad = betaRad;
    tm.thetaConeRad = thetaRad(1);
    tm.thetaConeDeg = rad2deg(thetaRad(1));
    tm.shock = shock;
    tm.thetaRad = thetaRad;
    tm.thetaDeg = rad2deg(thetaRad);
    tm.ur = ur;
    tm.utheta = utheta;
    tm.ux = ux;
    tm.ucyl = ucyl;
    tm.slopeRatio = ucyl ./ ux; % dr/dx = ucyl(theta) / ux(theta)
    tm.Mach = M_local;
    tm.p_p1 = p_p1;
    tm.Cp = (p_p1 - 1.0) ./ (0.5 * gamma * M1^2);
    tm.CpCone = tm.Cp(1);
end

function dydt = tmODE(theta, y, gamma)
    ur = y(1);
    ut = y(2);
    a_bar_sq = 0.5 * (gamma - 1.0) * (1.0 - ur^2 - ut^2);
    num = a_bar_sq * (2.0 * ur + ut * cot(theta)) - ur * ut^2;
    den = ut^2 - a_bar_sq;
    if abs(den) < 1e-12
        den = sign(den + eps) * 1e-12;
    end
    dydt = [ut; num / den];
end

function [value, isterminal, direction] = coneSurfaceEvent(~, y)
    value = y(2);      % Stop when u_theta = 0 (cone surface reached)
    isterminal = 1;    % Terminate integration
    direction = 1;     % u_theta approaches 0 from negative values
end
