function surfData = appendConformalFrustum(surfData, cfg)
%APPENDCONFORMALFRUSTUM Append a conformal transition frustum and conical rocket body.
%   surfData = APPENDCONFORMALFRUSTUM(surfData, cfg) smoothly lofts the
%   non-axisymmetric waverider trailing edge at X = L into an axisymmetric
%   circular booster mount of radius cfg.BoosterRadius over length
%   cfg.FrustumLength, optionally followed by a conical rocket stage of
%   length cfg.RocketStageLength and half-angle cfg.RocketConeAngleDeg.

    [NxWave, Ny] = size(surfData.X_low);
    surfData.NxWaverider = NxWave;
    surfData.HasAftFrustum = false;
    surfData.HasRocketStage = false;
    surfData.NxFrustumEnd = NxWave;
    surfData.NxTotal = NxWave;

    if ~cfg.EnableAftFrustum || cfg.FrustumLength <= 0
        return;
    end

    surfData.HasAftFrustum = true;
    L = cfg.Length;
    L_frust = cfg.FrustumLength;
    R_boost = cfg.BoosterRadius;

    % Waverider trailing edge profiles at X = L
    yTE_up  = surfData.Y_up(NxWave, :);
    zTE_up  = surfData.Z_up(NxWave, :);
    yTE_low = surfData.Y_low(NxWave, :);
    zTE_low = surfData.Z_low(NxWave, :);

    % Compute area-weighted vertical centroid of waverider base if FrustumZOffset is NaN
    if isnan(cfg.FrustumZOffset)
        dzBase = max(zTE_up - zTE_low, 0.0);
        baseArea = trapz(yTE_up, dzBase);
        if baseArea > 1e-9
            zCenter = trapz(yTE_up, 0.5 .* (zTE_up.^2 - zTE_low.^2)) / baseArea;
        else
            zCenter = 0.5 * (mean(zTE_up) + mean(zTE_low));
        end
    else
        zCenter = cfg.FrustumZOffset;
    end
    surfData.BoosterZCenter = zCenter;
    surfData.BoosterRadius = R_boost;

    % Parameterize circular booster interface ring at X = L + L_frust
    % Maps spanwise index k = 1..Ny (eta in [-1, 1]) to upper & lower semicircles
    eta = linspace(-1.0, 1.0, Ny);
    thetaUp  = 0.5 * pi .* (1.0 - eta);  % pi (left) -> pi/2 (top) -> 0 (right)
    thetaLow = -0.5 * pi .* (1.0 - eta); % -pi (left) -> -pi/2 (bottom) -> 0 (right)

    yRing_up  = R_boost .* cos(thetaUp);
    zRing_up  = zCenter + R_boost .* sin(thetaUp);
    yRing_low = R_boost .* cos(thetaLow);
    zRing_low = zCenter + R_boost .* sin(thetaLow);

    % Enforce exact left/right seam identity at k = 1 and k = Ny
    yRing_low(1)  = yRing_up(1);
    zRing_low(1)  = zRing_up(1);
    yRing_low(Ny) = yRing_up(Ny);
    zRing_low(Ny) = zRing_up(Ny);

    % Longitudinal stations along conformal frustum (excluding tau = 0 which is row NxWave)
    Nf = max(cfg.NumFrustumStations, 4);
    tauFull = linspace(0.0, 1.0, Nf).';
    tau = tauFull(2:end); % [(Nf-1) x 1]

    switch lower(cfg.FrustumBlendLaw)
        case 'cubic'
            sigma = 3.0 .* (tau.^2) - 2.0 .* (tau.^3);
        case 'quintic'
            sigma = 10.0 .* (tau.^3) - 15.0 .* (tau.^4) + 6.0 .* (tau.^5);
        case 'linear'
            sigma = tau;
        case 'conic'
            pExp = 1.45;
            sigma = (tau.^pExp) ./ ((tau.^pExp) + ((1.0 - tau).^pExp));
        otherwise
            sigma = 3.0 .* (tau.^2) - 2.0 .* (tau.^3);
    end

    X_f = repmat(L + tau .* L_frust, 1, Ny);
    Y_f_up  = (1.0 - sigma) * yTE_up  + sigma * yRing_up;
    Z_f_up  = (1.0 - sigma) * zTE_up  + sigma * zRing_up;
    Y_f_low = (1.0 - sigma) * yTE_low + sigma * yRing_low;
    Z_f_low = (1.0 - sigma) * zTE_low + sigma * zRing_low;

    % Enforce exact left/right seam coincidence along frustum
    Y_f_low(:, 1)  = Y_f_up(:, 1);
    Z_f_low(:, 1)  = Z_f_up(:, 1);
    Y_f_low(:, Ny) = Y_f_up(:, Ny);
    Z_f_low(:, Ny) = Z_f_up(:, Ny);

    % Append frustum rows to surface matrices
    surfData.X_up  = [surfData.X_up;  X_f];
    surfData.Y_up  = [surfData.Y_up;  Y_f_up];
    surfData.Z_up  = [surfData.Z_up;  Z_f_up];
    surfData.X_low = [surfData.X_low; X_f];
    surfData.Y_low = [surfData.Y_low; Y_f_low];
    surfData.Z_low = [surfData.Z_low; Z_f_low];

    % Estimate local Cp on frustum lower surface (weak expansion or flare compression)
    Cp_last = surfData.Cp_low(end, :);
    Cp_f = (1.0 - 0.65 .* sigma) * Cp_last;
    surfData.Cp_low = [surfData.Cp_low; Cp_f];

    surfData.NxFrustumEnd = size(surfData.X_low, 1);
    surfData.NxTotal = surfData.NxFrustumEnd;

    % Optional Traditional Conical / Cylindrical Rocket Booster Section
    if cfg.EnableRocketStage && cfg.RocketStageLength > 0
        surfData.HasRocketStage = true;
        L_rock = cfg.RocketStageLength;
        Nr = max(cfg.NumRocketStations, 4);
        tauRFull = linspace(0.0, 1.0, Nr).';
        tauR = tauRFull(2:end); % [(Nr-1) x 1]
        dxR = tauR .* L_rock;

        tanCone = tan(deg2rad(cfg.RocketConeAngleDeg));
        R_local = max(R_boost + dxR .* tanCone, 0.05 * R_boost); % [(Nr-1) x 1]

        X_r = repmat(L + L_frust + dxR, 1, Ny);
        Y_r_up  = R_local * cos(thetaUp);
        Z_r_up  = zCenter + R_local * sin(thetaUp);
        Y_r_low = R_local * cos(thetaLow);
        Z_r_low = zCenter + R_local * sin(thetaLow);

        Y_r_low(:, 1)  = Y_r_up(:, 1);
        Z_r_low(:, 1)  = Z_r_up(:, 1);
        Y_r_low(:, Ny) = Y_r_up(:, Ny);
        Z_r_low(:, Ny) = Z_r_up(:, Ny);

        % Newtonian Cp on conical rocket flare
        sinCone = max(sin(deg2rad(cfg.RocketConeAngleDeg)), 0.0);
        CpCone = repmat(2.0 * (sinCone^2), numel(tauR), Ny);

        surfData.X_up   = [surfData.X_up;   X_r];
        surfData.Y_up   = [surfData.Y_up;   Y_r_up];
        surfData.Z_up   = [surfData.Z_up;   Z_r_up];
        surfData.X_low  = [surfData.X_low;  X_r];
        surfData.Y_low  = [surfData.Y_low;  Y_r_low];
        surfData.Z_low  = [surfData.Z_low;  Z_r_low];
        surfData.Cp_low = [surfData.Cp_low; CpCone];

        surfData.NxTotal = size(surfData.X_low, 1);
        surfData.RocketEndRadius = R_local(end);
    else
        surfData.RocketEndRadius = R_boost;
    end
end
