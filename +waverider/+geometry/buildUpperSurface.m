function [X_up, Y_up, Z_up] = buildUpperSurface(X_low, Y_low, Z_low, cfg)
%BUILDUPPERSURFACE Construct upper surface from leading-edge curve and config.
%   [X_up, Y_up, Z_up] = BUILDUPPERSURFACE(X_low, Y_low, Z_low, cfg) builds
%   either a pure freestream-traced upper surface ('freestream') or a
%   volumetric cambered dorsal canopy ('canopy') that meets the lower
%   compression surface strictly along the entire leading edge (row 1) and
%   at both wingtips (columns 1 and end).

    [Nx, Ny] = size(X_low);

    % Leading edge coordinates (row 1 of lower surface)
    xLE = X_low(1, :);
    yLE = Y_low(1, :);
    zLE = Z_low(1, :);

    % Streamwise coordinate matches lower surface grid
    X_up = X_low;

    % Freestream tracing keeps Y constant from the leading edge, but to
    % ensure smooth spanwise spacing relative to any slight radial flare of
    % the lower surface, Y_up starts at yLE and remains at yLE (since at
    % the wingtips xLE = L so Y_low(end, tip) == yLE(tip)).
    Y_up = repmat(yLE, Nx, 1);

    % Baseline freestream Z surface: constant z = zLE along each streamline
    Z_fs = repmat(zLE, Nx, 1);

    if strcmpi(cfg.UpperSurfaceType, 'freestream') || cfg.CanopyHeightRatio <= 0
        Z_up = Z_fs;
    else
        % Volumetric dorsal canopy lofted smoothly above the freestream surface:
        % - Zero thickness at the leading edge (xi = 0) with smooth ramp
        % - Zero thickness at the wingtips (|y| / HalfSpan = 1)
        % - Maximum dorsal spine height near symmetry plane (y = 0)
        L = cfg.Length;
        bHalf = max(abs(yLE));
        if bHalf < 1e-9
            bHalf = cfg.HalfSpan;
        end

        % Normalized streamwise distance from local leading edge: s in [0, 1]
        dxStream = max(X_up - repmat(xLE, Nx, 1), 0.0);
        dxMax = max(L - xLE, 1e-9);
        sNorm = dxStream ./ repmat(dxMax, Nx, 1);
        % Zero out at tips where xLE == L
        isTip = (L - xLE) < 1e-7;
        sNorm(:, isTip) = 0.0;

        % Spanwise bell distribution: (1 - (|y|/bHalf)^p)^2
        etaSpan = min(abs(Y_up) ./ bHalf, 1.0);
        spanProfile = (1.0 - etaSpan.^cfg.CanopySpanPower).^2;

        % Longitudinal profile: smooth rise from LE (sNorm^0.85 * (1 - 0.25*sNorm))
        % and scaled by local chord fraction (1 - xLE/L) so apex & tips are sharp
        chordScale = repmat(max(1.0 - xLE ./ L, 0.0).^0.75, Nx, 1);
        longProfile = (sNorm.^0.85) .* (1.0 - 0.25 .* sNorm);

        hCanopy = (cfg.CanopyHeightRatio * L) .* chordScale .* spanProfile .* longProfile;
        Z_up = Z_fs + hCanopy;
    end

    % Enforce exact coincidence at the leading edge (row 1) and wingtips (col 1, Ny)
    X_up(1, :) = X_low(1, :);
    Y_up(1, :) = Y_low(1, :);
    Z_up(1, :) = Z_low(1, :);

    X_up(:, 1) = X_low(:, 1);
    Y_up(:, 1) = Y_low(:, 1);
    Z_up(:, 1) = Z_low(:, 1);

    X_up(:, Ny) = X_low(:, Ny);
    Y_up(:, Ny) = Y_low(:, Ny);
    Z_up(:, Ny) = Z_low(:, Ny);
end
