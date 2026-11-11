function [xGrid, rGrid, thetaGrid] = traceStreamline(tm, xStartLocal, xEndLocal, numStreamwise)
%TRACESTREAMLINE Vectorized streamline tracing through a Taylor-Maccoll flowfield.
%   [xGrid, rGrid, thetaGrid] = TRACESTREAMLINE(tm, xStartLocal, xEndLocal, numStreamwise)
%   integrates dr/dx = ucyl(theta) / ux(theta) from the conical shock wave
%   r(xStartLocal) = xStartLocal * tan(beta) downstream to xEndLocal across
%   numStreamwise stations for an arbitrary number of spanwise stations.
%
%   Inputs:
%     tm            - Struct returned by waverider.flow.taylorMaccoll
%     xStartLocal   - [1 x Ny] local axial coordinates where streamlines pierce shock
%     xEndLocal     - [1 x Ny] local axial coordinates at the aft base plane
%     numStreamwise - Number of streamwise grid stations (Nx)
%
%   Outputs:
%     xGrid         - [Nx x Ny] local axial coordinate along each streamline
%     rGrid         - [Nx x Ny] radial distance from local cone axis
%     thetaGrid     - [Nx x Ny] local conical ray angle theta [rad]

    xStartRow = xStartLocal(:).';
    xEndRow = xEndLocal(:).';
    Ny = numel(xStartRow);
    Nx = numStreamwise;

    xi = linspace(0.0, 1.0, Nx).';
    dxTotal = max(xEndRow - xStartRow, 0.0); % [1 x Ny]
    xGrid = xStartRow + xi .* dxTotal;       % [Nx x Ny]

    % Initial radial position on the conical shock wave: r0 = xStart * tan(beta)
    rGrid = zeros(Nx, Ny);
    thetaGrid = zeros(Nx, Ny);

    rCurr = xStartRow .* tan(tm.betaRad);
    rGrid(1, :) = rCurr;
    thetaGrid(1, :) = tm.betaRad;

    % Build fast PCHIP interpolant for slopeRatio(theta) = dr/dx
    slopeInterp = griddedInterpolant(tm.thetaRad, tm.slopeRatio, 'pchip', 'nearest');
    thetaMin = tm.thetaConeRad;
    thetaMax = tm.betaRad;

    % Use RK4 with 4 sub-steps per grid interval for high accuracy
    numSubSteps = 4;
    for i = 1:(Nx - 1)
        dxStep = ( dxTotal ./ (Nx - 1) ) / numSubSteps;
        xSub = xGrid(i, :);
        rSub = rCurr;

        for s = 1:numSubSteps
            th1 = clampTheta(atan2(rSub, max(xSub, 1e-9)), thetaMin, thetaMax);
            k1 = slopeInterp(th1);

            x2 = xSub + 0.5 * dxStep;
            r2 = rSub + 0.5 * dxStep .* k1;
            th2 = clampTheta(atan2(r2, max(x2, 1e-9)), thetaMin, thetaMax);
            k2 = slopeInterp(th2);

            r3 = rSub + 0.5 * dxStep .* k2;
            th3 = clampTheta(atan2(r3, max(x2, 1e-9)), thetaMin, thetaMax);
            k3 = slopeInterp(th3);

            x4 = xSub + dxStep;
            r4 = rSub + dxStep .* k3;
            th4 = clampTheta(atan2(r4, max(x4, 1e-9)), thetaMin, thetaMax);
            k4 = slopeInterp(th4);

            rSub = rSub + (dxStep / 6.0) .* (k1 + 2.0*k2 + 2.0*k3 + k4);
            xSub = x4;
        end

        rCurr = rSub;
        rGrid(i + 1, :) = rCurr;
        thetaGrid(i + 1, :) = clampTheta(atan2(rCurr, max(xGrid(i + 1, :), 1e-9)), thetaMin, thetaMax);
    end
end

function th = clampTheta(th, thMin, thMax)
    th = min(max(th, thMin), thMax);
end
