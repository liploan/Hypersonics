function rtiData = computeRangeProfile(body, freqGHz, bandwidthMHz, azVecDeg, elDeg, varargin)
%COMPUTERANGEPROFILE Compute Range-Angle-Intensity (RTI / Sinogram) and High-Resolution Range Profiles.
%   rtiData = COMPUTERANGEPROFILE(body, freqGHz, bandwidthMHz, azVecDeg, elDeg, Name, Value)
%   evaluates the pulse-compressed down-range scattering profile s(R, psi)
%   of a watertight hypersonic waverider across aspect angles azVecDeg for
%   a radar waveform with center frequency freqGHz and bandwidth bandwidthMHz.
%
%   Waveform & Range Resolution Physics:
%     - Speed of light: c0 = 299,792,458 m/s
%     - Unwindowed range resolution: deltaR_m = c0 / (2 * B)
%     - Down-range coordinate R [m] is referenced to the vehicle axial midpoint
%       r_ref = [0.5 * L_total, 0, 0]:
%         * Single-bounce facet f:
%             R_f(psi) = -(r_f - r_ref) . u_radar(psi, el)
%         * Concave underside double-bounce pair (P, Q):
%             R_PQ(psi) = -0.5 * ((r_P + r_Q) - 2*r_ref) . u_radar(psi, el) + 0.5 * d_PQ
%       Notice the +0.5 * d_PQ extra path delay produced by the internal
%       reflection across the concave lower compression cavity.

    p = inputParser;
    addParameter(p, 'IncludeMultiBounce', true, @islogical);
    addParameter(p, 'WindowType', 'hann', @(x) ischar(x) || isstring(x));
    addParameter(p, 'NumRangeBins', 241, @(x) isnumeric(x) && isscalar(x) && x >= 32);
    addParameter(p, 'PrecomputedConcavePairs', [], @(x) isstruct(x) || isempty(x));
    parse(p, varargin{:});
    opts = p.Results;

    c0 = 299792458.0; % [m/s]
    B_Hz = max(bandwidthMHz, 1.0) * 1e6;
    f0_Hz = freqGHz * 1e9;
    lambda0 = c0 / f0_Hz;
    k0 = 2.0 * pi / lambda0;
    deltaR_m = c0 / (2.0 * B_Hz);

    % Window broadening factor and sidelobe envelope
    winType = lower(char(opts.WindowType));
    switch winType
        case 'hann'
            winBroadening = 1.44;
            sidelobeDecay = 2.6;
        case 'hamming'
            winBroadening = 1.30;
            sidelobeDecay = 2.1;
        otherwise % 'rect'
            winBroadening = 1.00;
            sidelobeDecay = 0.85;
    end
    effectiveDeltaR_m = winBroadening * deltaR_m;

    % Vehicle geometry & reference center
    TR = body.MeshData.TR;
    V = TR.Points;
    F = TR.ConnectivityList;
    xMin = min(V(:, 1));
    xMax = max(V(:, 1));
    totalLength = xMax - xMin;
    rRef = [0.5 * (xMin + xMax), 0.0, 0.0];

    % Down-range grid [m] relative to vehicle midpoint
    maxSpan = max(abs(V(:, 2)));
    rLim = max(0.68 * totalLength, maxSpan * 1.15) + 1.2;
    nR = round(opts.NumRangeBins);
    rangeVec_m = linspace(-rLim, rLim, nR);

    v1 = V(F(:, 1), :);
    v2 = V(F(:, 2), :);
    v3 = V(F(:, 3), :);
    faceCenters = (v1 + v2 + v3) ./ 3.0;
    crossVec = 0.5 .* cross(v2 - v1, v3 - v1, 2);
    faceAreas = sqrt(sum(crossVec.^2, 2)) + 1e-15;
    faceNormals = crossVec ./ faceAreas;
    relCenters = faceCenters - rRef;

    % Concave lower-surface pairs for double-bounce range delay
    if opts.IncludeMultiBounce
        if isempty(opts.PrecomputedConcavePairs)
            concavePairs = waverider.rcs.buildConcavePairList(body);
        else
            concavePairs = opts.PrecomputedConcavePairs;
        end
    else
        concavePairs = struct('numPairs', 0);
    end

    azVecDeg = azVecDeg(:);
    nAz = numel(azVecDeg);
    elRad = deg2rad(elDeg);

    RTI_total_m2  = zeros(nAz, nR);
    RTI_single_m2 = zeros(nAz, nR);
    RTI_double_m2 = zeros(nAz, nR);

    % Landmark range tracks across azimuth [nAz x 1]
    R_nose_m = zeros(nAz, 1);
    R_waveriderBase_m = zeros(nAz, 1);
    R_aftBase_m = zeros(nAz, 1);

    prefactor = k0 / sqrt(pi);
    lobeSharpness = min(max(0.35 * k0, 3.5), 35.0);
    elecScale = min((1.0 / (0.4 * lambda0)).^0.6, 2.8);

    for iAz = 1:nAz
        azRad = deg2rad(azVecDeg(iAz));
        u = [-cos(elRad) * cos(azRad); ...
              cos(elRad) * sin(azRad); ...
              sin(elRad)];

        % Geometric landmark projections along line of sight
        R_nose_m(iAz)          = -([0.0, 0.0, 0.0] - rRef) * u;
        R_waveriderBase_m(iAz) = -([body.Config.Length, 0.0, 0.0] - rRef) * u;
        R_aftBase_m(iAz)       = -([xMax, 0.0, 0.0] - rRef) * u;

        % 1. Single-bounce lit facets
        cosTh = faceNormals * u;
        litIdx = find(cosTh > 1e-4);
        cosLit = cosTh(litIdx);
        areaLit = faceAreas(litIdx);
        rDotU_lit = relCenters(litIdx, :) * u;
        R_lit = -rDotU_lit; % Down-range [m]

        damp = 1.0 ./ (1.0 + 0.08 .* k0 .* sqrt(areaLit) .* (1.0 - cosLit));
        ampSingle = (2.0 * prefactor) .* areaLit .* cosLit .* damp .* exp(1i .* (2.0 * k0 .* rDotU_lit));

        % Evaluate band-limited PSF across range bins [nLit x nR]
        dR_single = (rangeVec_m - R_lit) ./ effectiveDeltaR_m;
        psfSingle = sinc_rad(pi .* dR_single) ./ (1.0 + (abs(dR_single) ./ sidelobeDecay).^2.5);
        profSingle = sum(ampSingle .* psfSingle, 1); % [1 x nR]

        % 2. Concave underside double-bounce pairs with +0.5*d_PQ cavity delay
        profDouble = zeros(1, nR);
        if opts.IncludeMultiBounce && concavePairs.numPairs > 0
            cosP = concavePairs.nP * u;
            cosQ = concavePairs.nQ * u;
            dDotU = concavePairs.dPQ * u;
            alignPQ = -dDotU + 2.0 .* cosP .* concavePairs.nP_dot_dPQ;
            alignQP =  dDotU + 2.0 .* cosQ .* concavePairs.nQ_dot_dQP;
            alignMax = max(alignPQ, alignQP);

            validMB = find((cosP > 0.02) & (cosQ > 0.02) & (alignMax > 0.62));
            if ~isempty(validMB)
                wLobe = exp(-lobeSharpness .* (1.0 - alignMax(validMB)));
                obliq = sqrt(cosP(validMB) .* cosQ(validMB));
                rMidRel = 0.5 .* (concavePairs.rP(validMB, :) + concavePairs.rQ(validMB, :)) - rRef;
                rMidDotU = rMidRel * u;
                distMB = concavePairs.distPQ(validMB);

                % Apparent monostatic down-range includes +0.5 * distMB cavity bounce delay
                R_mb = -rMidDotU + 0.5 .* distMB;
                phaseMB = k0 .* (2.0 .* rMidDotU - distMB);
                ampMB = (3.2 * prefactor) .* concavePairs.effArea(validMB) .* elecScale .* obliq .* wLobe .* exp(1i .* phaseMB);

                dR_mb = (rangeVec_m - R_mb) ./ effectiveDeltaR_m;
                psfMB = sinc_rad(pi .* dR_mb) ./ (1.0 + (abs(dR_mb) ./ sidelobeDecay).^2.5);
                profDouble = sum(ampMB .* psfMB, 1);
            end
        end

        noiseFloor = 1e-5;
        RTI_single_m2(iAz, :) = abs(profSingle).^2 + noiseFloor;
        RTI_double_m2(iAz, :) = abs(profDouble).^2 + noiseFloor;
        RTI_total_m2(iAz, :)  = abs(profSingle + profDouble).^2 + noiseFloor;
    end

    rtiData = struct();
    rtiData.freqGHz            = freqGHz;
    rtiData.bandwidthMHz       = bandwidthMHz;
    rtiData.windowType         = winType;
    rtiData.deltaR_m           = deltaR_m;
    rtiData.effectiveDeltaR_m  = effectiveDeltaR_m;
    rtiData.azVecDeg           = azVecDeg;
    rtiData.elDeg              = elDeg;
    rtiData.rangeVec_m         = rangeVec_m;
    rtiData.RTI_total_m2       = RTI_total_m2;
    rtiData.RTI_single_m2      = RTI_single_m2;
    rtiData.RTI_double_m2      = RTI_double_m2;
    rtiData.RTI_total_dBsm     = 10.0 .* log10(RTI_total_m2);
    rtiData.RTI_single_dBsm    = 10.0 .* log10(RTI_single_m2);
    rtiData.RTI_double_dBsm    = 10.0 .* log10(RTI_double_m2);
    rtiData.R_nose_m           = R_nose_m;
    rtiData.R_waveriderBase_m  = R_waveriderBase_m;
    rtiData.R_aftBase_m        = R_aftBase_m;
    rtiData.rRef               = rRef;
end

function y = sinc_rad(x)
%SINC_RAD Unnormalized sinc function sin(x)/x.
    y = ones(size(x));
    nz = abs(x) > 1e-8;
    y(nz) = sin(x(nz)) ./ x(nz);
end
