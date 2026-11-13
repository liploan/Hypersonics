function rtiData = computeRangeProfile(body, freqGHz, bandwidthMHz, azVecDeg, elDeg, varargin)
%COMPUTERANGEPROFILE Range-Angle-Intensity (RTI / sinogram) & High-Resolution Range Profiles.
%   rtiData = COMPUTERANGEPROFILE(body, freqGHz, bandwidthMHz, azVecDeg, elDeg, Name, Value)
%   returns the pulse-compressed down-range profile s(R, psi) of a waverider
%   for a waveform of centre frequency freqGHz and bandwidth bandwidthMHz.
%
%   Model (narrowband stepped-frequency / matched-filter equivalent):
%     Each scattering centre j (every lit PO facet, and every GO-PO
%     double-bounce overlap patch) has a complex amplitude a_j evaluated at
%     f_c with exactly the same physics as computeRCS, and an apparent
%     down-range R_j = (one-way path length) - u . r_ref. The compressed
%     profile is
%         s(R) = sum_j a_j * h(R - R_j)
%     where h is the analytic transform of the spectral window over the band:
%         rect    : h(x) = sinc(x)                               -3 dB 0.886 dR, PSL -13.3 dB
%         hann    : h(x) = [0.5 sinc(x) + 0.25 sinc(x-1) + 0.25 sinc(x+1)] / 0.5
%                                                                  -3 dB 1.44  dR, PSL -31.5 dB
%         hamming : h(x) = [0.54 sinc(x) + 0.23 sinc(x-1) + 0.23 sinc(x+1)] / 0.54
%                                                                  -3 dB 1.30  dR, PSL -42.7 dB
%     with x = (R - R_j) / dR, dR = c / (2 B) and sinc(x) = sin(pi x)/(pi x).
%     Normalising by the coherent gain makes an isolated scatterer peak at
%     its own RCS. For point scatterers this is identical to windowing a
%     stepped-frequency sweep and inverse-Fourier transforming it; it
%     assumes scattering amplitudes are frequency-independent across the
%     band, so fractional bandwidths B / f_c above ~20-30% are approximate
%     and B > f_c is rejected.
%
%   Down-range R [m] is referenced to the vehicle axial midpoint r_ref;
%   negative R is toward the radar. Double-bounce returns appear delayed by
%   half the extra internal path across the concave underside.

    p = inputParser;
    addParameter(p, 'IncludeMultiBounce', true, @islogical);
    addParameter(p, 'Polarization', 'Mean', @(x) ischar(x) || isstring(x));
    addParameter(p, 'WindowType', 'hann', @(x) ischar(x) || isstring(x));
    addParameter(p, 'NumRangeBins', 241, @(x) isnumeric(x) && isscalar(x) && x >= 32);
    addParameter(p, 'PrecomputedConcavePairs', [], @(x) isstruct(x) || isempty(x));
    parse(p, varargin{:});
    opts = p.Results;
    polTag = normalizePolarizationTag(opts.Polarization);
    hasSingleBouncePol = ismember(polTag, {'Mean', 'HH', 'VV', 'PP'});

    c0 = 299792458.0; % [m/s]
    if bandwidthMHz <= 0
        error('computeRangeProfile:Bandwidth', 'Bandwidth must be positive.');
    end
    if bandwidthMHz > 1000 * freqGHz
        error('computeRangeProfile:Bandwidth', ...
            'Bandwidth %.0f MHz exceeds the centre frequency %.2f GHz (band would include f <= 0).', ...
            bandwidthMHz, freqGHz);
    end
    B_Hz = bandwidthMHz * 1e6;
    lambda0 = c0 / (freqGHz * 1e9);
    k0 = 2.0 * pi / lambda0;
    deltaR_m = c0 / (2.0 * B_Hz);

    winType = lower(char(opts.WindowType));
    switch winType
        case 'hann',    a0 = 0.50; a1 = 0.25; w3 = 1.44;
        case 'hamming', a0 = 0.54; a1 = 0.23; w3 = 1.30;
        case 'rect',    a0 = 1.00; a1 = 0.00; w3 = 0.886;
        otherwise
            error('computeRangeProfile:Window', 'WindowType must be rect, hann or hamming.');
    end
    effectiveDeltaR_m = w3 * deltaR_m;   % -3 dB mainlobe width
    psf = @(x) (a0 .* nsinc(x) + a1 .* (nsinc(x - 1) + nsinc(x + 1))) ./ a0;

    % Geometry & reference centre
    TR = body.MeshData.TR;
    V = TR.Points;
    F = TR.ConnectivityList;
    xMin = min(V(:, 1));
    xMax = max(V(:, 1));
    totalLength = xMax - xMin;
    rRef = [0.5 * (xMin + xMax), 0.0, 0.0];

    v1 = V(F(:, 1), :);
    v2 = V(F(:, 2), :);
    v3 = V(F(:, 3), :);
    faceCenters = (v1 + v2 + v3) ./ 3.0;
    crossVec = 0.5 .* cross(v2 - v1, v3 - v1, 2);
    faceAreas = sqrt(sum(crossVec.^2, 2)) + 1e-15;
    faceNormals = crossVec ./ faceAreas;

    % Down-range grid: fine enough to sample the PSF (<= dR/8)
    maxSpan = max(abs(V(:, 2)));
    rLim = max(0.6 * totalLength, maxSpan * 1.1) + 4.0 * effectiveDeltaR_m + 0.3;
    nR = max(round(opts.NumRangeBins), ceil(2 * rLim / (deltaR_m / 8)) + 1);
    nR = min(nR, 20001);
    rangeVec_m = linspace(-rLim, rLim, nR);
    dBin = rangeVec_m(2) - rangeVec_m(1);
    halfK = ceil(6.0 * deltaR_m / dBin);
    kern = psf((-halfK:halfK) .* dBin ./ deltaR_m);

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
    R_nose_m = zeros(nAz, 1);
    R_waveriderBase_m = zeros(nAz, 1);
    R_aftBase_m = zeros(nAz, 1);

    for iAz = 1:nAz
        azRad = deg2rad(azVecDeg(iAz));
        u = [-cos(elRad) * cos(azRad), cos(elRad) * sin(azRad), sin(elRad)];
        eH = [sin(azRad), cos(azRad), 0.0];

        R_nose_m(iAz)          = -([xMin, 0, 0] - rRef) * u.';
        R_waveriderBase_m(iAz) = -([body.Config.Length, 0, 0] - rRef) * u.';
        R_aftBase_m(iAz)       = -([xMax, 0, 0] - rRef) * u.';

        % 1. Single-bounce PO scattering centres (zero in HV/VH/OP/LL channels)
        if hasSingleBouncePol
            cosTh = faceNormals * u.';
            lit = find(cosTh > 0);
            Ilit = waverider.rcs.facetPhaseIntegral(v1(lit, :), v2(lit, :), v3(lit, :), faceNormals(lit, :), 2.0 * k0 .* u);
            aS = (1i * k0 / sqrt(pi)) .* cosTh(lit) .* Ilit;
            RS = -(faceCenters(lit, :) - rRef) * u.';
            profSingle = splatAndCompress(aS, RS, rangeVec_m, kern);
        else
            profSingle = zeros(1, nR);
        end

        % 2. GO-PO concave double-bounce centres with true cavity path delay
        profDouble = zeros(1, nR);
        if opts.IncludeMultiBounce && concavePairs.numPairs > 0
            [~, con] = waverider.rcs.doubleBounceGOPO(concavePairs, u, k0, true, eH);
            rows = con{1};
            if ~isempty(rows)
                switch polTag
                    case 'Mean'
                        aD = rows(:, 1);
                    case 'HH'
                        aD = rows(:, 6);
                    case 'VV'
                        aD = rows(:, 7);
                    case {'HV', 'VH'}
                        aD = rows(:, 8);
                    case 'PP'
                        aD = 0.5 .* (rows(:, 6) + rows(:, 7)) + 0.5i .* (rows(:, 9) - rows(:, 8));
                    case 'OP'
                        aD = 0.5 .* (rows(:, 6) - rows(:, 7)) - 0.5i .* (rows(:, 8) + rows(:, 9));
                    case 'LL'
                        aD = 0.5 .* (rows(:, 6) - rows(:, 7)) + 0.5i .* (rows(:, 8) + rows(:, 9));
                end
                RD = real(rows(:, 2)) + u * rRef.';
                profDouble = splatAndCompress(aD, RD, rangeVec_m, kern);
            end
        end

        noiseFloor = 1e-5;   % -50 dBsm display floor
        RTI_single_m2(iAz, :) = abs(profSingle).^2 + noiseFloor;
        RTI_double_m2(iAz, :) = abs(profDouble).^2 + noiseFloor;
        RTI_total_m2(iAz, :)  = abs(profSingle + profDouble).^2 + noiseFloor;
    end

    rtiData = struct();
    rtiData.freqGHz            = freqGHz;
    rtiData.bandwidthMHz       = bandwidthMHz;
    rtiData.fractionalBW       = bandwidthMHz / (1000 * freqGHz);
    rtiData.polarization       = polTag;
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

function prof = splatAndCompress(a, R, rangeVec, kern)
%SPLATANDCOMPRESS Linear (2-tap) deposition of complex scatterers onto the
%   fine range grid followed by convolution with the sampled compressed PSF.
    nR = numel(rangeVec);
    dBin = rangeVec(2) - rangeVec(1);
    f = (R - rangeVec(1)) ./ dBin + 1;
    i0 = floor(f);
    t = f - i0;
    ok = i0 >= 1 & i0 < nR;
    acc = accumarray(i0(ok), a(ok) .* (1 - t(ok)), [nR, 1]) + ...
          accumarray(i0(ok) + 1, a(ok) .* t(ok), [nR, 1]);
    prof = conv(acc.', kern, 'same');
end

function y = nsinc(x)
%NSINC Normalised sinc sin(pi x)/(pi x).
    y = ones(size(x));
    nz = abs(x) > 1e-12;
    y(nz) = sin(pi .* x(nz)) ./ (pi .* x(nz));
end

function tag = normalizePolarizationTag(polIn)
    s = upper(strtrim(char(polIn)));
    switch s
        case {'MEAN', 'SCALAR', 'UNPOL'}
            tag = 'Mean';
        case 'HH'
            tag = 'HH';
        case 'VV'
            tag = 'VV';
        case 'HV'
            tag = 'HV';
        case 'VH'
            tag = 'VH';
        case {'PP', 'RL', 'LR'}
            tag = 'PP';
        case {'OP', 'RR'}
            tag = 'OP';
        case 'LL'
            tag = 'LL';
        otherwise
            error('computeRangeProfile:Polarization', ...
                'Unsupported Polarization "%s". Choose Mean, HH, VV, HV, VH, PP (RL/LR), OP (RR), or LL.', ...
                char(polIn));
    end
end
