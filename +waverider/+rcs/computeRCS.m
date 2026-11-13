function rcsData = computeRCS(body, freqGHz, azDeg, elDeg, varargin)
%COMPUTERCS Polarimetric monostatic RCS: exact facet PO + vector GO-PO double bounce.
%   rcsData = COMPUTERCS(body, freqGHz, azDeg, elDeg, Name, Value) computes
%   the monostatic Radar Cross Section (RCS) of a watertight waverider body
%   (UHF 0.3 GHz to X-band 12 GHz) across scalar, linear, and circular
%   polarizations (IEEE Backscatter Alignment convention).
%
%   Polarization ('Polarization' Name-Value, default 'Mean'):
%     - 'Mean'          : Scalar PEC baseline (Gamma = -1, polarization-averaged)
%     - 'HH', 'VV'      : Co-polarized linear (Horizontal-Horizontal / Vertical-Vertical)
%     - 'HV', 'VH'      : Cross-polarized linear (Lorentz-reciprocal S_HV = S_VH)
%     - 'PP' ('RL','LR'): Principal circular (opposite-sense, odd-bounce matched)
%     - 'OP' ('RR')     : Orthogonal circular (same-sense RHCP, even-bounce cavity matched)
%     - 'LL'            : Orthogonal circular (same-sense LHCP, even-bounce cavity matched)
%
%   All linear and circular channels are evaluated simultaneously in a single
%   geometric pass and returned in rcsData.Polarimetric with < 3% overhead.
%
%   Coordinate & Angle Convention:
%     - Vehicle nose points toward -X (nose at X=0, base at X=L).
%     - Unit vector FROM target TO radar:
%         u_radar = [-cos(el)*cos(az), cos(el)*sin(az), sin(el)]
%         e_H     = [ sin(az),         cos(az),         0      ]
%         e_V     = e_H x u_radar  (so e_H x e_V = k_i = -u_radar)
%
%   Scattering mechanisms:
%     1. Single-bounce Physical Optics over every lit facet using the EXACT
%        closed-form triangle integral (Gordon 1975; facetPhaseIntegral):
%            S_PO = (i k / sqrt(pi)) sum_f (n_f . u) \int_f exp(i 2k u . r) dA
%        For PEC single-bounce PO: S_HH = S_VV = S_PO, S_HV = S_VH = 0,
%        S_PP = S_PO, S_OP = S_LL = 0.
%     2. Concave-underside double bounce via full 3D vector GO-PO
%        (doubleBounceGOPO), validated against 0-deg and 45-deg canted
%        90-deg PEC dihedrals.
%     3. OPTIONAL leading-edge fringe term (IncludeEdgeDiffraction, default
%        false; uncalibrated heuristic projected onto edge tangent).

    p = inputParser;
    addParameter(p, 'IncludeMultiBounce', true, @islogical);
    addParameter(p, 'IncludeEdgeDiffraction', false, @islogical);
    addParameter(p, 'Polarization', 'Mean', @(x) ischar(x) || isstring(x));
    addParameter(p, 'PrecomputedConcavePairs', [], @(x) isstruct(x) || isempty(x));
    parse(p, varargin{:});
    opts = p.Results;
    polTag = normalizePolarizationTag(opts.Polarization);

    c0 = 299792458.0; % Speed of light [m/s]
    freqGHz = freqGHz(:).';
    nF = numel(freqGHz);
    lambdaVec = c0 ./ (freqGHz .* 1e9);
    kVec = 2.0 * pi ./ lambdaVec;

    % Support 1D vector of (azDeg, elDeg) directions
    azRad = deg2rad(azDeg(:));
    elRad = deg2rad(elDeg(:));
    if isscalar(azRad) && ~isscalar(elRad)
        azRad = repmat(azRad, size(elRad));
    elseif isscalar(elRad) && ~isscalar(azRad)
        elRad = repmat(elRad, size(azRad));
    end
    nDirs = numel(azRad);

    % Unit vector u_hat pointing from vehicle origin to radar & singularity-free BSA basis [nDirs x 3]
    uRadar = [-cos(elRad) .* cos(azRad), ...
               cos(elRad) .* sin(azRad), ...
               sin(elRad)];
    eH = [sin(azRad), cos(azRad), zeros(nDirs, 1)];
    eV = cross(eH, uRadar, 2);

    % Extract watertight mesh geometry
    TR = body.MeshData.TR;
    V = TR.Points;
    F = TR.ConnectivityList;
    patchID = body.MeshData.FacePatchID;

    v1 = V(F(:, 1), :);
    v2 = V(F(:, 2), :);
    v3 = V(F(:, 3), :);
    crossVec = 0.5 .* cross(v2 - v1, v3 - v1, 2);
    faceAreas = sqrt(sum(crossVec.^2, 2)) + 1e-15;
    faceNormals = crossVec ./ faceAreas; % Outward unit normals [nFaces x 3]

    isUpper = (patchID == 1);
    isLower = (patchID == 2);

    % Mutually visible concave super-cell pairs on the lower surface
    if opts.IncludeMultiBounce
        if isempty(opts.PrecomputedConcavePairs)
            concavePairs = waverider.rcs.buildConcavePairList(body);
        else
            concavePairs = opts.PrecomputedConcavePairs;
        end
    else
        concavePairs = struct('numPairs', 0);
    end

    % Leading-edge segments for the optional fringe heuristic
    if opts.IncludeEdgeDiffraction
        P_le = [body.SurfaceData.X_low(1, :).', body.SurfaceData.Y_low(1, :).', body.SurfaceData.Z_low(1, :).'];
        leMid = 0.5 .* (P_le(1:end-1, :) + P_le(2:end, :));
        dLE = P_le(2:end, :) - P_le(1:end-1, :);
        leLen = sqrt(sum(dLE.^2, 2)) + 1e-12;
        leTan = dLE ./ leLen;
    end

    % Allocate primary channel outputs [nDirs x nF]
    sigmaTotal_m2  = zeros(nDirs, nF);
    sigmaSingle_m2 = zeros(nDirs, nF);
    sigmaMB_m2     = zeros(nDirs, nF);
    sigmaUpper_m2  = zeros(nDirs, nF);
    sigmaLower_m2  = zeros(nDirs, nF);
    sigmaEdge_m2   = zeros(nDirs, nF);

    % Allocate full polarimetric complex scattering matrices [nDirs x nF]
    S_Mean_all = zeros(nDirs, nF);
    S_HH_all   = zeros(nDirs, nF);
    S_VV_all   = zeros(nDirs, nF);
    S_HV_all   = zeros(nDirs, nF);
    S_VH_all   = zeros(nDirs, nF);
    S_PP_all   = zeros(nDirs, nF);
    S_OP_all   = zeros(nDirs, nF);
    S_LL_all   = zeros(nDirs, nF);

    % Process directions in blocks to bound memory ([nFaces x block] complex)
    blk = max(1, floor(4e6 / size(F, 1)));

    for fIdx = 1:nF
        k0 = kVec(fIdx);
        S_single_all   = zeros(nDirs, 1);
        S_single_upper = zeros(nDirs, 1);
        S_single_lower = zeros(nDirs, 1);
        for d0 = 1:blk:nDirs
            dI = d0:min(d0 + blk - 1, nDirs);
            U = uRadar(dI, :);
            cosTheta = faceNormals * U.';                       % [nFaces x nb]
            lit = cosTheta > 0;
            Iface = waverider.rcs.facetPhaseIntegral(v1, v2, v3, faceNormals, 2.0 * k0 .* U);
            facetComplex = (1i * k0 / sqrt(pi)) .* (cosTheta .* lit) .* Iface;
            S_single_all(dI)   = sum(facetComplex, 1).';
            S_single_upper(dI) = sum(facetComplex(isUpper, :), 1).';
            S_single_lower(dI) = sum(facetComplex(isLower, :), 1).';
        end

        % Optional uncalibrated leading-edge fringe heuristic
        S_edge_mean = zeros(nDirs, 1);
        S_edge_HH   = zeros(nDirs, 1);
        S_edge_VV   = zeros(nDirs, 1);
        S_edge_HV   = zeros(nDirs, 1);
        if opts.IncludeEdgeDiffraction
            tDotU = leTan * uRadar.';
            sinBetaSq = max(1.0 - tDotU.^2, 1e-4);
            sinBetaEdge = sqrt(sinBetaSq);
            phiEdge = k0 .* leLen .* tDotU;
            sincEdge = ones(size(phiEdge));
            nz = abs(phiEdge) > 1e-12;
            sincEdge(nz) = sin(phiEdge(nz)) ./ phiEdge(nz);
            rLeDotU = leMid * uRadar.';
            edgeElem = (0.22 ./ sqrt(pi) .* sinBetaEdge .* leLen) .* sincEdge .* exp(1i * 2.0 * k0 .* rLeDotU);
            S_edge_mean = sum(edgeElem, 1).';
            tDotH = leTan * eH.';
            tDotV = leTan * eV.';
            S_edge_HH = sum(edgeElem .* ((tDotH.^2) ./ sinBetaSq), 1).';
            S_edge_VV = sum(edgeElem .* ((tDotV.^2) ./ sinBetaSq), 1).';
            S_edge_HV = sum(edgeElem .* ((tDotH .* tDotV) ./ sinBetaSq), 1).';
        end

        % GO-PO concave double bounce (scalar + full 3D vector Sinclair dyadic)
        S_mb_mean = zeros(nDirs, 1);
        S_mb_pol = struct('HH', zeros(nDirs, 1), 'VV', zeros(nDirs, 1), ...
                          'HV', zeros(nDirs, 1), 'VH', zeros(nDirs, 1));
        if opts.IncludeMultiBounce && concavePairs.numPairs > 0
            [S_mb_mean, ~, S_mb_pol] = waverider.rcs.doubleBounceGOPO(concavePairs, uRadar, k0, false, eH);
        end

        % Exact complex scattering amplitudes for all polarizations
        S_Mean_f = S_single_all + S_edge_mean + S_mb_mean;
        S_HH_f   = S_single_all + S_edge_HH   + S_mb_pol.HH;
        S_VV_f   = S_single_all + S_edge_VV   + S_mb_pol.VV;
        S_HV_f   =                S_edge_HV   + S_mb_pol.HV;
        S_VH_f   =                S_edge_HV   + S_mb_pol.VH;
        S_PP_f   = 0.5 .* (S_HH_f + S_VV_f) + 0.5i .* (S_VH_f - S_HV_f);
        S_OP_f   = 0.5 .* (S_HH_f - S_VV_f) - 0.5i .* (S_HV_f + S_VH_f);
        S_LL_f   = 0.5 .* (S_HH_f - S_VV_f) + 0.5i .* (S_HV_f + S_VH_f);

        S_Mean_all(:, fIdx) = S_Mean_f;
        S_HH_all(:, fIdx)   = S_HH_f;
        S_VV_all(:, fIdx)   = S_VV_f;
        S_HV_all(:, fIdx)   = S_HV_f;
        S_VH_all(:, fIdx)   = S_VH_f;
        S_PP_all(:, fIdx)   = S_PP_f;
        S_OP_all(:, fIdx)   = S_OP_f;
        S_LL_all(:, fIdx)   = S_LL_f;

        % Select active polarization channel decomposition
        switch polTag
            case 'Mean'
                S_sng_ch = S_single_all;   S_up_ch = S_single_upper; S_low_ch = S_single_lower;
                S_edg_ch = S_edge_mean;    S_mb_ch = S_mb_mean;
            case 'HH'
                S_sng_ch = S_single_all;   S_up_ch = S_single_upper; S_low_ch = S_single_lower;
                S_edg_ch = S_edge_HH;      S_mb_ch = S_mb_pol.HH;
            case 'VV'
                S_sng_ch = S_single_all;   S_up_ch = S_single_upper; S_low_ch = S_single_lower;
                S_edg_ch = S_edge_VV;      S_mb_ch = S_mb_pol.VV;
            case {'HV', 'VH'}
                S_sng_ch = zeros(nDirs, 1); S_up_ch = zeros(nDirs, 1); S_low_ch = zeros(nDirs, 1);
                S_edg_ch = S_edge_HV;       S_mb_ch = S_mb_pol.HV;
            case 'PP'
                S_sng_ch = S_single_all;   S_up_ch = S_single_upper; S_low_ch = S_single_lower;
                S_edg_ch = 0.5 .* (S_edge_HH + S_edge_VV);
                S_mb_ch  = 0.5 .* (S_mb_pol.HH + S_mb_pol.VV) + 0.5i .* (S_mb_pol.VH - S_mb_pol.HV);
            case 'OP'
                S_sng_ch = zeros(nDirs, 1); S_up_ch = zeros(nDirs, 1); S_low_ch = zeros(nDirs, 1);
                S_edg_ch = 0.5 .* (S_edge_HH - S_edge_VV) - 1i .* S_edge_HV;
                S_mb_ch  = 0.5 .* (S_mb_pol.HH - S_mb_pol.VV) - 0.5i .* (S_mb_pol.HV + S_mb_pol.VH);
            case 'LL'
                S_sng_ch = zeros(nDirs, 1); S_up_ch = zeros(nDirs, 1); S_low_ch = zeros(nDirs, 1);
                S_edg_ch = 0.5 .* (S_edge_HH - S_edge_VV) + 1i .* S_edge_HV;
                S_mb_ch  = 0.5 .* (S_mb_pol.HH - S_mb_pol.VV) + 0.5i .* (S_mb_pol.HV + S_mb_pol.VH);
        end

        % Small incoherent display floor (-50 dBsm) to avoid -Inf in nulls
        noiseFloor_m2 = 1e-5;
        sigmaSingle_m2(:, fIdx) = abs(S_sng_ch + S_edg_ch).^2 + noiseFloor_m2;
        sigmaMB_m2(:, fIdx)     = abs(S_mb_ch).^2;
        sigmaTotal_m2(:, fIdx)  = abs(S_sng_ch + S_edg_ch + S_mb_ch).^2 + noiseFloor_m2;
        sigmaUpper_m2(:, fIdx)  = abs(S_up_ch).^2 + 0.5 * noiseFloor_m2;
        sigmaLower_m2(:, fIdx)  = abs(S_low_ch + S_mb_ch).^2 + 0.5 * noiseFloor_m2;
        sigmaEdge_m2(:, fIdx)   = abs(S_edg_ch).^2;
    end

    noiseFloor_m2 = 1e-5;
    toDB = @(m2) 10.0 .* log10(max(m2, 1e-10));
    span_m2 = abs(S_HH_all).^2 + abs(S_VV_all).^2 + abs(S_HV_all).^2 + abs(S_VH_all).^2;

    polStruct = struct();
    polStruct.S_Mean    = S_Mean_all;
    polStruct.S_HH      = S_HH_all;
    polStruct.S_VV      = S_VV_all;
    polStruct.S_HV      = S_HV_all;
    polStruct.S_VH      = S_VH_all;
    polStruct.S_PP      = S_PP_all;
    polStruct.S_OP      = S_OP_all;
    polStruct.S_LL      = S_LL_all;
    polStruct.Mean_dBsm = toDB(abs(S_Mean_all).^2 + noiseFloor_m2);
    polStruct.HH_dBsm   = toDB(abs(S_HH_all).^2   + noiseFloor_m2);
    polStruct.VV_dBsm   = toDB(abs(S_VV_all).^2   + noiseFloor_m2);
    polStruct.HV_dBsm   = toDB(abs(S_HV_all).^2   + noiseFloor_m2);
    polStruct.VH_dBsm   = toDB(abs(S_VH_all).^2   + noiseFloor_m2);
    polStruct.PP_dBsm   = toDB(abs(S_PP_all).^2   + noiseFloor_m2);
    polStruct.OP_dBsm   = toDB(abs(S_OP_all).^2   + noiseFloor_m2);
    polStruct.LL_dBsm   = toDB(abs(S_LL_all).^2   + noiseFloor_m2);
    polStruct.Span_m2   = span_m2;
    polStruct.Span_dBsm = toDB(span_m2 + noiseFloor_m2);

    rcsData = struct();
    rcsData.FreqGHz = freqGHz;
    rcsData.Lambda_m = lambdaVec;
    rcsData.AzDeg = azDeg(:);
    rcsData.ElDeg = elDeg(:);
    rcsData.Polarization = polTag;
    rcsData.SigmaTotal_m2 = sigmaTotal_m2;
    rcsData.SigmaTotal_dBsm = toDB(sigmaTotal_m2);
    rcsData.SigmaSingle_dBsm = toDB(sigmaSingle_m2);
    rcsData.SigmaMultiBounce_dBsm = toDB(sigmaMB_m2);
    rcsData.SigmaUpper_dBsm = toDB(sigmaUpper_m2);
    rcsData.SigmaLower_dBsm = toDB(sigmaLower_m2);
    rcsData.SigmaEdge_dBsm = toDB(sigmaEdge_m2);
    rcsData.MultiBounceGain_dB = rcsData.SigmaTotal_dBsm - rcsData.SigmaSingle_dBsm;
    rcsData.Polarimetric = polStruct;
    rcsData.ConcavePairs = concavePairs;
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
            error('computeRCS:Polarization', ...
                'Unsupported Polarization "%s". Choose Mean, HH, VV, HV, VH, PP (RL/LR), OP (RR), or LL.', ...
                char(polIn));
    end
end
