function rcsData = computeRCS(body, freqGHz, azDeg, elDeg, varargin)
%COMPUTERCS Monostatic RCS: exact facet Physical Optics + GO-PO concave double bounce.
%   rcsData = COMPUTERCS(body, freqGHz, azDeg, elDeg, Name, Value) computes
%   the monostatic Radar Cross Section (RCS) of a watertight waverider body
%   (UHF 0.3 GHz to X-band 12 GHz).
%
%   Coordinate & Angle Convention:
%     - Vehicle nose points toward -X (nose at X=0, base at X=L).
%     - Unit vector FROM target TO radar:
%         u_radar = [-cos(el)*cos(az), cos(el)*sin(az), sin(el)]
%         az = 0 deg  -> nose-on,  az = 90 deg -> starboard broadside,
%         az = 180    -> tail-on,  el > 0 dorsal, el < 0 ventral (look-up).
%
%   Scattering mechanisms:
%     1. Single-bounce Physical Optics over every lit facet using the EXACT
%        closed-form triangle integral (Gordon 1975; facetPhaseIntegral):
%            S_PO = (i k / sqrt(pi)) sum_f (n_f . u) \int_f exp(i 2k u . r) dA
%        sigma = |S|^2 reproduces 4 pi A^2 / lambda^2 for a flat plate and
%        pi a^2 for a sphere (tests/test_physics_validation.m).
%        Lit test is n . u > 0 only (no self-shadowing between parts).
%     2. Concave-underside double bounce via GO-PO (doubleBounceGOPO),
%        validated against the closed-form 90-deg dihedral.
%     3. OPTIONAL leading-edge fringe term (IncludeEdgeDiffraction, default
%        false). This is an UNCALIBRATED heuristic, not a true MEC/PTD
%        solution; PO already contains the physical-optics edge contribution.
%
%   Scalar (polarisation-averaged) model: no HH/VV distinction, no creeping
%   waves, no travelling waves, no shadowing/blockage, PEC surfaces. Results
%   in the UHF resonance regime (body only a few wavelengths) are indicative.

    p = inputParser;
    addParameter(p, 'IncludeMultiBounce', true, @islogical);
    addParameter(p, 'IncludeEdgeDiffraction', false, @islogical);
    addParameter(p, 'Polarization', 'Mean', @(x) ischar(x) || isstring(x));
    addParameter(p, 'PrecomputedConcavePairs', [], @(x) isstruct(x) || isempty(x));
    parse(p, varargin{:});
    opts = p.Results;

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

    % Unit vector u_hat pointing from vehicle origin to radar [nDirs x 3]
    uRadar = [-cos(elRad) .* cos(azRad), ...
               cos(elRad) .* sin(azRad), ...
               sin(elRad)];

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

    % Allocate outputs [nDirs x nF]
    sigmaTotal_m2  = zeros(nDirs, nF);
    sigmaSingle_m2 = zeros(nDirs, nF);
    sigmaMB_m2     = zeros(nDirs, nF);
    sigmaUpper_m2  = zeros(nDirs, nF);
    sigmaLower_m2  = zeros(nDirs, nF);
    sigmaEdge_m2   = zeros(nDirs, nF);

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
        S_edge = zeros(nDirs, 1);
        if opts.IncludeEdgeDiffraction
            tDotU = leTan * uRadar.';
            sinBetaEdge = sqrt(max(1.0 - tDotU.^2, 1e-4));
            phiEdge = k0 .* leLen .* tDotU;
            sincEdge = ones(size(phiEdge));
            nz = abs(phiEdge) > 1e-12;
            sincEdge(nz) = sin(phiEdge(nz)) ./ phiEdge(nz);
            rLeDotU = leMid * uRadar.';
            edgeDiffCoeff = 0.22 ./ sqrt(pi) .* sinBetaEdge;
            S_edge = sum((edgeDiffCoeff .* leLen) .* sincEdge .* exp(1i * 2.0 * k0 .* rLeDotU), 1).';
        end

        % GO-PO concave double bounce
        S_mb = zeros(nDirs, 1);
        if opts.IncludeMultiBounce && concavePairs.numPairs > 0
            S_mb = waverider.rcs.doubleBounceGOPO(concavePairs, uRadar, k0);
        end

        % Small incoherent display floor (-50 dBsm) to avoid -Inf in nulls
        noiseFloor_m2 = 1e-5;
        sigmaSingle_m2(:, fIdx) = abs(S_single_all + S_edge).^2 + noiseFloor_m2;
        sigmaMB_m2(:, fIdx)     = abs(S_mb).^2;
        sigmaTotal_m2(:, fIdx)  = abs(S_single_all + S_edge + S_mb).^2 + noiseFloor_m2;
        sigmaUpper_m2(:, fIdx)  = abs(S_single_upper).^2 + 0.5 * noiseFloor_m2;
        sigmaLower_m2(:, fIdx)  = abs(S_single_lower + S_mb).^2 + 0.5 * noiseFloor_m2;
        sigmaEdge_m2(:, fIdx)   = abs(S_edge).^2;
    end

    rcsData = struct();
    rcsData.FreqGHz = freqGHz;
    rcsData.Lambda_m = lambdaVec;
    rcsData.AzDeg = azDeg(:);
    rcsData.ElDeg = elDeg(:);
    rcsData.SigmaTotal_m2 = sigmaTotal_m2;
    rcsData.SigmaTotal_dBsm = 10.0 .* log10(max(sigmaTotal_m2, 1e-10));
    rcsData.SigmaSingle_dBsm = 10.0 .* log10(max(sigmaSingle_m2, 1e-10));
    rcsData.SigmaMultiBounce_dBsm = 10.0 .* log10(max(sigmaMB_m2, 1e-10));
    rcsData.SigmaUpper_dBsm = 10.0 .* log10(max(sigmaUpper_m2, 1e-10));
    rcsData.SigmaLower_dBsm = 10.0 .* log10(max(sigmaLower_m2, 1e-10));
    rcsData.SigmaEdge_dBsm = 10.0 .* log10(max(sigmaEdge_m2, 1e-10));
    rcsData.MultiBounceGain_dB = rcsData.SigmaTotal_dBsm - rcsData.SigmaSingle_dBsm;
    rcsData.ConcavePairs = concavePairs;
end
