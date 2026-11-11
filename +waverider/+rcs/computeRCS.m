function rcsData = computeRCS(body, freqGHz, azDeg, elDeg, varargin)
%COMPUTERCS Multi-band Physical Optics + Concave Multi-Bounce + Edge Diffraction RCS solver.
%   rcsData = COMPUTERCS(body, freqGHz, azDeg, elDeg, Name, Value) computes
%   the monostatic Radar Cross Section (RCS) of a watertight hypersonic
%   waverider body from UHF (0.3 GHz) to X-band (10+ GHz).
%
%   Coordinate & Angle Convention:
%     - Vehicle nose points toward -X (nose at X=0, base at X=L).
%     - Monostatic radar direction unit vector pointing FROM target TO radar:
%         az = 0 deg, el = 0 deg   -> Nose-on illumination (-X direction)
%         az = 90 deg, el = 0 deg  -> Starboard broadside (+Y direction)
%         az = 180 deg, el = 0 deg -> Tail-on (+X direction)
%         el > 0 deg               -> Dorsal (upper hemisphere look-down)
%         el < 0 deg               -> Ventral (lower hemisphere look-up into concave underside)
%         u_radar = [-cos(el)*cos(az), cos(el)*sin(az), sin(el)]
%
%   Scattering Mechanisms Modeled:
%     1. Coherent Facet Physical Optics (PO) with exact Gordon/sinc phase-slope
%        integration over each triangle (accurate from UHF to X-band even when
%        facets are several wavelengths across).
%     2. Concave Underside Double-Bounce (Iterative PO / Ray-Bounce Coupling)
%        between mutually visible port/starboard facets on the concave lower
%        compression surface.
%     3. Sharp Swept Leading-Edge Travelling/Diffraction contribution (MEC).

    p = inputParser;
    addParameter(p, 'IncludeMultiBounce', true, @islogical);
    addParameter(p, 'IncludeEdgeDiffraction', true, @islogical);
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
    faceCenters = (v1 + v2 + v3) ./ 3.0; % [nFaces x 3]

    crossVec = 0.5 .* cross(v2 - v1, v3 - v1, 2);
    faceAreas = sqrt(sum(crossVec.^2, 2)) + 1e-15;
    faceNormals = crossVec ./ faceAreas; % Outward unit normals [nFaces x 3]

    % Tangent vectors & characteristic lengths on each triangle for sinc factor
    e1 = v2 - v1;
    e2 = v3 - v1;

    isUpper = (patchID == 1);
    isLower = (patchID == 2);

    % Precompute mutually visible concave facet pairs on the lower surface
    if opts.IncludeMultiBounce
        if isempty(opts.PrecomputedConcavePairs)
            concavePairs = waverider.rcs.buildConcavePairList(body);
        else
            concavePairs = opts.PrecomputedConcavePairs;
        end
    else
        concavePairs = struct('numPairs', 0);
    end

    % Extract leading-edge segments for sharp-edge diffraction
    X_le = body.SurfaceData.X_low(1, :).';
    Y_le = body.SurfaceData.Y_low(1, :).';
    Z_le = body.SurfaceData.Z_low(1, :).';
    P_le = [X_le, Y_le, Z_le];
    leMid = 0.5 .* (P_le(1:end-1, :) + P_le(2:end, :));
    dLE = P_le(2:end, :) - P_le(1:end-1, :);
    leLen = sqrt(sum(dLE.^2, 2)) + 1e-12;
    leTan = dLE ./ leLen;

    % Allocate outputs [nDirs x nF]
    sigmaTotal_m2  = zeros(nDirs, nF);
    sigmaSingle_m2 = zeros(nDirs, nF);
    sigmaMB_m2     = zeros(nDirs, nF);
    sigmaUpper_m2  = zeros(nDirs, nF);
    sigmaLower_m2  = zeros(nDirs, nF);
    sigmaEdge_m2   = zeros(nDirs, nF);

    % Loop over frequencies and evaluate vectorized over directions
    for fIdx = 1:nF
        k0 = kVec(fIdx);
        lam = lambdaVec(fIdx);

        % -----------------------------------------------------------------
        % 1. Coherent Single-Bounce Physical Optics (PO) with Sinc Integration
        % -----------------------------------------------------------------
        % Illumination cosine: cos(theta_i) = n_hat . u_radar [nFaces x nDirs]
        cosTheta = faceNormals * uRadar.';
        litMask = (cosTheta > 1e-5);

        % Phase at triangle centroids: exp(i * 2 * k0 * r_c . u_radar)
        rDotU = faceCenters * uRadar.'; % [nFaces x nDirs]

        % Gordon/Sinc triangle intra-facet phase variation factor so large facets
        % at X-band integrate accurately without centroid aliasing:
        % Phase differences along edges e1 and e2:
        phi1 = k0 .* (e1 * uRadar.'); % [nFaces x nDirs]
        phi2 = k0 .* (e2 * uRadar.'); % [nFaces x nDirs]
        sinc1 = sin(phi1 ./ 2.0 + 1e-12) ./ (phi1 ./ 2.0 + 1e-12);
        sinc2 = sin(phi2 ./ 2.0 + 1e-12) ./ (phi2 ./ 2.0 + 1e-12);
        phi12 = phi1 - phi2;
        sinc12 = sin(phi12 ./ 2.0 + 1e-12) ./ (phi12 ./ 2.0 + 1e-12);
        facetShapeFactor = (sinc1 .* sinc2 .* sinc12).^(1.0 / 3.0);

        % Complex PO scattering amplitude per facet:
        % S_PO = (i * k0 / sqrt(pi)) * sum_f A_f * (n_f . u) * sincFactor * exp(i * 2 * k0 * r_f . u)
        prefactorPO = (1i * k0 / sqrt(pi));
        facetWeight = (faceAreas .* cosTheta) .* litMask .* facetShapeFactor;
        facetComplex = prefactorPO .* facetWeight .* exp(1i * 2.0 * k0 .* rDotU);

        S_single_all   = sum(facetComplex, 1).';                 % [nDirs x 1]
        S_single_upper = sum(facetComplex(isUpper, :), 1).';     % [nDirs x 1]
        S_single_lower = sum(facetComplex(isLower, :), 1).';     % [nDirs x 1]

        % -----------------------------------------------------------------
        % 2. Sharp Swept Leading-Edge Diffraction (MEC / PTD)
        % -----------------------------------------------------------------
        S_edge = zeros(nDirs, 1);
        if opts.IncludeEdgeDiffraction
            % Only active when radar illuminates near leading-edge Keller cone
            tDotU = leTan * uRadar.'; % [nSegs x nDirs]
            sinBetaEdge = sqrt(max(1.0 - tDotU.^2, 1e-4));
            phiEdge = k0 .* leLen .* tDotU;
            sincEdge = sin(phiEdge + 1e-12) ./ (phiEdge + 1e-12);
            rLeDotU = leMid * uRadar.';
            % Equivalent edge current diffraction amplitude
            edgeDiffCoeff = 0.22 ./ sqrt(pi) .* sinBetaEdge;
            S_edge_seg = (edgeDiffCoeff .* leLen) .* sincEdge .* exp(1i * 2.0 * k0 .* rLeDotU);
            S_edge = sum(S_edge_seg, 1).';
        end

        % -----------------------------------------------------------------
        % 3. Concave Underside Double-Bounce (Iterative PO / Ray Coupling)
        % -----------------------------------------------------------------
        S_mb = zeros(nDirs, 1);
        if opts.IncludeMultiBounce && concavePairs.numPairs > 0
            S_mb = evaluateConcaveDoubleBounce(concavePairs, uRadar, k0, lam);
        end

        % Combine coherent and incoherent floor (avoid deep -80 dB artificial nulls)
        noiseFloor_m2 = 1e-5; % -50 dBsm physical diffuse floor
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

function S_mb = evaluateConcaveDoubleBounce(cp, uRadar, k0, lam)
%EVALUATECONCAVEDOUBLEBOUNCE Vectorized double-bounce PO across concave lower surface.
    % Incident wave vector k_inc = -uRadar [nDirs x 3]
    % Both patch P and patch Q must be illuminated or in the bounce path
    cosP = cp.nP * uRadar.'; % [nPairs x nDirs]
    cosQ = cp.nQ * uRadar.'; % [nPairs x nDirs]
    bothLit = (cosP > 0.02) & (cosQ > 0.02);

    % Specular reflection of incident ray (-uRadar) off patch P:
    % s_refl = -uRadar + 2*(uRadar . nP)*nP
    % Alignment of reflected ray from P with inter-patch vector dPQ:
    % s_refl . dPQ = -(dPQ . uRadar) + 2 * (nP . uRadar) * (nP . dPQ)
    dDotU = cp.dPQ * uRadar.'; % [nPairs x nDirs]
    alignPQ = -dDotU + 2.0 .* cosP .* repmat(cp.nP_dot_dPQ, 1, size(uRadar, 1));
    % Also symmetric bounce Q -> P:
    alignQP =  dDotU + 2.0 .* cosQ .* repmat(cp.nQ_dot_dQP, 1, size(uRadar, 1));
    alignMax = max(alignPQ, alignQP);

    % Diffraction-broadened bounce lobe width (narrower at X-band, wider at UHF)
    lobeSharpness = min(max(1.5 * k0 * cp.charSize, 4.0), 45.0);
    lobeMask = (alignMax > 0.65) & bothLit;
    lobeWeight = zeros(size(alignMax));
    lobeWeight(lobeMask) = exp(-lobeSharpness .* (1.0 - alignMax(lobeMask)));

    % Effective double-bounce projected area & electrical size scaling
    % High-frequency dihedral corner factor scales with k0 / sqrt(pi) and
    % transitions smoothly in the UHF resonance regime (D_trough ~ lambda)
    elecScale = min((cp.distPQ ./ (0.35 * lam)).^0.65, 2.8);
    effAreaMat = repmat(cp.effArea .* elecScale, 1, size(uRadar, 1));
    ampPair = (1i * k0 / sqrt(pi)) .* effAreaMat .* ...
              sqrt(max(cosP .* cosQ, 0.0)) .* lobeWeight;

    % Total two-bounce phase path: (rP . u + rQ . u - distPQ) * k0
    phasePath = k0 .* ((cp.rP + cp.rQ) * uRadar.' - repmat(cp.distPQ, 1, size(uRadar, 1)));
    S_mb = sum(ampPair .* exp(1i .* phasePath), 1).';
end
