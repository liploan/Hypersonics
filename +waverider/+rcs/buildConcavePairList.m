function cp = buildConcavePairList(body)
%BUILDCONCAVEPAIRLIST Identify mutually visible concave facet pairs on the lower surface.
%   cp = BUILDCONCAVEPAIRLIST(body) coarse-grains the lower compression
%   surface into super-facets and finds all pairs (P, Q) whose outward
%   normals face partially toward one another across the concave underside:
%     (rQ - rP) . nP > 0  AND  (rP - rQ) . nQ > 0

    surfData = body.SurfaceData;
    if isfield(surfData, 'IsAxisymmetric') && surfData.IsAxisymmetric
        cp = buildAxisymmetricConcavePairs(surfData);
        return;
    end

    NxWave = surfData.NxWaverider;
    X = surfData.X_low(1:NxWave, :);
    Y = surfData.Y_low(1:NxWave, :);
    Z = surfData.Z_low(1:NxWave, :);

    % Subsample structured lower surface into ~14 x 21 super-cells for fast pair coupling
    nI = min(14, size(X, 1));
    nJ = min(21, size(X, 2));
    iIdx = unique(round(linspace(1, size(X, 1), nI)));
    jIdx = unique(round(linspace(1, size(X, 2), nJ)));

    Xs = X(iIdx, jIdx);
    Ys = Y(iIdx, jIdx);
    Zs = Z(iIdx, jIdx);

    [Nr, Nc] = size(Xs);
    maxCells = (Nr - 1) * (Nc - 1);
    rCell = zeros(maxCells, 3);
    nCell = zeros(maxCells, 3);
    aCell = zeros(maxCells, 1);
    vCell = zeros(maxCells, 4, 3);   % planar quad polygon (CCW about outward normal)
    nCells = 0;

    for r = 1:(Nr - 1)
        for c = 1:(Nc - 1)
            p1 = [Xs(r, c),     Ys(r, c),     Zs(r, c)];
            p2 = [Xs(r+1, c),   Ys(r+1, c),   Zs(r+1, c)];
            p3 = [Xs(r+1, c+1), Ys(r+1, c+1), Zs(r+1, c+1)];
            p4 = [Xs(r, c+1),   Ys(r, c+1),   Zs(r, c+1)];

            % Outward normal of lower surface points downward (-Z):
            % (p2 - p1) is +X, (p4 - p1) is +Y -> cross(p4 - p1, p2 - p1) is -Z
            cr1 = cross(p3 - p1, p2 - p1);
            cr2 = cross(p4 - p1, p3 - p1);
            nVec = 0.5 .* (cr1 + cr2);
            areaVal = norm(nVec);
            if areaVal > 1e-6
                nCells = nCells + 1;
                ctr = 0.25 .* (p1 + p2 + p3 + p4);
                nHat = nVec ./ areaVal;
                rCell(nCells, :) = ctr;
                nCell(nCells, :) = nHat;
                aCell(nCells) = areaVal;
                % Project the (generally non-planar) quad onto its mean plane,
                % ordered CCW about the outward normal: p1 -> p4 -> p3 -> p2
                quad = [p1; p4; p3; p2];
                quad = quad - ((quad - ctr) * nHat.') .* nHat;
                vCell(nCells, :, :) = reshape(quad, 1, 4, 3);
            end
        end
    end

    rCell = rCell(1:nCells, :);
    nCell = nCell(1:nCells, :);
    aCell = aCell(1:nCells);
    vCell = vCell(1:nCells, :, :);

    % Find mutually visible concave pairs (P < Q)
    maxPairs = nCells * (nCells - 1) / 2;
    pList = zeros(maxPairs, 1);
    qList = zeros(maxPairs, 1);
    nP_dot_d = zeros(maxPairs, 1);
    nQ_dot_d = zeros(maxPairs, 1);
    distVec = zeros(maxPairs, 1);
    dVec = zeros(maxPairs, 3);
    nPairs = 0;

    for p = 1:(nCells - 1)
        rP = rCell(p, :);
        nP = nCell(p, :);
        for q = (p + 1):nCells
            vecPQ = rCell(q, :) - rP;
            dist = norm(vecPQ);
            if dist < 1e-4
                continue;
            end
            dHat = vecPQ ./ dist;
            % Mutual concavity test: patch P sees patch Q and patch Q sees patch P
            cP = dot(nP, dHat);
            cQ = dot(nCell(q, :), -dHat);
            if (cP > 0.015) && (cQ > 0.015)
                nPairs = nPairs + 1;
                pList(nPairs) = p;
                qList(nPairs) = q;
                nP_dot_d(nPairs) = cP;
                nQ_dot_d(nPairs) = cQ;
                distVec(nPairs) = dist;
                dVec(nPairs, :) = dHat;
            end
        end
    end

    % Keep ALL mutually visible pairs: the GO-PO double-bounce evaluator culls
    % per look direction using the exact specular-reflection geometry, so no
    % ad-hoc truncation by a heuristic "pair strength" is needed.
    if nPairs > 0
        pList = pList(1:nPairs);
        qList = qList(1:nPairs);
        nP_dot_d = nP_dot_d(1:nPairs);
        nQ_dot_d = nQ_dot_d(1:nPairs);
        distVec = distVec(1:nPairs);
        dVec = dVec(1:nPairs, :);
    end

    cp = struct();
    cp.numPairs = nPairs;
    % Cell (super-facet) geometry used by the GO-PO double-bounce evaluator
    cp.cellCenter = rCell;
    cp.cellNormal = nCell;
    cp.cellArea = aCell;
    cp.cellVerts = vCell;                 % [nCells x 4 x 3] planar CCW quads
    cp.cellRadius = zeros(nCells, 1);
    for c = 1:nCells
        cp.cellRadius(c) = max(vecnorm(squeeze(vCell(c, :, :)) - rCell(c, :), 2, 2));
    end
    if nPairs > 0
        cp.pIdx = pList;
        cp.qIdx = qList;
        cp.rP = rCell(pList, :);
        cp.rQ = rCell(qList, :);
        cp.nP = nCell(pList, :);
        cp.nQ = nCell(qList, :);
        cp.aP = aCell(pList);
        cp.aQ = aCell(qList);
        cp.dPQ = dVec;
        cp.distPQ = distVec;
        cp.nP_dot_dPQ = nP_dot_d;
        cp.nQ_dot_dQP = nQ_dot_d;
        cp.effArea = sqrt(cp.aP .* cp.aQ) .* sqrt(nP_dot_d .* nQ_dot_d);
    else
        cp.pIdx = zeros(0, 1); cp.qIdx = zeros(0, 1);
        cp.rP = zeros(0, 3); cp.rQ = zeros(0, 3);
        cp.nP = zeros(0, 3); cp.nQ = zeros(0, 3);
        cp.aP = zeros(0, 1); cp.aQ = zeros(0, 1);
        cp.dPQ = zeros(0, 3); cp.distPQ = zeros(0, 1);
        cp.nP_dot_dPQ = zeros(0, 1); cp.nQ_dot_dQP = zeros(0, 1);
        cp.effArea = zeros(0, 1);
    end
end

function cp = buildAxisymmetricConcavePairs(surfData)
%BUILDAXISYMMETRICCONCAVEPAIRS Build GO-PO super-cells for RV aft cavities & flared conic corners.
    hasConcaveAft = (~isfield(surfData, 'HasAftFrustum') || ~surfData.HasAftFrustum) && ...
                    isfield(surfData, 'AftCavityDepth_m') && (surfData.AftCavityDepth_m > 1e-5);
    hasFlaredCones = isfield(surfData, 'ConeAnglesDeg') && (numel(surfData.ConeAnglesDeg) >= 2) && ...
                     any(diff(surfData.ConeAnglesDeg) > 0.25);

    nPhi = 16; % Multiple of 4 -> exact Y <-> -Y and Z <-> -Z symmetry
    phiEdges = linspace(-pi, pi, nPhi + 1);

    maxCells = 256;
    rCell = zeros(maxCells, 3);
    nCell = zeros(maxCells, 3);
    aCell = zeros(maxCells, 1);
    vCell = zeros(maxCells, 4, 3);
    groupCell = zeros(maxCells, 1); % 1 = recessed aft cavity, 2 = flared exterior skin
    erCell = zeros(maxCells, 3);
    nCells = 0;

    % 1. Recessed Concave Aft Cavity Super-Cells (unobstructed convex interior void)
    if hasConcaveAft
        R_b = surfData.RBase_m;
        L_aft = surfData.XProfile(end);
        dAft = surfData.AftCavityDepth_m;
        pAft = surfData.AftCavityPower;
        nRad = 6;
        uEdges = linspace(1.0, 0.06, nRad + 1);

        for r = 1:nRad
            u1 = uEdges(r);
            u2 = uEdges(r + 1);
            x1 = L_aft - dAft * (1.0 - u1^pAft);
            x2 = L_aft - dAft * (1.0 - u2^pAft);
            for c = 1:nPhi
                ph1 = phiEdges(c);
                ph2 = phiEdges(c + 1);
                phMid = 0.5 * (ph1 + ph2);

                p1 = [x1, R_b * u1 * sin(ph1), R_b * u1 * cos(ph1)];
                p2 = [x2, R_b * u2 * sin(ph1), R_b * u2 * cos(ph1)];
                p3 = [x2, R_b * u2 * sin(ph2), R_b * u2 * cos(ph2)];
                p4 = [x1, R_b * u1 * sin(ph2), R_b * u1 * cos(ph2)];

                nVec = 0.5 .* (cross(p2 - p1, p3 - p1) + cross(p3 - p1, p4 - p1));
                areaVal = norm(nVec);
                if areaVal > 1e-7
                    nCells = nCells + 1;
                    ctr = 0.25 .* (p1 + p2 + p3 + p4);
                    nHat = nVec ./ areaVal;
                    rCell(nCells, :) = ctr;
                    nCell(nCells, :) = nHat;
                    aCell(nCells) = areaVal;
                    groupCell(nCells) = 1;
                    erCell(nCells, :) = [0.0, sin(phMid), cos(phMid)];
                    quad = [p1; p2; p3; p4];
                    quad = quad - ((quad - ctr) * nHat.') .* nHat;
                    vCell(nCells, :, :) = reshape(quad, 1, 4, 3);
                end
            end
        end
    end

    % 2. Flared Multi-Conic Exterior Skin Super-Cells (same-side corner reflections)
    if hasFlaredCones
        xProf = surfData.XProfile(:);
        rProf = surfData.RProfile(:);
        nP = numel(xProf);
        baseIdx = round(linspace(max(2, find(rProf > 1e-5, 1, 'first')), nP, 12));
        breakIdx = [];
        if isfield(surfData, 'StageBreaksX')
            for b = 1:numel(surfData.StageBreaksX)
                [~, idxB] = min(abs(xProf - surfData.StageBreaksX(b)));
                if rProf(idxB) > 1e-5
                    breakIdx(end + 1) = idxB; %#ok<AGROW>
                end
            end
        end
        iIdx = unique([baseIdx(:); breakIdx(:)]);

        for r = 1:(numel(iIdx) - 1)
            x1 = xProf(iIdx(r));     rad1 = rProf(iIdx(r));
            x2 = xProf(iIdx(r + 1)); rad2 = rProf(iIdx(r + 1));
            for c = 1:nPhi
                ph1 = phiEdges(c);
                ph2 = phiEdges(c + 1);
                phMid = 0.5 * (ph1 + ph2);

                p1 = [x1, rad1 * sin(ph1), rad1 * cos(ph1)];
                p2 = [x2, rad2 * sin(ph1), rad2 * cos(ph1)];
                p3 = [x2, rad2 * sin(ph2), rad2 * cos(ph2)];
                p4 = [x1, rad1 * sin(ph2), rad1 * cos(ph2)];

                nVec = 0.5 .* (cross(p2 - p1, p3 - p1) + cross(p3 - p1, p4 - p1));
                areaVal = norm(nVec);
                if areaVal > 1e-7
                    nCells = nCells + 1;
                    ctr = 0.25 .* (p1 + p2 + p3 + p4);
                    nHat = nVec ./ areaVal;
                    rCell(nCells, :) = ctr;
                    nCell(nCells, :) = nHat;
                    aCell(nCells) = areaVal;
                    groupCell(nCells) = 2;
                    erCell(nCells, :) = [0.0, sin(phMid), cos(phMid)];
                    quad = [p1; p2; p3; p4];
                    quad = quad - ((quad - ctr) * nHat.') .* nHat;
                    vCell(nCells, :, :) = reshape(quad, 1, 4, 3);
                end
            end
        end
    end

    rCell = rCell(1:nCells, :);
    nCell = nCell(1:nCells, :);
    aCell = aCell(1:nCells);
    vCell = vCell(1:nCells, :, :);
    groupCell = groupCell(1:nCells);
    erCell = erCell(1:nCells, :);

    maxPairs = max(nCells * (nCells - 1) / 2, 0);
    pList = zeros(maxPairs, 1);
    qList = zeros(maxPairs, 1);
    nP_dot_d = zeros(maxPairs, 1);
    nQ_dot_d = zeros(maxPairs, 1);
    distVec = zeros(maxPairs, 1);
    dVec = zeros(maxPairs, 3);
    nPairs = 0;

    for p = 1:(nCells - 1)
        rP = rCell(p, :);
        nP = nCell(p, :);
        gP = groupCell(p);
        for q = (p + 1):nCells
            if groupCell(q) ~= gP
                continue;
            end
            if (gP == 2) && (dot(erCell(p, :), erCell(q, :)) <= 0.25)
                % Prevent exterior rays from tunneling through the solid RV body
                continue;
            end
            vecPQ = rCell(q, :) - rP;
            dist = norm(vecPQ);
            if dist < 1e-4
                continue;
            end
            dHat = vecPQ ./ dist;
            cP = dot(nP, dHat);
            cQ = dot(nCell(q, :), -dHat);
            if (cP > 0.015) && (cQ > 0.015)
                nPairs = nPairs + 1;
                pList(nPairs) = p;
                qList(nPairs) = q;
                nP_dot_d(nPairs) = cP;
                nQ_dot_d(nPairs) = cQ;
                distVec(nPairs) = dist;
                dVec(nPairs, :) = dHat;
            end
        end
    end

    cp = struct();
    cp.numPairs = nPairs;
    cp.cellCenter = rCell;
    cp.cellNormal = nCell;
    cp.cellArea = aCell;
    cp.cellVerts = vCell;
    cp.cellRadius = zeros(nCells, 1);
    for c = 1:nCells
        cp.cellRadius(c) = max(vecnorm(squeeze(vCell(c, :, :)) - rCell(c, :), 2, 2));
    end
    if nPairs > 0
        pList = pList(1:nPairs);
        qList = qList(1:nPairs);
        nP_dot_d = nP_dot_d(1:nPairs);
        nQ_dot_d = nQ_dot_d(1:nPairs);
        distVec = distVec(1:nPairs);
        dVec = dVec(1:nPairs, :);

        cp.pIdx = pList;
        cp.qIdx = qList;
        cp.rP = rCell(pList, :);
        cp.rQ = rCell(qList, :);
        cp.nP = nCell(pList, :);
        cp.nQ = nCell(qList, :);
        cp.aP = aCell(pList);
        cp.aQ = aCell(qList);
        cp.dPQ = dVec;
        cp.distPQ = distVec;
        cp.nP_dot_dPQ = nP_dot_d;
        cp.nQ_dot_dQP = nQ_dot_d;
        cp.effArea = sqrt(cp.aP .* cp.aQ) .* sqrt(nP_dot_d .* nQ_dot_d);
    else
        cp.pIdx = zeros(0, 1); cp.qIdx = zeros(0, 1);
        cp.rP = zeros(0, 3); cp.rQ = zeros(0, 3);
        cp.nP = zeros(0, 3); cp.nQ = zeros(0, 3);
        cp.aP = zeros(0, 1); cp.aQ = zeros(0, 1);
        cp.dPQ = zeros(0, 3); cp.distPQ = zeros(0, 1);
        cp.nP_dot_dPQ = zeros(0, 1); cp.nQ_dot_dQP = zeros(0, 1);
        cp.effArea = zeros(0, 1);
    end
end
