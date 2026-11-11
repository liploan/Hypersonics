function cp = buildConcavePairList(body)
%BUILDCONCAVEPAIRLIST Identify mutually visible concave facet pairs on the lower surface.
%   cp = BUILDCONCAVEPAIRLIST(body) coarse-grains the lower compression
%   surface into super-facets and finds all pairs (P, Q) whose outward
%   normals face partially toward one another across the concave underside:
%     (rQ - rP) . nP > 0  AND  (rP - rQ) . nQ > 0

    surfData = body.SurfaceData;
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
                rCell(nCells, :) = 0.25 .* (p1 + p2 + p3 + p4);
                nCell(nCells, :) = nVec ./ areaVal;
                aCell(nCells) = areaVal;
            end
        end
    end

    rCell = rCell(1:nCells, :);
    nCell = nCell(1:nCells, :);
    aCell = aCell(1:nCells);

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

    % Keep top 350 strongest concave coupling pairs for speed & accuracy
    if nPairs > 0
        pList = pList(1:nPairs);
        qList = qList(1:nPairs);
        nP_dot_d = nP_dot_d(1:nPairs);
        nQ_dot_d = nQ_dot_d(1:nPairs);
        distVec = distVec(1:nPairs);
        dVec = dVec(1:nPairs, :);

        pairStrength = (aCell(pList) .* aCell(qList)) .* (nP_dot_d .* nQ_dot_d) ./ (distVec + 0.25);
        if nPairs > 350
            [~, sortIdx] = sort(pairStrength, 'descend');
            keep = sortIdx(1:350);
            pList = pList(keep);
            qList = qList(keep);
            nP_dot_d = nP_dot_d(keep);
            nQ_dot_d = nQ_dot_d(keep);
            distVec = distVec(keep);
            dVec = dVec(keep, :);
            nPairs = 350;
        end
    end

    cp = struct();
    cp.numPairs = nPairs;
    if nPairs > 0
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
        % Effective mutual capture area of the concave pair
        cp.effArea = sqrt(cp.aP .* cp.aQ) .* sqrt(nP_dot_d .* nQ_dot_d);
        cp.charSize = repmat(sqrt(mean([cp.aP; cp.aQ])), 1, 1);
    else
        cp.rP = zeros(0, 3); cp.rQ = zeros(0, 3);
        cp.nP = zeros(0, 3); cp.nQ = zeros(0, 3);
        cp.dPQ = zeros(0, 3); cp.distPQ = zeros(0, 1);
        cp.nP_dot_dPQ = zeros(0, 1); cp.nQ_dot_dQP = zeros(0, 1);
        cp.effArea = zeros(0, 1); cp.charSize = 0.2;
    end
end
