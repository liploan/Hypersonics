function rimData = buildRimPTDEdges(V, F, patchID, faceNormals, hasConcaveAft)
%BUILDRIMPTDEDGES Precompute exterior sharp aft base-rim wedge geometry for Ufimtsev PTD.
%   rimData = BUILDRIMPTDEDGES(V, F, patchID, faceNormals, hasConcaveAft)
%   extracts all sharp convex exterior wedge segments along the topological
%   seam between outer vehicle skin facets (patchID ~= 3) and aft base
%   facets (patchID == 3), and precomputes each segment's right-handed wedge
%   frame (tau1, n1, tHat), exterior wedge index n in (1, 2], and Ufimtsev
%   diffraction constants.

    if nargin < 5
        hasConcaveAft = false;
    end
    rimData = struct('count', 0, 'hasConcaveAft', logical(hasConcaveAft));
    if isempty(F) || isempty(patchID)
        return;
    end

    nF = size(F, 1);
    % Collect all half-edges [u, v, f, ic, isBase]
    heU = [F(:, 1); F(:, 2); F(:, 3)];
    heV = [F(:, 2); F(:, 3); F(:, 1)];
    heC = [F(:, 3); F(:, 1); F(:, 2)];
    heF = repmat((1:nF).', 3, 1);
    heB = (patchID(heF) == 3);

    eMin = min(heU, heV);
    eMax = max(heU, heV);
    [~, ~, grp] = unique([eMin, eMax], 'rows');
    nUnique = max(grp);

    skinFace = zeros(nUnique, 1);
    skinOpp  = zeros(nUnique, 1);
    skinA    = zeros(nUnique, 1);
    skinB    = zeros(nUnique, 1);
    baseFace = zeros(nUnique, 1);
    baseOpp  = zeros(nUnique, 1);

    for k = 1:numel(grp)
        g = grp(k);
        if heB(k)
            baseFace(g) = heF(k);
            baseOpp(g)  = heC(k);
        else
            skinFace(g) = heF(k);
            skinOpp(g)  = heC(k);
            skinA(g)    = heU(k);
            skinB(g)    = heV(k);
        end
    end

    cand = find((skinFace > 0) & (baseFace > 0));
    if isempty(cand)
        return;
    end

    mMax = numel(cand);
    re   = zeros(mMax, 3);
    len  = zeros(mMax, 1);
    tHat = zeros(mMax, 3);
    tau1 = zeros(mMax, 3);
    n1   = zeros(mMax, 3);
    n2   = zeros(mMax, 3);
    nW   = zeros(mMax, 1);
    cnt  = 0;

    for idx = 1:mMax
        g = cand(idx);
        va  = V(skinA(g), :);
        vb  = V(skinB(g), :);
        vc1 = V(skinOpp(g), :);
        vc2 = V(baseOpp(g), :);
        rmid = 0.5 .* (va + vb);
        de = vb - va;
        le = norm(de);
        if le < 1e-9
            continue;
        end

        nf1 = faceNormals(skinFace(g), :);
        nf2 = faceNormals(baseFace(g), :);
        % Strict convex exterior wedge test
        if dot(vc2 - rmid, nf1) >= -1e-9
            continue;
        end

        cos12 = max(-1.0, min(1.0, dot(nf1, nf2)));
        nwVal = 1.0 + acos(cos12) / pi;
        if nwVal <= 1.02
            continue;
        end

        dHat = de ./ le;
        q1 = cross(nf1, dHat);
        if dot(vc1 - rmid, q1) >= 0.0
            sgn = 1.0;
        else
            sgn = -1.0;
        end

        cnt = cnt + 1;
        re(cnt, :)   = rmid;
        len(cnt)     = le;
        tau1(cnt, :) = sgn .* q1;
        tHat(cnt, :) = sgn .* dHat; % tau1 x n1 = tHat
        n1(cnt, :)   = nf1;
        n2(cnt, :)   = nf2;
        nW(cnt)      = nwVal;
    end

    if cnt == 0
        return;
    end

    rimData.count         = cnt;
    rimData.re            = re(1:cnt, :);
    rimData.len           = len(1:cnt);
    rimData.tHat          = tHat(1:cnt, :);
    rimData.tau1          = tau1(1:cnt, :);
    rimData.n1            = n1(1:cnt, :);
    rimData.n2            = n2(1:cnt, :);
    rimData.nWedge        = nW(1:cnt);
    rimData.nPi           = rimData.nWedge .* pi;
    rimData.Xwedge        = -1.0 ./ (rimData.nWedge .* tan(pi ./ (2.0 .* rimData.nWedge)));
    rimData.sinPiOverN    = sin(pi ./ rimData.nWedge) ./ rimData.nWedge;
    rimData.cosPiOverN    = cos(pi ./ rimData.nWedge);
    rimData.poleRes       = -0.5 ./ (rimData.nWedge .* tan(pi ./ rimData.nWedge));
end
