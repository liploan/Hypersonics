function [S_mb, contrib] = doubleBounceGOPO(cp, uRadar, k0, wantContrib)
%DOUBLEBOUNCEGOPO Monostatic two-bounce scattering of the concave underside (GO-PO).
%   [S_mb, contrib] = DOUBLEBOUNCEGOPO(cp, uRadar, k0, wantContrib) evaluates
%   the complex two-bounce scattering amplitude S_mb [nDirs x 1] (sigma =
%   |S|^2, same normalisation as single-bounce PO) for every ordered pair of
%   mutually visible lower-surface super-cells (P -> Q) in cp (from
%   buildConcavePairList), using the standard GO-PO (SBR-style) model:
%
%     1. Plane wave k_i = -u reflects specularly off planar cell P (GO):
%            k_r = k_i - 2 (k_i . n_P) n_P
%     2. The reflected beam is the projection of P along k_r onto the plane
%        of Q; it illuminates the overlap polygon O = Q  n  proj_kr(P).
%     3. PO on O radiates back toward the radar with the bistatic weighting
%            (i k / (2 sqrt(pi))) * n_Q . (u - k_r) * Gamma
%        and the exact phase integral over O (Gordon, facetPhaseIntegral),
%        with total path  L(r) = k_r . r + c - u . r,  c = (k_i - k_r) . r_P.
%
%   Gamma = -1 (scalar soft-surface reflection; polarisation-averaged).
%   Validated against the closed-form 90-deg dihedral RCS
%   sigma = 16 pi a^2 b^2 sin^2(phi) / lambda^2 (tests/test_physics_validation).
%
%   Limitations: no blockage/shadowing by other parts of the body, flat
%   super-cells approximate the curved surface, and GO reflection is only
%   meaningful when cells are electrically large (roughly S-band and up).
%
%   If wantContrib is true, contrib{d} returns [amp, pathL/2, rCentroid(1:3)]
%   rows for every non-zero contribution in direction d, used to place the
%   delayed cavity returns in range profiles.

    if nargin < 4
        wantContrib = false;
    end
    nD = size(uRadar, 1);
    S_mb = zeros(nD, 1);
    contrib = cell(nD, 1);
    if cp.numPairs == 0
        return;
    end

    % Ordered pairs (P -> Q) and (Q -> P)
    P = [cp.pIdx; cp.qIdx];
    Q = [cp.qIdx; cp.pIdx];
    rP = cp.cellCenter(P, :);  nP = cp.cellNormal(P, :);
    rQ = cp.cellCenter(Q, :);  nQ = cp.cellNormal(Q, :);
    radP = cp.cellRadius(P);   radQ = cp.cellRadius(Q);
    Gamma = -1.0;
    pref = 1i * k0 / (2.0 * sqrt(pi));

    for d = 1:nD
        u = uRadar(d, :);
        ki = -u;
        cosP = nP * u.';                                   % P lit by radar
        kr = ki - 2.0 .* (nP * ki.') .* nP;                % [nPairs x 3]
        nQkr = sum(nQ .* kr, 2);                           % < 0 : Q lit by reflected beam
        tHit = sum(nQ .* (rQ - rP), 2) ./ min(nQkr, -1e-9);
        hit = rP + tHit .* kr;
        cand = find(cosP > 1e-6 & nQkr < -1e-6 & tHit > 0 & ...
                    vecnorm(hit - rQ, 2, 2) < (radQ + radP ./ max(-nQkr, 0.05)));
        if isempty(cand)
            continue;
        end
        Sd = 0.0;
        rows = zeros(numel(cand), 5);
        nRows = 0;
        for m = cand.'
            pv = squeeze(cp.cellVerts(P(m), :, :));        % 4 x 3
            qv = squeeze(cp.cellVerts(Q(m), :, :));
            nq = nQ(m, :); kvec = kr(m, :);
            % Project P's polygon along k_r onto the plane of Q
            tv = ((rQ(m, :) - pv) * nq.') ./ nQkr(m);
            if any(tv < -1e-9)
                continue;                                  % part of P lies beyond Q: not a valid bounce
            end
            proj = pv + tv .* kvec;
            poly = clipConvexPolygon(proj, qv, nq);
            if size(poly, 1) < 3
                continue;
            end
            % Ensure CCW orientation about n_Q (projection can mirror the winding)
            areaVec = sum(cross(poly - poly(1, :), circshift(poly, -1, 1) - poly(1, :), 2), 1);
            if dot(areaVec, nq) < 0
                poly = flipud(poly);
            end
            if abs(dot(areaVec, nq)) < 1e-12
                continue;
            end
            w = k0 .* (u - kvec);
            c = (ki - kvec) * rP(m, :).';
            % Fan-triangulate the convex overlap polygon and integrate exactly
            nT = size(poly, 1) - 2;
            v1 = repmat(poly(1, :), nT, 1);
            v2 = poly(2:end - 1, :);
            v3 = poly(3:end, :);
            Ipoly = sum(waverider.rcs.facetPhaseIntegral(v1, v2, v3, repmat(nq, nT, 1), w));
            amp = pref * (nq * (u - kvec).') * Gamma * Ipoly * exp(-1i * k0 * c);
            Sd = Sd + amp;
            if wantContrib
                triA = 0.5 .* vecnorm(cross(v2 - v1, v3 - v1, 2), 2, 2);
                rO = sum(((v1 + v2 + v3) ./ 3) .* triA, 1) ./ max(sum(triA), 1e-15);
                Lhalf = 0.5 * ((kvec - u) * rO.' + c);
                nRows = nRows + 1;
                rows(nRows, :) = [amp, Lhalf, rO];
            end
        end
        S_mb(d) = Sd;
        if wantContrib
            contrib{d} = rows(1:nRows, :);
        end
    end
end

function out = clipConvexPolygon(subject, clip, nHat)
%CLIPCONVEXPOLYGON Sutherland-Hodgman clipping of a planar polygon by a convex
%   polygon (both in the same plane, vertices CCW about nHat).
    out = subject;
    nC = size(clip, 1);
    for e = 1:nC
        if size(out, 1) == 0
            return;
        end
        a = clip(e, :);
        b = clip(mod(e, nC) + 1, :);
        inward = cross(nHat, b - a);               % points into the CCW polygon
        dist = (out - a) * inward.';
        inp = out;
        out = zeros(0, 3);
        nI = size(inp, 1);
        for i = 1:nI
            j = mod(i, nI) + 1;
            di = dist(i); dj = dist(j);
            if di >= 0
                out(end + 1, :) = inp(i, :); %#ok<AGROW>
            end
            if (di >= 0) ~= (dj >= 0)
                t = di / (di - dj);
                out(end + 1, :) = inp(i, :) + t .* (inp(j, :) - inp(i, :)); %#ok<AGROW>
            end
        end
    end
end
