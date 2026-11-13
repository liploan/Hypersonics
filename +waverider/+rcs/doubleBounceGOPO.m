function [S_mb, contrib, S_mb_pol] = doubleBounceGOPO(cp, uRadar, k0, wantContrib, eHOpt)
%DOUBLEBOUNCEGOPO Monostatic two-bounce scattering of the concave underside (GO-PO).
%   [S_mb, contrib, S_mb_pol] = DOUBLEBOUNCEGOPO(cp, uRadar, k0, wantContrib, eHOpt)
%   evaluates both the scalar baseline two-bounce scattering amplitude S_mb
%   [nDirs x 1] (Gamma = -1) and the full 3D vector dyadic Sinclair scattering
%   matrix S_mb_pol (fields HH, VV, HV, VH, each [nDirs x 1]; sigma = |S|^2,
%   IEEE Backscatter Alignment convention) for every ordered pair of mutually
%   visible lower-surface super-cells (P -> Q) in cp (from buildConcavePairList):
%
%     1. Plane wave k_i = -u reflects specularly off planar PEC cell P (GO):
%            k_r   = k_i - 2 (k_i . n_P) n_P
%            e_{r,q} = 2 (n_P . e_q) n_P - e_q,   q in {H, V}
%     2. The reflected beam projects P along k_r onto the plane of Q and
%        illuminates the exact clipped polygon O = Q \cap proj_{k_r}(P).
%     3. Secondary PO current J_Q = 2 n_Q x (k_r x e_{r,q}) / eta_0 radiates
%        back toward the radar with real geometric dyadic weight:
%            W_{pq} = 2 [ (n_Q . e_{r,q})(e_p . k_r) - (n_Q . k_r)(e_p . e_{r,q}) ]
%        and the exact Gordon contour integral over O (facetPhaseIntegral),
%        sharing ray culling, polygon clipping, and phase integration across
%        all polarizations with near-zero (< 3%) overhead.
%
%   Validated against the closed-form 90-deg dihedral (0-deg co-pol S_HH = -S_VV,
%   45-deg canted pure HV/VH depolarization, and roll-invariant circular OP) in
%   tests/test_physics_validation.m.
%
%   If wantContrib is true, contrib{d} returns [ampMean, pathL/2, rCentroid(1:3),
%   ampHH, ampVV, ampHV, ampVH] rows [nHits x 9] for down-range HRRP/RTI profiles.

    if nargin < 4 || isempty(wantContrib)
        wantContrib = false;
    end
    if nargin < 5
        eHOpt = [];
    end
    nD = size(uRadar, 1);
    S_mb = zeros(nD, 1);
    S_mb_pol = struct( ...
        'HH', zeros(nD, 1), ...
        'VV', zeros(nD, 1), ...
        'HV', zeros(nD, 1), ...
        'VH', zeros(nD, 1));
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
        if ~isempty(eHOpt)
            eH = eHOpt(d, :);
        else
            rhoXY = hypot(u(1), u(2));
            if rhoXY > 1e-15
                eH = [u(2), -u(1), 0.0] ./ rhoXY;
            else
                eH = [0.0, 1.0, 0.0];
            end
        end
        eV = cross(eH, u);

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
        Sd_HH = 0.0;
        Sd_VV = 0.0;
        Sd_HV = 0.0;
        rows = zeros(numel(cand), 9);
        nRows = 0;
        for m = cand.'
            pv = squeeze(cp.cellVerts(P(m), :, :));        % 4 x 3
            qv = squeeze(cp.cellVerts(Q(m), :, :));
            np = nP(m, :);
            nq = nQ(m, :);
            kvec = kr(m, :);
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
            baseAmp = pref * Ipoly * exp(-1i * k0 * c);

            % Scalar PEC baseline (Gamma = -1)
            amp = baseAmp * (nq * (u - kvec).') * Gamma;

            % Full 3D vector PEC reflection off P and PO radiation off Q (IEEE BSA)
            erH = 2.0 * (np * eH.') .* np - eH;
            erV = 2.0 * (np * eV.') .* np - eV;
            eH_kr = eH * kvec.';
            eV_kr = eV * kvec.';
            nq_erH = nq * erH.';
            nq_erV = nq * erV.';
            nq_kr  = nQkr(m);

            wHH = 2.0 * (nq_erH * eH_kr - nq_kr * (eH * erH.'));
            wVV = 2.0 * (nq_erV * eV_kr - nq_kr * (eV * erV.'));
            wHV = 2.0 * (nq_erV * eH_kr - nq_kr * (eH * erV.'));
            wVH = 2.0 * (nq_erH * eV_kr - nq_kr * (eV * erH.'));
            wCross = 0.5 * (wHV + wVH);                    % Lorentz-reciprocal GO-PO symmetrization

            ampHH = baseAmp * wHH;
            ampVV = baseAmp * wVV;
            ampHV = baseAmp * wCross;

            Sd    = Sd + amp;
            Sd_HH = Sd_HH + ampHH;
            Sd_VV = Sd_VV + ampVV;
            Sd_HV = Sd_HV + ampHV;
            if wantContrib
                triA = 0.5 .* vecnorm(cross(v2 - v1, v3 - v1, 2), 2, 2);
                rO = sum(((v1 + v2 + v3) ./ 3) .* triA, 1) ./ max(sum(triA), 1e-15);
                Lhalf = 0.5 * ((kvec - u) * rO.' + c);
                nRows = nRows + 1;
                rows(nRows, :) = [amp, Lhalf, rO, ampHH, ampVV, ampHV, ampHV];
            end
        end
        S_mb(d)        = Sd;
        S_mb_pol.HH(d) = Sd_HH;
        S_mb_pol.VV(d) = Sd_VV;
        S_mb_pol.HV(d) = Sd_HV;
        S_mb_pol.VH(d) = Sd_HV;
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
