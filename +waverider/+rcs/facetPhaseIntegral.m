function I = facetPhaseIntegral(v1, v2, v3, nHat, w)
%FACETPHASEINTEGRAL Exact closed-form integral of exp(i w.r) over flat triangles.
%   I = FACETPHASEINTEGRAL(v1, v2, v3, nHat, w) returns the complex surface
%   integral  I(f, d) = \int_{T_f} exp(i w_d . r) dA  over each triangle T_f
%   (vertices v1, v2, v3 [nF x 3], unit normal nHat [nF x 3]) for each phase
%   gradient vector w_d [nD x 3]. Output I is [nF x nD].
%
%   Uses Gordon's (1975) edge-sum form obtained from the in-plane divergence
%   theorem with G = w_t exp(i w.r) / (i |w_t|^2):
%
%     I = sum_edges [ w_t . (dr_e x n) ] / (i |w_t|^2) * exp(i w . r_mid,e)
%                  * sinc(w . dr_e / 2),            sinc(x) = sin(x)/x
%
%   where w_t is the component of w tangent to the facet. When |w_t| * L_max
%   is tiny the edge sum cancels catastrophically, so the second-order-
%   accurate centroid limit  I = A exp(i w . r_c)  is used instead.
%
%   For monostatic physical optics w = 2 k u_radar.

    nF = size(v1, 1);
    nD = size(w, 1);

    E = {v2 - v1, v3 - v2, v1 - v3};
    M = {0.5 .* (v1 + v2), 0.5 .* (v2 + v3), 0.5 .* (v3 + v1)};

    crossV = cross(v2 - v1, v3 - v1, 2);
    area = 0.5 .* sqrt(sum(crossV.^2, 2));
    rc = (v1 + v2 + v3) ./ 3.0;
    Lmax = max([vecnorm(E{1}, 2, 2), vecnorm(E{2}, 2, 2), vecnorm(E{3}, 2, 2)], [], 2);

    nDotW = nHat * w.';                                  % [nF x nD]
    wt2 = repmat(sum(w.^2, 2).', nF, 1) - nDotW.^2;      % |w_t|^2
    wt2 = max(wt2, 0.0);

    I = zeros(nF, nD);
    for e = 1:3
        m = cross(E{e}, nHat, 2);                        % outward in-plane edge normal * length
        wm = m * w.';                                    % w_t . m  (m lies in the plane)
        wd = E{e} * w.';                                 % w . dr_e
        wr = M{e} * w.';                                 % w . r_mid
        x = 0.5 .* wd;
        s = ones(size(x));
        nz = abs(x) > 1e-12;
        s(nz) = sin(x(nz)) ./ x(nz);
        I = I + wm .* exp(1i .* wr) .* s;
    end
    small = sqrt(wt2) .* Lmax < 1e-3;
    I(~small) = I(~small) ./ (1i .* wt2(~small));
    centroidLimit = area .* exp(1i .* (rc * w.'));       % [nF x nD]
    I(small) = centroidLimit(small);
end
