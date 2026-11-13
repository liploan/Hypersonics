function [S_mean, S_pol, segContrib] = evaluateRimPTD(rimData, uRadar, k0, eH, wantSeg)
%EVALUATERIMPTD Polarimetric Ufimtsev/Mitzner ILDC monostatic PTD base-rim diffraction.
%   [S_mean, S_pol, segContrib] = EVALUATERIMPTD(rimData, uRadar, k0, eH, wantSeg)
%   evaluates the singularity-free Ufimtsev/Mitzner Physical Theory of
%   Diffraction (PTD) fringe current contribution across all sharp aft base
%   perimeter wedge segments for look vectors uRadar [nDirs x 3] at
%   wavenumber k0 [rad/m] and horizontal BSA polarization vectors eH [nDirs x 3].
%
%   Outputs:
%     S_mean     : [nDirs x 1] polarization-mean PTD fringe complex amplitude
%     S_pol      : struct with fields HH, VV, HV, VH ([nDirs x 1], HV == VH)
%     segContrib : optional [nDirs x 1] cell of [re_dot_u, aMean, aHH, aVV, aHV]
%                  for wideband HRRP / RTI range-bin deposition.

    if nargin < 5
        wantSeg = false;
    end
    nD = size(uRadar, 1);
    S_mean = zeros(nD, 1);
    S_pol = struct('HH', zeros(nD, 1), 'VV', zeros(nD, 1), ...
                   'HV', zeros(nD, 1), 'VH', zeros(nD, 1));
    segContrib = cell(nD, 1);
    if isempty(rimData) || rimData.count == 0
        return;
    end

    eV = cross(eH, uRadar, 2);
    hasConcaveAft = rimData.hasConcaveAft;
    halfPi = 0.5 * pi;
    sqrtPi = sqrt(pi);

    re   = rimData.re;
    len  = rimData.len;
    tHat = rimData.tHat;
    tau1 = rimData.tau1;
    n1   = rimData.n1;
    n2   = rimData.n2;
    nW   = rimData.nWedge;
    nPi  = rimData.nPi;
    Xw   = rimData.Xwedge;
    sPN  = rimData.sinPiOverN;
    cPN  = rimData.cosPiOverN;
    pRes = rimData.poleRes;

    for d = 1:nD
        u = uRadar(d, :);
        eh = eH(d, :);
        ev = eV(d, :);

        lit1 = (n1 * u.') > 1e-9;
        lit2 = ((n2 * u.') > 1e-9) & (~hasConcaveAft || (u(1) > 1e-7));
        active = lit1 | lit2;
        if ~any(active)
            continue;
        end

        tu = tHat * u.';
        sinB2 = 1.0 - tu.^2;
        active = active & (sinB2 >= 1e-6);
        idx = find(active);
        if isempty(idx)
            continue;
        end

        tuA   = tu(idx);
        sinBA = sqrt(sinB2(idx));
        invSB = 1.0 ./ sinBA;
        tA    = tHat(idx, :);

        uPerp = (u - tuA .* tA) .* invSB;
        xw = sum(uPerp .* tau1(idx, :), 2);
        yw = sum(uPerp .* n1(idx, :), 2);
        phi = atan2(yw, xw);
        neg = phi < 0.0;
        phi(neg) = phi(neg) + 2.0 * pi;

        validPhi = (phi > 1e-5) & (phi < (nPi(idx) - 1e-5));
        if ~any(validPhi)
            continue;
        end
        idx   = idx(validPhi);
        phi   = phi(validPhi);
        tuA   = tuA(validPhi);
        sinBA = sinBA(validPhi);
        invSB = invSB(validPhi);
        tA    = tHat(idx, :);
        lit2A = lit2(idx);

        phi2 = nPi(idx) - phi;
        Yf = zeros(size(phi));

        near1 = abs(phi - halfPi) < 1e-5;
        near2 = (~near1) & lit2A & (abs(phi2 - halfPi) < 1e-5);
        reg   = ~(near1 | near2);

        if any(near1)
            po2_n1 = zeros(sum(near1), 1);
            m2 = lit2A(near1) & (phi2(near1) > 0.0) & (phi2(near1) < pi);
            p2Sub = phi2(near1);
            po2_n1(m2) = -0.5 .* tan(p2Sub(m2));
            Yf(near1) = pRes(idx(near1)) - po2_n1;
        end
        if any(near2)
            po1_n2 = zeros(sum(near2), 1);
            p1Sub = phi(near2);
            m1 = (p1Sub > 0.0) & (p1Sub < pi);
            po1_n2(m1) = -0.5 .* tan(p1Sub(m1));
            Yf(near2) = pRes(idx(near2)) - po1_n2;
        end
        if any(reg)
            ir = idx(reg);
            pr = phi(reg);
            p2r = phi2(reg);
            symSum = sPN(ir) ./ (cPN(ir) - cos(2.0 .* pr ./ nW(ir)));
            po1 = zeros(size(pr));
            m1 = (pr > 0.0) & (pr < pi);
            po1(m1) = -0.5 .* tan(pr(m1));
            po2 = zeros(size(p2r));
            m2 = lit2A(reg) & (p2r > 0.0) & (p2r < pi);
            po2(m2) = -0.5 .* tan(p2r(m2));
            Yf(reg) = symSum - po1 - po2;
        end

        fs = Yf - Xw(idx);
        fh = Yf + Xw(idx);

        hPar = (tA * eh.') .* invSB;
        vPar = (tA * ev.') .* invSB;
        hPar2 = hPar.^2;
        vPar2 = vPar.^2;

        wHH = fs .* hPar2 + fh .* vPar2;
        wVV = fs .* vPar2 + fh .* hPar2;
        wHV = (fs - fh) .* hPar .* vPar;

        lenA = len(idx);
        xSinc = k0 .* lenA .* tuA;
        sincVal = ones(size(xSinc));
        nz = abs(xSinc) > 1e-8;
        sincVal(nz) = sin(xSinc(nz)) ./ xSinc(nz);
        sincVal(~nz) = 1.0 - (xSinc(~nz).^2) ./ 6.0;

        uDotRe = re(idx, :) * u.';
        elem = ((lenA .* sinBA) ./ sqrtPi) .* sincVal .* exp(1i * 2.0 * k0 .* uDotRe);

        aMean = elem .* Yf;
        aHH   = elem .* wHH;
        aVV   = elem .* wVV;
        aHV   = elem .* wHV;

        S_mean(d)   = sum(aMean);
        S_pol.HH(d) = sum(aHH);
        S_pol.VV(d) = sum(aVV);
        S_pol.HV(d) = sum(aHV);
        S_pol.VH(d) = S_pol.HV(d);

        if wantSeg
            segContrib{d} = [uDotRe, aMean, aHH, aVV, aHV];
        end
    end
end
