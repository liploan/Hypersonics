classdef WaveriderBody
%WAVERIDERBODY Parametric watertight hypersonic waveriding glide body.
%   body = waverider.WaveriderBody(cfg) generates the 3D lower compression
%   surface, upper freestream/canopy surface, optional conformal transition
%   frustum and conical rocket booster stage, and aft base cap, stitches
%   them into a verified 2-manifold watertight triangulation, and computes
%   volumetric, sectional, and preliminary inviscid hypersonic aero metrics.

    properties
        Config waverider.WaveriderConfig
        SurfaceData struct
        MeshData struct
        Report struct
        AeroMetrics struct
        SectionalAttributes struct
    end

    methods
        function obj = WaveriderBody(varargin)
            %WAVERIDERBODY Construct and generate a watertight waverider body.
            if nargin == 0
                cfg = waverider.WaveriderConfig();
            elseif nargin == 1 && isa(varargin{1}, 'waverider.WaveriderConfig')
                cfg = varargin{1};
            else
                cfg = waverider.WaveriderConfig(varargin{:});
            end
            obj.Config = cfg;
            obj = obj.generate();
        end

        function obj = generate(obj)
            %GENERATE Run surface generation, frustum lofting, stitching, and verification.
            cfg = obj.Config;
            cfg.validate();
            switch lower(cfg.Method)
                case {'conical', 'osculating_cone'}
                    obj.SurfaceData = waverider.geometry.generateOsculatingCone(cfg);
                case 'caret'
                    obj.SurfaceData = waverider.geometry.generateCaretWedge(cfg);
                case 'cst'
                    obj.SurfaceData = waverider.geometry.generateCSTBody(cfg);
                case 'rv'
                    obj.SurfaceData = waverider.geometry.generateReentryVehicle(cfg);
                otherwise
                    error('WaveriderBody:UnsupportedMethod', ...
                        'Unsupported method "%s".', cfg.Method);
            end

            % Append optional conformal transition frustum & conical rocket stage
            obj.SurfaceData = waverider.geometry.appendConformalFrustum(obj.SurfaceData, cfg);

            obj.MeshData = waverider.geometry.stitchWatertightMesh(obj.SurfaceData, cfg);
            obj.Report = waverider.mesh.verifyWatertight(obj.MeshData, cfg, obj.SurfaceData);
            obj.AeroMetrics = obj.computePreliminaryAero();
            obj.SectionalAttributes = obj.computeSectionalAttributes();
        end

        function fig = plot(obj, varargin)
            %PLOT Render a 4-panel 3D engineering plot of the waverider.
            fig = waverider.viz.plotWaverider(obj, varargin{:});
        end

        function fig = plotGeometryLines(obj, varargin)
            %PLOTGEOMETRYLINES Render 3-view lofting lines, rib slices, and area rule distribution.
            fig = waverider.viz.plotGeometryLines(obj, varargin{:});
        end

        function fig = plotKeyAttributes(obj, varargin)
            %PLOTKEYATTRIBUTES Render aerodynamic polars, shock flowfield, and key attributes.
            fig = waverider.viz.plotKeyAttributes(obj, varargin{:});
        end

        function fig = plotBoosterIntegration(obj, varargin)
            %PLOTBOOSTERINTEGRATION Render conformal aft frustum & conical rocket adapter analysis.
            fig = waverider.viz.plotBoosterIntegration(obj, varargin{:});
        end

        function rcsData = computeRCS(obj, freqGHz, azDeg, elDeg, varargin)
            %COMPUTERCS Compute monostatic RCS (UHF to X-band) with concave multi-bounce & polarimetry.
            rcsData = waverider.rcs.computeRCS(obj, freqGHz, azDeg, elDeg, varargin{:});
        end

        function rtiData = computeRangeProfile(obj, freqGHz, bandwidthMHz, azVecDeg, elDeg, varargin)
            %COMPUTERANGEPROFILE Compute wideband range-angle-intensity (RTI / HRRP) sinogram.
            rtiData = waverider.rcs.computeRangeProfile(obj, freqGHz, bandwidthMHz, azVecDeg, elDeg, varargin{:});
        end

        function [fig, summaryMetrics] = plotRCS(obj, varargin)
            %PLOTRCS Render 4-panel UHF-to-X-band RCS & concave underside scattering dashboard.
            [fig, summaryMetrics] = waverider.viz.plotRCSAnalysis(obj, varargin{:});
        end

        function [fig, rtiData] = plotRangeAngleIntensity(obj, varargin)
            %PLOTRANGEANGLEINTENSITY Render 4-panel wideband RTI sinogram & HRRP cuts.
            [fig, rtiData] = waverider.viz.plotRangeAngleIntensity(obj, varargin{:});
        end

        function filePath = exportSTL(obj, filePath, mode)
            %EXPORTSTL Export the watertight mesh to binary or ASCII STL.
            if nargin < 3
                mode = 'binary';
            end
            filePath = waverider.mesh.exportSTL(obj.MeshData, filePath, mode);
        end

        function filePath = exportOBJ(obj, filePath)
            %EXPORTOBJ Export the watertight mesh to Wavefront OBJ format.
            filePath = waverider.mesh.exportOBJ(obj.MeshData, filePath);
        end

        function polar = computeAeroPolar(obj, alphaDegVec)
            %COMPUTEAEROPOLAR Compute CL, CD_wave, L/D, viscous L/D, and Cm across AoA alpha [deg].
            %   Uses the exact same calibrated inviscid + Prandtl-Meyer surface pressure
            %   model as AeroMetrics so polar(alpha=0) matches AeroMetrics to machine precision.
            if nargin < 2 || isempty(alphaDegVec)
                alphaDegVec = linspace(-4.0, 16.0, 41);
            end
            alphaDegVec = alphaDegVec(:);
            nA = numel(alphaDegVec);

            CL = zeros(nA, 1);
            CD = zeros(nA, 1);
            CY = zeros(nA, 1);
            LD = zeros(nA, 1);
            CD_visc = zeros(nA, 1);
            LD_visc = zeros(nA, 1);
            Cm = zeros(nA, 1);

            geomCache = obj.buildAeroGeometryCache();

            for k = 1:nA
                res = obj.evaluateAeroState(alphaDegVec(k), geomCache);
                CL(k) = res.CL;
                CD(k) = res.CD_wave;
                CY(k) = res.CY;
                LD(k) = res.LD_inviscid;
                CD_visc(k) = res.CD_total_viscous_est;
                LD_visc(k) = res.LD_viscous_est;
                Cm(k) = res.Cm;
            end

            polar = struct();
            polar.AlphaDeg = alphaDegVec;
            polar.CL = CL;
            polar.CD_wave = CD;
            polar.CY = CY;
            polar.LD_inviscid = LD;
            polar.CD_friction_est = geomCache.CD_friction_est;
            polar.CD_viscous_est = CD_visc;
            polar.LD_viscous_est = LD_visc;
            polar.Cm = Cm;
        end

        function sec = computeSectionalAttributes(obj)
            %COMPUTESECTIONALATTRIBUTES Compute axial cross-sectional & spanwise distributions.
            cfg = obj.Config;
            surfData = obj.SurfaceData;
            L = cfg.Length;
            NxWave = surfData.NxWaverider;

            % Total vehicle axial extent (including frustum & rocket stage if enabled)
            xMax = max(surfData.X_up(:));
            nStations = 32;
            xStations = linspace(0.05 * L, xMax, nStations).';
            areaX = zeros(nStations, 1);
            halfSpanX = zeros(nStations, 1);
            thicknessX = zeros(nStations, 1);

            yLE = surfData.Y_up(1, :);
            xLE = surfData.X_up(1, :);
            Ny = numel(yLE);
            isAxisym = isfield(surfData, 'IsAxisymmetric') && surfData.IsAxisymmetric;

            for s = 1:nStations
                xs = xStations(s);
                if isAxisym
                    xAxisVec = surfData.X_up(:, ceil(Ny / 2));
                    yUpSlice  = interp1(xAxisVec, surfData.Y_up, xs, 'linear', 'extrap');
                    zUpSlice  = interp1(xAxisVec, surfData.Z_up, xs, 'linear', 'extrap');
                    yLowSlice = interp1(xAxisVec, surfData.Y_low, xs, 'linear', 'extrap');
                    zLowSlice = interp1(xAxisVec, surfData.Z_low, xs, 'linear', 'extrap');

                    yPoly = [yUpSlice, fliplr(yLowSlice)];
                    zPoly = [zUpSlice, fliplr(zLowSlice)];
                    areaX(s) = polyarea(yPoly, zPoly);
                    halfSpanX(s) = max(abs(yPoly));
                    thicknessX(s) = max(zUpSlice - zLowSlice);
                elseif xs <= L + 1e-7
                    ySlice = [];
                    zUpSlice = [];
                    zLowSlice = [];
                    for k = 1:Ny
                        if xLE(k) <= xs
                            xCol = surfData.X_up(1:NxWave, k);
                            if xCol(end) - xCol(1) < 1e-8
                                zu = surfData.Z_up(NxWave, k);
                                zl = surfData.Z_low(NxWave, k);
                            else
                                zu = interp1(xCol, surfData.Z_up(1:NxWave, k), xs, 'linear', 'extrap');
                                zl = interp1(xCol, surfData.Z_low(1:NxWave, k), xs, 'linear', 'extrap');
                            end
                            ySlice(end + 1) = yLE(k); %#ok<AGROW>
                            zUpSlice(end + 1) = zu; %#ok<AGROW>
                            zLowSlice(end + 1) = zl; %#ok<AGROW>
                        end
                    end
                    if numel(ySlice) >= 3
                        dzSlice = max(zUpSlice - zLowSlice, 0.0);
                        areaX(s) = trapz(ySlice, dzSlice);
                        halfSpanX(s) = max(abs(ySlice));
                        thicknessX(s) = max(dzSlice);
                    end
                else
                    % Station lies in the conformal frustum or conical rocket section (xs > L)
                    xAftVec = surfData.X_up(NxWave:end, round(Ny/2));
                    yUpSlice  = interp1(xAftVec, surfData.Y_up(NxWave:end, :), xs, 'linear', 'extrap');
                    zUpSlice  = interp1(xAftVec, surfData.Z_up(NxWave:end, :), xs, 'linear', 'extrap');
                    yLowSlice = interp1(xAftVec, surfData.Y_low(NxWave:end, :), xs, 'linear', 'extrap');
                    zLowSlice = interp1(xAftVec, surfData.Z_low(NxWave:end, :), xs, 'linear', 'extrap');

                    yPoly = [yUpSlice, fliplr(yLowSlice)];
                    zPoly = [zUpSlice, fliplr(zLowSlice)];
                    areaX(s) = polyarea(yPoly, zPoly);
                    halfSpanX(s) = max(abs(yPoly));
                    thicknessX(s) = max(zUpSlice - zLowSlice);
                end
            end

            cumVolumeX = cumtrapz(xStations, areaX);

            % Spanwise distributions on forebody (X <= L)
            if isAxisym
                spanY = surfData.Y_up(NxWave, :).';
                bHalfRef = max(surfData.RBase_m, 1e-12);
                rProf = surfData.RProfile;
                xProf = surfData.XProfile;
                xEntry = zeros(Ny, 1);
                for k = 1:Ny
                    rTgt = abs(spanY(k));
                    idx = find(rProf >= rTgt - 1e-12, 1, 'first');
                    if isempty(idx)
                        xEntry(k) = L;
                    elseif idx == 1
                        xEntry(k) = xProf(1);
                    else
                        dr = rProf(idx) - rProf(idx - 1);
                        if abs(dr) < 1e-12
                            xEntry(k) = xProf(idx - 1);
                        else
                            xEntry(k) = xProf(idx - 1) + ((rTgt - rProf(idx - 1)) / dr) * (xProf(idx) - xProf(idx - 1));
                        end
                    end
                end
                chordY = max(L - xEntry, 0.0);
            else
                spanY = yLE(:);
                bHalfRef = cfg.HalfSpan;
                chordY = max(L - xLE(:), 0.0);
            end
            baseHeightY = max(surfData.Z_up(NxWave, :).' - surfData.Z_low(NxWave, :).', 0.0);
            meanCpSpan = mean(surfData.Cp_low(1:NxWave, :), 1).';
            spanLoadNorm = meanCpSpan .* (chordY ./ L);

            sec = struct();
            sec.xStations = xStations;
            sec.xNorm = xStations ./ L;
            sec.AreaX = areaX;
            sec.CumulativeVolumeX = cumVolumeX;
            sec.HalfSpanX = halfSpanX;
            sec.MaxThicknessX = thicknessX;
            sec.SpanY = spanY;
            sec.SpanEta = spanY ./ bHalfRef;
            sec.ChordY = chordY;
            sec.BaseHeightY = baseHeightY;
            sec.MeanCpSpan = meanCpSpan;
            sec.SpanLoadNorm = spanLoadNorm;
        end

        function printSummary(obj)
            %PRINTSUMMARY Display a formatted summary of geometry, mesh, and aero metrics.
            cfg = obj.Config;
            rep = obj.Report;
            aero = obj.AeroMetrics;

            fprintf('\n----------------------------------------------------------------\n');
            if strcmpi(cfg.Method, 'rv')
                fprintf('  REENTRY VEHICLE (RV) SUMMARY (%s | %s NOSE)\n', ...
                    upper(strrep(cfg.RVBodyType, '_', '-')), upper(cfg.RVNoseType));
                fprintf('----------------------------------------------------------------\n');
                fprintf('  Design Mach / Conic Angles: M = %.2f | theta = [%s] deg\n', ...
                    cfg.Mach, num2str(obj.SurfaceData.ConeAnglesDeg, '%.2f '));
                fprintf('  Forebody Dimensions (LxD) : L = %.2f m, D_base = %.3f m (R_base = %.3f m)\n', ...
                    cfg.Length, 2.0 * obj.SurfaceData.RBase_m, obj.SurfaceData.RBase_m);
                fprintf('  Aft Plane Bulkhead        : %s (Cavity Depth = %.3f m)\n', ...
                    upper(obj.SurfaceData.RVAftType), obj.SurfaceData.AftCavityDepth_m);
            else
                fprintf('  WAVERIDER GLIDE BODY SUMMARY (%s)\n', upper(cfg.Method));
                fprintf('----------------------------------------------------------------\n');
                fprintf('  Design Mach / Shock Angle : M = %.2f, beta = %.2f deg\n', cfg.Mach, cfg.ShockAngleDeg);
                fprintf('  Forebody Dimensions (Lxb) : %.2f m x %.2f m (HalfSpan = %.2f m)\n', ...
                    cfg.Length, 2.0 * cfg.HalfSpan, cfg.HalfSpan);
                fprintf('  Upper Surface Mode        : %s (Canopy ratio = %.3f)\n', ...
                    cfg.UpperSurfaceType, cfg.CanopyHeightRatio);
            end
            if obj.SurfaceData.HasAftFrustum
                fprintf('  Aft Conformal Frustum     : L_frust = %.2f m -> Booster R = %.3f m (%s blend)\n', ...
                    cfg.FrustumLength, cfg.BoosterRadius, cfg.FrustumBlendLaw);
                if obj.SurfaceData.HasRocketStage
                    fprintf('  Conical Rocket Stage      : L_rocket = %.2f m, Half-Angle = %.2f deg (R_end = %.3f m)\n', ...
                        cfg.RocketStageLength, cfg.RocketConeAngleDeg, obj.SurfaceData.RocketEndRadius);
                end
            end
            fprintf('  Mesh Topology             : %d Vertices | %d Edges | %d Faces\n', ...
                rep.NumVertices, rep.NumEdges, rep.NumFaces);
            fprintf('  Watertight Verification   : %s (FreeBoundary = %d, NonManifold = %d, Euler = %d)\n', ...
                 ternaryStr(rep.IsWatertight, 'PASS [MANIFOLD WATERTIGHT]', 'FAIL'), ...
                rep.NumBoundaryEdges, rep.NumNonManifoldEdges, rep.EulerCharacteristic);
            fprintf('  Enclosed Volume           : %.4f m^3 (Centroid X = %.2f m, Z = %.3f m)\n', ...
                rep.EnclosedVolume, rep.VolumeCentroid(1), rep.VolumeCentroid(3));
            fprintf('  Planform / Wetted Area    : S_plan = %.3f m^2 | S_wet = %.3f m^2 | S_base = %.3f m^2\n', ...
                rep.PlanformArea, rep.WettedArea, rep.BaseArea);
            fprintf('  Volumetric Efficiency     : eta_V = V^(2/3)/S_plan = %.4f\n', rep.VolumetricEfficiency);
            fprintf('  Inviscid Wave Aero (AoA=0): C_L = %.4f | C_D,wave = %.4f | (L/D)_inviscid = %.2f\n', ...
                aero.CL, aero.CD_wave, aero.LD_inviscid);
            fprintf('  Viscous Est (30 km StdAtm): C_D,f = %.4f | C_D,base = %.4f | (L/D)_viscous = %.2f\n', ...
                aero.CD_friction_est, aero.CD_base_est, aero.LD_viscous_est);
            fprintf('----------------------------------------------------------------\n');
        end
    end

    methods (Access = private)
        function cache = buildAeroGeometryCache(obj)
            %BUILDAEROGEOMETRYCACHE Precompute facet geometry, calibrated lower-surface factors, & Eckert Cf.
            cfg = obj.Config;
            TR = obj.MeshData.TR;
            V = TR.Points;
            F = TR.ConnectivityList;

            v1 = V(F(:, 1), :);
            v2 = V(F(:, 2), :);
            v3 = V(F(:, 3), :);
            faceCenters = (v1 + v2 + v3) ./ 3.0;

            crossVec = 0.5 .* cross(v2 - v1, v3 - v1, 2); % Outward dS * n_hat
            faceAreas = sqrt(sum(crossVec.^2, 2)) + 1e-15;
            unitNormals = crossVec ./ faceAreas;

            M = cfg.Mach;
            g = cfg.Gamma;
            CpMax = (2.0 / (g * M^2)) * ...
                (((g + 1.0)^2 * M^2 / (4.0 * g * M^2 - 2.0 * (g - 1.0)))^(g / (g - 1.0)) * ...
                 ((1.0 - g + 2.0 * g * M^2) / (g + 1.0)) - 1.0);

            isUpper  = (obj.MeshData.FacePatchID == 1);
            isLower  = (obj.MeshData.FacePatchID == 2);
            isBase   = (obj.MeshData.FacePatchID == 3);
            isAftExt = (obj.MeshData.FacePatchID == 4) | (obj.MeshData.FacePatchID == 5);
            activeMask = isUpper | isLower | isAftExt;

            % Calibrate embedded Newtonian factor K_cal on facets carrying exact design Cp
            % (lower surface on waveriders; both dorsal & ventral surfaces on sharp TM cones).
            if isfield(obj.MeshData, 'FaceCpDesign') && numel(obj.MeshData.FaceCpDesign) == size(F, 1)
                faceCpDesignVec = obj.MeshData.FaceCpDesign(:);
            else
                NxWave = obj.SurfaceData.NxWaverider;
                faceCpDesignVec = nan(size(F, 1), 1);
                faceCpDesignVec(isLower) = mean(obj.SurfaceData.Cp_low(1:NxWave, :), 'all');
            end
            isCal = (isUpper | isLower) & isfinite(faceCpDesignVec);
            isUncal = activeMask & ~isCal;

            sinImpact0Cal = max(-unitNormals(isCal, 1), 0.0);
            cpDesignCal = faceCpDesignVec(isCal);
            cpMN0Cal = CpMax .* (sinImpact0Cal.^2);
            kCal = ones(sum(isCal), 1);
            validWindward0 = sinImpact0Cal > 1e-4;
            kCal(validWindward0) = cpDesignCal(validWindward0) ./ cpMN0Cal(validWindward0);
            kCal(~validWindward0) = 1.0;

            % Eckert (1955) reference-temperature turbulent skin-friction estimate at h = 30 km (US Std Atm)
            T_inf = 226.51;
            rho_inf = 0.01841;
            a_inf = sqrt(g * 287.05 * T_inf);
            U_inf = M * a_inf;
            T_wall = 1000.0;
            Pr = 0.72;
            rRec = Pr^(1.0 / 3.0);
            T_aw = T_inf * (1.0 + 0.5 * rRec * (g - 1.0) * M^2);
            T_star = T_inf + 0.5 * (T_wall - T_inf) + 0.22 * (T_aw - T_inf);
            mu_star = 1.458e-6 * (T_star^1.5) / (T_star + 110.4);
            rho_star = rho_inf * (T_inf / T_star);
            L_char = max(0.67 * cfg.Length, 0.1);
            Re_star = max(rho_star * U_inf * L_char / mu_star, 1e4);
            Cf_star = 0.074 / (Re_star^0.2);
            Cf_inf = Cf_star * (rho_star / rho_inf);

            Sref = max(obj.Report.PlanformArea, 1e-12);
            CD_friction_est = Cf_inf * (obj.Report.WettedArea / Sref);
            Cp_base = -1.0 / (g * M^2);
            projBaseArea = sum(max(crossVec(isBase, 1), 0.0));
            CD_base_est = -Cp_base * (projBaseArea / Sref);

            cache = struct();
            cache.crossVec = crossVec;
            cache.faceAreas = faceAreas;
            cache.unitNormals = unitNormals;
            cache.rArm = faceCenters - obj.Report.VolumeCentroid;
            cache.M = M;
            cache.g = g;
            cache.CpMax = CpMax;
            cache.isUpper = isUpper;
            cache.isLower = isLower;
            cache.isAftExt = isAftExt;
            cache.isCal = isCal;
            cache.isUncal = isUncal;
            cache.activeMask = activeMask;
            cache.cpDesignCal = cpDesignCal;
            cache.validWindward0 = validWindward0;
            cache.kCal = kCal;
            cache.Sref = Sref;
            cache.Lref = cfg.Length;
            cache.CD_friction_est = CD_friction_est;
            cache.CD_base_est = CD_base_est;
        end

        function res = evaluateAeroState(obj, alphaDeg, cache)
            %EVALUATEAEROSTATE Evaluate surface Cp, forces, and moments at angle of attack alphaDeg.
            if nargin < 3
                cache = obj.buildAeroGeometryCache();
            end
            aRad = deg2rad(alphaDeg);
            uInf = [cos(aRad), 0.0, sin(aRad)];
            liftDir = [-sin(aRad), 0.0, cos(aRad)];

            nDotU = cache.unitNormals * uInf.';
            sinImpact = max(-nDotU, 0.0);
            sinLee = max(nDotU, 0.0);

            M = cache.M;
            g = cache.g;
            pmFactor = 2.0 / (g * M^2);
            pmExp = 2.0 * g / (g - 1.0);
            CpPM = pmFactor .* (max(1.0 - 0.5 * (g - 1.0) * M .* sinLee, 0.0).^pmExp - 1.0);

            CpFace = zeros(size(cache.unitNormals, 1), 1);

            if any(cache.isCal)
                sinImpCal = sinImpact(cache.isCal);
                if abs(alphaDeg) < 1e-13
                    cpCal = cache.cpDesignCal;
                    leeMask0 = ~cache.validWindward0 & (sinLee(cache.isCal) > 0);
                    cpCalArr = CpPM(cache.isCal);
                    cpCal(leeMask0) = cpCalArr(leeMask0);
                    CpFace(cache.isCal) = cpCal;
                else
                    cpWindCal = cache.kCal .* cache.CpMax .* (sinImpCal.^2);
                    cpLeeCal = CpPM(cache.isCal);
                    CpFace(cache.isCal) = cpWindCal .* (sinImpCal > 0) + cpLeeCal .* (sinImpCal == 0);
                end
            end

            if any(cache.isUncal)
                sinImpUA = sinImpact(cache.isUncal);
                CpFace(cache.isUncal) = cache.CpMax .* (sinImpUA.^2) .* (sinImpUA > 0) + CpPM(cache.isUncal) .* (sinImpUA == 0);
            end

            dF = -CpFace(cache.activeMask) .* cache.crossVec(cache.activeMask, :);
            F_total = sum(dF, 1) ./ cache.Sref;

            CD_wave = dot(F_total, uInf);
            CL = dot(F_total, liftDir);
            CY = F_total(2);
            LD_inviscid = CL / max(CD_wave, 1e-9);

            CD_total_visc = CD_wave + cache.CD_friction_est;
            LD_viscous_est = CL / max(CD_total_visc, 1e-9);

            dM = cross(cache.rArm(cache.activeMask, :), dF, 2);
            Cm = sum(dM(:, 2)) / (cache.Sref * cache.Lref);

            res = struct();
            res.CL = CL;
            res.CD_wave = CD_wave;
            res.CY = CY;
            res.LD_inviscid = LD_inviscid;
            res.CD_friction_est = cache.CD_friction_est;
            res.CD_base_est = cache.CD_base_est;
            res.CD_total_viscous_est = CD_total_visc;
            res.LD_viscous_est = LD_viscous_est;
            res.Cm = Cm;
            res.CpMax = cache.CpMax;
            res.CpFace = CpFace;
            res.FaceNormals = cache.unitNormals;
            res.FaceAreas = cache.faceAreas;
        end

        function aero = computePreliminaryAero(obj)
            %COMPUTEPRELIMINARYAERO Evaluate design-point (alpha = 0 deg) aerodynamic metrics.
            aero = obj.evaluateAeroState(0.0);
        end
    end
end

function out = ternaryStr(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end
