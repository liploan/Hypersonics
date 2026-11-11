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
            switch lower(cfg.Method)
                case {'conical', 'osculating_cone'}
                    obj.SurfaceData = waverider.geometry.generateOsculatingCone(cfg);
                case 'caret'
                    obj.SurfaceData = waverider.geometry.generateCaretWedge(cfg);
                case 'cst'
                    obj.SurfaceData = waverider.geometry.generateCSTBody(cfg);
                otherwise
                    error('WaveriderBody:UnsupportedMethod', ...
                        'Unsupported method "%s".', cfg.Method);
            end

            % Append optional conformal transition frustum & conical rocket stage
            obj.SurfaceData = waverider.geometry.appendConformalFrustum(obj.SurfaceData, cfg);

            obj.MeshData = waverider.geometry.stitchWatertightMesh(obj.SurfaceData, cfg);
            obj.Report = waverider.mesh.verifyWatertight(obj.MeshData, cfg);
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
            %COMPUTERCS Compute monostatic RCS (UHF to X-band) with concave multi-bounce.
            rcsData = waverider.rcs.computeRCS(obj, freqGHz, azDeg, elDeg, varargin{:});
        end

        function [fig, summaryMetrics] = plotRCS(obj, varargin)
            %PLOTRCS Render 4-panel UHF-to-X-band RCS & concave underside scattering dashboard.
            [fig, summaryMetrics] = waverider.viz.plotRCSAnalysis(obj, varargin{:});
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
            %COMPUTEAEROPOLAR Compute CL, CD_wave, L/D, and Cm across angle of attack alpha [deg].
            if nargin < 2 || isempty(alphaDegVec)
                alphaDegVec = linspace(-4.0, 16.0, 41);
            end
            alphaDegVec = alphaDegVec(:);
            nA = numel(alphaDegVec);

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
            isAftExt = (obj.MeshData.FacePatchID == 4) | (obj.MeshData.FacePatchID == 5);
            activeMask = isUpper | isLower | isAftExt;

            Sref = max(obj.Report.PlanformArea, 1e-12);
            Lref = cfg.Length;
            xCG = obj.Report.VolumeCentroid;
            rArm = faceCenters - xCG;

            NxWave = obj.SurfaceData.NxWaverider;
            meanFlowCp0 = mean( reshape(obj.SurfaceData.Cp_low(1:NxWave, :), [], 1) );
            sinDelta0 = max(mean(-unitNormals(isLower, 1)), 1e-4);

            CL = zeros(nA, 1);
            CD = zeros(nA, 1);
            LD = zeros(nA, 1);
            Cm = zeros(nA, 1);

            for k = 1:nA
                aRad = deg2rad(alphaDegVec(k));
                uInf = [cos(aRad), 0.0, sin(aRad)];
                liftDir = [-sin(aRad), 0.0, cos(aRad)];

                nDotU = unitNormals * uInf.';
                sinImpact = max(-nDotU, 0.0);
                sinLee = max(nDotU, 0.0);

                CpFace = zeros(size(F, 1), 1);
                impactScale = (sinImpact(isLower) ./ sinDelta0).^2;
                CpFace(isLower) = (0.55 * meanFlowCp0 .* min(impactScale, 4.0) + ...
                                   0.45 * CpMax .* (sinImpact(isLower).^2)) .* (sinImpact(isLower) > 0) - ...
                                  min(1.0 / (M^2), (1.5 / M) .* sinLee(isLower)) .* (sinImpact(isLower) == 0);

                CpExpansionUp = -min(1.2 / (M^2), (2.0 / M) .* sinLee(isUpper));
                CpFace(isUpper) = CpMax .* (sinImpact(isUpper).^2) + CpExpansionUp .* (sinImpact(isUpper) == 0);

                if any(isAftExt)
                    CpExpAft = -min(1.0 / (M^2), (1.5 / M) .* sinLee(isAftExt));
                    CpFace(isAftExt) = CpMax .* (sinImpact(isAftExt).^2) + CpExpAft .* (sinImpact(isAftExt) == 0);
                end

                dF = -CpFace(activeMask) .* crossVec(activeMask, :);
                F_total = sum(dF, 1) ./ Sref;

                CD(k) = dot(F_total, uInf);
                CL(k) = dot(F_total, liftDir);
                LD(k) = CL(k) / max(CD(k), 1e-6);

                dM = cross(rArm(activeMask, :), dF, 2);
                Cm(k) = sum(dM(:, 2)) / (Sref * Lref);
            end

            polar = struct();
            polar.AlphaDeg = alphaDegVec;
            polar.CL = CL;
            polar.CD_wave = CD;
            polar.LD_inviscid = LD;
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

            for s = 1:nStations
                xs = xStations(s);
                if xs <= L + 1e-7
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

            % Spanwise distributions on waverider forebody (X <= L)
            spanY = yLE(:);
            chordY = max(L - xLE(:), 0.0);
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
            sec.SpanEta = spanY ./ cfg.HalfSpan;
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
            fprintf('  WAVERIDER GLIDE BODY SUMMARY (%s)\n', upper(cfg.Method));
            fprintf('----------------------------------------------------------------\n');
            fprintf('  Design Mach / Shock Angle : M = %.2f, beta = %.2f deg\n', cfg.Mach, cfg.ShockAngleDeg);
            fprintf('  Forebody Dimensions (Lxb) : %.2f m x %.2f m (HalfSpan = %.2f m)\n', ...
                cfg.Length, 2.0 * cfg.HalfSpan, cfg.HalfSpan);
            fprintf('  Upper Surface Mode        : %s (Canopy ratio = %.3f)\n', ...
                cfg.UpperSurfaceType, cfg.CanopyHeightRatio);
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
            fprintf('----------------------------------------------------------------\n');
        end
    end

    methods (Access = private)
        function aero = computePreliminaryAero(obj)
            %COMPUTEPRELIMINARYAERO Integrate pressure forces over surface patches.
            cfg = obj.Config;
            TR = obj.MeshData.TR;
            V = TR.Points;
            F = TR.ConnectivityList;

            v1 = V(F(:, 1), :);
            v2 = V(F(:, 2), :);
            v3 = V(F(:, 3), :);

            crossVec = 0.5 .* cross(v2 - v1, v3 - v1, 2);
            faceAreas = sqrt(sum(crossVec.^2, 2)) + 1e-15;
            unitNormals = crossVec ./ faceAreas;

            M = cfg.Mach;
            g = cfg.Gamma;
            CpMax = (2.0 / (g * M^2)) * ...
                (((g + 1.0)^2 * M^2 / (4.0 * g * M^2 - 2.0 * (g - 1.0)))^(g / (g - 1.0)) * ...
                 ((1.0 - g + 2.0 * g * M^2) / (g + 1.0)) - 1.0);

            sinImpact = max(-unitNormals(:, 1), 0.0);

            CpFace = zeros(size(F, 1), 1);
            isUpper  = (obj.MeshData.FacePatchID == 1);
            isLower  = (obj.MeshData.FacePatchID == 2);
            isAftExt = (obj.MeshData.FacePatchID == 4) | (obj.MeshData.FacePatchID == 5);

            NxWave = obj.SurfaceData.NxWaverider;
            meanFlowCp = mean( reshape(obj.SurfaceData.Cp_low(1:NxWave, :), [], 1) );
            CpFace(isLower) = 0.65 * meanFlowCp + 0.35 * CpMax .* (sinImpact(isLower).^2);

            sinLee = max(unitNormals(isUpper, 1), 0.0);
            CpExpansion = -min(1.2 / (M^2), (2.0 / M) .* sinLee);
            CpFace(isUpper) = CpMax .* (sinImpact(isUpper).^2) + CpExpansion .* (sinImpact(isUpper) == 0);

            if any(isAftExt)
                sinLeeAft = max(unitNormals(isAftExt, 1), 0.0);
                CpExpAft = -min(1.0 / (M^2), (1.5 / M) .* sinLeeAft);
                CpFace(isAftExt) = CpMax .* (sinImpact(isAftExt).^2) + CpExpAft .* (sinImpact(isAftExt) == 0);
            end

            Sref = max(obj.Report.PlanformArea, 1e-12);
            activeMask = isUpper | isLower | isAftExt;
            F_vec = sum(-CpFace(activeMask) .* crossVec(activeMask, :), 1) ./ Sref;

            CD_wave = F_vec(1);
            CL = F_vec(3);

            aero = struct();
            aero.CL = CL;
            aero.CD_wave = CD_wave;
            aero.LD_inviscid = CL / max(CD_wave, 1e-9);
            aero.CpMax = CpMax;
            aero.CpFace = CpFace;
            aero.FaceNormals = unitNormals;
            aero.FaceAreas = faceAreas;
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
