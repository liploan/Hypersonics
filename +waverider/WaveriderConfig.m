classdef WaveriderConfig
%WAVERIDERCONFIG Parametric configuration for hypersonic waveriding glide bodies.
%   Supports four primary waverider generation methods:
%     1. 'conical'         - Axisymmetric cone-derived waverider (Taylor-Maccoll flow)
%     2. 'osculating_cone' - Sobieczky Osculating Cones waverider (variable curvature shock)
%     3. 'caret'           - Nonweiler Caret / 2D oblique-shock wedge waverider
%     4. 'cst'             - Parametric CST (Class-Shape Transformation) boost-glide body
%
%   Supports configurable upper surface models:
%     - 'freestream'       - Pure freestream-traced expansion-free upper surface
%     - 'canopy'           - Volumetric dorsal canopy blended with freestream leading edge
%
%   Supports an aft conformal transition frustum & conical rocket booster adapter:
%     - Smoothly lofts the non-axisymmetric waverider base at X = L into a
%       circular rocket stage interface (radius BoosterRadius) over length
%       FrustumLength, followed by an optional conical/cylindrical rocket stage.

    properties
        % Primary Method Selection ('conical', 'osculating_cone', 'caret', 'cst', or 'rv')
        Method (1,:) char {mustBeMember(Method, {'conical', 'osculating_cone', 'caret', 'cst', 'rv'})} = 'osculating_cone'

        % Flight Condition & Gas Properties
        Mach (1,1) double {mustBeGreaterThan(Mach, 1.05)} = 6.0
        Gamma (1,1) double {mustBeGreaterThan(Gamma, 1.0)} = 1.4

        % Global Planform Dimensions
        Length (1,1) double {mustBePositive} = 5.0          % Forebody waverider length L [m]
        HalfSpan (1,1) double {mustBePositive} = 1.6        % Semi-span b/2 at X = L [m]

        % Shock & Flowfield Design Parameters
        ShockAngleDeg (1,1) double {mustBeInRange(ShockAngleDeg, 8, 70)} = 14.0  % Design shock angle beta [deg]
        VertexOffsetRatio (1,1) double {mustBeInRange(VertexOffsetRatio, 0.05, 0.95)} = 0.38 % Normalized apex depth in shock cone

        % Leading-Edge Planform Parameterization: x_LE / L = (|y| / HalfSpan)^(1 / PlanformPower)
        PlanformPower (1,1) double {mustBeInRange(PlanformPower, 0.25, 2.5)} = 0.65
        LeadingEdgeAnhedralDeg (1,1) double {mustBeInRange(LeadingEdgeAnhedralDeg, -35, 35)} = 8.0

        % Osculating Cones (ICCA / FCT) Parameters
        FlatShockFraction (1,1) double {mustBeInRange(FlatShockFraction, 0.0, 0.85)} = 0.25 % Fraction of span with flat shock segment
        ShockCurvaturePower (1,1) double {mustBeInRange(ShockCurvaturePower, 1.2, 4.0)} = 2.0 % Outer shock curve exponent

        % Upper Surface Configuration
        UpperSurfaceType (1,:) char {mustBeMember(UpperSurfaceType, {'freestream', 'canopy'})} = 'canopy'
        CanopyHeightRatio (1,1) double {mustBeInRange(CanopyHeightRatio, 0.0, 0.35)} = 0.06  % Dorsal canopy height / Length
        CanopySpanPower (1,1) double {mustBeInRange(CanopySpanPower, 1.0, 6.0)} = 2.2        % Transverse bell-curve exponent

        % CST Boost-Glide Parameters (used when Method == 'cst')
        CSTNoseExponent (1,1) double {mustBeInRange(CSTNoseExponent, 0.2, 1.5)} = 0.55       % Class function N1 longitudinal
        CSTAftExponent (1,1) double {mustBeInRange(CSTAftExponent, 0.0, 1.5)} = 0.15         % Class function N2 longitudinal
        CSTUpperWeights (1,:) double {mustBeNonempty, mustBeFinite, mustBeNonnegative} = [0.12, 0.16, 0.14, 0.11] % Bernstein weights for upper ridge
        CSTLowerArchRatio (1,1) double {mustBeInRange(CSTLowerArchRatio, -0.15, 0.25)} = 0.025 % Concave compression arch on lower surface

        % Axisymmetric Reentry Vehicle (RV) Parameters (used when Method == 'rv')
        RVBodyType (1,:) char {mustBeMember(RVBodyType, {'mono_conic', 'bi_conic', 'tri_conic'})} = 'bi_conic'
        RVNoseType (1,:) char {mustBeMember(RVNoseType, {'conical', 'spherical', 'ogive'})} = 'spherical'
        RVNoseBluntnessRatio (1,1) double {mustBeInRange(RVNoseBluntnessRatio, 0.0, 0.80)} = 0.18 % R_nose / R_base (spherical) or ogive fullness
        RVConeAnglesDeg (1,3) double {mustBeFinite} = [11.5, 6.0, 12.0]                            % Stage half-angles [theta_1, theta_2, theta_3] [deg]
        RVStageFractions (1,3) double {mustBeNonnegative, mustBeFinite} = [0.42, 0.58, 0.0]        % Relative axial stage weights
        RVAftType (1,:) char {mustBeMember(RVAftType, {'flat', 'concave'})} = 'flat'               % Aft bulkhead geometry ('flat' or 'concave')
        RVAftCavityDepthRatio (1,1) double {mustBeInRange(RVAftCavityDepthRatio, 0.0, 0.85)} = 0.35 % Concave aft recess depth d_aft / R_base
        RVAftCavityPower (1,1) double {mustBeInRange(RVAftCavityPower, 1.0, 4.0)} = 2.0            % Concave aft dome exponent (2=paraboloid, 1=conical)

        % Aft Conformal Frustum Transition & Conical Rocket Booster Mount
        EnableAftFrustum (1,1) logical = false                                               % Enable conformal transition frustum aft of X = L
        FrustumLength (1,1) double {mustBeNonnegative} = 1.5                                 % Axial length of conformal transition section [m]
        BoosterRadius (1,1) double {mustBePositive} = 0.65                                   % Radius of circular rocket interface at frustum end [m]
        FrustumBlendLaw (1,:) char {mustBeMember(FrustumBlendLaw, {'cubic', 'quintic', 'linear', 'conic'})} = 'cubic'
        FrustumZOffset (1,1) double = NaN                                                    % Vertical center Z of rocket ring (NaN = auto-align with base)
        EnableRocketStage (1,1) logical = false                                              % Append conical/cylindrical rocket stage aft of frustum
        RocketStageLength (1,1) double {mustBeNonnegative} = 1.8                             % Length of attached conical rocket booster section [m]
        RocketConeAngleDeg (1,1) double {mustBeInRange(RocketConeAngleDeg, -15, 25)} = 4.5   % Conical half-angle of rocket body [deg] (0 = cylinder)

        % Discretization / Mesh Resolution
        NumStreamwise (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(NumStreamwise, 10)} = 65
        NumHalfSpan (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(NumHalfSpan, 10)} = 45
        NumBaseLayers (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(NumBaseLayers, 2)} = 12
        NumFrustumStations (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(NumFrustumStations, 4)} = 25
        NumRocketStations (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(NumRocketStations, 4)} = 18
    end

    methods
        function obj = WaveriderConfig(varargin)
            %WAVERIDERCONFIG Construct a configuration with optional name-value pairs.
            if nargin == 1 && (ischar(varargin{1}) || isstring(varargin{1}))
                obj = waverider.WaveriderConfig.fromPreset(char(varargin{1}));
                return;
            end
            if mod(numel(varargin), 2) ~= 0
                error('WaveriderConfig:OddNameValueArgs', ...
                    'Expected an even number of Name-Value arguments (got %d).', numel(varargin));
            end
            for k = 1:2:numel(varargin)
                propName = varargin{k};
                if isprop(obj, propName)
                    obj.(propName) = varargin{k+1};
                else
                    error('WaveriderConfig:UnknownProperty', ...
                        'Unknown configuration property "%s".', propName);
                end
            end
        end

        function validate(obj)
            %VALIDATE Check coupled physical constraints across configuration properties.
            if strcmpi(obj.Method, 'rv')
                if obj.RVConeAnglesDeg(1) < 1.0 || obj.RVConeAnglesDeg(1) > 55.0
                    error('WaveriderConfig:RVConeAngle', ...
                        'RVConeAnglesDeg(1) (%.2f deg) must lie in [1.0, 55.0] deg.', obj.RVConeAnglesDeg(1));
                end
                if any(obj.RVConeAnglesDeg(2:3) < -12.0 | obj.RVConeAnglesDeg(2:3) > 55.0)
                    error('WaveriderConfig:RVConeAngle', ...
                        'Aft conic half-angles RVConeAnglesDeg(2:3) must lie in [-12.0, 55.0] deg.');
                end
                nStages = 1 + strcmpi(obj.RVBodyType, 'bi_conic') + 2 * strcmpi(obj.RVBodyType, 'tri_conic');
                if sum(obj.RVStageFractions(1:nStages)) <= 0
                    error('WaveriderConfig:RVStageFractions', ...
                        'Active RVStageFractions(1:%d) must have strictly positive sum.', nStages);
                end
                return;
            end
            muDeg = asind(1.0 / obj.Mach);
            if obj.ShockAngleDeg <= muDeg + 0.15
                error('WaveriderConfig:ShockBelowMachAngle', ...
                    'ShockAngleDeg (%.2f deg) must exceed the freestream Mach wave angle asin(1/M) = %.2f deg at Mach %.2f.', ...
                    obj.ShockAngleDeg, muDeg, obj.Mach);
            end
        end
    end

    methods (Static)
        function cfg = fromPreset(presetName)
            %FROMPRESET Create a WaveriderConfig from a canonical preset name.
            cfg = waverider.WaveriderConfig();
            switch lower(strtrim(presetName))
                case {'conical', 'cone', 'taylor_maccoll'}
                    cfg.Method = 'conical';
                    cfg.Mach = 6.0;
                    cfg.Length = 5.0;
                    cfg.HalfSpan = 1.55;
                    cfg.ShockAngleDeg = 14.5;
                    cfg.VertexOffsetRatio = 0.42;
                    cfg.PlanformPower = 0.68;
                    cfg.UpperSurfaceType = 'freestream';
                    cfg.CanopyHeightRatio = 0.0;

                case {'osculating', 'osculating_cone', 'oc'}
                    cfg.Method = 'osculating_cone';
                    cfg.Mach = 7.0;
                    cfg.Length = 5.5;
                    cfg.HalfSpan = 1.75;
                    cfg.ShockAngleDeg = 13.5;
                    cfg.VertexOffsetRatio = 0.36;
                    cfg.FlatShockFraction = 0.28;
                    cfg.ShockCurvaturePower = 2.2;
                    cfg.PlanformPower = 0.62;
                    cfg.UpperSurfaceType = 'canopy';
                    cfg.CanopyHeightRatio = 0.045;

                case {'caret', 'wedge', 'nonweiler'}
                    cfg.Method = 'caret';
                    cfg.Mach = 6.5;
                    cfg.Length = 5.0;
                    cfg.HalfSpan = 1.50;
                    cfg.ShockAngleDeg = 15.0;
                    cfg.LeadingEdgeAnhedralDeg = 11.0;
                    cfg.PlanformPower = 0.85;
                    cfg.UpperSurfaceType = 'canopy';
                    cfg.CanopyHeightRatio = 0.035;

                case {'cst', 'htv2', 'boost_glide'}
                    cfg.Method = 'cst';
                    cfg.Mach = 8.0;
                    cfg.Length = 6.0;
                    cfg.HalfSpan = 1.80;
                    cfg.ShockAngleDeg = 12.5;
                    cfg.PlanformPower = 0.60;
                    cfg.LeadingEdgeAnhedralDeg = 6.0;
                    cfg.UpperSurfaceType = 'canopy';
                    cfg.CanopyHeightRatio = 0.042;
                    cfg.CSTUpperWeights = [0.14, 0.18, 0.15, 0.11];
                    cfg.CSTLowerArchRatio = 0.025;

                case {'mounted', 'booster_mounted', 'frustum_mount'}
                    cfg = waverider.WaveriderConfig.fromPreset('osculating_cone');
                    cfg.EnableAftFrustum = true;
                    cfg.FrustumLength = 1.65;
                    cfg.BoosterRadius = 0.72;
                    cfg.FrustumBlendLaw = 'cubic';
                    cfg.EnableRocketStage = true;
                    cfg.RocketStageLength = 1.85;
                    cfg.RocketConeAngleDeg = 4.5;

                case {'rv_monoconic', 'rv_mono', 'mk21'}
                    cfg.Method = 'rv';
                    cfg.RVBodyType = 'mono_conic';
                    cfg.RVNoseType = 'spherical';
                    cfg.RVNoseBluntnessRatio = 0.16;
                    cfg.Mach = 12.0;
                    cfg.Length = 2.60;
                    cfg.RVConeAnglesDeg = [8.0, 8.0, 8.0];
                    cfg.RVStageFractions = [1.0, 0.0, 0.0];
                    cfg.RVAftType = 'flat';
                    cfg.RVAftCavityDepthRatio = 0.0;

                case {'rv_biconic', 'rv_bi', 'biconic'}
                    cfg.Method = 'rv';
                    cfg.RVBodyType = 'bi_conic';
                    cfg.RVNoseType = 'spherical';
                    cfg.RVNoseBluntnessRatio = 0.22;
                    cfg.Mach = 14.0;
                    cfg.Length = 3.10;
                    cfg.RVConeAnglesDeg = [11.5, 6.0, 6.0];
                    cfg.RVStageFractions = [0.42, 0.58, 0.0];
                    cfg.RVAftType = 'concave';
                    cfg.RVAftCavityDepthRatio = 0.38;
                    cfg.RVAftCavityPower = 2.0;

                case {'rv_triconic', 'rv_triconic_flare', 'triconic'}
                    cfg.Method = 'rv';
                    cfg.RVBodyType = 'tri_conic';
                    cfg.RVNoseType = 'ogive';
                    cfg.RVNoseBluntnessRatio = 0.25;
                    cfg.Mach = 10.0;
                    cfg.Length = 3.40;
                    cfg.RVConeAnglesDeg = [9.5, 3.5, 13.5];
                    cfg.RVStageFractions = [0.36, 0.34, 0.30];
                    cfg.RVAftType = 'concave';
                    cfg.RVAftCavityDepthRatio = 0.42;
                    cfg.RVAftCavityPower = 2.6;

                case {'rv_sharp_conic', 'rv_sharp'}
                    cfg.Method = 'rv';
                    cfg.RVBodyType = 'mono_conic';
                    cfg.RVNoseType = 'conical';
                    cfg.RVNoseBluntnessRatio = 0.0;
                    cfg.Mach = 15.0;
                    cfg.Length = 2.80;
                    cfg.RVConeAnglesDeg = [7.5, 7.5, 7.5];
                    cfg.RVStageFractions = [1.0, 0.0, 0.0];
                    cfg.RVAftType = 'flat';
                    cfg.RVAftCavityDepthRatio = 0.0;

                otherwise
                    error('WaveriderConfig:UnknownPreset', ...
                        ['Unknown preset "%s". Valid presets: conical, osculating_cone, caret, cst, ', ...
                         'booster_mounted, rv_monoconic, rv_biconic, rv_triconic_flare, rv_sharp_conic.'], presetName);
            end
        end
    end
end
