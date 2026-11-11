function sweepTable = demo_parametric_sweep()
%DEMO_PARAMETRIC_SWEEP Parametric trade study across shock angle and canopy height.
%   sweepTable = DEMO_PARAMETRIC_SWEEP() evaluates how shock angle beta and
%   upper canopy height ratio affect volumetric efficiency eta_V, enclosed
%   volume, and inviscid L/D while maintaining 100% watertightness.

    rootDir = fileparts(fileparts(mfilename('fullpath')));
    addpath(rootDir);

    betaValues = [12.0, 14.0, 16.0, 18.0];
    canopyValues = [0.00, 0.04, 0.08];

    rows = [];
    fprintf('\n================================================================================\n');
    fprintf('  PARAMETRIC DESIGN SWEEP (Osculating Cones, M = 7.0, L = 5.5 m)\n');
    fprintf('================================================================================\n');
    fprintf('  %-10s %-12s %-12s %-12s %-10s %-10s %-12s\n', ...
        'Beta[deg]', 'CanopyRatio', 'Volume[m^3]', 'S_plan[m^2]', 'eta_V', 'L/D_inv', 'Watertight');
    fprintf('--------------------------------------------------------------------------------\n');

    for b = 1:numel(betaValues)
        for c = 1:numel(canopyValues)
            cfg = waverider.WaveriderConfig.fromPreset('osculating_cone');
            cfg.ShockAngleDeg = betaValues(b);
            cfg.CanopyHeightRatio = canopyValues(c);
            if canopyValues(c) == 0
                cfg.UpperSurfaceType = 'freestream';
            else
                cfg.UpperSurfaceType = 'canopy';
            end

            body = waverider.WaveriderBody(cfg);
            rep = body.Report;
            aero = body.AeroMetrics;

            fprintf('  %-10.1f %-12.3f %-12.3f %-12.3f %-10.4f %-10.2f %-12d\n', ...
                cfg.ShockAngleDeg, cfg.CanopyHeightRatio, rep.EnclosedVolume, ...
                rep.PlanformArea, rep.VolumetricEfficiency, aero.LD_inviscid, rep.IsWatertight);

            rows = [rows; cfg.ShockAngleDeg, cfg.CanopyHeightRatio, ...
                    rep.EnclosedVolume, rep.PlanformArea, rep.VolumetricEfficiency, ...
                    aero.CL, aero.CD_wave, aero.LD_inviscid, double(rep.IsWatertight)]; %#ok<AGROW>
        end
    end

    sweepTable = array2table(rows, 'VariableNames', ...
        {'ShockAngleDeg', 'CanopyHeightRatio', 'EnclosedVolume_m3', 'PlanformArea_m2', ...
         'VolumetricEfficiency', 'CL', 'CD_wave', 'LD_inviscid', 'IsWatertight'});
    fprintf('================================================================================\n');
end
