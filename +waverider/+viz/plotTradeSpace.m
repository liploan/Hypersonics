function [fig, tradeData] = plotTradeSpace(basePreset, varargin)
%PLOTTRADESPACE Render 2D design trade-space contours and Pareto front.
%   [fig, tradeData] = PLOTTRADESPACE(basePreset, Name, Value) sweeps
%   ShockAngleDeg and CanopyHeightRatio (plus PlanformPower) and renders a
%   4-panel design trade-space figure:
%     1. Contour map of Inviscid L/D overlaid with Volumetric Efficiency eta_V iso-lines
%     2. Contour map of Enclosed Volume [m^3] vs (Shock Angle, Canopy Ratio)
%     3. Pareto Trade-Off Scatter: Inviscid L/D vs Volumetric Efficiency eta_V
%     4. Planform Power Exponent P sensitivity on L/D, Volume, and Aspect Ratio

    if nargin < 1 || isempty(basePreset)
        basePreset = 'osculating_cone';
    end

    p = inputParser;
    addParameter(p, 'Visible', 'on', @(x) ischar(x) || isstring(x));
    addParameter(p, 'BetaGrid', linspace(11.5, 18.5, 8), @isnumeric);
    addParameter(p, 'CanopyGrid', linspace(0.0, 0.08, 7), @isnumeric);
    parse(p, varargin{:});
    opts = p.Results;

    betaVec = opts.BetaGrid;
    canopyVec = opts.CanopyGrid;
    nB = numel(betaVec);
    nC = numel(canopyVec);

    [BETA, CANOPY] = meshgrid(betaVec, canopyVec);
    LD_grid = zeros(nC, nB);
    ETAV_grid = zeros(nC, nB);
    VOL_grid = zeros(nC, nB);

    baseCfg = waverider.WaveriderConfig.fromPreset(basePreset);
    % Use fast grid resolution for smooth contour sweep
    baseCfg.NumStreamwise = 35;
    baseCfg.NumHalfSpan = 25;

    for i = 1:nC
        for j = 1:nB
            cfg = baseCfg;
            cfg.ShockAngleDeg = BETA(i, j);
            cfg.CanopyHeightRatio = CANOPY(i, j);
            if cfg.CanopyHeightRatio <= 1e-6
                cfg.UpperSurfaceType = 'freestream';
            else
                cfg.UpperSurfaceType = 'canopy';
            end
            b = waverider.WaveriderBody(cfg);
            LD_grid(i, j) = b.AeroMetrics.LD_inviscid;
            ETAV_grid(i, j) = b.Report.VolumetricEfficiency;
            VOL_grid(i, j) = b.Report.EnclosedVolume;
        end
    end

    % 1D sweep over PlanformPower P
    pVec = linspace(0.35, 1.25, 12);
    LD_p = zeros(size(pVec));
    ETAV_p = zeros(size(pVec));
    AR_p = zeros(size(pVec));
    for k = 1:numel(pVec)
        cfg = baseCfg;
        cfg.PlanformPower = pVec(k);
        b = waverider.WaveriderBody(cfg);
        LD_p(k) = b.AeroMetrics.LD_inviscid;
        ETAV_p(k) = b.Report.VolumetricEfficiency;
        AR_p(k) = b.Report.AspectRatio;
    end

    fig = figure('Color', 'w', 'Position', [70, 70, 1400, 880], 'Visible', opts.Visible);
    waverider.viz.stylePublicationFigure(fig, 'init');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
    title(tl, sprintf('Parametric Design Space & Pareto Trade-Offs [%s | M = %.1f, L = %.1f m]', ...
        upper(strrep(baseCfg.Method, '_', ' ')), baseCfg.Mach, baseCfg.Length), ...
        'FontSize', 14, 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Tile 1: Inviscid L/D Contour Map + Volumetric Efficiency Iso-Lines
    % ---------------------------------------------------------------------
    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    contourf(ax1, BETA, CANOPY, LD_grid, 14, 'LineColor', 'none', 'HandleVisibility', 'off');
    colormap(ax1, 'turbo');
    cb1 = colorbar(ax1);
    cb1.Label.String = 'Inviscid L/D (\alpha = 0^\circ)';

    [C_eta, h_eta] = contour(ax1, BETA, CANOPY, ETAV_grid, 7, '--', ...
        'LineColor', [0.06, 0.09, 0.16], 'LineWidth', 1.6, 'DisplayName', '\eta_V Iso-Contours');
    clabel(C_eta, h_eta, 'Color', [0.06, 0.09, 0.16], 'FontSize', 8.5, 'FontWeight', 'bold');
    plot(ax1, baseCfg.ShockAngleDeg, baseCfg.CanopyHeightRatio, 'kp', ...
        'MarkerSize', 13, 'MarkerFaceColor', 'y', 'DisplayName', 'Baseline Preset');

    xlabel(ax1, 'Design Shock Angle \beta [deg]');
    ylabel(ax1, 'Canopy Height Ratio h_{canopy} / L [-]');
    title(ax1, 'Inviscid L/D Map with Dashed \eta_V Iso-Contours', 'FontWeight', 'bold');
    legend(ax1, 'Location', 'northeast', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 2: Enclosed Internal Volume Contour Map [m^3]
    % ---------------------------------------------------------------------
    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    [C_vol, h_vol] = contourf(ax2, BETA, CANOPY, VOL_grid, 14, 'LineColor', [0.15, 0.20, 0.28]);
    clabel(C_vol, h_vol, 'Color', [0.06, 0.09, 0.16], 'FontSize', 8.5, 'FontWeight', 'bold');
    cb2 = colorbar(ax2);
    cb2.Label.String = 'Enclosed Volume V [m^3]';
    plot(ax2, baseCfg.ShockAngleDeg, baseCfg.CanopyHeightRatio, 'kp', ...
        'MarkerSize', 13, 'MarkerFaceColor', 'y', 'DisplayName', 'Baseline Preset');

    xlabel(ax2, 'Design Shock Angle \beta [deg]');
    ylabel(ax2, 'Canopy Height Ratio h_{canopy} / L [-]');
    title(ax2, 'Enclosed Packaging Volume V(\beta, h_{canopy}) [m^3]', 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Tile 3: Pareto Front Scatter (L/D vs Volumetric Efficiency eta_V)
    % ---------------------------------------------------------------------
    ax3 = nexttile(tl, 3);
    hold(ax3, 'on');
    scatter(ax3, ETAV_grid(:), LD_grid(:), 55, BETA(:), 'filled', 'MarkerEdgeColor', 'k');
    cb3 = colorbar(ax3);
    cb3.Label.String = 'Shock Angle \beta [deg]';
    grid(ax3, 'on');
    xlabel(ax3, 'Volumetric Efficiency \eta_V = V^{2/3} / S_{plan} [-]');
    ylabel(ax3, 'Inviscid Lift-to-Drag Ratio L/D [-]');
    title(ax3, 'Aerodynamic vs. Volumetric Packaging Pareto Front', 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Tile 4: Planform Exponent Sensitivity (PlanformPower P)
    % ---------------------------------------------------------------------
    ax4 = nexttile(tl, 4);
    yyaxis(ax4, 'left');
    hold(ax4, 'on');
    plot(ax4, pVec, LD_p, '-o', 'Color', [0.12, 0.46, 0.80], 'LineWidth', 2.2, ...
        'MarkerSize', 5, 'DisplayName', 'Inviscid L/D');
    plot(ax4, pVec, AR_p, '--s', 'Color', [0.20, 0.65, 0.35], 'LineWidth', 2.0, ...
        'MarkerSize', 5, 'DisplayName', 'Aspect Ratio AR');
    ylabel(ax4, 'Inviscid L/D & Aspect Ratio [-]');

    yyaxis(ax4, 'right');
    plot(ax4, pVec, ETAV_p, '-^', 'Color', [0.86, 0.32, 0.14], 'LineWidth', 2.2, ...
        'MarkerSize', 5, 'DisplayName', 'Volumetric Efficiency \eta_V');
    ylabel(ax4, 'Volumetric Efficiency \eta_V [-]');

    grid(ax4, 'on');
    xlabel(ax4, 'Planform Power-Law Exponent P (y/b = (x/L)^P) [-]');
    title(ax4, 'Planform Leading-Edge Exponent Sensitivity', 'FontWeight', 'bold');
    legend(ax4, 'Location', 'best', 'FontSize', 8);

    waverider.viz.stylePublicationFigure(fig, 'finalize');

    tradeData = struct('BETA', BETA, 'CANOPY', CANOPY, ...
        'LD', LD_grid, 'EtaV', ETAV_grid, 'Volume', VOL_grid, ...
        'PlanformPowerVec', pVec, 'LD_p', LD_p, 'EtaV_p', ETAV_p, 'AR_p', AR_p);
end
