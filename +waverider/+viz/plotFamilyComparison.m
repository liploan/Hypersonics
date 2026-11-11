function fig = plotFamilyComparison(bodies, varargin)
%PLOTFAMILYCOMPARISON Compare geometry and key attributes across waverider families.
%   fig = PLOTFAMILYCOMPARISON(bodies, Name, Value) generates a 4-panel
%   comparative analysis across a cell array of WaveriderBody objects
%   (defaults to {'conical', 'osculating_cone', 'caret', 'cst'}).

    p = inputParser;
    addParameter(p, 'Visible', 'on', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    opts = p.Results;

    if nargin < 1 || isempty(bodies)
        presets = {'conical', 'osculating_cone', 'caret', 'cst'};
        bodies = cell(1, numel(presets));
        for k = 1:numel(presets)
            bodies{k} = waverider.WaveriderBody(presets{k});
        end
    end

    nB = numel(bodies);
    colors = [0.12, 0.46, 0.78;
              0.88, 0.36, 0.14;
              0.18, 0.66, 0.34;
              0.58, 0.26, 0.72];
    if nB > size(colors, 1)
        colors = lines(nB);
    end

    fig = figure('Color', 'w', 'Position', [60, 60, 1360, 840], 'Visible', opts.Visible);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Comparative Geometry & Key Attributes Across Hypersonic Waverider Families', ...
        'FontSize', 14, 'FontWeight', 'bold');

    labels = cell(1, nB);
    for k = 1:nB
        labels{k} = upper(strrep(bodies{k}.Config.Method, '_', ' '));
    end

    % ---------------------------------------------------------------------
    % Tile 1: Normalized Aft Base Cross-Section Comparison (Y/b_half vs Z/b_half)
    % ---------------------------------------------------------------------
    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    for k = 1:nB
        b = bodies{k};
        bHalf = b.Config.HalfSpan;
        yUp = b.SurfaceData.Y_up(end, :) ./ bHalf;
        zUp = b.SurfaceData.Z_up(end, :) ./ bHalf;
        yLow = b.SurfaceData.Y_low(end, :) ./ bHalf;
        zLow = b.SurfaceData.Z_low(end, :) ./ bHalf;

        yLoop = [yUp, fliplr(yLow), yUp(1)];
        zLoop = [zUp, fliplr(zLow), zUp(1)];
        plot(ax1, yLoop, zLoop, '-', 'Color', colors(k, :), 'LineWidth', 2.2, ...
            'DisplayName', sprintf('%s (\\eta_V = %.3f)', labels{k}, b.Report.VolumetricEfficiency));
    end
    axis(ax1, 'equal');
    grid(ax1, 'on');
    xlabel(ax1, 'Normalized Span Y / (b/2) [-]');
    ylabel(ax1, 'Normalized Vertical Z / (b/2) [-]');
    title(ax1, 'Normalized Aft Base Cross-Sections at X / L = 1.0', 'FontWeight', 'bold');
    legend(ax1, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 2: Normalized Longitudinal Area Progression A(x) / S_plan
    % ---------------------------------------------------------------------
    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    for k = 1:nB
        b = bodies{k};
        sec = b.SectionalAttributes;
        plot(ax2, sec.xNorm, sec.AreaX ./ b.Report.PlanformArea, '-', ...
            'Color', colors(k, :), 'LineWidth', 2.3, ...
            'DisplayName', sprintf('%s (Vol = %.2f m^3)', labels{k}, b.Report.EnclosedVolume));
    end
    grid(ax2, 'on');
    xlabel(ax2, 'Normalized Axial Coordinate x / L [-]');
    ylabel(ax2, 'Normalized Cross-Section Area A(x) / S_{plan} [-]');
    title(ax2, 'Longitudinal Cross-Sectional Area Growth', 'FontWeight', 'bold');
    legend(ax2, 'Location', 'northwest', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 3: Angle-of-Attack Inviscid Lift-to-Drag Polar Comparison
    % ---------------------------------------------------------------------
    ax3 = nexttile(tl, 3);
    hold(ax3, 'on');
    alphaVec = linspace(-2.0, 14.0, 33);
    peakLD = zeros(1, nB);
    for k = 1:nB
        b = bodies{k};
        pol = b.computeAeroPolar(alphaVec);
        peakLD(k) = max(pol.LD_inviscid);
        plot(ax3, pol.AlphaDeg, pol.LD_inviscid, '-', 'Color', colors(k, :), 'LineWidth', 2.3, ...
            'DisplayName', sprintf('%s (M=%.1f, (L/D)_{max}=%.2f)', labels{k}, b.Config.Mach, peakLD(k)));
    end
    xline(ax3, 0.0, 'k:', 'Design \alpha = 0^\circ', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    grid(ax3, 'on');
    xlabel(ax3, 'Angle of Attack \alpha [deg]');
    ylabel(ax3, 'Inviscid Lift-to-Drag Ratio L/D [-]');
    title(ax3, 'Lift-to-Drag Ratio Polar Across Angle of Attack', 'FontWeight', 'bold');
    legend(ax3, 'Location', 'northeast', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 4: Multi-Attribute Bar Chart (Volumetric Efficiency, L/D, Wetted Ratio)
    % ---------------------------------------------------------------------
    ax4 = nexttile(tl, 4);
    etaV_scaled = zeros(nB, 1);
    LD_zero = zeros(nB, 1);
    LD_max = zeros(nB, 1);
    wetRatio = zeros(nB, 1);

    for k = 1:nB
        b = bodies{k};
        etaV_scaled(k) = 30.0 * b.Report.VolumetricEfficiency; % Scaled x30 for visual comparison
        LD_zero(k) = b.AeroMetrics.LD_inviscid;
        LD_max(k) = peakLD(k);
        wetRatio(k) = b.Report.WettedArea / b.Report.PlanformArea;
    end

    barMatrix = [etaV_scaled, LD_zero, LD_max, wetRatio];
    bChart = bar(ax4, categorical(labels, labels), barMatrix, 'grouped');
    bChart(1).FaceColor = [0.20, 0.55, 0.85];
    bChart(2).FaceColor = [0.90, 0.42, 0.18];
    bChart(3).FaceColor = [0.95, 0.76, 0.15];
    bChart(4).FaceColor = [0.28, 0.70, 0.42];

    grid(ax4, 'on');
    ylabel(ax4, 'Metric Value [-]');
    title(ax4, 'Key Performance & Geometric Attribute Comparison', 'FontWeight', 'bold');
    legend(ax4, {'30 \times \eta_V (Vol. Eff.)', '(L/D)_{\alpha=0^\circ}', '(L/D)_{max}', 'S_{wet} / S_{plan}'}, ...
        'Location', 'northeast', 'FontSize', 8);
end
