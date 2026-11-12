function fig = plotBoosterIntegration(body, varargin)
%PLOTBOOSTERINTEGRATION Visualize conformal aft transition frustum & rocket booster mount.
%   fig = PLOTBOOSTERINTEGRATION(body, Name, Value) renders a 4-panel
%   engineering analysis of the conformal transition section connecting the
%   waverider aft base (X = L) to a traditional conical rocket body:
%     1. 3D Integrated Stack (Waverider + Conformal Frustum + Conical Booster)
%     2. Cross-Sectional Morphing Slices along the Frustum (tau = 0 -> 1)
%     3. Comparison of Frustum Blending Laws (Linear, Cubic, Quintic, Conic)
%     4. Full-Stack Longitudinal Area A(x), Half-Span Y_max(x), & Volume V(x)

    p = inputParser;
    addParameter(p, 'Visible', 'on', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    opts = p.Results;

    % Ensure conformal frustum is enabled on the body for this visualization
    cfg = body.Config;
    if ~cfg.EnableAftFrustum || cfg.FrustumLength <= 0
        cfg.EnableAftFrustum = true;
        cfg.FrustumLength = 0.30 * cfg.Length;
        cfg.BoosterRadius = 0.42 * cfg.HalfSpan;
        cfg.EnableRocketStage = true;
        cfg.RocketStageLength = 0.35 * cfg.Length;
        body = waverider.WaveriderBody(cfg);
    end

    surfData = body.SurfaceData;
    meshData = body.MeshData;
    rep = body.Report;
    sec = body.SectionalAttributes;

    L = cfg.Length;
    Lf = cfg.FrustumLength;
    NxWave = surfData.NxWaverider;
    NxFrust = surfData.NxFrustumEnd;
    NxTot = surfData.NxTotal;

    fig = figure('Color', 'w', 'Position', [65, 65, 1400, 880], 'Visible', opts.Visible);
    waverider.viz.stylePublicationFigure(fig, 'init');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'loose', 'Padding', 'compact');

    title(tl, sprintf('Conformal Aft Transition Frustum & Conical Rocket Booster Mount [%s | R_{booster} = %.2f m]', ...
        upper(strrep(cfg.Method, '_', ' ')), cfg.BoosterRadius), ...
        'FontSize', 14, 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Tile 1: 3D Integrated Waverider + Conformal Frustum + Conical Booster
    % ---------------------------------------------------------------------
    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    V = meshData.Vertices;

    patch(ax1, 'Faces', meshData.FacesUpper, 'Vertices', V, ...
        'FaceColor', [0.18, 0.50, 0.82], 'EdgeColor', 'none', 'FaceAlpha', 0.94, ...
        'DisplayName', 'Waverider Upper Surface');
    patch(ax1, 'Faces', meshData.FacesLower, 'Vertices', V, ...
        'FaceColor', [0.88, 0.42, 0.20], 'EdgeColor', 'none', 'FaceAlpha', 0.94, ...
        'DisplayName', 'Waverider Lower Surface');
    patch(ax1, 'Faces', meshData.FacesFrustum, 'Vertices', V, ...
        'FaceColor', [0.58, 0.64, 0.72], 'EdgeColor', 'none', 'FaceAlpha', 0.92, ...
        'DisplayName', sprintf('Conformal Frustum (L_f=%.2fm)', Lf));
    if ~isempty(meshData.FacesRocket)
        patch(ax1, 'Faces', meshData.FacesRocket, 'Vertices', V, ...
            'FaceColor', [0.84, 0.85, 0.88], 'EdgeColor', 'none', 'FaceAlpha', 0.95, ...
            'DisplayName', sprintf('Conical Rocket Body (\\theta_r=%.1f^\\circ)', cfg.RocketConeAngleDeg));
    end
    patch(ax1, 'Faces', meshData.FacesBase, 'Vertices', V, ...
        'FaceColor', [0.32, 0.70, 0.44], 'EdgeColor', 'none', 'FaceAlpha', 0.95, ...
        'DisplayName', 'Watertight Booster Base');

    % Highlight morphing station rings along the conformal frustum
    ringRows = unique(round(linspace(NxWave, NxFrust, 5)));
    ringColors = turbo(numel(ringRows));
    for rIdx = 1:numel(ringRows)
        r = ringRows(rIdx);
        xR = [surfData.X_up(r, :), fliplr(surfData.X_low(r, :))];
        yR = [surfData.Y_up(r, :), fliplr(surfData.Y_low(r, :))];
        zR = [surfData.Z_up(r, :), fliplr(surfData.Z_low(r, :))];
        plot3(ax1, xR, yR, zR, '-', 'Color', ringColors(rIdx, :), 'LineWidth', 2.0, ...
            'HandleVisibility', 'off');
    end

    axis(ax1, 'equal');
    grid(ax1, 'on');
    view(ax1, [-134, 22]);
    camlight(ax1, 'headlight');
    lighting(ax1, 'gouraud');
    xlabel(ax1, 'Axial X [m]'); ylabel(ax1, 'Span Y [m]'); zlabel(ax1, 'Vertical Z [m]');
    title(ax1, sprintf('3D Integrated Stack (%d Faces, FreeBoundary = %d, \\chi = %d)', ...
        rep.NumFaces, rep.NumBoundaryEdges, rep.EulerCharacteristic), 'FontWeight', 'bold');
    legend(ax1, 'Location', 'southwest', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 2: Progressive Cross-Sectional Morphing along Conformal Frustum
    % ---------------------------------------------------------------------
    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    tauList = [0.00, 0.25, 0.50, 0.75, 1.00];
    for idx = 1:numel(ringRows)
        r = ringRows(idx);
        yR = [surfData.Y_up(r, :), fliplr(surfData.Y_low(r, :)), surfData.Y_up(r, 1)];
        zR = [surfData.Z_up(r, :), fliplr(surfData.Z_low(r, :)), surfData.Z_up(r, 1)];
        lbl = sprintf('\\tau = %.2f (X = %.2f m)', tauList(min(idx, end)), surfData.X_up(r, 1));
        plot(ax2, yR, zR, '-', 'Color', ringColors(idx, :), 'LineWidth', 2.3, 'DisplayName', lbl);
    end
    if surfData.HasRocketStage
        yEnd = [surfData.Y_up(NxTot, :), fliplr(surfData.Y_low(NxTot, :)), surfData.Y_up(NxTot, 1)];
        zEnd = [surfData.Z_up(NxTot, :), fliplr(surfData.Z_low(NxTot, :)), surfData.Z_up(NxTot, 1)];
        plot(ax2, yEnd, zEnd, 'k--', 'LineWidth', 2.0, ...
            'DisplayName', sprintf('Rocket Stage Exit (R = %.2f m)', surfData.RocketEndRadius));
    end
    plot(ax2, 0.0, surfData.BoosterZCenter, 'kp', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'y', 'DisplayName', 'Booster Axis Center');

    axis(ax2, 'equal');
    grid(ax2, 'on');
    xlabel(ax2, 'Spanwise Coordinate Y [m]');
    ylabel(ax2, 'Vertical Coordinate Z [m]');
    title(ax2, 'Conformal Frustum Cross-Sectional Morphing (\tau = 0 \rightarrow 1)', 'FontWeight', 'bold');
    legend(ax2, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 7.5);

    % ---------------------------------------------------------------------
    % Tile 3: Comparison of Frustum Blending Laws (Span Boat-Tail & Dorsal Loft)
    % ---------------------------------------------------------------------
    ax3 = nexttile(tl, 3);
    hold(ax3, 'on');
    laws = {'linear', 'cubic', 'quintic', 'conic'};
    lawColors = lines(numel(laws));
    tauFine = linspace(0.0, 1.0, 60);
    xFine = L + tauFine .* Lf;

    yTip0 = abs(surfData.Y_up(NxWave, end));
    yTip1 = cfg.BoosterRadius;
    midCol = cfg.NumHalfSpan;
    zTop0 = surfData.Z_up(NxWave, midCol);
    zTop1 = surfData.BoosterZCenter + cfg.BoosterRadius;

    for m = 1:numel(laws)
        switch laws{m}
            case 'linear'
                sig = tauFine;
            case 'cubic'
                sig = 3.0 .* (tauFine.^2) - 2.0 .* (tauFine.^3);
            case 'quintic'
                sig = 10.0 .* (tauFine.^3) - 15.0 .* (tauFine.^4) + 6.0 .* (tauFine.^5);
            case 'conic'
                pE = 1.45;
                sig = (tauFine.^pE) ./ ((tauFine.^pE) + ((1.0 - tauFine).^pE));
        end
        yBoat = (1.0 - sig) .* yTip0 + sig .* yTip1;
        zCrown = (1.0 - sig) .* zTop0 + sig .* zTop1;

        plot(ax3, xFine, yBoat, '-', 'Color', lawColors(m, :), 'LineWidth', 2.2, ...
            'DisplayName', sprintf('%s: Half-Span Y(x)', upper(laws{m})));
        plot(ax3, xFine, zCrown, '--', 'Color', lawColors(m, :), 'LineWidth', 1.8, ...
            'HandleVisibility', 'off');
    end
    grid(ax3, 'on');
    xlabel(ax3, 'Frustum Axial Coordinate X [m]');
    ylabel(ax3, 'Half-Span Y(x) [Solid] & Crown Z(x) [Dashed] [m]');
    title(ax3, 'Frustum Blending Law Contours (Solid = Span Boat-Tail, Dashed = Crown)', 'FontWeight', 'bold');
    legend(ax3, 'Location', 'northeast', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 4: Full-Stack Longitudinal Area, Half-Span, & Cumulative Volume
    % ---------------------------------------------------------------------
    ax4 = nexttile(tl, 4);
    yyaxis(ax4, 'left');
    hold(ax4, 'on');
    plot(ax4, sec.xStations, sec.AreaX, '-', 'Color', [0.12, 0.46, 0.80], 'LineWidth', 2.4, ...
        'DisplayName', 'Cross-Section Area A(x) [m^2]');
    plot(ax4, sec.xStations, sec.HalfSpanX, '--', 'Color', [0.20, 0.65, 0.35], 'LineWidth', 2.0, ...
        'DisplayName', 'Local Half-Span Y_{max}(x) [m]');
    ylabel(ax4, 'Area [m^2] / Half-Span [m]');

    yyaxis(ax4, 'right');
    plot(ax4, sec.xStations, sec.CumulativeVolumeX, '-', 'Color', [0.86, 0.32, 0.14], 'LineWidth', 2.4, ...
        'DisplayName', 'Cumulative Volume V(x) [m^3]');
    ylabel(ax4, 'Cumulative Volume [m^3]');

    xline(ax4, L, 'k:', 'Waverider Base (X=L)', 'LineWidth', 1.4, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    xline(ax4, L + Lf, 'k-.', 'Booster Ring', 'LineWidth', 1.4, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    grid(ax4, 'on');
    xlabel(ax4, 'Axial Station X [m]');
    title(ax4, sprintf('Full-Stack Area & Volume Distribution (Total Vol = %.2f m^3)', ...
        rep.EnclosedVolume), 'FontWeight', 'bold');
    legend(ax4, 'Location', 'northwest', 'FontSize', 8);

    waverider.viz.stylePublicationFigure(fig, 'finalize');
end
