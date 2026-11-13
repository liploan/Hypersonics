function fig = plotGeometryLines(body, varargin)
%PLOTGEOMETRYLINES Render 3-view lofting blueprint, station slices, and area distribution.
%   fig = PLOTGEOMETRYLINES(body, Name, Value) generates a 4-panel
%   geometric engineering lines plan:
%     1. 3D Parametric Lofting Frame (streamlines, cross-station ribs, LE, CG)
%     2. Multi-Station Transverse Cross-Sections at x/L = [0.25, 0.50, 0.75, 1.00]
%     3. Top Planform & Streamline Blueprint with Chord/Span Contours
%     4. Longitudinal Area-Rule Distribution A(x/L), Cumulative Volume V(x/L), & Thickness

    p = inputParser;
    addParameter(p, 'Visible', 'on', @(x) ischar(x) || isstring(x));
    addParameter(p, 'StationFractions', [0.25, 0.50, 0.75, 1.00], @isnumeric);
    parse(p, varargin{:});
    opts = p.Results;

    cfg = body.Config;
    surfData = body.SurfaceData;
    rep = body.Report;
    sec = body.SectionalAttributes;
    L = cfg.Length;

    fig = figure('Color', 'w', 'Position', [70, 70, 1400, 880], 'Visible', opts.Visible);
    waverider.viz.stylePublicationFigure(fig, 'init');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'loose', 'Padding', 'compact');

    title(tl, sprintf('Geometric Lofting Blueprint & Sectional Lines [%s | L = %.1f m, b = %.2f m]', ...
        upper(strrep(cfg.Method, '_', ' ')), cfg.Length, 2.0 * cfg.HalfSpan), ...
        'FontSize', 14, 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Tile 1: 3D Parametric Lofting Rib & Streamline Wireframe
    % ---------------------------------------------------------------------
    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');

    % Translucent upper and lower surfaces
    surf(ax1, surfData.X_up, surfData.Y_up, surfData.Z_up, ...
        'FaceColor', [0.22, 0.54, 0.84], 'EdgeColor', 'none', 'FaceAlpha', 0.35, ...
        'DisplayName', 'Upper Loft Surface');
    surf(ax1, surfData.X_low, surfData.Y_low, surfData.Z_low, ...
        'FaceColor', [0.88, 0.44, 0.20], 'EdgeColor', 'none', 'FaceAlpha', 0.45, ...
        'DisplayName', 'Lower Compression Surface');

    % Plot sparse longitudinal streamlines
    [Nx, Ny] = size(surfData.X_low);
    colIdx = unique(round(linspace(1, Ny, 15)));
    for k = colIdx
        plot3(ax1, surfData.X_low(:, k), surfData.Y_low(:, k), surfData.Z_low(:, k), ...
            '-', 'Color', [0.65, 0.22, 0.08], 'LineWidth', 0.9, 'HandleVisibility', 'off');
        plot3(ax1, surfData.X_up(:, k), surfData.Y_up(:, k), surfData.Z_up(:, k), ...
            '-', 'Color', [0.10, 0.32, 0.58], 'LineWidth', 0.9, 'HandleVisibility', 'off');
    end

    % Plot transverse rib stations in 3D
    stFracs = opts.StationFractions;
    stColors = lines(numel(stFracs));
    for s = 1:numel(stFracs)
        xs = stFracs(s) * L;
        [yLoop, zLoop] = extractStationLoop(surfData, xs);
        if ~isempty(yLoop)
            xLoop = repmat(xs, size(yLoop));
            plot3(ax1, xLoop, yLoop, zLoop, '-', 'Color', stColors(s, :), ...
                'LineWidth', 2.2, 'DisplayName', sprintf('Station x/L = %.2f', stFracs(s)));
        end
    end

    % Leading edge & volume centroid
    plot3(ax1, surfData.X_low(1, :), surfData.Y_low(1, :), surfData.Z_low(1, :), ...
        'k-', 'LineWidth', 2.4, 'DisplayName', 'Leading Edge Curve');
    cg = rep.VolumeCentroid;
    plot3(ax1, cg(1), cg(2), cg(3), 'p', 'MarkerSize', 12, ...
        'MarkerFaceColor', [0.95, 0.80, 0.10], 'MarkerEdgeColor', 'k', ...
        'DisplayName', sprintf('Volume Centroid (X=%.2fm)', cg(1)));

    axis(ax1, 'equal');
    grid(ax1, 'on');
    view(ax1, [-132, 24]);
    xlabel(ax1, 'Axial X [m]'); ylabel(ax1, 'Span Y [m]'); zlabel(ax1, 'Vertical Z [m]');
    title(ax1, '3D Parametric Rib Lofting & Streamline Network', 'FontWeight', 'bold');
    legend(ax1, 'Location', 'southwest', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 2: Multi-Station Transverse Cross-Section Evolution ("Body Plan")
    % ---------------------------------------------------------------------
    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    for s = 1:numel(stFracs)
        xs = stFracs(s) * L;
        [yLoop, zLoop, yShk, zShk] = extractStationLoop(surfData, xs);
        if ~isempty(yLoop)
            plot(ax2, yLoop, zLoop, '-', 'Color', stColors(s, :), 'LineWidth', 2.2, ...
                'DisplayName', sprintf('Body x/L = %.2f (X = %.2f m)', stFracs(s), xs));
            plot(ax2, yShk, zShk, '--', 'Color', stColors(s, :), 'LineWidth', 1.2, ...
                'HandleVisibility', 'off');
        end
    end
    axis(ax2, 'equal');
    grid(ax2, 'on');
    xlabel(ax2, 'Spanwise Coordinate Y [m]');
    ylabel(ax2, 'Vertical Coordinate Z [m]');
    title(ax2, 'Transverse Station Cross-Sections (Solid = Body, Dashed = Shock)', 'FontWeight', 'bold');
    legend(ax2, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 3: Top Planform View with Local Thickness Contours & Streamlines
    % ---------------------------------------------------------------------
    ax3 = nexttile(tl, 3);
    hold(ax3, 'on');
    dZ_thickness = max(surfData.Z_up - surfData.Z_low, 0.0);
    contourf(ax3, surfData.X_up, surfData.Y_up, dZ_thickness, 14, ...
        'LineColor', 'none', 'HandleVisibility', 'off');
    colormap(ax3, 'parula');
    cb3 = colorbar(ax3);
    cb3.Label.String = 'Local Thickness \DeltaZ(x,y) [m]';

    % Overlay leading edge, trailing edge, and streamlines
    plot(ax3, surfData.X_low(1, :), surfData.Y_low(1, :), 'k-', 'LineWidth', 2.4, ...
        'DisplayName', 'Leading Edge');
    plot(ax3, surfData.X_low(end, :), surfData.Y_low(end, :), 'k--', 'LineWidth', 1.8, ...
        'DisplayName', 'Aft Base Plane (X = L)');
    for k = colIdx
        plot(ax3, surfData.X_low(:, k), surfData.Y_low(:, k), '-', ...
            'Color', [1.0, 1.0, 1.0, 0.65], 'LineWidth', 0.9, 'HandleVisibility', 'off');
    end
    plot(ax3, cg(1), cg(2), 'kp', 'MarkerSize', 11, 'MarkerFaceColor', 'y', ...
        'DisplayName', 'Volume Centroid');

    axis(ax3, 'equal');
    grid(ax3, 'on');
    xlabel(ax3, 'Axial X [m]');
    ylabel(ax3, 'Span Y [m]');
    title(ax3, sprintf('Planform Thickness Map & Streamlines (S_{plan} = %.2f m^2, P = %.2f)', ...
        rep.PlanformArea, cfg.PlanformPower), 'FontWeight', 'bold');
    legend(ax3, 'Location', 'northwest', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 4: Longitudinal Area-Rule A(x/L), Cumulative Volume V(x/L), & Thickness
    % ---------------------------------------------------------------------
    ax4 = nexttile(tl, 4);
    yyaxis(ax4, 'left');
    hold(ax4, 'on');
    plot(ax4, sec.xNorm, sec.AreaX, '-', 'Color', [0.12, 0.45, 0.78], 'LineWidth', 2.4, ...
        'DisplayName', 'Cross-Section Area A(x) [m^2]');
    plot(ax4, sec.xNorm, sec.MaxThicknessX, '--', 'Color', [0.18, 0.65, 0.35], 'LineWidth', 2.0, ...
        'DisplayName', 'Max Thickness t_{max}(x) [m]');
    ylabel(ax4, 'Area [m^2] / Max Thickness [m]');

    yyaxis(ax4, 'right');
    plot(ax4, sec.xNorm, sec.CumulativeVolumeX, '-', 'Color', [0.85, 0.32, 0.15], 'LineWidth', 2.4, ...
        'DisplayName', 'Cumulative Volume V(x) [m^3]');
    ylabel(ax4, 'Cumulative Enclosed Volume [m^3]');

    grid(ax4, 'on');
    xlabel(ax4, 'Normalized Axial Station x / L [-]');
    title(ax4, sprintf('Longitudinal Area Progression (Total Vol = %.3f m^3, Base Area = %.3f m^2)', ...
        rep.EnclosedVolume, rep.BaseArea), 'FontWeight', 'bold');
    legend(ax4, 'Location', 'northwest', 'FontSize', 8);

    waverider.viz.stylePublicationFigure(fig, 'finalize');
end

function [yLoop, zLoop, yShk, zShk] = extractStationLoop(surfData, xs)
    if isfield(surfData, 'IsAxisymmetric') && surfData.IsAxisymmetric
        midCol = ceil(size(surfData.X_up, 2) / 2);
        xAxis = surfData.X_up(:, midCol);
        if xs < xAxis(1) - 1e-7 || xs > xAxis(end) + 1e-7
            yLoop = []; zLoop = []; yShk = []; zShk = [];
            return;
        end
        yUpSlice  = interp1(xAxis, surfData.Y_up, xs, 'linear', 'extrap');
        zUpSlice  = interp1(xAxis, surfData.Z_up, xs, 'linear', 'extrap');
        yLowSlice = interp1(xAxis, surfData.Y_low, xs, 'linear', 'extrap');
        zLowSlice = interp1(xAxis, surfData.Z_low, xs, 'linear', 'extrap');
        xShkAxis  = surfData.X_shock(:, midCol);
        yShkSlice = interp1(xShkAxis, surfData.Y_shock, xs, 'linear', 'extrap');
        zShkSlice = interp1(xShkAxis, surfData.Z_shock, xs, 'linear', 'extrap');

        yLoop = [yUpSlice, fliplr(yLowSlice), yUpSlice(1)];
        zLoop = [zUpSlice, fliplr(zLowSlice), zUpSlice(1)];
        yShk = [yShkSlice, fliplr(yShkSlice), yShkSlice(1)];
        zShk = [-zShkSlice, fliplr(zShkSlice), -zShkSlice(1)];
        return;
    end

    xLE = surfData.X_up(1, :);
    yLE = surfData.Y_up(1, :);
    Ny = numel(yLE);

    ySlice = [];
    zUpSlice = [];
    zLowSlice = [];
    zShkSlice = [];

    for k = 1:Ny
        if xLE(k) <= xs + 1e-7
            xCol = surfData.X_up(:, k);
            if abs(xCol(end) - xCol(1)) < 1e-8
                zu = surfData.Z_up(end, k);
                zl = surfData.Z_low(end, k);
                zs = surfData.Z_shock(end, k);
            else
                zu = interp1(xCol, surfData.Z_up(:, k), xs, 'linear', 'extrap');
                zl = interp1(xCol, surfData.Z_low(:, k), xs, 'linear', 'extrap');
                zs = interp1(xCol, surfData.Z_shock(:, k), xs, 'linear', 'extrap');
            end
            ySlice(end + 1) = yLE(k); %#ok<AGROW>
            zUpSlice(end + 1) = zu; %#ok<AGROW>
            zLowSlice(end + 1) = zl; %#ok<AGROW>
            zShkSlice(end + 1) = zs; %#ok<AGROW>
        end
    end

    if numel(ySlice) < 3
        yLoop = []; zLoop = []; yShk = []; zShk = [];
        return;
    end

    yLoop = [ySlice, fliplr(ySlice), ySlice(1)];
    zLoop = [zUpSlice, fliplr(zLowSlice), zUpSlice(1)];
    yShk = ySlice;
    zShk = zShkSlice;
end
