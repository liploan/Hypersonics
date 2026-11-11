function fig = plotWaverider(body, varargin)
%PLOTWAVERIDER Multi-view 3D visualization of a watertight waverider body.
%   fig = PLOTWAVERIDER(body, Name, Value) renders a 4-panel engineering
%   overview of the waverider geometry (including optional conformal aft
%   transition frustum and conical rocket booster stage):
%     1. 3D Isometric view of the watertight mesh patches
%     2. 3D Compression surface colored by pressure coefficient Cp + shock sheet
%     3. Aft cross-section contour (and circular booster interface ring if enabled)
%     4. Longitudinal symmetry profile view (Y = 0)

    p = inputParser;
    addParameter(p, 'ShowShock', true, @islogical);
    addParameter(p, 'ShowEdges', false, @islogical);
    addParameter(p, 'Visible', 'on', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    opts = p.Results;

    cfg = body.Config;
    surfData = body.SurfaceData;
    meshData = body.MeshData;
    rep = body.Report;
    NxWave = surfData.NxWaverider;

    fig = figure('Color', 'w', 'Position', [80, 80, 1280, 820], 'Visible', opts.Visible);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    if strlength(string(opts.Title)) > 0
        titleStr = char(opts.Title);
    else
        mountTag = '';
        if surfData.HasAftFrustum
            mountTag = ' + Conformal Booster Frustum';
        end
        titleStr = sprintf('Watertight Hypersonic Waverider [%s%s | M = %.1f | \\beta = %.1f^\\circ]', ...
            upper(strrep(cfg.Method, '_', ' ')), mountTag, cfg.Mach, cfg.ShockAngleDeg);
    end
    title(tl, titleStr, 'FontSize', 14, 'FontWeight', 'bold');

    edgeColor = 'none';
    if opts.ShowEdges
        edgeColor = [0.25, 0.25, 0.25];
    end

    % ---------------------------------------------------------------------
    % Tile 1: 3D Watertight Triangulated Mesh by Surface Patch
    % ---------------------------------------------------------------------
    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    V = meshData.Vertices;

    patch(ax1, 'Faces', meshData.FacesUpper, 'Vertices', V, ...
        'FaceColor', [0.20, 0.52, 0.82], 'EdgeColor', edgeColor, 'FaceAlpha', 0.95, ...
        'DisplayName', 'Waverider Upper Surface');
    patch(ax1, 'Faces', meshData.FacesLower, 'Vertices', V, ...
        'FaceColor', [0.88, 0.42, 0.22], 'EdgeColor', edgeColor, 'FaceAlpha', 0.95, ...
        'DisplayName', 'Waverider Lower Surface');
    if ~isempty(meshData.FacesFrustum)
        patch(ax1, 'Faces', meshData.FacesFrustum, 'Vertices', V, ...
            'FaceColor', [0.62, 0.66, 0.72], 'EdgeColor', edgeColor, 'FaceAlpha', 0.95, ...
            'DisplayName', 'Conformal Transition Frustum');
    end
    if ~isempty(meshData.FacesRocket)
        patch(ax1, 'Faces', meshData.FacesRocket, 'Vertices', V, ...
            'FaceColor', [0.85, 0.86, 0.88], 'EdgeColor', edgeColor, 'FaceAlpha', 0.95, ...
            'DisplayName', 'Conical Rocket Booster');
    end
    patch(ax1, 'Faces', meshData.FacesBase, 'Vertices', V, ...
        'FaceColor', [0.35, 0.72, 0.45], 'EdgeColor', edgeColor, 'FaceAlpha', 0.95, ...
        'DisplayName', 'Aft Watertight Base Cap');

    plot3(ax1, surfData.X_low(1, :), surfData.Y_low(1, :), surfData.Z_low(1, :), ...
        'k-', 'LineWidth', 2.0, 'DisplayName', 'Sharp Leading Edge');

    axis(ax1, 'equal');
    grid(ax1, 'on');
    view(ax1, [-138, 22]);
    camlight(ax1, 'headlight');
    lighting(ax1, 'gouraud');
    xlabel(ax1, 'X [m]'); ylabel(ax1, 'Y [m]'); zlabel(ax1, 'Z [m]');
    title(ax1, sprintf('Watertight Mesh (%d Faces, Boundary Edges = %d)', ...
        rep.NumFaces, rep.NumBoundaryEdges), 'FontWeight', 'bold');
    legend(ax1, 'Location', 'southwest', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 2: Lower Surface Cp Distribution & Attached Shock Wave Sheet
    % ---------------------------------------------------------------------
    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    surf(ax2, surfData.X_low, surfData.Y_low, surfData.Z_low, surfData.Cp_low, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.96, 'DisplayName', 'Lower Surface C_p');
    surf(ax2, surfData.X_up, surfData.Y_up, surfData.Z_up, ...
        'FaceColor', [0.75, 0.80, 0.86], 'EdgeColor', 'none', 'FaceAlpha', 0.30, ...
        'DisplayName', 'Upper Envelope');
    if opts.ShowShock
        surf(ax2, surfData.X_shock, surfData.Y_shock, surfData.Z_shock, ...
            'FaceColor', [0.92, 0.18, 0.25], 'EdgeColor', 'none', 'FaceAlpha', 0.22, ...
            'DisplayName', 'Attached Shock Sheet');
    end
    colormap(ax2, 'turbo');
    cb = colorbar(ax2);
    cb.Label.String = 'Compression C_p';
    axis(ax2, 'equal');
    grid(ax2, 'on');
    view(ax2, [-125, -18]);
    xlabel(ax2, 'X [m]'); ylabel(ax2, 'Y [m]'); zlabel(ax2, 'Z [m]');
    title(ax2, 'Shock Capture & Lower Surface Pressure Coefficient', 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Tile 3: Aft Cross-Sections (Waverider X = L & Booster Ring if enabled)
    % ---------------------------------------------------------------------
    ax3 = nexttile(tl, 3);
    hold(ax3, 'on');
    yBaseUp = surfData.Y_up(NxWave, :);
    zBaseUp = surfData.Z_up(NxWave, :);
    yBaseLow = surfData.Y_low(NxWave, :);
    zBaseLow = surfData.Z_low(NxWave, :);
    yShockBase = surfData.Y_shock(end, :);
    zShockBase = surfData.Z_shock(end, :);

    fill(ax3, [yBaseUp, fliplr(yBaseLow)], [zBaseUp, fliplr(zBaseLow)], ...
        [0.82, 0.90, 0.96], 'EdgeColor', 'none', 'FaceAlpha', 0.65, ...
        'DisplayName', 'Waverider Base (X = L)');
    plot(ax3, yBaseUp, zBaseUp, '-', 'Color', [0.15, 0.45, 0.78], 'LineWidth', 2.2, ...
        'DisplayName', 'Upper TE (X = L)');
    plot(ax3, yBaseLow, zBaseLow, '-', 'Color', [0.85, 0.35, 0.15], 'LineWidth', 2.2, ...
        'DisplayName', 'Lower TE (X = L)');
    plot(ax3, yShockBase, zShockBase, '--', 'Color', [0.85, 0.12, 0.20], 'LineWidth', 1.8, ...
        'DisplayName', 'Shock Wave (X = L)');

    if surfData.HasAftFrustum
        NxF = surfData.NxFrustumEnd;
        yRing = [surfData.Y_up(NxF, :), fliplr(surfData.Y_low(NxF, :))];
        zRing = [surfData.Z_up(NxF, :), fliplr(surfData.Z_low(NxF, :))];
        plot(ax3, yRing, zRing, '-', 'Color', [0.25, 0.25, 0.28], 'LineWidth', 2.4, ...
            'DisplayName', sprintf('Booster Interface Ring (R=%.2fm)', cfg.BoosterRadius));
    end

    axis(ax3, 'equal');
    grid(ax3, 'on');
    xlabel(ax3, 'Span Y [m]'); ylabel(ax3, 'Vertical Z [m]');
    title(ax3, sprintf('Cross-Section at X = %.1f m (Vol = %.3f m^3, \\eta_V = %.3f)', ...
        cfg.Length, rep.EnclosedVolume, rep.VolumetricEfficiency), 'FontWeight', 'bold');
    legend(ax3, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 4: Symmetry Plane Profile (Y = 0)
    % ---------------------------------------------------------------------
    ax4 = nexttile(tl, 4);
    hold(ax4, 'on');
    midCol = cfg.NumHalfSpan;
    xSym = surfData.X_up(:, midCol);
    zSymUp = surfData.Z_up(:, midCol);
    zSymLow = surfData.Z_low(:, midCol);

    xSymShock = surfData.X_shock(:, midCol);
    zSymShock = surfData.Z_shock(:, midCol);

    fill(ax4, [xSym; flipud(xSym)], [zSymUp; flipud(zSymLow)], ...
        [0.92, 0.86, 0.78], 'EdgeColor', 'none', 'FaceAlpha', 0.65, ...
        'DisplayName', 'Centerline Slice (Y = 0)');
    plot(ax4, xSym, zSymUp, '-', 'Color', [0.15, 0.45, 0.78], 'LineWidth', 2.2, ...
        'DisplayName', 'Upper Centerline');
    plot(ax4, xSym, zSymLow, '-', 'Color', [0.85, 0.35, 0.15], 'LineWidth', 2.2, ...
        'DisplayName', 'Lower Centerline');
    plot(ax4, xSymShock, zSymShock, '--', 'Color', [0.85, 0.12, 0.20], 'LineWidth', 1.8, ...
        'DisplayName', 'Centerline Shock Ray');

    if surfData.HasAftFrustum
        xline(ax4, cfg.Length, 'k:', 'Waverider Base', 'LineWidth', 1.3, 'HandleVisibility', 'off');
        xline(ax4, cfg.Length + cfg.FrustumLength, 'k-.', 'Booster Mount', 'LineWidth', 1.3, 'HandleVisibility', 'off');
    end

    axis(ax4, 'equal');
    grid(ax4, 'on');
    xlabel(ax4, 'Axial X [m]'); ylabel(ax4, 'Vertical Z [m]');
    title(ax4, sprintf('Centerline Profile (Y = 0) | S_{plan} = %.2f m^2, AR = %.2f', ...
        rep.PlanformArea, rep.AspectRatio), 'FontWeight', 'bold');
    legend(ax4, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);
end
