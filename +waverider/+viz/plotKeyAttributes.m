function fig = plotKeyAttributes(body, varargin)
%PLOTKEYATTRIBUTES Render key aerodynamic, flowfield, spanwise, and mesh attributes.
%   fig = PLOTKEYATTRIBUTES(body, Name, Value) generates a 4-panel dashboard:
%     1. Shock-Layer Flowfield Profile (Taylor-Maccoll Cp(theta) & Mach(theta))
%     2. Angle-of-Attack Aerodynamic Performance Polar (CL, CD_wave, and L/D vs alpha)
%     3. Spanwise Geometric & Aerodynamic Loading Distribution across y/(b/2)
%     4. 3D Surface Pressure Coefficient & Outward Normal Vectors on Watertight Mesh

    p = inputParser;
    addParameter(p, 'Visible', 'on', @(x) ischar(x) || isstring(x));
    addParameter(p, 'AlphaRangeDeg', linspace(-4.0, 16.0, 41), @isnumeric);
    parse(p, varargin{:});
    opts = p.Results;

    cfg = body.Config;
    sec = body.SectionalAttributes;
    polar = body.computeAeroPolar(opts.AlphaRangeDeg);

    % Always compute exact Taylor-Maccoll profile at design (Mach, ShockAngleDeg)
    tm = waverider.flow.taylorMaccoll(cfg.Mach, cfg.ShockAngleDeg, cfg.Gamma);

    fig = figure('Color', 'w', 'Position', [80, 80, 1340, 840], 'Visible', opts.Visible);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    title(tl, sprintf('Waverider Key Aerodynamic, Flowfield & Geometric Attributes [%s | M = %.1f, \\beta = %.1f^\\circ]', ...
        upper(strrep(cfg.Method, '_', ' ')), cfg.Mach, cfg.ShockAngleDeg), ...
        'FontSize', 14, 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Tile 1: Taylor-Maccoll Conical Shock-Layer Profile (theta_c to beta)
    % ---------------------------------------------------------------------
    ax1 = nexttile(tl, 1);
    yyaxis(ax1, 'left');
    hold(ax1, 'on');
    plot(ax1, tm.thetaDeg, tm.Cp, '-', 'Color', [0.85, 0.25, 0.15], 'LineWidth', 2.4, ...
        'DisplayName', 'Pressure Coeff C_p(\theta)');
    plot(ax1, tm.thetaDeg, abs(tm.utheta), '--', 'Color', [0.65, 0.18, 0.55], 'LineWidth', 1.8, ...
        'DisplayName', 'Normal Velocity |u_\theta(\theta)| / V_{max}');
    ylabel(ax1, 'C_p(\theta) / Normalized |u_\theta|');

    yyaxis(ax1, 'right');
    plot(ax1, tm.thetaDeg, tm.Mach, '-', 'Color', [0.12, 0.48, 0.82], 'LineWidth', 2.4, ...
        'DisplayName', 'Local Mach M(\theta)');
    ylabel(ax1, 'Shock-Layer Local Mach Number [-]');

    xline(ax1, tm.thetaConeDeg, 'k:', sprintf('Cone \\theta_c = %.2f^\\circ', tm.thetaConeDeg), ...
        'LineWidth', 1.4, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    xline(ax1, tm.betaDeg, 'r:', sprintf('Shock \\beta = %.2f^\\circ', tm.betaDeg), ...
        'LineWidth', 1.4, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    grid(ax1, 'on');
    xlabel(ax1, 'Conical Polar Ray Angle \theta [deg]');
    title(ax1, sprintf('Taylor-Maccoll Shock Layer (M_\\infty = %.1f, \\theta_c = %.2f^\\circ \\rightarrow \\beta = %.2f^\\circ)', ...
        cfg.Mach, tm.thetaConeDeg, tm.betaDeg), 'FontWeight', 'bold');
    legend(ax1, 'Location', 'north', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 2: Angle-of-Attack Aerodynamic Polar (CL, CD, and L/D vs Alpha)
    % ---------------------------------------------------------------------
    ax2 = nexttile(tl, 2);
    yyaxis(ax2, 'left');
    hold(ax2, 'on');
    plot(ax2, polar.AlphaDeg, polar.CL, '-', 'Color', [0.12, 0.48, 0.82], 'LineWidth', 2.4, ...
        'DisplayName', 'Lift Coefficient C_L(\alpha)');
    plot(ax2, polar.AlphaDeg, polar.CD_wave, '--', 'Color', [0.20, 0.68, 0.35], 'LineWidth', 2.2, ...
        'DisplayName', 'Wave Drag C_{D,wave}(\alpha)');
    ylabel(ax2, 'Force Coefficients C_L, C_{D,wave} [-]');

    yyaxis(ax2, 'right');
    plot(ax2, polar.AlphaDeg, polar.LD_inviscid, '-', 'Color', [0.88, 0.35, 0.12], 'LineWidth', 2.5, ...
        'DisplayName', 'Inviscid (L/D)(\alpha)');
    [maxLD, idxMax] = max(polar.LD_inviscid);
    alphaOpt = polar.AlphaDeg(idxMax);
    plot(ax2, alphaOpt, maxLD, 'kp', 'MarkerSize', 11, 'MarkerFaceColor', [0.98, 0.82, 0.12], ...
        'DisplayName', sprintf('Peak L/D = %.2f at \\alpha = %.1f^\\circ', maxLD, alphaOpt));
    ylabel(ax2, 'Inviscid Lift-to-Drag Ratio L/D [-]');

    xline(ax2, 0.0, 'k:', 'Design \alpha = 0^\circ', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    grid(ax2, 'on');
    xlabel(ax2, 'Angle of Attack \alpha [deg]');
    title(ax2, sprintf('Hypersonic Inviscid Aerodynamic Polar (L/D|_{\\alpha=0} = %.2f)', ...
        body.AeroMetrics.LD_inviscid), 'FontWeight', 'bold');
    legend(ax2, 'Location', 'northwest', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 3: Spanwise Geometric & Aerodynamic Loading Distribution
    % ---------------------------------------------------------------------
    ax3 = nexttile(tl, 3);
    yyaxis(ax3, 'left');
    hold(ax3, 'on');
    plot(ax3, sec.SpanEta, sec.ChordY ./ cfg.Length, '-', 'Color', [0.15, 0.45, 0.78], 'LineWidth', 2.2, ...
        'DisplayName', 'Normalized Chord c(y) / L');
    plot(ax3, sec.SpanEta, (sec.BaseHeightY ./ cfg.Length) * 5.0, '--', 'Color', [0.22, 0.65, 0.42], 'LineWidth', 2.0, ...
        'DisplayName', 'Base Height 5 \times h_{base}(y) / L');
    ylabel(ax3, 'Normalized Geometric Distribution [-]');

    yyaxis(ax3, 'right');
    plot(ax3, sec.SpanEta, sec.SpanLoadNorm, '-', 'Color', [0.85, 0.30, 0.15], 'LineWidth', 2.4, ...
        'DisplayName', 'Sectional Lift Index C_p(y) \cdot c(y)/L');
    ylabel(ax3, 'Sectional Compression Load Index [-]');

    grid(ax3, 'on');
    xlabel(ax3, 'Normalized Spanwise Coordinate \eta = y / (b/2) [-]');
    title(ax3, sprintf('Spanwise Chord, Base Thickness & Lift Distribution (AR = %.2f)', ...
        body.Report.AspectRatio), 'FontWeight', 'bold');
    legend(ax3, 'Location', 'south', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 4: 3D Mesh Face Pressure Cp & Outward Normal Verification Quiver
    % ---------------------------------------------------------------------
    ax4 = nexttile(tl, 4);
    hold(ax4, 'on');
    TR = body.MeshData.TR;
    V = TR.Points;
    F = TR.ConnectivityList;
    CpFace = body.AeroMetrics.CpFace;

    patch(ax4, 'Faces', F, 'Vertices', V, 'FaceVertexCData', CpFace, ...
        'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 0.95, ...
        'DisplayName', 'Surface Face C_p');
    colormap(ax4, 'turbo');
    cb4 = colorbar(ax4);
    cb4.Label.String = 'Surface Pressure Coefficient C_p (\alpha = 0^\circ)';

    % Subsample face centers and outward unit normals for quiver display
    fc = incenter(TR);
    fn = faceNormal(TR);
    step = max(floor(size(F, 1) / 220), 1);
    idxQ = 1:step:size(F, 1);
    quiver3(ax4, fc(idxQ, 1), fc(idxQ, 2), fc(idxQ, 3), ...
        fn(idxQ, 1), fn(idxQ, 2), fn(idxQ, 3), 0.45, ...
        'Color', [0.15, 0.15, 0.15], 'LineWidth', 0.8, ...
        'DisplayName', 'Outward Unit Normals');

    axis(ax4, 'equal');
    grid(ax4, 'on');
    view(ax4, [-135, 22]);
    xlabel(ax4, 'X [m]'); ylabel(ax4, 'Y [m]'); zlabel(ax4, 'Z [m]');
    title(ax4, sprintf('3D Face C_p & Verified Outward Manifold Normals (\\chi = %d)', ...
        body.Report.EulerCharacteristic), 'FontWeight', 'bold');
    legend(ax4, 'Location', 'southwest', 'FontSize', 8);
end
