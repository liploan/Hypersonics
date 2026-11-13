function [fig, summaryMetrics] = plotRCSAnalysis(body, varargin)
%PLOTRCSANALYSIS Multi-band (UHF to X-band) RCS & concave underside scattering analysis.
%   [fig, summaryMetrics] = PLOTRCSANALYSIS(body, Name, Value) generates a
%   4-panel electromagnetic scattering dashboard from UHF (0.3 GHz) to
%   X-band (10+ GHz):
%     1. Multi-Band Azimuthal RCS Pattern (UHF, L, S, C, X bands) at Ventral El = -20 deg
%     2. Ventral vs. Dorsal Elevation Cut (el = -90 to +90 deg) showing Concave Double-Bounce
%     3. Continuous Frequency Response (0.3 to 12 GHz) across Key Threat Aspects
%     4. 3D Underside Concave Cavity Hotspot & Port-to-Starboard Double-Bounce Ray Paths

    p = inputParser;
    addParameter(p, 'Visible', 'on', @(x) ischar(x) || isstring(x));
    addParameter(p, 'VentralElevationDeg', -20.0, @isnumeric);
    addParameter(p, 'Polarization', 'Mean', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    opts = p.Results;

    cfg = body.Config;
    cp = waverider.rcs.buildConcavePairList(body);

    % Radar bands of interest: UHF (0.5 GHz), L (1.5 GHz), S (3.0 GHz), C (6.0 GHz), X (10.0 GHz)
    bandFreqsGHz = [0.5, 1.5, 3.0, 6.0, 10.0];
    bandNames = {'UHF (0.5 GHz, \lambda=60cm)', ...
                 'L-Band (1.5 GHz, \lambda=20cm)', ...
                 'S-Band (3.0 GHz, \lambda=10cm)', ...
                 'C-Band (6.0 GHz, \lambda=5cm)', ...
                 'X-Band (10.0 GHz, \lambda=3cm)'};
    bandColors = [0.15, 0.48, 0.82;
                  0.18, 0.68, 0.38;
                  0.92, 0.62, 0.12;
                  0.85, 0.28, 0.15;
                  0.58, 0.20, 0.72];

    % 1. Compute Azimuthal Sweep (0 to 180 deg) at Ventral Elevation
    azVec = linspace(0.0, 180.0, 181).';
    rcsAz = waverider.rcs.computeRCS(body, bandFreqsGHz, azVec, opts.VentralElevationDeg, ...
        'Polarization', opts.Polarization, 'PrecomputedConcavePairs', cp);

    % 2. Compute Elevation Cut (-90 deg Ventral to +90 deg Dorsal) at Oblique Az = 45 deg
    elVec = linspace(-90.0, 90.0, 181).';
    rcsEl = waverider.rcs.computeRCS(body, [0.5, 3.0, 10.0], 45.0, elVec, ...
        'Polarization', opts.Polarization, 'PrecomputedConcavePairs', cp);

    % 3. Compute Continuous Frequency Sweep (0.3 GHz UHF to 12.0 GHz X-band)
    freqSweepGHz = linspace(0.3, 12.0, 60);
    aspectAz = [0.0;   45.0;  45.0;  85.0];
    aspectEl = [-10.0; -22.0; +22.0; -20.0];
    rcsFreq = waverider.rcs.computeRCS(body, freqSweepGHz, aspectAz, aspectEl, ...
        'Polarization', opts.Polarization, 'PrecomputedConcavePairs', cp);

    fig = figure('Color', 'w', 'Position', [55, 55, 1420, 890], 'Visible', opts.Visible);
    waverider.viz.stylePublicationFigure(fig, 'init');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'loose', 'Padding', 'compact');

    title(tl, sprintf('Multi-Band RCS & Concave Underside Scattering Analysis (UHF to X-Band) [%s | L = %.1f m | Pol: %s]', ...
        upper(strrep(cfg.Method, '_', ' ')), cfg.Length, rcsAz.Polarization), ...
        'FontSize', 14, 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Tile 1: Multi-Band Azimuthal RCS Sweep (UHF, L, S, C, X)
    % ---------------------------------------------------------------------
    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    for b = 1:numel(bandFreqsGHz)
        % Smooth over 3 deg running window to show envelope clearly alongside coherentlobes
        sigSmooth = movmean(rcsAz.SigmaTotal_dBsm(:, b), 3);
        plot(ax1, azVec, sigSmooth, '-', 'Color', bandColors(b, :), 'LineWidth', 2.0, ...
            'DisplayName', bandNames{b});
    end
    xline(ax1, 68.0, 'k:', 'LE Normal Specular', 'LineWidth', 1.3, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    grid(ax1, 'on');
    xlim(ax1, [0, 180]);
    xlabel(ax1, 'Azimuth Angle \psi [deg] (0^\circ = Nose-On, 90^\circ = Broadside, 180^\circ = Tail)');
    ylabel(ax1, 'Monostatic RCS \sigma [dBsm]');
    title(ax1, sprintf('Ventral Azimuthal RCS Pattern (Elevation \\theta_{el} = %.0f^\\circ Look-Up)', ...
        opts.VentralElevationDeg), 'FontWeight', 'bold');
    legend(ax1, 'Location', 'northwest', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 2: Ventral vs. Dorsal Elevation Cut (Single-Bounce vs Concave Multi-Bounce)
    % ---------------------------------------------------------------------
    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    elColors = [0.15, 0.48, 0.82; 0.90, 0.55, 0.10; 0.58, 0.20, 0.72];
    elFreqNames = {'UHF (0.5 GHz)', 'S-Band (3.0 GHz)', 'X-Band (10.0 GHz)'};
    for k = 1:3
        totCurve = movmean(rcsEl.SigmaTotal_dBsm(:, k), 3);
        sngCurve = movmean(rcsEl.SigmaSingle_dBsm(:, k), 3);
        plot(ax2, elVec, totCurve, '-', 'Color', elColors(k, :), 'LineWidth', 2.2, ...
            'DisplayName', sprintf('%s Total (PO + Multi-Bounce)', elFreqNames{k}));
        plot(ax2, elVec, sngCurve, '--', 'Color', elColors(k, :), 'LineWidth', 1.4, ...
            'DisplayName', sprintf('%s Single-Bounce Only', elFreqNames{k}));
    end
    xline(ax2, 0.0, 'k-', 'Horizon (\theta_{el} = 0^\circ)', 'LineWidth', 1.4, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    grid(ax2, 'on');
    xlim(ax2, [-90, 90]);
    xlabel(ax2, 'Elevation Angle \theta_{el} [deg] (< 0^\circ Ventral Concave Underside | > 0^\circ Dorsal Canopy)');
    ylabel(ax2, 'Monostatic RCS \sigma [dBsm]');
    title(ax2, 'Ventral vs. Dorsal Elevation Cut at Oblique Azimuth \psi = 45^\circ', 'FontWeight', 'bold');
    legend(ax2, 'Location', 'northeast', 'FontSize', 7.5, 'NumColumns', 2);

    % ---------------------------------------------------------------------
    % Tile 3: Continuous Frequency Sweep (0.3 GHz UHF to 12.0 GHz X-Band)
    % ---------------------------------------------------------------------
    ax3 = nexttile(tl, 3);
    hold(ax3, 'on');
    aspectLabels = {'Nose-On Ventral (\psi=0^\circ, \theta_{el}=-10^\circ)', ...
                    'Oblique Ventral Concave Trough (\psi=45^\circ, \theta_{el}=-22^\circ)', ...
                    'Oblique Dorsal Convex Canopy (\psi=45^\circ, \theta_{el}=+22^\circ)', ...
                    'Broadside Ventral (\psi=85^\circ, \theta_{el}=-20^\circ)'};
    aspectStyles = {'-o', '-s', '--^', '-d'};
    aspectCols = [0.15, 0.48, 0.82; 0.85, 0.25, 0.15; 0.20, 0.68, 0.35; 0.58, 0.20, 0.72];

    for a = 1:4
        yCurve = movmean(rcsFreq.SigmaTotal_dBsm(a, :), 3);
        plot(ax3, freqSweepGHz, yCurve, aspectStyles{a}, 'Color', aspectCols(a, :), ...
            'LineWidth', 2.0, 'MarkerSize', 4, 'MarkerIndices', 1:6:numel(freqSweepGHz), ...
            'DisplayName', aspectLabels{a});
    end
    % Mark IEEE radar band boundaries
    xline(ax3, 1.0, 'k:', 'UHF|L', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    xline(ax3, 2.0, 'k:', 'L|S',   'LineWidth', 1.0, 'HandleVisibility', 'off');
    xline(ax3, 4.0, 'k:', 'S|C',   'LineWidth', 1.0, 'HandleVisibility', 'off');
    xline(ax3, 8.0, 'k:', 'C|X',   'LineWidth', 1.0, 'HandleVisibility', 'off');

    grid(ax3, 'on');
    xlim(ax3, [0.3, 12.0]);
    xlabel(ax3, 'Radar Frequency f [GHz] (UHF 0.3-1 | L 1-2 | S 2-4 | C 4-8 | X 8-12 GHz)');
    ylabel(ax3, 'Monostatic RCS \sigma(f) [dBsm]');
    title(ax3, 'Frequency Scaling Across UHF, L, S, C, and X Bands', 'FontWeight', 'bold');
    legend(ax3, 'Location', 'southwest', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 4: 3D Underside Concave Cavity Hotspot & Double-Bounce Ray Paths
    % ---------------------------------------------------------------------
    ax4 = nexttile(tl, 4);
    hold(ax4, 'on');
    surfData = body.SurfaceData;
    NxW = surfData.NxWaverider;
    X_low = surfData.X_low(1:NxW, :);
    Y_low = surfData.Y_low(1:NxW, :);
    Z_low = surfData.Z_low(1:NxW, :);

    % Color lower surface by transverse concave depth relative to straight chord between wingtips
    zLine = repmat(Z_low(:, 1), 1, size(Z_low, 2));
    concavityDepth = max(Z_low - zLine, 0.0);
    surf(ax4, X_low, Y_low, Z_low, concavityDepth, 'EdgeColor', 'none', 'FaceAlpha', 0.92, ...
        'DisplayName', 'Concave Trough Depth [m]');
    colormap(ax4, 'turbo');
    cb4 = colorbar(ax4);
    cb4.Label.String = 'Underside Concave Cavity Depth \DeltaZ_{arch} [m]';

    % Plot top mutually visible port-to-starboard double-bounce coupling rays
    if cp.numPairs > 0
        nShow = min(40, cp.numPairs);
        for m = 1:nShow
            pts = [cp.rP(m, :); cp.rQ(m, :)];
            if m == 1
                plot3(ax4, pts(:, 1), pts(:, 2), pts(:, 3), 'w-', 'LineWidth', 1.4, ...
                    'DisplayName', 'Concave Double-Bounce Coupling Rays');
            else
                plot3(ax4, pts(:, 1), pts(:, 2), pts(:, 3), 'w-', 'LineWidth', 1.1, ...
                    'HandleVisibility', 'off');
            end
        end
        plot3(ax4, cp.rP(1:nShow, 1), cp.rP(1:nShow, 2), cp.rP(1:nShow, 3), 'yo', ...
            'MarkerSize', 4, 'MarkerFaceColor', 'y', 'DisplayName', 'Cavity Bounce Centres');
    end

    axis(ax4, 'equal');
    grid(ax4, 'on');
    view(ax4, [-128, -26]); % Bottom-up ventral view into concave underside
    xlabel(ax4, 'X [m]'); ylabel(ax4, 'Y [m]'); zlabel(ax4, 'Z [m]');
    title(ax4, sprintf('3D Ventral View: Concave Cavity & %d Double-Bounce Pairs', cp.numPairs), ...
        'FontWeight', 'bold');
    legend(ax4, 'Location', 'southwest', 'FontSize', 8);

    waverider.viz.stylePublicationFigure(fig, 'finalize');

    summaryMetrics = struct();
    summaryMetrics.BandFreqsGHz = bandFreqsGHz;
    summaryMetrics.AzSweep = rcsAz;
    summaryMetrics.ElSweep = rcsEl;
    summaryMetrics.FreqSweep = rcsFreq;
    summaryMetrics.NumConcavePairs = cp.numPairs;
end
