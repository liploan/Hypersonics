function fig = plotRangeAngleIntensity(body, varargin)
%PLOTRANGEANGLEINTENSITY Generate 4-panel Range-Angle-Intensity (RTI) & Bandwidth Resolution analysis.
%   fig = PLOTRANGEANGLEINTENSITY(body, Name, Value) evaluates and plots:
%     Panel 1: 2D Range-Angle-Intensity Sinogram (Total: Single + Concave Double-Bounce)
%     Panel 2: Concave Cavity Double-Bounce Range-Angle Enhancement (Delayed Range Lobes)
%     Panel 3: Waveform Bandwidth Comparison (B = 150, 600, 2000 MHz -> Delta R = 100, 25, 7.5 cm)
%     Panel 4: 1D High-Resolution Range Profile (HRRP) Decomposition at Oblique Ventral Look

    p = inputParser;
    addParameter(p, 'CenterFreqGHz', 10.0, @isnumeric);
    addParameter(p, 'BandwidthMHz', 1000.0, @isnumeric);
    addParameter(p, 'ElevationDeg', -22.0, @isnumeric);
    addParameter(p, 'CutAzimuthDeg', 45.0, @isnumeric);
    addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    opts = p.Results;

    fc = opts.CenterFreqGHz;
    bw = opts.BandwidthMHz;
    elLook = opts.ElevationDeg;
    azCut = opts.CutAzimuthDeg;

    azVec = linspace(0, 180, 121);
    pairs = waverider.rcs.buildConcavePairList(body);

    % Compute primary RTI at configured bandwidth
    rtiMain = waverider.rcs.computeRangeProfile(body, fc, bw, azVec, elLook, ...
        'IncludeMultiBounce', true, 'WindowType', 'hann', 'PrecomputedConcavePairs', pairs);

    % Compute multi-bandwidth HRRP cuts at azCut
    bwList = [150, 600, 2000];
    rtiLowBW  = waverider.rcs.computeRangeProfile(body, fc, bwList(1), azCut, elLook, ...
        'IncludeMultiBounce', true, 'WindowType', 'hann', 'PrecomputedConcavePairs', pairs);
    rtiMedBW  = waverider.rcs.computeRangeProfile(body, fc, bwList(2), azCut, elLook, ...
        'IncludeMultiBounce', true, 'WindowType', 'hann', 'PrecomputedConcavePairs', pairs);
    rtiHighBW = waverider.rcs.computeRangeProfile(body, fc, bwList(3), azCut, elLook, ...
        'IncludeMultiBounce', true, 'WindowType', 'hann', 'PrecomputedConcavePairs', pairs);

    fig = figure('Name', 'Waverider Range-Angle-Intensity & Waveform Bandwidth Analysis', ...
        'Color', [0.05, 0.07, 0.12], 'Position', [60, 60, 1520, 920]);

    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Range-Angle-Intensity (Sinogram) & Waveform Bandwidth Analysis  (f_c = %.1f GHz, B = %.0f MHz, \\DeltaR = %.1f cm, \\theta_{el} = %.0f^\\circ)', ...
        fc, bw, rtiMain.deltaR_m * 100, elLook), ...
        'Color', [0.95, 0.97, 1.0], 'FontSize', 14, 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Panel 1: Total Range-Angle-Intensity Sinogram (Single + Double Bounce)
    % ---------------------------------------------------------------------
    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    imagesc(ax1, rtiMain.rangeVec_m, rtiMain.azVecDeg, rtiMain.RTI_total_dBsm);
    set(ax1, 'YDir', 'normal');
    colormap(ax1, turbo);
    clim(ax1, [-38, 18]);
    cb1 = colorbar(ax1, 'Color', [0.85, 0.90, 0.96]);
    cb1.Label.String = 'Scattering Intensity \sigma(R, \psi) [dBsm]';
    cb1.Label.FontSize = 10;

    % Overlay geometric Nose and Base trajectories
    plot(ax1, rtiMain.R_nose_m, rtiMain.azVecDeg, '--', 'Color', [0.22, 0.95, 0.98], ...
        'LineWidth', 1.6, 'DisplayName', 'Nose Tip Locus R_{nose}(\psi)');
    plot(ax1, rtiMain.R_waveriderBase_m, rtiMain.azVecDeg, '-.', 'Color', [0.98, 0.75, 0.20], ...
        'LineWidth', 1.5, 'DisplayName', 'Waverider Base (X = L)');
    if max(abs(rtiMain.R_aftBase_m - rtiMain.R_waveriderBase_m)) > 0.05
        plot(ax1, rtiMain.R_aftBase_m, rtiMain.azVecDeg, ':', 'Color', [0.35, 0.95, 0.55], ...
            'LineWidth', 1.6, 'DisplayName', 'Booster Aft Base');
    end
    yline(ax1, azCut, '--w', sprintf('HRRP Cut (\\psi = %.0f^\\circ)', azCut), ...
        'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left', ...
        'Color', [0.95, 0.97, 1.0], 'HandleVisibility', 'off');

    xlim(ax1, [min(rtiMain.rangeVec_m), max(rtiMain.rangeVec_m)]);
    ylim(ax1, [0, 180]);
    xlabel(ax1, 'Down-Range Coordinate R [m]  (Relative to Vehicle Center)', 'Color', [0.92, 0.95, 1.0], 'FontWeight', 'bold');
    ylabel(ax1, 'Azimuth Aspect Angle \psi [deg]  (0^\circ=Nose, 90^\circ=Broadside, 180^\circ=Tail)', 'Color', [0.92, 0.95, 1.0], 'FontWeight', 'bold');
    title(ax1, sprintf('1. Total Range-Angle-Intensity Sinogram (B = %.0f MHz, \\DeltaR = %.1f cm)', ...
        bw, rtiMain.deltaR_m * 100), 'Color', [0.95, 0.97, 1.0], 'FontSize', 11.5);
    legend(ax1, 'Location', 'southwest', 'TextColor', [0.92, 0.95, 1.0], 'Color', [0.08, 0.11, 0.18]);
    styleDarkAxes(ax1);

    % ---------------------------------------------------------------------
    % Panel 2: Concave Cavity Double-Bounce Range-Angle Isolation
    % ---------------------------------------------------------------------
    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    imagesc(ax2, rtiMain.rangeVec_m, rtiMain.azVecDeg, rtiMain.RTI_double_dBsm);
    set(ax2, 'YDir', 'normal');
    colormap(ax2, turbo);
    clim(ax2, [-38, 18]);
    cb2 = colorbar(ax2, 'Color', [0.85, 0.90, 0.96]);
    cb2.Label.String = 'Concave 2-Bounce Intensity [dBsm]';
    cb2.Label.FontSize = 10;

    plot(ax2, rtiMain.R_nose_m, rtiMain.azVecDeg, '--', 'Color', [0.22, 0.95, 0.98], ...
        'LineWidth', 1.4, 'DisplayName', 'Nose Tip Locus');
    plot(ax2, rtiMain.R_waveriderBase_m, rtiMain.azVecDeg, '-.', 'Color', [0.98, 0.75, 0.20], ...
        'LineWidth', 1.4, 'DisplayName', 'Waverider Base (X = L)');

    xlim(ax2, [min(rtiMain.rangeVec_m), max(rtiMain.rangeVec_m)]);
    ylim(ax2, [0, 180]);
    xlabel(ax2, 'Down-Range Coordinate R [m]  (+0.5 d_{PQ} Cavity Delay Visible)', 'Color', [0.92, 0.95, 1.0], 'FontWeight', 'bold');
    ylabel(ax2, 'Azimuth Aspect Angle \psi [deg]', 'Color', [0.92, 0.95, 1.0], 'FontWeight', 'bold');
    title(ax2, '2. Isolated Concave Underside Double-Bounce Sinogram (Delayed Cavity Returns)', ...
        'Color', [0.95, 0.97, 1.0], 'FontSize', 11.5);
    legend(ax2, 'Location', 'southwest', 'TextColor', [0.92, 0.95, 1.0], 'Color', [0.08, 0.11, 0.18]);
    styleDarkAxes(ax2);

    % ---------------------------------------------------------------------
    % Panel 3: Effect of Waveform Bandwidth B on 1D Range Resolution
    % ---------------------------------------------------------------------
    ax3 = nexttile(tl, 3);
    hold(ax3, 'on');
    plot(ax3, rtiLowBW.rangeVec_m, rtiLowBW.RTI_total_dBsm(1, :), '-.', ...
        'Color', [0.98, 0.75, 0.20], 'LineWidth', 1.9, ...
        'DisplayName', sprintf('B = %d MHz (\\DeltaR = %.1f cm)', bwList(1), rtiLowBW.deltaR_m * 100));
    plot(ax3, rtiMedBW.rangeVec_m, rtiMedBW.RTI_total_dBsm(1, :), '-', ...
        'Color', [0.22, 0.74, 0.97], 'LineWidth', 2.0, ...
        'DisplayName', sprintf('B = %d MHz (\\DeltaR = %.1f cm)', bwList(2), rtiMedBW.deltaR_m * 100));
    plot(ax3, rtiHighBW.rangeVec_m, rtiHighBW.RTI_total_dBsm(1, :), '-', ...
        'Color', [0.75, 0.52, 0.98], 'LineWidth', 1.6, ...
        'DisplayName', sprintf('B = %d MHz (\\DeltaR = %.1f cm)', bwList(3), rtiHighBW.deltaR_m * 100));

    xlim(ax3, [min(rtiMain.rangeVec_m), max(rtiMain.rangeVec_m)]);
    ylim(ax3, [-42, 22]);
    xlabel(ax3, 'Down-Range Coordinate R [m]  (Relative to Vehicle Center)', 'Color', [0.92, 0.95, 1.0], 'FontWeight', 'bold');
    ylabel(ax3, 'Compressed Range Profile \sigma(R) [dBsm]', 'Color', [0.92, 0.95, 1.0], 'FontWeight', 'bold');
    title(ax3, sprintf('3. Waveform Bandwidth Impact on Down-Range Scatterer Separation (\\psi = %.0f^\\circ, \\theta_{el} = %.0f^\\circ)', ...
        azCut, elLook), 'Color', [0.95, 0.97, 1.0], 'FontSize', 11.5);
    legend(ax3, 'Location', 'northeast', 'TextColor', [0.92, 0.95, 1.0], 'Color', [0.08, 0.11, 0.18]);
    styleDarkAxes(ax3);

    % ---------------------------------------------------------------------
    % Panel 4: 1D HRRP Single-Bounce vs. Concave Double-Bounce Decomposition
    % ---------------------------------------------------------------------
    ax4 = nexttile(tl, 4);
    hold(ax4, 'on');
    [~, idxCut] = min(abs(rtiMain.azVecDeg - azCut));
    plot(ax4, rtiMain.rangeVec_m, rtiMain.RTI_total_dBsm(idxCut, :), '-', ...
        'Color', [0.22, 0.74, 0.97], 'LineWidth', 2.3, ...
        'DisplayName', 'Total Compressed Return (PO + Concave 2-Bounce)');
    plot(ax4, rtiMain.rangeVec_m, rtiMain.RTI_single_dBsm(idxCut, :), '--', ...
        'Color', [0.98, 0.57, 0.24], 'LineWidth', 1.7, ...
        'DisplayName', 'Single-Bounce PO Only');
    plot(ax4, rtiMain.rangeVec_m, rtiMain.RTI_double_dBsm(idxCut, :), '-', ...
        'Color', [0.75, 0.52, 0.98], 'LineWidth', 1.9, ...
        'DisplayName', 'Concave Underside 2-Bounce (+0.5 d_{PQ} Delayed)');

    xline(ax4, rtiMain.R_nose_m(idxCut), ':c', 'Nose Tip', ...
        'LineWidth', 1.3, 'Color', [0.22, 0.95, 0.98], 'HandleVisibility', 'off');
    xline(ax4, rtiMain.R_waveriderBase_m(idxCut), ':y', 'Waverider Base (X=L)', ...
        'LineWidth', 1.3, 'Color', [0.98, 0.75, 0.20], 'HandleVisibility', 'off');

    xlim(ax4, [min(rtiMain.rangeVec_m), max(rtiMain.rangeVec_m)]);
    ylim(ax4, [-42, 22]);
    xlabel(ax4, 'Down-Range Coordinate R [m]  (Relative to Vehicle Center)', 'Color', [0.92, 0.95, 1.0], 'FontWeight', 'bold');
    ylabel(ax4, 'Range Profile Intensity \sigma(R) [dBsm]', 'Color', [0.92, 0.95, 1.0], 'FontWeight', 'bold');
    title(ax4, sprintf('4. 1D High-Resolution Range Profile Decomposition (B = %.0f MHz, \\psi = %.0f^\\circ)', ...
        bw, azCut), 'Color', [0.95, 0.97, 1.0], 'FontSize', 11.5);
    legend(ax4, 'Location', 'northeast', 'TextColor', [0.92, 0.95, 1.0], 'Color', [0.08, 0.11, 0.18]);
    styleDarkAxes(ax4);

    if ~isempty(opts.SavePath)
        exportgraphics(fig, opts.SavePath, 'Resolution', 200, 'BackgroundColor', [0.05, 0.07, 0.12]);
    end
end

function styleDarkAxes(ax)
    set(ax, 'Color', [0.07, 0.10, 0.17], ...
        'XColor', [0.80, 0.85, 0.93], ...
        'YColor', [0.80, 0.85, 0.93], ...
        'GridColor', [0.35, 0.45, 0.60], ...
        'GridAlpha', 0.28, ...
        'Box', 'on', ...
        'FontSize', 10);
    grid(ax, 'on');
end
