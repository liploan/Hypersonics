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
    addParameter(p, 'WindowType', 'hann', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Polarization', 'Mean', @(x) ischar(x) || isstring(x));
    addParameter(p, 'IncludeEdgeDiffraction', false, @islogical);
    addParameter(p, 'Visible', 'on', @(x) ischar(x) || isstring(x));
    addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    opts = p.Results;

    fc = opts.CenterFreqGHz;
    bw = opts.BandwidthMHz;
    elLook = opts.ElevationDeg;
    azCut = opts.CutAzimuthDeg;
    winType = opts.WindowType;
    polIn = opts.Polarization;
    extraArgs = {};
    if opts.IncludeEdgeDiffraction
        extraArgs = {'IncludeEdgeDiffraction', true};
    end

    azVec = unique([0:1.5:171, 172:0.5:180]);
    pairs = waverider.rcs.buildConcavePairList(body);

    % Compute primary RTI at configured bandwidth
    rtiMain = waverider.rcs.computeRangeProfile(body, fc, bw, azVec, elLook, ...
        'IncludeMultiBounce', true, 'WindowType', winType, 'Polarization', polIn, ...
        'PrecomputedConcavePairs', pairs, extraArgs{:});

    % Compute multi-bandwidth HRRP cuts at azCut
    bwList = [150, 600, 2000];
    rtiLowBW  = waverider.rcs.computeRangeProfile(body, fc, bwList(1), azCut, elLook, ...
        'IncludeMultiBounce', true, 'WindowType', winType, 'Polarization', polIn, ...
        'PrecomputedConcavePairs', pairs, extraArgs{:});
    rtiMedBW  = waverider.rcs.computeRangeProfile(body, fc, bwList(2), azCut, elLook, ...
        'IncludeMultiBounce', true, 'WindowType', winType, 'Polarization', polIn, ...
        'PrecomputedConcavePairs', pairs, extraArgs{:});
    rtiHighBW = waverider.rcs.computeRangeProfile(body, fc, bwList(3), azCut, elLook, ...
        'IncludeMultiBounce', true, 'WindowType', winType, 'Polarization', polIn, ...
        'PrecomputedConcavePairs', pairs, extraArgs{:});
    peakMainDB = max(rtiMain.RTI_total_dBsm(:));
    cMaxRTI = max(18, min(55, ceil((peakMainDB + 2) / 5) * 5));
    peakCutDB = max([rtiLowBW.RTI_total_dBsm(:); rtiMedBW.RTI_total_dBsm(:); rtiHighBW.RTI_total_dBsm(:)]);
    yMaxHRRP = max(22, min(55, ceil((peakCutDB + 3) / 5) * 5));

    fig = figure('Name', 'Waverider Range-Angle-Intensity & Waveform Bandwidth Analysis', ...
        'Color', 'w', 'Position', [60, 60, 1460, 900], 'Visible', opts.Visible);
    waverider.viz.stylePublicationFigure(fig, 'init');

    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
    title(tl, sprintf('Range-Angle-Intensity (Sinogram) & Waveform Bandwidth Analysis  (f_c = %.1f GHz, B = %.0f MHz, \\DeltaR = %.1f cm, \\theta_{el} = %.0f^\\circ, Pol: %s)', ...
        fc, bw, rtiMain.deltaR_m * 100, elLook, rtiMain.polarization), ...
        'FontSize', 14, 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Panel 1: Total Range-Angle-Intensity Sinogram (Single + Double Bounce)
    % ---------------------------------------------------------------------
    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    imagesc(ax1, rtiMain.rangeVec_m, rtiMain.azVecDeg, rtiMain.RTI_total_dBsm);
    set(ax1, 'YDir', 'normal', 'Box', 'on');
    colormap(ax1, turbo);
    clim(ax1, [-38, cMaxRTI]);
    cb1 = colorbar(ax1);
    cb1.Label.String = 'Scattering Intensity \sigma(R, \psi) [dBsm]';

    % Overlay geometric Nose and Base trajectories
    plot(ax1, rtiMain.R_nose_m, rtiMain.azVecDeg, '--', 'Color', [0.88, 0.98, 1.00], ...
        'LineWidth', 1.8, 'DisplayName', 'Nose Tip Locus R_{nose}(\psi)');
    plot(ax1, rtiMain.R_waveriderBase_m, rtiMain.azVecDeg, '-.', 'Color', [1.00, 0.84, 0.25], ...
        'LineWidth', 1.7, 'DisplayName', 'Waverider Base (X = L)');
    if max(abs(rtiMain.R_aftBase_m - rtiMain.R_waveriderBase_m)) > 0.05
        plot(ax1, rtiMain.R_aftBase_m, rtiMain.azVecDeg, ':', 'Color', [0.55, 1.00, 0.65], ...
            'LineWidth', 1.8, 'DisplayName', 'Booster Aft Base');
    end
    yline(ax1, azCut, '--w', sprintf('HRRP Cut (\\psi = %.0f^\\circ)', azCut), ...
        'LineWidth', 1.4, 'LabelHorizontalAlignment', 'left', ...
        'Color', [1.0, 1.0, 1.0], 'HandleVisibility', 'off');

    xlim(ax1, [min(rtiMain.rangeVec_m), max(rtiMain.rangeVec_m)]);
    ylim(ax1, [0, 180]);
    xlabel(ax1, 'Down-Range Coordinate R [m]  (Relative to Vehicle Center)');
    ylabel(ax1, 'Azimuth Aspect Angle \psi [deg]  (0^\circ=Nose, 90^\circ=Broadside, 180^\circ=Tail)');
    title(ax1, sprintf('1. Total Range-Angle-Intensity Sinogram (B = %.0f MHz, \\DeltaR = %.1f cm)', ...
        bw, rtiMain.deltaR_m * 100), 'FontWeight', 'bold');
    legend(ax1, 'Location', 'southwest', 'FontSize', 8.5);

    % ---------------------------------------------------------------------
    % Panel 2: Concave Cavity Double-Bounce Range-Angle Isolation
    % ---------------------------------------------------------------------
    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    imagesc(ax2, rtiMain.rangeVec_m, rtiMain.azVecDeg, rtiMain.RTI_double_dBsm);
    set(ax2, 'YDir', 'normal', 'Box', 'on');
    colormap(ax2, turbo);
    clim(ax2, [-38, cMaxRTI]);
    cb2 = colorbar(ax2);
    cb2.Label.String = 'Concave 2-Bounce Intensity [dBsm]';

    plot(ax2, rtiMain.R_nose_m, rtiMain.azVecDeg, '--', 'Color', [0.88, 0.98, 1.00], ...
        'LineWidth', 1.7, 'DisplayName', 'Nose Tip Locus');
    plot(ax2, rtiMain.R_waveriderBase_m, rtiMain.azVecDeg, '-.', 'Color', [1.00, 0.84, 0.25], ...
        'LineWidth', 1.7, 'DisplayName', 'Waverider Base (X = L)');

    xlim(ax2, [min(rtiMain.rangeVec_m), max(rtiMain.rangeVec_m)]);
    ylim(ax2, [0, 180]);
    xlabel(ax2, 'Down-Range Coordinate R [m]  (+0.5 d_{PQ} Cavity Delay Visible)');
    ylabel(ax2, 'Azimuth Aspect Angle \psi [deg]');
    title(ax2, '2. Isolated Concave Underside Double-Bounce Sinogram (Delayed Cavity Returns)', ...
        'FontWeight', 'bold');
    legend(ax2, 'Location', 'southwest', 'FontSize', 8.5);

    % ---------------------------------------------------------------------
    % Panel 3: Effect of Waveform Bandwidth B on 1D Range Resolution
    % ---------------------------------------------------------------------
    ax3 = nexttile(tl, 3);
    hold(ax3, 'on');
    plot(ax3, rtiLowBW.rangeVec_m, rtiLowBW.RTI_total_dBsm(1, :), '-.', ...
        'Color', [0.85, 0.36, 0.08], 'LineWidth', 2.2, ...
        'DisplayName', sprintf('B = %d MHz (\\DeltaR = %.1f cm)', bwList(1), rtiLowBW.deltaR_m * 100));
    plot(ax3, rtiMedBW.rangeVec_m, rtiMedBW.RTI_total_dBsm(1, :), '-', ...
        'Color', [0.06, 0.42, 0.78], 'LineWidth', 2.2, ...
        'DisplayName', sprintf('B = %d MHz (\\DeltaR = %.1f cm)', bwList(2), rtiMedBW.deltaR_m * 100));
    plot(ax3, rtiHighBW.rangeVec_m, rtiHighBW.RTI_total_dBsm(1, :), '-', ...
        'Color', [0.52, 0.16, 0.72], 'LineWidth', 1.8, ...
        'DisplayName', sprintf('B = %d MHz (\\DeltaR = %.1f cm)', bwList(3), rtiHighBW.deltaR_m * 100));

    grid(ax3, 'on');
    set(ax3, 'Box', 'on');
    xlim(ax3, [min(rtiMain.rangeVec_m), max(rtiMain.rangeVec_m)]);
    ylim(ax3, [-42, yMaxHRRP]);
    xlabel(ax3, 'Down-Range Coordinate R [m]  (Relative to Vehicle Center)');
    ylabel(ax3, 'Compressed Range Profile \sigma(R) [dBsm]');
    title(ax3, sprintf('3. Waveform Bandwidth Impact on Down-Range Scatterer Separation (\\psi = %.0f^\\circ, \\theta_{el} = %.0f^\\circ)', ...
        azCut, elLook), 'FontWeight', 'bold');
    legend(ax3, 'Location', 'northeast', 'FontSize', 8.5);

    % ---------------------------------------------------------------------
    % Panel 4: 1D HRRP Single-Bounce vs. Concave Double-Bounce Decomposition
    % ---------------------------------------------------------------------
    ax4 = nexttile(tl, 4);
    hold(ax4, 'on');
    [~, idxCut] = min(abs(rtiMain.azVecDeg - azCut));
    plot(ax4, rtiMain.rangeVec_m, rtiMain.RTI_total_dBsm(idxCut, :), '-', ...
        'Color', [0.06, 0.42, 0.78], 'LineWidth', 2.4, ...
        'DisplayName', 'Total Compressed Return (PO + Concave 2-Bounce)');
    plot(ax4, rtiMain.rangeVec_m, rtiMain.RTI_single_dBsm(idxCut, :), '--', ...
        'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.9, ...
        'DisplayName', 'Single-Bounce PO Only');
    plot(ax4, rtiMain.rangeVec_m, rtiMain.RTI_double_dBsm(idxCut, :), '-', ...
        'Color', [0.52, 0.16, 0.72], 'LineWidth', 2.0, ...
        'DisplayName', 'Concave Underside 2-Bounce (+0.5 d_{PQ} Delayed)');

    xline(ax4, rtiMain.R_nose_m(idxCut), ':', 'Nose Tip', ...
        'LineWidth', 1.4, 'Color', [0.05, 0.36, 0.56], ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    xline(ax4, rtiMain.R_waveriderBase_m(idxCut), ':', 'Waverider Base (X=L)', ...
        'LineWidth', 1.4, 'Color', [0.68, 0.28, 0.05], ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    grid(ax4, 'on');
    set(ax4, 'Box', 'on');
    xlim(ax4, [min(rtiMain.rangeVec_m), max(rtiMain.rangeVec_m)]);
    ylim(ax4, [-42, yMaxHRRP]);
    xlabel(ax4, 'Down-Range Coordinate R [m]  (Relative to Vehicle Center)');
    ylabel(ax4, 'Range Profile Intensity \sigma(R) [dBsm]');
    title(ax4, sprintf('4. 1D High-Resolution Range Profile Decomposition (B = %.0f MHz, \\psi = %.0f^\\circ)', ...
        bw, azCut), 'FontWeight', 'bold');
    legend(ax4, 'Location', 'northeast', 'FontSize', 8.5);

    waverider.viz.stylePublicationFigure(fig, 'finalize');

    if ~isempty(opts.SavePath)
        waverider.viz.exportFigurePNG(fig, opts.SavePath, 200);
    end
end
