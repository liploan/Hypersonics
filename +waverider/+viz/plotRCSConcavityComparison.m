function [fig, compData] = plotRCSConcavityComparison(baseBody, varargin)
%PLOTRCSCONCAVITYCOMPARISON 2D Az-El RCS maps, family comparison & concavity trade-off.
%   [fig, compData] = PLOTRCSCONCAVITYCOMPARISON(baseBody, Name, Value)
%   generates a 4-panel analysis focused on the electromagnetic impact of
%   underside concavity:
%     1. 2D Azimuth x Elevation Total RCS Heatmap at X-Band (10 GHz)
%     2. 2D Azimuth x Elevation Concave Double-Bounce Gain Map (dB)
%     3. 4-Family Ventral RCS Comparison (Conical, Osculating Cone, Caret, CST)
%     4. Parametric Sweep of Underside Anhedral/Concavity Depth vs RCS & L/D

    p = inputParser;
    addParameter(p, 'Visible', 'on', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    opts = p.Results;

    if nargin < 1 || isempty(baseBody)
        baseBody = waverider.WaveriderBody('caret');
    end

    cpBase = waverider.rcs.buildConcavePairList(baseBody);

    % 1. 2D Azimuth x Elevation Grid at X-Band (10 GHz)
    azGridVec = linspace(0.0, 180.0, 46);
    elGridVec = linspace(-75.0, 75.0, 39);
    [AZ, EL] = meshgrid(azGridVec, elGridVec);
    rcs2D = waverider.rcs.computeRCS(baseBody, 10.0, AZ(:), EL(:), ...
        'PrecomputedConcavePairs', cpBase);

    sigma2D_total = reshape(rcs2D.SigmaTotal_dBsm, size(AZ));
    gain2D_mb     = reshape(rcs2D.MultiBounceGain_dB, size(AZ));

    % 2. Compare all 4 Waverider Families in Ventral Sector (el = -22 deg, X-band 10 GHz)
    presets = {'conical', 'osculating_cone', 'caret', 'cst'};
    famLabels = {'CONICAL', 'OSCULATING CONE', 'CARET (Inverted-V)', 'CST (Concave Arch)'};
    azFine = linspace(0.0, 180.0, 121).';
    famRCS_X = zeros(numel(azFine), numel(presets));
    famPairs = zeros(1, numel(presets));

    for k = 1:numel(presets)
        bFam = waverider.WaveriderBody(presets{k});
        cpFam = waverider.rcs.buildConcavePairList(bFam);
        famPairs(k) = cpFam.numPairs;
        resFam = waverider.rcs.computeRCS(bFam, 10.0, azFine, -22.0, ...
            'PrecomputedConcavePairs', cpFam);
        famRCS_X(:, k) = movmean(resFam.SigmaTotal_dBsm, 3);
    end

    % 3. Sweep Underside Anhedral / Concavity Depth (0 deg flat -> 18 deg deep trough)
    anhedralVec = linspace(0.0, 18.0, 8);
    meanVentralX_dBsm   = zeros(size(anhedralVec));
    meanVentralUHF_dBsm = zeros(size(anhedralVec));
    mbGainPeak_dB       = zeros(size(anhedralVec));
    aeroLD_vec          = zeros(size(anhedralVec));

    azSector = linspace(25.0, 75.0, 15).';
    elSector = repmat(-25.0, size(azSector));

    for aIdx = 1:numel(anhedralVec)
        cfgSweep = waverider.WaveriderConfig.fromPreset('osculating_cone');
        cfgSweep.LeadingEdgeAnhedralDeg = anhedralVec(aIdx);
        cfgSweep.NumStreamwise = 40;
        cfgSweep.NumHalfSpan = 28;
        bSwp = waverider.WaveriderBody(cfgSweep);
        cpSwp = waverider.rcs.buildConcavePairList(bSwp);

        rcsSwp = waverider.rcs.computeRCS(bSwp, [0.5, 10.0], azSector, elSector, ...
            'PrecomputedConcavePairs', cpSwp);
        meanVentralUHF_dBsm(aIdx) = 10.0 * log10(mean(rcsSwp.SigmaTotal_m2(:, 1)));
        meanVentralX_dBsm(aIdx)   = 10.0 * log10(mean(rcsSwp.SigmaTotal_m2(:, 2)));
        mbGainPeak_dB(aIdx)       = max(rcsSwp.MultiBounceGain_dB(:, 2));
        aeroLD_vec(aIdx)          = bSwp.AeroMetrics.LD_inviscid;
    end

    fig = figure('Color', 'w', 'Position', [65, 65, 1380, 860], 'Visible', opts.Visible);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Underside Concavity Electromagnetic Impact: 2D Az-El Maps, Family Comparison & Trade-Offs', ...
        'FontSize', 14, 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Tile 1: 2D Azimuth x Elevation Total RCS Heatmap (X-Band 10 GHz)
    % ---------------------------------------------------------------------
    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    contourf(ax1, AZ, EL, sigma2D_total, 18, 'LineColor', 'none', 'HandleVisibility', 'off');
    colormap(ax1, 'turbo');
    cb1 = colorbar(ax1);
    cb1.Label.String = 'Total Monostatic RCS \sigma [dBsm] (10 GHz)';
    yline(ax1, 0.0, 'w--', 'Horizon (\theta_{el} = 0^\circ)', 'LineWidth', 1.6);
    xlim(ax1, [0, 180]); ylim(ax1, [-75, 75]);
    xlabel(ax1, 'Azimuth Angle \psi [deg] (0^\circ = Nose, 90^\circ = Broadside, 180^\circ = Tail)');
    ylabel(ax1, 'Elevation \theta_{el} [deg] (<0^\circ Ventral | >0^\circ Dorsal)');
    title(ax1, sprintf('2D Azimuth-Elevation RCS Sphere [%s, X-Band 10 GHz]', ...
        upper(strrep(baseBody.Config.Method, '_', ' '))), 'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Tile 2: 2D Azimuth x Elevation Concave Double-Bounce Gain Map (dB)
    % ---------------------------------------------------------------------
    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    contourf(ax2, AZ, EL, gain2D_mb, 16, 'LineColor', 'none', 'HandleVisibility', 'off');
    cb2 = colorbar(ax2);
    cb2.Label.String = 'Concave Multi-Bounce Amplification \Delta\sigma_{MB} [dB]';
    yline(ax2, 0.0, 'w--', 'Horizon (\theta_{el} = 0^\circ)', 'LineWidth', 1.6);
    xlim(ax2, [0, 180]); ylim(ax2, [-75, 75]);
    xlabel(ax2, 'Azimuth Angle \psi [deg]');
    ylabel(ax2, 'Elevation \theta_{el} [deg] (Note Ventral Amplification for \theta_{el} < 0^\circ)');
    title(ax2, 'Concave Underside Double-Bounce Enhancement (\sigma_{total} - \sigma_{single})', ...
        'FontWeight', 'bold');

    % ---------------------------------------------------------------------
    % Tile 3: 4-Family Ventral RCS Comparison (X-Band 10 GHz, el = -22 deg)
    % ---------------------------------------------------------------------
    ax3 = nexttile(tl, 3);
    hold(ax3, 'on');
    famCols = [0.12, 0.46, 0.78; 0.88, 0.36, 0.14; 0.18, 0.66, 0.34; 0.58, 0.26, 0.72];
    for k = 1:numel(presets)
        plot(ax3, azFine, famRCS_X(:, k), '-', 'Color', famCols(k, :), 'LineWidth', 2.2, ...
            'DisplayName', sprintf('%s (%d concave pairs)', famLabels{k}, famPairs(k)));
    end
    grid(ax3, 'on');
    xlim(ax3, [0, 180]);
    xlabel(ax3, 'Azimuth Angle \psi [deg] (Ventral Look-Up \theta_{el} = -22^\circ)');
    ylabel(ax3, 'X-Band Monostatic RCS \sigma [dBsm]');
    title(ax3, 'Cross-Family Ventral RCS Comparison at X-Band (10 GHz)', 'FontWeight', 'bold');
    legend(ax3, 'Location', 'northwest', 'FontSize', 8);

    % ---------------------------------------------------------------------
    % Tile 4: Underside Concavity Anhedral Depth vs. Ventral RCS & Aero L/D
    % ---------------------------------------------------------------------
    ax4 = nexttile(tl, 4);
    yyaxis(ax4, 'left');
    hold(ax4, 'on');
    plot(ax4, anhedralVec, meanVentralX_dBsm, '-s', 'Color', [0.85, 0.22, 0.15], ...
        'LineWidth', 2.3, 'MarkerSize', 6, 'DisplayName', 'Ventral Sector Mean RCS (X-Band 10 GHz)');
    plot(ax4, anhedralVec, meanVentralUHF_dBsm, '--o', 'Color', [0.15, 0.48, 0.82], ...
        'LineWidth', 2.1, 'MarkerSize', 6, 'DisplayName', 'Ventral Sector Mean RCS (UHF 0.5 GHz)');
    plot(ax4, anhedralVec, mbGainPeak_dB, ':^', 'Color', [0.58, 0.20, 0.72], ...
        'LineWidth', 2.0, 'MarkerSize', 6, 'DisplayName', 'Peak Multi-Bounce Gain \Delta\sigma_{MB} [dB]');
    ylabel(ax4, 'Sector RCS [dBsm] / Multi-Bounce Gain [dB]');

    yyaxis(ax4, 'right');
    plot(ax4, anhedralVec, aeroLD_vec, '-d', 'Color', [0.18, 0.65, 0.35], ...
        'LineWidth', 2.3, 'MarkerSize', 6, 'DisplayName', 'Inviscid Aerodynamic L/D');
    ylabel(ax4, 'Inviscid Aerodynamic L/D [-]');

    grid(ax4, 'on');
    xlabel(ax4, 'Underside Leading-Edge Anhedral / Concave Trough Angle \Gamma [deg]');
    title(ax4, 'Concavity Trade-Off: Ventral Cavity RCS vs. Hypersonic L/D', 'FontWeight', 'bold');
    legend(ax4, 'Location', 'northwest', 'FontSize', 8);

    compData = struct();
    compData.AZ = AZ;
    compData.EL = EL;
    compData.Sigma2D_dBsm = sigma2D_total;
    compData.Gain2D_dB = gain2D_mb;
    compData.AnhedralDegVec = anhedralVec;
    compData.MeanVentralX_dBsm = meanVentralX_dBsm;
    compData.MeanVentralUHF_dBsm = meanVentralUHF_dBsm;
    compData.PeakMBGain_dB = mbGainPeak_dB;
    compData.AeroLD = aeroLD_vec;
end
