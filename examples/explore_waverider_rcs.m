function [figBands, figConcavity, rcsReport] = explore_waverider_rcs(presetName, visibleMode)
%EXPLORE_WAVERIDER_RCS Explore UHF to X-band RCS & underside concave cavity scattering.
%   [figBands, figConcavity, rcsReport] = EXPLORE_WAVERIDER_RCS(presetName, visibleMode)
%   evaluates monostatic RCS across UHF (0.5 GHz), L (1.5 GHz), S (3.0 GHz),
%   C (6.0 GHz), and X (10.0 GHz) bands, quantifies the ventral vs. dorsal
%   asymmetry and double-bounce amplification caused by the concave lower
%   compression surface, and exports two 4-panel figures to ./output/.

    if nargin < 1 || isempty(presetName)
        presetName = 'caret';
    end
    if nargin < 2 || isempty(visibleMode)
        visibleMode = 'on';
    end

    rootDir = fileparts(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    outDir = fullfile(rootDir, 'output');
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    body = waverider.WaveriderBody(presetName);
    cp = waverider.rcs.buildConcavePairList(body);

    bandFreqs = [0.5, 1.5, 3.0, 6.0, 10.0];
    bandNames = {'UHF (0.5 GHz)', 'L-Band (1.5 GHz)', 'S-Band (3.0 GHz)', 'C-Band (6.0 GHz)', 'X-Band (10.0 GHz)'};

    % Evaluate Ventral (el = -22 deg) vs Dorsal (el = +22 deg) across forward/oblique sector (az = 20..70 deg)
    azSec = linspace(20.0, 70.0, 26).';
    rcsVentral = waverider.rcs.computeRCS(body, bandFreqs, azSec, -22.0, 'PrecomputedConcavePairs', cp);
    rcsDorsal  = waverider.rcs.computeRCS(body, bandFreqs, azSec, +22.0, 'PrecomputedConcavePairs', cp);
    rcsNose    = waverider.rcs.computeRCS(body, bandFreqs, 0.0,   -10.0, 'PrecomputedConcavePairs', cp);

    fprintf('\n====================================================================================================\n');
    fprintf('  UHF TO X-BAND RCS & CONCAVE UNDERSIDE SCATTERING REPORT [%s | %d Concave Bounce Pairs]\n', ...
        upper(presetName), cp.numPairs);
    fprintf('====================================================================================================\n');
    fprintf('  %-18s %-12s %-14s %-15s %-15s %-15s\n', ...
        'Radar Band', 'Lambda [m]', 'Nose-On[dBsm]', 'Dorsal(+22el)', 'Ventral(-22el)', 'Max MB Gain[dB]');
    fprintf('----------------------------------------------------------------------------------------------------\n');

    for b = 1:numel(bandFreqs)
        nose_dB    = rcsNose.SigmaTotal_dBsm(1, b);
        dorsal_dB  = 10.0 * log10(mean(rcsDorsal.SigmaTotal_m2(:, b)));
        ventral_dB = 10.0 * log10(mean(rcsVentral.SigmaTotal_m2(:, b)));
        mbGain_dB  = max(rcsVentral.MultiBounceGain_dB(:, b));

        fprintf('  %-18s %-12.3f %+12.2f   %+13.2f   %+13.2f   %+13.2f\n', ...
            bandNames{b}, rcsVentral.Lambda_m(b), nose_dB, dorsal_dB, ventral_dB, mbGain_dB);
    end
    fprintf('====================================================================================================\n');

    % Render & Export Figure 1: Multi-Band UHF-to-X-Band RCS & 3D Concave Ray-Bounce Dashboard
    [figBands, summaryMetrics] = waverider.viz.plotRCSAnalysis(body, 'Visible', visibleMode);
    f1Path = fullfile(outDir, sprintf('waverider_rcs_uhf_to_xband_%s.png', presetName));
    exportgraphics(figBands, f1Path, 'Resolution', 180);
    fprintf('Saved Multi-Band RCS Analysis figure to    : %s\n', f1Path);

    % Render & Export Figure 2: 2D Az-El Sphere, Multi-Bounce Gain, 4-Family & Anhedral Trade-Off
    [figConcavity, compData] = waverider.viz.plotRCSConcavityComparison(body, 'Visible', visibleMode);
    f2Path = fullfile(outDir, 'waverider_rcs_family_concavity.png');
    exportgraphics(figConcavity, f2Path, 'Resolution', 180);
    fprintf('Saved Concavity Comparison & Trade figure  : %s\n', f2Path);

    rcsReport = struct('SummaryMetrics', summaryMetrics, 'ComparisonData', compData);
end
