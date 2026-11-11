function [figComp, figTrade] = plot_family_comparison(visibleMode)
%PLOT_FAMILY_COMPARISON Compare all 4 waverider families and plot 2D design trade space.
%   [figComp, figTrade] = PLOT_FAMILY_COMPARISON(visibleMode) generates:
%     1. Side-by-side comparison across 'conical', 'osculating_cone', 'caret',
%        and 'cst' waverider families (base cross-sections, area progression,
%        AoA L/D polars, and multi-attribute bar chart)
%     2. 2D Parametric Design Trade-Space Contours & Pareto Front

    if nargin < 1 || isempty(visibleMode)
        visibleMode = 'on';
    end

    rootDir = fileparts(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    outDir = fullfile(rootDir, 'output');
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    presets = {'conical', 'osculating_cone', 'caret', 'cst'};
    bodies = cell(1, numel(presets));
    for k = 1:numel(presets)
        bodies{k} = waverider.WaveriderBody(presets{k});
    end

    figComp = waverider.viz.plotFamilyComparison(bodies, 'Visible', visibleMode);
    compPath = fullfile(outDir, 'waverider_family_comparison.png');
    exportgraphics(figComp, compPath, 'Resolution', 180);
    fprintf('Saved 4-Family Comparison figure to: %s\n', compPath);

    figTrade = waverider.viz.plotTradeSpace('osculating_cone', 'Visible', visibleMode);
    tradePath = fullfile(outDir, 'waverider_trade_space.png');
    exportgraphics(figTrade, tradePath, 'Resolution', 180);
    fprintf('Saved Parametric Trade Space figure to: %s\n', tradePath);
end
