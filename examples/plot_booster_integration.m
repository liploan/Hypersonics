function [figBoost, figOverview, body] = plot_booster_integration(visibleMode)
%PLOT_BOOSTER_INTEGRATION Generate, verify, export, and plot a booster-mounted waverider.
%   [figBoost, figOverview, body] = PLOT_BOOSTER_INTEGRATION(visibleMode)
%   generates an Osculating Cones waverider fitted with a conformal aft
%   transition frustum and conical rocket booster stage ('booster_mounted'
%   preset), verifies 2-manifold watertightness, exports STL/OBJ meshes,
%   and saves the booster integration figures to ./output/.

    if nargin < 1 || isempty(visibleMode)
        visibleMode = 'on';
    end

    rootDir = fileparts(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    outDir = fullfile(rootDir, 'output');
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    body = waverider.WaveriderBody('booster_mounted');
    body.printSummary();

    % Export watertight STL and grouped OBJ for the integrated stack
    stlPath = fullfile(outDir, 'waverider_booster_mounted.stl');
    objPath = fullfile(outDir, 'waverider_booster_mounted.obj');
    body.exportSTL(stlPath);
    body.exportOBJ(objPath);
    fprintf('Exported Booster-Mounted Watertight STL: %s\n', stlPath);
    fprintf('Exported Booster-Mounted Grouped OBJ   : %s\n', objPath);

    % 1. Dedicated 4-panel Booster Integration Figure
    figBoost = body.plotBoosterIntegration('Visible', visibleMode);
    boostPng = fullfile(outDir, 'waverider_booster_integration.png');
    exportgraphics(figBoost, boostPng, 'Resolution', 180);
    fprintf('Saved Booster Integration figure to    : %s\n', boostPng);

    % 2. Standard 4-panel Overview Figure with Frustum & Rocket Stage
    figOverview = body.plot('Visible', visibleMode);
    ovPng = fullfile(outDir, 'waverider_booster_mounted_overview.png');
    exportgraphics(figOverview, ovPng, 'Resolution', 180);
    fprintf('Saved Booster Overview figure to       : %s\n', ovPng);
end
