function fig = plot_geometry_blueprint(presetName, visibleMode)
%PLOT_GEOMETRY_BLUEPRINT Generate and save the 4-panel geometric lofting blueprint.
%   fig = PLOT_GEOMETRY_BLUEPRINT(presetName, visibleMode) builds a
%   watertight waverider body (default: 'osculating_cone') and plots:
%     - 3D Parametric Lofting Frame (streamlines, rib stations, LE, volume centroid)
%     - Multi-Station Transverse Cross-Sections at x/L = [0.25, 0.50, 0.75, 1.00]
%     - Planform Local Thickness Contour Map & Streamlines
%     - Longitudinal Area-Rule A(x/L), Cumulative Volume V(x/L), & Max Thickness

    if nargin < 1 || isempty(presetName)
        presetName = 'osculating_cone';
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
    fig = body.plotGeometryLines('Visible', visibleMode);

    outPath = fullfile(outDir, sprintf('geometry_blueprint_%s.png', presetName));
    exportgraphics(fig, outPath, 'Resolution', 180);
    fprintf('Saved Geometry Lofting Blueprint figure to: %s\n', outPath);
end
