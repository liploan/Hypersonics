function bodies = demo_all_waveriders(saveOutputs)
%DEMO_ALL_WAVERIDERS Generate, verify, visualize, and export all 4 waverider families.
%   bodies = DEMO_ALL_WAVERIDERS(saveOutputs) generates canonical watertight
%   hypersonic waverider bodies for:
%     1. Axisymmetric Cone-Derived ('conical')
%     2. Sobieczky Osculating Cones ('osculating_cone')
%     3. Nonweiler Caret / 2D Wedge ('caret')
%     4. Kulfan CST Boost-Glide Waverider ('cst')
%   and exports STL/OBJ meshes and overview figures into ./output/.

    if nargin < 1
        saveOutputs = true;
    end

    rootDir = fileparts(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    outDir = fullfile(rootDir, 'output');
    if saveOutputs && ~isfolder(outDir)
        mkdir(outDir);
    end

    presets = {'conical', 'osculating_cone', 'caret', 'cst'};
    bodies = cell(1, numel(presets));

    for k = 1:numel(presets)
        name = presets{k};
        body = waverider.WaveriderBody(name);
        body.printSummary();
        bodies{k} = body;

        if saveOutputs
            stlPath = fullfile(outDir, sprintf('waverider_%s.stl', name));
            objPath = fullfile(outDir, sprintf('waverider_%s.obj', name));
            body.exportSTL(stlPath);
            body.exportOBJ(objPath);
            fprintf('  -> Exported STL: %s\n', stlPath);
            fprintf('  -> Exported OBJ: %s\n', objPath);

            fig = body.plot('Visible', 'off');
            pngPath = fullfile(outDir, sprintf('waverider_%s_overview.png', name));
            exportgraphics(fig, pngPath, 'Resolution', 180);
            close(fig);
            fprintf('  -> Saved Figure: %s\n', pngPath);
        end
    end
end
