function fig = plot_key_attributes(presetName, visibleMode)
%PLOT_KEY_ATTRIBUTES Generate and save the 4-panel aerodynamic & flowfield attribute dashboard.
%   fig = PLOT_KEY_ATTRIBUTES(presetName, visibleMode) builds a watertight
%   waverider body (default: 'osculating_cone') and plots:
%     - Taylor-Maccoll Conical Shock-Layer Profile (Cp, Mach, |u_theta| vs ray angle theta)
%     - Angle-of-Attack Aerodynamic Polar (CL, CD_wave, and Inviscid L/D vs alpha)
%     - Spanwise Chord, Base Height & Sectional Lift Loading across y/(b/2)
%     - 3D Surface Pressure Cp & Verified Outward Manifold Normals Quiver Plot

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
    fig = body.plotKeyAttributes('Visible', visibleMode);

    outPath = fullfile(outDir, sprintf('key_attributes_%s.png', presetName));
    waverider.viz.exportFigurePNG(fig, outPath, 200);
    fprintf('Saved Key Attributes Dashboard figure to: %s\n', outPath);
end
