function rootDir = setup_waverider()
%SETUP_WAVERIDER Initialize the Hypersonic Waverider MATLAB environment.
%   rootDir = SETUP_WAVERIDER() adds the project root, examples, and tests
%   directories to the MATLAB path and displays available presets.

    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir);
    if isfolder(fullfile(rootDir, 'examples'))
        addpath(fullfile(rootDir, 'examples'));
    end
    if isfolder(fullfile(rootDir, 'tests'))
        addpath(fullfile(rootDir, 'tests'));
    end

    fprintf('================================================================\n');
    fprintf('  Hypersonic Parametric Watertight Waverider Generator (MATLAB) \n');
    fprintf('  Root directory added to path: %s\n', rootDir);
    fprintf('================================================================\n');
end
