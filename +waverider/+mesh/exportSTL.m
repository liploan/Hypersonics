function filePath = exportSTL(meshData, filePath, mode)
%EXPORTSTL Export watertight waverider triangulation to an STL file.
%   filePath = EXPORTSTL(meshData, filePath, mode) writes the MATLAB
%   triangulation to disk in 'binary' (default) or 'text' (ASCII) STL format.

    if nargin < 3 || isempty(mode)
        mode = 'binary';
    end

    [folder, ~, ~] = fileparts(filePath);
    if ~isempty(folder) && ~isfolder(folder)
        mkdir(folder);
    end

    stlwrite(meshData.TR, filePath, mode);
end
