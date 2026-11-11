function filePath = exportOBJ(meshData, filePath)
%EXPORTOBJ Export watertight waverider triangulation to a Wavefront OBJ file.
%   filePath = EXPORTOBJ(meshData, filePath) writes the global vertex list,
%   vertex normals, and grouped faces (UpperSurface, LowerCompressionSurface,
%   ConformalTransitionFrustum, ConicalRocketBooster, AftBaseCap) to a .obj file.

    [folder, ~, ~] = fileparts(filePath);
    if ~isempty(folder) && ~isfolder(folder)
        mkdir(folder);
    end

    TR = meshData.TR;
    V = TR.Points;
    VN = vertexNormal(TR);

    fid = fopen(filePath, 'w');
    if fid < 0
        error('exportOBJ:FileOpenError', 'Could not open "%s" for writing.', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid));

    fprintf(fid, '# Watertight Hypersonic Waverider Glide Body\n');
    fprintf(fid, '# Vertices: %d | Faces: %d\n\n', size(V, 1), size(TR.ConnectivityList, 1));

    fprintf(fid, 'v %.8f %.8f %.8f\n', V.');
    fprintf(fid, '\n');

    fprintf(fid, 'vn %.6f %.6f %.6f\n', VN.');
    fprintf(fid, '\n');

    writeGroup(fid, 'UpperSurface', meshData.FacesUpper);
    writeGroup(fid, 'LowerCompressionSurface', meshData.FacesLower);
    if isfield(meshData, 'FacesFrustum') && ~isempty(meshData.FacesFrustum)
        writeGroup(fid, 'ConformalTransitionFrustum', meshData.FacesFrustum);
    end
    if isfield(meshData, 'FacesRocket') && ~isempty(meshData.FacesRocket)
        writeGroup(fid, 'ConicalRocketBooster', meshData.FacesRocket);
    end
    writeGroup(fid, 'AftBaseCap', meshData.FacesBase);
end

function writeGroup(fid, groupName, F)
    fprintf(fid, 'g %s\n', groupName);
    for i = 1:size(F, 1)
        v1 = F(i, 1);
        v2 = F(i, 2);
        v3 = F(i, 3);
        fprintf(fid, 'f %d//%d %d//%d %d//%d\n', v1, v1, v2, v2, v3, v3);
    end
    fprintf(fid, '\n');
end
