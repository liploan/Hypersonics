function meshData = stitchWatertightMesh(surfData, cfg)
%STITCHWATERTIGHTMESH Assemble forebody, optional aft frustum/rocket, and base into a watertight mesh.
%   meshData = STITCHWATERTIGHTMESH(surfData, cfg) constructs a closed,
%   2-manifold, outward-oriented triangular surface mesh (MATLAB triangulation)
%   by sharing exact topological vertex IDs along the sharp leading edge,
%   wingtip corners, lateral frustum/rocket seams, and aft base perimeter.

    X_low = surfData.X_low;
    Y_low = surfData.Y_low;
    Z_low = surfData.Z_low;

    X_up = surfData.X_up;
    Y_up = surfData.Y_up;
    Z_up = surfData.Z_up;

    [NxTotal, Ny] = size(X_low);
    if isfield(surfData, 'NxWaverider')
        NxWave = surfData.NxWaverider;
        NxFrustEnd = surfData.NxFrustumEnd;
    else
        NxWave = NxTotal;
        NxFrustEnd = NxTotal;
    end
    Nb = max(cfg.NumBaseLayers, 2);

    % Construct structured aft base cap grid [Nb x Ny] at the final axial station (row NxTotal)
    tau = linspace(0.0, 1.0, Nb).';
    X_base = repmat(X_up(NxTotal, :), Nb, 1) .* (1.0 - tau) + repmat(X_low(NxTotal, :), Nb, 1) .* tau;
    Y_base = repmat(Y_up(NxTotal, :), Nb, 1) .* (1.0 - tau) + repmat(Y_low(NxTotal, :), Nb, 1) .* tau;
    Z_base = repmat(Z_up(NxTotal, :), Nb, 1) .* (1.0 - tau) + repmat(Z_low(NxTotal, :), Nb, 1) .* tau;

    % Preallocate vertex list and topological index maps
    maxVerts = 2 * NxTotal * Ny + Nb * Ny;
    vertices = zeros(maxVerts, 3);
    nVerts = 0;

    function id = addVertex(x, y, z)
        nVerts = nVerts + 1;
        vertices(nVerts, :) = [x, y, z];
        id = nVerts;
    end

    vid_up = zeros(NxTotal, Ny);
    vid_low = zeros(NxTotal, Ny);
    vid_base = zeros(Nb, Ny);

    % 1. Left and Right Wingtip Singular Vertices along Waverider Forebody (i = 1 .. NxWave)
    vLeftTip = addVertex(X_up(1, 1), Y_up(1, 1), Z_up(1, 1));
    vRightTip = addVertex(X_up(1, Ny), Y_up(1, Ny), Z_up(1, Ny));

    vid_up(1:NxWave, 1)  = vLeftTip;
    vid_low(1:NxWave, 1) = vLeftTip;
    vid_up(1:NxWave, Ny)  = vRightTip;
    vid_low(1:NxWave, Ny) = vRightTip;

    % 2. Shared Left/Right Lateral Seam Vertices along Conformal Frustum & Rocket Stage
    for i = (NxWave + 1):NxTotal
        vL = addVertex(X_up(i, 1), Y_up(i, 1), Z_up(i, 1));
        vR = addVertex(X_up(i, Ny), Y_up(i, Ny), Z_up(i, Ny));
        vid_up(i, 1)   = vL;
        vid_low(i, 1)  = vL;
        vid_up(i, Ny)  = vR;
        vid_low(i, Ny) = vR;
    end

    % Aft base cap left/right corners share the final station's left/right seam vertices
    vid_base(:, 1)  = vid_up(NxTotal, 1);
    vid_base(:, Ny) = vid_up(NxTotal, Ny);

    % 3. Interior Spanwise Stations (k = 2 .. Ny-1)
    for k = 2:(Ny - 1)
        % Shared Leading Edge vertex (i = 1)
        vLE = addVertex(X_up(1, k), Y_up(1, k), Z_up(1, k));
        vid_up(1, k) = vLE;
        vid_low(1, k) = vLE;

        % Upper surface interior & aft stations (i = 2 .. NxTotal)
        for i = 2:NxTotal
            vid_up(i, k) = addVertex(X_up(i, k), Y_up(i, k), Z_up(i, k));
        end

        % Lower surface interior & aft stations (i = 2 .. NxTotal)
        for i = 2:NxTotal
            vid_low(i, k) = addVertex(X_low(i, k), Y_low(i, k), Z_low(i, k));
        end

        % Base cap top edge shares aft-most Upper Trailing Edge (m = 1 -> i = NxTotal)
        vid_base(1, k) = vid_up(NxTotal, k);

        % Base cap interior layers (m = 2 .. Nb-1)
        for m = 2:(Nb - 1)
            vid_base(m, k) = addVertex(X_base(m, k), Y_base(m, k), Z_base(m, k));
        end

        % Base cap bottom edge shares aft-most Lower Trailing Edge (m = Nb -> i = NxTotal)
        vid_base(Nb, k) = vid_low(NxTotal, k);
    end

    vertices = vertices(1:nVerts, :);

    % 4. Triangulate Waverider Forebody Upper & Lower (rows 1 .. NxWave)
    facesUpWave = triangulateStructuredGrid(vid_up(1:NxWave, :), false);
    [facesLowWave, cpLowWave] = triangulateStructuredGrid(vid_low(1:NxWave, :), true, surfData.Cp_low(1:NxWave, :));

    % 5. Triangulate Conformal Transition Frustum (rows NxWave .. NxFrustEnd)
    facesFrust = zeros(0, 3);
    if NxFrustEnd > NxWave
        fUpFrust  = triangulateStructuredGrid(vid_up(NxWave:NxFrustEnd, :), false);
        fLowFrust = triangulateStructuredGrid(vid_low(NxWave:NxFrustEnd, :), true);
        facesFrust = [fUpFrust; fLowFrust];
    end

    % 6. Triangulate Conical Rocket Body (rows NxFrustEnd .. NxTotal)
    facesRocket = zeros(0, 3);
    if NxTotal > NxFrustEnd
        fUpRock  = triangulateStructuredGrid(vid_up(NxFrustEnd:NxTotal, :), false);
        fLowRock = triangulateStructuredGrid(vid_low(NxFrustEnd:NxTotal, :), true);
        facesRocket = [fUpRock; fLowRock];
    end

    % 7. Triangulate Aft Base Cap (outward normal: +X)
    facesBase = triangulateStructuredGrid(vid_base, false);

    faces = [facesUpWave; facesLowWave; facesFrust; facesRocket; facesBase];
    facePatchID = [1 * ones(size(facesUpWave, 1), 1); ...
                   2 * ones(size(facesLowWave, 1), 1); ...
                   4 * ones(size(facesFrust, 1), 1); ...
                   5 * ones(size(facesRocket, 1), 1); ...
                   3 * ones(size(facesBase, 1), 1)];
    faceCpDesign = [nan(size(facesUpWave, 1), 1); ...
                    cpLowWave; ...
                    nan(size(facesFrust, 1), 1); ...
                    nan(size(facesRocket, 1), 1); ...
                    nan(size(facesBase, 1), 1)];

    TR = triangulation(faces, vertices);

    meshData = struct();
    meshData.TR = TR;
    meshData.Vertices = vertices;
    meshData.Faces = faces;
    meshData.FacesUpper = facesUpWave;
    meshData.FacesLower = facesLowWave;
    meshData.FacesFrustum = facesFrust;
    meshData.FacesRocket = facesRocket;
    meshData.FacesBase = facesBase;
    meshData.FacePatchID = facePatchID; % 1=Upper, 2=Lower, 3=Base, 4=Frustum, 5=Rocket
    meshData.FaceCpDesign = faceCpDesign;
    meshData.X_base = X_base;
    meshData.Y_base = Y_base;
    meshData.Z_base = Z_base;
end

function [faces, faceVals] = triangulateStructuredGrid(vidGrid, reverseWinding, valGrid)
%TRIANGULATESTRUCTUREDGRID Mirror-symmetric triangulation of a structured vertex-ID grid.
%   Uses diagonal A-C on the port half (c <= floor((Nc-1)/2)) and mirror
%   diagonal B-D on the starboard half (c > floor((Nc-1)/2)) so symmetric
%   waveriders have exact port/starboard aerodynamic and RCS symmetry.
    [Nr, Nc] = size(vidGrid);
    maxF = 2 * (Nr - 1) * (Nc - 1);
    faces = zeros(maxF, 3);
    hasVal = (nargin >= 3) && ~isempty(valGrid);
    if hasVal
        faceVals = zeros(maxF, 1);
    else
        faceVals = [];
    end
    nFaces = 0;
    halfCols = floor((Nc - 1) / 2);

    for r = 1:(Nr - 1)
        for c = 1:(Nc - 1)
            A = vidGrid(r, c);
            B = vidGrid(r + 1, c);
            C = vidGrid(r + 1, c + 1);
            D = vidGrid(r, c + 1);

            if c <= halfCols
                % Port half: split quad along diagonal A-C
                if (A ~= B) && (B ~= C) && (C ~= A)
                    nFaces = nFaces + 1;
                    if reverseWinding
                        faces(nFaces, :) = [A, C, B];
                    else
                        faces(nFaces, :) = [A, B, C];
                    end
                    if hasVal
                        faceVals(nFaces) = (valGrid(r, c) + valGrid(r + 1, c) + valGrid(r + 1, c + 1)) / 3.0;
                    end
                end
                if (A ~= C) && (C ~= D) && (D ~= A)
                    nFaces = nFaces + 1;
                    if reverseWinding
                        faces(nFaces, :) = [A, D, C];
                    else
                        faces(nFaces, :) = [A, C, D];
                    end
                    if hasVal
                        faceVals(nFaces) = (valGrid(r, c) + valGrid(r + 1, c + 1) + valGrid(r, c + 1)) / 3.0;
                    end
                end
            else
                % Starboard half: split quad along mirror diagonal B-D
                if (A ~= B) && (B ~= D) && (D ~= A)
                    nFaces = nFaces + 1;
                    if reverseWinding
                        faces(nFaces, :) = [A, D, B];
                    else
                        faces(nFaces, :) = [A, B, D];
                    end
                    if hasVal
                        faceVals(nFaces) = (valGrid(r, c) + valGrid(r + 1, c) + valGrid(r, c + 1)) / 3.0;
                    end
                end
                if (B ~= C) && (C ~= D) && (D ~= B)
                    nFaces = nFaces + 1;
                    if reverseWinding
                        faces(nFaces, :) = [B, D, C];
                    else
                        faces(nFaces, :) = [B, C, D];
                    end
                    if hasVal
                        faceVals(nFaces) = (valGrid(r + 1, c) + valGrid(r + 1, c + 1) + valGrid(r, c + 1)) / 3.0;
                    end
                end
            end
        end
    end

    faces = faces(1:nFaces, :);
    if hasVal
        faceVals = faceVals(1:nFaces);
    end
end
