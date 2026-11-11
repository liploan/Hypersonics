function report = verifyWatertight(meshData, cfg)
%VERIFYWATERTIGHT Rigorous topological & volumetric verification of a 3D mesh.
%   report = VERIFYWATERTIGHT(meshData, cfg) checks boundary edges, edge
%   manifoldness, half-edge orientation consistency, Euler characteristic,
%   face areas, enclosed volume, planform area, and volumetric efficiency.

    TR = meshData.TR;
    V = TR.Points;
    F = TR.ConnectivityList;

    numVertices = size(V, 1);
    numFaces = size(F, 1);

    % 1. Free boundary edges via MATLAB triangulation
    FB = freeBoundary(TR);
    numBoundaryEdges = size(FB, 1);

    % 2. Directed half-edge and undirected edge manifoldness audit
    he1 = F(:, [1, 2]);
    he2 = F(:, [2, 3]);
    he3 = F(:, [3, 1]);
    halfEdges = [he1; he2; he3]; % [3*numFaces x 2]

    sortedEdges = sort(halfEdges, 2);
    [uniqueEdges, ~, ic] = unique(sortedEdges, 'rows');
    numEdges = size(uniqueEdges, 1);
    edgeCounts = accumarray(ic, 1);
    numNonManifoldEdges = sum(edgeCounts ~= 2);

    % Check half-edge cancellation: +1 if he(:,1) < he(:,2), -1 if he(:,1) > he(:,2)
    heSign = sign(halfEdges(:, 2) - halfEdges(:, 1));
    netEdgeOrientation = accumarray(ic, heSign);
    numConflictingHalfEdges = sum(netEdgeOrientation ~= 0);

    % 3. Euler characteristic V - E + F (equals 2 for a closed genus-0 polyhedron)
    eulerCharacteristic = numVertices - numEdges + numFaces;

    % 4. Face areas, normals, and signed enclosed volume via divergence theorem
    v1 = V(F(:, 1), :);
    v2 = V(F(:, 2), :);
    v3 = V(F(:, 3), :);

    cross23 = cross(v2 - v1, v3 - v1, 2);
    faceAreas = 0.5 .* sqrt(sum(cross23.^2, 2));
    numDegenerateFaces = sum(faceAreas < 1e-12);

    tetSignedVols = (1.0 / 6.0) .* sum(v1 .* cross(v2, v3, 2), 2);
    enclosedVolume = sum(tetSignedVols);
    tetCentroids = 0.25 .* (v1 + v2 + v3);
    if abs(enclosedVolume) > 1e-12
        volumeCentroid = sum(tetCentroids .* tetSignedVols, 1) ./ enclosedVolume;
    else
        volumeCentroid = mean(V, 1);
    end

    % 5. Surface patch areas (Upper, Lower, Frustum, Rocket, Base) and Planform Area
    isUpper   = (meshData.FacePatchID == 1);
    isLower   = (meshData.FacePatchID == 2);
    isBase    = (meshData.FacePatchID == 3);
    isFrustum = (meshData.FacePatchID == 4);
    isRocket  = (meshData.FacePatchID == 5);

    upperArea   = sum(faceAreas(isUpper));
    lowerArea   = sum(faceAreas(isLower));
    frustumArea = sum(faceAreas(isFrustum));
    rocketArea  = sum(faceAreas(isRocket));
    baseArea    = sum(faceAreas(isBase));
    totalArea   = sum(faceAreas);

    % Planform area = projection of Upper forebody surface onto horizontal (X-Y) plane
    planformArea = 0.5 * sum(abs(cross23(isUpper, 3)));
    if planformArea < 1e-9
        planformArea = 0.5 * sum(abs(cross23(isLower, 3)));
    end

    % Volumetric efficiency: eta_V = V^(2/3) / S_planform
    volEfficiency = (max(enclosedVolume, 0.0)^(2.0 / 3.0)) / max(planformArea, 1e-12);

    isWatertight = (numBoundaryEdges == 0) && ...
                   (numNonManifoldEdges == 0) && ...
                   (numConflictingHalfEdges == 0) && ...
                   (numDegenerateFaces == 0) && ...
                   (enclosedVolume > 0);

    report = struct();
    report.IsWatertight = isWatertight;
    report.NumVertices = numVertices;
    report.NumEdges = numEdges;
    report.NumFaces = numFaces;
    report.NumBoundaryEdges = numBoundaryEdges;
    report.NumNonManifoldEdges = numNonManifoldEdges;
    report.NumConflictingHalfEdges = numConflictingHalfEdges;
    report.NumDegenerateFaces = numDegenerateFaces;
    report.EulerCharacteristic = eulerCharacteristic;
    report.EnclosedVolume = enclosedVolume;
    report.VolumeCentroid = volumeCentroid;
    report.PlanformArea = planformArea;
    report.WettedArea = upperArea + lowerArea + frustumArea + rocketArea;
    report.UpperSurfaceArea = upperArea;
    report.LowerSurfaceArea = lowerArea;
    report.FrustumArea = frustumArea;
    report.RocketStageArea = rocketArea;
    report.BaseArea = baseArea;
    report.TotalSurfaceArea = totalArea;
    report.VolumetricEfficiency = volEfficiency;
    if nargin >= 2 && ~isempty(cfg)
        report.AspectRatio = (2.0 * cfg.HalfSpan)^2 / max(planformArea, 1e-12);
    end
end
