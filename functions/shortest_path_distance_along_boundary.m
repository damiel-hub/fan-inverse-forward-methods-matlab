function xysBoundary = shortest_path_distance_along_boundary(xMesh_crop, yMesh_crop, zMesh_crop)

    % Calculate diagonal length of the mesh grid
    diagonal_length = sqrt((xMesh_crop(1,1) - xMesh_crop(1,end))^2 + (yMesh_crop(1,1) - yMesh_crop(end,1))^2);
    
    % Find the apex point (highest point) in the mesh
    [~, iApex] = max(zMesh_crop(:));
    xApex = xMesh_crop(iApex);
    yApex = yMesh_crop(iApex);
    
    % Create a wall mesh with boundary values set to a high value
    wallMesh = zeros(size(zMesh_crop));
    wallMesh(isnan(zMesh_crop)) = diagonal_length * 10;
    
    % Set the apex height for shortest path calculation
    zApex_s = diagonal_length * 10;
    
    % Calculate visibility polygon and apex information
    [~,~,xyzkApexAll,xyzVisPolygon,~,~] = FanTopo(xMesh_crop, yMesh_crop, wallMesh, xApex, yApex, zApex_s, ...
                                                           'tanAlphaM', 1, 'saveVisPolygon', 1);

    % Extract 3D outline from the visibility polygon and apex information
    poly_xyz_all = get_3D_outline_from_FanTopo_cone(xyzkApexAll, xyzVisPolygon, 1);
    
    % Initialize the boundary array
    xysBoundary = [];
    
    % Calculate the shortest path distance along the boundary
    for i = 1:length(poly_xyz_all)
        sOutline = zApex_s - poly_xyz_all{i}(:,3);
        xysBoundary = [xysBoundary; poly_xyz_all{i}(:,1:2), sOutline;nan nan nan];
    end
    
    
end
