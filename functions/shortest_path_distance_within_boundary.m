function sMap = shortest_path_distance_within_boundary(xMesh_crop, yMesh_crop, zMesh_crop, pltFlag, xApex, yApex)
% shortest_path_distance_within_boundary
% Now accepts optional xApex and yApex arguments.
% If they are not provided, the function finds the highest point automatically.

    % Calculate the diagonal length of the mesh grid
    diagonal_length = sqrt((xMesh_crop(1,1) - xMesh_crop(1,end))^2 + (yMesh_crop(1,1) - yMesh_crop(end,1))^2);
    
    % --- MODIFICATION START ---
    % Check if xApex and yApex are provided (total inputs should be 6)
    % If fewer than 6 inputs, we calculate them automatically (Backward Compatibility)
    if nargin < 6
        [~, iApex] = max(zMesh_crop(:));
        xApex = xMesh_crop(iApex);
        yApex = yMesh_crop(iApex);
    end
    % --- MODIFICATION END ---
    
    % Create a wall mesh with boundary values set to a high value
    wallMesh = zeros(size(zMesh_crop));
    wallMesh(isnan(zMesh_crop)) = diagonal_length * 10;
    
    % Set the apex height for the shortest path calculation
    zApex_s = diagonal_length * 10;
    
    % writeGeoTiff(wallMesh, 'zWall.tif', 3826, min(xMesh_crop(:)), max(xMesh_crop(:)), min(yMesh_crop(:)), max(yMesh_crop(:)), 'north', 'west')
    
    % Compute the shortest path topography
    % (Uses either the passed Apex or the calculated Apex from above)
    [sTopo, ~, ~, ~, ~, ~] = FanTopo(xMesh_crop, yMesh_crop, wallMesh, xApex, yApex, zApex_s, 'tanAlphaM', 1);
    
    % Calculate the shortest path distance map
    sMap = zApex_s - sTopo;
    
    % Plot the shortest path distance map if pltFlag is true
    if pltFlag
        figure
        pcolor(xMesh_crop, yMesh_crop, wallMesh)
        shading flat
        axis equal
        
        figure
        imagesc(xMesh_crop(1,:), yMesh_crop(:,1), sMap)
        hold on
        plot(xApex, yApex, 'r.', 'MarkerSize', 6)
        contour(xMesh_crop, yMesh_crop, sMap, 0:100:max(sMap(:)), 'k')
        axis xy
        axis equal
        axis tight
        c = colorbar;
        ylabel(c,'Shortest Path Distance (m)')
        title('Shortest Path Distance Map')
        xlabel('Easting (m)')
        ylabel('Northing (m)')
    end
end