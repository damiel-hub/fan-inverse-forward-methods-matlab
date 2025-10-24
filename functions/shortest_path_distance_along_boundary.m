function xyzsBoundary = shortest_path_distance_along_boundary(inputGeoTiffPath, shapeFilePath)
    %% 1. SETUP AND DATA LOADING
    epsilon = 1e-9;

    % Read georeferenced raster and shapefile data
    [zMesh, R] = readgeoraster(inputGeoTiffPath);
    [xMesh, yMesh] = worldGrid(R);
    fan_boundary = shaperead(shapeFilePath);
    fan_boundary_x = fan_boundary.X(:); % Ensure column vector
    fan_boundary_y = fan_boundary.Y(:);

    %% 2. PARSE AND VALIDATE ENVIRONMENT POLYGONS
    environment = {};
    nan_indices = [0; find(isnan(fan_boundary_x))];

    for k = 1:length(nan_indices)
        if k == length(nan_indices)
            start_idx = nan_indices(k) + 1;
            end_idx = length(fan_boundary_x);
        else
            start_idx = nan_indices(k) + 1;
            end_idx = nan_indices(k+1) - 1;
        end
        
        if end_idx >= start_idx
            poly_segment = [fan_boundary_x(start_idx:end_idx), fan_boundary_y(start_idx:end_idx)];
            environment{end+1} = poly_segment;
        end
    end
    
    % Reorder to ensure the largest polygon (assumed outer boundary) is first
    if numel(environment) > 1
        poly_areas = cellfun(@(p) polyarea(p(:,1), p(:,2)), environment);
        [~, max_area_idx] = max(poly_areas);
        environment = [environment(max_area_idx), environment(1:max_area_idx-1), environment(max_area_idx+1:end)];
    end

    % --- **CRITICAL: ENFORCE CORRECT WINDING ORDER** ---
    % Outer boundary (first cell) must be Counter-Clockwise (CCW)
    if ispolycw(environment{1}(:,1), environment{1}(:,2))
        fprintf('Info: Outer boundary was clockwise, flipping to CCW.\n');
        environment{1} = flipud(environment{1});
    end

    % Inner holes (all other cells) must be Clockwise (CW)
    for k = 2:numel(environment)
        if ~ispolycw(environment{k}(:,1), environment{k}(:,2))
            fprintf('Info: Hole %d was CCW, flipping to CW.\n', k-1);
            environment{k} = flipud(environment{k});
        end
    end
    
    %% 3. PREPARE FOR CALCULATION
    dxdy = xMesh(1,2) - xMesh(1,1);
    snap_distance = dxdy;
    fan_boundary_xy = environment{1};
    
    % Resample the boundary
    fan_boundary_s = [0; cumsum(sqrt(sum(diff(fan_boundary_xy).^2,2)))];
    fan_boundary_s_resample = (0:dxdy:fan_boundary_s(end))';
    fan_boundary_x_resample = interp1(fan_boundary_s, fan_boundary_xy(:,1), fan_boundary_s_resample);
    fan_boundary_y_resample = interp1(fan_boundary_s, fan_boundary_xy(:,2), fan_boundary_s_resample);
    
    % Find the Apex
    mask = inpolygon(xMesh, yMesh, fan_boundary_xy(:,1), fan_boundary_xy(:,2));
    zMesh_original = zMesh;
    zMesh(~mask) = nan;
    [~,index_max] = max(zMesh(:));
    [row, col] = ind2sub(size(zMesh), index_max);
    xApex = xMesh(row, col);
    yApex = yMesh(row, col);
    
    %% 4. ROBUST SHORTEST PATH CALCULATION
    num_points = length(fan_boundary_y_resample);
    distance_all = NaN(num_points, 1);


    % Ensure the Apex itself is valid
    if ~in_environment([xApex yApex], environment, epsilon)
        error('The calculated Apex at (%.2f, %.2f) is outside the fan boundary.', xApex, yApex);
    end

    parfor i = 1:num_points
        start_point = [xApex, yApex];
        end_point_original = [fan_boundary_x_resample(i), fan_boundary_y_resample(i)];
        
        % Call the MEX function to get path vertices
        path_vertices = shortest_path(start_point, end_point_original, environment, epsilon, snap_distance);
        
        % Calculate total path length from the vertices
        if size(path_vertices, 1) > 1
            distance_all(i) = sum(sqrt(sum(diff(path_vertices, 1, 1).^2, 2)));
        else
            distance_all(i) = 0;
        end
    end
    
    %% 5. FINALIZE OUTPUT
    fan_boundary_z_resample = interp2(xMesh, yMesh, zMesh_original, fan_boundary_x_resample, fan_boundary_y_resample);
    xyzsBoundary = [fan_boundary_x_resample(:) fan_boundary_y_resample(:) fan_boundary_z_resample(:) distance_all(:)];
end