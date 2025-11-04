function [inMesh, onedgeMesh] = inpolygon_optimized(xMesh, yMesh, xv, yv)
% INPOLYGON_OPTIMIZED A highly efficient, vectorized scan-line algorithm.
%
% DESCRIPTION:
%   This function is an optimized version of the corrected scan-line
%   algorithm. It vectorizes the scan-line fill using a cumsum technique
%   and optimizes the on-edge detection by operating on sub-grids.

    %% --- 1. Initialization ---
    grid_size = size(xMesh);
    value_x = xMesh(1, :);
    value_y = yMesh(:, 1)';
    
    diffCount = zeros(grid_size); % Temporary matrix for the cumsum trick
    onedgeMesh = false(grid_size);
    dist_tolerance = 1e-9;

    %% --- 2. Process Polygons ---
    nan_indices = [0; find(isnan(xv(:))); length(xv(:)) + 1];
    
    for p = 1:length(nan_indices) - 1
        start_idx = nan_indices(p) + 1;
        end_idx = nan_indices(p+1) - 1;
        
        polygon_x = xv(start_idx:end_idx);
        polygon_y = yv(start_idx:end_idx);
        
        if isempty(polygon_x); continue; end
        if polygon_x(1) ~= polygon_x(end) || polygon_y(1) ~= polygon_y(end)
            polygon_x(end+1) = polygon_x(1);
            polygon_y(end+1) = polygon_y(1);
        end

        % --- Part A: Optimized On-Edge Detection ---
        for edge_i = 1:length(polygon_x) - 1
            ax = polygon_x(edge_i);     ay = polygon_y(edge_i);
            bx = polygon_x(edge_i+1);   by = polygon_y(edge_i+1);
            
            edge_len_sq = (bx - ax)^2 + (by - ay)^2;
            if edge_len_sq < dist_tolerance^2; continue; end
            
            % OPTIMIZATION: Find sub-grid indices instead of a full logical mask
            minX = min(ax, bx); maxX = max(ax, bx);
            minY = min(ay, by); maxY = max(ay, by);
            
            col_indices = find(value_x >= minX - dist_tolerance & value_x <= maxX + dist_tolerance);
            row_indices = find(value_y >= minY - dist_tolerance & value_y <= maxY + dist_tolerance);
            
            if isempty(col_indices) || isempty(row_indices); continue; end

            % Perform calculations only on the small sub-grid
            [sub_xMesh, sub_yMesh] = meshgrid(value_x(col_indices), value_y(row_indices));
            
            cross_product = (sub_yMesh - ay) .* (bx - ax) - (sub_xMesh - ax) .* (by - ay);
            distance_to_line = abs(cross_product) / sqrt(edge_len_sq);

            % Map results back to the full onedgeMesh
            on_this_edge_sub = distance_to_line < dist_tolerance;
            onedgeMesh(row_indices, col_indices) = onedgeMesh(row_indices, col_indices) | on_this_edge_sub;
        end
        
        % --- Part B: Vectorized Scan-Line Fill ---
        value_x_index = 1:length(value_x);
        value_y_index = 1:length(value_y);
        polygon_x_index = interp1(value_x, value_x_index, polygon_x, 'linear', 'extrap');
        polygon_y_index = interp1(value_y, value_y_index, polygon_y, 'linear', 'extrap');

        for edge_i = 1:length(polygon_x_index) - 1
            x1 = polygon_x_index(edge_i);       y1 = polygon_y_index(edge_i);
            x2 = polygon_x_index(edge_i + 1);   y2 = polygon_y_index(edge_i + 1);
            
            if abs(y1 - y2) < 1e-9; continue; end
            
            y_upper = max(y1, y2); y_lower = min(y1, y2);
            y_min_index = ceil(y_lower);
            y_max_index = floor(y_upper);
            if abs(y_upper - y_max_index) < 1e-9
                y_max_index = y_max_index - 1;
            end
            if y_min_index > y_max_index; continue; end
            
            y_seg = y_min_index:y_max_index;
            x_intersections = floor(interp1([y1 y2], [x1 x2], y_seg, 'linear'));
            
            % OPTIMIZATION: Vectorized update to the diffCount matrix
            valid_mask = (y_seg >= 1) & (y_seg <= grid_size(1));
            y_seg = y_seg(valid_mask);
            x_intersections = x_intersections(valid_mask);

            start_cols = ones(size(y_seg));
            end_cols = x_intersections + 1;
            
            valid_starts = x_intersections >= 1;
            valid_ends = (end_cols >= 1) & (end_cols <= grid_size(2));

            start_indices = sub2ind(grid_size, y_seg(valid_starts), start_cols(valid_starts));
            end_indices = sub2ind(grid_size, y_seg(valid_ends), end_cols(valid_ends));

            diffCount(start_indices) = diffCount(start_indices) + 1;
            diffCount(end_indices) = diffCount(end_indices) - 1;
        end
    end
    
    %% --- 3. Finalize Outputs ---
    crossingCount = cumsum(diffCount, 2);
    totalInclusionMask = mod(crossingCount, 2) == 1;
    inMesh = totalInclusionMask & ~onedgeMesh;
end