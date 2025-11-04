function [inMesh, onedgeMesh] = inpolygon_optimized_fixed(xMesh, yMesh, xv, yv)
% INPOLYGON_OPTIMIZED A highly efficient, vectorized scan-line algorithm.
%
% DESCRIPTION:
%   This function is an optimized version of the corrected scan-line
%   algorithm. It vectorizes the scan-line fill using a cumsum technique
%   and optimizes the on-edge detection by operating on sub-grids.
%
%   FIXED: 
%   1. Replaced meshgrid() call in loop with matrix slicing.
%   2. Corrected on-edge detection to use dot-product (segment vs. line).
%   3. Corrected scan-line fill logic to increment count at intersections.
%   4. Replaced interp1() call in loop with manual interpolation.
    
    %% --- 1. Initialization ---
    grid_size = size(xMesh);
    value_x = xMesh(1, :);
    value_y = yMesh(:, 1)';
    
    diffCount = zeros(grid_size); % Temporary matrix for the cumsum trick
    onedgeMesh = false(grid_size);
    dist_tolerance = 1e-9;
    dist_tolerance_sq = dist_tolerance^2; % --- FIX: Pre-calculate squared tolerance

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
            if edge_len_sq < dist_tolerance_sq; continue; end
            
            % OPTIMIZATION: Find sub-grid indices
            minX = min(ax, bx); maxX = max(ax, bx);
            minY = min(ay, by); maxY = max(ay, by);
            
            col_indices = find(value_x >= minX - dist_tolerance & value_x <= maxX + dist_tolerance);
            row_indices = find(value_y >= minY - dist_tolerance & value_y <= maxY + dist_tolerance);
            
            if isempty(col_indices) || isempty(row_indices); continue; end
            
            % --- FIX 1: Use matrix slicing, NOT meshgrid() ---
            sub_xMesh = xMesh(row_indices, col_indices);
            sub_yMesh = yMesh(row_indices, col_indices);
            
            % --- FIX 2: Correct on-segment logic (dot product) ---
            cross_product = (sub_yMesh - ay) .* (bx - ax) - (sub_xMesh - ax) .* (by - ay);
            distance_sq = cross_product.^2 / edge_len_sq;

            dot_product = (sub_xMesh - ax) .* (bx - ax) + (sub_yMesh - ay) .* (by - ay);
            on_segment = (dot_product >= -dist_tolerance_sq) & (dot_product <= edge_len_sq + dist_tolerance_sq);

            on_this_edge_sub = (distance_sq < dist_tolerance_sq) & on_segment;
            
            % Map results back to the full onedgeMesh
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
            
            if abs(y1 - y2) < 1e-9; continue; end % Skip horizontal edges
            
            y_upper = max(y1, y2); y_lower = min(y1, y2);
            
            % Standard scan-line: process rows from ceil(lower) to floor(upper)
            y_min_index = ceil(y_lower);
            y_max_index = floor(y_upper);

            % Standard scan-line: don't process top vertex pixel
            if abs(y_upper - y_max_index) < 1e-9
                y_max_index = y_max_index - 1;
            end
            
            if y_min_index > y_max_index; continue; end
            
            y_seg = (y_min_index:y_max_index)'; % Must be column vector
            
            % --- FIX 4: Replace interp1 with manual calculation ---
            inv_slope = (x2 - x1) / (y2 - y1);
            x_intersections = x1 + (y_seg - y1) .* inv_slope;
            
            % --- FIX 3: Correct cumsum update logic ---
            % We increment the count at the pixel column just AFTER the intersection
            x_cols_to_update = floor(x_intersections) + 1;

            % Find valid indices that are inside the grid
            valid_mask = (y_seg >= 1) & (y_seg <= grid_size(1)) & ...
                         (x_cols_to_update >= 1) & (x_cols_to_update <= grid_size(2));
            
            if ~any(valid_mask); continue; end

            y_seg_valid = y_seg(valid_mask);
            x_cols_valid = x_cols_to_update(valid_mask);
            
            % Convert to linear indices for fast matrix update
            indices = sub2ind(grid_size, y_seg_valid, x_cols_valid);
            
            % Add 1 to the diffCount at these locations
            diffCount(indices) = diffCount(indices) + 1;
        end
    end
    
    %% --- 3. Finalize Outputs ---
    crossingCount = cumsum(diffCount, 2);
    totalInclusionMask = mod(crossingCount, 2) == 1;
    
    % 'in' = (Inside based on fill) AND (Not on the edge)
    inMesh = totalInclusionMask & ~onedgeMesh;
end