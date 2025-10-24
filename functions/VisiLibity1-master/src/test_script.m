% SCRIPT TO COMPARE PERFORMANCE OF visibility_polygon_v2 AND visiPolygon
% This version uses a simple polygon with NO holes.

clear; clc; close all;

%% 1. DEFINE TEST ENVIRONMENT AND OBSERVER
% =========================================================================
% Outer boundary (a non-convex polygon, counter-clockwise)
outer_boundary = [0, 0; 100, -10; 60, 50; 100, 100; 0, 100; -10, 50; 0, 0];

% Observer (viewpoint) position inside the polygon
observer_x = 80;
observer_y = 0;


% Parameters for the functions
epsilon = 1e-9;
snap_distance = 0.01;
threshold = 0.1; % For visiPolygon
dispflag = 0;      % Don't display figures inside the loop

fprintf('Starting performance comparison (Outer Boundary Only)...\n');

%% 2. PREPARE INPUTS FOR EACH FUNCTION
% =========================================================================

% --- For visibility_polygon_v2 (cell array with one polygon) ---
% Note: We remove the last point because the function closes the polygon automatically
environment_v2 = {outer_boundary(1:end-1,:)};

% --- For visiPolygon (simple x and y vectors) ---
xBoundary = outer_boundary(:,1);
yBoundary = outer_boundary(:,2);


%% 3. RUN BENCHMARK FOR visibility_polygon_v2
% =========================================================================
num_iterations = 1000;
time_v2 = 0;

for i = 1:num_iterations
    tic;
    [vis_poly_v2, ~] = visibility_polygon_v2([observer_x observer_y], environment_v2, epsilon, snap_distance);
    time_v2 = time_v2 + toc;
end
avg_time_v2 = time_v2 / num_iterations;


%% 4. RUN BENCHMARK FOR visiPolygon
% =========================================================================
time_visiPolygon = 0;

for i = 1:num_iterations
    tic;
    [xVisi, yVisi, ~, ~] = visiPolygon(xBoundary, yBoundary, observer_x, observer_y, threshold, dispflag);
    time_visiPolygon = time_visiPolygon + toc;
end
avg_time_visiPolygon = time_visiPolygon / num_iterations;


%% 5. DISPLAY RESULTS
% =========================================================================
fprintf('\n--- Performance Results ---\n');
fprintf('Average time for visibility_polygon_v2 (MEX): %.6f seconds\n', avg_time_v2);
fprintf('Average time for visiPolygon (M-file):      %.6f seconds\n', avg_time_visiPolygon);
fprintf('\n');

if avg_time_v2 < avg_time_visiPolygon
    speedup = avg_time_visiPolygon / avg_time_v2;
    fprintf('The MEX function (visibility_polygon_v2) is %.2f times faster.\n', speedup);
else
    speedup = avg_time_v2 / avg_time_visiPolygon;
    fprintf('The M-file (visiPolygon) is %.2f times faster.\n', speedup);
end


%% 6. VISUALIZE AND VALIDATE THE OUTPUTS
% =========================================================================
% Recalculate one last time to get the polygon data for plotting
[vis_poly_v2, growing_verts] = visibility_polygon_v2([observer_x observer_y], environment_v2, epsilon, snap_distance);
[xVisi, yVisi, xEff, yEff] = visiPolygon(xBoundary, yBoundary, observer_x, observer_y, threshold, 0);

figure('Name', 'Visibility Polygon Comparison (Outer Boundary Only)', 'NumberTitle', 'off');
hold on;
grid on;
axis equal;

% Plot Environment
plot(outer_boundary(:,1), outer_boundary(:,2), 'k-', 'LineWidth', 2, 'DisplayName', 'Environment');

% Plot Observer
plot(observer_x, observer_y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Observer');

% Plot results from visibility_polygon_v2
plot([vis_poly_v2(:,1); vis_poly_v2(1,1)], [vis_poly_v2(:,2); vis_poly_v2(1,2)], 'b-', 'LineWidth', 2, 'DisplayName', 'visibility_polygon_v2');
plot(growing_verts(:,1), growing_verts(:,2), 'g*', 'MarkerSize', 8, 'DisplayName', 'Growing Vertices (v2)');

% Plot results from visiPolygon
plot(xVisi, yVisi, 'm--', 'LineWidth', 2, 'DisplayName', 'visiPolygon');
plot(xEff, yEff, 'c^', 'MarkerSize', 8, 'DisplayName', 'Effective Corners (visiPolygon)');

title('Comparison of Visibility Polygon Functions');
legend('show', 'Location', 'best');
xlabel('X Coordinate');
ylabel('Y Coordinate');
hold off;