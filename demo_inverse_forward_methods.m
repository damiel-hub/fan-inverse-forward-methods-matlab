% Add necessary paths for function dependencies
addpath('functions', 'functions/inpoly-master');

%% Set DEM resolution and file paths
DEM_resolution = 10; % Available resolutions: 10, 20 meters

% Estimated Running Time for Different DEM Resolutions:
% System 1: AMD Ryzen 9 5900X, NVIDIA GeForce RTX 3080, MATLAB R2023b
% 10-meter: 73.2 seconds, 20-meter: 28.5 seconds
%
% System 2: Intel Core i7-7700HQ, NVIDIA GeForce GTX 1050, MATLAB R2022a
% 10-meter: 176.1 seconds, 20-meter: 90.5 seconds

tic

% File paths based on resolution
topo_pre_event = sprintf('data/topo_PuTunPuNas_min_before2014_riverbed2014_%dm.tif', DEM_resolution);
topo_post_event = sprintf('data/topo_PuTunPuNas_2014_%dm.tif', DEM_resolution);
shape_fan_boundary = 'data/shape/PT2014.shp';

%% Inverse Method: Extract the fan profile from the field data

% Step 1:
% Load and process DEM data
[xMesh, yMesh, zMesh_post] = readGeoTiff(topo_post_event);
[xMesh_crop, yMesh_crop, zMesh_post_crop] = clipGeoTiff(topo_post_event, shape_fan_boundary);

% Calculate shortest path distance within and along boundary
sMap = shortest_path_distance_within_boundary(xMesh_crop, yMesh_crop, zMesh_post_crop, 1);
xysBoundary = shortest_path_distance_along_boundary(xMesh_crop, yMesh_crop, zMesh_post_crop);
zBoundary = interp2(xMesh, yMesh, zMesh_post, xysBoundary(:,1), xysBoundary(:,2));

% Step 2:
% Fit quadratic elevation-distance relationship
bin_size = 100;
ds = 5;
outlength = 500;

% Plot and fit within boundary
figure
scatter(sMap(:), zMesh_post_crop(:), 'k.')
hold on
fitting_s_z_within_boundary = process_s_z_relationship(sMap, zMesh_post_crop, bin_size, ds, outlength, 1);
xlabel('Shortest path distance to all data points, s (m)')
ylabel('Elevation, z (m)')

% Plot and fit along boundary
figure
plot(xysBoundary(:,3), zBoundary, 'k-')
hold on
fitting_s_z_along_boundary = process_s_z_relationship(xysBoundary(:,3), zBoundary, bin_size, ds, outlength, 1, 'medianFilter', 0);
xlabel('Shortest path distance to boundary points, s (m)')
ylabel('Elevation, z (m)')

%% Forward Method: Reconstruct the debris and alluvial fan topography

% Step 3
% Load initial topography
[xMesh, yMesh, zMesh_pre] = readGeoTiff(topo_pre_event);
[~, ~, zMesh_pre_crop] = clipGeoTiff(topo_pre_event, shape_fan_boundary);

% Calculate fan volume
fanSimVolume = sum(zMesh_post_crop - zMesh_pre_crop, 'all', 'omitnan') * (xMesh_crop(1,2) - xMesh_crop(1,1))^2;
fprintf('The expected simulated fan volume is %.2f cubic meters (L^3) within the defined boundary.\n', fanSimVolume);

% Determine the apex location
[zApex, iApex] = max(zMesh_post_crop(:));
xApex = xMesh_crop(iApex);
yApex = yMesh_crop(iApex);

% Set initial guesses for height above ground
guessHeightAboveGround_bottom = 1;
guessHeightAboveGround_top = 10;

% Reconstruct topography
[zTopo_sim, heightAG_Volume_All] = reconstruct_fan_surface(xMesh, yMesh, zMesh_pre, xApex, yApex, fanSimVolume, guessHeightAboveGround_top, guessHeightAboveGround_bottom, fitting_s_z_within_boundary, "fanBoundarySHP", shape_fan_boundary);

% Plot results
plotFanTopoResults(xMesh, yMesh, zTopo_sim, zMesh_pre, xApex, yApex); % Volume within simulation area
plotFanTopoResults(xMesh, yMesh, zTopo_sim, zMesh_pre, xApex, yApex, shape_fan_boundary); % Volume within given boundary
plotFanTopoResults(xMesh, yMesh, zMesh_post, zMesh_pre, xApex, yApex, shape_fan_boundary); % Elevation difference before and after event

toc
