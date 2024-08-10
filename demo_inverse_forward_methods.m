path(path,'functions')
path(path,'functions\inpoly-master')

%% Set data path name
DEM_resolution = 5; % In the demo-code, we provide 5, 10 and 20 meter resolution tif file.

% Estimated Running Time for Different DEM Resolutions:
% 5-meter: 349.4 seconds
% 10-meter: 73.2 seconds
% 20-meter: 28.5 seconds
% 
% System Specifications:
% CPU: AMD Ryzen 9 5900X
% GPU: NVIDIA GeForce RTX 3080
% MATLAB Version: R2023b

tic
%% DEM_resolution
topo_pre_event = ['data\topo_PuTunPuNas_min_before2014_riverbed2014_' num2str(DEM_resolution) 'm.tif'];
topo_post_event = ['data\topo_PuTunPuNas_2014_' num2str(DEM_resolution) 'm.tif'];
shape_fan_boundary = 'data\shape\PT2014.shp';

%% Inverse method: Extract the fan profile from the field data
% Step 1: Calculate shortest path distance within or along boundary
[xMesh, yMesh, zMesh_post] = readGeoTiff(topo_post_event);
[xMesh_crop, yMesh_crop, zMesh_post_crop] = clipGeoTiff(topo_post_event, shape_fan_boundary);

% Within boundary
sMap = shortest_path_distance_within_boundary(xMesh_crop, yMesh_crop, zMesh_post_crop, 1);

% Along boundary
xysBoundary = shortest_path_distance_along_boundary(xMesh_crop, yMesh_crop, zMesh_post_crop);
zBoundary = interp2(xMesh, yMesh, zMesh_post, xysBoundary(:,1), xysBoundary(:,2));

%% Step 2: fitting with quadratic elevation-distance relationship
bin_size = 100;
ds = 5;
outlength = 500;

% Within boundary
figure
scatter(sMap(:), zMesh_post_crop(:), 'k.');
hold on
fitting_s_z_within_boundary = process_s_z_relationship(sMap, zMesh_post_crop, bin_size, ds, outlength, 1);
xlabel('Shortest in polygon distance, s (m)')
ylabel('Elevation, z (m)')

% Along boundary
figure
plot(xysBoundary(:,3), zBoundary, 'k-');
hold on
fitting_s_z_along_boundary = process_s_z_relationship(xysBoundary(:,3), zBoundary, bin_size, ds, outlength, 1, 'medianFilter', 0);
xlabel('Shortest distance on boundary of polygon, s (m)')
ylabel('Elevation, z (m)')
%% Forward method: Reconstruct the debris and alluvial fan topography using the elevation-distance profile

% Read initial topography
[xMesh, yMesh, zMesh_pre] = readGeoTiff(topo_pre_event);
[~, ~, zMesh_pre_crop] = clipGeoTiff(topo_pre_event, shape_fan_boundary);

% Calculate fan volume
fanSimVolume = sum(zMesh_post_crop - zMesh_pre_crop, 'all', 'omitnan') * (xMesh_crop(1,2) - xMesh_crop(1,1)).^2;
fprintf('The expected simulated fan volume is %.2f cubic meters (L^3) within the defined boundary.\n', fanSimVolume);

% Use the highest point in the boundary as apex location
[zApex, iApex] = max(zMesh_post_crop(:));
xApex = xMesh_crop(iApex);
yApex = yMesh_crop(iApex);

guessHeightAboveGround_bottom = 1;
guessHeightAboveGround_top = 10;

% Reconstruct the topography
[zTopo_sim, heightAG_Volume_All] = reconstruct_fan_surface(xMesh, yMesh, zMesh_pre, xApex, yApex, fanSimVolume, guessHeightAboveGround_top, guessHeightAboveGround_bottom, fitting_s_z_within_boundary,"fanBoundarySHP",shape_fan_boundary);

% Plot the result (calculate volume within simulation area)
plotFanTopoResults(xMesh, yMesh, zTopo_sim, zMesh_pre, xApex, yApex)
% Plot the result (calculate volume within within given boundary)
plotFanTopoResults(xMesh, yMesh, zTopo_sim, zMesh_pre, xApex, yApex, shape_fan_boundary)
% Plot the field elevation difference before and after event
plotFanTopoResults(xMesh, yMesh, zMesh_post, zMesh_pre, xApex, yApex, shape_fan_boundary)

toc
