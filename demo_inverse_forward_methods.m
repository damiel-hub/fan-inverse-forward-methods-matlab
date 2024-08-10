path(path,'functions')
path(path,'functions\inpoly-master')
%% Set data path name
topo_pre_event = 'data\topo_PuTunPuNas_min_before2014_riverbed2014.tif';
topo_post_event = 'data\topo_PuTunPuNas_2014.tif';
shape_fan_boundary = 'data\shape\PT2014.shp';

%% Inverse method: Extract the fan profile from the field data
% Step 1: Calculate shortest path distance within or along boundary
[xMesh, yMesh, zMesh_post] = readGeoTiff(topo_post_event);
[xMesh_crop, yMesh_crop, zMesh_post_crop] = clipGeoTiff(topo_post_event, shape_fan_boundary);

% Within boundary
sMap = shortest_path_distance_within_boundary(xMesh_crop, yMesh_crop, zMesh_post_crop, 0);

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

% Along boundary
figure
plot(xysBoundary(:,3), zBoundary, 'k.-');
hold on
fitting_s_z_along_boundary = process_s_z_relationship(xysBoundary(:,3), zBoundary, bin_size, ds, outlength, 1);

%% Forward method: Reconstruct the debris and alluvial fan topography using the elevation-distance profile

% Read initial topography
[xMesh, yMesh, zMesh_pre] = readGeoTiff(topo_pre_event);
[~, ~, zMesh_pre_crop] = clipGeoTiff(topo_pre_event, shape_fan_boundary);

% Calculate fan volume
fanSimVolume = sum(zMesh_post_crop - zMesh_pre_crop, 'all', 'omitnan') * (xMesh_crop(1,2) - xMesh_crop(1,1)).^2;

% Use the highest point in the boundary as apex location
[zApex, iApex] = max(zMesh_post_crop(:));
xApex = xMesh_crop(iApex);
yApex = yMesh_crop(iApex);

guessHeightAboveGround_top = 10;
guessHeightAboveGround_bottom = 1;

% Reconstruct the topography
[zTopo_sim, heightAG_Volume_All] = reconstruct_fan_surface(xMesh, yMesh, zMesh_pre, xApex, yApex, fanSimVolume, guessHeightAboveGround_top, guessHeightAboveGround_bottom, fitting_s_z_within_boundary,"fanBoundarySHP",shape_fan_boundary);
plotFanTopoResults(xMesh, yMesh, zTopo_sim, zMesh_pre, xApex, yApex)
plotFanTopoResults(xMesh, yMesh, zTopo_sim, zMesh_pre, xApex, yApex, shape_fan_boundary)

