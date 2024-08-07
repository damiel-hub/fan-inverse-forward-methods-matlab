path(path,'functions')
path(path,'functions\inpoly-master')

%% Inverse method: Extract the fan profile from the field data
% Step 1: Calculate shortest path distance within or along boundary
[xMesh_post, yMesh_post, zMesh_post] = readGeoTiff('data\topo_fan.tif');
[xMesh_post_crop, yMesh_post_crop, zMesh_post_crop] = clipGeoTiff('data\topo_fan.tif', 'data\shape\fan_extent.shp');

% Within boundary
sMap = shortest_path_distance_within_boundary(xMesh_post_crop, yMesh_post_crop, zMesh_post_crop, 0);

% Along boundary
xysBoundary = shortest_path_distance_along_boundary(xMesh_post_crop, yMesh_post_crop, zMesh_post_crop);
zBoundary = interp2(xMesh_post, yMesh_post, zMesh_post, xysBoundary(:,1), xysBoundary(:,2));

%% Step 2: fitting with quadratic elevation-distance relationship
bin_size = 10;
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
[xMesh_pre, yMesh_pre, zMesh_pre] = readGeoTiff('data\topo_initial.tif');

% Use the highest point in the boundary as apex location
[zApex, iApex] = max(zMesh_post_crop(:));
xApex = xMesh_post_crop(iApex);
yApex = yMesh_post_crop(iApex);

% Reconstruct the topography
[zTopo, ~, ~, ~, ~, ~] = FanTopo_slope_bd(xMesh_pre, yMesh_pre, zMesh_pre, xApex, yApex, zApex, 'caseName', 'myProfile','dz_interpM', {fitting_s_z_within_boundary});

zMap = zTopo;
zMap(isnan(zMap)) = zMesh_pre(isnan(zMap));

% Plot the result
figure
pcolor(xMesh_pre, yMesh_pre, zTopo - zMesh_pre)
shading flat
axis equal
axis tight
hold on
contour(xMesh_pre, yMesh_pre, zMap, 50, 'k')
clim([0 max(zTopo - zMesh_pre,[],"all")])
c = colorbar;
ylabel(c,'\Delta z (m)')
xlabel('Easting (m)')
ylabel('Northing (m)')