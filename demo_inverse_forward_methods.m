close all; clc; clear;

% Add necessary paths for function dependencies
addpath(genpath('functions'))

% INPOLY: A fast point-in-polygon function developed by Darren Engwirda
% 
% This function provides a faster alternative to MATLAB's inpolygon() using 
% a crossing-number test with improved computational efficiency.
% 
% Source: https://github.com/dengwirda/inpoly
% License: Free for private, research, and institutional use with restrictions 
%          (see inpoly.m for details). Commercial use requires permission 
%          from the author.
% 
% Reference:
%     J. Kepner, D. Engwirda, V. Gadepally, C. Hill, T. Kraska, M. Jones, 
%     A. Kipf, L. Milechin, N. Vembar: Fast Mapping onto Census Blocks, 
%     IEEE HPEC, 2020.
% 
% Adding inpoly to the MATLAB path:
addpath('functions/inpoly-master')

% freezeColors: Enable multiple colormaps within a single figure or axis
%
% This function allows multiple colormaps to coexist in a single figure.
%
% Developed by John Iversen (2005-2022).
%
% Source: https://github.com/jiversen/freezeColors
% Licensed under the MIT License.
%
% Adding freezeColors to the MATLAB path:
addpath('functions/freezeColors-main')

%% Set DEM resolution and file paths
DEM_resolution = 10; % Available resolutions: 10, 20 meters

% Estimated Running Time for Different DEM Resolutions:
% System 1: AMD Ryzen 9 5900X, NVIDIA GeForce RTX 3080, MATLAB R2023b
% 10-meter: 73.2 seconds, 20-meter: 28.5 seconds
%
% System 2: Intel Core i7-7700HQ, NVIDIA GeForce GTX 1050, MATLAB R2022a
% 10-meter: 176.1 seconds, 20-meter: 90.5 seconds
%
% System 3: Intel Utral 7 155H, MATLAB R2024a
% 10-meter: 91.6 seconds, 20-meter: 41.2 seconds


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
xyzsBoundary = shortest_path_distance_along_boundary(topo_post_event, shape_fan_boundary);
zBoundary = xyzsBoundary(:,3);

% Step 2:
% Fit quadratic elevation-distance relationship
bin_num = 50;
ds = 5;
outlength = 500;

% Plot and fit within boundary
figure
scatter(sMap(:), zMesh_post_crop(:), 'k.')
hold on
fitting_s_z_within_boundary = process_s_z_relationship(sMap, zMesh_post_crop, bin_num, ds, outlength, 1);
xlabel('Shortest path distance to all data points, s (m)')
ylabel('Elevation, z (m)')

% Plot and fit along boundary
figure
plot(xyzsBoundary(:,4), zBoundary, 'k-')
hold on
fitting_s_z_along_boundary = process_s_z_relationship(xyzsBoundary(:,4), zBoundary, bin_num, ds, outlength, 1, 'medianFilter', 0);
xlabel('Shortest path distance to boundary points, s (m)')
ylabel('Elevation, z (m)')

%% Forward Method: Reconstruct the debris and alluvial fan topography

% Step 3
% Load initial topography
[~, ~, zMesh_pre] = readGeoTiff(topo_pre_event);
[~, ~, zMesh_pre_crop] = clipGeoTiff(topo_pre_event, shape_fan_boundary);


% Determine the apex location
[~, iApex] = max(zMesh_post_crop(:));
xApex = xMesh_crop(iApex);
yApex = yMesh_crop(iApex);


zApex = interp1(fitting_s_z_within_boundary(:,1), fitting_s_z_within_boundary(:,2), 0);
[zTopo,~,~,~,~,~] = FanTopo(xMesh, yMesh, zMesh_pre, xApex, yApex, zApex, "caseName", 'myProfile', 'dz_interpM', {fitting_s_z_within_boundary});

%%
% Plot results
plotFanTopoResults(xMesh, yMesh, zTopo, zMesh_pre, xApex, yApex); % Volume within simulation area
plotFanTopoResults(xMesh, yMesh, zTopo, zMesh_pre, xApex, yApex, shape_fan_boundary); % Volume within given boundary
plotFanTopoResults(xMesh, yMesh, zMesh_post, zMesh_pre, xApex, yApex, shape_fan_boundary); % Elevation difference before and after event

toc