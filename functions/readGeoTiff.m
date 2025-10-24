function [xMesh,yMesh,zMesh] = readGeoTiff(filename)
% read terrain mesh data from GeoTiff file (*.tif)
% edit by Yuan Hung Damiel Chiu 2022/02/07

try
    % readgeorastre is Better function to read geoTiff, but it is not avaible in old version Matlab
    [Map,R]   = readgeoraster(filename); 
catch
    [Map,R]   = geotiffread(filename);
end
dx        = R.CellExtentInWorldX;
dy        = R.CellExtentInWorldY;
xmin      = R.XWorldLimits(1);
xmax      = R.XWorldLimits(2);
ymin      = R.YWorldLimits(1);
ymax      = R.YWorldLimits(2);
xgrid     = xmin+dx/2:dx:xmax-dx/2;
ygrid     = ymin+dy/2:dy:ymax-dy/2;
[xMesh,yMesh] = meshgrid(xgrid,ygrid);
% zMesh     = double(Map);
zMesh = Map;
clear Map
zMesh(zMesh <= -9999) = nan;
zMesh(zMesh == 0) = nan;
zMesh(zMesh > 1000000) = nan;
if strcmp(R.ColumnsStartFrom,'north')
    zMesh=single(flipud(zMesh));
end
if strcmp(R.RowsStartFrom,'east')
    zMesh=single(fliplr(zMesh));
end
end

% function [xMesh,yMesh,zMesh] = readGeoTiff(inputGeoTiffPath)
%     % last edit on 2024/8/6 by Yuan-Hung, Chiu (Damiel)
% 
%     % Read the georeferenced raster data
%     [zMesh, R] = readgeoraster(inputGeoTiffPath);
%     % Generate grid coordinates based on the spatial reference
%     [xMesh, yMesh] = worldGrid(R);
% end