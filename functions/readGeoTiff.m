function [xMesh,yMesh,zMesh] = readGeoTiff(inputGeoTiffPath)
    % last edit on 2024/8/6 by Yuan-Hung, Chiu (Damiel)
    
    % Read the georeferenced raster data
    [zMesh, R] = readgeoraster(inputGeoTiffPath);
    % Generate grid coordinates based on the spatial reference
    [xMesh, yMesh] = worldGrid(R);
end