function [xMesh_crop, yMesh_crop, zMesh_crop] = clipGeoGrid(zMesh, R, key, shapeFilePath, outputGeoTiffPath)
    % clipGeoTiff Clips a geotiff based on a polygon shapefile.
    %
    % Inputs:
    %   xMesh, yMesh, zMesh, R, key: obtain from getGrid function
    %   shapeFilePath (string): Path to the shapefile defining the clipping fan boundary polygon.
    %   outputGeotiffPath (string, optional): Path to save the output clipped geotiff.
    %                                          If not provided, the output will not be saved.
    
    % last edit on 2025/1/2 by Yuan-Hung, Chiu (Damiel)

    % Read the shapefile containing the fan boundary polygon
    fan_boundary = shaperead(shapeFilePath);
    fan_boundary_x = fan_boundary.X(1:end-1);
    fan_boundary_y = fan_boundary.Y(1:end-1);

    % Crop the raster to the bounding box of the fan extent
    [zMesh_crop, R_crop] = mapcrop(zMesh, R, [min(fan_boundary_x) max(fan_boundary_x)], [min(fan_boundary_y) max(fan_boundary_y)]);
    
    % Generate grid coordinates for the cropped data
    [xMesh_crop, yMesh_crop] = worldGrid(R_crop);

    % Create a logical mask for points inside the polygon
    mask = inpolygon(xMesh_crop, yMesh_crop, fan_boundary_x, fan_boundary_y);
    zMesh_crop(~mask) = nan;

    % If outputGeotiffPath is provided, write the clipped and cropped data to a new GeoTIFF
    if nargin > 2 && ~isempty(outputGeoTiffPath)
        geotiffwrite(outputGeoTiffPath, zMesh_crop, R_crop, GeoKeyDirectoryTag=key);
    end
end
