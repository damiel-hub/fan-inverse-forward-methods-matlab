function [xMesh_crop, yMesh_crop, zMesh_crop] = clipGeoTiff(inputGeoTiffPath, shapeFilePath, outputGeoTiffPath)
    % clipGeoTiff Clips a geotiff based on a polygon shapefile.
    %
    % Inputs:
    %   inputGeotiffPath (string): Path to the input geotiff file.
    %   shapeFilePath (string): Path to the shapefile defining the clipping fan boundary polygon.
    %   outputGeotiffPath (string, optional): Path to save the output clipped geotiff.
    %                                          If not provided, the output will not be saved.
    
    % last edit on 2024/8/6 by Yuan-Hung, Chiu (Damiel)

    % Read the georeferenced raster data
    [zMesh, R] = readgeoraster(inputGeoTiffPath);

    % Get the GeoTIFF metadata
    info = geotiffinfo(inputGeoTiffPath);
    key = info.GeoTIFFTags.GeoKeyDirectoryTag;

    % Generate grid coordinates based on the spatial reference
    [xMesh, yMesh] = worldGrid(R);

    % Read the shapefile containing the fan boundary polygon
    fan_boundary = shaperead(shapeFilePath);
    fan_boundary_x = fan_boundary.X(1:end-1);
    fan_boundary_y = fan_boundary.Y(1:end-1);

    % Create a logical mask for points inside the polygon
    mask = inpolygon(xMesh, yMesh, fan_boundary_x, fan_boundary_y);

    % Clip the elevation data to the fan extent
    zMesh_clip = zMesh;
    zMesh_clip(~mask) = nan;

    % Crop the raster to the bounding box of the fan extent
    [zMesh_crop, R_crop] = mapcrop(zMesh_clip, R, [min(fan_boundary_x) max(fan_boundary_x)], [min(fan_boundary_y) max(fan_boundary_y)]);
    
    % Generate grid coordinates for the cropped data
    [xMesh_crop, yMesh_crop] = worldGrid(R_crop);

    % If outputGeotiffPath is provided, write the clipped and cropped data to a new GeoTIFF
    if nargin > 2 && ~isempty(outputGeoTiffPath)
        geotiffwrite(outputGeoTiffPath, zMesh_crop, R_crop, GeoKeyDirectoryTag=key);
    end
end
