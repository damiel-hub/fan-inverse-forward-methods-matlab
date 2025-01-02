function [xMesh, yMesh, zMesh, R, key] = getGrid(inputGeoTiffPath)
    % inputGeotiffPath (string): Path to the input geotiff file.
    % Read the georeferenced raster data
    [zMesh, R] = readgeoraster(inputGeoTiffPath);

    % Get the GeoTIFF metadata
    info = geotiffinfo(inputGeoTiffPath);
    key = info.GeoTIFFTags.GeoKeyDirectoryTag;

    % Generate grid coordinates based on the spatial reference
    [xMesh, yMesh] = worldGrid(R);
end