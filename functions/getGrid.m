function [zMesh, R, key] = getGrid(inputGeoTiffPath)
    % inputGeotiffPath (string): Path to the input geotiff file.
    % Read the georeferenced raster data
    [zMesh, R] = readgeoraster(inputGeoTiffPath);

    % Get the GeoTIFF metadata
    info = geotiffinfo(inputGeoTiffPath);
    key = info.GeoTIFFTags.GeoKeyDirectoryTag;

end