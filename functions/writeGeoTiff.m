function writeGeoTiff(outputData,outputTiffName,coordRefSysCode,xCorMin,xCorMax,yCorMin,yCorMax,colStart,rowStart)

% colStart='south' or 'north'
% rowStart='west' or 'east'
% coordRefSysCode = 32720 (PCS_WGS84_UTM_zone_20S)
% coordRefSysCode = 3826 (TWD97)
% last edit on 24/08/05 by Damiel

% TWD97 geoTags
%{
geoTags.GTModelTypeGeoKey = 1;
geoTags.GTRasterTypeGeoKey = 1;
geoTags.GTCitationGeoKey = 'PCS Name = TWD97_TM2_zone_121';
geoTags.GeographicTypeGeoKey = 32767;
geoTags.GeogCitationGeoKey = 'GCS Name = GCS_TWD97|Datum = D_TWD_1997|Ellipsoid = GRS_1980|Primem = Greenwich|';
geoTags.GeogGeodeticDatumGeoKey = 1026;
geoTags.GeogPrimeMeridianGeoKey = 8901;
geoTags.GeogAngularUnitsGeoKey = 9102;
geoTags.GeogAngularUnitSizeGeoKey = 0.0175;
geoTags.GeogEllipsoidGeoKey = 7019;
geoTags.GeogSemiMajorAxisGeoKey = 6378137;
geoTags.GeogInvFlatteningGeoKey = 298.2572;
geoTags.GeogPrimeMeridianLongGeoKey = 0;
geoTags.ProjectedCSTypeGeoKey = 32767;
geoTags.PCSCitationGeoKey = 'ESRI PE String = PROJCS["TWD97_TM2_zone_121",GEOGCS["GCS_TWD97",DATUM["D_TWD_1997",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["false_easting",250000.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",121.0],PARAMETER["scale_factor",0.9999],PARAMETER["latitude_of_origin",0.0],UNIT["m",1.0]]';
geoTags.ProjectionGeoKey = 32767;
geoTags.ProjCoordTransGeoKey = 1;
geoTags.ProjLinearUnitsGeoKey = 9001;
geoTags.ProjNatOriginLongGeoKey = 121;
geoTags.ProjNatOriginLatGeoKey = 0;
geoTags.ProjFalseEastingGeoKey = 250000;
geoTags.ProjFalseNorthingGeoKey = 0;
geoTags.ProjScaleAtNatOriginGeoKey = 0.9999;
%}

    R = maprefcells();
    dx = (xCorMax-xCorMin)/(size(outputData,2)-1);
    dy = (yCorMax-yCorMin)/(size(outputData,1)-1);
    xmin = xCorMin-dx/2;
    xmax = xCorMax+dx/2;
    ymin = yCorMin-dy/2;
    ymax = yCorMax+dy/2;
    R.CellExtentInWorldX=dx;
    R.CellExtentInWorldY=dy;
    R.XWorldLimits=[xmin,xmax];
    R.YWorldLimits=[ymin,ymax];
    R.RasterSize=size(outputData);
    R.ColumnsStartFrom = colStart;
    R.RowsStartFrom = rowStart;
    R.ProjectedCRS = projcrs(coordRefSysCode);
    geotiffwrite(outputTiffName, outputData, R,'CoordRefSysCode', coordRefSysCode);

end