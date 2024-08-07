function poly_xyz_all = get_3D_outline_from_FanTopo_slope_bd_cone(xyzkApex,xyzVisi,tanAlpha)

warning('off', 'all');
polyUnion = polyshape(xyzVisi{1}(:,1),xyzVisi{1}(:,2));
xVisiUnion = polyUnion.Vertices(:,1);
yVisiUnion = polyUnion.Vertices(:,2);
D = sqrt((xVisiUnion-xyzkApex(1,1)).^2 + (yVisiUnion-xyzkApex(1,2)).^2);
zVisiUnion = coneFunction(xyzkApex(1,3),D, 'caseName', 'cone', 'tanAlpha', tanAlpha);

for i = 2:length(xyzVisi)
    polyUnion = polyshape(xVisiUnion,yVisiUnion);
    poly1 = polyshape(xyzVisi{i}(:,1),xyzVisi{i}(:,2));
    polyX = poly1.Vertices(:,1);
    polyY = poly1.Vertices(:,2);
    D = sqrt((polyX-xyzkApex(i,1)).^2 + (polyY-xyzkApex(i,2)).^2);
    polyZ = coneFunction(xyzkApex(i,3),D, 'caseName', 'cone', 'tanAlpha', tanAlpha);
    
    [polyUnion,shapeID,vertexID] = union([polyUnion poly1]);

    xVisiUnion = polyUnion.Vertices(:,1);
    yVisiUnion = polyUnion.Vertices(:,2);
    
    max_length = max(length(zVisiUnion),length(polyZ));
    zVisi_All = nan(max_length,2);
    zVisi_All(1:length(zVisiUnion),1) = zVisiUnion;
    zVisi_All(1:length(polyZ),2) = polyZ;
    
    vertexID(shapeID==0) = nan;
    shapeID(shapeID==0) = nan;
    
    zVisi_index = sub2ind(size(zVisi_All),vertexID,shapeID);
    nanMask = isnan(zVisi_index);
    zVisi_index(nanMask) = 1;
    zVisi_All = zVisi_All(zVisi_index);
    zVisi_All(nanMask) = nan;
    zVisiUnion = zVisi_All;
end
warning('on', 'all');


nan_index = find(isnan(polyUnion.Vertices(:,1)));
nan_index = [0;nan_index;length(polyUnion.Vertices(:,1))+1];
start_index = nan_index + 1;
end_index = nan_index - 1;
start_index(end) = [];
end_index(1) = [];

poly_xyz_all = {};
for i = 1:length(start_index)
    poly_x = polyUnion.Vertices(start_index(i):end_index(i),1);
    poly_y = polyUnion.Vertices(start_index(i):end_index(i),2);
    poly_z = zVisiUnion(start_index(i):end_index(i));
    nanMask = isnan(poly_z);
    poly_x(nanMask) = [];
    poly_y(nanMask) = [];
    poly_z(nanMask) = [];

    poly_xyz_all{end+1} = [poly_x,poly_y,poly_z;poly_x(1),poly_y(1),poly_z(1)];
    hold on
end