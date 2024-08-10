function [zTopo,kTopoAll,xyzkApexAll,xyzVisPolygon,xyVisPolygonAll,thetaMesh] = FanTopo(xMesh,yMesh,zMesh,xApexM,yApexM,zApexM,options)
%FANTOPO_SLOPE constructs the constant-slope or concave fan morphology, the apexes
%positions, and source provenance

% >> [zTopo,kTopoAll,xyzkApexAll,xyzVisPolygon,xyVisPolygonAll,thetaMesh] = reconstruct_fan_surface(xMesh,yMesh,zMesh,xApexM,yApexM,zApexM,options)

% Inputs:
% xMesh - 2D matrix of x-coordinates for mesh grid points.
% yMesh - 2D matrix of y-coordinates for mesh grid points.
% zMesh - 2D matrix of initial elevation values before fan aggradation.
% xApexM - Vector of x-coordinates for fan apex(es).
% yApexM - Vector of y-coordinates for fan apex(es).
% zApexM - Vector of z-coordinates (elevations) for fan apex(es).
% options - Structure containing optional parameters:
%   caseName - (string) Type of fan morphology to generate (e.g., 'cone', 'concave', 'infinite', 'myProfile').
%   caseName = 'cone'
%       tanAlphaM - (vector) Slope angles (tangents) for each apex, defining fan steepness.
%   caseName = 'concave'
%       tanAlphaM - (vector) Slope angles (tangents) for each apex, defining fan steepness.
%       KM - (vector) Concavity factors for each apex, controlling the curvature of the fan.
%   caseName = 'infinite'
%       tanAlphaM - (vector) Slope angles (tangents) for each apex, defining fan steepness.
%       KM - (vector) Concavity factors for each apex, controlling the curvature of the fan.    
%       tanInfiniteM - (vector) Slope values for cases where the tangent approaches infinity.
%   caseName = 'myProfile'
%       dz_interpM - (cell array) Interpolation values for elevation, used in spline-based morphologies.
%   dispflag - (scalar) Flag to display the generated topography (1 for on, 0 for off).
%   saveVisPolygon - (scalar) Flag to save visibility polygons (1 for yes, 0 for no).

% Outputs:
% zTopo - 2D matrix of final fan topography (elevation after aggradation).
% kTopoAll - 2D matrix with indices of the apex dominating each mesh grid point.
% xyzkApexAll - Matrix of apex coordinates and indices (including child apexes).
% xyzVisPolygon - Cell array of 3D coordinates (x, y, z) for visibility polygons.
% xyVisPolygonAll - Matrix of x and y coordinates for all visibility polygons.
% thetaMesh - 2D matrix of angular distribution relative to apex(es).

% Tzu-Yin Kasha Chen, March 2019; modified Aug 2022
% Yuan-Hung Chiu modified Aug 2024

arguments
    xMesh double
    yMesh double
    zMesh double
    xApexM double
    yApexM double
    zApexM double
    options.caseName = 'cone'
    options.tanAlphaM = nan(1,length(zApexM))
    options.KM = nan(1,length(zApexM))
    options.tanInfiniteM = nan(1,length(zApexM))
    options.dz_interpM = arrayfun(@(x) nan(1, x), 1:length(zApexM), 'UniformOutput', false);
    options.dispflag = 0
    options.saveVisPolygon = 0
end


% Expand the x and y coordinates while preserving the order
flip_lr = xMesh(1,1) > xMesh(1,end);
flip_ud = yMesh(1,1) > yMesh(end,1);

if flip_lr
    xMesh = fliplr(xMesh);
    zMesh = fliplr(zMesh);
end
if flip_ud
    yMesh = flipud(yMesh);
    zMesh = flipud(zMesh);
end

% add a high wall around the domain
xMin = min(min(xMesh)); xMax = max(max(xMesh));
yMin = min(min(yMesh)); yMax = max(max(yMesh));
zMin = min(min(zMesh)); zMax = max(max(zMesh));
zMax = max(max(zApexM),zMax);
zMesh0 = zMesh;
dxMesh = (xMax-xMin)/(size(xMesh,2)-1);
dyMesh = (yMax-yMin)/(size(yMesh,1)-1);
[xMesh,yMesh] = meshgrid((xMin-dxMesh):dxMesh:(xMax+dxMesh),(yMin-dyMesh):dyMesh:(yMax+dyMesh));
zMesh = ones(size(xMesh))*zMax;
zMesh(2:end-1,2:end-1) = zMesh0;
zMesh(isnan(zMesh)) = zMax;

[nr, nc] = size(zMesh);

% initialize topography and sorted apex list:
xyzkApexAll = [];
kTopoAll = nan(size(zMesh));
xyzVisPolygon = {};
xyVisPolygonAll = [];

if ~options.saveVisPolygon
    disp('xyzVisPolygon: off');
    disp('xyVisPolygonAll: off');
end

if ~options.dispflag
    figure
end

zTopo = nan(size(zMesh));
thetaMesh = nan(size(zMesh));

for jj = 1:length(zApexM)
    kTopo = zeros(size(zMesh));
    xyzkApex = [xApexM(jj), yApexM(jj), zApexM(jj), nan];
    kApex = 1;
    
    % loop over apexes:
    while kApex<= size(xyzkApex,1)
        %%
        % select active apex:
        xApex = xyzkApex(kApex,1);
        yApex = xyzkApex(kApex,2);
        zApex = xyzkApex(kApex,3);
        
        % find intersection polygon of cone surface and boundary surface:
        D = sqrt( (xMesh-xApex).^2 + (yMesh-yApex).^2 );
        zCone = coneFunction(zApex,D, 'caseName', options.caseName,'tanAlpha', options.tanAlphaM(jj), 'K', options.KM(jj), 'zApex0', zApexM(jj), 'tanInfinite', options.tanInfiniteM(jj), 'dz_interp', options.dz_interpM{jj});
        C = contour(xMesh,yMesh,zCone-zMesh,[0,0],'Visible','off');
        n_nodes = C(2,C(1,:)==0);
        if max(n_nodes)>5 % ignore the apex whose impact is too small
            kNan = 1;
            while(kNan(end)+C(2,kNan(end))<size(C,2))
                kNan = [kNan;kNan(end)+C(2,kNan(end))+1];
            end
            xContour = C(1,:)'; xContour(kNan) = nan; xContour(1) = [];
            yContour = C(2,:)'; yContour(kNan) = nan; yContour(1) = [];
            % find visibility polygon and children apexes:
            [xVisi,yVisi,xChildApex,yChildApex] = visiPolygon(xContour,yContour,xApex,yApex,min(dxMesh,dyMesh)/10,0);
            
            if length(xVisi)>5 % ignore the apex whose impact is too small
                
                if options.saveVisPolygon
                    D_Visi = sqrt( (xVisi-xApex).^2 + (yVisi-yApex).^2 );
                    zVisi = coneFunction(zApex,D_Visi, 'caseName', options.caseName,'tanAlpha', options.tanAlphaM(jj), 'K', options.KM(jj), 'zApex0', zApexM(jj), 'tanInfinite', options.tanInfiniteM(jj), 'dz_interp', options.dz_interpM{jj});
                    xyzVisPolygon{end+1} = [xVisi yVisi zVisi];
                end
                
                % update fan surface to the visible sector occluded by boundar surface and other fan sectors:
                [NODE, EDGE] = getNodeAndEdge(xVisi, yVisi);
                isVisible = inpoly2([xMesh(:),yMesh(:)],NODE,EDGE);
                isVisible = reshape(isVisible,nr,nc);

                thetaMesh_temp = atan2(xMesh - xApex, yMesh - yApex);
                mask = isVisible & (zCone > zTopo | isnan(zTopo));
                thetaMesh(mask) = thetaMesh_temp(mask);
                zTopo(mask) = zCone(mask);
                kTopo(zCone==zTopo) = kApex;
                
                if options.saveVisPolygon
                    if isempty(xyVisPolygonAll)
                        xyVisPolygonAll = [xVisi, yVisi];
                    else
                        isVisible = inpoly2([xyVisPolygonAll(:,1), xyVisPolygonAll(:,2)], NODE, EDGE);
                        xyVisPolygonAll(isVisible, :) = [];
                        xyVisPolygonAll = [xyVisPolygonAll; [xVisi, yVisi]];
                    end
                end

                % add effective children apexes into the apex list
                CTopo = contour(xMesh,yMesh,kTopo,[1e-6,1e-6],'Visible','off'); % the boundary of previous visibility polygons
                CTopo(:,CTopo(1,:)==1e-6) = [];


                min_d_xyVisi = min(sqrt((diff(xVisi).^2+diff(yVisi).^2))); % threshold for finding the semi-apexes that are too close
                D = sqrt( (xChildApex-xApex).^2 + (yChildApex-yApex).^2 );
                zConeChildApex = coneFunction(zApex,D, 'caseName', options.caseName,'tanAlpha', options.tanAlphaM(jj), 'K', options.KM(jj), 'zApex0', zApexM(jj), 'tanInfinite', options.tanInfiniteM(jj), 'dz_interp', options.dz_interpM{jj});
                for i = 1:length(xChildApex)
                    dist_CTopo = min(sqrt((CTopo(1,:)-xChildApex(i)).^2+(CTopo(2,:)-yChildApex(i)).^2));
                    if dist_CTopo < sqrt(2)*dxMesh*2 % only keep the semi-apex that close to the boundary of previous visibility polygons
                        dx_existChildApex = xyzkApex(:,1)-xChildApex(i);
                        dy_existChildApex = xyzkApex(:,2)-yChildApex(i);
                        ds_existChildApex = sqrt(dx_existChildApex.^2+(dy_existChildApex).^2); % distance to existing semi-apexes
                        isTooClose = find(ds_existChildApex < min_d_xyVisi/4); % find the existing semi-apexes that are too close to the new semi-apex
                        isTooClose(isTooClose == kApex) = [];
                        isSameXorY = find((dx_existChildApex==0 & abs(dy_existChildApex)<dyMesh/8) | (abs(dx_existChildApex)<dxMesh/8 & dy_existChildApex==0)); % find the existing semi-apexes that have the same x or y cordination as the new semi-apex
                        isSameXorY(isSameXorY == kApex) = [];
                        if ~isempty(isTooClose) || ~isempty(isSameXorY)
                            % update the z value of the too-close/sameXorY existing semi-apexes
                            D = sqrt((xyzkApex(isTooClose,1)-xApex).^2+(xyzkApex(isTooClose,2)-yApex).^2);
                            xyzkApex(isTooClose,3) = max(xyzkApex(isTooClose,3), coneFunction(zApex,D, 'caseName', options.caseName,'tanAlpha', options.tanAlphaM(jj), 'K', options.KM(jj), 'zApex0', zApexM(jj), 'tanInfinite', options.tanInfiniteM(jj), 'dz_interp', options.dz_interpM{jj}));
                            D = sqrt((xyzkApex(isSameXorY,1)-xApex).^2+(xyzkApex(isSameXorY,2)-yApex).^2);
                            xyzkApex(isSameXorY,3) = max(xyzkApex(isSameXorY,3), coneFunction(zApex,D, 'caseName', options.caseName,'tanAlpha', options.tanAlphaM(jj), 'K', options.KM(jj), 'zApex0', zApexM(jj), 'tanInfinite', options.tanInfiniteM(jj), 'dz_interp', options.dz_interpM{jj}));
                        else
                            % add new semi-apex
                            zChildApex = zConeChildApex(i);
                            xyzkApex = [xyzkApex; xChildApex(i) yChildApex(i) zChildApex kApex];
                        end
                    end
                end
                % remove buried apexes
                zAtopo = interp2(xMesh,yMesh,zTopo,xyzkApex(:,1),xyzkApex(:,2));
                zAtopo_vale = coneFunction(zAtopo,sqrt(2)*dxMesh*2, 'caseName', options.caseName,'tanAlpha', options.tanAlphaM(jj), 'K', options.KM(jj), 'zApex0', zApexM(jj), 'tanInfinite', options.tanInfiniteM(jj), 'dz_interp', options.dz_interpM{jj});
                xyzkApex(xyzkApex(:,3)<zAtopo_vale,:) = [];

                % sort the apexes by elevation
                if size(xyzkApex,1)>kApex
                    xyzkApex(kApex+1:end,:) = sortrows(xyzkApex(kApex+1:end,:),3,'descend');
                end
            else
                if options.saveVisPolygon
                    xyzVisPolygon{end+1} = [nan nan nan];
                end
            end
        else
            if options.saveVisPolygon
                xyzVisPolygon{end+1} = [nan nan nan];
            end
        end
        % show topography, active fan contour, and visible sector and apexes:
        if options.dispflag
            axis equal
            axis([xMin, xMax, yMin, yMax])
            clim([zMin, zMax])
            hold on
            plot(xVisi, yVisi, 'g-')
            plot(xApex, yApex, 'ko')
            plot(xyzkApex(:, 1), xyzkApex(:, 2), 'k.')
            plot(xChildApex, yChildApex, 'kv')
            title(['Apex no. ', int2str(kApex)])
            hold off
            drawnow
        end

        % proceed to next apex on the list:
        kApex = kApex + 1;
    end
    if jj>1
        xyzkApex(:,4) = xyzkApex(:,4)+size(xyzkApexAll,1);
        kTopoAll(kTopo>0)=kTopo(kTopo>0)+size(xyzkApexAll,1);
    else
        kTopoAll=kTopo;
    end
    xyzkApexAll = [xyzkApexAll;xyzkApex];
    
    zMesh(~isnan(zTopo)) = zTopo(~isnan(zTopo));
end
if ~options.dispflag
    close
end

% Remove the high wall
zTopo = zTopo(2:end-1, 2:end-1);
thetaMesh = thetaMesh(2:end-1, 2:end-1);
kTopoAll = kTopoAll(2:end-1, 2:end-1);
kTopoAll(kTopoAll == 0) = nan;

% Flip back if necessary
if flip_ud
    zTopo = flipud(zTopo);
    thetaMesh = flipud(thetaMesh);
    kTopoAll = flipud(kTopoAll);
end
if flip_lr
    zTopo = fliplr(zTopo);
    thetaMesh = fliplr(thetaMesh);
    kTopoAll = fliplr(kTopoAll);
end


function [NODE, EDGE] = getNodeAndEdge(x, y)
    % getNodeAndEdge creates the NODE and EDGE arrays from x and y coordinates.
    %
    % Inputs:
    %   x - A vector of x coordinates of the polygon's vertices
    %   y - A vector of y coordinates of the polygon's vertices
    %
    % Outputs:
    %   NODE - An Mx2 array of the polygon's vertices
    %   EDGE - A Px2 array of edge indexing

    % Combine x and y into NODE array
    NODE = [x(:), y(:)];

    % Create EDGE array
    numVertices = length(x);
    EDGE = [1:numVertices; 2:numVertices+1]';
    EDGE(end) = 1;

end

end
