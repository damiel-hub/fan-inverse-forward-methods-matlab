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
% xyzVisPolygon - Cell array of 3D coordinates (`x`, `y`, `z`) for visibility polygons. Only generated if `saveVisPolygon` is set to `1`.
% xyVisPolygonAll - Matrix of `x` and `y` coordinates for all visibility polygons. Only generated if `saveVisPolygon` is set to `1`.
% thetaMesh - 2D matrix of angular distribution relative to apex(es).
% Tzu-Yin Kasha Chen, March 2019; modified Aug 2022
% Yuan-Hung Chiu modified Aug 2024; modified Dec 2025

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
    options.runthetaMesh = 0
end

% 1. Grid Pre-processing
flip_lr = xMesh(1,1) > xMesh(1,end);
flip_ud = yMesh(1,1) > yMesh(end,1);
if flip_lr; xMesh = fliplr(xMesh); zMesh = fliplr(zMesh); end
if flip_ud; yMesh = flipud(yMesh); zMesh = flipud(zMesh); end

% Add a high wall around the domain
xMin = min(min(xMesh)); xMax = max(max(xMesh));
yMin = min(min(yMesh)); yMax = max(max(yMesh));
zMin = min(min(zMesh)); zMax = max(max(zMesh));
zMax = max(max(zApexM),zMax);
zMesh0 = zMesh;

% Calculate Grid Spacing
dxMesh = (xMax-xMin)/(size(xMesh,2)-1);
dyMesh = (yMax-yMin)/(size(yMesh,1)-1);
[xMesh,yMesh] = meshgrid((xMin-dxMesh):dxMesh:(xMax+dxMesh),(yMin-dyMesh):dyMesh:(yMax+dyMesh));
zMesh = ones(size(xMesh))*zMax;
zMesh(2:end-1,2:end-1) = zMesh0;
zMesh(isnan(zMesh)) = zMax;
F = griddedInterpolant(xMesh', yMesh', zMesh');

% =========================================================================
% DYNAMIC CONFIGURATION BLOCK
% =========================================================================
avg_dMesh = (dxMesh + dyMesh) / 2;

config = struct();
% Level to determine "intersection" with the ground (was 0 or 1e-6)
config.contourLevel = 1e-6; 

% Minimum nodes to consider a polygon valid (was 5)
config.minNodeCount = 5; 

% Threshold for visibility polygon simplification (was min(dx,dy)/2)
config.visiThreshold = avg_dMesh / 2;

% Minimum distance a new child apex must be from the visibility boundary.
% Previously: sqrt(2)*dxMesh*2 (~2.8*dx). 
% If this is too small, infinite tiny steps occur.
config.minDistFromBoundary = 3.0 * avg_dMesh; 

% Minimum separation between a new apex and existing apexes to merge them.
% Previously: derived dynamically from min_d_xyVisi/4
config.minApexSeparation = avg_dMesh / 2; 

% Tolerance for checking if x or y coordinates align (Grid snapping).
% Previously: dxMesh/8
config.alignmentTol = avg_dMesh / 10;
% =========================================================================

% Initialize Output
xyzkApexAll = [];
kTopoAll = nan(size(zMesh));
xyzVisPolygon = {};
xyVisPolygonAll = [];
zTopo = nan(size(zMesh));
thetaMesh = nan(size(zMesh));

if ~options.dispflag
    figure; axis equal; axis([xMin, xMax, yMin, yMax]); clim([zMin, zMax]);
end

for jj = 1:length(zApexM)
    kTopo = zeros(size(zMesh));
    xyzkApex = [xApexM(jj), yApexM(jj), zApexM(jj), nan];
    kApex = 1;
    
    % Loop over apexes (Main Simulation Loop)
    while kApex <= size(xyzkApex,1)
        
        % Select active apex
        xApex = xyzkApex(kApex,1);
        yApex = xyzkApex(kApex,2);
        zApex = xyzkApex(kApex,3);
        
        % Calculate Cone Surface
        D = sqrt( (xMesh-xApex).^2 + (yMesh-yApex).^2 );
        zCone = coneFunction(zApex,D, 'caseName', options.caseName,'tanAlpha', options.tanAlphaM(jj), 'K', options.KM(jj), 'zApex0', zApexM(jj), 'tanInfinite', options.tanInfiniteM(jj), 'dz_interp', options.dz_interpM{jj});
        
        % Find intersection polygon (Cone vs Boundary)
        % Using config.contourLevel instead of 0 for consistency
        C = contour(xMesh,yMesh,zCone-zMesh,[config.contourLevel, config.contourLevel],'Visible','off');
        
        hasNodes = ~isempty(C) && size(C,2) > 1;
        
        if hasNodes
            % Parse contour matrix C
            n_nodes_vec = [];
            idx = 1; 
            kNan = [];
            while idx < size(C,2)
                len = C(2, idx);
                n_nodes_vec = [n_nodes_vec, len];
                kNan = [kNan; idx + len + 1];
                idx = idx + len + 1;
            end
            kNan(kNan > size(C,2)) = [];
            
            % Check if the contour is large enough to matter
            if max(n_nodes_vec) > config.minNodeCount
                xContour = C(1,:)'; xContour(kNan) = nan; xContour(1) = [];
                yContour = C(2,:)'; yContour(kNan) = nan; yContour(1) = [];
                
                % Compute Visibility Polygon using dynamic threshold
                [xVisi,yVisi,xChildApex,yChildApex] = visiPolygon_optimized(xContour,yContour,xApex,yApex, config.visiThreshold, 0);
                
                if length(xVisi) > config.minNodeCount
                    
                    if options.saveVisPolygon
                        D_Visi = sqrt( (xVisi-xApex).^2 + (yVisi-yApex).^2 );
                        zVisi = coneFunction(zApex,D_Visi, 'caseName', options.caseName,'tanAlpha', options.tanAlphaM(jj), 'K', options.KM(jj), 'zApex0', zApexM(jj), 'tanInfinite', options.tanInfiniteM(jj), 'dz_interp', options.dz_interpM{jj});
                        xyzVisPolygon{end+1} = [xVisi yVisi zVisi];
                    end
                    
                    % Update fan surface
                    [NODE, EDGE] = getNodeAndEdge(xVisi, yVisi);
                    [isVisible, onVisible] = inpolygon_mex_v23(xMesh(1,:), yMesh(:,1), xVisi, yVisi);
                    isVisible = isVisible | onVisible;
                    
                    mask = isVisible & (zCone > zTopo | isnan(zTopo));                
                    zTopo(mask) = zCone(mask);
                    kTopo(zCone==zTopo) = kApex;
                    
                    if options.runthetaMesh
                        thetaMesh_temp = atan2(xMesh - xApex, yMesh - yApex);
                        thetaMesh(mask) = thetaMesh_temp(mask);
                    end
                    
                    if options.saveVisPolygon
                        if isempty(xyVisPolygonAll)
                            xyVisPolygonAll = [xVisi, yVisi];
                        else
                            isVisible = inpoly2([xyVisPolygonAll(:,1), xyVisPolygonAll(:,2)], NODE, EDGE);
                            xyVisPolygonAll(isVisible, :) = [];
                            xyVisPolygonAll = [xyVisPolygonAll; [xVisi, yVisi]];
                        end
                    end
                    
                    % -----------------------------------------------------
                    % ADD CHILD APEXES (Using Dynamic Tolerances)
                    % -----------------------------------------------------
                    % Get boundary of current visibility
                    CTopo = contourc(xMesh(1,:),yMesh(:,1),kTopo,[config.contourLevel, config.contourLevel]); 
                    CTopo(:,CTopo(1,:)==config.contourLevel) = [];
                    
                    if ~isempty(CTopo) && size(CTopo, 2) > 0
                        
                        D = sqrt( (xChildApex-xApex).^2 + (yChildApex-yApex).^2 );
                        zConeChildApex = coneFunction(zApex,D, 'caseName', options.caseName,'tanAlpha', options.tanAlphaM(jj), 'K', options.KM(jj), 'zApex0', zApexM(jj), 'tanInfinite', options.tanInfiniteM(jj), 'dz_interp', options.dz_interpM{jj});
                        
                        for i = 1:length(xChildApex)
                            % Distance to the previous boundary
                            dist_CTopo = min(sqrt((CTopo(1,:)-xChildApex(i)).^2+(CTopo(2,:)-yChildApex(i)).^2));
                            
                            % STABILITY FIX: Use config.minDistFromBoundary
                            % Ensures we don't create children too close to the edge (micro-stepping)
                            if dist_CTopo < config.minDistFromBoundary 
                                
                                dx_exist = xyzkApex(:,1)-xChildApex(i);
                                dy_exist = xyzkApex(:,2)-yChildApex(i);
                                ds_exist = sqrt(dx_exist.^2 + dy_exist.^2); 
                                
                                % STABILITY FIX: Use config.minApexSeparation
                                isTooClose = find(ds_exist < config.minApexSeparation); 
                                isTooClose(isTooClose == kApex) = [];
                                
                                % STABILITY FIX: Use config.alignmentTol for robust coordinate comparison
                                isSameXorY = find((abs(dx_exist) < eps & abs(dy_exist) < config.alignmentTol) | ...
                                                  (abs(dx_exist) < config.alignmentTol & abs(dy_exist) < eps));
                                isSameXorY(isSameXorY == kApex) = [];
                                
                                if ~isempty(isTooClose) || ~isempty(isSameXorY)
                                    % Update existing apexes (merge)
                                    idx_update = unique([isTooClose; isSameXorY]);
                                    D_exist = sqrt((xyzkApex(idx_update,1)-xApex).^2+(xyzkApex(idx_update,2)-yApex).^2);
                                    z_new = coneFunction(zApex,D_exist, 'caseName', options.caseName,'tanAlpha', options.tanAlphaM(jj), 'K', options.KM(jj), 'zApex0', zApexM(jj), 'tanInfinite', options.tanInfiniteM(jj), 'dz_interp', options.dz_interpM{jj});
                                    xyzkApex(idx_update,3) = max(xyzkApex(idx_update,3), z_new);
                                else
                                    % Add new semi-apex
                                    xyzkApex = [xyzkApex; xChildApex(i) yChildApex(i) zConeChildApex(i) kApex];
                                end
                            end
                        end
                    end
                    
                    % Remove buried apexes
                    if sum(~isnan(zTopo(:))) <= 4 
                        zAtopo = mean(zTopo(:), 'omitmissing');
                    else
                        F.Values = zTopo';
                        zAtopo = F(xyzkApex(:,1),xyzkApex(:,2));
                    end
                    
                    if isnan(zAtopo); zAtopo = NaN; end
                    
                    % STABILITY FIX: Use config.minDistFromBoundary for height check
                    zAtopo_vale = coneFunction(zAtopo, config.minDistFromBoundary, 'caseName', options.caseName,'tanAlpha', options.tanAlphaM(jj), 'K', options.KM(jj), 'zApex0', zApexM(jj), 'tanInfinite', options.tanInfiniteM(jj), 'dz_interp', options.dz_interpM{jj});
                    
                    toRemove = xyzkApex(:,3) < zAtopo_vale;
                    toRemove(isnan(toRemove)) = 0;
                    xyzkApex(toRemove,:) = [];
                    
                    % Sort apexes by elevation
                    if size(xyzkApex,1)>kApex
                        xyzkApex(kApex+1:end,:) = sortrows(xyzkApex(kApex+1:end,:),3,'descend');
                    end
                else
                    if options.saveVisPolygon; xyzVisPolygon{end+1} = [nan nan nan]; end
                end
            else
                if options.saveVisPolygon; xyzVisPolygon{end+1} = [nan nan nan]; end
            end
        else
            if options.saveVisPolygon; xyzVisPolygon{end+1} = [nan nan nan]; end
        end
        
        % Plotting
        if options.dispflag
            plot(xVisi, yVisi, 'g-'); plot(xApex, yApex, 'ko');
            plot(xyzkApex(:, 1), xyzkApex(:, 2), 'k.');
            plot(xChildApex, yChildApex, 'kv');
            title(['Apex no. ', int2str(kApex)]);
            drawnow; 
        end
        
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

if ~options.dispflag; close; end

% Post-processing (remove walls and flip back)
zTopo = zTopo(2:end-1, 2:end-1);
thetaMesh = thetaMesh(2:end-1, 2:end-1);
kTopoAll = kTopoAll(2:end-1, 2:end-1);
kTopoAll(kTopoAll == 0) = nan;

if flip_ud; zTopo = flipud(zTopo); thetaMesh = flipud(thetaMesh); kTopoAll = flipud(kTopoAll); end
if flip_lr; zTopo = fliplr(zTopo); thetaMesh = fliplr(thetaMesh); kTopoAll = fliplr(kTopoAll); end

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