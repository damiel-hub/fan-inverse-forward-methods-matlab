function [xVisi,yVisi,xEffCorner,yEffCorner] = visiPolygon_optimized(xBoundary,yBoundary,xApex,yApex,threshold,dispflag)
%VISIPOLYGON calculates the visibility polygon of an apex in the given
% polygon
%
% >> [xVisi,yVisi,xEffCorner,yEffCorner] = visiPolygon(xBoundary,yBoundary,xApex,yApex,threshold,dispflag)
%
% This is an optimized version that avoids growing arrays in loops.
%
% Inputs:
% [xBoundary,yBoundary] = x, y coordinates of boundary vertices, vertices
% of various cycles of edges are seperated by NaN. A closed region is
% bounded by a cycle of edges linked by vertices in counterclockwise order.
% In contrast, the interior obstacles are bounded by a cycles of edges
% linked by vertices in clockwise order.
% [xApex,yApex] = x, y coordinates of the apex (view point)
% threshold = referenced distance threshold
% dispflag = 1 to request view of image; = 0 for no view
% Outputs:
% [xVisi,yVisi] = x, y coordinates of the visibility polygon of the apex
% [xEffCorner,yEffCorner] = x, y coordinates of the effective corners
%
% H. Capart, 2018/10/26; revised by T.Y.K. Chen, 2020/05/16
% Optimized by pre-allocation, 2025

%% construct doubly linked list of vertices:
xBoundary = xBoundary(:);
yBoundary = yBoundary(:);
kNan = find(isnan(xBoundary));
nBoundary = length(xBoundary);
nVertex = nBoundary - sum(isnan(xBoundary));
mVertex = 1;
xVertex = nan*ones(nVertex,1);
yVertex = nan*ones(nVertex,1);
kNext = nan*ones(nVertex,1);
kPrev = nan*ones(nVertex,1);
mBoundary = [1; kNan+1];
nBoundary = [kNan-1; nBoundary];
cwBoundary = zeros(1,length(mBoundary));
xyBoundary = cell(1,length(mBoundary));
for kCycle = 1:length(mBoundary)
    xyBoundary{kCycle} = [xBoundary(mBoundary(kCycle):nBoundary(kCycle)),yBoundary(mBoundary(kCycle):nBoundary(kCycle))];
    xyBoundary{kCycle} = unique(xyBoundary{kCycle},'rows','stable');
    cwBoundary(kCycle) = ispolycw(xyBoundary{kCycle}(:,1),xyBoundary{kCycle}(:,2));
end
ccwID = find(~cwBoundary);
cwID = find(cwBoundary);
closeBoundaryID = [];
for kCycle = ccwID
    if length(xyBoundary{kCycle})>=3
        isInApex = inpolygon(xApex,yApex,xyBoundary{kCycle}(:,1),xyBoundary{kCycle}(:,2));
        if isInApex
            closeBoundaryID = kCycle;
        end
    end
end
if isempty(closeBoundaryID)
    ccwdist = nan(1,length(mBoundary));
    for kCycle = ccwID
        distVertex = sqrt( (xApex-xyBoundary{kCycle}(:,1)).^2 + (yApex-xyBoundary{kCycle}(:,2)).^2 ); 
        ccwdist(kCycle) = min(distVertex);
    end
    if min(ccwdist)<3*threshold
        [~,closeBoundaryID]=min(ccwdist);
    else
        xVisi=[];yVisi=[];xEffCorner = [];yEffCorner = [];
        return;
    end
end
nVertex = mVertex + length(xyBoundary{closeBoundaryID}) - 1;
xVertex(mVertex:nVertex) = xyBoundary{closeBoundaryID}(:,1);
yVertex(mVertex:nVertex) = xyBoundary{closeBoundaryID}(:,2);
kPrev(mVertex:nVertex) = [nVertex, (mVertex:nVertex-1)]';
kNext(mVertex:nVertex) = [(mVertex+1:nVertex), mVertex]';
mBoundary(closeBoundaryID) = mVertex;
nBoundary(closeBoundaryID) = nVertex;
mVertex = nVertex+1;
for kCycle = cwID
    if length(xyBoundary{kCycle})>=3
        isInPoly = inpolygon(xyBoundary{kCycle}(1,1),xyBoundary{kCycle}(1,2),xyBoundary{closeBoundaryID}(:,1),xyBoundary{closeBoundaryID}(:,2));
        if isInPoly
            nVertex = mVertex + length(xyBoundary{kCycle}) - 1;
            xVertex(mVertex:nVertex) = xyBoundary{kCycle}(:,1);
            yVertex(mVertex:nVertex) = xyBoundary{kCycle}(:,2);
            kPrev(mVertex:nVertex) = [nVertex, (mVertex:nVertex-1)]';
            kNext(mVertex:nVertex) = [(mVertex+1:nVertex), mVertex]';
            mBoundary(kCycle) = mVertex;
            nBoundary(kCycle) = nVertex;
            mVertex = nVertex+1;
        end     
    end
end
mBoundary(isnan(mBoundary)) = [];
xVertex(mVertex:end) = [];
yVertex(mVertex:end) = [];
kNext(mVertex:end) = [];
kPrev(mVertex:end) = [];
if isempty(xVertex)
    xVisi=[];yVisi=[];xEffCorner = [];yEffCorner = [];
    return;
end
%% modify the boundary to pass through the apex
distVertex = sqrt( (xApex-xVertex).^2 + (yApex-yVertex).^2 ); % distance from the Apex to each vertex
[distMinVertex,kNearestVertex] = min(distVertex);
uEdge = xVertex(kNext)-xVertex;
vEdge = yVertex(kNext)-yVertex;
uNearest = -vEdge;
vNearest = uEdge;
detNearest = -uNearest.*vEdge + vNearest.*uEdge;
lambda = (-(xVertex-xApex).*vEdge + (yVertex-yApex).*uEdge)./detNearest;
mu = (uNearest.*(yVertex-yApex) - vNearest.*(xVertex-xApex))./detNearest;
xEdge = xApex + lambda.*uNearest; % position of the perpendicular foot
yEdge = yApex + lambda.*vNearest;
distEdge = sqrt( (xApex-xEdge).^2 + (yApex-yEdge).^2 ); % distance from the Apex to each edge
isOnEdge = (mu>0)&(mu<1);
distEdge(~isOnEdge) = Inf;
[distMinEdge,kNearestEdge] = min(distEdge);
isInApex = inpolygon(xApex,yApex,xBoundary,yBoundary);
if isInApex % Apex in polygon
    if distMinVertex < 10^-1*threshold
        kApex = kNearestVertex;
        xVertex(kApex,1) = xApex;
        yVertex(kApex,1) = yApex;
    elseif distMinEdge < 10^-1*threshold
        kApex = length(xVertex) + 1;
        xVertex(kApex,1) = xApex;
        yVertex(kApex,1) = yApex;
        kPrev(kApex) = kNearestEdge;
        kPrev(kNext(kNearestEdge)) = kApex;
        kNext(kApex) = kNext(kNearestEdge);
        kNext(kNearestEdge) = kApex;
    elseif distMinVertex < distMinEdge
        kApex = length(xVertex) + 1;
        xVertex(kApex,1) = xApex;
        yVertex(kApex,1) = yApex;
        xVertex(kApex+1,1) = xVertex(kNearestVertex);
        yVertex(kApex+1,1) = yVertex(kNearestVertex);
        kPrev(kNext(kNearestVertex)) = kApex + 1;
        kPrev(kApex+1) = kApex;
        kPrev(kApex) = kNearestVertex;
        kNext(kApex+1) = kNext(kNearestVertex);
        kNext(kApex) = kApex + 1;
        kNext(kNearestVertex) = kApex;
    else
        kApex = length(xVertex) + 1;
        xVertex(kApex,1) = xApex;
        yVertex(kApex,1) = yApex;
        kIntersect1 = kApex + 1;
        xVertex(kIntersect1,1) = xEdge(kNearestEdge);
        yVertex(kIntersect1,1) = yEdge(kNearestEdge);
        kIntersect2 = kApex + 2;
        xVertex(kIntersect2,1) = xEdge(kNearestEdge);
        yVertex(kIntersect2,1) = yEdge(kNearestEdge);
        kPrev(kNext(kNearestEdge)) = kIntersect2;
        kNext(kIntersect2) = kNext(kNearestEdge);
        kPrev(kIntersect2) = kApex;
        kNext(kApex) = kIntersect2;
        kPrev(kApex) = kIntersect1;
        kNext(kIntersect1) = kApex;
        kPrev(kIntersect1) = kNearestEdge;
        kNext(kNearestEdge) = kIntersect1;
    end
else % Apex outside polygon
    uRay = xVertex(kNearestVertex) - xApex;
    vRay = yVertex(kNearestVertex) - yApex;
    uNext = xVertex(kNext(kNearestEdge)) - xVertex(kNearestEdge);
    vNext = yVertex(kNext(kNearestEdge)) - yVertex(kNearestEdge);
    detIntersect = -uRay*vNext + vRay*uNext;
    lambda = ( -(xVertex(kNearestEdge)-xApex)*vNext + (yVertex(kNearestEdge)-yApex)*uNext )/detIntersect;
    mu = ( uRay*(yVertex(kNearestEdge)-yApex) - vRay*(xVertex(kNearestEdge)-xApex) )/detIntersect;
    if lambda > 0  && lambda < 1 && mu > 0  && mu < 1 % link between the apex and the nearest vertex intersects the nearest edge
        kApex = length(xVertex) + 1;
        xVertex(kApex,1) = xApex;
        yVertex(kApex,1) = yApex;
        kPrev(kNext(kNearestEdge)) = kApex;
        kPrev(kApex) = kNearestEdge;
        kNext(kApex) = kNext(kNearestEdge);
        kNext(kNearestEdge) = kApex;
    elseif distMinVertex < threshold/2 || (distMinVertex <= distMinEdge && distMinVertex<=3*threshold)
        kApex = kNearestVertex;
        xVertex(kApex,1) = xApex;
        yVertex(kApex,1) = yApex;
    elseif distMinVertex<=3*threshold
        kApex = length(xVertex) + 1;
        xVertex(kApex,1) = xApex;
        yVertex(kApex,1) = yApex;
        kPrev(kNext(kNearestEdge)) = kApex;
        kNext(kApex) = kNext(kNearestEdge);
        kPrev(kApex) = kNearestEdge;
        kNext(kNearestEdge) = kApex;
    else
        xVisi=[];yVisi=[];xEffCorner = [];yEffCorner = [];
        return;
    end
end
if dispflag
    figure
    plot([xVertex,xVertex(kNext)]',[yVertex,yVertex(kNext)]','k-')
    hold on
    plot(xApex,yApex,'ko')
    axis equal
    axis tight
end

% --- OPTIMIZATION 1: Vertex-to-Boundary Map ---
% Create a lookup map to avoid repeated 'find' calls inside loops
% This maps each vertex index to its boundary cycle index
nAllVertices = length(xVertex); % Get size *before* pre-allocation
kVertexToBoundary = zeros(nAllVertices, 1);
for k = 1:length(mBoundary)
    kVertexToBoundary(mBoundary(k):nBoundary(k)) = k;
end
% --- END OPTIMIZATION 1 ---

%% find corners
tol_parallel = -10^-10; %tolerance for parallel edge and sightline
uRay = xVertex-xApex;
vRay = yVertex-yApex;
uPrev = xVertex-xVertex(kPrev);
vPrev = yVertex-yVertex(kPrev);
wPrev = uRay.*vPrev - vRay.*uPrev; % >=0 if vectorPrev is visible to the Apex
uNext = xVertex(kNext)-xVertex;
vNext = yVertex(kNext)-yVertex;
wNext = uRay.*vNext - vRay.*uNext; % >=0 if vectorNext is visible to the Apex
acute = uPrev.*vNext - vPrev.*uNext; % a corner must have an obtuse angle
isRightCorner = (wPrev>=tol_parallel)&(wNext<tol_parallel)&(acute<0);
isLeftCorner = (wPrev<tol_parallel)&(wNext>=tol_parallel)&(acute<0);
kRightCorner = find(isRightCorner);
kLeftCorner = find(isLeftCorner);
dRC_Apex = sqrt(uRay(isRightCorner).^2+vRay(isRightCorner).^2);
dLC_Apex = sqrt(uRay(isLeftCorner).^2+vRay(isLeftCorner).^2);
[~,RCsort] = sort(dRC_Apex);
[~,LCsort] = sort(dLC_Apex);
kRightCorner = kRightCorner(RCsort);
kLeftCorner = kLeftCorner(LCsort);

% --- OPTIMIZATION 2: Pre-allocation ---
% Pre-allocate arrays to avoid growing them inside loops.
% This is the most significant performance improvement.
% Each corner-to-intersection link can add at most 3 new vertices.
nMaxNewVertices = 3 * (length(kRightCorner) + length(kLeftCorner));
nOldVertices = length(xVertex);

% --- FIX 4: Only resize if new vertices are actually possible ---
if nMaxNewVertices > 0
    nTotalVertices = nOldVertices + nMaxNewVertices;

    % Resize arrays once by assigning to the last projected element
    xVertex(nTotalVertices, 1) = 0;
    yVertex(nTotalVertices, 1) = 0;
    kNext(nTotalVertices, 1) = 0;
    kPrev(nTotalVertices, 1) = 0;
    
    % Update vertex-to-boundary map for the *full* pre-allocated size
    % New vertices will inherit the boundary ID of the corner they came from.
    % This is a simplification; we only *really* need the map for kCorner.
    % The original map is sufficient.
    kVertexToBoundary(nTotalVertices, 1) = 0;
    kVertexToBoundary(1:nOldVertices) = kVertexToBoundary(1:nOldVertices);
else
    nTotalVertices = nOldVertices; % Unused, but explicit
end
% --- END FIX 4 ---


% Use a counter for the next available vertex slot
mVertex = nOldVertices + 1;
% --- END OPTIMIZATION 2 ---

isleapRC = false(size(kRightCorner)); % if corner and intersection point lie on various cycles of edge
isleapLC = false(size(kLeftCorner));
isPassedRC = false(size(kRightCorner)); % if the right corner is passed by sightline of apex to other right corner
isPassedLC = false(size(kLeftCorner)); % if the left corner is passed by sightline of apex to other left corner
kRCAdded = zeros(size(kRightCorner)); % index of the added vertex at the corner position when connecting corner to intersection point
kLCAdded = zeros(size(kLeftCorner));
%% link right corners to intersection points
tolmu = -10^-10; %tolerance for passing through a vertex

% --- FIX 2 & 3 ---
% Define ALL valid-range vectors *before* either loop runs,
% so they are not contaminated by new vertices.
kNext_valid = kNext(1:nOldVertices);
kPrev_valid = kPrev(1:nOldVertices);
xVertex_valid = xVertex(1:nOldVertices);
yVertex_valid = yVertex(1:nOldVertices);

% Pre-calculate vectors for Right Corner loop
uNext_R = xVertex_valid(kNext_valid) - xVertex_valid;
vNext_R = yVertex_valid(kNext_valid) - yVertex_valid;
dNext_R = sqrt(uNext_R.^2+vNext_R.^2);

% Pre-calculate vectors for Left Corner loop
uPrev_L = xVertex_valid(kPrev_valid) - xVertex_valid;
vPrev_L = yVertex_valid(kPrev_valid) - yVertex_valid;
dPrev_L = sqrt(uPrev_L.^2+vPrev_L.^2);
% --- END FIX 2 & 3 ---

for i = 1:length(kRightCorner)
    if ~isPassedRC(i)
        kCorner = kRightCorner(i);
        % Use O(1) lookup map instead of O(logN) find
        kBoundary = kVertexToBoundary(kCorner); 
        
        uRay = xVertex(kCorner) - xApex;
        vRay = yVertex(kCorner) - yApex;
        
        % --- FIX for Pre-allocation Error ---
        % Use the vectors defined *before* the loop
        detIntersect = -uRay*vNext_R + vRay*uNext_R;
        detIntersect(abs(detIntersect)<abs(tol_parallel)) = 0;
        
        lambda = ( -(xVertex_valid-xApex).*vNext_R + (yVertex_valid-yApex).*uNext_R )./detIntersect;
        mu = ( uRay.*(yVertex_valid-yApex) - vRay.*(xVertex_valid-xApex) )./detIntersect;
        
        % Vectorized intersection check (on all vertices up to original size)
        isIntersect = (mu>=tolmu) & (mu<1+tolmu) & (lambda>1) & (dNext_R>0);
        % --- END FIX ---
        
        isIntersect([kApex,kPrev(kApex),kCorner,kPrev(kCorner)])=0;
        
        kIntersectIdx = find(isIntersect);
        if isempty(kIntersectIdx)
            lambdaIntersect = [];
        else
            lambdaIntersect = lambda(kIntersectIdx);
        end
        
        if dispflag
            plot(xVertex(kCorner),yVertex(kCorner),'k>','markerfacecolor','k','markersize',3)
        end
        
        if ~isempty(lambdaIntersect)
            [lambdaMin,jMin] = min(lambdaIntersect);
            kMin = kIntersectIdx(jMin);
            alignedRCi = find(kRightCorner==kMin);
            while ~isempty(alignedRCi) % intersection point is another right corner
                kIntersectIdx(jMin)=[];
                lambdaIntersect(jMin)=[];
                isPassedRC(alignedRCi) = 1;
                if isempty(lambdaIntersect)
                    kMin = [];
                    break;
                end
                [lambdaMin,jMin] = min(lambdaIntersect);
                kMin = kIntersectIdx(jMin);
                alignedRCi = find(kRightCorner==kMin);
            end
            
            if isempty(kMin)
                continue; % All intersections were other corners
            end

            % Use O(1) lookup map
            isleapRC(i) = kBoundary ~= kVertexToBoundary(kMin);
            
            xIntersect = xApex + lambdaMin*uRay;
            yIntersect = yApex + lambdaMin*vRay;
            
            if dispflag
                plot([xVertex(kCorner),xIntersect],[yVertex(kCorner),yIntersect],'k:')
                plot(xIntersect,yIntersect,'k>','markerfacecolor','w','markersize',4)
            end
            
            % --- OPTIMIZATION: Use pre-allocated slots ---
            kIntersect = mVertex; % Get next available slot
            kRCAdded(i) = mVertex + 1;
            
            if abs(mu(kMin))<=abs(tolmu) % intersection point is a vertex (2 new vertices)
                xVertex(kIntersect) = xIntersect;
                yVertex(kIntersect) = yIntersect;
                xVertex(kIntersect+1) = xVertex(kCorner);
                yVertex(kIntersect+1) = yVertex(kCorner);
                
                kPrev(kNext(kCorner)) = kIntersect + 1;
                kNext(kIntersect+1) = kNext(kCorner);
                kPrev(kIntersect+1) = kMin;
                
                kNext(kCorner) = kIntersect;
                kPrev(kIntersect) = kCorner;
                kNext(kIntersect) = kNext(kMin);
                kPrev(kNext(kMin)) = kIntersect;
                
                kNext(kMin) = kIntersect + 1;
                
                mVertex = mVertex + 2; % Increment counter
            else % intersection point lies on an edge (3 new vertices)
                xVertex(kIntersect) = xIntersect;
                yVertex(kIntersect) = yIntersect;
                xVertex(kIntersect+1) = xVertex(kCorner);
                yVertex(kIntersect+1) = yVertex(kCorner);
                xVertex(kIntersect+2) = xIntersect;
                yVertex(kIntersect+2) = yIntersect;
                
                kPrev(kNext(kCorner)) = kIntersect + 1;
                kNext(kIntersect+1) = kNext(kCorner);
                kPrev(kIntersect+1) = kIntersect + 2;
                
                kNext(kCorner) = kIntersect;
                kPrev(kIntersect) = kCorner;
                kNext(kIntersect) = kNext(kMin);
                kPrev(kNext(kMin)) = kIntersect;
                
                kNext(kMin) = kIntersect + 2;
                kPrev(kIntersect+2) = kMin;
                kNext(kIntersect+2) = kIntersect + 1;
                
                mVertex = mVertex + 3; % Increment counter
            end
            % --- END OPTIMIZATION ---
        end
    end
end
%% link left corners to intersection points

% --- FIX 2 ---
% (Definitions moved before the right-corner loop)
% --- END FIX 2 ---

for i = 1:length(kLeftCorner)
    if ~isPassedLC(i)
        kCorner = kLeftCorner(i);
        % Use O(1) lookup map
        kBoundary = kVertexToBoundary(kCorner);
        
        uRay = xVertex(kCorner) - xApex;
        vRay = yVertex(kCorner) - yApex;
        
        % --- FIX for Pre-allocation Error ---
        % Use the vectors defined *before* the loop
        detIntersect = -uRay*vPrev_L + vRay*uPrev_L;
        detIntersect(abs(detIntersect)<10^-10) = 0;

        lambda = ( -(xVertex_valid-xApex).*vPrev_L + (yVertex_valid-yApex).*uPrev_L )./detIntersect;
        mu = ( uRay.*(yVertex_valid-yApex) - vRay.*(xVertex_valid-xApex) )./detIntersect;
        
        % Vectorized intersection check (on all vertices up to original size)
        isIntersect = (mu>=tolmu) & (mu<1+tolmu) & (lambda>1) & (dPrev_L>0);
        % --- END FIX ---

        isIntersect([kApex,kNext(kApex),kCorner,kNext(kCorner)])=0;

        kIntersectIdx = find(isIntersect);
        if isempty(kIntersectIdx)
            lambdaIntersect = [];
        else
            lambdaIntersect = lambda(kIntersectIdx);
        end

        if dispflag
            plot(xVertex(kCorner),yVertex(kCorner),'k<','markerfacecolor','k','markersize',3)
        end
        
        if ~isempty(lambdaIntersect)
            [lambdaMin,jMin] = min(lambdaIntersect);
            kMin = kIntersectIdx(jMin);
            alignedLCi = find(kLeftCorner==kMin);
            while ~isempty(alignedLCi) % intersection point is another left corner
                kIntersectIdx(jMin)=[];
                lambdaIntersect(jMin)=[];
                isPassedLC(alignedLCi) = 1;
                if isempty(lambdaIntersect)
                    kMin = [];
                    break;
                end
                [lambdaMin,jMin] = min(lambdaIntersect);
                kMin = kIntersectIdx(jMin);
                alignedLCi = find(kLeftCorner==kMin);
            end
            
            if isempty(kMin)
                continue; % All intersections were other corners
            end
            
            % Use O(1) lookup map
            isleapLC(i) = kBoundary ~= kVertexToBoundary(kMin);
            
            xIntersect = xApex + lambdaMin*uRay;
            yIntersect = yApex + lambdaMin*vRay;
            
            if dispflag
                plot([xVertex(kCorner),xIntersect],[yVertex(kCorner),yIntersect],'k:')
                plot(xIntersect,yIntersect,'k<','markerfacecolor','w','markersize',4)
            end
            
            % --- OPTIMIZATION: Use pre-allocated slots ---
            kIntersect = mVertex; % Get next available slot
            kLCAdded(i) = mVertex + 1;
            
            if abs(mu(kMin))<=abs(tolmu) % intersection point is a vertex (2 new vertices)
                xVertex(kIntersect) = xIntersect;
                yVertex(kIntersect) = yIntersect;
                xVertex(kIntersect+1) = xVertex(kCorner);
                yVertex(kIntersect+1) = yVertex(kCorner);
                
                kNext(kPrev(kCorner)) = kIntersect + 1;
                kPrev(kIntersect+1) = kPrev(kCorner);
                kNext(kIntersect+1) = kMin;
                
                kPrev(kCorner) = kIntersect;
                kNext(kIntersect) = kCorner;
                kNext(kPrev(kMin)) = kIntersect;
                kPrev(kIntersect) = kPrev(kMin);
                
                kPrev(kMin) = kIntersect + 1;
                
                mVertex = mVertex + 2; % Increment counter
            else % intersection point lies on an edge (3 new vertices)
                xVertex(kIntersect) = xIntersect;
                yVertex(kIntersect) = yIntersect;
                xVertex(kIntersect+1) = xVertex(kCorner);
                yVertex(kIntersect+1) = yVertex(kCorner);
                xVertex(kIntersect+2) = xIntersect;
                yVertex(kIntersect+2) = yIntersect;

                kNext(kPrev(kCorner)) = kIntersect + 1;
                kPrev(kIntersect+1) = kPrev(kCorner);
                kNext(kIntersect+1) = kIntersect + 2;
                
                kPrev(kCorner) = kIntersect;
                kNext(kIntersect) = kCorner;
                kNext(kPrev(kMin)) = kIntersect;
                kPrev(kIntersect) = kPrev(kMin);
                
                kPrev(kMin) = kIntersect + 2;
                kNext(kIntersect+2) = kMin;
                kPrev(kIntersect+2) = kIntersect + 1;
                
                mVertex = mVertex + 3; % Increment counter
            end
            % --- END OPTIMIZATION ---
        end
    end
end

% --- OPTIMIZATION 3: Trim pre-allocated arrays ---
% Remove the unused, pre-allocated space
% --- FIX 4: Only trim if we actually resized ---
if nMaxNewVertices > 0
    xVertex(mVertex:end) = [];
    yVertex(mVertex:end) = [];
    kNext(mVertex:end) = [];
    kPrev(mVertex:end) = [];
end
% --- END FIX 4 & OPTIMIZATION 3 ---

%% obtain vertices of the visibility polygon
% loop around starting from Apex:
kVisi = kApex;
while kNext(kVisi(end))~=kApex
    kVisi = [kVisi; kNext(kVisi(end))];
end
kVisi = [kVisi; kNext(kVisi(end))];
xVisi = xVertex(kVisi);
yVisi = yVertex(kVisi);
%% obtain effective corners (children apexes)
% if 0
xRightCorner = xVertex(kRightCorner);
xLeftCorner = xVertex(kLeftCorner);
yRightCorner = yVertex(kRightCorner);
yLeftCorner = yVertex(kLeftCorner);
% find corners on the visibility polygon:
isRCOnVisi = inpolygon(xRightCorner,yRightCorner,xVisi,yVisi);
isLCOnVisi = inpolygon(xLeftCorner,yLeftCorner,xVisi,yVisi);
% remove corners with small influences:
kPendingRC = kRCAdded(~isleapRC & isRCOnVisi & kRCAdded~=0 & ~isPassedRC);
kPendingLC = kLCAdded(~isleapLC & isLCOnVisi & kLCAdded~=0 & ~isPassedLC);
for i = 1:length(kPendingRC)
    kRemovedPolygon = kPendingRC(i);
    while kNext(kRemovedPolygon(end))~=kPendingRC(i)
        kRemovedPolygon = [kRemovedPolygon; kNext(kRemovedPolygon(end))];
    end
    %     plot(xVertex(kRemovedPolygon),yVertex(kRemovedPolygon),'c-')
    if (length(kRemovedPolygon)<=5 && polyarea(xVertex(kRemovedPolygon),yVertex(kRemovedPolygon))<=threshold^2/50)
        kPendingRC(i) = nan;
    end
end
for i = 1:length(kPendingLC)
    kRemovedPolygon = kPendingLC(i);
    while kNext(kRemovedPolygon(end))~=kPendingLC(i)
        kRemovedPolygon = [kRemovedPolygon; kNext(kRemovedPolygon(end))];
    end
    %     plot(xVertex(kRemovedPolygon),yVertex(kRemovedPolygon),'c-')
    if (length(kRemovedPolygon)<=5 && polyarea(xVertex(kRemovedPolygon),yVertex(kRemovedPolygon))<=threshold^2/50)
        kPendingLC(i) = nan;
    end
end
xEffCorner = [xVertex(kRightCorner(isleapRC & isRCOnVisi & ~isPassedRC));xVertex(kPendingRC(~isnan(kPendingRC)));...
    xVertex(kLeftCorner(isleapLC & isLCOnVisi & ~isPassedLC));xVertex(kPendingLC(~isnan(kPendingLC)))];
yEffCorner = [yVertex(kRightCorner(isleapRC & isRCOnVisi & ~isPassedRC));yVertex(kPendingRC(~isnan(kPendingRC)));...
    yVertex(kLeftCorner(isleapLC & isLCOnVisi & ~isPassedLC));yVertex(kPendingLC(~isnan(kPendingLC)))];
if dispflag
    plot(xEffCorner,yEffCorner,'r^')
    plot(xVisi,yVisi,'r-')
end
% else
%     xEffCorner = [xVertex(kRightCorner);xVertex(kLeftCorner)];
%     yEffCorner = [yVertex(kRightCorner);yVertex(kLeftCorner)];
% end
%% validation
% uRay = xVisi-xApex;
% vRay = yVisi-yApex;
% uNext = xVisi([2:end,1]')-xVisi;
% vNext = yVisi([2:end,1]')-yVisi;
% wNext = uRay.*vNext - vRay.*uNext;
% if ~isempty(find(wNext<-10^-4,1))
%     wNext(wNext<-10^-4)
%     disp('error in visibility polygon')
%     save errordata xBoundary yBoundary xApex yApex threshold dispflag
% end
end




