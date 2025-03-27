function filledImg = scanlineFillWithTopPatch(xMesh, yMesh, polygon)
    %SCANLINEFILLWITHTOPPATCH Fill a polygon using a scanline algorithm
    %   and then manually patch horizontal edges.
    %
    %   INPUTS:
    %       xMesh   - 2D array (or vector) defining the X-coordinates of the mesh.
    %                 In this code, xMesh(1,:) is used as the horizontal axis.
    %       yMesh   - 2D array (or vector) defining the Y-coordinates of the mesh.
    %                 Here, yMesh(:,1)' is treated as the vertical axis.
    %       polygon - Nx2 matrix of real-world polygon vertices specified as (x,y).
    %
    %   The function maps the polygon's (x,y) coordinates to discrete pixel indices
    %   (imgWidth, imgHeight) via 'interp1'. Then, it performs a standard scanline
    %   fill (using "lower inclusive, upper exclusive" to avoid double-counting).
    %   Finally, it manually patches any horizontal edges (y1 == y2) to ensure they
    %   are filled as well.
    %
    %   OUTPUT:
    %       filledImg - A logical array of size [imgHeight x imgWidth], where 'true'
    %                   marks the filled polygon area and 'false' is the background.
    %
    %   Example usage:
    %       >> filledImg = scanlineFillWithTopPatch(xMesh, yMesh, polygon);
    %
    %   Note:
    %   - xMesh(1,:) defines the grid's X-axis. We compute imgWidth = length(x_interp).
    %   - yMesh(:,1)' defines the grid's Y-axis. We compute imgHeight = length(y_interp).
    %   - The polygon is first interpolated into discrete image coordinates,
    %     then filled row by row. Horizontal edges are added in a final pass.
    x_interp = xMesh(1,:);
    y_interp = yMesh(:,1)';
    imgWidth = length(x_interp);
    imgHeight = length(y_interp);
    polygonX = interp1(x_interp, 1:length(x_interp),  polygon(:,1),'linear','extrap');
    polygonY = interp1(y_interp, 1:length(y_interp), polygon(:,2),'linear','extrap');
    polygon = [polygonX polygonY];

    % Ensure polygon is closed
    if ~isequal(polygon(1,:), polygon(end,:))
        polygon(end+1,:) = polygon(1,:);
    end

    % 1) STANDARD SCANLINE FILL (Lower Inclusive, Upper Exclusive)

    % Create empty logical image
    filledImg = false(imgHeight, imgWidth);

    % Determine bounding box in Y
    yMin = floor(min(polygon(:,2)));
    yMax = ceil(max(polygon(:,2)));
    yMin = max(yMin, 1);
    yMax = min(yMax, imgHeight);

    % Number of edges (last vertex repeated, so -1)
    numEdges = size(polygon,1) - 1;

    for y = yMin : yMax

        xIntersections = [];

        for e = 1:numEdges
            x1 = polygon(e,1);  y1 = polygon(e,2);
            x2 = polygon(e+1,1);y2 = polygon(e+1,2);

            % Order the edge so that y1 <= y2
            if y2 < y1
                [x1,x2] = deal(x2, x1);
                [y1,y2] = deal(y2, y1);
            end

            % Check if this horizontal scan line intersects edge
            % "Lower inclusive, upper exclusive" => y1 <= y < y2
            if (y >= y1) && (y < y2)
                % Edge is not horizontal => safe to compute intersection
                if (y1 ~= y2)
                    dx = (x2 - x1);
                    dy = (y2 - y1);
                    intersectX = x1 + dx * (y - y1) / dy;
                    xIntersections(end+1) = intersectX; %#ok<AGROW>
                end
            end
        end

        % Sort intersections
        xIntersections = sort(xIntersections);

        % Fill between pairs of intersections
        for i = 1 : 2 : length(xIntersections)
            if (i+1 <= length(xIntersections))
                xStart = ceil(xIntersections(i));
                xEnd   = floor(xIntersections(i+1));

                % Clip
                xStart = max(xStart, 1);
                xEnd   = min(xEnd, imgWidth);

                if xStart <= xEnd
                    filledImg(y, xStart:xEnd) = true;
                end
            end
        end
    end

    % 2) MANUAL PATCH FOR HORIZONTAL EDGES
    % If you want horizontal edges to be filled, do it here.

    for e = 1:numEdges
        x1 = polygon(e,1);  y1 = polygon(e,2);
        x2 = polygon(e+1,1);y2 = polygon(e+1,2);

        % Check if the edge is horizontal: y1 == y2
        if abs(y2 - y1) < 1e-13    % or just (y1 == y2) if integer coords
            % We want to fill that entire horizontal span
            yH = round(y1);  % same as y2
            if (yH >= 1) && (yH <= imgHeight)

                % Determine left and right x
                if x2 < x1
                    xLeft = x2; 
                    xRight = x1;
                else
                    xLeft = x1;
                    xRight = x2;
                end

                xStart = max(1, ceil(xLeft));
                xEnd   = min(imgWidth, floor(xRight));

                if xStart <= xEnd
                    filledImg(yH, xStart:xEnd) = true;
                end
            end
        end
    end

    for i_vertex = 1:length(polygon(:,2))
        if abs(polygon(i_vertex,1) - round(polygon(i_vertex,1))) < 1e-13 && abs(polygon(i_vertex,2) - round(polygon(i_vertex,2))) < 1e-13
            pointXY = polygon(i_vertex,:);
            filledImg(round(pointXY(2)), round(pointXY(1)))= true;
        end
    end

end