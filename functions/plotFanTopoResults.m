function plotFanTopoResults(xMesh, yMesh, zTopo, zMesh, xApex, yApex, fanBoundarySHP, contourInterval)

    if nargin < 7 || isempty(fanBoundarySHP)
        fanBoundarySHP = nan;
    end
    if nargin < 8 || isempty(contourInterval)
        contourInterval = 10;
    end


    % Fill NaNs in zTopo using zMesh
    zTopoFill = zTopo;
    zTopoFill(isnan(zTopo)) = zMesh(isnan(zTopo));

    % Calculate difference between zTopo and zMesh
    zDiff = zTopo - zMesh;

    % Plotting
    figure
    lightterrain2D_imagesc(xMesh, yMesh, zTopoFill)
    hold on
    freezeColors

    % Plot color difference map
    pcolor(xMesh, yMesh, zDiff)
    shading flat

    if ~isnan(fanBoundarySHP)
        % Read fan boundary shapefile
        fan_boundary = shaperead(fanBoundarySHP);
        fan_boundary_x = fan_boundary.X(1:end-1);
        fan_boundary_y = fan_boundary.Y(1:end-1);
        volumeCalculateExtentXY = [fan_boundary_x', fan_boundary_y'];
    
        % Plot volume calculation extent if provided
        plot(volumeCalculateExtentXY(:,1), volumeCalculateExtentXY(:,2), 'r-')

        inMask = inpolygon(xMesh,yMesh,volumeCalculateExtentXY(:,1), volumeCalculateExtentXY(:,2));
    else
        inMask = true(size(zTopo));
    end

    dzMesh = zTopo - zMesh;
    dzMesh(~inMask) = nan;
    % Calculate fan volume
    fanVolume = sum(dzMesh, 'all', 'omitnan') * (xMesh(1,2) - xMesh(1,1)).^2;

    colormap("turbo")
    % Plot apex point
    plot(xApex, yApex, 'r.', 'MarkerSize', 8)

    % Set color limits and colorbar
    clim([min(zDiff, [], 'all'), max(zDiff, [], 'all')])
    c = colorbar;
    ylabel(c,'\Delta z')

    % Plot contour
    contour(xMesh, yMesh, zTopoFill, min(zTopoFill(:)):contourInterval:max(zTopoFill(:)), 'k')
    axis equal

    if ~isnan(fanBoundarySHP)
        % Add title
        title(['sim volume = ' num2str(fanVolume, '%.0f') ' [L^3] in boundary'])
    else
        title(['sim volume = ' num2str(fanVolume, '%.0f') ' [L^3]'])
    end
    hold off
end