function plotFanTopoResults(xMesh, yMesh, zTopo, zMesh, xApex, yApex, fanBoundarySHP, contourInterval, clim_value)

    if nargin < 7 || isempty(fanBoundarySHP)
        fanBoundarySHP = nan;
    end
    if nargin < 8 || isempty(contourInterval)
        contourInterval = 10;
    end


    % Fill NaNs in zTopo using zMesh
    zTopoFill = zTopo;
    zTopoFill(isnan(zTopo)) = zMesh(isnan(zTopo));

    % Plotting
    figure
    lightterrain2D_imagesc(xMesh, yMesh, zTopoFill)
    hold on
    freezeColors

    if ~isnan(fanBoundarySHP)
        % Read fan boundary shapefile
        fan_boundary = shaperead(fanBoundarySHP);
        fan_boundary_x = fan_boundary.X(1:end-1);
        fan_boundary_y = fan_boundary.Y(1:end-1);
        volumeCalculateExtentXY = [fan_boundary_x', fan_boundary_y'];
        inMask = inpolygon(xMesh,yMesh,volumeCalculateExtentXY(:,1), volumeCalculateExtentXY(:,2));
    else
        inMask = true(size(zTopo));
    end

    % Calculate difference between zTopo and zMesh
    zDiff = zTopo - zMesh;
    zDiff(~inMask) = nan;

    % Plot color difference map
    pcolor(xMesh, yMesh, zDiff)
    shading flat
    colormap("turbo")

    % Calculate fan volume
    fanVolume = sum(zDiff, 'all', 'omitnan') * (xMesh(1,2) - xMesh(1,1)).^2;

    if ~isnan(fanBoundarySHP)
        % Plot volume calculation extent if provided
        plot(volumeCalculateExtentXY(:,1), volumeCalculateExtentXY(:,2), 'r-')
    end
    
    % Plot apex point
    plot(xApex, yApex, 'r.', 'MarkerSize', 8)

    % Set color limits and colorbar
    % clim([min(clim_value), max(clim_value)])
    c = colorbar;
    ylabel(c,'\Delta z (m)')

    % Plot contour
    contour(xMesh, yMesh, zTopoFill, min(zTopoFill(:)):contourInterval:max(zTopoFill(:)), 'k', 'LineWidth',0.1)
    axis equal

    if ~isnan(fanBoundarySHP)
        % Add title
        title(['sim volume = ' num2str(fanVolume, '%.0f') ' (m^3) in boundary'])
    else
        title(['sim volume = ' num2str(fanVolume, '%.0f') ' (m^3)'])
    end
    xlabel('Easting (m)')
    ylabel('Northing (m)')
    hold off
end