function [zTopo_sim, heightAG_Volume_All] = reconstruct_fan_surface(xMesh, yMesh, zMesh, xApex, yApex, volume_expected, guessHeightAboveGround_top, guessHeightAboveGround_bottom, dz_profile, options)

arguments
    xMesh double
    yMesh double
    zMesh double
    xApex double
    yApex double
    volume_expected double
    guessHeightAboveGround_top
    guessHeightAboveGround_bottom
    dz_profile
    options.fanBoundarySHP = nan;
    options.tol = 0.03
end

fan_boundary_shp = options.fanBoundarySHP;

% Handle fan boundary if provided
if ~isnan(fan_boundary_shp)
    fan_boundary = shaperead(fan_boundary_shp);
    fan_boundary_x = fan_boundary.X(1:end-1);
    fan_boundary_y = fan_boundary.Y(1:end-1);
    volumeCalculateExtentXY = [fan_boundary_x', fan_boundary_y'];
    volumeCalculateMask = inpolygon(xMesh, yMesh, fan_boundary_x, fan_boundary_y);
else
    volumeCalculateExtentXY = nan;
    volumeCalculateMask = true(size(zMesh));
end

% Simulation Process
heights = [guessHeightAboveGround_bottom, guessHeightAboveGround_top];
minMaxInitialGuessHeightVolume = zeros(2, 2);




parfor i = 1:2
    guessHeight = heights(i);
    zApex = interp2(xMesh, yMesh, zMesh, xApex, yApex) + guessHeight;
    [zTopo, ~, ~] = FanTopo(xMesh, yMesh, zMesh, xApex, yApex, zApex, 'caseName', 'myProfile', 'dz_interpM', {dz_profile});
    
    dod = zTopo - zMesh;
    dod(~volumeCalculateMask) = nan;

    fanVolume = sum(dod, 'all', 'omitnan') * (xMesh(1,2) - xMesh(1,1))^2;
    
    minMaxInitialGuessHeightVolume(i, :) = [guessHeight, fanVolume];
    if any(~volumeCalculateMask(:))
        fprintf('HAG = %.2f [L], Volume = %.2f [L^3], within given boundary\n', guessHeight, fanVolume);
    else
        fprintf('HAG = %.2f [L], Volume = %.2f [L^3], within simulation area\n', guessHeight, fanVolume);
    end
end

[zTopo_sim, heightAG_Volume_All] = fanTopoSimVolumeMask(volume_expected,minMaxInitialGuessHeightVolume,xMesh,yMesh,zMesh,xApex,yApex,volumeCalculateMask,dz_profile,options.tol);

end