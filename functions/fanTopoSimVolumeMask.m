function [zTopo, heightAG_Volume_All] = fanTopoSimVolumeMask(interpV, minMaxInitialGuessHeightVolume, xMesh, yMesh, zMesh, xApex, yApex, volMask, dz_profile, tol)
    % Set default tolerance value if not provided
    if nargin < 10
        tol = 0.03;
    end

    % Initial error check for apex heights
    minMaxVolumeDiff = minMaxInitialGuessHeightVolume(:,2) - min(interpV);
    if all(minMaxVolumeDiff >= 0)
        error('Bottom apex too high');
    elseif all(minMaxInitialGuessHeightVolume(:,2) - max(interpV) >= 0)
        error('Top apex too low');
    end
    
    zApexGround = interp2(xMesh, yMesh, zMesh, xApex, yApex);

    % Use the previously calculated height above ground and volume relationship
    % to run the code
    if size(minMaxInitialGuessHeightVolume, 1) > 2
        initialFlag = false;
        heightAG_Volume_All = minMaxInitialGuessHeightVolume;
    else
        initialFlag = true;
    end

    while ~isempty(interpV)
        interpV(interpV < 0) = []; % Remove negative volumes
        if isempty(interpV), break; end
        
        if ~initialFlag
            [~, uniqueIdx] = unique(heightAG_Volume_All(:,2));
            heightAG_Volume_All = heightAG_Volume_All(uniqueIdx,:);
        else
            heightAG_Volume_All = minMaxInitialGuessHeightVolume;
            initialFlag = false;
        end

        interpHAG = interp1(heightAG_Volume_All(:,2), heightAG_Volume_All(:,1), interpV);
        
        for i = 1:length(interpHAG)
            % Run fan topo simulation process
            zApex = zApexGround + interpHAG(i);
            [zTopo,~,~,~,~,~] = FanTopo(xMesh, yMesh, zMesh, xApex, yApex, zApex, 'caseName', 'myProfile', 'dz_interpM', {dz_profile});
            
            zMesh_sim = zTopo;
            zMesh_sim(isnan(zMesh_sim)) = zMesh(isnan(zMesh_sim));
            dod = zTopo - zMesh;
            dod(~volMask) = nan;
            
            fanVolume = sum(dod, 'all', 'omitnan') * (xMesh(1,2) - xMesh(1,1))^2;
            if any(~volMask(:))
                fprintf('HAG = %.2f [L], Volume = %.2f [L^3], within given boundary\n', interpHAG(i), fanVolume);
            else
                fprintf('HAG = %.2f [L], Volume = %.2f [L^3], within simulation area\n', interpHAG(i), fanVolume);
            end
            heightAG_Volume_All = [heightAG_Volume_All; interpHAG(i), fanVolume];
            if abs(fanVolume - interpV(i)) <= interpV(i) * tol
                interpV(i) = -1;
            end
        end
    end
end
