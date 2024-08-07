function zCone = coneFunction(zApex,D,options)
    arguments
        zApex double
        D double
        options.caseName = 'cone'
        options.tanAlpha = nan
        options.K = nan
        options.zApex0 = nan
        options.tanInfinite = nan
        options.dz_interp = nan

    end
    switch options.caseName
        case 'cone'
            zCone = zApex - options.tanAlpha*D; % cone surface
        case 'concavity'
            S = options.tanAlpha-options.K*(options.zApex0-zApex);
            zCone = zApex-S/options.K*(1-exp(-options.K*D));
        case 'infinite'
            S = options.tanAlpha-options.K*(options.zApex0-zApex);
            zCone = zApex-(S-options.tanInfinite)/options.K*(1-exp(-options.K*D))-options.tanInfinite*D;
        case 'myProfile'
            out_linear = 1;
            if ~out_linear
                dOffset = interp1_extrap(options.dz_interp(:,2),options.dz_interp(:,1),zApex);
            else
                dOffset = interp1(options.dz_interp(:,2),options.dz_interp(:,1),zApex,'linear','extrap');
            end
            if numel(zApex) ==1
                if ~out_linear
                    zCone = interp1_extrap(options.dz_interp(:,1)-dOffset,options.dz_interp(:,2),D(:));
                else
                    zCone = interp1(options.dz_interp(:,1)-dOffset,options.dz_interp(:,2),D(:),'linear','extrap');                
                end
                zCone = reshape(zCone,size(D));
            else
                zCone = nan(size(zApex));
                for i = 1:length(zCone)
                    if isnan(zApex(i))
                        zCone(i) = nan;
                    else
                        if ~out_linear
                            zCone(i) = interp1_extrap(options.dz_interp(:,1)-dOffset(i),options.dz_interp(:,2),D(:));
                        else
                            zCone(i) = interp1(options.dz_interp(:,1)-dOffset(i),options.dz_interp(:,2),D(:),'linear','extrap');
                        end
                    end
                end
            end
    end
end

function vq = interp1_extrap(x,v,xq)
    [~,minIndex] = min(x);
    [~,maxIndex] = max(x);
    xq_small = xq<x(minIndex);
    xq_large = xq>x(maxIndex);
    vq = interp1(x,v,xq);
    vq(xq_small) = v(minIndex);
    vq(xq_large) = v(maxIndex);
end
