function fitting_s_z = process_s_z_relationship(sMap, zMap, bin_size, ds, outlength, pltFlag)

    % Flatten the input matrices
    s_flaten = sMap(:);
    z_flaten = zMap(:);

    % Define bins for s values
    bin_pts_s = min(s_flaten):bin_size:max(s_flaten);
    bin_pts_mid_s = [];
    Q2 = [];

    % Calculate the median z-value for each s bin, removing outliers
    for i = 1:length(bin_pts_s) - 1
        logic = s_flaten > bin_pts_s(i) & s_flaten < bin_pts_s(i + 1);
        Q2(end+1) = quantile(z_flaten(logic), 0.5);
        bin_pts_mid_s(end+1) = (bin_pts_s(i) + bin_pts_s(i + 1)) / 2;
    end

    % Remove NaN values
    bin_pts_mid_s(isnan(Q2)) = [];
    Q2(isnan(Q2)) = [];

    % Polynomial fitting and extrapolation
    dd_max = max(bin_pts_mid_s);
    dd_in = 0:ds:dd_max;
    
    if ~isempty(bin_pts_mid_s)
        p = polyfit(bin_pts_mid_s, Q2, 2);
        a = p(1);
        b = p(2);
        c = p(3);
        z_in = polyval(p, dd_in);

        % Extrapolate for the upper and lower ranges
        dd_up = -outlength:ds:-ds;
        z_up = b * dd_up + c;
        dd_do = dd_in(end) + ds + (0:ds:outlength);
        z_do = (2 * a * dd_max + b) * dd_do - a * dd_max^2 + c;
        
        if pltFlag
        % Plot the data
        plot(bin_pts_mid_s, Q2, 'b.')
        % Plot the fit and extrapolations
        plot(dd_in, z_in, 'b-');
        plot(dd_up, z_up, 'b--');
        plot(dd_do, z_do, 'b--');
        plot([0 dd_max], z_in([1 end]), 'bo', 'MarkerSize', 6);
        axis padded
        daspect([5 1 1])
        grid on
        box on
        end
        ss = [dd_up dd_in dd_do];
        zz = [z_up z_in z_do];
        fitting_s_z = [ss' zz'];
    end

end
