function [fitting_s_z, bin_pts_mid_s, Q2, ss_in, z_in, ss_up, z_up, ss_do, z_do, s_T, L, S, P, RMSE, dimensionless_drop, slope, numberPoints_in_bin] = process_s_z_relationship_return_info_width(sMap, zMap, bin_size, ds, outlength, pltFlag, options)

arguments
    sMap
    zMap
    bin_size
    ds
    outlength
    pltFlag
    options.medianFilter = true;
end

% Flatten the input matrices
s_flaten = sMap(:);
z_flaten = zMap(:);
nan_index = isnan(s_flaten) | isnan(z_flaten);
s_flaten(nan_index) = [];
z_flaten(nan_index) = [];

if options.medianFilter

    % Define bins for s values
    bin_pts_s = min(s_flaten):bin_size:max(s_flaten);
    bin_pts_mid_s = [];
    Q2 = [];
    numberPoints_in_bin = [];
    % Calculate the median z-value for each s bin, removing outliers
    for i = 1:length(bin_pts_s) - 1
        logic = s_flaten > bin_pts_s(i) & s_flaten < bin_pts_s(i + 1);
        Q2(end+1) = quantile(z_flaten(logic), 0.5);
        bin_pts_mid_s(end+1) = (bin_pts_s(i) + bin_pts_s(i + 1)) / 2;
        numberPoints_in_bin(end+1) = numel(z_flaten(logic));
    end

    % Remove NaN values
    bin_pts_mid_s(isnan(Q2)) = [];
    Q2(isnan(Q2)) = [];

    % Polynomial fitting and extrapolation
    s_T = max(s_flaten);
    ss_in = 0:ds:s_T;
    
    if ~isempty(bin_pts_mid_s)
        p = polyfit(bin_pts_mid_s, Q2, 2);
        L = p(1);
        S = p(2);
        P = p(3);
        z_flaten_fit = polyval(p, s_flaten);
        RMSE = sqrt(mean((z_flaten_fit - z_flaten).^2));
        % Extrapolate for the upper and lower ranges
        ss_up = -outlength:ds:-ds;
        z_up = S * ss_up + P;
        ss_do = ss_in(end) + ds + (0:ds:outlength);
        z_do = (2 * L * s_T + S) * ss_do - L * s_T^2 + P;
        
        % Calculate RMSE
        z_in = polyval(p, ss_in);
        if pltFlag
        % Plot the data
        plot(bin_pts_mid_s, Q2, 'b.')
        % Plot the fit and extrapolations
        plot(ss_in, z_in, 'b-');
        plot(ss_up, z_up, 'b--');
        plot(ss_do, z_do, 'b--');
        plot([0 s_T], z_in([1 end]), 'bo', 'MarkerSize', 6);
        axis padded
        daspect([5 1 1])
        grid on
        box on
        end
        ss = [ss_up ss_in ss_do];
        zz = [z_up z_in z_do];
        fitting_s_z = [ss' zz'];
    end
else
    bin_pts_mid_s = [];
    Q2 = [];
    s_T = max(s_flaten);
    ss_in = 0:ds:s_T;

    if ~isempty(s_flaten)
        p = polyfit(s_flaten, z_flaten, 2);
        L = p(1);
        S = p(2);
        P = p(3);
        z_in = polyval(p, ss_in);
        z_flaten_fit = polyval(p, s_flaten);
        RMSE = sqrt(mean((z_flaten_fit - z_flaten).^2));
        % Extrapolate for the upper and lower ranges
        ss_up = -outlength:ds:-ds;
        z_up = S * ss_up + P;
        ss_do = ss_in(end) + ds + (0:ds:outlength);
        z_do = (2 * L * s_T + S) * ss_do - L * s_T^2 + P;
        
        if pltFlag
        % Plot the fit and extrapolations
        plot(ss_in, z_in, 'b-');
        plot(ss_up, z_up, 'b--');
        plot(ss_do, z_do, 'b--');
        plot([0 s_T], z_in([1 end]), 'bo', 'MarkerSize', 6);
        axis padded
        daspect([5 1 1])
        grid on
        box on
        end
        ss = [ss_up ss_in ss_do];
        zz = [z_up z_in z_do];
        fitting_s_z = [ss' zz'];
    end
end

dimensionless_drop = L*s_T/4;
slope = -L*s_T-S;

end
