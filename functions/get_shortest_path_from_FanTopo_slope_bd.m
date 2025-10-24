function semiApexsAll = get_shortest_path_from_FanTopo_slope_bd(xMesh,yMesh,sigmaTopo,kTopo_sigma,xyzkApex_sigma,xyzVisi_sigma,ex_pt,pltflag)
%---Output---
% semiApexsAll: The shorest path from the given points saved in cell format. 
% In each cell, it is the path x,y,z save in n-3 array 

%---Input---
% xMesh,yMesh,sigmaTopo,kTopo_sigma,xyzkApex_sigma,xyzVisi_sigma: The information get 
% from FanTopo_slope_bd function.

% ex_pt: Given start points for finding the shortest path.


    k_apex_useful = unique(kTopo_sigma(~isnan(kTopo_sigma)));
    
    xyzkApex_useful = xyzkApex_sigma(k_apex_useful,1:3);
    
    semiApexsAll = {};
    ex_pt_k = interp2(xMesh,yMesh,kTopo_sigma,ex_pt(:,1),ex_pt(:,2));
    ex_pt_z = interp2(xMesh,yMesh,sigmaTopo,ex_pt(:,1),ex_pt(:,2));
    for i = 1:length(ex_pt_k)
       ex_pt_k(i) = find(k_apex_useful==ex_pt_k(i));
    end
    
    
    
    for ii = 1:length(ex_pt_k)
        semi_x = xyzkApex_useful(ex_pt_k(ii),1);
        semi_y = xyzkApex_useful(ex_pt_k(ii),2);
        semi_z = xyzkApex_useful(ex_pt_k(ii),3);
        semiApexs = [ex_pt(ii,1), ex_pt(ii,2), ex_pt_z(ii,1);...
                     semi_x, semi_y, semi_z];
    
        while 1
        dd_from_kApex = sum((repmat([semi_x semi_y],length(xyzkApex_useful),1) - xyzkApex_useful(:,1:2)).^2,2).^0.5;
        dz_from_kApex = repmat(semi_z,length(xyzkApex_useful),1) - xyzkApex_useful(:,3);
        [~,slope_from_kApex_index] = sort(dz_from_kApex./dd_from_kApex);
        
        i = 1;
        while i<=length(slope_from_kApex_index)
            useful_index = slope_from_kApex_index(i);
            main_index = k_apex_useful(useful_index);
            in = inpolygon(semi_x,semi_y,xyzVisi_sigma{main_index}(:,1),xyzVisi_sigma{main_index}(:,2));
            if in
                semiApexs(end+1,1:3) = xyzkApex_sigma(main_index,1:3);
                semi_x = semiApexs(end,1);
                semi_y = semiApexs(end,2);
                semi_z = semiApexs(end,3);
                break
            else
                i = i + 1;
            end
        
        end
        if main_index == 1
            semiApexsAll{end+1} = semiApexs;
            break
        end
        end
    end
    
    if pltflag
        figure
        pcolor(xMesh,yMesh,kTopo_sigma)
        shading flat
        axis equal
        hold on
        plot(xyzkApex_sigma(:,1),xyzkApex_sigma(:,2),'r.')
        for i = 1:length(xyzkApex_sigma)
            %text(xyzkApex_sigma(i,1),xyzkApex_sigma(i,2),[num2str(i) ',' num2str(xyzkApex_sigma(i,4))],'fontsize',8)
        end
        
        k_apex_useful = unique(kTopo_sigma(~isnan(kTopo_sigma)));
        plot(xyzkApex_sigma(k_apex_useful,1),xyzkApex_sigma(k_apex_useful,2),'go')
        for i = 1:length(ex_pt_k)
            plot(semiApexsAll{i}(:,1),semiApexsAll{i}(:,2),'r-')
        end
        plot(ex_pt(:,1),ex_pt(:,2),'ro')
    end

end