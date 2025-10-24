function lightterrain2D_imagesc(xMap,yMap,zMesh)
[nx, ny, nz] = surfnorm(xMap,yMap,zMesh);
lightVector = [0,-1,-1];
Nx=reshape(nx,[size(nx,1)*size(nx,2),1]);
Ny=reshape(ny,[size(ny,1)*size(ny,2),1]);
Nz=reshape(nz,[size(nx,1)*size(nz,2),1]);
terrainN = [Nx,Ny,Nz];
lightVector = repmat(lightVector,length(terrainN),1);
angle = acos(dot((lightVector./vecnorm(lightVector,2,2))', (terrainN./vecnorm(terrainN,2,2))'));
angle = reshape(angle,[size(nx,1),size(nx,2)]);
angleIndex=pi:-pi/100:0;
colorIndex=1:-1/100:0;
color=interp1(angleIndex,colorIndex,angle);
%blackwhitemap=repmat(color',1,3);
C(:,:,1)=reshape(color,[size(nx,1),size(nx,2)]);
C(:,:,2)=C(:,:,1);
C(:,:,3)=C(:,:,1);
%surf(xMap,yMap,zMesh,C);
imagesc(xMap(1,:),yMap(:,1),C(:,:,1))
shading flat
axis equal
axis xy
colormap(flipud(repmat(colorIndex',1,3)))
hold on
freezeColors
%contour3(xMap,yMap,zMesh,50,'LineColor','k');
%clabel(c,h)
end
