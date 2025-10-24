
xContour = [0 10 10 50 50 0 0];
yContour = [0 0 10 10 40 40 0];
figure
plot(xContour, yContour, 'k-')

%%
xApex = 9;
yApex = 2;
threshold = 1;

[xVisi,yVisi,xChildApex,yChildApex] = visiPolygon(xContour,yContour,xApex,yApex,threshold,1);