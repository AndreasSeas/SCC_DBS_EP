function [x,y,xi,yi,tri,X,Y,dt]=setup_plot()
% set up the plot

[X,Y] = meshgrid(-1:0.01:1, -1:0.01:1);
xq = X(:);
yq = Y(:);
rTest = sqrt(xq.^2 + yq.^2);
ckUnit = rTest <= 1;
x = xq(ckUnit);
y = yq(ckUnit);
dt = delaunayTriangulation(x,y) ;
tri = dt.ConnectivityList ;
xi = dt.Points(:,1) ;
yi = dt.Points(:,2) ;
