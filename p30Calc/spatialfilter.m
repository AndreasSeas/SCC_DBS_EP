function [sfe] = spatialfilter(V,x,y,dscale)
% - V: voltage in indexes for spatial filtering (electrodes)
    % for simplicity, just have one timepoint rn
% - x, y: x and y coordinates of these data
% - dscale: scale of the distance value
% https://www.frontiersin.org/articles/10.3389/fneur.2019.00325/full
% • For each electrode, the values of the 6 closest neighbors are 
%       determined, plus the central electrode value itself.
% • The 7 data points are sorted.
% • The minimal and maximal values are removed by dropping the first and 
%       last items of this list.
% • The remaining 5 values are then averaged, with weights proportional to 
%       the inverse distance to the central electrode. The central 
%       electrode is given a weight of 1.

% modification, will use the direct closest neighbors based on delauney
% triangulation instead

%%%%%%%

V(isnan(V))=[];

[~,I]=sort(V);% sort V
Ikeep=I(2:end-1);
Vkeep=V(Ikeep);% get the ones to average
dist=sqrt((x(1)-x).^2+(y(1)-y).^2)*dscale;
dist(1)=1;% set this weight to 1

sfe=sum(Vkeep./dist(Ikeep))/sum(1./dist(Ikeep));

% future iterations
% https://arxiv.org/abs/1702.02914
% https://journals-physiology-org.proxy.lib.duke.edu/doi/full/10.1152/jn.00560.2019