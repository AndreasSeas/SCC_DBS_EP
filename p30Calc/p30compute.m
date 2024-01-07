% Author    : Andreas Seas
% Created   : 22-10-28
% Edited    : 22-10-28
% Descr     : generate figure 3 and its subfigures

%% clear slate
close all;clear;clc

%% Inputs
CBcolor = [255, 215, 0]/255;
FMcolor = [255, 0, 0]/255;
CCcolor = [0, 0, 213/255];
% save figs?
savefigs=false;
% set up fonts
s = settings;
s.matlab.fonts.editor.code.Name.TemporaryValue = 'Arial';
s.matlab.fonts.editor.code.Name.TemporaryValue = 'Arial';

%% set up the data
[DB,V,~,chanmap,homedir]=setmeup();
[num_chan,num_t,num_exp]=size(V);
clc

%% get full data
FullData=readtable('FullData.xlsx');
pause(0.1)
clc

%% Get connectivity list
xy=[chanmap.X,chanmap.Y];

DT_eeg = delaunayTriangulation(xy(:,1),xy(:,2));
avgidx=cell(num_chan,1);
for i=1:num_chan
    avgidx{i}=unique(DT_eeg.ConnectivityList(find(sum(DT_eeg.ConnectivityList==i,2)),:));
end

figure('Position',[834 42 830 983]);hold on;
triplot(DT_eeg,xy(:,1),xy(:,2));
daspect([1,1,1]);

scatter(xy(1,1),xy(1,2),100,"red",'filled','o')
scatter(xy(avgidx{1},1),xy(avgidx{1},2),50,"green",'filled','o')

%% do spatial filtration
V_spat=zeros(num_chan,num_exp);
t=30;
t_idx=t+50;
dscale=0.1;
w = waitbar(0,'spatial filtration');

for j=1:num_exp
    waitbar(j/num_exp,w,'working');
    for k=1:num_chan
        V_spat(k,j)=spatialfilter(V(avgidx{k},t_idx,j),xy(avgidx{k},1), ...
            xy(avgidx{k},2),dscale);
    end
end
close(w)

%% find location of peak + max voltage therein
[x,y,xi,yi,tri,X,Y,dt]=setup_plot();
w = waitbar(0,'spatial filtration');

Xpct=0.1;

x_peak=zeros(num_exp,1);
y_peak=zeros(num_exp,1);

x_trough=zeros(num_exp,1);
y_trough=zeros(num_exp,1);

mag_peak=zeros(num_exp,1);
mag_trough=zeros(num_exp,1);

idx_max_cell=cell(num_exp,1);
idx_min_cell=cell(num_exp,1);

outline_peak=cell(num_exp,1);
outline_trough=cell(num_exp,1);



for j=1:num_exp
    waitbar(j/num_exp,w,'working');
    
    z = griddata(xy(:,1), xy(:,2),V_spat(:,j), x, y);
    F = scatteredInterpolant(x,y,z);
    zi = F(xi,yi) ;

    [N,edges]=histcounts(z,'NumBins',100);
    numbins=0;
    chansadded=0;
    while chansadded<numel(z)*Xpct
        chansadded=chansadded+N(end-numbins);
        numbins=numbins+1;
    end

    lowerbd=edges(end-numbins);
    idx_max=find(zi>lowerbd);
    idx_max_cell{j}=idx_max;

    numbins=0;
    chansadded=0;
    while chansadded<numel(z)*Xpct
        chansadded=chansadded+N(numbins+1);
        numbins=numbins+1;
    end
    
    upperbd = edges(numbins+1);
    idx_min=find(zi<upperbd);
    idx_min_cell{j}=idx_min;
    
    x_peak(j)=mean(xi(idx_max));
    y_peak(j)=mean(yi(idx_max));
    mag_peak(j)=max(z);

    x_trough(j)=mean(xi(idx_min));
    y_trough(j)=mean(yi(idx_min));
    mag_trough(j)=min(z);

    x_max=xi(idx_max);
    y_max=yi(idx_max);
    dt_max=delaunayTriangulation(x_max,y_max);
    [F,P] = freeBoundary(dt_max);
    outline_peak{j}=[P;P(1,:)];

    x_min=xi(idx_min);
    y_min=yi(idx_min);
    dt_min=delaunayTriangulation(x_min,y_min);
    [F,P] = freeBoundary(dt_min);
    outline_trough{j}=[P;P(1,:)];
    
end
close(w)

% 
% %% export data
% ExportData=FullData;
% ExportData=addvars(ExportData,x_peak);
% ExportData=addvars(ExportData,y_peak);
% ExportData=addvars(ExportData,mag_peak);
% DEPtype=cell(num_exp,1);
% for i=1:num_exp
%     if FullData.ActFM(i)>0
%         DEPtype{i}='B-DEP';
%     else
%         DEPtype{i}='U-DEP';
%     end
% end
% 
% ExportData=addvars(ExportData,DEPtype);


% writetable(ExportData,'ExportData_p40.xlsx','Sheet',1)

% T = addvars(T,LastName,'Before',"Age");
% FullData=readtable('FullData.xlsx');

%% four column figure 6

% f=figure('Position',[100,100,500,500]); hold on; grid off;


idx1=find(FullData.ActFM>0 & FullData.ActCB>0);
idx2=find(FullData.ActFM>0 & FullData.ActCB==0);
idx3=find(FullData.ActFM==0 & FullData.ActCB>0 & ~strcmp(FullData.Type_Real,"EP1"));
idx4=find(FullData.ActFM==0 & FullData.ActCB>0 & strcmp(FullData.Type_Real,"EP1"));

mag_scatter={mag_peak(idx1),mag_peak(idx2),mag_peak(idx3),mag_peak(idx4)};
mag_scatter{1}=rmoutliers(mag_scatter{1},'quartiles');
mag_scatter{2}=rmoutliers(mag_scatter{2},'quartiles');
mag_scatter{3}=rmoutliers(mag_scatter{3},'quartiles');
mag_scatter{4}=rmoutliers(mag_scatter{4},'quartiles');

mag_min_scatter={mag_trough(idx1),mag_trough(idx2),mag_trough(idx3),mag_trough(idx4)};
mag_min_scatter{1}=rmoutliers(mag_min_scatter{1},'quartiles');
mag_min_scatter{2}=rmoutliers(mag_min_scatter{2},'quartiles');
mag_min_scatter{3}=rmoutliers(mag_min_scatter{3},'quartiles');
mag_min_scatter{4}=rmoutliers(mag_min_scatter{4},'quartiles');

% mag_scatter={FullData.P1_FP1_mag(idx1),FullData.P1_FP1_mag(idx2),FullData.P1_FP1_mag(idx3)};
% boxchart(mag_peak(idx),'MarkerStyle','none')
% scatter(ones(size(mag_peak(idx)))+(rand(size(mag_peak(idx)))*0.5)-0.25,mag_peak(idx),[],'k')

minmax=[-25,25];

f=figure('Position',[50,50 1550 800]); hold on; 

t=tiledlayout(2,4,'TileSpacing','compact');
nexttile;hold on;
% boxchart(mag_scatter{1},'MarkerStyle','none');


boxchart(mag_scatter{2},'MarkerStyle','none','BoxFaceColor','red');
scatter(ones(size(mag_scatter{2}))+(rand(size(mag_scatter{2}))*0.5)-0.25,mag_scatter{2},[],'r','filled')

boxchart(mag_min_scatter{2},'MarkerStyle','none','BoxFaceColor','blue');
scatter(ones(size(mag_min_scatter{2}))+(rand(size(mag_min_scatter{2}))*0.5)-0.25,mag_min_scatter{2},[],'b','filled')
xticklabels("")
ylabel('Magnitude of peak, \muV', 'Interpreter','Tex')
ylim(minmax)
title('FM > 0, CB == 0')
subtitle('(bilateral)')

nexttile;hold on;
boxchart(mag_scatter{1},'MarkerStyle','none','BoxFaceColor','red');
scatter(ones(size(mag_scatter{1}))+(rand(size(mag_scatter{1}))*0.5)-0.25,mag_scatter{1},[],'r','filled')

boxchart(mag_min_scatter{1},'MarkerStyle','none','BoxFaceColor','blue');
scatter(ones(size(mag_min_scatter{1}))+(rand(size(mag_min_scatter{1}))*0.5)-0.25,mag_min_scatter{1},[],'b','filled')
xticklabels("")


title('FM > 0, CB > 0')
subtitle('(bilateral)')
ylim(minmax)

nexttile;hold on;
boxchart(mag_scatter{4},'MarkerStyle','none','BoxFaceColor','red');
scatter(ones(size(mag_scatter{4}))+(rand(size(mag_scatter{4}))*0.5)-0.25,mag_scatter{4},[],'r','filled')

boxchart(mag_min_scatter{4},'MarkerStyle','none','BoxFaceColor','blue');
scatter(ones(size(mag_min_scatter{4}))+(rand(size(mag_min_scatter{4}))*0.5)-0.25,mag_min_scatter{4},[],'b','filled')
xticklabels("")
ylim(minmax)
title('FM == 0, CB > 0')
subtitle('(bilateral)')

nexttile;hold on;
boxchart(mag_scatter{3},'MarkerStyle','none','BoxFaceColor','red');
scatter(ones(size(mag_scatter{3}))+(rand(size(mag_scatter{3}))*0.5)-0.25,mag_scatter{3},[],'r','filled')

boxchart(mag_min_scatter{3},'MarkerStyle','none','BoxFaceColor','blue');
scatter(ones(size(mag_min_scatter{3}))+(rand(size(mag_min_scatter{3}))*0.5)-0.25,mag_min_scatter{3},[],'b','filled')
xticklabels("")
ylim(minmax)
title('FM == 0, CB > 0')
subtitle('(unilateral)')

%%%%%%%%%%%%%%

cbdmax=20;

nexttile;hold on;
daspect([1,1,1]);
plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
scatter(x_peak(idx2),y_peak(idx2),100,mag_peak(idx2),'filled','o','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.5,0.5,0.5])
scatter(x_trough(idx2),y_trough(idx2),100,mag_trough(idx2),'filled','o','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.5,0.5,0.5])
crb=redblue();
colormap(crb)
caxis([-cbdmax,cbdmax]);
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');

nexttile;hold on;
daspect([1,1,1]);
plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
scatter(x_peak(idx1),y_peak(idx1),100,mag_peak(idx1),'filled','o','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.5,0.5,0.5])
scatter(x_trough(idx1),y_trough(idx1),100,mag_trough(idx1),'filled','o','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.5,0.5,0.5])
crb=redblue();
colormap(crb)
caxis([-cbdmax,cbdmax]);
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');



nexttile;hold on;
daspect([1,1,1]);
plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
scatter(x_peak(idx4),y_peak(idx4),100,mag_peak(idx4),'filled','o','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.5,0.5,0.5])
scatter(x_trough(idx4),y_trough(idx4),100,mag_trough(idx4),'filled','o','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.5,0.5,0.5])
crb=redblue();
colormap(crb)
caxis([-cbdmax,cbdmax]);
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');
% a=colorbar();
% a.Label.String='Magnitude of peak, \muV';


nexttile;hold on;
daspect([1,1,1]);
plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
scatter(x_peak(idx3),y_peak(idx3),100,mag_peak(idx3),'filled','o','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.5,0.5,0.5])
scatter(x_trough(idx3),y_trough(idx3),100,mag_trough(idx3),'filled','o','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.5,0.5,0.5])
crb=redblue();
colormap(crb)
caxis([-cbdmax,cbdmax]);
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');
a=colorbar();
a.Label.String='Magnitude of peak, \muV';



% plot top 10 and bottom 10
% hold on;daspect([1,1,1]);
% plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
% for i=1:numel(idx3)
%     %     fill(outline_trough{idx3(i)}(:,1),outline_trough{idx3(i)}(:,2),'r', ...
%     %         'FaceAlpha',0.1,'EdgeColor','none')
%     %     fill(outline_trough{idx3(i)}(:,1),outline_trough{idx3(i)}(:,2),'r', ...
%     %         'FaceAlpha',0.1,'EdgeColor','none')
% 
%     scatter(x(idx_max_cell{idx3(i)}), ...
%         y(idx_max_cell{idx3(i)}),10,'blue','filled','MarkerFaceAlpha',0.01)
% 
%     scatter(x(idx_min_cell{idx3(i)}), ...
%         y(idx_min_cell{idx3(i)}),10,'red','filled','MarkerFaceAlpha',0.01)
%     return
% end
% set(gca,'XColor', 'none','YColor','none')
% set(gca, 'color', 'none');

set(findall(gcf,'-property','FontSize'),'FontSize',17)
saveas(f,"fig6_0517_v2.svg");

return
% get indices



B_DEP_idx=FullData.ActFM>0;
B_DEP_idx=FullData.ActFM>0;
U_DEP_idx=FullData.ActFM==0;
Response_idx=cellfun('length', strfind(FullData.Type_Real,'NR'))==0;
Left_idx=cellfun('length', strfind(FullData.Side,'L'));
Right_idx=cellfun('length', strfind(FullData.Side,'R'));

%% 

case1_idx=FullData.ActFM>0 & FullData.ActCB>0;
case2_idx=FullData.ActFM==0 & FullData.ActCB>0;
case3_idx=FullData.ActFM>0 & FullData.ActCB==0;
%     nexttile
%     hold on;
f=figure('Position',[100,100,500,500]); hold on
boxchart(mag_peak(case3_idx),'MarkerStyle','none')
% boxchart(mag_peak(case2_idx),'MarkerStyle','none')
% boxchart(mag_peak(case3_idx),'MarkerStyle','none')
scatter(ones(size(mag_peak(case3_idx)))+(rand(size(mag_peak(case3_idx)))*0.5)-0.25,mag_peak(case3_idx),[],'k')


%% plot big plot
% f=figure('Position',[1 41 1664 992]);hold on;
sidename={'Left','Left','Right','Right'};
depname={'B-DEP','U-DEP','B-DEP','U-DEP'};
sidefilt={Left_idx,Left_idx,Right_idx,Right_idx};
depfilt={B_DEP_idx,U_DEP_idx,B_DEP_idx,U_DEP_idx};

% tld=tiledlayout(3,4,'Padding', 'none', 'TileSpacing', 'compact'); 
cbdmax=max(mag_peak(Response_idx));

rng(1);
groups=cell(size(mag_peak));
sets=cell(4,1);
for j=1:4
    filt=find((depfilt{j}==1).*(sidefilt{j}==1).*(Response_idx==1));
    for i=1:numel(filt)
        groups{filt(i)}=[sidename{j},',',depname{j}];
    end

%     nexttile
%     hold on;
    f=figure('Position',[100,100,500,500]); hold on
    boxchart(mag_peak(filt),'MarkerStyle','none')
    scatter(ones(size(mag_peak(filt)))+(rand(size(mag_peak(filt)))*0.5)-0.25,mag_peak(filt),[],'k')
    sets{j}=mag_peak(filt);
    ylim([0,30]);
    title(sidename{j} + ", " + depname{j})

    if j==1
        ylabel('P30 Magnitude')
    end

    pd = fitdist(mag_peak(filt),'Normal');
    disp(pd);
    ci = paramci(pd);
    disp(ci)
    fignamebase="maxbox_"+sidename{j} + "_" + depname{j};
    exportgraphics(f,fignamebase+".eps");
    exportgraphics(f,fignamebase+".png",'Resolution',600);
    saveas(f,fignamebase+".fig");
    saveas(f,fignamebase+".svg");


end
[p,tbl,stats] = anova1(mag_peak(Response_idx),groups(Response_idx));

[h,p,ci,stats] = ttest2(sets{1},sets{2});
p
[h,p,ci,stats] = ttest2(sets{3},sets{4});
p

return
%%
for j=1:4
    filt=find((depfilt{j}==1).*(sidefilt{j}==1).*(Response_idx==1));
    % scatter
    % subplot(4,3,3*j-2);
    nexttile
    hold on;daspect([1,1,1]);
    plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
    scatter(x_peak(filt),y_peak(filt),100,mag_peak(filt),'filled','o','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.5,0.5,0.5])
    crb=redblue();
    colormap(crb)
    caxis([-cbdmax,cbdmax]);
    set(gca,'XColor', 'none','YColor','none')
    set(gca, 'color', 'none');
    if j==4
        colorbar()
    end
    
    if j==1
        text(-1.2,0,'Centroid of p30','HorizontalAlignment','center','rotation',90)
    end
end

for j=1:4
    filt=find((depfilt{j}==1).*(sidefilt{j}==1).*(Response_idx==1));
   % points
    % subplot(4,3,3*j-1);
    nexttile
    hold on;daspect([1,1,1]);
    plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
    for i=1:numel(filt)
        scatter(x(idx_max_cell{filt(i)}), ...
            y(idx_max_cell{filt(i)}),10,'blue','filled','MarkerFaceAlpha',0.01)
    end
    set(gca,'XColor', 'none','YColor','none')
    set(gca, 'color', 'none');
    
    if j==1
        text(-1.2,0,'top 10% using spatial filter','HorizontalAlignment','center','rotation',90)
    end

end
% saveas(f,'compositefigure_p40.png');

return

%% plot big plot
figure('Position',[834 42 830 983]);hold on;
sidename={'Left','Left','Right','Right'};
depname={'B-DEP','U-DEP','B-DEP','U-DEP'};
sidefilt={Left_idx,Left_idx,Right_idx,Right_idx};
depfilt={B_DEP_idx,U_DEP_idx,B_DEP_idx,U_DEP_idx};

tld=tiledlayout(4,3, 'Padding', 'none', 'TileSpacing', 'compact'); 
cbdmax=max(mag_peak(Response_idx));

for j=1:4
    
    filt=find((depfilt{j}==1).*(sidefilt{j}==1).*(Response_idx==1));
    % scatter
    % subplot(4,3,3*j-2);
    nexttile
    hold on;daspect([1,1,1]);
    plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
    scatter(x_peak(filt),y_peak(filt),100,mag_peak(filt),'filled','o','MarkerFaceAlpha',0.3)
    crb=redblue();
    colormap(crb)
    caxis([-cbdmax,cbdmax]);

    %     scatter(x_peak(filt),y_peak(filt),100,'blue','filled','o','MarkerFaceAlpha',0.3)
    
    set(gca,'XColor', 'none','YColor','none')
    set(gca, 'color', 'none');
    

    % points
    % subplot(4,3,3*j-1);
    nexttile
    hold on;daspect([1,1,1]);
    plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
    for i=1:numel(filt)
        scatter(x(idx_max_cell{filt(i)}), ...
            y(idx_max_cell{filt(i)}),10,'blue','filled','MarkerFaceAlpha',0.01)
    end
    set(gca,'XColor', 'none','YColor','none')
    set(gca, 'color', 'none');
    title(sidename{j} + ", " + depname{j})
    
    % contours
    % subplot(4,3,3*j); 
    nexttile
    hold on;daspect([1,1,1]);
    plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
    for i=1:numel(filt)
        fill(outline_peak{filt(i)}(:,1),outline_peak{filt(i)}(:,2),'b', ...
            'FaceAlpha',0.1,'EdgeColor','none')
    end
    set(gca,'XColor', 'none','YColor','none')
    set(gca, 'color', 'none');
end

return
%% plot big plot
figure('Position',[834 42 830 983]);hold on;
sidename={'Left','Left','Right','Right'};
depname={'B-DEP','U-DEP','B-DEP','U-DEP'};
sidefilt={Left_idx,Left_idx,Right_idx,Right_idx};
depfilt={B_DEP_idx,U_DEP_idx,B_DEP_idx,U_DEP_idx};

tiledlayout(4,3, 'Padding', 'none', 'TileSpacing', 'compact'); 

for j=1:4

    filt=find((depfilt{j}==1).*(sidefilt{j}==1).*(Response_idx==1));
    % scatter
    % subplot(4,3,3*j-2);
    nexttile
    hold on;daspect([1,1,1]);
    plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
    scatter(x_peak(filt),y_peak(filt),100,'blue','filled','o','MarkerFaceAlpha',0.3)
    
    set(gca,'XColor', 'none','YColor','none')
    set(gca, 'color', 'none');
    

    % points
    % subplot(4,3,3*j-1);
    nexttile
    hold on;daspect([1,1,1]);
    plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
    for i=1:numel(filt)
        scatter(x(idx_max_cell{filt(i)}), ...
            y(idx_max_cell{filt(i)}),10,'blue','filled','MarkerFaceAlpha',0.01)
    end
    set(gca,'XColor', 'none','YColor','none')
    set(gca, 'color', 'none');
    title(sidename{j} + ", " + depname{j})
    
    % contours
    % subplot(4,3,3*j); 
    nexttile
    hold on;daspect([1,1,1]);
    plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
    for i=1:numel(filt)
        fill(outline_peak{filt(i)}(:,1),outline_peak{filt(i)}(:,2),'b', ...
            'FaceAlpha',0.1,'EdgeColor','none')
    end
    set(gca,'XColor', 'none','YColor','none')
    set(gca, 'color', 'none');
end

%% plot the points
figure('Position',[834 42 830 983]);hold on;
plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
filt=find((U_DEP_idx==1).*(Left_idx==1).*(Response_idx==1));
for j=1:numel(filt)
    
    scatter(x(idx_max_cell{filt(j)}), ...
        y(idx_max_cell{filt(j)}),10,'blue','filled','MarkerFaceAlpha',0.1)

end

daspect([1,1,1]);
%% plot the contours
figure('Position',[834 42 830 983]);hold on;
plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
filt=find((U_DEP_idx==1).*(Left_idx==1).*(Response_idx==1));
for j=1:numel(filt)
    
    fill(outline_peak{filt(j)}(:,1),outline_peak{filt(j)}(:,2),'b','FaceAlpha',0.1)

end

daspect([1,1,1]);
%% plot the points
figure('Position',[834 42 830 983]);hold on;
plot(xy(:,1),xy(:,2),'.','Color',[0.7,0.7,0.7])
filt=find((B_DEP_idx==1).*(Right_idx==1).*(Response_idx==1));
% scatter(x_peak(filt),y_peak(filt),mag_peak(filt)*10,'blue','filled','o','MarkerFaceAlpha',0.3)
scatter(x_peak(filt),y_peak(filt),100,'blue','filled','o','MarkerFaceAlpha',0.3)
filt=find((U_DEP_idx==1).*(Right_idx==1).*(Response_idx==1));
% scatter(x_peak(filt),y_peak(filt),mag_peak(filt)*10,'red','filled','o','MarkerFaceAlpha',0.3)
scatter(x_peak(filt),y_peak(filt),100,'red','filled','o','MarkerFaceAlpha',0.3)
daspect([1,1,1]);
% figure;
% plot(FullData.ActCB,mag_peak,'ok')
% boxchart(FullData.Side,mag_peak)

% scatter(x_peak,y_peak,[],B_DEP_idx,'filled');

return
%% plot sample
[x,y,xi,yi,tri,X,Y,dt]=setup_plot();
% VV=V(:,t_idx,4);
VV=V_spat(:,50);
z = griddata(xy(:,1), xy(:,2),VV, x, y);
F = scatteredInterpolant(x,y,z);
zi = F(xi,yi) ;
figure;
h=trisurf(tri,xi,yi,zi) ;
view(2)
shading interp
cRB=redblue();
colormap(cRB);
colorbar;
[zmin,zmax]=caxis;
cbdmax=max(abs([zmin,zmax]));
caxis([-cbdmax,cbdmax]);
daspect([1,1,1])
% set(gca,'visible','off')

Z=zeros(size(X));
xidx=int64(x*100+101);
yidx=int64(y*100+101);
for j=1:numel(x)
    Z(xidx(j),yidx(j))=z(j);
end
hold on
% plot3(xy(idx_max,1),xy(idx_max,2),100*ones(size(xy(idx_max,2))),'ok')

%% get X% of max
Xpct=0.1;
[N,edges]=histcounts(z,'NumBins',100);
numbins=0;
chansadded=0;
while chansadded<numel(z)*Xpct
    chansadded=chansadded+N(end-numbins);
    numbins=numbins+1;

end

lowerbd=edges(end-numbins);

idx_max=find(zi>lowerbd);
x_max=xi(idx_max);
y_max=yi(idx_max);
dt_max=delaunayTriangulation(x_max,y_max);
[F,P] = freeBoundary(dt_max);
plot3(P(:,1),P(:,2),ones(numel(P(:,1)))*100,'-r','LineWidth',2)

%% plot dt
plot3(xi(idx_max),yi(idx_max),10*ones(size(xi(idx_max))),'.k')

return

%% plot sample
VV=V_spat(:,4);
z = griddata(xy(:,1), xy(:,2),VV, x, y);
F = scatteredInterpolant(x,y,z);
zi = F(xi,yi) ;
figure;
h=trisurf(tri,xi,yi,zi) ;
view(2)
shading interp
cRB=redblue();
colormap(cRB);
colorbar;
[zmin,zmax]=caxis;
cbdmax=max(abs([zmin,zmax]));
caxis([-cbdmax,cbdmax]);
daspect([1,1,1])
% set(gca,'visible','off')

% Z=zeros(size(X));
% xidx=int64(x*100+101);
% yidx=int64(y*100+101);
% for j=1:numel(x)
%     Z(xidx(j),yidx(j))=z(j);
% end
% hold on
hold on

idx_max=find(zi>lowerbd);
plot3(xi(idx_max),yi(idx_max),100*ones(size(xi(idx_max))),'ok')



return

% hold off
% figure;
% contourf(-X,-Y,Z,10,'LineWidth',1)

%%
growth_frac=0.2;
for j=4%1:num_exp
    [M,I]=max(zi);
    
    idx_max=I;% set initial index

    allconn=avgidx{I};
%     allconn=allconn(2:end);% outer boundary
    while numel(idx_max)<numel(zi)/100
        
        allconn=setdiff(allconn,idx_max);% remove internal and boundary ones
        vset=zi(allconn);
        allconn=allconn(isnan(vset)==0);% remove NaN
        vset=vset(isnan(vset)==0);

        [~,i_max]=sort(vset,'descend');
        
        ikeep=allconn(1:floor(numel(i_max)*growth_frac));
        
        idx_max=[idx_max;ikeep];
        
%         disp(idx_max);

        allconn=setdiff(unique(cat(1,avgidx{idx_max})),idx_max);
%         pause;
%         
% 
%         [~,min_i]=min(vset);
%         nbhd=avgidx{idx_max(min_i)};
%         nbhd_diff=setdiff(nbhd,idx_max);
%         % remove existing stuff
%         [~,max_i]=max(V_spat(nbhd_diff));
%         
%         idx_max=[idx_max,nbhd_diff(max_i)];

    end


end


%% plot sample
[x,y,xi,yi,tri,X,Y]=setup_plot();
% VV=V(:,t_idx,4);
VV=V_spat(:,4);
z = griddata(xy(:,1), xy(:,2),VV, x, y);
F = scatteredInterpolant(x,y,z);
zi = F(xi,yi) ;
figure;
h=trisurf(tri,xi,yi,zi) ;
view(2)
shading interp
cRB=redblue();
colormap(cRB);
colorbar;
[zmin,zmax]=caxis;
cbdmax=max(abs([zmin,zmax]));
caxis([-cbdmax,cbdmax]);
daspect([1,1,1])
% set(gca,'visible','off')

Z=zeros(size(X));
xidx=int64(x*100+101);
yidx=int64(y*100+101);
for j=1:numel(x)
    Z(xidx(j),yidx(j))=z(j);
end
hold on
plot3(xi(idx_max),yi(idx_max),100*ones(size(xi(idx_max))),'ok')



return
%% find maximum voltage contour
growth_frac=0.3;

avgidx=cell(numel(z),1);
for i=1:numel(z)
    avgidx{i}=unique(dt.ConnectivityList(find(sum(dt.ConnectivityList==i,2)),:));
end

%%
for j=4%1:num_exp
    [M,I]=max(zi);
    
    idx_max=I;% set initial index

    allconn=avgidx{I};
%     allconn=allconn(2:end);% outer boundary
    while numel(idx_max)<numel(zi)/100
        
        allconn=setdiff(allconn,idx_max);% remove internal and boundary ones
        vset=zi(allconn);
        allconn=allconn(isnan(vset)==0);% remove NaN
        vset=vset(isnan(vset)==0);

        [~,i_max]=sort(vset,'descend');
        
        ikeep=allconn(1:floor(numel(i_max)*growth_frac));
        
        idx_max=[idx_max;ikeep];
        
%         disp(idx_max);

        allconn=setdiff(unique(cat(1,avgidx{idx_max})),idx_max);
%         pause;
%         
% 
%         [~,min_i]=min(vset);
%         nbhd=avgidx{idx_max(min_i)};
%         nbhd_diff=setdiff(nbhd,idx_max);
%         % remove existing stuff
%         [~,max_i]=max(V_spat(nbhd_diff));
%         
%         idx_max=[idx_max,nbhd_diff(max_i)];

    end


end


%% plot sample
[x,y,xi,yi,tri,X,Y]=setup_plot();
% VV=V(:,t_idx,4);
VV=V_spat(:,4);
z = griddata(xy(:,1), xy(:,2),VV, x, y);
F = scatteredInterpolant(x,y,z);
zi = F(xi,yi) ;
figure;
h=trisurf(tri,xi,yi,zi) ;
view(2)
shading interp
cRB=redblue();
colormap(cRB);
colorbar;
[zmin,zmax]=caxis;
cbdmax=max(abs([zmin,zmax]));
caxis([-cbdmax,cbdmax]);
daspect([1,1,1])
% set(gca,'visible','off')

Z=zeros(size(X));
xidx=int64(x*100+101);
yidx=int64(y*100+101);
for j=1:numel(x)
    Z(xidx(j),yidx(j))=z(j);
end
hold on
plot3(xi(idx_max),yi(idx_max),100*ones(size(xi(idx_max))),'ok')
