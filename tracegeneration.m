% Author    : Andreas Seas
% Descr     : Generate figure 5 for SCC DBS DEP Manuscript
%               - Adapted from M2_Code_raters_V2_seastrytomakethebox

%% clear slate
close all;clear;clc

%% Inputs
% colors
CBcolor = [0.8, 0.6, 0.1];%[0.4000    0.3647    0.1059];
FMcolor = "r";
TCcolor = "b";

% alphas
facealph = 0.60;
edgealph = 1;

% labels
CBlabel = "Cingulum Bundle Selective";
FMlabel = "Forceps Minor Selective";
TClabel = "Non-Selective";

% Example EPs
i_EP1=110;%38;
i_EP2=89;%4;
i_EP1_ex=9;
i_EP2_ex=29;

% tiled layout param
tightness='tight';
fpos=[0,0,1000,1500];

% figure name
% mkdir("figure5_v4");
fignamebase="figure5_v5";

% save figs?
flag_figsave=true;
flag_figdisp=true;


%% run setup file
[DB,V,xy,chanmap,homedir]=setmeup();

fulldata=readtable('C:\Users\as822\Box\collab-Sinai\manuscripts\M2 - Spatiotemporal correlates of SCC DBS\Figures\Source Files and Code\Background Calculations\FullData.xlsx');

%% 

return

% %% load ratings
% R1="C:\Users\as822\Box\collab-Sinai\manuscripts\M2 - Spatiotemporal correlates of SCC DBS\raterCagResults\IndividualRunDB_modifiedlabeledtree_Bryan.xlsx";
% R1=readtable(R1);
% 
% R2="C:\Users\as822\Box\collab-Sinai\manuscripts\M2 - Spatiotemporal correlates of SCC DBS\raterCagResults\IndividualRunDB_modifiedlabeledtree_Andreas.xlsx";
% R2=readtable(R2);
% 
% R3="C:\Users\as822\Box\collab-Sinai\manuscripts\M2 - Spatiotemporal correlates of SCC DBS\raterCagResults\IndividualRunDB_modifiedlabeledtree_Sohail.xlsx";
% R3=readtable(R3);

%% do LR equivalence stuff
% get LRequivalency set
[LRequivalency] = LReq(chanmap);

% get coherence metric for timespan
[LRcoh] = LRcoherence(105:135, V, DB,LRequivalency);

R1 = addvars(R1,LRcoh);
R2 = addvars(R2,LRcoh);
R3 = addvars(R3,LRcoh);

%% do boxplots
% saveas


%% perform Fleiss' Kappa, and identify Noise/Complex v Signal
% https://www.statology.org/fleiss-kappa-excel/#:~:text=The%20actual%20formula%20used%20to,1%20%E2%80%93%200.2128)%20%3D%200.2099.
% https://en.wikipedia.org/wiki/Fleiss%27_kappa
classID_R123=[R1.ManualClass,R2.ManualClass,R3.ManualClass];
classID=cell(numel(classID_R123)/3,1);
agr_mat=zeros(numel(classID_R123)/3,3);
for i=1:numel(classID_R123)/3
    agr_mat(i,:)=[sum(strcmp(classID_R123(i,:),'EP1')),...
        sum(strcmp(classID_R123(i,:),'EP2')),...
        sum(strcmp(classID_R123(i,:),'NR'))];

    if agr_mat(i,3)>1
        classID{i}='NR';
    elseif agr_mat(i,1)>1
        classID{i}='EP1';
    elseif agr_mat(i,2)>1
        classID{i}='EP2';
    else
        classID{i}='NoAgr';
    end

end

N = numel(classID_R123)/3;% number of samples
n = 3;% number of raters
k = 3;% number of options

pj=sum(agr_mat,1)./(N*n);
sumsq=sum(agr_mat.^2-agr_mat,2);%=(J2^2-J2)+(K2^2-K2)+(L2^2-L2);
Pi=(1/(n*(n-1))).*sumsq;

Pbar=sum(Pi)/N;
Pbar_e=sum(pj.^2);
kappa = (Pbar-Pbar_e)/(1-Pbar_e);

%% Figure 5A
fpos=[100,100,1400,500];
fid="A";
switch flag_figdisp
    case true
        f=figure('visible','on','Position',fpos);hold on;
    case false
        f=figure('visible','off','Position',fpos);hold on;
end
i_signal = ~strcmp(classID,"NoAgr") & ~strcmp(classID,"NR"); % get idx
i_CB = strcmp(DB.Type_real,"CB") & i_signal;
i_FM = strcmp(DB.Type_real,"FM") & i_signal;
i_TC = strcmp(DB.Type_real,"OT") & i_signal;

DB_=DB;% plot scatter
C=scatter(DB_.ActCB(i_CB),DB_.ActFM(i_CB),100,'Filled', ...
    'MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);
F=scatter(DB_.ActCB(i_FM),DB_.ActFM(i_FM),100,'Filled', ...
    'MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);
T=scatter(DB_.ActCB(i_TC),DB_.ActFM(i_TC),100,'Filled', ...
    'MarkerFaceColor',TCcolor ,'MarkerEdgeColor',TCcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);

XLIM=xlim;YLIM=ylim;% add line of unity
U=plot([0,100],[0,100],'--k','DisplayName','Line of Unity');
xlim(XLIM);ylim(YLIM);

xlabel('% Activation of CB'); % Adjust Formatting
ylabel('% Activation of FM')
C.DisplayName=CBlabel;
F.DisplayName=FMlabel;
T.DisplayName=TClabel;
legend([C,F,T,U],'location','northeast');
set(gca,'FontSize',16);

switch flag_figsave
    case true
        exportgraphics(f,fignamebase+fid+".eps");
        exportgraphics(f,fignamebase+fid+".png",'Resolution',600);
        saveas(f,fignamebase+fid+".fig");
        saveas(f,fignamebase+fid+".svg");
        saveas(f,fignamebase+fid+".pdf");
end
return
%% Figure 5B
fpos=[100,100,700,500];
fid="B";
switch flag_figdisp
    case true
        f=figure('visible','on','Position',fpos);hold on;
    case false
        f=figure('visible','off','Position',fpos);hold on;
end

i_EP = strcmp(classID,"EP1"); % get EP1 idx
set(gca,'FontSize',16);
title('EP1',FontSize=20)

xlim(XLIM);ylim(YLIM);


scatter(DB_.ActCB(i_CB & i_EP),DB_.ActFM(i_CB & i_EP),100,'Filled', ...
    'MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
    'Marker',"s");
scatter(DB_.ActCB(i_FM & i_EP),DB_.ActFM(i_FM & i_EP),100,'Filled', ...
    'MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
    'Marker',"s");
scatter(DB_.ActCB(i_TC & i_EP),DB_.ActFM(i_TC & i_EP),100,'Filled', ...
    'MarkerFaceColor',TCcolor ,'MarkerEdgeColor',TCcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
    'Marker',"s");

XLIM=xlim;YLIM=ylim;% add line of unity
U=plot([0,100],[0,100],'--k','DisplayName','Line of Unity');


EPobj=scatter(-100,-100,100,'Filled', ...
    'MarkerFaceColor',"k",'MarkerEdgeColor',"k", ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
    'Marker',"s","DisplayName","EP1");

xlim(XLIM);ylim(YLIM);
% 
% UEP=scatter(DB_.ActCB(i_EP1),DB_.ActFM(i_EP1),300, ...
%     'MarkerEdgeColor',"k", "linewidth",2,...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"o","DisplayName","Characteristic EP1");
% 
% scatter(DB_.ActCB(i_EP1_ex),DB_.ActFM(i_EP1_ex),300, ...
%     'MarkerEdgeColor',"k", "linewidth",2,...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"o","DisplayName","EP1 exception");

xlabel('% Activation of CB',FontSize=16); % Adjust Formatting
ylabel('% Activation of FM',FontSize=16)

text(50,23,num2str(sum(i_CB & i_EP)*100/sum(i_CB),3) + "% of CB",Color=CBcolor,FontSize=20);
text(50,20,num2str(sum(i_FM & i_EP)*100/sum(i_FM),3) + "% of FM",Color=FMcolor,FontSize=20);
text(50,17,num2str(sum(i_TC & i_EP)*100/sum(i_TC),3) + "% of TC",Color=TCcolor,FontSize=20);

switch flag_figsave
    case true
        exportgraphics(f,fignamebase+fid+".eps");
        exportgraphics(f,fignamebase+fid+".png",'Resolution',600);
        saveas(f,fignamebase+fid+".fig");
        saveas(f,fignamebase+fid+".svg");
end
exportgraphics(f,"C:\Users\as822\Downloads\EP1.png",'Resolution',600);

%% Figure 5C
fpos=[100,100,700,500];
fid="C";
switch flag_figdisp
    case true
        f=figure('visible','on','Position',fpos);hold on;
    case false
        f=figure('visible','off','Position',fpos);hold on;
end

i_EP = strcmp(classID,"EP2"); % get EP1 idx
set(gca,'FontSize',16);
title('EP2',FontSize=20)

xlim(XLIM);ylim(YLIM);


scatter(DB_.ActCB(i_CB & i_EP),DB_.ActFM(i_CB & i_EP),100,'Filled', ...
    'MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
    'Marker',"v");
scatter(DB_.ActCB(i_FM & i_EP),DB_.ActFM(i_FM & i_EP),100,'Filled', ...
    'MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
    'Marker',"v");
scatter(DB_.ActCB(i_TC & i_EP),DB_.ActFM(i_TC & i_EP),100,'Filled', ...
    'MarkerFaceColor',TCcolor ,'MarkerEdgeColor',TCcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
    'Marker',"v"); 

XLIM=xlim;YLIM=ylim;% add line of unity
U=plot([0,100],[0,100],'--k','DisplayName','Line of Unity');


EPobj=scatter(-100,-100,100,'Filled', ...
    'MarkerFaceColor',"k",'MarkerEdgeColor',"k", ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
    'Marker',"s","DisplayName","EP2");

xlim(XLIM);ylim(YLIM);

% UEP=scatter(DB_.ActCB(i_EP2),DB_.ActFM(i_EP2),300, ...
%     'MarkerEdgeColor',"k", "linewidth",2,...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"o","DisplayName","Characteristic EP2");
% 
% scatter(DB_.ActCB(i_EP2_ex),DB_.ActFM(i_EP2_ex),300, ...
%     'MarkerEdgeColor',"k", "linewidth",2,...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"o","DisplayName","EP2 exception");

xlabel('% Activation of CB',FontSize=16); % Adjust Formatting
ylabel('% Activation of FM',FontSize=16)

text(50,23,num2str(sum(i_CB & i_EP)*100/sum(i_CB),3) + "% of CB",Color=CBcolor,FontSize=20);
text(50,20,num2str(sum(i_FM & i_EP)*100/sum(i_FM),3) + "% of FM",Color=FMcolor,FontSize=20);
text(50,17,num2str(sum(i_TC & i_EP)*100/sum(i_TC),3) + "% of TC",Color=TCcolor,FontSize=20);

switch flag_figsave
    case true
        exportgraphics(f,fignamebase+fid+".eps");
        exportgraphics(f,fignamebase+fid+".png",'Resolution',600);
        saveas(f,fignamebase+fid+".fig");
        saveas(f,fignamebase+fid+".svg");
end

exportgraphics(f,"C:\Users\as822\Downloads\EP2.png",'Resolution',600);
% return
% figure;i_EP = strcmp(classID,"EP2"); histogram(DB_.ActCB(i_EP),BinEdges=0:5:100)

%% remake critical topograms
% get indices
topo_idx=[i_EP1,i_EP2,i_EP1_ex,i_EP2_ex];
savename={'EP1','EP2','EP1_exception','EP2_exception'};
% identify idxrem
idxrem={[];
    [165,94,111,145,241];
    [90,101,89,129,128,130,241,31,37,36,256,35,140];
    [130,129,131,100,89,142,110,90,101,144,252,253,254,73]};

% reference code for running to identify which channels to remove
% meanvolts=V(:,:,i_EP2_ex);
% temp=meanvolts(:,70:end)';
% figure;
% s=plot(temp);
% 
% for i=1:257
% row = dataTipTextRow('idx',repmat(i,432,1));
% s(i).DataTipTemplate.DataTipRows(end+1) = row;
% end

% get delauney triangulation of EEG
DT_eeg = delaunayTriangulation(xy(:,1),xy(:,2));

% set up plotting
tval=[30,60,90,120] + 50;
[X,Y] = meshgrid(-1:0.01:1, -1:0.01:1);
xq = X(:);
yq = Y(:);
rTest = sqrt(xq.^2 + yq.^2);
ckUnit = rTest <= 1;
xq = xq(ckUnit);
yq = yq(ckUnit);
x=xq;
y=yq;
dt = delaunayTriangulation(x,y) ;
nInt = length(xq);
tri = dt.ConnectivityList ;
xi = dt.Points(:,1) ;
yi = dt.Points(:,2) ;

for i = 1:numel(topo_idx)
    f=figure(); hold on; grid on; box off;
    tiledlayout(1,4,'TileSpacing','none');
    
    % get voltages for this pt
    meanvolts=V(:,:,topo_idx(i));
    
    % find nan rows and average in triangulation matrix:

    % for this set of idxrem, come up with averaging indices to average it out
    
    if numel(idxrem{i})>0
        meanvolts(idxrem{i},:)=NaN;
        avgidx=cell(numel(idxrem{i}),1);
        for j=1:numel(avgidx)
            idxlist=unique(DT_eeg.ConnectivityList(find(sum(DT_eeg.ConnectivityList==idxrem{i}(j),2)),:));
            avgidx{j}=idxlist;
            meanvolts(idxrem{i}(j),:)=mean(meanvolts(avgidx{j},:),'omitnan');
    
        end
    end

    [ntimes,nchan]=size(meanvolts);
    EEG_int = zeros(numel(xq), ntimes);
    maxval=max(max(abs(meanvolts(:,tval))));
    colorbds=[-maxval',maxval'];

    for col=1:4
        nexttile
%         maxval=max(max(abs(meanvolts(:,tval(col)))));
%         colorbds=[-maxval',maxval'];

        z = griddata(xy(:,1), xy(:,2),meanvolts(1:end-1,tval(col)), x, y);
%         subplot(4,4,[row,col]); %hold on; grid off; box off; axis off;
        F = scatteredInterpolant(x,y,z);
        zi = F(xi,yi) ;
        h=trisurf(tri,xi,yi,zi) ;
        view(2)
        shading interp
        cRB=redblue();
        colormap(cRB);
        caxis(colorbds);
        daspect([1,1,1])
        set(gca,'visible','off')
%         set(gca,'Color','None')
%         set(gcf,'Color','None')

    end
    exportgraphics(f,"figure_5_"+savename{i}+".png",'Resolution',600);
    saveas(f,"figure_5_"+savename{i}+".fig")
end

return


%% create base figure
% % switch flag_figdisp
% %     case true
% %         f=figure('visible','on','Position',fpos);%[-2494 373 830 983]); 
% %     case false
% %         f=figure('visible','off','Position',fpos);%[-2494 373 830 983]);
% % end
% % t = tiledlayout(8,4,'TileSpacing',tightness);
% % 
% % % PANEL 1: plot FM/CB/TC distinguishing
% % nexttile([2,4]); grid off; hold on; % set up tiled layout for this panel
% % % daspect([1,1,1])
% i_signal = ~strcmp(classID,"NoAgr") & ~strcmp(classID,"NR"); % get idx
% i_CB = strcmp(DB.Type_real,"CB") & i_signal;
% i_FM = strcmp(DB.Type_real,"FM") & i_signal;
% i_TC = strcmp(DB.Type_real,"OT") & i_signal;
% 
% DB_=DB;% plot scatter
% C=scatter(DB_.ActCB(i_CB),DB_.ActFM(i_CB),100,'Filled', ...
%     'MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor, ...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);
% F=scatter(DB_.ActCB(i_FM),DB_.ActFM(i_FM),100,'Filled', ...
%     'MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor, ...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);
% T=scatter(DB_.ActCB(i_TC),DB_.ActFM(i_TC),100,'Filled', ...
%     'MarkerFaceColor',TCcolor ,'MarkerEdgeColor',TCcolor, ...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);
% 
% XLIM=xlim;YLIM=ylim;% add line of unity
% U=plot([0,100],[0,100],'--k','DisplayName','Line of Unity');
% xlim(XLIM);ylim(YLIM);
% 
% xlabel('% Activation of CB'); % Adjust Formatting
% ylabel('% Activation of FM')
% C.DisplayName=CBlabel;
% F.DisplayName=FMlabel;
% T.DisplayName=TClabel;
% legend([C,F,T,U],'location','northeast');
% % set(gca,'FontSize',20);
% 
% % PANEL 2: EP1 example
% nexttile([1,2]); % set up tiled layout for this panel
% nuname2read=DB.NUNAME{i_EP1};
% cd(homedir);cd ..\Data\'EEG Processing Pipeline Output_HPF'\topograms\OT\
% I=imread(nuname2read+"topo_hpf_only.png");
% cd(homedir)
% imshow(I)
% 
% % PANEL 3: EP2 example
% nexttile([1,2]); % set up tiled layout for this panel
% nuname2read=DB.NUNAME{i_EP2};
% cd(homedir);cd ..\Data\'EEG Processing Pipeline Output_HPF'\topograms\OT\
% I=imread(nuname2read+"topo_hpf_only.png");
% cd(homedir)
% imshow(I)
% 
% % PANEL 4: EP1 MRI example
% nexttile([2,2]); % set up tiled layout for this panel
% fill([0,0,1,1],[0,1,1,0],'w');
% text(0.5,0.5, "EP1 MRI Example",'HorizontalAlignment','center')
% set(gca,'Visible','off');
% 
% % PANEL 5: EP2 MRI example
% nexttile([2,2]); % set up tiled layout for this panel
% fill([0,0,1,1],[0,1,1,0],'w');
% text(0.5,0.5, "EP2 MRI Example",'HorizontalAlignment','center')
% set(gca,'Visible','off');
% 
% % PANEL 6: EP1 Overlay
% nexttile([2,2]); hold on; % set up tiled layout for this panel
% grid off
% % daspect([1,1,1])
% i_EP = strcmp(classID,"EP1"); % get EP1 idx
% title('EP1')
% % create KDE background
% % gridx=-10:1:100;
% % gridy=-10:1:100;
% % [x,y]=meshgrid(gridx,gridy);
% % xi=[x(:),y(:)];
% % [z,ep] = ksdensity([DB_.ActCB(i_EP),DB_.ActFM(i_EP)],xi);
% % X = reshape(ep(:,1),length(gridx),length(gridy));
% % Y = reshape(ep(:,2),length(gridx),length(gridy));
% % Z = reshape(z,length(gridx),length(gridy));
% % [M,c]=contourf(X,Y,Z,5,"DisplayName","KDE Contour");
% % set(findobj(gca,'Type','patch','UserData',2),'EdgeColor',[1 1 1])
% % colormap('pink')
% % colormap(gray)
% % colormap([1,1,1; ...
% %     0.95,0.95,0.95; ...
% %     0.9,0.9,0.9; ...
% %     0.85,0.85,0.85;...
% %     0.8,0.8,0.8; ...
% %     0.75,0.75,0.75])
% xlim(XLIM);ylim(YLIM);
% 
% % xg=0:10:100;yg=0:10:100;
% % xgridLines = [xg, repelem(xg(1),1,numel(yg)); xg, repelem(xg(end),1,numel(xg))];
% % ygridLines = [repelem(yg(1),1,numel(xg)), yg; repelem(yg(end),1,numel(xg)), yg];
% % plot(xgridLines, ygridLines, 'Color', [.8 .8 .8]);
% 
% scatter(DB_.ActCB(i_CB & i_EP),DB_.ActFM(i_CB & i_EP),100,'Filled', ...
%     'MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor, ...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"s");
% scatter(DB_.ActCB(i_FM & i_EP),DB_.ActFM(i_FM & i_EP),100,'Filled', ...
%     'MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor, ...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"s");
% scatter(DB_.ActCB(i_TC & i_EP),DB_.ActFM(i_TC & i_EP),100,'Filled', ...
%     'MarkerFaceColor',TCcolor ,'MarkerEdgeColor',TCcolor, ...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"s");
% 
% % plot(mean(DB_.ActCB(i_EP)),mean(DB_.ActFM(i_EP)),'xk',MarkerSize=30);
% 
% % 
% % scatter(DB_.ActCB(i_CB & ~i_EP),DB_.ActFM(i_CB & ~i_EP),100,'Filled', ...
% %     'MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor, ...
% %     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
% %     'Marker',".");
% % scatter(DB_.ActCB(i_FM & ~i_EP),DB_.ActFM(i_FM & ~i_EP),100,'Filled', ...
% %     'MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor, ...
% %     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
% %     'Marker',".");
% % scatter(DB_.ActCB(i_TC & ~i_EP),DB_.ActFM(i_TC & ~i_EP),100,'Filled', ...
% %     'MarkerFaceColor',TCcolor ,'MarkerEdgeColor',TCcolor, ...
% %     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
% %     'Marker',".");
% 
% XLIM=xlim;YLIM=ylim;% add line of unity
% U=plot([0,100],[0,100],'--k','DisplayName','Line of Unity');
% 
% 
% EPobj=scatter(-100,-100,100,'Filled', ...
%     'MarkerFaceColor',"k",'MarkerEdgeColor',"k", ...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"s","DisplayName","EP1");
% % notEPobj=scatter(-100,-100,100,'Filled', ...
% %     'MarkerFaceColor',"k" ,'MarkerEdgeColor',"k", ...
% %     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
% %     'Marker',".","DisplayName","not EP1");
% xlim(XLIM);ylim(YLIM);
% 
% UEP=scatter(DB_.ActCB(i_EP1),DB_.ActFM(i_EP1),300, ...
%     'MarkerEdgeColor',"k", "linewidth",3,...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"o","DisplayName","Characteristic EP1");
% 
% scatter(DB_.ActCB(i_EP1_ex),DB_.ActFM(i_EP1_ex),300, ...
%     'MarkerEdgeColor',"k", "linewidth",3,...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"x","DisplayName","EP1 exception");
% 
% xlabel('% Activation of CB'); % Adjust Formatting
% ylabel('% Activation of FM')
% % legend([EPobj,notEPobj,U,UEP,c],'location','northeast');
% % legend([EPobj,U,UEP],'location','northeast');
% % set(gca,'FontSize',20);
% 
% text(50,25,num2str(sum(i_CB & i_EP)*100/sum(i_CB),3) + "% of CB == EP1",Color=CBcolor);
% text(50,20,num2str(sum(i_FM & i_EP)*100/sum(i_FM),3) + "% of FM == EP1",Color=FMcolor);
% text(50,15,num2str(sum(i_TC & i_EP)*100/sum(i_TC),3) + "% of TC == EP1",Color=TCcolor);
% 
% % PANEL 7: EP2 Overlay
% nexttile([2,2]); hold on; % set up tiled layout for this panel
% grid off
% % daspect([1,1,1])
% i_EP = strcmp(classID,"EP2"); % get EP2 idx
% 
% % create KDE background
% % gridx=-10:1:100;
% % gridy=-10:1:100;
% % [x,y]=meshgrid(gridx,gridy);
% % xi=[x(:),y(:)];
% % [z,ep] = ksdensity([DB_.ActCB(i_EP),DB_.ActFM(i_EP)],xi);
% % X = reshape(ep(:,1),length(gridx),length(gridy));
% % Y = reshape(ep(:,2),length(gridx),length(gridy));
% % Z = reshape(z,length(gridx),length(gridy));
% % [M,c]=contourf(X,Y,Z,5,"DisplayName","KDE Contour");
% % set(findobj(gca,'Type','patch','UserData',2),'EdgeColor',[1 1 1])
% % colormap('pink')
% % colormap(gray)
% % colormap([1,1,1; ...
% %     0.95,0.95,0.95; ...
% %     0.9,0.9,0.9; ...
% %     0.85,0.85,0.85;...
% %     0.8,0.8,0.8; ...
% %     0.75,0.75,0.75])
% xlim(XLIM);ylim(YLIM);
% 
% % xg=0:10:100;yg=0:10:100;
% % xgridLines = [xg, repelem(xg(1),1,numel(yg)); xg, repelem(xg(end),1,numel(xg))];
% % ygridLines = [repelem(yg(1),1,numel(xg)), yg; repelem(yg(end),1,numel(xg)), yg];
% % plot(xgridLines, ygridLines, 'Color', [.8 .8 .8]);
% 
% scatter(DB_.ActCB(i_CB & i_EP),DB_.ActFM(i_CB & i_EP),100,'Filled', ...
%     'MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor, ...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"v");
% scatter(DB_.ActCB(i_FM & i_EP),DB_.ActFM(i_FM & i_EP),100,'Filled', ...
%     'MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor, ...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"v");
% scatter(DB_.ActCB(i_TC & i_EP),DB_.ActFM(i_TC & i_EP),100,'Filled', ...
%     'MarkerFaceColor',TCcolor ,'MarkerEdgeColor',TCcolor, ...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"v");
% 
% % plot(mean(DB_.ActCB(i_EP)),mean(DB_.ActFM(i_EP)),'xk',MarkerSize=30);
% % 
% % scatter(DB_.ActCB(i_CB & ~i_EP),DB_.ActFM(i_CB & ~i_EP),100,'Filled', ...
% %     'MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor, ...
% %     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
% %     'Marker',".");
% % scatter(DB_.ActCB(i_FM & ~i_EP),DB_.ActFM(i_FM & ~i_EP),100,'Filled', ...
% %     'MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor, ...
% %     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
% %     'Marker',".");
% % scatter(DB_.ActCB(i_TC & ~i_EP),DB_.ActFM(i_TC & ~i_EP),100,'Filled', ...
% %     'MarkerFaceColor',TCcolor ,'MarkerEdgeColor',TCcolor, ...
% %     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
% %     'Marker',".");
% 
% XLIM=xlim;YLIM=ylim;% add line of unity
% U=plot([0,100],[0,100],'--k','DisplayName','Line of Unity');
% 
% 
% EPobj=scatter(-100,-100,100,'Filled', ...
%     'MarkerFaceColor',"k",'MarkerEdgeColor',"k", ...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"v","DisplayName","EP2");
% % notEPobj=scatter(-100,-100,100,'Filled', ...
% %     'MarkerFaceColor',"k" ,'MarkerEdgeColor',"k", ...
% %     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
% %     'Marker',".","DisplayName","not EP2");
% xlim(XLIM);ylim(YLIM);
% 
% UEP=scatter(DB_.ActCB(i_EP2),DB_.ActFM(i_EP2),300, ...
%     'MarkerEdgeColor',"k", "linewidth",3,...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"o","DisplayName","Characteristic EP2");
% 
% scatter(DB_.ActCB(i_EP2_ex),DB_.ActFM(i_EP2_ex),300, ...
%     'MarkerEdgeColor',"k", "linewidth",3,...
%     'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph,...
%     'Marker',"x","DisplayName","EP2 exception");
% 
% xlabel('% Activation of CB'); % Adjust Formatting
% ylabel('% Activation of FM')
% % legend([EPobj,notEPobj,U,UEP, c],'location','northeast');
% % legend([EPobj,U,UEP],'location','northeast');
% % set(gca,'FontSize',20);
% title('EP2')
% text(50,25,num2str(sum(i_CB & i_EP)*100/sum(i_CB),3) + "% of CB == EP2",Color=CBcolor);
% text(50,20,num2str(sum(i_FM & i_EP)*100/sum(i_FM),3) + "% of FM == EP2",Color=FMcolor);
% text(50,15,num2str(sum(i_TC & i_EP)*100/sum(i_TC),3) + "% of TC == EP2",Color=TCcolor);

% PANEL 8: EP1 example
nexttile([1,2]); % set up tiled layout for this panel
nuname2read=DB.NUNAME{i_EP1_ex};
cd(homedir);cd ..\Data\'EEG Processing Pipeline Output_HPF'\topograms\
I=imread(nuname2read+"topo_hpf_only.png");
cd(homedir)
imshow(I)

% PANEL 9: EP2 example
nexttile([1,2]); % set up tiled layout for this panel
nuname2read=DB.NUNAME{i_EP2_ex};
cd(homedir);cd ..\Data\'EEG Processing Pipeline Output_HPF'\topograms\
I=imread(nuname2read+"topo_hpf_only.png");
cd(homedir)
imshow(I)

switch flag_figsave
    case true
        exportgraphics(f,figname,'Resolution',600);
end
return
%% perform "voting"
% EPtype=


%% plot scatterplot
f=figure;hold on; grid on;
gscatter(R1.ActCB,R1.ActFM,R1.ManualClass)

f=figure;hold on; grid on;
gscatter(R2.ActCB,R2.ActFM,R2.ManualClass)

%% plot 3d scatterplot
f=figure;hold on; grid on;
scatter3(R1.ActCB(strcmp(R1.ManualClass,"EP1")),R1.ActFM(strcmp(R1.ManualClass,"EP1")),R1.LRcoh(strcmp(R1.ManualClass,"EP1")),50,'black','filled')
scatter3(R1.ActCB(strcmp(R1.ManualClass,"EP2")),R1.ActFM(strcmp(R1.ManualClass,"EP2")),R1.LRcoh(strcmp(R1.ManualClass,"EP2")),50,'green','filled')

f=figure;hold on; grid on;
scatter3(R2.ActCB(strcmp(R2.ManualClass,"EP1")),R2.ActFM(strcmp(R2.ManualClass,"EP1")),R2.LRcoh(strcmp(R2.ManualClass,"EP1")),50,'black','filled')
scatter3(R2.ActCB(strcmp(R2.ManualClass,"EP2")),R2.ActFM(strcmp(R2.ManualClass,"EP2")),R2.LRcoh(strcmp(R2.ManualClass,"EP2")),50,'green','filled')

%%
f=figure;hold on; grid on;
scatter3(R2.ActCB,R2.ActFM,R2.LRcoh,R2.ManualClass)

% [mdl]=scatterSVM(DB,i_CB,i_FM,CBcolor,FMcolor,facealph, edgealph);

return
%
%% figure 6 - SVM maker
% make model using only those indices that are non-noise based on visual
% inspection
i_signal = ~strcmp(DB.ManualClass,"Noise");
i_CB = strcmp(DB.Type_real,"CB") & i_signal;
i_FM = strcmp(DB.Type_real,"FM") & i_signal;
i_OT = strcmp(DB.Type_real,"OT") & i_signal;
i_symm = strcmp(DB.ManualClass,"Symmetric") & i_signal;
i_asymm = strcmp(DB.ManualClass,"Asymmetric") & i_signal;

% get LRequivalency set
[LRequivalency] = LReq(chanmap);

% get coherence metric for timespan
[LRcoh] = LRcoherence(15:45, V, DB,LRequivalency);

%% fig
% basic scatter
% simplescatter(DB,i_CB,i_FM,CBcolor,FMcolor,facealph, edgealph)

% SVM scatter in 2d
[mdl]=scatterSVM(DB,i_CB,i_FM,CBcolor,FMcolor,facealph, edgealph);

% naive bayes for fun
% [mdl]=scatterNB(DB_,i_CB,i_FM,CBcolor,FMcolor,facealph, edgealph);

% classification of controls
classifycontrols(DB,i_CB,i_FM,i_OT,CBcolor,FMcolor,TCcolor,0.7, edgealph,mdl)

% classscatter(DB,i_CB,i_FM,i_symm,i_asymm,CBcolor,FMcolor,facealph, edgealph)

allclass(DB,i_CB,i_FM,i_OT,i_signal,CBcolor,FMcolor,TCcolor,facealph, edgealph,LRcoh)

[mdl]=scatterSVM3(DB,i_CB,i_FM,CBcolor,FMcolor,facealph, edgealph,LRcoh);
% no need, makes a vertical plane

%% try tsne
%
% LRcohmat=[LRcoherence(15:45, V, DB,LRequivalency),...
%     LRcoherence(45:75, V, DB,LRequivalency),...
%     LRcoherence(75:105, V, DB,LRequivalency),...
%     LRcoherence(105:135, V, DB,LRequivalency)];
%
% Xall=[DB.ActCB,DB.ActFM,LRcohmat];
%
% X_=Xall(i_signal,:);
% label_=DB.ManualClass(i_signal);
% %DB.Type_real(i_signal);
%
% Y = tsne(X_,'');
% figure;
% gscatter(Y(:,1),Y(:,2),label_);

%% functions

% setmeup
function [DB,V,xy,chanmap,homedir]=setmeup()
%(independent of your naming of directories)
disp("running: setmeup");
homedir=pwd;
addpath(genpath(homedir));% add the path of all folders here
% cd('..');
% cd('Data\EEG Processing Pipeline Output_HPF\eegdata_mean\');
cd('C:\Users\as822\Box\collab-Sinai\project-selectEPs\EEG-dataSumm\eegdata_mean');
meandir=pwd;
cd(homedir);
cd('..');
cd('C:\Users\as822\Box\collab-Sinai\Code\EEG Processing Pipeline\referencefiles');
refdir=pwd;
addpath(genpath(refdir))
cd(homedir);

%% constants for all patients
eventloc=50;
triallength=501;

%% load metadata
% raw=readtable('IndividualRunDB_220523.xlsx');
% raw=readtable('IndividualRunDB_220617(adjusteddescriptions).xlsx');
cd(refdir)
raw=readtable('C:\Users\as822\Box\collab-Sinai\manuscripts\M2 - Spatiotemporal correlates of SCC DBS\IndividualRunDB_220617_modifiedlabeledtree.xlsx');
xy = load('adultAvg_xy256_unitCircle.txt');% load 2d projection
rSens = load('adultAvg_xyz256_unitSphere.txt');% load 3d projection
rInt = load('unitSphere_1000pts.txt');% load 1000 point interpolation
chanmap=readtable('net256_channelMap_edited.xlsx');
cd(homedir)
%% get only the specific trials you want
F=zeros(size(raw.ID));
F=F | strcmp(raw.Type_real,"CB");% get CB ones
F=F | strcmp(raw.Type_real,"FM");% get FM ones
F=F | strcmp(raw.Type_real,"OT");% get OT ones
% F=F & raw.AmpIDX==4;% get max amplitude only
DB=removevars(raw(F,:),{'NOTES','EPNotes','DEP_type'});
DB(any(ismissing(DB),2), :) = [];% only get the ones with all the data we need
% size(raw)

%% create the db for the ones you want to test
[L,~]=size(DB);
cd(meandir)
for i=1:L
    load(DB.NUNAME{i}+".mat");
    V(:,:,i)=meanvolts;
end
cd(homedir)

% rmpath(genpath(homedir));% remove path
% rmpath(genpath(refdir))
end

% LReq
function [LRequivalency] = LReq(chanmap)
% identify left-right equivalent electrodes
% make adjacency matrix

disp("running: LReq");

xyadj=nan(numel(chanmap.x256_Index),numel(chanmap.x256_Index));

for i=1:numel(chanmap.x256_Index)
    xyadj(:,i)=sqrt((chanmap.X-chanmap.X(i)).^2+(chanmap.Y-chanmap.Y(i)).^2);
end

% find matching XY pairs
LRequivalency=[];% unknown rows, 2 columns (left and right)
for i=1:numel(chanmap.x256_Index)
    if sum(LRequivalency==i)==0
        if chanmap.X(i)>0
            R=i;
            L=find(chanmap.X==-chanmap.X(i));

        elseif chanmap.X(i)<0
            L=i;
            R=find(chanmap.X==-chanmap.X(i));
        end
    end
    if numel(L)>0 && numel(R)>0
        LRequivalency=[LRequivalency;[L,R]];
    end
end
LRequivalency=unique(LRequivalency,'rows');
end

% LRcoherence
function [LRcoh] = LRcoherence(rng, V, DB,LRequivalency)

disp("running: LRcoherence");

[n,~]=size(DB);
LRcoh=zeros(n,1);
% rng=ctr(1)-radius:ctr(1)+radius;% only get the first range
for i=1:n

    tempcohi=zeros(numel(LRequivalency)/2,1);
    for j=1:numel(LRequivalency)/2
        tempmat=corrcoef(V(LRequivalency(j,1),rng,i),V(LRequivalency(j,2),rng,i));
        tempcohi(j)=tempmat(2,1);
    end
    LRcoh(i)=mean(tempcohi,'omitnan');
end


end

% simplescatter
function simplescatter(DB_,i_CB,i_FM,CBcolor,FMcolor,facealph, edgealph)

ff = figure; hold on; grid on; box off;
c=scatter(DB_.ActCB(i_CB),DB_.ActFM(i_CB),100,'Filled','MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
f=scatter(DB_.ActCB(i_FM),DB_.ActFM(i_FM),100,'Filled','MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
xlabel('% Activation of CB'); ylabel('% Activation of FM')
c.DisplayName="Cingulum Bundle Selective";
f.DisplayName="Forceps Minor Selective";
legend([c,f],'location','northeast');
daspect([1,1,1])
set(gca,'FontSize',20);

end

% % simplescatter
% function simplescatter(DB_,i_CB,i_FM,CBcolor,FMcolor,facealph, edgealph)
%
% ff = figure; hold on; grid on; box off;
% c=scatter(DB_.ActCB(i_CB),DB_.ActFM(i_CB),100,'Filled','MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
% f=scatter(DB_.ActCB(i_FM),DB_.ActFM(i_FM),100,'Filled','MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
% xlabel('% Activation of CB'); ylabel('% Activation of FM')
% c.DisplayName="Cingulum Bundle Selective";
% f.DisplayName="Forceps Minor Selective";
% legend([c,f],'location','northeast');
% daspect([1,1,1])
% set(gca,'FontSize',20);
%
% end


% classscatter
function classscatter(DB_,i_CB,i_FM,i_symm,i_asymm,CBcolor,FMcolor,facealph, edgealph)

i_symm = i_symm & (i_CB | i_FM);
i_asymm = i_asymm & (i_CB | i_FM);

ff = figure; hold on; grid on; box off;
S=scatter(DB_.ActCB(i_symm),DB_.ActFM(i_symm),300,'Filled','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
A=scatter(DB_.ActCB(i_asymm),DB_.ActFM(i_asymm),300,'Filled','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');

c=scatter(DB_.ActCB(i_CB),DB_.ActFM(i_CB),100,'Filled','MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
f=scatter(DB_.ActCB(i_FM),DB_.ActFM(i_FM),100,'Filled','MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');

xlabel('% Activation of CB'); ylabel('% Activation of FM')
S.DisplayName = "Symmetric";
A.DisplayName = "Asymmetric";
c.DisplayName="Cingulum Bundle Selective";
f.DisplayName="Forceps Minor Selective";
legend([S,A,c,f],'location','northeast');
daspect([1,1,1])
set(gca,'FontSize',20);
xlim([0,100]);ylim([0,100]);
end

% allclass
function allclass(DB_,i_CB,i_FM,i_OT,i_signal,CBcolor,FMcolor,TCcolor,facealph, edgealph,LRcoh)
i_slight=strcmp(DB_.Subgroup,"SlightAsymm") & i_signal;
i_symm=strcmp(DB_.ManualClass,"Symmetric") & i_signal & ~i_slight;
i_asymm=strcmp(DB_.ManualClass,"Asymmetric") & i_signal & ~i_slight;

ff = figure; hold on; grid on; box off;
S=scatter3(DB_.ActCB(i_symm),DB_.ActFM(i_symm),LRcoh(i_symm),300,'MarkerEdgeColor','g','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',edgealph,'LineWidth',2);%,'MarkerEdgeColor','none');
A=scatter3(DB_.ActCB(i_asymm),DB_.ActFM(i_asymm),LRcoh(i_asymm),300,'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',edgealph,'LineWidth',2);%,'MarkerEdgeColor','none');
SA=scatter3(DB_.ActCB(i_slight),DB_.ActFM(i_slight),LRcoh(i_slight),300,'MarkerEdgeColor','m','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',edgealph,'LineWidth',2);%,'MarkerEdgeColor','none');

% S=scatter3(DB_.ActCB(i_symm),DB_.ActFM(i_symm),LRcoh(i_symm),300,'Filled','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
% A=scatter3(DB_.ActCB(i_asymm),DB_.ActFM(i_asymm),LRcoh(i_asymm),300,'Filled','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
% SA=scatter3(DB_.ActCB(i_slight),DB_.ActFM(i_slight),LRcoh(i_slight),300,'Filled','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
c=scatter3(DB_.ActCB(i_CB),DB_.ActFM(i_CB),LRcoh(i_CB),100,'Filled','MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
f=scatter3(DB_.ActCB(i_FM),DB_.ActFM(i_FM),LRcoh(i_FM),100,'Filled','MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
t=scatter3(DB_.ActCB(i_OT),DB_.ActFM(i_OT),LRcoh(i_OT),100,'Filled','MarkerFaceColor',TCcolor ,'MarkerEdgeColor',TCcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');

xlabel('% Activation of CB'); ylabel('% Activation of FM')
S.DisplayName = "Symmetric";
A.DisplayName = "Asymmetric";
SA.DisplayName = "Slight Asymmetry";
c.DisplayName="Cingulum Bundle Selective";
f.DisplayName="Forceps Minor Selective";
t.DisplayName="Clinical Contact";
legend([S,A,SA,c,f,t],'location','northeast');
% daspect([1,1,1])
view(3)
set(gca,'FontSize',20);
xlim([0,100]);ylim([0,100]);
end

% scatter SVM
function [mdl]=scatterSVM(DB_,i_CB,i_FM,CBcolor,FMcolor,facealph, edgealph)

ff = figure; hold on; grid on; %box off;
c=scatter(DB_.ActCB(i_CB),DB_.ActFM(i_CB),100,'Filled','MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
f=scatter(DB_.ActCB(i_FM),DB_.ActFM(i_FM),100,'Filled','MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
xlabel('% Activation of CB'); ylabel('% Activation of FM')
c.DisplayName="Cingulum Bundle Selective";
f.DisplayName="Forceps Minor Selective";
legend([c,f],'location','northeast');
% daspect([1,1,1])
set(gca,'FontSize',20);

% get SVM
i_svm = i_CB | i_FM; % only include CB and FM in segmentation
xtest = [DB_.ActCB(i_svm),DB_.ActFM(i_svm)];
ytest = i_CB(i_svm);
mdl = fitcsvm(xtest,ytest);
% xlim([0,100]);ylim([0,100]);
yline=-([0,100]*mdl.Beta(1) + mdl.Bias)/mdl.Beta(2);
l = plot([0,100],yline,'-k','linewidth',2);
l.DisplayName="SVM boundary";
ybound1=-([0,100]*mdl.Beta(1) + mdl.Bias-1)/mdl.Beta(2);
ybound2=-([100,0]*mdl.Beta(1) + mdl.Bias+1)/mdl.Beta(2);
b = plot([0,100],ybound1,'--k');
plot([100,0],ybound2,'--k'); % top
b.DisplayName = "Margin";
% fill([0,100,100,0],[ybound1,ybound2],'k','FaceAlpha',0.1,'LineStyle','none');

legend([c,f,l,b],'location','northwest');
xlim([0,100]);ylim([0,100]);
daspect([1,1,1])


end

% scatter SVM
function [mdl]=scatterSVM3(DB_,i_CB,i_FM,CBcolor,FMcolor,facealph, edgealph,LRcoh)

ff = figure; hold on; grid on; %box off;
c=scatter3(DB_.ActCB(i_CB),DB_.ActFM(i_CB),LRcoh(i_CB),100,'Filled','MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
f=scatter3(DB_.ActCB(i_FM),DB_.ActFM(i_FM),LRcoh(i_FM),100,'Filled','MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
xlabel('% Activation of CB'); ylabel('% Activation of FM')
c.DisplayName="Cingulum Bundle Selective";
f.DisplayName="Forceps Minor Selective";
legend([c,f],'location','northeast');
% daspect([1,1,1])
set(gca,'FontSize',20);

% get SVM
i_svm = i_CB | i_FM; % only include CB and FM in segmentation
xtest = [DB_.ActCB(i_svm),DB_.ActFM(i_svm), LRcoh(i_svm)];
ytest = i_CB(i_svm);
mdl = fitcsvm(xtest,ytest);
X=xtest;Y=ytest;
sv =  mdl.SupportVectors;
d =1;% step size for fine sampling
[x, y, z] = meshgrid(min(X(:,1)):d:max(X(:,1)), ...
    min(X(:,2)):d:max(X(:,2)), ...
    min(X(:,3)):d:max(X(:,3))); % generate a simple grid for fine sampling
xGrid = [x(:),y(:),z(:)];
%get scores, f
[ temp , F] = predict(mdl,xGrid);
%reshape to same grid size as the input
F = reshape(F(:,2), size(x));
% Assume class labels are 1 and 0 and convert to logical
t = logical(Y);
%plot data points, color by class label
[faces,verts,~] = isosurface(x, y, z, F, 0, x);
xx = verts(:,1);
yy = verts(:,2);
zz = verts(:,3);
N = length(xx);
O = ones(N,1);
C = [xx yy O]\zz;

XL=xlim;YL=ylim;
d=5;
[X,Y] = meshgrid(XL(1):d:XL(2),YL(1):d:YL(2));
% [X,Y] = meshgrid(XL,YL);
ZZ= C(1)*X+C(2)*Y + C(3);

s = surf(X,Y,ZZ);%,'EdgeColor','k','facecolor','k','facealpha',0.2);
% % s.FaceColor = 'flat';
s.FaceColor = 'flat';
colormap('gray');
% caxis([0,1])

legend([c,f],'location','northwest');
xlim([0,100]);ylim([0,100]);
% daspect([1,1,1])


end


% scatter NB
function [mdl]=scatterNB(DB_,i_CB,i_FM,CBcolor,FMcolor,facealph, edgealph)

ff = figure; hold on; grid on; box off;
c=scatter3(DB_.ActCB(i_CB),DB_.ActFM(i_CB),ones(size(DB_.ActFM(i_CB))),100,'Filled','MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
f=scatter3(DB_.ActCB(i_FM),DB_.ActFM(i_FM),ones(size(DB_.ActFM(i_FM))),100,'Filled','MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor,'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);%,'MarkerEdgeColor','none');
xlabel('% Activation of CB'); ylabel('% Activation of FM')
c.DisplayName="Cingulum Bundle Selective";
f.DisplayName="Forceps Minor Selective";
legend([c,f],'location','northeast');
% daspect([1,1,1])
set(gca,'FontSize',20);

% get naive bayes
i_svm = i_CB | i_FM; % only include CB and FM in segmentation
xtest = [DB_.ActCB(i_svm),DB_.ActFM(i_svm)];
ytest = i_CB(i_svm);
mdl = fitcnb(xtest,ytest);
xLim=xlim;yLim=ylim;
[xx1, xx2] = meshgrid(0:0.1:100,0:0.1:100);
sz = size(xx1);
XGrid = [xx1(:) xx2(:)];
[predictedspecies,Posterior,~] = predict(mdl,XGrid);
% surf(xx1,xx2,reshape(Posterior(:,1),sz),'EdgeColor','none')
contourf(xx1,xx2,reshape(Posterior(:,2),sz))
colormap('cool')
colorbar
view(2)

xlim([0,100]);ylim([0,100]);
daspect([1,1,1])


end

% classify control (TC) cases
function classifycontrols(DB_,i_CB,i_FM,i_OT,CBcolor,FMcolor,TCcolor,facealph, edgealph,mdl)

xval=DB_.ActCB(i_OT);yval=DB_.ActFM(i_OT);
ybound=-(xval*mdl.Beta(1) +mdl.Bias)/mdl.Beta(2);
control_CB=(ybound-yval)>abs(1/mdl.Beta(2));
sum(control_CB)
control_FM=(yval-ybound)>abs(1/mdl.Beta(2));
sum(control_FM)
control_undecided=abs(yval-ybound)<abs(1/mdl.Beta(2));
sum(control_undecided)
ff = figure; hold on; grid on;% box off;
c=scatter(DB_.ActCB(i_CB),DB_.ActFM(i_CB),100,'Filled', ...
    'MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor, ...
    'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
f=scatter(DB_.ActCB(i_FM),DB_.ActFM(i_FM),100,'Filled', ...
    'MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor, ...
    'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
t_c=scatter(xval(control_CB),yval(control_CB),100,'d','Filled', ...
    'MarkerFaceColor',CBcolor ,'MarkerEdgeColor',CBcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);
t_f=scatter(xval(control_FM),yval(control_FM),100,'d','Filled', ...
    'MarkerFaceColor',FMcolor ,'MarkerEdgeColor',FMcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);
t_u=scatter(xval(control_undecided),yval(control_undecided),100,'d', ...
    'Filled','MarkerFaceColor',TCcolor ,'MarkerEdgeColor',TCcolor, ...
    'MarkerFaceAlpha',facealph,'MarkerEdgeAlpha',edgealph);
yline=-([0,100]*mdl.Beta(1) + mdl.Bias)/mdl.Beta(2);
l = plot([0,100],yline,'-k','linewidth',2);
l.DisplayName="SVM boundary";
ybound1=-([0,100]*mdl.Beta(1) + mdl.Bias-1)/mdl.Beta(2);
ybound2=-([100,0]*mdl.Beta(1) + mdl.Bias+1)/mdl.Beta(2);
b = plot([0,100],ybound1,'--k');
plot([100,0],ybound2,'--k'); % top
b.DisplayName = "Margin";
% fill([0,100,100,0],[ybound1,ybound2],'k','FaceAlpha',0.1,'LineStyle','none');
xlabel('% Activation of CB'); ylabel('% Activation of FM')
c.DisplayName="Cingulum Bundle Selective";
f.DisplayName="Forceps Minor Selective";
t_c.DisplayName = "Control - Cingulum Predominant";
t_f.DisplayName = "Control - Forceps Predominant";
t_u.DisplayName = "Control - Unclassified";
legend([c,f,t_c,t_f,t_u],'location','northwest');
daspect([1,1,1])
xlim([0,100]);ylim([0,100]);
set(gca,'FontSize',20);

end




%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%