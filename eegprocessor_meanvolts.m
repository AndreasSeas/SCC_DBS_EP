% Filename  : eegprocessor.m
% Author    : Andreas Seas
% Created   : 17-11-2021
% Edited    : 21-11-2021
% Descr     : ***
%             *** means fix me

%% clear slate
close all;clear;clc


%% Set up Workspace
disp("Set Up Workspace");
homedir=pwd;
addpath(genpath('C:\Users\as822\Box\PST Project'));
addpath(genpath(homedir))
topodir='C:\Users\as822\Box\SCC_DBS_Sinai\Data\EEG Processing Pipeline Output_HPF\topograms';
meandir='c:\Users\as822\Box\SCC_DBS_Sinai\Data\EEG Processing Pipeline Output_HPF\eegdata_mean\';
% %% Set up Workspace
% disp("Set Up Workspace");
% homedir=pwd;
% % addpath(genpath('/Users/as822/Box/PST Project'));% add the path of all eeg 
% % addpath(genpath(homedir));% add the path of all folders here
% % raweegdir='/Users/as822/Box/SCC_DBS_Sinai/Data/***';
% % tempdir='/Users/as822/Box/SCC_DBS_Sinai/Data/EEG Processing Pipeline Output_HPF/topograms';
% % eegdatadir='/Users/as822/Box/SCC_DBS_Sinai/Data/EEG Processing Pipeline Output_HPF/eegdatadir';

%% TEMP - definition of run info for testing (before function)
bool_plot=1;
bool_save=0;

%% constants for all patients
fSamp=1000;% Hz
fStim=2; % Hz
fNyquist = fSamp / 2; % Hz
crit_isi=400;% for finding artifacts
trialtrim=5;% number of trials to trim from beginning and end
eventloc=50;
triallength=501;


%% load metadata
DB=readtable('IndividualRunDB.xlsx');
xy = load('adultAvg_xy256_unitCircle.txt');% load 2d projection
rSens = load('adultAvg_xyz256_unitSphere.txt');% load 3d projection
rInt = load('unitSphere_1000pts.txt');% load 1000 point interpolation
chanmap=readtable('net256_channelMap_edited.xlsx');
% %% firpm optimal FIR filter
% order=91; % minimum frequency to resolve is 1 Hz... therefore need to 
% % be at least 1000 long
% f=[0, 1, 5, 50, 55,65, 70, 150, 160, fNyquist]./fNyquist;
% a=[0, 0, 1, 1, 0,0, 1, 1, 0, 0];
% b = firpm(order,f,a);
% [h,w] = freqz(b,1,512);
% plot(f,a,w/pi,abs(h))

% %% fir1 bandpass
% order=3000;
% Wn=[5, 50]./fNyquist;
% b=fir1(order,Wn);
% figure;
% freqz(b,1,512)

% %% design filter using firls
% order=3000; % minimum frequency to resolve is 1 Hz... therefore need to 
% % be at least 1000 long
% freqval=[0, 1, 5, 50, 55,65, 70, 150, 160, fNyquist]./fNyquist;
% ampval=[0, 0, 1, 1, 0,0, 1, 1, 0, 0];
% 
% b=firls(order,freqval,ampval);
% % fvtool(b);
% % this filter design removes 60 Hz noise, and allows for passage for most
% % if not all neural signal... use firls in order to minimize ringing

%% pick what type and side you wanna look at
mytype='CB';
myside='L';
db_idx=(strcmp(DB.Type,mytype)) & (strcmp(DB.Side,myside));
% DB2=DB(find(db_idx),:);
DB2=DB;
%% identify individual files
files2load=unique(DB2.Filename);

%% prepare outputs
sensor_data=cell(0,0);
patientid=cell(0,0);
% FP1=find((strcmp(chanmap.Sensor_Name,'FP1')));
% FP2=find((strcmp(chanmap.Sensor_Name,'FP2')));
% P3=find((strcmp(chanmap.Sensor_Name,'P3')));
% P4=find((strcmp(chanmap.Sensor_Name,'P4')));
% sensor_idx=[FP1,FP2,P3,P4];
% sensor_name={'FP1','FP2','P3','P4'};
% sensor_data=cell(1,4);
% patientid=cell(0,0);

%% run through files
for k=1:numel(files2load)%K%15%24
    idx=find(strcmp(DB2.Filename,files2load{k}));
    loadme=[DB2.Directory{idx(1)},'/', DB2.Filename{idx(1)}];
    dataStruct = load(loadme);
    dataCell = struct2cell(dataStruct);
    rawEEG = dataCell{1};
    clear dataStruct;
    clear dataCell;
    
    [nchan,nsamp]=size(rawEEG);
    tval=1/fSamp:1/fSamp:nsamp/fSamp;
    
    [~,L]=max(DB2.Amp(idx));

    for l=1:numel(idx)%L
        disp(['working on ' DB2.NUNAME{idx(l)}])
        
        idx1=DB2.trial_start(idx(l));
        idxend=DB2.trial_end(idx(l));
        
        tempEEG=rawEEG(:,idx1:idxend);
%         filtEEG=filtermat(b,tempEEG);
        filtEEG=filtermat2(tempEEG);
        
        remEEG=rembadchann(filtEEG);
        
        EEG=rereference(remEEG);
        
        % detect artifacts
        MV=movvar(EEG(15,:)',10);
        SNR=hpme(EEG(15,:))';
        [~,pks_MV]=findpeaks(MV, 'MinPeakDistance', crit_isi);
        [~,pks_SNR]=findpeaks(MV, 'MinPeakDistance', crit_isi);
        
        if numel(pks_MV)>numel(pks_SNR)
            pks_i=pks_MV;
        elseif numel(pks_MV)<numel(pks_SNR)
            pks_i=pks_SNR;
        else
            pks_i=pks_MV;
        end
        pks_i=pks_i(trialtrim:end-trialtrim+1);
        
        voltages=zeros(nchan,triallength,numel(pks_i)-1);
        for i=1:numel(pks_i)-1
           voltages(:,:,i)=EEG(:,pks_i(i)-eventloc+1:pks_i(i)-eventloc+triallength); 
        end
        
        [voltages_good,eff]=removebadtrials(voltages);
        
        meanvolts=mean(voltages_good,3);
        cd(meandir);
        save([DB.NUNAME{idx(l)},'.mat'],'meanvolts');
        cd(homedir);
    end
    
end


return
%% testbed
% qrcode(xy,mean(voltages,3),2,eventloc,[0,1])
% qrcode(xy,mean(voltages,3),2,eventloc,[-10,120],20,DB.NUNAME{idx(l)})
% qrcode(xy,mean(voltages_good,3),2,eventloc,[-10,120],20,DB.NUNAME{idx(l)},0)
% topobutter([80,120,150,250],meanvolts,xy,DB.NUNAME{idx(l)})

topogram([30,60,90,120]+50,meanvolts,xy,DB.NUNAME{idx(l)})

% Vprime=removebadtrials(voltages);
%% function database

function topogram(snaptime,meanvolts,xy,titname)

[X,Y] = meshgrid(-1:0.025:1, -1:0.025:1);
xq = X(:);
yq = Y(:);
rTest = sqrt(xq.^2 + yq.^2);
ckUnit = rTest <= 1;
xq = xq(ckUnit);
yq = yq(ckUnit);
nInt = length(xq);

[ntimes,nchan]=size(meanvolts);

EEG_int = zeros(5025, ntimes); 
% interpolate at each time time point per case
for jj = 1:ntimes
    EEG_int(:,jj) = griddata(xy(:,1), xy(:,2), ...
        meanvolts(1:end-1,jj), xq, yq);
end

f=figure('visible','off');
f.Position=[26 1 1267 500];
maxval=max(max(abs(EEG_int(:,snaptime))));
minmax=[min(min(meanvolts)),max(max(meanvolts))];
colorbds=[-maxval,maxval];

t=tiledlayout(1,numel(snaptime),'TileSpacing','tight');

for i=1:numel(snaptime)
%     subplot(1,numel(snaptime),i);
    nexttile
    scatter(xq,yq, 100, EEG_int(:,snaptime(i)), 'filled');
    
    cRB=redblue();
    colormap(cRB);
    caxis(colorbds);
    daspect([1,1,1])
    set(gca,'visible','off')
    if i==numel(snaptime)
       colorbar() 
    end
    %     title(num2str(snaptime(i)-50)+" ms")
    text(0,0.9,num2str(snaptime(i)-50)+" ms",'FontSize',20,'horizontalalignment','center')
    
    
end
title(t,titname,'FontSize',30,'Interpreter','None')
exportgraphics(f,[titname,'topo_hpf_only.png'],'Resolution',300);

end

function topobutter(snaptime,meanvolts,xy,titname)

[X,Y] = meshgrid(-1:0.025:1, -1:0.025:1);
xq = X(:);
yq = Y(:);
rTest = sqrt(xq.^2 + yq.^2);
ckUnit = rTest <= 1;
xq = xq(ckUnit);
yq = yq(ckUnit);
nInt = length(xq);

[ntimes,nchan]=size(meanvolts);

EEG_int = zeros(5025, ntimes); 
% interpolate at each time time point per case
for jj = 1:ntimes
    EEG_int(:,jj) = griddata(xy(:,1), xy(:,2), ...
        meanvolts(1:end-1,jj), xq, yq);
end

f=figure;
f.Position=[26 1 1267 984];
subplot(2,numel(snaptime),numel(snaptime)+1:2*numel(snaptime));
grid on;hold on;
% plot([0:ntimes-1]'-50,meanvolts');
plot(meanvolts');
xlim([0,max(snaptime)*1.2])
% get colorbounds
maxval=max(max(abs(EEG_int(:,snaptime))));
minmax=[min(min(meanvolts)),max(max(meanvolts))];
for i=1:numel(snaptime)
    plot([snaptime(i),snaptime(i)],minmax,'--','linewidth',3)
end
title(titname,'Interpreter','none');

colorbds=[-maxval,maxval];

for i=1:numel(snaptime)
    subplot(2,numel(snaptime),i);
    scatter(xq,yq, 100, EEG_int(:,snaptime(i)), 'filled');
    colorbar()
    cRB=redblue();
    colormap(cRB);
    caxis(colorbds);
    daspect([1,1,1])
    
    title(num2str(snaptime(i)-50)+" ms")
    
end
exportgraphics(f,[titname,'topo.png'],'Resolution',300);

end

function qrcode(xy,eegdata,edges,eventloc,bounds,eventblock,ID,normorno)
% get important info
[~,nsamp]=size(eegdata);
t_orig=-eventloc:nsamp-eventloc;


% find regions
dist=sqrt(xy(:,1).^2+xy(:,2).^2);
distthresh=0.3;
URi=find(xy(:,1)>0.0001 & xy(:,2)>0.0001 & dist>=distthresh);
LRi=find(xy(:,1)>0.0001 & xy(:,2)<-0.0001 & dist>=distthresh);
ULi=find(xy(:,1)<-0.0001 & xy(:,2)>0.0001 & dist>=distthresh);
LLi=find(xy(:,1)<-0.0001 & xy(:,2)<-0.0001 & dist>=distthresh);
Ci=find(dist<distthresh);
regionset={ULi,URi,Ci,LLi,LRi};

% get barcode values
barcodevals=zeros(numel(regionset),nsamp);

for i=1:nsamp
    
    indices2avg=i-edges:i+edges;
    indices2avg(indices2avg<1)=[];% remove negative vals
    indices2avg(indices2avg>nsamp)=[];% remove values out of range
    
    for k=1:numel(regionset)
       sett=eegdata(regionset{k},indices2avg);
       barcodevals(k,i)=mean(mean(((sett>0).*sett),1,'omitnan'),'omitnan');
    end
    
end

t_vals=bounds(1):bounds(2);

% normorno=1;
if normorno==1
    normbarcode=barcodevals./repmat(max(barcodevals')',1,nsamp);
    plotme=normbarcode(:,find(t_orig==bounds(1)):find(t_orig==bounds(2)));
else
    plotme=barcodevals(:,find(t_orig==bounds(1)):find(t_orig==bounds(2)));
    plotme(:,1:find(t_vals==eventblock))=0;
end

figure
h=heatmap(t_vals,{'UL','UR','C','LL','LR'},...
        plotme,'CellLabelColor','none');
idxrem=(rem(t_vals,20)==0);
h.ColorbarVisible = 'off';
h.GridVisible='off';
h.XDisplayLabels(~idxrem)={''};
% s=struct(h);
% s.XAxis.TickLabelRotation = 60;  % angled
colormap('hot');
set(gcf,'Position',[1 27 960 778]);
title(strrep(ID,'_',' '));
    
end

function Y=filtermat(b,X)

[nchan,nsamp]=size(X);
sampadd=ceil(nsamp/10);

Xpad=[fliplr(X(:,1:sampadd)),X,fliplr(X(:,end-sampadd:end))];
Ypad=zeros(size(Xpad));
for i=1:nchan
   
    Ypad(i,:)=filter(b,1,Xpad(i,:));
    
end

Y = Ypad(:,sampadd+1:end-sampadd-1);

end

function Y=filtermat2(X)

[nchan,nsamp]=size(X);
sampadd=ceil(nsamp/10);

Xpad=[fliplr(X(:,1:sampadd)),X,fliplr(X(:,end-sampadd:end))];
Ypad=zeros(size(Xpad));
for i=1:nchan
    Ypad(i,:)=highpass(Xpad(i,:),0.4, 1000, 'ImpulseResponse', 'iir', 'StopbandAttenuation', 40);
%     Ypad(i,:)=lowpass(temp,55, 1000, 'ImpulseResponse', 'iir', 'StopbandAttenuation', 40);
    
%     Ypad(i,:)=bandpass(Xpad(i,:),[0.4,50], 1000, 'ImpulseResponse', 'iir', 'StopbandAttenuation', 40);
    
    % gives same as original result... try making notch really steep?
%     =bandstop(temp,[55,65],1000,'ImpulseResponse','iir', 'StopbandAttenuation', 40);
    
end

Y = Ypad(:,sampadd+1:end-sampadd-1);

end

function y=hpme(x)
% add padding left and right
sampadd=ceil(numel(x)/10);
x=[fliplr(x(1:sampadd)), x, fliplr(x(end-sampadd:end))];
ytemp=highpass(x, 0.95*500, 1000);
y=ytemp(sampadd+1:end-sampadd-1);
% adds padding to left and right
end

function Xclean=rembadchann(X)
% load X in as row = nchan,col = samples, fs = 1000;
% https://eeglab.org/tutorials/06_RejectArtifacts/cleanrawdata.html

% (1) remove bad channels based off high standard deviation
S=std(X')';
pd=fitdist(S,'Normal');
idxrem=find(S>pd.sigma*5);
X(idxrem,:)=NaN;
Xclean=X;
% (2) remove bad channels based off poor correlation with others
% ***

end

function Xreref=rereference(X)

[nchan,~]=size(X);

avg_sig=mean(X,'omitnan');
M_avg_sig=repmat(avg_sig,nchan,1);
Xreref=X-M_avg_sig;

end

function [Vprime,eff]=removebadtrials(V)

%thresholds
thBad = 200;
thBlink = 90;
% thMvt = 40;

%eyechans
eyeChan = [18, 37, 226, 238, 241, 252];
window=80;% 80 ms moving average

[~,~,ntrial]=size(V);

trialrem=[];

for k=1:ntrial
    
    Vi=V(:,:,k);
    if sum(range(movmean(Vi,window,2),2) > thBad)>0
        trailrem=[trialrem;k];
    elseif sum(range(movmean(Vi(eyeChan,:),window,2),2) > thBlink)>0
        trailrem=[trialrem;k];
    end
    
end

Vprime=V;
Vprime(:,:,trialrem)=[];

[~,~,num]=size(Vprime);

eff=num/ntrial;

end
