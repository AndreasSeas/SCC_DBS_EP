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
addpath(genpath('/Users/as822/Box/PST Project'));% add the path of all eeg 
addpath(genpath(homedir));% add the path of all folders here
raweegdir='/Users/as822/Box/SCC_DBS_Sinai/Data/***';


%% TEMP - definition of run info for testing (before function)
bool_plot=1;
bool_save=0;

%% constants for all patients
fSamp=1000;% Hz
fStim=2; % Hz
fNyquist = fSamp / 2; % Hz

%% load metadata
% load ID database to go through
IDDB=readtable('ID_database.xlsx');
IDDB_=IDDB;
xy = load('adultAvg_xy256_unitCircle.txt');% load 2d projection
rSens = load('adultAvg_xyz256_unitSphere.txt');% load 3d projection
rInt = load('unitSphere_1000pts.txt');% load 1000 point interpolation

%% filter IDDB to only the types you want (here doing all)
idxrun=find(strcmp(IDDB.Type,'OT') |...
    strcmp(IDDB.Type,'CB') | ...
    strcmp(IDDB.Type,'FM'));
IDDB=IDDB(idxrun,:);
[numruns,~]=size(IDDB);

%% run through each dataset, getting the reference channel
refchann=cell(numruns,1);% all have fSamp = 1000 Hz = 1kHz
nChan=zeros(size(refchann));
nSamp=zeros(size(refchann));
for k=1:numruns
%     chro
    loadme=[IDDB.Directory{k},'/', IDDB.Filename{k}];
    dataStruct = load(loadme);
    dataCell = struct2cell(dataStruct);
    rawEEG = dataCell{1};
    [nChan(k), nSamp(k)] = size(rawEEG);
    refchann{k} = rawEEG(15,:);
    
    disp(['done with ', IDDB.Filename{k}]);
    clear dataStruct;
    clear dataCell;
end
return
%% filter and plot reference channels
bpref=refchann;


for k=1:numruns
    figure; hold on; grid on;
    x=refchann{k};
    tval=1/fSamp:1/fSamp:numel(x)/fSamp;
%     plot(abs(diff(x))./max(abs(diff(x))));
    gtmp = input('');
%     [x_,lowidx,hiidx]=transientremoval(x);
%     figure; hold on; grid on;
%     plot(tval,refchann{k});
%     plot(tval(lowidx:hiidx),x_);
%     title(IDDB.Filename{k})
%    % y = bandpass(x,fpass,fs)
%     bpref{k}=bandpass(refchann{k},[0.4,50],fSamp);
%     plot(bpref{k});
    
%     pause(0.5)
end

% remember, the role of this step is to basically go through each dataset
% (in this example) and find the epoch start/stop ubiquitously... batch
% processing rather than serial processing for now
%   Delta = 1 - 4 Hz
%   Theta = 4 - 8 Hz
%   Alpha = 8 - 12 Hz
%   Beta  = 12 - 30 Hz
%   Gamma = 30 - 70, 70 - 150 Hz

%% function database
function [x_,lowidx,hiidx]=transientremoval(x)
dx=abs(diff(x));
midpt=round(numel(x)/2);% find midpoint
idx_hidx=find(dx>0.5*max(dx));
lowidx=max(idx_hidx(idx_hidx<midpt));
hiidx=min(idx_hidx(idx_hidx>midpt));

if numel(lowidx)==0
    lowidx=1;
end

if numel(hiidx)==0
    hiidx=numel(x);
end

lowidx=lowidx+5;
hiidx=hiidx-5;

x_=x(lowidx:hiidx);
end
