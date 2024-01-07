function [DB,V,xy,chanmap,homedir]=setmeup()
%(independent of your naming of directories)
disp("Set Up Workspace");
homedir=pwd;
addpath(genpath(homedir));% add the path of all folders here
% cd('..');
cd('C:\Users\as822\Box\collab-Sinai\project-selectEPs\EEG-dataSumm\eegdata_mean');
meandir=pwd;
cd(homedir);
cd('C:\Users\as822\Box\collab-Sinai\Code\EEG Processing Pipeline\referencefiles');
refdir=pwd;
addpath(genpath(refdir))
cd(homedir);

%% constants for all patients
eventloc=50;
triallength=501;

%% load metadata
% raw=readtable('IndividualRunDB_220523.xlsx');
raw=readtable('IndividualRunDB_220617(adjusteddescriptions).xlsx');
% raw=readtable('IndividualRunDB_220617_modifiedlabeledtree.xlsx');
xy = load('adultAvg_xy256_unitCircle.txt');% load 2d projection
rSens = load('adultAvg_xyz256_unitSphere.txt');% load 3d projection
rInt = load('unitSphere_1000pts.txt');% load 1000 point interpolation
chanmap=readtable('net256_channelMap_edited.xlsx');

%% get only the specific trials you want
F=zeros(size(raw.ID));
F=F | strcmp(raw.Type_real,"CB");% get CB ones
F=F | strcmp(raw.Type_real,"FM");% get FM ones
F=F | strcmp(raw.Type_real,"OT");% get OT ones
% F=F & raw.AmpIDX==4;% get max amplitude only
DB=removevars(raw(F,:),{'NOTES'});
DB(any(ismissing(DB),2), :) = [];% only get the ones with all the data we need

%% create the db for the ones you want to test
[L,~]=size(DB);
cd(meandir)
for i=1:L
    load(DB.NUNAME{i}+".mat");
    V(:,:,i)=meanvolts;
end
cd(homedir)
end