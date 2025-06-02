
% The purpose of this code is to analyze breathing data
% by sorting trials based on reversal locations and TMS conditions

%%
clear
clc
close all

studydir = '/Users/qingfangliu/Experiment';
ResDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ExptRes';
ScriptDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ForExperiment';
LabChartPath = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/LabChartData'; % directory saved all preprocessed labchart mat file of all subjects

addpath(ScriptDir)
parentDir = fileparts(ResDir); % the entire Analysis folder
addpath(parentDir)

SubInfoFile = fullfile(parentDir,'SubjectConds.XLSX');
SubInfo = readtable(SubInfoFile) ;
Subs = SubInfo.Sub; % subject index
Conds = SubInfo.Cond; % TMS condition
OdorDay1 = SubInfo.OdorDay1; % whether use day1 odor ratings
UseDay4 = SubInfo.UseDay4; % if use data from day3 and day4
Excluded = SubInfo.Excluded; % whether excluding from analysis

Subs = Subs(Excluded==0); % excluding subs from analysis
Conds = Conds(Excluded==0);

load(fullfile(parentDir,'Alld.mat'))  % load d matrix

nruns = 3; 
nSess = 2;

subctr = 0;

%%
% subj = 4;

for subj = Subs'
    
    %%
subctr = subctr + 1;

subno = sprintf('Sub%d',subj);

% to save processed lab chart output
LabChartdir = fullfile(studydir,subno,'ProcessedLC');

% load breathing data
load(fullfile(LabChartdir,sprintf('Agg_resp_data_%s.mat',subno)));

rctr = 0; % count runs

% also need to compute baseline breath signals for each session

Breath1 = [agg_resp_data{1,1};agg_resp_data{1,2};agg_resp_data{1,3}];
BreathBase(1) = mean(Breath1(:,2001:3001),'all');

Breath2 = [agg_resp_data{2,1};agg_resp_data{2,2};agg_resp_data{2,3}];
BreathBase(2) = mean(Breath2(:,2001:3001),'all');

%%
for ss = 1:nSess

for r = 1:nruns

    rctr = rctr + 1;  % count up run idx
    
    % load behavioral d matrix for current run
    d = Alld{subj==Subs,rctr};
    
    % resp data for current run
    resp_data = agg_resp_data{ss,r};
    
    % split d matrix to get rev loc indices
    cue1d = d(d(:,1)==1,:);
    cue2d = d(d(:,1)==2,:);
    
    rev_loc_1 = find(cue1d(:,3)==1);
    rev_loc_2 = find(cue2d(:,3)==1);
    
    cue1d(rev_loc_1,14) = 2;
    cue1d(rev_loc_1-1,14) = 1;
    cue1d(rev_loc_1+1,14) = 3;
    cue1d(rev_loc_1+2,14) = 4;
    
    cue2d(rev_loc_2,14) = 2;
    cue2d(rev_loc_2-1,14) = 1;
    cue2d(rev_loc_2+1,14) = 3;
    cue2d(rev_loc_2+2,14) = 4;
    
    % combine two d matrices back and sort by col 14
    Newd = [cue1d;cue2d];
    Newd = sortrows(Newd,7);
    
    if subj==1 && rctr==4
        Newd = [nan(4,14);Newd];
    end
    
    id1 = find(Newd(:,14)==1);
    id2 = find(Newd(:,14)==2);
    id3 = find(Newd(:,14)==3);
    id4 = find(Newd(:,14)==4);
    
    % correct baseline and save
    resp_id1(subctr,rctr,:,:) = resp_data(id1,:) - BreathBase(ss);
    resp_id2(subctr,rctr,:,:) = resp_data(id2,:) - BreathBase(ss);
    resp_id3(subctr,rctr,:,:) = resp_data(id3,:) - BreathBase(ss);
    resp_id4(subctr,rctr,:,:) = resp_data(id4,:) - BreathBase(ss);

end % end of run

end % end of session

end % end of subject

%%
Avg_resp_id1 = squeeze(nanmean(resp_id1,[2,3]));
M(:,1) = max(Avg_resp_id1,[],2,'omitnan');

Avg_resp_id2 = squeeze(nanmean(resp_id2,[2,3]));
M(:,2) = max(Avg_resp_id2,[],2,'omitnan');

Avg_resp_id3 = squeeze(nanmean(resp_id3,[2,3]));
M(:,3) = max(Avg_resp_id3,[],2,'omitnan');

Avg_resp_id4 = squeeze(nanmean(resp_id4,[2,3]));
M(:,4) = max(Avg_resp_id4,[],2,'omitnan');

figure
subplot(1,2,1)
boxplot(M)
title('Peak Amplitude')

Durations(:,1) = FindDuration(Avg_resp_id1);
Durations(:,2) = FindDuration(Avg_resp_id2);
Durations(:,3) = FindDuration(Avg_resp_id3);
Durations(:,4) = FindDuration(Avg_resp_id4);
subplot(1,2,2)
boxplot(Durations)
title('Inhale Duration')

%% function to output individual breathing durations

function Dur = FindDuration(m)
    % m is a matrix
    Dur = zeros(size(m,1),1);
    for i = 1:size(m,1)
        x = m(i,:);
        out = find(x<1); % use 1 instead of 0 as threshold to accomodate sub30
        out = out(out>4000);
        Dur(i) = (out(1) - 3000)/1000;
    end
end


