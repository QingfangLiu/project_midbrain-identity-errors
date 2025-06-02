
% This script organizes the behavioral data matrix d into 'Alld.mat' and
% saves it. 
% Later analysis can simply load this mat file

%%
clear; clc; close all
studydir = '/Users/qingfangliu/Experiment';
ResDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ExptRes';
ScriptDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ForExperiment';

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
UseDay4 = UseDay4(Excluded==0);

nSubs = length(Subs); % number of subjects

nruns = 3; 
nSess = 2;
nallruns = nruns * nSess;

subctr = 0;
Alld = cell(nSubs,nallruns);

%%
for subj = Subs'

subno = sprintf('Sub%d',subj);
subdir = fullfile(studydir,subno);
fprintf('%s\n',subno)
subctr = subctr + 1;

% prepare to organize betas into runs
data = cell(nallruns,1);
OdorLabels = cell(nallruns,1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load behavioral data
if UseDay4(subj==Subs)==1 
    dat = dir(fullfile(subdir,'DAY4','comp*')); % dat from day4 
    load(fullfile(subdir,'DAY4',dat.name))
else
    dat = dir(fullfile(subdir,'DAY3','comp*')); % dat from day3
    load(fullfile(subdir,'DAY3',dat.name))
end
sessdirs = dir(fullfile(subdir,'func','Day*'));


%% organize data for each run

rctr = 0;
for ss = 1:nSess % sessions
    sess_name = sessdirs(ss).name; % current session name
    
for r = 1:nruns % runs
    run_name = sprintf('Run%d',r);    
    funcRunDir = fullfile(studydir,subno,'func',sess_name,run_name);
    fprintf('Run%d %s\n',r,sess_name)
    rctr = rctr + 1;
    
    % load behavioral data   
    if ss == 1 && UseDay4(subj==Subs)==0
       d = res.reversal_learning_task_DAY2{1,r};
    elseif ss == 2 && UseDay4(subj==Subs)==0
       d = res.reversal_learning_task_DAY3{1,r};
    elseif ss == 1 && UseDay4(subj==Subs)==1 
       d = res.reversal_learning_task_DAY3{1,r};
    elseif ss == 2 && UseDay4(subj==Subs)==1 
       d = res.reversal_learning_task_DAY4{1,r};
    end
    
    % handle sub1-day4-run1 differently by removing first four trials
    % (trials before scan starts)
    if subj==1 && ss==2 && r==1
        d = d(5:end,:);
    end
    
    ntrials = size(d,1);
    Alld{subctr,rctr} = d;
    
end
end

end

%%
% save all behavioral data in a 31*6 cell 
% each entry is a matrix of 64 (or 60) * 13
% only for sub1 run4, only 60 trials were useable. 

dname = fullfile(parentDir,'Alld.mat');
save(dname,'Alld')



