
% This script takes part from the main code, only to count trial number
% gaps

%%
clear; clc; close all
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'; 
SubInfoDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/SubInfo';
maskdir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs';

%%
SubInfoFile = fullfile(SubInfoDir,'SubjectConds.XLSX');
SubInfo = readtable(SubInfoFile);

Subs = SubInfo.Sub;                     % subject index
Conds = SubInfo.Cond;                   % TMS condition
OdorDay1 = SubInfo.OdorDay1;            % whether use day1 odor ratings
UseDay4 = SubInfo.UseDay4;              % if use data from day3 and day4
Excluded = SubInfo.Excluded;            % whether excluding from analysis

Subs = Subs(Excluded==0);               % excluding subs from analysis
Conds = Conds(Excluded==0);
UseDay4 = UseDay4(Excluded==0);
SubInfo = SubInfo(Excluded==0,:);
nsubs = length(Subs);

nvols = 430;                            % volume number per run
nruns = 3; 
nSess = 2;
nallruns = nruns * nSess;

% load beh data Alld.mat
load(fullfile(ResDir,'Behavior','Alld.mat'))

%% PATHS AND VARIABLES TO CHANGE

nrevs_per_run = 12;   % # of revs per run
nlocs_per_rev = 4;    % # of locs per rev

subctr = 0;
all_dist = zeros(nsubs,nallruns,2);

for subj = Subs'
    subno = sprintf('Sub%d',subj);
    fprintf('Working on %s\n',subno)
    subctr = subctr + 1;
    rctr = 0;                           % count number of runs
    for ss = 1:nSess 
    for r = 1:nruns    
        rctr = rctr + 1;
        fprintf('Data organization run %d\n',rctr)
        d = Alld{subj==Subs,rctr};
        ntrials = size(d,1);
        d = [d(:,[1:3,7]),(1:ntrials)'];
        id_mat = zeros(nrevs_per_run,nlocs_per_rev);                        % matrix of trial id for this run: 12 reversals * 4 rev locs
        for cue = 1:2
            sep_d = d(d(:,1)==cue,:);
            revloc = find(sep_d(:,3)==1); 
            id_mat((1:6) + 6*(cue-1),2) = sep_d(revloc,4);    % cue-onset of each trial
            id_mat((1:6) + 6*(cue-1),1) = sep_d(revloc-1,4);
            id_mat((1:6) + 6*(cue-1),3) = sep_d(revloc+1,4);
            id_mat((1:6) + 6*(cue-1),4) = sep_d(revloc+2,4);
        end
        
        % distance between trials pre and post reversal
        pre_dist = mean(id_mat(:,2) - id_mat(:,1));
        post_dist = mean(id_mat(:,3) - id_mat(:,2));
        
        all_dist(subctr,rctr,1) = pre_dist;
        all_dist(subctr,rctr,2) = post_dist;
        
    end % run
    end % sess
end


%% average trial distance 

% across all sham and cTBS runs
avg_dist = squeeze(mean(all_dist,2));
mean(avg_dist)
std(avg_dist)
boxplot(avg_dist)
[h, p] = ttest(avg_dist(:, 1), avg_dist(:, 2));

%%
% across all sham runs
all_dist_sham = zeros(nsubs, nruns, 2);
for i = 1:nsubs
    if strcmp(Conds{i}, 'PA')
        all_dist_sham(i, :, :) = all_dist(i, 1:3, :);
    else
        all_dist_sham(i, :, :) = all_dist(i, 4:6, :);
    end
end

avg_dist_sham = squeeze(mean(all_dist_sham,2));
mean(avg_dist_sham)
std(avg_dist_sham)
boxplot(avg_dist_sham)
[h_sham, p_sham] = ttest(avg_dist_sham(:, 1), avg_dist_sham(:, 2));

%%
% across all cTBS runs
all_dist_cTBS = zeros(nsubs, nruns, 2);
for i = 1:nsubs
    if strcmp(Conds{i}, 'AP')
        all_dist_cTBS(i, :, :) = all_dist(i, 1:3, :);
    else
        all_dist_cTBS(i, :, :) = all_dist(i, 4:6, :);
    end
end

avg_dist_cTBS = squeeze(mean(all_dist_cTBS,2));
mean(avg_dist_cTBS)
std(avg_dist_cTBS)
boxplot(avg_dist_cTBS)
[h_cTBS, p_cTBS] = ttest(avg_dist_cTBS(:, 1), avg_dist_cTBS(:, 2));

%%
vec_avg_dist = [avg_dist_sham(:);avg_dist_cTBS(:)];

% factor 1: sham or cTBS
TMS = sort(repmat(1:2,1,62));
% factor 2: pre or post
PrePost = [repmat(1,1,31),repmat(2,1,31),repmat(1,1,31),repmat(2,1,31)];

table_data = table(repmat(Subs,4,1),TMS',PrePost',vec_avg_dist, ...
    'VariableNames', {'Subs','TMS','PrePost','Delay'});
file_name = fullfile(ResDir,'MVPA_Searchlight','TrialOnsetDelay.xlsx');
writetable(table_data, file_name, 'WriteRowNames', false);


