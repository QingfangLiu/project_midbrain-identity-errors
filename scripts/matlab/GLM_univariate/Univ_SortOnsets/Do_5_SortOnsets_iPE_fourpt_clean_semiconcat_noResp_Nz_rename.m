
% This script is used to find the time onsets of task events aligned with
% the scan timing.
% The events include: cue onset, response onset, reversal prediction,
% non-reversal prediction (-1,+1,+2).

%% initialize things
nvols = 430; % number of volumes each run

% set up model directory
modelname = 'iPE_fourpt_clean_semiconcat_noResp_Nz_rename';
fprintf('Working on sorting onsets of %s for %s\n',modelname,subno);
modeldir = fullfile(subdir,'fxUnivariate',modelname);
if ~exist(modeldir,'dir')
    mkdir(modeldir)
else
end

% load time correction mat file
load(fullfile(subdir,sprintf('AllTimeCorr_%s',subno)))

% load complete behavioral data
if UseDay4(subj==Subs)==0 % for all subjects not involving day4
    bfile = dir(fullfile(subdir,'DAY3','complete*'));
else 
    bfile = dir(fullfile(subdir,'DAY4','complete*'));
end

load(fullfile(bfile.folder,bfile.name));

sessdirs = dir(fullfile(studydir,subno,'func','Day*'));

clear tmpons durations NRs badvols onsets

% condition names
%names = [];
%name_tmp = [{'Outcomes_Reversal'} {'Outcomes_Minus1'} {'Outcomes_Plus1'} {'Outcomes_Plus2'}]; 

for i = 1:8 % 4 for each session, 2 sessions in total
    onsets{i} = [];
end

rctr = 0; % count runs
ssctr = 0; % count sessions
cues_ons = []; % all cue onsets will be in one column

%% loop over sessions and runs
for ss = 1:nSess
    sess_name = sessdirs(ss).name;

for r = 1:nruns

    fprintf('Run%d %s\n',r,sess_name)
    run_name = sprintf('Run%d',r);   
    rundir = fullfile(subdir,'func',sess_name,sprintf('Run%d',r));
    
    % find NR of current run
    NRdir = fullfile(SubResDir,'NuisanceRegressor');
    NRfile = fullfile(NRdir,sprintf('Nz_nuisance_regressors_%s_%s_%s.txt',subno,sess_name,run_name));
    
    NRtmp = load(NRfile);
    NRs{rctr+1} = NRtmp(:,1:34); % the first 34 NR, to concatenate across runs
    badvols{rctr+1} = NRtmp(:,35:size(NRtmp,2)); % other NR about bad volumes (number varies), empty if no bad vols at this run
    

    if UseDay4(subj==Subs)==0 % for all subjects except for sub1 and sub35
        if ss==1
             d = res.reversal_learning_task_DAY2{r}; % load behavioral data from day2
        else
             d = res.reversal_learning_task_DAY3{r}; % load behavioral data from day3
        end
    end

    if UseDay4(subj==Subs)==1 % set an exception for sub1 and sub35, when using day4
         if ss==1
            d = res.reversal_learning_task_DAY3{r}; % load behavioral data from day3
         else
            d = res.reversal_learning_task_DAY4{r}; % load behavioral data from day4
         end
    end

    % handle sub1-day4-run1 differently by removing first four trials
    % (trials before scan starts)
    if subj==1 && ss==2 && r==1
        d = d(5:end,:);
    end
    
    timecorr = alltimecorr(ss,r); % time adjust of current run in ms
    
    choice_acc = d(:,2)==d(:,13); % find acc of each trial
    
    tmpons{1} = (d(d(:,3)==1,11) + timecorr)./1000; % reversal
    tmpons{5} = (d(:,7) + timecorr)./1000; % cue presentation timing
    
    % need some difference in defining the regressors
    cue1d = d(d(:,1)==1,:);
    cue2d = d(d(:,1)==2,:);
    
    choice_acc_1 = choice_acc(d(:,1)==1);
    choice_acc_2 = choice_acc(d(:,1)==2);

    rev_loc_1 = find(cue1d(:,3)==1);
    rev_loc_2 = find(cue2d(:,3)==1);
    
    Minus1_outcomes_cue1 = (cue1d(rev_loc_1 - 1,11) + timecorr)./1000;
    Minus1_outcomes_cue1 = Minus1_outcomes_cue1(choice_acc_1(rev_loc_1 - 1)); % correct acc only
    
    Plus1_outcomes_cue1 = (cue1d(rev_loc_1 + 1,11) + timecorr)./1000;
    Plus1_outcomes_cue1 = Plus1_outcomes_cue1(choice_acc_1(rev_loc_1 + 1)); % correct acc only
    
    Plus2_outcomes_cue1 = (cue1d(rev_loc_1 + 2,11) + timecorr)./1000;
    Plus2_outcomes_cue1 = Plus2_outcomes_cue1(choice_acc_1(rev_loc_1 + 2)); % correct acc only
    
    Minus1_outcomes_cue2 = (cue2d(rev_loc_2 - 1,11) + timecorr)./1000;
    Minus1_outcomes_cue2 = Minus1_outcomes_cue2(choice_acc_2(rev_loc_2 - 1)); % correct acc only
    
    Plus1_outcomes_cue2 = (cue2d(rev_loc_2 + 1,11) + timecorr)./1000;
    Plus1_outcomes_cue2 = Plus1_outcomes_cue2(choice_acc_2(rev_loc_2 + 1)); % correct acc only
    
    Plus2_outcomes_cue2 = (cue2d(rev_loc_2 + 2,11) + timecorr)./1000;
    Plus2_outcomes_cue2 = Plus2_outcomes_cue2(choice_acc_2(rev_loc_2 + 2)); % correct acc only

    tmpons{2} = [Minus1_outcomes_cue1;Minus1_outcomes_cue2];
    tmpons{3} = [Plus1_outcomes_cue1;Plus1_outcomes_cue2];
    tmpons{4} = [Plus2_outcomes_cue1;Plus2_outcomes_cue2];

    % adjust onset timings for runs from the second 
    for i = 1:length(tmpons)
        tmpons{i} = tmpons{i} + nvols * TR * rctr;
    end
            
    % keep reversal and non-reversal onsets for current run
    onsets{ssctr*4+1} = [onsets{ssctr*4+1}; tmpons{1}];  % reversal onsets
    onsets{ssctr*4+2} = [onsets{ssctr*4+2}; tmpons{2}];  % non-reversal minus1 onsets
    onsets{ssctr*4+3} = [onsets{ssctr*4+3}; tmpons{3}];  % non-reversal plus1 onsets
    onsets{ssctr*4+4} = [onsets{ssctr*4+4}; tmpons{4}];  % non-reversal plus2 onsets
    
    % concatenate cues and resps onsets across all runs
    cues_ons = [cues_ons; tmpons{5}];    
    rctr = rctr + 1;  % count up run idx
    
end 

    %names = [names name_tmp]; % add condition names after each run
    ssctr = ssctr + 1; % count up session idx
    fprintf('\n');
end


% add names of other two conditions
% names = [names {'cues'}];

names = [{'Outcomes_Reversal_S1'} {'Outcomes_Minus1_S1'} {'Outcomes_Plus1_S1'} {'Outcomes_Plus2_S1'}...
    {'Outcomes_Reversal_S2'} {'Outcomes_Minus1_S2'} {'Outcomes_Plus1_S2'} {'Outcomes_Plus2_S2'}...
    {'cues'}]; 

% add other two onsets one by one after reversal trial onsets
onsets{length(onsets)+1} = cues_ons;

% define durations 
for i = 1:length(onsets)
    durations{i} = zeros(size(onsets{i}));
end

% concatenate all first 34 NRs across 6 runs
ConcatNR = vertcat(NRs{:});

% concatenate bad volume nuisance regressors diagonally (revised from
% previous lab code)
allbvols = blkdiag(badvols{:});
ConcatNR = [ConcatNR allbvols];

%% save and write

save(fullfile(modeldir,'Onsets.mat'),'names','onsets','durations');
dlmwrite(fullfile(modeldir,'NuisanceRegressors.txt'),ConcatNR);

clear durations onsets names cues_ons


 
            
        
        