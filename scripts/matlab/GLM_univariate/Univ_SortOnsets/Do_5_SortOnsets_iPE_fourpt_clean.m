
% This script is used to find the time onsets of task events aligned with
% the scan timing.
% The events include: cue onset, response onset, reversal prediction,
% non-reversal prediction (-1,+1,+2).

%% initialize things
fprintf('Working on sorting onsets iPE fourpt clean for %s\n',subno);

% set up model directory
modelname = 'iPE_fourpt_clean';
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

%% loop over sessions and runs
for ss = 1:nSess
    sess_name = sessdirs(ss).name;

for r = 1:nruns

    fprintf('Run%d %s\n',r,sess_name)
    rundir = fullfile(subdir,'func',sess_name,sprintf('Run%d',r));

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
    
    iPE_outcomes = (d(d(:,3)==1 & ~choice_acc,11) + timecorr)./1000;  % wrong acc only
    cues = (d(:,7) + timecorr)./1000;
    responses = (d(~isnan(d(:,6)),8) + d(~isnan(d(:,6)),6) + timecorr)./1000;

    % need some difference in defining the regressors
    cue1d = d(d(:,1)==1,:);
    cue2d = d(d(:,1)==2,:);

    choice_acc_1 = choice_acc(d(:,1)==1);
    choice_acc_2 = choice_acc(d(:,1)==2);
     
    % loc of reversal for each cue
    rev_loc_1 = find(cue1d(:,3)==1);
    rev_loc_2 = find(cue2d(:,3)==1);
    
    Minus1_outcomes_cue1 = (cue1d(rev_loc_1 - 1,11) + timecorr)./1000;
    Minus1_outcomes_cue1 = Minus1_outcomes_cue1(choice_acc_1(rev_loc_1 - 1)); % correct acc only
    
    Plus1_outcomes_cue1 = (cue1d(rev_loc_1 + 1,11) + timecorr)./1000;
    Plus1_outcomes_cue1 = Plus1_outcomes_cue1(choice_acc_1(rev_loc_1 + 1)); % correct acc only
    
    Plus2_outcomes_cue1 = (cue1d(rev_loc_1 + 2,11) + timecorr)./1000;
    Plus2_outcomes_cue1 = Plus2_outcomes_cue1(choice_acc_1(rev_loc_1 + 2)); % correct acc only

    Minus1_outcomes_cue2 = (cue1d(rev_loc_2 - 1,11) + timecorr)./1000;
    Minus1_outcomes_cue2 = Minus1_outcomes_cue2(choice_acc_2(rev_loc_2 - 1)); % correct acc only
        
    Plus1_outcomes_cue2 = (cue1d(rev_loc_2 + 1,11) + timecorr)./1000;
    Plus1_outcomes_cue2 = Plus1_outcomes_cue2(choice_acc_2(rev_loc_2 + 1)); % correct acc only
       
    Plus2_outcomes_cue2 = (cue1d(rev_loc_2 + 2,11) + timecorr)./1000;
    Plus2_outcomes_cue2 = Plus2_outcomes_cue2(choice_acc_2(rev_loc_2 + 2)); % correct acc only

    Minus1_outcomes = [Minus1_outcomes_cue1;Minus1_outcomes_cue2];
    Plus1_outcomes = [Plus1_outcomes_cue1;Plus1_outcomes_cue2];
    Plus2_outcomes = [Plus2_outcomes_cue1;Plus2_outcomes_cue2];

    names{1} = 'Outcomes_Reversal';
    names{2} = 'Outcomes_Minus1';
    names{3} = 'Outcomes_Plus1';
    names{4} = 'Outcomes_Plus2';
    names{5} = 'Cues';
    names{6} = 'Responses';

    onsets{1} = iPE_outcomes;
    onsets{2} = Minus1_outcomes; 
    onsets{3} = Plus1_outcomes; 
    onsets{4} = Plus2_outcomes; 
    onsets{5} = cues;
    onsets{6} = responses;

    for i = 1:length(onsets)
        durations{i} = zeros(size(onsets{i}));
    end

    save(fullfile(modeldir,sprintf('%s_Run%d_Onsets.mat',sess_name,r)),'names','onsets','durations');

end 

fprintf('\n');

end
 
            
        
        