
% This script is used to find the time onsets of task events aligned with
% the scan timing.
% The events include: cue onset, response onset, reversal prediction,
% non-reversal prediction (-1,+1,+2).

% This script was modifed based on the iPE_fourpt code to be a concatenate
% model. 

%% initialize things
nvols = 430; % number of volumes each run
fprintf('Working on sorting onsets iPE fourpt concat for %s\n',subno);

% set up model directory
modelname = 'iPE_fourpt_concat';
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

clear onsets names tmpons durations

% condition names
names{1} = 'Outcomes_Reversal';
names{2} = 'Outcomes_Minus1';
names{3} = 'Outcomes_Plus1';
names{4} = 'Outcomes_Plus2';
names{5} = 'Cues';
names{6} = 'Responses';
    
for i = 1:length(names)
    onsets{i} = [];
end

rctr = 0; % count runs

%% loop over sessions and runs
for ss = 1:nSess
    
    sess_name = sessdirs(ss).name;

for r = 1:nruns

    run_name = sprintf('Run%d',r);   
    fprintf('Run%d %s\n',r,sess_name)
    rundir = fullfile(subdir,'func',sess_name,sprintf('Run%d',r));
    
    % find NR of current run
    NRdir = fullfile(subdir,'NuisanceRegressor');
    NRfile = fullfile(NRdir,sprintf('nuisance_regressors_%s_%s_%s.txt',subno,sess_name,run_name));
    
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
    
    iPE_outcomes = (d(d(:,3)==1,11) + timecorr)./1000;
    cues = (d(:,7) + timecorr)./1000;
    responses = (d(~isnan(d(:,6)),8) + d(~isnan(d(:,6)),6) + timecorr)./1000;

    % need some difference in defining the regressors
    cue1d = d(d(:,1)==1,:);
    cue2d = d(d(:,1)==2,:);

    rev_loc_1 = find(cue1d(:,3)==1);
    Minus1_outcomes_cue1 = (cue1d(rev_loc_1 - 1,11) + timecorr)./1000;
    Plus1_outcomes_cue1 = (cue1d(rev_loc_1 + 1,11) + timecorr)./1000;
    Plus2_outcomes_cue1 = (cue1d(rev_loc_1 + 2,11) + timecorr)./1000;

    rev_loc_2 = find(cue2d(:,3)==1);
    Minus1_outcomes_cue2 = (cue2d(rev_loc_2 - 1,11) + timecorr)./1000;
    Plus1_outcomes_cue2 = (cue2d(rev_loc_2 + 1,11) + timecorr)./1000;
    Plus2_outcomes_cue2 = (cue2d(rev_loc_2 + 2,11) + timecorr)./1000;

    tmpons{1} = iPE_outcomes;
    tmpons{2} = [Minus1_outcomes_cue1;Minus1_outcomes_cue2]; 
    tmpons{3} = [Plus1_outcomes_cue1;Plus1_outcomes_cue2];
    tmpons{4} = [Plus2_outcomes_cue1;Plus2_outcomes_cue2];
    tmpons{5} = cues;
    tmpons{6} = responses;
    
    % adjust onset timings for runs from the second 
    for i = 1:length(tmpons)
        tmpons{i} = tmpons{i} + nvols * TR * rctr;
    end
    
    rctr = rctr + 1;  % count up run idx
    
    % concat at each new run
    onsets{1} = [onsets{1}; tmpons{1}]; 
    onsets{2} = [onsets{2}; tmpons{2}];
    onsets{3} = [onsets{3}; tmpons{3}];
    onsets{4} = [onsets{4}; tmpons{4}];
    onsets{5} = [onsets{5}; tmpons{5}];
    onsets{6} = [onsets{6}; tmpons{6}];
        
end 

fprintf('\n');

end
 
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



        
        