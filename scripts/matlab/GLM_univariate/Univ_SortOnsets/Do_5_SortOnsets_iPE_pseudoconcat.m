
% This script is used to find the time onsets of task events aligned with
% the scan timing.
% The events include: cue onset, response onset, reversal prediction,
% all non-reversal prediction.

% This is for a pseudoconcatenate model.
% To estimate two regressors of interest (reversal, non-reversal) for each
% run, but concatenate other two regressors (cues, response), and 34
% nuisance regressors.

%% initialize things

nvols = 430; % number of volumes each run
fprintf('Working on sorting onsets iPE_concat for %s\n',subno);

% set up model directory
modelname = 'iPE_pseudoconcat';
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

% load behavioral mat file
load(fullfile(bfile.folder,bfile.name));

sessdirs = dir(fullfile(studydir,subno,'func','Day*'));

clear tmpons durations NRs badvols onsets

% condition names
names = [];
name_tmp = [{'Outcomes_Reversal'} {'Outcomes_NonRev'}]; 

for i = 1:12 % 2 for each run, 6 runs in total
    onsets{i} = [];
end

rctr = 0; % count runs
cues_ons = []; % all cue onsets will be in one column
resp_ons = []; % all response onsets will be in one column

%%
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
    
    tmpons{1} = (d(d(:,3)==1,11) + timecorr)./1000; % reversal
    tmpons{2} = (d(d(:,3)==0,11) + timecorr)./1000; % non-reversal
    tmpons{3} = (d(:,7) + timecorr)./1000; % cue presentation timing
    tmpons{4} = (d(~isnan(d(:,6)),8) + d(~isnan(d(:,6)),6) + timecorr)./1000; % response timing

    % adjust onset timings for runs from the second 
    for i = 1:length(tmpons)
        tmpons{i} = tmpons{i} + nvols * TR * rctr;
    end
            
    % keep reversal and non-reversal onsets for current run
    onsets{rctr*2+1} = tmpons{1};  % reversal onsets
    onsets{rctr*2+2} = tmpons{2};  % non-reversal onsets
    
    % concatenate cues and resps onsets across all runs
    cues_ons = [cues_ons; tmpons{3}];
    resp_ons = [resp_ons; tmpons{4}];
    
    names = [names name_tmp]; % add condition names after each run
    rctr = rctr + 1;  % count up run idx

end  % end of the run
end  % end of session

% add names of other two conditions
names = [names {'cues'} {'resp'}];

% add other two onsets one by one after reversal trial onsets
onsets{length(onsets)+1} = cues_ons;
onsets{length(onsets)+1} = resp_ons;

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

clear durations onsets names cues_ons resp_ons
    
        
        