
% The purpose of this script is (1) to specify time onsets of conditions
% for later decoding analysis: we care about onset times of odor
% expectations. Three conditions of interest are 3 expected odors.
% (2) to concatenate nuisance regressors across sessions and runs

%% set up
nvols = 430; %volumes per run
fprintf('Sorting onsets for decoding on %s\n',subno);
DecodingSubDir = fullfile(parentDir,'DecodingGLM',subno);

% set a new model directory
modelname = 'Outcome_pseudoconcat_Nz';
modeldir = fullfile(DecodingSubDir,'fxMultivariate',modelname);
if ~exist(modeldir,'dir')
    mkdir(modeldir)
else
end

% load behavioral data
if UseDay4(subj==Subs)==1 
    dat = dir(fullfile(subdir,'DAY4','comp*')); % dat from day4 
    load(fullfile(subdir,'DAY4',dat.name))
else
    dat = dir(fullfile(subdir,'DAY3','comp*')); % dat from day3
    load(fullfile(subdir,'DAY3',dat.name))
end

% load time correction mat file
load(fullfile(subdir,sprintf('AllTimeCorr_%s',subno)))

% name specifies reversal types
names = [];
name_tmp = [{'O1'} {'O2'} {'O3'}]; 
onsets = [];
durations = [];
clear NRs badvols 

rctr = 0; % count runs
odorctr = 0; % count reversal types (cumulative across and within runs)

cues_ons = []; % all cue will be in one column
sessdirs = dir(fullfile(subdir,'func','Day*'));

%% loop over sessions and runs

for ss = 1:nSess % sessions
    sess_name = sessdirs(ss).name; % current session name
    
for r = 1:nruns % runs
    run_name = sprintf('Run%d',r);    
    funcRunDir = fullfile(studydir,subno,'func',sess_name,run_name);
    fprintf('Run%d %s\n',r,sess_name)

    % find NR of current run
    NRdir = fullfile(SubResDir,'NuisanceRegressor');
    NRfile = fullfile(NRdir,sprintf('Nz_nuisance_regressors_%s_%s_%s.txt',subno,sess_name,run_name));

    NRtmp = load(NRfile);
    NRs{rctr+1} = NRtmp(:,1:34); % the first 34 NR, to concatenate across runs
    badvols{rctr+1} = NRtmp(:,35:size(NRtmp,2)); % other NR about bad volumes (number varies), empty if no bad vols at this run

    % time correction of current run in ms
    timecorr = alltimecorr(ss,r); 

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
    
    % find received odor labels
    Outcomes = d(:,2);
    
    for odor = 1:3
        odorctr = odorctr + 1; % revctr increases at each run so that we have separate conditions for each run
        onsets{odorctr} = (d(Outcomes==odor,11) + timecorr)./1000 + nvols * TR * rctr;
    end
    
    % cue presentation
    cues_ons = [cues_ons; ((d(:,7) + timecorr)./1000 + nvols * TR * rctr)];
        
    names = [names name_tmp]; % add odor condition names after each run
    rctr = rctr + 1;  % count up run idx

end % end of run loop
end % end of session loop

%%
% add names of other three conditions
names = [names {'cues'}];

% add other three onsets one by one after reversal trial onsets
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

% save onset mat file
save(fullfile(modeldir,'Onsets.mat'),'names','onsets','durations');

% write nuisance regressors
dlmwrite(fullfile(modeldir,'NuisanceRegressors.txt'),ConcatNR);
    

        
        