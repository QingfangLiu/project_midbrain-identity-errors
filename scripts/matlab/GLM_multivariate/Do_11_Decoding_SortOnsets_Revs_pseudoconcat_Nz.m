
% The purpose of this script is (1) to specify time onsets of conditions
% for later decoding analysis: we care about onset times of reversal trials, 
% separated by reversal types, also conditions of non-rev trials, cues, and
% responses. (2) to concatenate nuisance regressors across sessions and runs

%% set up
nvols = 430; %volumes per run
fprintf('Sorting onsets for decoding on %s\n',subno);
DecodingSubDir = fullfile(parentDir,'DecodingGLM',subno);

% set a new model directory
modelname = 'Revs_pseudoconcat_Nz';
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
name_tmp = [{'1->2'} {'1->3'} {'2->1'} {'2->3'} {'3->1'} {'3->2'}]; 
onsets = [];
durations = [];
clear NRs badvols 

rctr = 0; % count runs
revctr = 0; % count reversal types (cumulative across and within runs)

nonPE_ons = []; % all non reversal outcomes will be in one column
cues_ons = []; % all cue onsets will be in one column
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
        
    %sort reversals by type, numbers as such:
    %12 = 1 -> 2
    %13 = 1 -> 3
    %21 = 2 -> 1
    %23 = 2 -> 3
    %31 = 3 -> 1
    %33 = 3 -> 2

    rtypes = [12 13 21 23 31 32]; %reversal types, coded as below
    
    for i = 1:2 %loop through 2 cue types
        dtmp = d(d(:,1)==i,:); % d matrix of current cue type
        ridx = find(dtmp(:,3)==1); %find index of reversals for that cue type
        rtype = dtmp(ridx-1,2)*10 + dtmp(ridx,2); %code reversals according to above numbering system
        d((d(:,1)==i) & d(:,3),3) = rtype; % put coded reversal types back into d matrix - replacing col 3
    end

    for rev = rtypes
        revctr = revctr + 1; % revctr increases at each run so that we have separate conditions for each run
        onsets{revctr} = (d(d(:,3)==rev,11) + timecorr)./1000 + nvols * TR * rctr;
    end

    % non-reversal trials
    nonPE_ons = [nonPE_ons; ((d(d(:,3)==0,11) + timecorr)./1000 + nvols * TR * rctr)];
    % cue presentation
    cues_ons = [cues_ons; ((d(:,7) + timecorr)./1000 + nvols * TR * rctr)];
    % response
    %resp_ons = [resp_ons; ((d(~isnan(d(:,6)),8) + d(~isnan(d(:,6)),6) + timecorr)./1000 + nvols * TR * rctr)];
    
    names = [names name_tmp]; % add rev condition names after each run
    rctr = rctr + 1;  % count up run idx

end % end of run loop
end % end of session loop

%%
% add names of other three conditions
names = [names {'nonPE'} {'cues'}];

% add other three onsets one by one after reversal trial onsets
onsets{length(onsets)+1} = nonPE_ons;
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
    

        
        