

%% set up
fprintf('Sorting onsets for decoding on %s\n',subno);

% set a new model directory
modelname = 'STCues_LSS'; 

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

rctr = 0; % count runs
sessdirs = dir(fullfile(subdir,'func','Day*'));

%% loop over sessions and runs

for ss = 1:nSess % sessions
    sess_name = sessdirs(ss).name; % current session name
for r = 1:nruns % runs
    run_name = sprintf('Run%d',r);    
    funcRunDir = fullfile(studydir,subno,'func',sess_name,run_name);
    fprintf('Run%d %s\n',r,sess_name)
    rctr = rctr + 1;
   
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
        
    % all reversal trials
    Rev_ons = (d(d(:,3)==1,11) + timecorr)./1000;
    
    % non-reversal trials
    nonRev_ons =  (d(d(:,3)~=1,11) + timecorr)./1000;
        
    % cue presentation
    cues_ons = (d(:,7) + timecorr)./1000;
    
    % all condition names
    names = [{'Cue'} {'OtherCues'} {'nonRev'} {'Revs'}];
                
    for t = 1:size(d,1) % loop over number of trials
        
        % specify model directory (one GLM for each trial)
        modeldir = fullfile(DecodingSubDir,'fxMultivariate',modelname,sprintf('Run%dTrial%d',rctr,t));
        if ~exist(modeldir,'dir')
            mkdir(modeldir)
        end
        
        % current cue
        tempCue_ons = (d(t,7) + timecorr)./1000;
        
        % all other reversals
        OtherCues_ons = setdiff(cues_ons,tempCue_ons);
            
        onsets{1} = tempCue_ons;
        onsets{2} = OtherCues_ons;
        onsets{3} = nonRev_ons;
        onsets{4} = Rev_ons;

        clear durations
        % define durations 
        for i = 1:length(onsets)
            durations{i} = zeros(size(onsets{i}));
        end

        % save onset mat file
        save(fullfile(modeldir,'Onsets.mat'),'names','onsets','durations');

    end % loop of trials
end % loop of runs

end % end of session loop

        
        