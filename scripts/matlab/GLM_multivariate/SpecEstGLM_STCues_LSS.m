
%% setup
modelname = 'STCues_LSS';
fprintf('Specify and estimate GLM for multivariate %s on %s\n',modelname,subno);

%load behavioral data
if UseDay4(subj==Subs)==1 
    dat = dir(fullfile(subdir,'DAY4','comp*')); % dat from day4 
    load(fullfile(subdir,'DAY4',dat.name))
else
    dat = dir(fullfile(subdir,'DAY3','comp*')); % dat from day3
    load(fullfile(subdir,'DAY3',dat.name))
end

subdirfunc = fullfile(studydir,subno,'func');
sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));
rctr = 0;

%%
for ss = 1:nSess 
    
    sess_name = sessdirsfunc(ss).name;
    
for r = 1:nruns
    rctr = rctr + 1;
    fprintf('Running run %d/6 %s\n',rctr,subno)
    
    run_name = sprintf('Run%d',r); 
    scandir = fullfile(subdir,'func',sess_name,sprintf('Run%d',r));
    filenames = dir(fullfile(scandir,'s2mmwf*.nii')); % using the 2mm images for decoding analysis (smaller smoothing size)
    scans = [];
    for i = 1:length(filenames)
        scans = [scans; {fullfile(scandir,sprintf('%s,1',filenames(i).name))}];
    end
        
    % find NR of current run
    NRdir = fullfile(subdir,'NuisanceRegressor');
    regfile = fullfile(NRdir,sprintf('nuisance_regressors_%s_%s_%s.txt',subno,sess_name,run_name));

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
    
    for t = 1:size(d,1)
                
        modeldir = fullfile(DecodingSubDir,'fxMultivariate',modelname,sprintf('Run%dTrial%d',rctr,t));
        onsfile = fullfile(modeldir,'Onsets.mat');

        clear matlabbatch
        matlabbatch{1}.spm.stats.fmri_spec.dir = {modeldir};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

        matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = scans;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {onsfile};
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {regfile};
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;

        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.1;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

        % Estimate the model
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(modeldir,'SPM.mat')};
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

        % Run model
        spm_jobman('run', matlabbatch)
        clear matlabbatch
        
    end
end    
end

