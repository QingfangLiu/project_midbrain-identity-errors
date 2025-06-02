

%% setup
fprintf('Specifying and estimating GLM for %s\n',subno);
    
modelname = 'iPE_pseudoconcat_noResp';
modeldir = fullfile(subdir,'fxUnivariate',modelname);
fprintf('%s\n\n',modeldir);

subdirfunc = fullfile(studydir,subno,'func');
sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));
    
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = {modeldir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

scans = [];
rctr = 0;
clear nvol
    
%% find all functional scans across sessions and runs   
    
for ss = 1:nSess 

    sess_name = sessdirsfunc(ss).name;
    
for r = 1:nruns
    rctr = rctr + 1;
    scandir = fullfile(subdir,'func',sess_name,sprintf('Run%d',r));
    filenames = dir(fullfile(scandir,'s6wf*.nii')); % using the 6mm images
    nvol(rctr) = length(filenames);
    for i = 1:nvol(rctr)
        scans = [scans; {fullfile(scandir,sprintf('%s,1',filenames(i).name))}];
    end    
end
end

    
% these two files are combined across runs
onsfile = fullfile(modeldir,'Onsets.mat');
regfile = fullfile(modeldir,'NuisanceRegressors.txt');

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

% Save model
jobname = fullfile(Jobdir, sprintf('SpecGLM_Univariate_%s_%s.mat',modelname,subno));    
save(jobname, 'matlabbatch');

% Run model
spm_jobman('run', matlabbatch)
clear matlabbatch

% To adjust session difference
% handle concatenate design with this SPM function
spm_fmri_concatenate(fullfile(modeldir,'SPM.mat'),nvol); 
    
% Estimate the model
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(modeldir,'SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
  
% Save model
jobname = fullfile(Jobdir, sprintf('EstGLM_Univariate_%s_%s.mat',modelname,subno));    
save(jobname, 'matlabbatch');

% Run model
spm_jobman('run', matlabbatch)
clear matlabbatch
    
    

