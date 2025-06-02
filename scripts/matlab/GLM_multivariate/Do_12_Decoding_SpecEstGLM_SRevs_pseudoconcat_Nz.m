
% The purpose of this script is to specify and estimate GLM
% of the concatenate models for decoding analysis.

% Each time rerunning this script will redo the GLM estimation and ask user
% to confirm.

%% setup
modelname = 'SRevs_pseudoconcat_Nz';
DecodingSubDir = fullfile(parentDir,'DecodingGLM',subno);

modeldir = fullfile(DecodingSubDir,'fxMultivariate',modelname);

fprintf('Specify and estimate GLM for multivariate %s on %s\n',modelname,subno);

subdirfunc = fullfile(studydir,subno,'func');
sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));

scans = [];
rctr = 0;
clear nvol

%%
for ss = 1:nSess 
    
    sess_name = sessdirsfunc(ss).name;
    
for r = 1:nruns
    rctr = rctr + 1;
    scandir = fullfile(subdir,'func',sess_name,sprintf('Run%d',r));
    filenames = dir(fullfile(scandir,'s2mmwf*.nii')); % using the 2mm images for decoding analysis (smaller smoothing size)
    nvol(rctr) = length(filenames);
    for i = 1:nvol(rctr)
        scans = [scans; {fullfile(scandir,sprintf('%s,1',filenames(i).name))}];
    end    
end
end

 % these two files are combined across runs and sessions
onsfile = fullfile(modeldir,'Onsets.mat');
regfile = fullfile(modeldir,'NuisanceRegressors.txt');

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

% Save model
jobname = fullfile(Jobdir, sprintf('SpecGLM_Multivariate_%s_%s.mat',modelname,subno));    
save(jobname, 'matlabbatch');

% Run model
spm_jobman('run', matlabbatch)
clear matlabbatch
 
% handle concatenate design with this SPM function
spm_fmri_concatenate(fullfile(modeldir,'SPM.mat'),nvol); 
 
% Estimate the model
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(modeldir,'SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

% Save model
jobname = fullfile(Jobdir, sprintf('EstGLM_Multivariate_%s_%s.mat',modelname,subno));    
save(jobname, 'matlabbatch');

% Run model
spm_jobman('run', matlabbatch)
clear matlabbatch




