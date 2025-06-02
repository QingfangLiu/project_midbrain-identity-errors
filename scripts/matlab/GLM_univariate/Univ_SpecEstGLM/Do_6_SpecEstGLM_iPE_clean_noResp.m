
% The purpose of this script is to specify and estimate GLM.

% Rerun this script for the same subject will ask to confirm and remove all
% previous GLM estimation results (including contrasts)

fprintf('Specifying and estimating GLM for %s\n',subno);

modelname = 'iPE_clean_noResp';
modeldir = fullfile(subdir,'fxUnivariate',modelname);
fprintf('%s\n\n',modeldir);

spm fmri % code has to have this line to use 'cfg_dep' function
clear matlabbatch
% Specify the model
matlabbatch{1}.spm.stats.fmri_spec.dir = {modeldir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

sessdirs = dir(fullfile(studydir,subno,'func','Day*'));

rctr = 0;
for ss = 1:nSess 
    sess_name = sessdirs(ss).name;
for r = 1:nruns
    rctr = rctr + 1;
    scandir = fullfile(subdir,'func',sess_name,sprintf('Run%d',r));
    onsfile = fullfile(modeldir,sprintf('%s_Run%d_Onsets.mat',sess_name,r));
    regfile = fullfile(subdir,'NuisanceRegressor',sprintf('nuisance_regressors_%s_%s_Run%d.txt',subno,sess_name,r)); % nuisance regressor

    tmp = [];
    filenames = dir(fullfile(scandir,'s6wf*.nii'));

    for i = 1:length(filenames)
        tmp = [tmp; {fullfile(scandir,sprintf('%s,1',filenames(i).name))}];
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess(rctr).scans = tmp;
    matlabbatch{1}.spm.stats.fmri_spec.sess(rctr).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(rctr).multi = {onsfile};
    matlabbatch{1}.spm.stats.fmri_spec.sess(rctr).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(rctr).multi_reg = {regfile};
    matlabbatch{1}.spm.stats.fmri_spec.sess(rctr).hpf = 128;

end
end

matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.1;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''}; %{fullfile(studydir,'ROIs','rmean_fullMB_point9.nii')};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

%Estimate the model
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% save the job
fname = fullfile(Jobdir, sprintf('SpecEstGLM_Univariate_%s_%s.mat',modelname,subno));
save(fname, 'matlabbatch');
spm_jobman('run', matlabbatch);   

clear matlabbatch
   



