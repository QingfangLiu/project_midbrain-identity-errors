
fprintf('Specifying and estimating GLM for %s\n',subno);

modelname = 'iPE_clean';
modeldir = fullfile(subdir,'fxUnivariate',modelname);
fprintf('%s\n\n',modeldir);

clear matlabbatch
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

    Scan_4D = dir(fullfile(scandir,'s6*_4D.nii'));
    im_name = fullfile(scandir,Scan_4D.name);
    
    ims = spm_select('expand', im_name);  % get 3D image names from 4D
    ims = cellstr(ims); % convert char to cell

    matlabbatch{1}.spm.stats.fmri_spec.sess(rctr).scans = ims;
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
matlabbatch{1}.spm.stats.fmri_spec.mask = {''}; 
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(modeldir,'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

fname = fullfile(Jobdir, sprintf('SpecEstGLM_Univariate_%s_%s.mat',modelname,subno));
save(fname, 'matlabbatch');
spm_jobman('run', matlabbatch);   
clear matlabbatch
   
