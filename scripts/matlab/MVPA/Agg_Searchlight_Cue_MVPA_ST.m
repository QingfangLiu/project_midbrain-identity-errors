
% This code aggregates the MVPSA results from each subject into group
% level: both saved in external drive

clear; clc; close all

parentDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis';
WBmask = fullfile(parentDir,'ROIs','SPEEDTMS_external_mask.nii');           % use as the explicit mask

External = '/Volumes/QF10TB/SPEEDTMS_results'; 
MVPADir = fullfile(External,'MVPA_Searchlight');

UseSubDirs = dir(fullfile(MVPADir,'Sub*'));
subno_list = {UseSubDirs(:).name};
nsubs = length(subno_list);

model_list = {'STCues_pseudoconcat_Nz','STCues_LSS'};
n_model_list = length(model_list);

mask_files = {'SPEEDTMS_external_mask.nii'};
mask_labels = {'WB'};
n_masks = length(mask_files);

rgets = [2,3,4];


%%

for md = 1:n_model_list
    modelname = model_list{md}; 
    
for m = 1:n_masks
    maskname = mask_labels{m};

for rget = rgets  
    
    % do one-sample t-tests at the group level, either containing both TMS conds,
    % or sham only, or cTBS only
    test_names = {'pre-post','pre-post_sham','pre-post_cTBS'};
    n_test = length(test_names);

    for t = 1:n_test
    
    % where to save aggregate searchlight results
    pathname = fullfile(MVPADir, 'Agg', modelname, maskname, sprintf('r%d',rget), test_names{t});
    if ~exist(pathname, 'dir')
        mkdir(pathname);
    end

    % check if t-test has been run in this folder
    if ~any(size(dir(fullfile(pathname,'beta*.nii')),1))
        CorrDiffImages = cell(nsubs,1); 
    
    for j = 1:nsubs
        subno = subno_list{j};
        filename = sprintf('s6NeuralCorr_avg_%s.nii',test_names{t});                           
        CorrDiffImages{j} = fullfile(MVPADir, subno, modelname, maskname, sprintf('r%d',rget),filename); 
    end

    clear matlabbatch        
    matlabbatch{1}.spm.stats.factorial_design.dir = {pathname}; 
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = CorrDiffImages;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {WBmask};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    matlabbatch{3}.spm.stats.con.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'positive';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'negative';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

    jobname = fullfile(pathname, 'JobGroupTest.mat');
    save(jobname, 'matlabbatch');
    spm_jobman('run', matlabbatch);  
    clear matlabbatch

    end % check existence
    end % end of tests

%%
% to compare the neural diff between sham and cTBS
pathname = fullfile(MVPADir,'Agg', modelname, maskname, sprintf('r%d',rget), 'pre-post_sham_vs_cTBS');
if ~exist(pathname, 'dir')
    mkdir(pathname);
end

% check if t-test has been run in this folder
if  ~any(size(dir(fullfile(pathname,'beta*.nii')),1))
    ConImages_corr_diff_sham = cell(nsubs,1); 
    ConImages_corr_diff_cTBS = cell(nsubs,1); 
    for j = 1:nsubs
        subno = subno_list{j};
        filename_sham = 's6NeuralCorr_avg_pre-post_sham.nii';          
        filename_cTBS = 's6NeuralCorr_avg_pre-post_cTBS.nii';                              
        ConImages_corr_diff_sham{j} = fullfile(MVPADir, subno, modelname, maskname, sprintf('r%d',rget),filename_sham); 
        ConImages_corr_diff_cTBS{j} = fullfile(MVPADir, subno, modelname, maskname, sprintf('r%d',rget),filename_cTBS); 
    end

    clear matlabbatch    
    matlabbatch{1}.spm.stats.factorial_design.dir = {pathname};
    for j = 1:nsubs
        matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(j).scans = {ConImages_corr_diff_sham{j}; ConImages_corr_diff_cTBS{j}};
    end
    matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {WBmask};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    matlabbatch{3}.spm.stats.con.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'sham_vs_cTBS';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1,-1];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'cTBS_vs_sham';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1,1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

    jobname = fullfile(pathname, 'JobGroupTest.mat');
    save(jobname, 'matlabbatch');
    spm_jobman('run', matlabbatch);  
    clear matlabbatch

end % check existence

end
end
end





