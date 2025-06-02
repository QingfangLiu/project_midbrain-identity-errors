

%%
modelname = 'iPE_fourpt_clean_pseudoconcat_noResp';
modeldir = fullfile(External,'UnivariateGLM',subno,'fxUnivariate',modelname); 
maskdir = fullfile(parentDir,'ROIs');

%% read ROI idx 
roi_files = {'WB_corrected_masks/Conj_p5_MB_b.nii'...
    'WB_corrected_masks/LPFC_b.img'};
roi_labels = {'MB','LPFC'};
nrois = length(roi_files);
roiidx = cell(nrois,1);

for i = 1:nrois
    roi_vol = spm_read_vols(spm_vol(fullfile(maskdir,roi_files{i})));
    roiidx{i} = find(roi_vol);
end
   
% find beta idx of interest  (1: reversal; 2: minus one; 3: plus one; 4: plus two)
bidx = reshape(1:24,4,6)';

%%
pathname = fullfile(parentDir,'NeuralAnalyzeRes',subno);
if ~exist(pathname,'dir')
    mkdir(pathname)
end

for roi = 1:nrois
for sess = 1:2 
    for run = 1:3
        tmp_run_idx = run + 3*(sess-1);
        for b = 1:4 % loop each beta of interest
            bvol = spm_read_vols(spm_vol_nifti(fullfile(modeldir,sprintf('beta_%04d.nii',bidx(tmp_run_idx,b)))));
            beta_vals(roi,sess,run,b) = nanmean(bvol(roiidx{roi}));
        end      
    end
end
end

beta_file_name = 'betas.mat';
save(fullfile(pathname,beta_file_name),'beta_vals','roi_labels')     