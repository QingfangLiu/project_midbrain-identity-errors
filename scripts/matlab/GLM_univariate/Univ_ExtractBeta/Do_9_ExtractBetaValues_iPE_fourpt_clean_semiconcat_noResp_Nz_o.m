
%% extract values in functional midbrain masks (WB corrected)

%%
modelname = 'iPE_fourpt_clean_semiconcat_noResp_Nz';
modeldir = fullfile(ResDir,'UnivariateGLM',subno,'fxUnivariate',modelname); 
fprintf('Extracting %s beta values for %s\n',modelname,subno);

%% read ROI idx 
roiname = {'Insula_b.img','mPFC.img','Striatum_b.img','Conj_Thalamus.nii','AMG_b.nii'};

WBdir = fullfile(maskdir,'WB_corrected_masks');
Indepdir = fullfile(maskdir,'Indep');
roidir = {WBdir,WBdir,WBdir,WBdir,Indepdir};

clear roiidx

for i = 1:length(roiname)
    roi_vol = spm_read_vols(spm_vol(fullfile(roidir{i},roiname{i})));
    roiidx{i} = find(roi_vol);
end
   
%%
% find beta idx of interest  (1: reversal; 2: minus one; 3: plus one; 4: plus two)
bidx = reshape(1:8,4,2)';

%%
for sess = 1:2 % loop over two sessions
    for r = 1:length(roiidx) % loop over each ROI
        for b = 1:4 % loop each beta of interest
            bvol = spm_read_vols(spm_vol_nifti(fullfile(modeldir,sprintf('beta_%04d.nii',bidx(sess,b)))));
            beta_vals{r}(b,sess,:) = bvol(roiidx{r});
            MeanBOLD(b,r,sess) = nanmean(bvol(roiidx{r}));
        end      
    end
end

beta_file_name = sprintf('Betas_%s_others.mat',modelname); % updated with 'others' version
save(fullfile(ResDir,'Neural',subno,beta_file_name),'beta_vals','MeanBOLD')

clear beta_vals MeanBOLD roiidx