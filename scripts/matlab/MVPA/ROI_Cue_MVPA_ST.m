
% This script uses single-trial beta estimates to do ROI-based MVPA
% This code runs under 'RunCode.m'

%%
subjpath = fullfile(ResDir,'MVPA_ROI',subno); 
if ~exist(subjpath,'dir')
    mkdir(subjpath)
end

overwrite = 0;

% load beh data Alld.mat
load(fullfile(ResDir,'Behavior','Alld.mat'))

%% PATHS AND VARIABLES TO CHANGE

% model_list = {'STCues_pseudoconcat_Nz','STCues_LSS'};
model_list = {'STCues_pseudoconcat_Nz','STCues_LSS'};

n_model_list = length(model_list);
          
% mask labels and mask file names
mask_files = {'OFC_cTBS','SeedRegion_OFC_bilateral.nii';...
              'lOFC_cTBS','SeedRegion_lOFC.nii';...
              'rOFC_cTBS','SeedRegion_rOFC.nii';...
              'LPFC_cTBS','TargetRegion_LPFC_bilateral.nii';...
              'lLPFC_cTBS','TargetRegion_lLPFC.nii';...
              'rLPFC_cTBS','TargetRegion_rLPFC.nii';...
              'OFC_func','OFC_func_1e-3.img';...
              'lOFC_med','aal_OFC/rOFC_med_L.img';...
              'rOFC_med','aal_OFC/rOFC_med_R.img';...
              'OFC_med','aal_OFC/rOFC_med_b.img';...
              'lOFC_SimiComp','FromSimiComp/Left_OFC_p5e-3.img';... 
              'rOFC_SimiComp','FromSimiComp/Right_OFC_p5e-3.img';...
              'OFC_SimiComp','FromSimiComp/Bilateral_OFC_p5e-3.img';...
              'lLPFC_SimiComp','FromSimiComp/Left_LPFC_p5e-3.img';... 
              'rLPFC_SimiComp','FromSimiComp/Right_LPFC_p5e-3.img';...
              'LPFC_SimiComp','FromSimiComp/Bilateral_LPFC_p5e-3.img';...
              'Func_OFC','FromSimiComp/Func_OFC.img';...
              'Func_LPFC','FromSimiComp/Func_LPFC.img';...
              'OFC_MVPSA','FromMVPSA_updated/Left_OFC_p5e-3.img';...   
              'LPFC_MVPSA','FromMVPSA_updated/Bilateral_LPFC_p5e-3.img';...
              'MB_b','WB_corrected_masks/Conj_p5_MB_b.nii';...
              'mPFC','WB_corrected_masks/mPFC.img';...
              'Striatum_b','WB_corrected_masks/Striatum_b.img';...
              'AMG_b','Indep/AMG_b.nii';...
              'Thalamus','WB_corrected_masks/Conj_Thalamus.nii';...
              'Func_DS','FromSimiComp/Func_DS.img';...
              'Func_VS','FromSimiComp/Func_VS.img';...
              'Func_MB','FromSimiComp/Func_MB.img';...
              'Func_mPFC','FromSimiComp/Func_mPFC.img';...
              'Func_HPC','FromSimiComp/Func_HPC.img';...
              'Func_Insula','FromSimiComp/Func_Insula.img';...
              'Func_Thalamus','FromSimiComp/Func_Thalamus.img';...
              'Func_Acc','FromSimiComp/Func_Acc.img';...
              };
                 
n_masks = length(mask_files);

nrevs_per_run = 12;   % # of revs per run
nlocs_per_rev = 4;    % # of locs per rev
UseNames = {'pre_rev','post_rev'};
ncorr_types = length(UseNames);      % # of correlation types

%%
for md = 1:n_model_list
    modelname = model_list{md}; 
    
for m = 1:n_masks
    maskmap = fullfile(maskdir,mask_files{m,2});
    maskname = mask_files{m,1};
    maskvol_vol = spm_read_vols(spm_vol(maskmap));
    lin_index = find(maskvol_vol);
    nvox = length(lin_index);
    sz = size(maskvol_vol);  % 79*95*79

    res_path = fullfile(subjpath, modelname, maskname);
    if ~exist(res_path,'dir')
        mkdir(res_path)
    end
    
    % how results are going to be saved
    save_name = fullfile(res_path,'corr_vars.mat');
    
% check if results already exist
if (~exist(save_name,'file'))
        
    fprintf('Neural Similarity calculation in: \n')
    fprintf('%s \n',res_path)

%% data organize (for all voxels)
rctr = 0;                           % count number of runs
trial_ctr = 0;                      % count number of trials
data = cell(nallruns,1);

for ss = 1:nSess 
    for r = 1:nruns    
        rctr = rctr + 1;
        fprintf('Data organization run %d\n',rctr)
        d = Alld{subj==Subs,rctr};
        ntrials = size(d,1);
        d = [d(:,1:3),(1:ntrials)'];
        id_mat = zeros(nrevs_per_run,nlocs_per_rev);                        % a matrix of trial id for this run: 12 reversals * 4 rev locs
        for cue = 1:2
            sep_d = d(d(:,1)==cue,:);
            revloc = find(sep_d(:,3)==1); 
            id_mat((1:6) + 6*(cue-1),2) = sep_d(revloc,4);
            id_mat((1:6) + 6*(cue-1),1) = sep_d(revloc-1,4);
            id_mat((1:6) + 6*(cue-1),3) = sep_d(revloc+1,4);
            id_mat((1:6) + 6*(cue-1),4) = sep_d(revloc+2,4);
        end
        beta_temp = zeros(nlocs_per_rev,nrevs_per_run,nvox); 
        for rev = 1:nrevs_per_run
            parfor c = 1:nlocs_per_rev
                switch modelname
                case 'STCues_pseudoconcat_Nz'
                    bidx = id_mat(rev,c) + trial_ctr; 
                    bname = fullfile(ResDir,'Decoding_GLM',subno,'fxMultivariate',modelname,sprintf('beta_%04d.nii', bidx));
                case 'STCues_LSS'
                    bidx = id_mat(rev,c);
                    bname = fullfile(ResDir,'Decoding_GLM',subno,'fxMultivariate',modelname,sprintf('Run%dTrial%d',rctr,bidx),'beta_0001.nii');
                end
                vol = spm_read_vols(spm_vol(bname));
                beta_temp(c,rev,:) = vol(lin_index); 
            end 
        end % rev
        trial_ctr = trial_ctr + ntrials;
        data{rctr,1} = beta_temp;
    end % run
end % sess

%%
corr_res = zeros(nallruns,ncorr_types);                                     % neural similarity in each run
corr_res_rev = zeros(nallruns,ncorr_types,nrevs_per_run);                   % also track each rev

for r = 1:nallruns      
    run_data = data{r,1};
    corrmat = zeros(nlocs_per_rev,nlocs_per_rev,nrevs_per_run);             % correlation matrix for each reversal
    for rev = 1:nrevs_per_run
        use_data = squeeze(run_data(:,rev,:));
        corrmat(:,:,rev) = corrcoef(use_data','Rows','complete');           % calculate pairwise correlation & ignore NAs 
        corrmat(:,:,rev) = atanh(corrmat(:,:,rev));                         % do Fisher's z-transformation on the corrmat before averaging
        corr_res_rev(r,1,rev) = corrmat(1,2,rev);
        corr_res_rev(r,2,rev) = corrmat(2,3,rev);
    end
    avg_corrmat = squeeze(mean(corrmat,3));                                 % average across revs
    corr_res(r,1) = avg_corrmat(1,2);
    corr_res(r,2) = avg_corrmat(2,3);
end

save(save_name,"corr_res","corr_res_rev")                                   % save both variables

else
    
    fprintf('Already done! \n')
    fprintf('%s \n',res_path)
    
end % end of checking exist

end % mask loop
end % model list

