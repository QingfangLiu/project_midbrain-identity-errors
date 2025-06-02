
% the purpose of this code is to extract the beta values for iPE index, for
% each reversal
% this has to be based on the single-trial models with onsets of outcome
% use the three functional ROIs in the iPE part of the analysis

%%
maskdir = fullfile(parentDir,'ROIs');
model_list = {'ST_Outcome_LSS','STOdorD_pseudoconcat_Nz'};
n_model_list = length(model_list);

mask_files = {'WB_corrected_masks/Conj_p5_MB_b.nii'...
              'WB_corrected_masks/LPFC_b.hdr'};
mask_labels = {'Midbrain','LPFC'};          
n_masks = length(mask_files);
assert(length(mask_labels) == n_masks) % ensure mask_files and mask_labels have equal # elements

nrevs_per_run = 12;   % # of revs per run
nlocs_per_rev = 4;    % # of locs per rev

% to save all beta values averaged across voxels in each mask
all_betas = zeros(n_model_list,n_masks,nSess,nruns,nrevs_per_run,nlocs_per_rev);

%%
for md = 1:n_model_list
    modelname = model_list{md}; 
    
for m = 1:n_masks
    maskmap = fullfile(maskdir,mask_files{m});
    maskname = mask_labels{m};
    maskvol_vol = spm_read_vols(spm_vol(maskmap));
    lin_index = find(maskvol_vol);
    
    %% data organize 
    rctr = 0;                           % count number of runs
    trial_ctr = 0;                      % count number of trials

    for ss = 1:nSess 
    for r = 1:nruns    
        rctr = rctr + 1;    
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
         
        for rev = 1:nrevs_per_run
            parfor c = 1:nlocs_per_rev
                switch modelname
                case 'STOdorD_pseudoconcat_Nz'
                    bidx = id_mat(rev,c) + trial_ctr; 
                    bname = fullfile(External,'UnivariateGLM',subno,'fxUnivariate',modelname,sprintf('beta_%04d.nii', bidx));
                case 'ST_Outcome_LSS'
                    bidx = id_mat(rev,c);
                    bname = fullfile(External,'UnivariateGLM',subno,'fxUnivariate',modelname,sprintf('Run%dTrial%d',rctr,bidx),'beta_0001.nii');
                end
                vol = spm_read_vols(spm_vol(bname));
                all_betas(md,m,ss,r,rev,c) = nanmean(vol(lin_index)); 
            end 
        end % rev
        trial_ctr = trial_ctr + ntrials;
    
    end % run
    end % sess
    
end
end

pathname = fullfile(parentDir,'NeuralAnalyzeRes',subno);
if ~exist(pathname,'dir')
    mkdir(pathname)
end

beta_file_name = 'betas_per_rev_from_ST.mat';
save(fullfile(pathname,beta_file_name),'all_betas')   



