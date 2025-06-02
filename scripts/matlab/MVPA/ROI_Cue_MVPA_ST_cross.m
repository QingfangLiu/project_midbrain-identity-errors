
% This script uses single-trial beta estimates to do ROI-based MVPA
% This code runs under 'RunCode.m'

%%
maskdir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs';
subjpath = fullfile(ResDir,'MVPA_ROI_cross',subno); 
if ~exist(subjpath,'dir')
    mkdir(subjpath)
end

overwrite = 0;

% load beh data Alld.mat
load(fullfile(ResDir,'Behavior','Alld.mat'))

%% PATHS AND VARIABLES TO CHANGE

% model_list = {'STCues_pseudoconcat_Nz','STCues_LSS'};
model_list = {'STCues_LSS'};

n_model_list = length(model_list);
          
% mask labels and mask file names
mask_files = {'Func_OFC','FromSimiComp/Func_OFC.img';...
              'Func_LPFC','FromSimiComp/Func_LPFC.img';...
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

% reversal types: 1-2, 1-3, 2-1, 2-3, 3-1, 3-2
rev_types = [1,2;1,3;2,1;2,3;3,1;3,2];
n_rev_types = 6;

n_locs = 3;  % number of rev locs
n_cues = 2;  % number of cues

data = cell(nallruns,1);

for ss = 1:nSess 
    for r = 1:nruns    
        rctr = rctr + 1;
        fprintf('Data organization run %d\n',rctr)
        d = Alld{subj==Subs,rctr};
        ntrials = size(d,1);
        d = [d(:,1:3),(1:ntrials)'];
        
        id_array = zeros(n_locs,n_cues,n_rev_types);
        
        for cue = 1:2
            sep_d = d(d(:,1)==cue,:);
            revloc = find(sep_d(:,3)==1);
            useloc = sort([revloc-1;revloc;revloc+1]);  % find all trials pre-rev, rev, post-rev
            yy = sep_d(useloc,2:4);
            split_matrices = reshape(yy.', 3, 3, 6);
            result_matrix = zeros(5, 6);
            for i = 1:6
                result_matrix(:, i) = [split_matrices(1, 1, i); split_matrices(1, 2, i);...
                    split_matrices(3, 1, i); split_matrices(3, 2, i); split_matrices(3, 3, i)];
            end
            for type = 1:6
                which_col = find(result_matrix(1,:)==rev_types(type,1) & result_matrix(2,:)==rev_types(type,2));
                id_array(1,cue,type) = result_matrix(3,which_col);
                id_array(2,cue,type) = result_matrix(4,which_col);
                id_array(3,cue,type) = result_matrix(5,which_col);
            end
        end
        
         beta_temp = zeros(n_locs,n_cues,n_rev_types,nvox); 
        for loc = 1:n_locs
        for cue = 1:n_cues
        for type = 1:n_rev_types
            bidx = id_array(loc,cue,type);
            bname = fullfile(ResDir,'Decoding_GLM',subno,'fxMultivariate',modelname,sprintf('Run%dTrial%d',rctr,bidx),'beta_0001.nii');
            vol = spm_read_vols(spm_vol(bname));
            beta_temp(loc,cue,type,:) = vol(lin_index);
        end 
        end
        end
   
        data{rctr,1} = beta_temp;
    end % run
end % sess

%%
corr_res = zeros(nallruns,ncorr_types);                                     % neural similarity in each run
corr_res_rev = zeros(nallruns,ncorr_types,nrevs_per_run);                   % also track each rev

for r = 1:nallruns      
    run_data = data{r,1};
    
    % compute correlation across cues for each odor
    corr_pre = zeros(2,n_rev_types);
    corr_post = zeros(2,n_rev_types);
        
    for type = 1:n_rev_types
            % rev-1 (cue1) & rev (cue2)
            vec1 = squeeze(run_data(1,1,type,:));
            vec2 = squeeze(run_data(2,2,type,:));
            corr_pre(1,type) = atanh(corr(vec1,vec2));

            % rev (cue1) & rev-1 (cue2)
            vec1 = squeeze(run_data(2,1,type,:));
            vec2 = squeeze(run_data(1,2,type,:));
            corr_pre(2,type) = atanh(corr(vec1,vec2));

            % rev (cue1) & rev+1 (cue2)
            vec1 = squeeze(run_data(2,1,type,:));
            vec2 = squeeze(run_data(3,2,type,:));
            corr_post(1,type) = atanh(corr(vec1,vec2));

            % rev+1 (cue1) & rev (cue2)
            vec1 = squeeze(run_data(3,1,type,:));
            vec2 = squeeze(run_data(2,2,type,:));
            corr_post(2,type) = atanh(corr(vec1,vec2));
    end
        
    corr_res(r,1) = mean(corr_pre,'all');  
    corr_res(r,2) = mean(corr_post,'all');   
end

save(save_name,"corr_res")                                   % save both variables

else
    
    fprintf('Already done! \n')
    fprintf('%s \n',res_path)
    
end % end of checking exist

end % mask loop
end % model list

