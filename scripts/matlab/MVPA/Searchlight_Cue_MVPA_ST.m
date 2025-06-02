
% This script uses single-trial beta estimates to do searchlight RSA type
% of analysis. 
% This analysis is very important, making sure every step is correct!!!
% This code runs under 'RunCode.m'

%%
maskdir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs';
subjpath = fullfile(ResDir,'MVPA_Searchlight',subno);
if ~exist(subjpath,'dir')
    mkdir(subjpath)
end

% load beh data Alld.mat
load(fullfile(ResDir,'Behavior','Alld.mat'))

overwrite = 0;
% in post-analyze, comparing two TMS conds
if strcmp(Conds(subj==Subs),'PA')
    runsP = 1:3;
    runsA = 4:6;
else 
    runsP = 4:6;
    runsA = 1:3;
end

%% PATHS AND VARIABLES TO CHANGE

model_list = {'STCues_pseudoconcat_Nz','STCues_LSS'};
n_model_list = length(model_list);

mask_files = {'SPEEDTMS_external_mask.nii'};
mask_labels = {'WB'};
n_masks = length(mask_files);
assert(length(mask_labels) == n_masks) % ensure mask_files and mask_labels have equal # elements

rgets = [2,3,4];
nrevs_per_run = 12;   % # of revs per run
nlocs_per_rev = 4;    % # of locs per rev
UseNames = {'pre_rev','post_rev'};
ncorr_types = length(UseNames);      % # of correlation types

%%
for md = 1:n_model_list
    modelname = model_list{md}; 
    
for m = 1:n_masks
    maskmap = fullfile(maskdir,mask_files{m});
    maskname = mask_labels{m};
    maskvol_vol = spm_read_vols(spm_vol(maskmap));
    lin_index = find(maskvol_vol);
    nvox = length(lin_index);
    sz = size(maskvol_vol);                                                 % 79*95*79
    ref_vox = round(sz/2);                                                  % define reference voxel
    [MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3));                         % calculate distance to ref voxel for every voxel (radii)
    radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2);

for rget = rgets  
    res_path = fullfile(subjpath, modelname, maskname, sprintf('r%d',rget));
    res_path_per_rev = fullfile(res_path,'per_rev');
    jobpath = fullfile(res_path,'jobs');
    
    if ~exist(res_path_per_rev,'dir')
        mkdir(res_path_per_rev)
    end
    if ~exist(jobpath,'dir')
        mkdir(jobpath)
    end
    
if (~exist(fullfile(res_path,'s6NeuralCorr_run3_pre_rev.nii'),'file')) || (overwrite == 1)
    
    fprintf('Neural Similarity calculation in: \n')
    fprintf('%s \n',res_path)

%% prepare for searchlight
radius_index = find(radii<rget) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3)); % prototype sphere index for radii<rget that can be used everywhere in the brain
linindexconv = zeros(sz);                                               % resultsvolume index
linindexconv(lin_index) = 1:length(lin_index);                          % consective integers at mask locations
    
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
        id_mat = zeros(nrevs_per_run,nlocs_per_rev);                        % matrix of trial id for this run: 12 reversals * 4 rev locs
        for cue = 1:2
            sep_d = d(d(:,1)==cue,:);
            revloc = find(sep_d(:,3)==1); 
            id_mat((1:6) + 6*(cue-1),2) = sep_d(revloc,4);
            id_mat((1:6) + 6*(cue-1),1) = sep_d(revloc-1,4);
            id_mat((1:6) + 6*(cue-1),3) = sep_d(revloc+1,4);
            id_mat((1:6) + 6*(cue-1),4) = sep_d(revloc+2,4);
        end
        beta_temp = zeros(nlocs_per_rev,nrevs_per_run,nvox); 
        for rev = 1:nrevs_per_run                                           % loop through 12 reversals
            parfor c = 1:nlocs_per_rev                                      % 4 rev locs each reversal
                switch modelname
                case 'STCues_pseudoconcat_Nz'
                    bidx = id_mat(rev,c) + trial_ctr; 
                    bname = fullfile(External,'DecodingGLM',subno,'fxMultivariate',modelname,sprintf('beta_%04d.nii', bidx));
                case 'STCues_LSS'
                    bidx = id_mat(rev,c);
                    bname = fullfile(External,'DecodingGLM',subno,'fxMultivariate',modelname,sprintf('Run%dTrial%d',rctr,bidx),'beta_0001.nii');
                end
                vol = spm_read_vols(spm_vol(bname));
                beta_temp(c,rev,:) = vol(lin_index); 
            end 
        end % rev
        trial_ctr = trial_ctr + ntrials;
        data{rctr,1} = beta_temp;
        
        % also saving id_mat to summarize
        
        
    end % run
end % sess

%%
% save searchlight results for each correlation type and each run
all_sl_res = zeros(nallruns,ncorr_types,nvox);   
% save searchlight results also for each reversal
all_sl_res_rev = zeros(nallruns,ncorr_types,nrevs_per_run,nvox);

parfor cnt = 1:nvox % searchlight loop
    if mod(cnt,1e4)==0
        fprintf('Searchlight - V: %d/%d - %s\n',cnt, nvox, subno);
    end

    % prepare voxel index
    indexindex = radius_index + lin_index(cnt);                             % move the prototype sphere index to position of current voxel
    Useindex = intersect(lin_index, indexindex);                            % eliminate indices outside the brain mask
    dat_index = linindexconv(Useindex);                                     % index of current searchlight voxels
    
    temp_sl_res = zeros(nallruns,ncorr_types);
    temp_sl_res_revs = zeros(nallruns,ncorr_types,nrevs_per_run);           % keep track of each reversal
    
    for r = 1:nallruns      
        run_data = data{r,1};                                               % data from this run
        corrmat = zeros(nlocs_per_rev,nlocs_per_rev,nrevs_per_run);         % correlation matrix for each reversal
        for rev = 1:nrevs_per_run
            use_data = squeeze(run_data(:,rev,dat_index));
            corrmat(:,:,rev) = corrcoef(use_data','Rows','complete');       % calculate pairwise correlation & ignore NAs 
            corrmat(:,:,rev) = atanh(corrmat(:,:,rev));                     % do Fisher's z-transformation on the corrmat before averaging
            temp_sl_res_revs(r,1,rev) = corrmat(1,2,rev);
            temp_sl_res_revs(r,2,rev) = corrmat(2,3,rev);
        end
        avg_corrmat = squeeze(mean(corrmat,3));                             % average across revs
        temp_sl_res(r,1) = avg_corrmat(1,2);                                % pre-rev pattern similarity
        temp_sl_res(r,2) = avg_corrmat(2,3);                                % post-rev pattern similarity
    end
    
    all_sl_res(:,:,cnt) = temp_sl_res;                               
    all_sl_res_rev(:,:,:,cnt) = temp_sl_res_revs;
 end % searchlight loop

%% to save to disk
% get a tmp beta nii to use its header
switch modelname
case 'STCues_pseudoconcat_Nz'
    tmpbname = fullfile(External,'DecodingGLM',subno,'fxMultivariate',modelname,'beta_0001.nii');
case 'STCues_LSS'
    tmpbname = fullfile(External,'DecodingGLM',subno,'fxMultivariate',modelname,'Run1Trial1','beta_0001.nii');
end

% write the searchlight results into disk
for r = 1:nallruns
    for k = 1:ncorr_types
        resultsvol_vol = zeros(sz);   
        resultsvol_vol(lin_index) = all_sl_res(r,k,:); 
        resultsvol_hdr = spm_vol(tmpbname);
        
        Resultname = sprintf('NeuralCorr_run%d_%s.nii',r,UseNames{k});
        resultsvol_hdr.fname = fullfile(res_path,Resultname);                       
        spm_write_vol(resultsvol_hdr, resultsvol_vol);
    end
end

% save the nifti files for each run, each rev, and two simi measures
% into subfolders for clarity
for r = 1:nallruns
    for rev = 1:nrevs_per_run
    for k = 1:ncorr_types
        resultsvol_vol = zeros(sz);   
        resultsvol_vol(lin_index) = all_sl_res_rev(r,k,rev,:); 
        resultsvol_hdr = spm_vol(tmpbname);
        
        Resultname = sprintf('NeuralCorr_run%d_rev%d_%s.nii',r,rev,UseNames{k});
        resultsvol_hdr.fname = fullfile(res_path_per_rev,Resultname);                       
        spm_write_vol(resultsvol_hdr, resultsvol_vol);
    end
    end
end

%% post-analyze
    % the above analyses calculated and saved two correlations for each run
    % now it calculates two averages across all six runs
    % and subtracts the two averaged correlation maps
    % also separates based on TMS conditions
    % all of these are done using unsmoothed correlation images

    fprintf('Working on post-calculation: \n')
       
    %% average pre-rev correlation
    Allfile = dir(fullfile(res_path,'NeuralCorr_run*_pre_rev.nii'));
    file_paths = arrayfun(@(x) fullfile(x.folder,x.name), Allfile, 'UniformOutput', false);

    clear matlabbatch  
    matlabbatch{1}.spm.util.imcalc.input = file_paths;
    matlabbatch{1}.spm.util.imcalc.output = 'NeuralCorr_avg_pre_rev';
    matlabbatch{1}.spm.util.imcalc.outdir = {res_path};
    matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5+i6)/6';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    % average pre-rev correlation (P only)
    matlabbatch{2}.spm.util.imcalc.input = file_paths(runsP);
    matlabbatch{2}.spm.util.imcalc.output = 'NeuralCorr_avg_pre_rev_sham';
    matlabbatch{2}.spm.util.imcalc.outdir = {res_path};
    matlabbatch{2}.spm.util.imcalc.expression = '(i1+i2+i3)/3';
    matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{2}.spm.util.imcalc.options.mask = 0;
    matlabbatch{2}.spm.util.imcalc.options.interp = 1;
    matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
    
    % average pre-rev correlation (A only)
    matlabbatch{3}.spm.util.imcalc.input = file_paths(runsA);
    matlabbatch{3}.spm.util.imcalc.output = 'NeuralCorr_avg_pre_rev_cTBS';
    matlabbatch{3}.spm.util.imcalc.outdir = {res_path};
    matlabbatch{3}.spm.util.imcalc.expression = '(i1+i2+i3)/3';
    matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{3}.spm.util.imcalc.options.mask = 0;
    matlabbatch{3}.spm.util.imcalc.options.interp = 1;
    matlabbatch{3}.spm.util.imcalc.options.dtype = 4;
    
    %% average post-rev correlation
    Allfile = dir(fullfile(res_path,'NeuralCorr_run*_post_rev.nii'));
    file_paths = arrayfun(@(x) fullfile(x.folder,x.name), Allfile, 'UniformOutput', false);

    matlabbatch{4}.spm.util.imcalc.input = file_paths;
    matlabbatch{4}.spm.util.imcalc.output = 'NeuralCorr_avg_post_rev';
    matlabbatch{4}.spm.util.imcalc.outdir = {res_path};
    matlabbatch{4}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5+i6)/6';
    matlabbatch{4}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{4}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{4}.spm.util.imcalc.options.mask = 0;
    matlabbatch{4}.spm.util.imcalc.options.interp = 1;
    matlabbatch{4}.spm.util.imcalc.options.dtype = 4;
    
    % average post-rev correlation (P only)
    matlabbatch{5}.spm.util.imcalc.input = file_paths(runsP);
    matlabbatch{5}.spm.util.imcalc.output = 'NeuralCorr_avg_post_rev_sham';
    matlabbatch{5}.spm.util.imcalc.outdir = {res_path};
    matlabbatch{5}.spm.util.imcalc.expression = '(i1+i2+i3)/3';
    matlabbatch{5}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{5}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{5}.spm.util.imcalc.options.mask = 0;
    matlabbatch{5}.spm.util.imcalc.options.interp = 1;
    matlabbatch{5}.spm.util.imcalc.options.dtype = 4;
    
    % average post-rev correlation (A only)
    matlabbatch{6}.spm.util.imcalc.input = file_paths(runsA);
    matlabbatch{6}.spm.util.imcalc.output = 'NeuralCorr_avg_post_rev_cTBS';
    matlabbatch{6}.spm.util.imcalc.outdir = {res_path};
    matlabbatch{6}.spm.util.imcalc.expression = '(i1+i2+i3)/3';
    matlabbatch{6}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{6}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{6}.spm.util.imcalc.options.mask = 0;
    matlabbatch{6}.spm.util.imcalc.options.interp = 1;
    matlabbatch{6}.spm.util.imcalc.options.dtype = 4;
    
    jobname = fullfile(jobpath, 'Job_average_corrs.mat');
    save(jobname, 'matlabbatch');
    spm_jobman('run', matlabbatch);
    clear matlabbatch
    
    %% get corr (pre - post), using the averaged image
    matlabbatch{1}.spm.util.imcalc.input = {fullfile(res_path,'NeuralCorr_avg_pre_rev.nii');...
                                            fullfile(res_path,'NeuralCorr_avg_post_rev.nii')};
    matlabbatch{1}.spm.util.imcalc.output = 'NeuralCorr_avg_pre-post';
    matlabbatch{1}.spm.util.imcalc.outdir = {res_path};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1 - i2';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    % corr (pre - post), P only
    matlabbatch{2}.spm.util.imcalc.input = {fullfile(res_path,'NeuralCorr_avg_pre_rev_sham.nii');...
                                            fullfile(res_path,'NeuralCorr_avg_post_rev_sham.nii')};
    matlabbatch{2}.spm.util.imcalc.output = 'NeuralCorr_avg_pre-post_sham';
    matlabbatch{2}.spm.util.imcalc.outdir = {res_path};
    matlabbatch{2}.spm.util.imcalc.expression = 'i1 - i2';
    matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{2}.spm.util.imcalc.options.mask = 0;
    matlabbatch{2}.spm.util.imcalc.options.interp = 1;
    matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
    
    % corr (pre - post), A only
    matlabbatch{3}.spm.util.imcalc.input = {fullfile(res_path,'NeuralCorr_avg_pre_rev_cTBS.nii');...
                                            fullfile(res_path,'NeuralCorr_avg_post_rev_cTBS.nii')};
    matlabbatch{3}.spm.util.imcalc.output = 'NeuralCorr_avg_pre-post_cTBS';
    matlabbatch{3}.spm.util.imcalc.outdir = {res_path};
    matlabbatch{3}.spm.util.imcalc.expression = 'i1 - i2';
    matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{3}.spm.util.imcalc.options.mask = 0;
    matlabbatch{3}.spm.util.imcalc.options.interp = 1;
    matlabbatch{3}.spm.util.imcalc.options.dtype = 4;
    
    jobname = fullfile(jobpath, 'Job_get_corr_diff.mat');
    save(jobname, 'matlabbatch');
    spm_jobman('run', matlabbatch);
    clear matlabbatch
    
    %% smooth all the images in the folder
    Allfile = dir(fullfile(res_path,'NeuralCorr*.nii'));
    file_paths = arrayfun(@(x) fullfile(x.folder,x.name), Allfile, 'UniformOutput', false);
    
    matlabbatch{1}.spm.spatial.smooth.data = file_paths;
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6,6,6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's6';
    
    Allfile = dir(fullfile(res_path_per_rev,'NeuralCorr*.nii'));
    file_paths = arrayfun(@(x) fullfile(x.folder,x.name), Allfile, 'UniformOutput', false);
    
    matlabbatch{2}.spm.spatial.smooth.data = file_paths;
    matlabbatch{2}.spm.spatial.smooth.fwhm = [6,6,6];
    matlabbatch{2}.spm.spatial.smooth.dtype = 0;
    matlabbatch{2}.spm.spatial.smooth.im = 0;
    matlabbatch{2}.spm.spatial.smooth.prefix = 's6';

    jobname = fullfile(jobpath, 'Job_smooth_images.mat');
    save(jobname, 'matlabbatch');
    spm_jobman('run', matlabbatch);
    clear matlabbatch
    
else
    fprintf('Already done! \n')
    fprintf('%s \n',res_path)
    
end % check existence
    
end % r_gets
end % mask loop
end % model list

