
% updated to do MVPA using cross cues representation

%%
maskdir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs';
subjpath = fullfile(ResDir,'MVPA_Searchlight_cross',subno);
if ~exist(subjpath,'dir')
    mkdir(subjpath)
end

overwrite = 0;

% load beh data Alld.mat
load(fullfile(ResDir,'Behavior','Alld.mat'))

% in post-analyze, comparing two TMS conds
if strcmp(Conds(subj==Subs),'PA')
    runsP = 1:3;
    runsA = 4:6;
else 
    runsP = 4:6;
    runsA = 1:3;
end

%% PATHS AND VARIABLES TO CHANGE

% reduced to this version, using just one model, one mask (WB), 
% one rget (r=2)
% but keep the same data organizing structure for consistency

modelname = 'STCues_LSS';
maskmap = fullfile(maskdir,'SPEEDTMS_external_mask.nii');
maskname = 'WB';
rget = 2;

UseNames = {'pre_rev','post_rev'};
ncorr_types = length(UseNames);      % # of correlation types

%%
maskvol_vol = spm_read_vols(spm_vol(maskmap));
lin_index = find(maskvol_vol);
nvox = length(lin_index);
sz = size(maskvol_vol);                                                 % 79*95*79
ref_vox = round(sz/2);                                                  % define reference voxel
[MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3));                         % calculate distance to ref voxel for every voxel (radii)
radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2);

res_path = fullfile(subjpath, modelname, maskname, sprintf('r%d',rget));
jobpath = fullfile(res_path,'jobs');
if ~exist(jobpath,'dir')
    mkdir(jobpath)
end

%%
    
if (~exist(fullfile(res_path,'s6NeuralCorr_run3_pre_rev.nii'),'file')) || (overwrite == 1)
    
    fprintf('Neural Similarity calculation in: \n')
    fprintf('%s \n',res_path)

%% prepare for searchlight
radius_index = find(radii<rget) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3)); % prototype sphere index for radii<rget that can be used everywhere in the brain
linindexconv = zeros(sz);                                               % resultsvolume index
linindexconv(lin_index) = 1:length(lin_index);                          % consective integers at mask locations
    
%% data organize (for all voxels)
rctr = 0;                           % count number of runs

% reversal types: 1-2, 1-3, 2-1, 2-3, 3-1, 3-2
rev_types = [1,2;1,3;2,1;2,3;3,1;3,2];
n_rev_types = 6;

n_locs = 3;  % number of rev locs
n_cues = 2;  % number of cues

data = cell(nallruns,1);
% data organizing updated: now also considering two cues and three odors

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
% save searchlight results for each correlation type and each run
all_sl_res = zeros(nallruns,ncorr_types,nvox);   

parfor cnt = 1:nvox % searchlight loop
    if mod(cnt,1e4)==0
        fprintf('Searchlight - V: %d/%d - %s\n',cnt, nvox, subno);
    end

    % prepare voxel index
    indexindex = radius_index + lin_index(cnt);                             % move the prototype sphere index to position of current voxel
    Useindex = intersect(lin_index, indexindex);                            % eliminate indices outside the brain mask
    dat_index = linindexconv(Useindex);                                     % index of current searchlight voxels
    temp_sl_res = zeros(nallruns,ncorr_types);
    
    for r = 1:nallruns      
        run_data = data{r,1};                                               % data from this run
        
        % compute correlation across cues for each odor
        corr_pre = zeros(2,n_rev_types);
        corr_post = zeros(2,n_rev_types);
        
        for type = 1:n_rev_types
            % rev-1 (cue1) & rev (cue2)
            vec1 = squeeze(run_data(1,1,type,dat_index));
            vec2 = squeeze(run_data(2,2,type,dat_index));
            corr_pre(1,type) = atanh(corr(vec1,vec2));

            % rev (cue1) & rev-1 (cue2)
            vec1 = squeeze(run_data(2,1,type,dat_index));
            vec2 = squeeze(run_data(1,2,type,dat_index));
            corr_pre(2,type) = atanh(corr(vec1,vec2));

            % rev (cue1) & rev+1 (cue2)
            vec1 = squeeze(run_data(2,1,type,dat_index));
            vec2 = squeeze(run_data(3,2,type,dat_index));
            corr_post(1,type) = atanh(corr(vec1,vec2));

            % rev+1 (cue1) & rev (cue2)
            vec1 = squeeze(run_data(3,1,type,dat_index));
            vec2 = squeeze(run_data(2,2,type,dat_index));
            corr_post(2,type) = atanh(corr(vec1,vec2));
        end
        
        temp_sl_res(r,1) = mean(corr_pre,'all');                                  % pre-rev pattern similarity
        temp_sl_res(r,2) = mean(corr_post,'all');                                 % post-rev pattern similarity
    end
    all_sl_res(:,:,cnt) = temp_sl_res;                               
 end % searchlight loop

%% to save to disk
% get a tmp beta nii to use its header
tmpbname = fullfile(ResDir,'Decoding_GLM',subno,'fxMultivariate',modelname,'Run1Trial1','beta_0001.nii');

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
    
    jobname = fullfile(jobpath, 'Job_smooth_images.mat');
    save(jobname, 'matlabbatch');
    spm_jobman('run', matlabbatch);
    clear matlabbatch
    
else
    fprintf('Already done! \n')
    fprintf('%s \n',res_path)
    
end % check existence
    
