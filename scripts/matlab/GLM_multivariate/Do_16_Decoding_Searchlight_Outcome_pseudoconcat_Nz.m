

% note: revisiting this analysis for my DS job interview as a
% classification problem (4/29/2025)
% folder structures have been reorganized since previous work so cleaned up
% the code to make it ready to run

% previous code is structured as a function to allow testing different ways
% of scaling the data, but now I removed this messy structure 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear
clc
close all

parentDir = '/Users/qingfangliu/Library/CloudStorage/Dropbox/KahntLab/Project_SPEEDTMS/Analysis';  % the entire Analysis folder
addpath(parentDir)
addpath('/Users/qingfangliu/libsvm-3.25/matlab') % add package for svm

External = '/Volumes/QF10TB/SPEEDTMS_results';  % external drive containing results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SubInfoFile = fullfile(parentDir,'SubInfo','SubjectConds.XLSX');
SubInfo = readtable(SubInfoFile) ;
Subs = SubInfo.Sub; % subject index
Excluded = SubInfo.Excluded; % whether excluding from analysis
Subs = Subs(Excluded==0); % excluding subs from analysis
nsubs = length(Subs);

load(fullfile(External,'Behavior','Alld.mat')); % load behavioral data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PATHS AND VARIABLES TO CHANGE
modelname = 'Outcome_pseudoconcat_Nz';
decodetype = 'LOOCV';

do_agg = 1; % aggregate across subjects
overwrite = 0; % whether to overwrite previous results

% smooth
sk = [6 6 6]; % smoothing kernel for resulting map
do_smooth = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nruns = 3; 
nSess = 2;
nallruns = nruns * nSess;


%% PATHS AND VARIABLES TO CHANGE

% use beta (default) or t-stat for decoding 
dat_use_strings = {'beta'};

% whether scale data by subtracting means
%dat_scale_strings = {'raw','Mean_across_conds','Mean_across_all_voxels','Mean_across_sl_voxels'};
dat_scale_strings = {'raw'};

maskmap = fullfile(parentDir,'Masks','Group_mask.nii');
maskname = 'WB';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GET DIMENSIONS
maskvol_hdr = spm_vol(maskmap);
maskvol_vol = spm_read_vols(maskvol_hdr);
sz = size(maskvol_vol);
lin_index = find(maskvol_vol); % the linear index of the mask in 79*95*79 space


%%

for n_dat_use = 1:length(dat_use_strings)
    dat_use = dat_use_strings{n_dat_use};
    
    for n_dat_scale = 1:length(dat_scale_strings)
        dat_scale = dat_scale_strings{n_dat_scale};
    
        for rget = 4 % searchlight radius to get in voxel
            RunSearchlight(maskmap,maskname,dat_use,dat_scale,rget,do_agg,overwrite)
        end
    
    end
end


%%
dat_use = 'beta';
dat_scale = 'raw';
rget = 4;


%% SEARCHLIGHT

% define reference voxel
ref_vox = round(sz/2);

% calculate distance to ref voxel for every voxel (radii)
[MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3));
radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2);

% prototype sphere index for radii<rget that can be used everywhere in the brain
% 251 voxels based on rget=4 (8mm sphere)
radius_index = find(radii<rget) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3));

% compute number of searchlights
nvox = length(lin_index);

% resultsvolume index
linindexconv = zeros(sz);
linindexconv(lin_index) = 1:length(lin_index); % consective integers at mask locations

    
%%
for i = 1:nsubs  % do decoding for each subject separately
    
subno = sprintf('Sub%d',Subs(i));
fprintf('%s\n',subno)
subjpath = fullfile(External,'Decoding_Searchlight',subno); % where to store the searchlight decoding results
    
% results path (where to save current searchlight results)
res_path = fullfile(subjpath, modelname, maskname, dat_use, dat_scale, sprintf('r%d',rget));
if ~exist(res_path,'dir')
    mkdir(res_path)
end

fprintf('%s - %s - %s\n',subno,dat_use,dat_scale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data
% prepare to organize betas into runs
data = cell(nallruns,1);
Labels = cell(nallruns,1);
nconds = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% organize data for each run

rctr = 0;
for ss = 1:nSess % sessions
for r = 1:nruns % runs   
    rctr = rctr + 1;
    fprintf('Data organization: %s-r%d\n',subno,rctr)
    temp = zeros(nconds,nvox); % dimension of each run's data: nconds * nvox
    for n = 1:nconds
        bidx = n + 3 * (rctr - 1);  % beta estimate index
        bname = fullfile(External,'DecodingGLM',subno,'fxMultivariate',modelname,sprintf('beta_%04d.nii',bidx));
        vol = spm_read_vols(spm_vol(bname));
        temp(n,:) = vol(lin_index); % data from ROI voxels
    end 
    data{rctr,1} = temp; % data organized for each run
    Labels{rctr,1} = (1:3)'; % labels for each run
    
end % run
end % session


%% subtract mean from each voxel

if strcmp(dat_scale,'Mean_across_all_voxels')
    data_avgd = cell(size(data));
    for r = 1:nruns
        cond_avg = mean(data{r,1},2,'omitnan'); % average of data across voxels for each cond
        cond_avg_mat = repmat(cond_avg,1,nvox); % organize averages into a matrix
        data_avgd{r} = data{r} - cond_avg_mat; % substract this average from the raw data
    end
    data = data_avgd; % use the scaled data for later decoding analysis
end


%%
% check if there is decoding result already
if ~any(size(dir(fullfile(res_path,'s*Acc*.nii')),1)) || (overwrite == 1)
  
    % define the labels for decoding for each run
    resultsvol_vol = zeros(sz); % used when saving the decoding accuracy (moved this to be clearer)
    temp_sl_res = zeros(nvox,1); % used for saving searchlight results
    
parfor cnt = 1:nvox
    
    if mod(cnt,1000)==0
        fprintf('Searchlight - V: %d/%d - %d\n',cnt, nvox,Subs(i));
    end
    
    % move the prototype sphere index to position of current voxel
    indexindex = radius_index + lin_index(cnt);
    
    % eliminate indices outside the brain mask
    Useindex = intersect(lin_index, indexindex);
    
    % searchlight indices in data space
    dat_index = linindexconv(Useindex);
     
    runvec = 1:nallruns; % LOOCV using each run as test run by loop
    rctr = 0;
    accvals = zeros(length(runvec),1); % pre-define before saving acc vals in each LOOCV
    
    for r = runvec
        rctr = rctr + 1;
        rtrain = runvec(runvec~=r); % leave one run cross-validation
        rtest = r;
        vectors_train = [];
        labels_train = [];
    
        for t = rtrain
            vectors_train = [vectors_train; data{t,1}(:,dat_index)];
            labels_train = [labels_train; Labels{t,1}];
        end

        vectors_test = data{rtest,1}(:,dat_index);
        labels_test = Labels{rtest,1};
                   
        if strcmp(dat_scale,'Mean_across_conds')
            voxel_avg_train = mean(vectors_train,1); % average of data across training conds for each voxel
            voxel_avg_mat_train = repmat(voxel_avg_train,half_length_train * 2, 1); % organize averages into a matrix
            vectors_train = vectors_train - voxel_avg_mat_train; % substract this average from the raw data
            
            voxel_avg_test = mean(vectors_test,1); % average of data across testing conds for each voxel
            voxel_avg_mat_test = repmat(voxel_avg_test,half_length_test * 2, 1); % organize averages into a matrix
            vectors_test = vectors_test - voxel_avg_mat_test; % substract this average from the raw data
        end
        
        if strcmp(dat_scale,'Mean_across_sl_voxels')
            cond_avg_train = mean(vectors_train,2,'omitnan'); % average of training data across voxels for each cond
            cond_avg_mat_train = repmat(cond_avg_train,1,length(dat_index)); 
            vectors_train = vectors_train - cond_avg_mat_train; 

            cond_avg_test = mean(vectors_test,2,'omitnan'); % average of test data across voxels for each cond
            cond_avg_mat_test = repmat(cond_avg_test,1,length(dat_index)); 
            vectors_test = vectors_test - cond_avg_mat_test; 
        end

        model = svmtrain(labels_train,vectors_train,'-s 0 -t 0 -q -c 0.0001'); % -s svm_type; -t kernel_type; -q quite mode
        [pred_label, accuracy, ~] = svmpredict(labels_test,vectors_test,model,'-q'); % pred_label: labels of testing data predicted from the trained model
        accvals(rctr) = accuracy(1);
        
    end
    
    outcome = mean(accvals); % average accuracies (one value saved in current voxel location)
    temp_sl_res(cnt) = outcome - 33.33;  % save sl results into 
    
end % cnt

resultsvol_vol(lin_index) = temp_sl_res;  % assign the acc results to 3d space

% results dimensions
resultsvol_hdr = spm_vol(bname);

% save result maps
Resultname = sprintf('Acc_%s_%s.nii', subno, decodetype);
resultsvol_hdr.fname = fullfile(res_path,Resultname); % save file path into hdr.fname
spm_write_vol(resultsvol_hdr, resultsvol_vol);

%% Smooth
if do_smooth
    clear matlabbatch
    matlabbatch{1}.spm.spatial.smooth.data = {fullfile(res_path,Resultname)};
    matlabbatch{1}.spm.spatial.smooth.fwhm = sk;
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = sprintf('s%01d', sk(1));
    
    spm_jobman('run', matlabbatch);
    clear matlabbatch
end

end % end of checking existence of decoding results

end % subject


%% aggregate decoding results across subjects

if do_agg

    % where to save aggregate searchlight results
    pathname = fullfile(parentDir,'Searchlight', 'Agg', modelname, maskname, dat_use, dat_scale, sprintf('r%d',rget));
    if ~exist(pathname, 'dir')
        mkdir(pathname);
    end
    
    % check if t-test has been run in this folder
    if ~any(size(dir(fullfile(pathname,'beta*.nii')),1)) || (overwrite == 1)
    
    subctr = 0;
    ffiles = cell(nsubs,1);  % read decoding accuracy images from all subjects
    for subj = Subs'
        subctr = subctr + 1;
        subno = sprintf('Sub%d',subj);
        subjpath = fullfile(parentDir,'Searchlight',subno);
        filename = sprintf('Acc_%s_%s.nii', subno, decodetype);
        filename = sprintf('s6%s,1',filename); % use smoothed images
        ffiles{subctr} = fullfile(subjpath, modelname, maskname, dat_use, dat_scale, sprintf('r%d',rget), filename); 
    end
            
    clear matlabbatch
    matlabbatch{1}.spm.stats.factorial_design.dir = {pathname}; 
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = ffiles;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {maskmap};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    % model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    % contrast
    matlabbatch{3}.spm.stats.con.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'DecodingAcc';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

    % save job batch
    jobname = fullfile(pathname, 'JobGroupTest.mat');
    save(jobname, 'matlabbatch');
        
    %run batch
    spm_jobman('run', matlabbatch);  
    clear matlabbatch
    
    end % end of check

end % if agg



