
clc; clear;

ResDir = '/Volumes/QF10TB/SPEEDTMS_results'; 
maskpath = fullfile(ResDir,'GlobalConn','Masks');

SubInfoFile = fullfile(ResDir,'Behavior','SubjectConds.XLSX');
SubInfo = readtable(SubInfoFile) ;
Subs = SubInfo.Sub; % subject index
Conds = SubInfo.Cond; % TMS condition
Excluded = SubInfo.Excluded; % whether excluding from analysis

Subs = Subs(Excluded==0); % excluding subs from analysis
Conds = Conds(Excluded==0);
nsubs = length(Subs);

nSess = 2; % number of sessions
nruns = 3; % number of runs
nscans = 430; % number of scans in each run

conn_name = 'ConnectednessMap_3mm'; % change where the whole result is saved
Option = 'Unsigned'; % Unsigned or Signed or Squared

do_smooth = 1; % whether to do smooth
sk = 6; % smoothing kernel size (for individual contrast map before aggregation)
agg_contrast_map = 1; % whether to aggregate contrast maps

%% load masks
hdr = spm_vol(fullfile(maskpath,'gm_0.1_3mm.nii'));
gm = spm_read_vols(hdr);
gm_idx = find(gm > 0);
nvox = length(gm_idx); % number of voxels

%%
for subj = Subs'   
    subno = sprintf('Sub%d',subj);
    subdir = fullfile(studydir,subno);
    subdirfunc = fullfile(subdir,'func');
    sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));
    
    % prepare where to save the connResults
    respath = fullfile(ResDir, 'GlobalConn', conn_name, subno, Option);
    if ~exist(respath,'dir')
        mkdir(respath)
    end
    
    for ss = 1:nSess
        
        sess_name = sessdirsfunc(ss).name;
    
    for r = 1:nruns
        
        fprintf('Working on %s Sess %d Run %d\n', subno, ss, r);
        nhdr = spm_vol(fullfile(maskpath,'wm_0.9_3mm.nii'));
        nhdr.fname = fullfile(respath, sprintf('%s_%s_Run%d.nii', conn_name, sess_name, r));
        
        if ~exist(nhdr.fname,'file') % check if it has been created
            % load extracted and filtered scan data
            load(fullfile(ResDir, 'GlobalConn', subno, sess_name, sprintf('Run%d',r), 'tc_filtered_3mm.mat')); 
            dat = zscore(dat); % zcore for convenient Pearson correlation calculation
            cors = zeros(1,nvox); % initialize a vector for each voxel
        
            parfor i = 1:nvox
                tmp = (dat(:,i)' * dat)./ nscans; % Pearson correlation on z-scored values (with each voxel)

                if strcmp(Option,'Unsigned')
                    cors(i) = sum(abs(atanh(tmp(tmp<1))))/nvox; %UNSIGNED recommended (average across voxels, transformed first)
                elseif strcmp(Option,'Signed')
                    cors(i) = sum(atanh(tmp(tmp<1)))/nvox; % SIGNED
                else
                    cors(i) = sum((atanh(tmp(tmp<1))).^2)/nvox; % Squared
                end

                if ~mod(i,5000)
                    fprintf('\t\t%d of %d\n', i, nvox);
                end

            end

            vol = zeros(hdr.dim);
            vol(gm_idx) = cors;
            spm_write_vol(nhdr, vol); 
        
        end
        
    end % end of runs
    end % end of sessions
    
end % end of subjects



%% Compute difference map & do smoothing
% average across 3 runs

for subj = Subs'   
    subno = sprintf('Sub%d',subj);
    subdir = fullfile(ResDir, 'GlobalConn', subno);
    sessdirsfunc = dir(fullfile(subdir,'Day*'));
    
    % prepare where to save the connResults
    respath = fullfile(ResDir, 'GlobalConn', conn_name, subno, Option);
    image_dim = size(spm_read_vols(spm_vol(fullfile(maskpath,'wm_0.9_3mm.nii'))));
    map = zeros([image_dim,nSess,nruns]);
    
    for ss = 1:nSess
        sess_name = sessdirsfunc(ss).name;
        for r = 1:nruns
            map(:,:,:,ss,r) = spm_read_vols(spm_vol(fullfile(respath,sprintf('%s_%s_Run%d.nii', conn_name, sess_name, r))));
        end
    end
    
    % save difference map
    out1 = (map(:,:,:,1,1) + map(:,:,:,1,2) + map(:,:,:,1,3))./3;
    out2 = (map(:,:,:,2,1) + map(:,:,:,2,2) + map(:,:,:,2,3))./3;
    
    if strcmp(Conds{subj==Subs}, 'PA')
        out_P = out1;
        out_A = out2;
    else
        out_A = out1;
        out_P = out2;
    end
       
    % 1: out1 - out2
    hdr.fname = fullfile(respath, sprintf('%s_1_vs_2.nii', conn_name));
    spm_write_vol(hdr, out1 - out2);
    
    % 2: out_A - out_P
    hdr.fname = fullfile(respath, sprintf('%s_A_vs_P.nii', conn_name));
    spm_write_vol(hdr, out_A - out_P);

    % Smooth
    if do_smooth
        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {fullfile(respath,sprintf('%s_1_vs_2.nii', conn_name))
                                                  fullfile(respath, sprintf('%s_A_vs_P.nii', conn_name))};
        matlabbatch{1}.spm.spatial.smooth.fwhm = ones(1,3)*sk;
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's6';

        jobname = fullfile(respath, 'Job_smooth.mat');
        save(jobname, 'matlabbatch')

        spm_jobman('run', matlabbatch);
        clear matlabbatch
    end
end % subject



%% Compute difference map & do smoothing
% only first run of each

for subj = Subs'   
    
    subno = sprintf('Sub%d',subj);
    subdir = fullfile(ResDir, 'GlobalConn', subno);
    sessdirsfunc = dir(fullfile(subdir,'Day*'));
    
    % prepare where to save the connResults
    respath = fullfile(ResDir, 'GlobalConn', conn_name, subno, Option);
    image_dim = size(spm_read_vols(spm_vol(fullfile(maskpath,'wm_0.9_3mm.nii'))));
    map = zeros([image_dim,nSess,nruns]);
    
    for ss = 1:nSess
        sess_name = sessdirsfunc(ss).name;
        for r = 1:nruns
            map(:,:,:,ss,r) = spm_read_vols(spm_vol(fullfile(respath,sprintf('%s_%s_Run%d.nii', conn_name, sess_name, r))));
        end
    end
    
    % save difference map
    out1 = map(:,:,:,1,1);
    out2 = map(:,:,:,2,1);
    
    if strcmp(Conds{subj==Subs}, 'PA')
        out_P = out1;
        out_A = out2;
    else
        out_A = out1;
        out_P = out2;
    end
       
    % 1: out1 - out2
    hdr.fname = fullfile(respath, sprintf('%s_1_vs_2_run1.nii', conn_name));
    spm_write_vol(hdr, out1 - out2);
    
    % 2: out_A - out_P
    hdr.fname = fullfile(respath, sprintf('%s_A_vs_P_run1.nii', conn_name));
    spm_write_vol(hdr, out_A - out_P);

    % Smooth
    if do_smooth
        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {fullfile(respath,sprintf('%s_1_vs_2_run1.nii', conn_name))
                                                  fullfile(respath, sprintf('%s_A_vs_P_run1.nii', conn_name))};
        matlabbatch{1}.spm.spatial.smooth.fwhm = ones(1,3)*sk;
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's6';

        spm_jobman('run', matlabbatch);
        clear matlabbatch
    end
end % subject


%% aggregate contrast maps across subjects
    
contrast_names = {'1_vs_2','A_vs_P','1_vs_2_run1','A_vs_P_run1'};
n_contrast = length(contrast_names);

if agg_contrast_map
    
    for c = 1:n_contrast
    
    % where to save aggregate contrast map results
    pathname = fullfile(ResDir, 'GlobalConn', conn_name, 'Agg', Option, contrast_names{c});
    if ~exist(pathname, 'dir')
        mkdir(pathname);
    end

    % check if t-test has been run in this folder
    if ~any(size(dir(fullfile(pathname,'beta*.nii')),1))
    
    subctr = 0;
    ffiles = cell(nsubs,1);  % read decoding accuracy images from all subjects
    for subj = Subs'
        subctr = subctr + 1;
        subno = sprintf('Sub%d',subj);
        filename = sprintf('%s_%s.nii', conn_name, contrast_names{c});
        filename = sprintf('s6%s,1',filename); % use smoothed images
        ffiles{subctr} = fullfile(ResDir, 'GlobalConn', conn_name, subno, Option, filename);
    end
            
    clear matlabbatch
    matlabbatch{1}.spm.stats.factorial_design.dir = {pathname}; 
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = ffiles;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {fullfile(maskpath,'gm_0.1_3mm.nii')}; % use the gray matter mask
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % contrast
    matlabbatch{3}.spm.stats.con.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = sprintf('positive_%s', contrast_names{c});
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = sprintf('negative_%s', contrast_names{c});
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 1; % delete existing contrasts
    
    % save job
    jobname = fullfile(pathname, 'Job_agg.mat');
    save(jobname, 'matlabbatch')
    
    %run batch
    spm_jobman('run', matlabbatch);  
    clear matlabbatch
    
    end % end of check
    
    end % end loop of contrasts

end


%% aggregate contrast maps across subjects
% with considering covariates
% covariates include: 
% (1) subject-level TMS ratings (2 columns)
% (2) motion parameters
% only focus on A_vs_P_run1 contrast

% first, create a matrix R
TMSratingFile = fullfile(ResDir,'Behavior','TMSratings.XLSX');
TMSrating = readtable(TMSratingFile);
uncomfort_diff = TMSrating.uncomfort_cTBS - TMSrating.uncomfort_sham;
strong_diff = TMSrating.strong_cTBS - TMSrating.strong_sham;
R = [uncomfort_diff,strong_diff];
save(fullfile(ResDir,'Behavior','TMSratings_diff.mat'), 'R');
    
%%
contrast_names = {'A_vs_P_run1'};
n_contrast = length(contrast_names);
cov_file = fullfile(ResDir,'Behavior','TMSratings_diff.mat');

for c = 1:n_contrast

% where to save aggregate contrast map results
pathname = fullfile(ResDir, 'GlobalConn', conn_name, 'Agg_w_cov', Option, contrast_names{c});
if ~exist(pathname, 'dir')
    mkdir(pathname);
end

% check if t-test has been run in this folder
if ~any(size(dir(fullfile(pathname,'beta*.nii')),1))

subctr = 0;
ffiles = cell(nsubs,1);  % read decoding accuracy images from all subjects
for subj = Subs'
    subctr = subctr + 1;
    subno = sprintf('Sub%d',subj);
    filename = sprintf('%s_%s.nii', conn_name, contrast_names{c});
    filename = sprintf('s6%s,1',filename); % use smoothed images
    ffiles{subctr} = fullfile(ResDir, 'GlobalConn', conn_name, subno, Option, filename);
end

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {pathname}; 
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = ffiles;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov.files = {cov_file};
matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {fullfile(maskpath,'gm_0.1_3mm.nii')}; % use the gray matter mask
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% model estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% contrast
matlabbatch{3}.spm.stats.con.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = sprintf('positive_%s', contrast_names{c});
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = sprintf('negative_%s', contrast_names{c});
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1; % delete existing contrasts

% save job
jobname = fullfile(pathname, 'Job_agg_w_cov.mat');
save(jobname, 'matlabbatch')

%run batch
spm_jobman('run', matlabbatch);  
clear matlabbatch

end % end of check

end % end loop of contrasts

