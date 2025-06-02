

% This script does global conn calculation
% (1) focused on ROIs resliced on 3mm space
% (2) split each run into 2 time bins

clc; clear;

%% folder path
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'; 
maskpath = fullfile(ResDir,'GlobalConn','Masks'); % where the 3mm space masks are saved
ROIpath = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs/GlobalConn';

%%
SubInfoFile = fullfile(ResDir,'Behavior','SubjectConds.XLSX');
SubInfo = readtable(SubInfoFile);

Subs = SubInfo.Sub; % subject index
Conds = SubInfo.Cond; % TMS condition
Excluded = SubInfo.Excluded; % whether excluding from analysis
Subs = Subs(Excluded==0); % excluding subs from analysis
Conds = Conds(Excluded==0);
nsubs = length(Subs);

nSess = 2; % number of sessions
nruns = 3; % number of runs
nscans = 430 * 3; % number of scans in each run

conn_name = 'ConnectednessMap_3mm_2bins'; % adding '_2bins'
Option = 'Unsigned'; % Unsigned or Signed or Squared

%% load masks

% I still need this gray matter mask, although when looping voxels, I only
% need to calculate those voxels falling in the ROIs

hdr = spm_vol(fullfile(maskpath,'gm_0.1_3mm.nii'));
gm = spm_read_vols(hdr);
gm_idx = find(gm > 0);
nvox = length(gm_idx); % number of voxels

% ROI labels and ROI file names
ROIfiles = {'OFC','r3mmSeedRegion_OFC_bilateral.nii';...
            'LPFC','r3mmTargetRegion_LPFC_bilateral.nii';...
            'func_OFC','func_OFC.img';...
            'func_LPFC','func_LPFC.img';...
            'Insula_l','r3mmInsula_l.img';...
            'Insula_r','r3mmInsula_r.img';...
            'Insula_b','r3mmInsula_b.img';...
            'MB_l','r3mmConj_p5_MB_l.nii';...
            'MB_r','r3mmConj_p5_MB_r.nii';...
            'MB_b','r3mmConj_p5_MB_b.nii';...
            'mPFC','r3mmmPFC.img';...
            'Striatum_l','r3mmStriatum_l.img';...
            'Striatum_r','r3mmStriatum_r.img';...
            'Striatum_b','r3mmStriatum_b.img';...
            'AMG_l','r3mmAMG_l.nii';...
            'AMG_r','r3mmAMG_r.nii';...
            'AMG_b','r3mmAMG_b.nii';...
            'Thalamus','r3mmConj_Thalamus.nii';...
            };
nROIs = length(ROIfiles);

% split the 430 volumes into half
bins = [215, 215];
windows = [1,215;216,430];

%%
tic

for R = 1:nROIs
    
    ROI_vols = spm_read_vols(spm_vol(fullfile(ROIpath,ROIfiles{R,2})));
    ROI_idx = find(ROI_vols > 0);
    [LIA,LOCB] = ismember(ROI_idx,gm_idx);
    % if a voxel doesn't exist, then ignore it
    LOCB = LOCB(LIA); 
    n_ROI_vox = sum(LIA);

%%
for subj = Subs'   
    subno = sprintf('Sub%d',subj);
    sessdirsfunc = dir(fullfile(ResDir,'GlobalConn',subno,'Day*'));
    
    % prepare where to save the connResults
    respath = fullfile(ResDir, 'GlobalConn_ROI', conn_name, subno, Option);
    if ~exist(respath,'dir')
        mkdir(respath)
    end
    
for ss = 1:nSess
    sess_name = sessdirsfunc(ss).name;
for r = 1:nruns       
    save_name = sprintf('%s_cors_Sess%d_Run%d.mat',ROIfiles{R,1},ss,r);
    
if ~exist(fullfile(respath, save_name),'file')
    fprintf('Working on %s %s\n', save_name, subno);
    % load extracted and filtered scan data
    load(fullfile(ResDir, 'GlobalConn', subno, sess_name, sprintf('Run%d',r), 'tc_filtered_3mm.mat')); 
    
    cors = zeros(n_ROI_vox,length(bins)); 
    for t = 1:length(bins) % go through time windows

        use_dat = dat(windows(t,1):windows(t,2),:); 
        bin_size = windows(t,2) - windows(t,1) + 1;
        use_dat = zscore(use_dat);

        parfor i = 1:n_ROI_vox % here the change is simple, now only loop through voxels inside the ROIs
            idx = LOCB(i);   % find the voxel index of the OFC voxels in the gray matter voxel space
            tmp = (use_dat(:,idx)' * use_dat)./ bin_size; % Pearson correlation on z-scored values (with each voxel)

            if strcmp(Option,'Unsigned')
                cors(i,t) = sum(abs(atanh(tmp(tmp<1))))/nvox; %UNSIGNED
            elseif strcmp(Option,'Signed')
                cors(i,t) = sum(atanh(tmp(tmp<1)))/nvox; % SIGNED
            else
                cors(i,t) = sum((atanh(tmp(tmp<1))).^2)/nvox; % Squared
            end
        end % end of voxel loop
    end % end of bins
    avg_cors = mean(cors,1); % average across voxels
    % save corr matrix for each sess and each subject
    save(fullfile(respath, save_name), 'cors','avg_cors')
else
    fprintf('%s for %s already exists! \n',save_name,subno)
end % check exist

end % end of runs
end % end of sessions
end % end of subjects
end % end of ROIs

toc

%% organize GC results into an array
% dimension: nROIs * nsubs * nSess * ntimes

dat = zeros(nROIs,nsubs,nSess,nruns*2);
subctr = 0;
for subj = Subs'  
    subctr = subctr + 1;
    subno = sprintf('Sub%d',subj);
    for roi = 1:nROIs
    for s = 1:nSess
        for r = 1:nruns
            load(fullfile(ResDir, 'GlobalConn_ROI', conn_name, subno, Option,sprintf('%s_cors_Sess%d_Run%d.mat',ROIfiles{roi,1},s,r)))        
            dat(roi,subctr,s,2*r-1) = avg_cors(1);
            dat(roi,subctr,s,2*r) = avg_cors(2);
        end
    end
    end
end

ROIlabels = ROIfiles(:,1);
save(fullfile(ResDir, 'GlobalConn_ROI','GC_2bins_dat.mat'), 'dat', 'ROIlabels'); 


