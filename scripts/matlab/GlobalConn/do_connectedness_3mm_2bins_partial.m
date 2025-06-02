
% This version was modified
% (1) focus on 2 ROIs to calculate global connectedness
% (2) split each run into 2 bins to get one more data point


clc; clear;

External = '/Volumes/QF10TB/SPEEDTMS_results';
studydir = '/Users/qingfangliu/Experiment';
ResDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ExptRes';
parentDir = fileparts(ResDir); % the entire Analysis folder

SubInfoFile = fullfile(parentDir,'SubjectConds.XLSX');
SubInfo = readtable(SubInfoFile) ;
Subs = SubInfo.Sub; % subject index
Conds = SubInfo.Cond; % TMS condition
Excluded = SubInfo.Excluded; % whether excluding from analysis

Subs = Subs(Excluded==0); % excluding subs from analysis
Conds = Conds(Excluded==0);
nsubs = length(Subs);

nSess = 2; % number of sessions
nruns = 3; % number of runs
nscans = 430 * 3; % number of scans in each run

conn_name = 'ConnectednessMap_3mm_2bins_partial'; % adding '_2bins_partial'
Option = 'Unsigned'; % Unsigned or Signed or Squared

%% load masks

% I still need this gray matter mask, although when looping voxels, I only
% need to calculate those voxels falling in the ROIs

maskpath = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/Global_connectedness';
hdr = spm_vol(fullfile(maskpath,'gm_0.1_3mm.nii'));
gm = spm_read_vols(hdr);
gm_idx = find(gm > 0);
nvox = length(gm_idx); % number of voxels

ROIfiles = {'r3mmSeedRegion_OFC_bilateral.nii','r3mmTargetRegion_LPFC_bilateral.nii'};
ROIlabels = {'OFC','LPFC'};

% split the 430 volumes into half
bins = [215, 215];
windows = [1,215;216,430];

%%
tic

for R = 1:length(ROIfiles)
    
    ROI_vols = spm_read_vols(spm_vol(fullfile(maskpath,'ROIs',ROIfiles{R})));
    ROI_idx = find(ROI_vols > 0);
    [LIA,LOCB] = ismember(ROI_idx,gm_idx);
    % if a voxel doesn't exist, then ignore it
    LOCB = LOCB(LIA); 
    n_ROI_vox = sum(LIA);
    
    % also get voxel idx of another ROI (use '_1' to differentiate)
    ROI_vols_1 = spm_read_vols(spm_vol(fullfile(maskpath,'ROIs',ROIfiles{3-R})));
    ROI_idx_1 = find(ROI_vols_1 > 0);
    [LIA_1,LOCB_1] = ismember(ROI_idx_1,gm_idx);
    LOCB_1 = LOCB_1(LIA_1); 
    n_ROI_vox_1 = sum(LIA_1);

%%
for subj = Subs'   
    subno = sprintf('Sub%d',subj);
    subdir = fullfile(studydir,subno);
    subdirfunc = fullfile(subdir,'func');
    sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));
    
    % prepare where to save the connResults
    respath = fullfile(External, 'GlobalConnResults', conn_name, subno, Option);
    if ~exist(respath,'dir')
        mkdir(respath)
    end
    
for ss = 1:nSess
    sess_name = sessdirsfunc(ss).name;
for r = 1:nruns       
    save_name = sprintf('%s_cors_Sess%d_Run%d.mat',ROIlabels{R},ss,r);
    
if ~exist(fullfile(respath, save_name),'file')
    fprintf('Working on %s %s\n', save_name, subno);
    % load extracted and filtered scan data
    load(fullfile(External, 'GlobalConnResults', subno, sess_name, sprintf('Run%d',r), 'tc_filtered_3mm.mat')); 
    
    cors = zeros(n_ROI_vox,length(bins)); 
    for t = 1:length(bins) % go through time windows

        use_dat = dat(windows(t,1):windows(t,2),:); 
        bin_size = windows(t,2) - windows(t,1) + 1;
        use_dat = zscore(use_dat);

        parfor i = 1:n_ROI_vox % here the change is simple, now only loop through voxels inside the ROIs
            idx = LOCB(i);   % find the voxel index of the OFC voxels in the gray matter voxel space
            tmp = (use_dat(:,idx)' * use_dat)./ bin_size; % Pearson correlation on z-scored values (with each voxel)

            % !! changes are made here to exclude voxels from another ROI
            % but pay attention not to rerun this line of code
            tmp(LOCB_1) = [];
            
            % updated this block of code (simply use "mean" so no need to
            % adjust number of voxels being divided)
            if strcmp(Option,'Unsigned')
                cors(i,t) = mean(abs(atanh(tmp(tmp<1)))); %UNSIGNED
            elseif strcmp(Option,'Signed')
                cors(i,t) = mean(atanh(tmp(tmp<1))); % SIGNED
            else
                cors(i,t) = mean((atanh(tmp(tmp<1))).^2); % Squared
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

%% updated data organization way
% organize this into an array
% dimension: nROIs * nsubs * nSess * ntimes
% so it's easier to load this into R and get both sess and TMS info
% together

dat = zeros(length(ROIfiles),nsubs,nSess,nruns*2);
subctr = 0;
for subj = Subs'  
    subctr = subctr + 1;
    subno = sprintf('Sub%d',subj);
    for roi = 1:length(ROIlabels)
    for s = 1:nSess
        for r = 1:nruns
            load(fullfile(External, 'GlobalConnResults', conn_name, subno, Option,sprintf('%s_cors_Sess%d_Run%d.mat',ROIlabels{roi},s,r)))        
            dat(roi,subctr,s,2*r-1) = avg_cors(1);
            dat(roi,subctr,s,2*r) = avg_cors(2);
        end
    end
    end
end

save(fullfile(maskpath,'ROI_analysis','GC_2bins_dat_partial.mat'), 'dat'); 


