
% extract global connectedness values using anatomical ROIs used for TMS

%% prepare
clc; clear;
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
nscans = 430; % number of scans in each run

conn_name = 'ConnectednessMap_3mm'; % change where the whole result is saved
Option = 'Squared'; % Unsigned or Signed or Squared
GCDir = fullfile(parentDir,'Global_connectedness','FinalResBackUp',conn_name);

%% find ROIs to use
% using TMS seed region OFC & target region LPFC
roiname = {'SeedRegion_OFC_bilateral.nii','SeedRegion_lOFC.nii','SeedRegion_rOFC.nii',...
    'TargetRegion_LPFC_bilateral.nii','TargetRegion_lLPFC.nii','TargetRegion_rLPFC.nii'};
roi_labels = {'OFC','lOFC','rOFC','LPFC','lLPFC','rLPFC'};

% copy those ROIs to the GC analysis folder
GCROI_dir = fullfile(parentDir,'Global_connectedness','ROIs');
if ~exist(GCROI_dir,'dir')
    mkdir(GCROI_dir)
end

% organize all these roi images
sourceImg = cell(length(roiname),1); % all ROIs in 2*2*2 space
for i = 1:length(roiname)
    roi_file = fullfile(parentDir,'ROIs',roiname{i});
    new_file = fullfile(GCROI_dir,roiname{i});
    copyfile(roi_file,new_file)
    sourceImg{i,1} = new_file;
end

%% reslice ROIs from 2*2*2 to the 3*3*3 voxel space
refImg = fullfile(parentDir,'Global_connectedness','FinalResBackUp',conn_name,'Sub1',...
    Option,'ConnectednessMap_3mm_Day3_Run1.nii'); % Normalized data in 3*3*3 size
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref = {refImg};
matlabbatch{1}.spm.spatial.coreg.write.source = sourceImg;
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0; % use nearst neighbor
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0]; % no wrap
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r3mm';
save(fullfile(GCROI_dir,'ResliceJob.mat'), 'matlabbatch');
spm_jobman('run', matlabbatch)

%% read ROI idx on 3*3*3 voxel space
roiidx = cell(length(roiname),1);
for i = 1:length(roiname)
    roi_vol = spm_read_vols(spm_vol(fullfile(GCROI_dir,sprintf('r3mm%s',roiname{i}))));
    roiidx{i} = find(roi_vol);
end

% test code using the functional ROI, but ignore this later
%roi_vol = spm_read_vols(spm_vol(fullfile(GCROI_dir,'func_lOFC.hdr')));
%roiidx{6} = find(roi_vol);
%roiname = [roiname,'func_lOFC'];

%% want to extract GC voxel values for each run and for each subject
clear GC_voxels MeanGC
subctr = 0;

for subj = Subs'   
    subctr = subctr + 1;
    subno = sprintf('Sub%d',subj);
    fprintf('Extracting GC values from %s\n',subno)
    subdir = fullfile(studydir,subno);
    subdirfunc = fullfile(subdir,'func');
    sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));
    for ss = 1:nSess 
        sess_name = sessdirsfunc(ss).name; % get session name
    for r = 1:nruns
        TempMap = fullfile(GCDir,subno,Option,sprintf('%s_%s_Run%d.nii', conn_name, sess_name, r));
        for i = 1:length(roiidx) % loop over each ROI
            bvol = spm_read_vols(spm_vol_nifti(TempMap));
            GC_voxels{subctr,i}(ss,r,:) = bvol(roiidx{i}); % make sure ROI and map are in the same voxel space
            MeanGC(i,ss,r,subctr) = nanmean(bvol(roiidx{i})); % keep the mean collapse across individual voxel values
        end
    end
    end
end

%% save GC voxel values
SaveDir = fullfile(parentDir,'Global_connectedness','ROI_analysis');
if ~exist(SaveDir,'dir')
    mkdir(SaveDir)
end

matname = fullfile(SaveDir,sprintf('Extracted_GC_values_%s_%s.mat',conn_name,Option));
save(matname,'GC_voxels','MeanGC','roi_labels')




         