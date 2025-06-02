
% updated code with 2mm voxel size

clear; clc;

%% First: reslice TPM to match current data size & dim

sourceImg = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/Global_connectedness/TPM.nii'; % copied from /Applications/spm12/tmp/TPM.nii
refImg = '/Users/qingfangliu/Experiment/Sub1/func/Day3/Run1/s6wfSPEEDTMS_01_D3-0002-00001-000001-01.nii,1'; % Normalized data in 2*2*2 space

% reslice
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref = {refImg};
matlabbatch{1}.spm.spatial.coreg.write.source = {sourceImg};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4; % use 4th degree spline
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0]; % no wrap
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r2mm';

% run job
spm_jobman('run', matlabbatch)

% a file 'rTPM.nii' should be generated 

%% Second: threshold rTMP to get masks

hdr = spm_vol('/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/Global_connectedness/r2mmTPM.nii');
nhdr = hdr(1);

% gray matter
x = spm_read_vols(hdr(1)); 
y = zeros(size(x));
y(x > 0.1) = 1; % thereshold at 10%
nhdr.fname = 'gm_0.1_2mm.nii';
spm_write_vol(nhdr,y);

%% Check size

% reference image
y1 = spm_read_vols(spm_vol(refImg));
size(y1)

% gray matter
y2 = spm_read_vols(spm_vol('gm_0.1_2mm.nii'));
size(y2)

% gray matter should have the same size as the reference image






