
% The purpose of this script is to use coregister-reslice function to
% match the voxel size/dim to the reference image

% This script should be run under its parent directory

clc; clear; close all
matlabbatch = [];

% SPM function - Coregister: Reslice 
% ref: Image defining space
% source: Images to reslice
matlabbatch{1}.spm.spatial.coreg.write.ref = {'/Volumes/QF10TB/SPEEDTMS_results/TempFuncCopied/Sub1/func/Day3/Run1/wfSPEEDTMS_01_D3-0002-00001-000001-01.nii'};
matlabbatch{1}.spm.spatial.coreg.write.source = {'/Users/qingfangliu/Dropbox/KahntLab/MakeTissueMasks/gm_0.1.nii,1'
                                                '/Users/qingfangliu/Dropbox/KahntLab/MakeTissueMasks/wm_0.9.nii,1'
                                                '/Users/qingfangliu/Dropbox/KahntLab/MakeTissueMasks/csf_0.9.nii,1'
                                                 };
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0; % use nearest neighbor interpolation
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0]; % no wrap
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

% save job
fname = 'Job.mat';
save(fname, 'matlabbatch');

% run job
spm_jobman('run', matlabbatch)

%% to check various things

% Check differences between images before and after 
Y = spm_read_vols(spm_vol('rTPM.nii'));
size(Y)
% 79*95*79*6

% view image in matlab
figure
imagesc(Y(:,:,10,1)) % display the axial slice where z=10 in the 1st volume

Y1 = spm_read_vols(spm_vol('gm_0.1.nii')); % load mask before resliced
size(Y1)   % 79*95*79
yy1 = reshape(Y1,1,[]);
tabulate(yy1) % should contain only 1 and 0

Y2 = spm_read_vols(spm_vol('rgm_0.1.nii')); % load mask after resliced
size(Y2)   % 104*96*58
yy2 = reshape(Y2,1,[]);
tabulate(yy2) % should contain only 1 and 0

V3 = spm_vol('/Users/qingfangliu/Experiment/Sub1/func/Day3/Run1/meanfSPEEDTMS_01_D3-0002-00001-000001-01.nii,1');
Y3 = spm_read_vols(V3);  
size(Y3) % 104*96*58

V4 = spm_vol('/Users/qingfangliu/Experiment/Sub1/func/Day3/Run1/s2mmwfSPEEDTMS_01_D3-0002-00001-000001-01.nii');
Y4 = spm_read_vols(V4);
size(Y4) % 79*95*79

% size of Y1: 79*95*79
% size of Y2,Y3: 104*96*58
% so reslicing alters the size of the mask to match the reference image






