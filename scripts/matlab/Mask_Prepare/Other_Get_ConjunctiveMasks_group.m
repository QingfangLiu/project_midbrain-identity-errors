

clear; clc
parentDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis';

%%
name1 = fullfile(parentDir,'Masks','Group_mask.nii');
% name2 = '/Users/qingfangliu/Dropbox/KahntLab/MakeTissueMasks/gm_0.1.nii';

name2 = fullfile(parentDir,'ROIs','rTPM_point1_NoCerebellum.nii');
name3 = fullfile(parentDir,'ROIs','rmean_fullMB_point5.nii');

% all three ROIs are in the standard MNI space
% we want the mask to include both gray matter and MB voxels, and then
% conjuncted with the group mask from subjects to take care of any missing
% voxels

%newname = 'SPEEDTMS_group_mask.nii'; % define new name for the conjunctive mask
newname = 'SPEEDTMS_external_mask.nii'; % use a new name

% read mask1
tpmhdr1 = spm_vol(name1);
[Y1,~] = spm_read_vols(tpmhdr1);
sz = size(Y1);
idx1 = find(Y1 == 1); % 1 for included 

% read mask2
tpmhdr2 = spm_vol(name2);
[Y2,~] = spm_read_vols(tpmhdr2);
sz2 = size(Y2);
idx2 = find(Y2 == 255); % 255 for included

% read mask3
tpmhdr3 = spm_vol(name3);
[Y3,~] = spm_read_vols(tpmhdr3);
sz3 = size(Y3);
idx3 = find(Y3 == 1); 

% take union of gray matter and midbrain, then intersect with group mask
idx = intersect(idx1,union(idx2,idx3));

% create new mask
newhdr = tpmhdr1;
newvol = zeros(sz);
newhdr.fname = fullfile(parentDir,'ROIs',newname);
newvol(idx) = 1;
spm_write_vol(newhdr,newvol);

fprintf('Conjunction successfully! \n')






