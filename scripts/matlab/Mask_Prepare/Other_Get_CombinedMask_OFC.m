

clear; clc
% combine two masks

parentDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis';

name1 = fullfile(parentDir,'ROIs','SeedRegion_lOFC.nii');
name2 = fullfile(parentDir,'ROIs','SeedRegion_rOFC.nii');
name3 = fullfile(parentDir,'ROIs','SPEEDTMS_group_mask.nii');
  
newname = 'SeedRegion_OFC_bilateral.nii'; % define new name for the conjunctive mask

% read mask1
tpmhdr1 = spm_vol(name1);
[Y1,~] = spm_read_vols(tpmhdr1);
sz = size(Y1);
idx1 = find(Y1 == 1); % 1 for included 

% read mask2
tpmhdr2 = spm_vol(name2);
[Y2,~] = spm_read_vols(tpmhdr2);
idx2 = find(Y2 == 1); % 1 for included

% read mask3
tpmhdr3 = spm_vol(name3);
[Y3,~] = spm_read_vols(tpmhdr3);
idx3 = find(Y3 == 1); % 1 for included

% union of two idx for the new mask
idx = union(idx1,idx2);
idx = intersect(idx,idx3);

% create new mask
newhdr = tpmhdr1;
newvol = zeros(sz);
newhdr.fname = fullfile(parentDir,'ROIs',newname);
newvol(idx) = 1;
spm_write_vol(newhdr,newvol);

fprintf('Conjunction successfully! \n')


