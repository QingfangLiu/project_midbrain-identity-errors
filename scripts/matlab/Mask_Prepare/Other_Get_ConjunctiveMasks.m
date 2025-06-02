

clear; clc
% conjunction of two masks

P = mfilename('fullpath');
Path = fileparts(fileparts(P));
ROIdir = fullfile(Path,'ROIs'); % path of ROIs directory

%% Change this part
%name1 = fullfile(parentDir,'Masks','Group_mask.nii');
%name2 = fullfile(parentDir,'ROIs','rTPM_point1_NoCerebellum.nii');
    
%name1 = 'LeftMidbrain_1e-7.hdr';
%name1 = 'New_Midbrain_func_1e-7.hdr';
name1 = 'New_Right_Midbrain_func_1e-7.hdr';

name2 = 'rmean_fullMB_point8.nii'; % anat mask for midbrain

newname = sprintf('Conj_%s',name1); % add 'Conj' to the old functional mask name
newname = strrep(newname,'hdr','nii'); % replace 'hdr' with 'nii'

%%
% read mask1
tpmhdr1 = spm_vol(fullfile(ROIdir,name1));
[Y1,~] = spm_read_vols(tpmhdr1);
sz = size(Y1);
idx1 = find(Y1 == 1);

% read mask2
tpmhdr2 = spm_vol_nifti(fullfile(ROIdir,name2));
[Y2,~] = spm_read_vols(tpmhdr2);
idx2 = find(Y2 == 1);

% union of two idx for the new mask
idx = intersect(idx1,idx2);

% create new mask
newhdr = tpmhdr1;
newvol = zeros(sz);
newhdr.fname = fullfile(ROIdir,newname);
newvol(idx) = 1;
spm_write_vol(newhdr,newvol);

fprintf('Conjunction successfully! \n')


