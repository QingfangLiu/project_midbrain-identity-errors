
% this script gets conjunction of functional MB masks (WB FWE corrected) and anat MB masks

clear; clc
% conjunction of two masks

P = mfilename('fullpath');
Path = fileparts(fileparts(P));
ROIdir = fullfile(Path,'ROIs'); % path of ROIs directory

%% Change this part

func_mb = {'MB_b.hdr','MB_l.hdr','MB_r.hdr'};
n_func_mb = length(func_mb);
anat_mb = {'rmean_fullMB_point5.nii','rmean_fullMB_point8.nii'};
anat_labels = {'p5','p8'};
n_anat_mb = length(anat_mb);

for p = 1:n_anat_mb
for x = 1:n_func_mb

newname = sprintf('Conj_%s_%s',anat_labels{p},func_mb{x}); % add 'Conj' to the old functional mask name
newname = strrep(newname,'hdr','nii'); % replace 'hdr' with 'nii'

%%
% read mask1
tpmhdr1 = spm_vol(fullfile(ROIdir,'WB_corrected_masks',func_mb{x}));
[Y1,~] = spm_read_vols(tpmhdr1);
sz = size(Y1);
idx1 = find(Y1 == 1);

% read mask2
tpmhdr2 = spm_vol_nifti(fullfile(ROIdir,anat_mb{p}));
[Y2,~] = spm_read_vols(tpmhdr2);
idx2 = find(Y2 == 1);

% union of two idx for the new mask
idx = intersect(idx1,idx2);

% create new mask
newhdr = tpmhdr1;
newvol = zeros(sz);
newhdr.fname = fullfile(ROIdir,'WB_corrected_masks',newname);
newvol(idx) = 1;
spm_write_vol(newhdr,newvol);

fprintf('Conjunction successfully! \n')

end
end


