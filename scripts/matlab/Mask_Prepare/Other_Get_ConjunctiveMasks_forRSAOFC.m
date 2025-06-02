
% updated 3/29/23
% the locations where these func masks were saved have been changed and
% renamed.


clear; clc
% conjunction of two masks

P = mfilename('fullpath');
Path = fileparts(fileparts(P));
ROIdir = fullfile(Path,'ROIs'); % path of ROIs directory

%% Change this part

ffiles = {'LOFC_med_1e-3.hdr','rOFC_med_1e-3.hdr'};
afiles = {'rOFC_med_L.hdr','rOFC_med_R.hdr'};

for i = 1:2
    
    name1 = ffiles{i};
    name2 = afiles{i};
    
    newname = sprintf('Conj_%s',name1); % add 'Conj' to the old functional mask name
    newname = strrep(newname,'hdr','nii'); % replace 'hdr' with 'nii'

%%
% read mask1
tpmhdr1 = spm_vol(fullfile(ROIdir,'FromRSA',name1));
[Y1,~] = spm_read_vols(tpmhdr1);
sz1 = size(Y1);
idx1 = find(Y1 == 1);
assert(~isempty(idx1))

% read mask2
tpmhdr2 = spm_vol_nifti(fullfile(ROIdir,'aal_OFC',name2));
[Y2,~] = spm_read_vols(tpmhdr2);
sz2 = size(Y2);
Y2 = double(Y2);
idx2 = find(Y2);
%idx2 = find(Y2 ~= 0); % don't know why, but ==1 doesn't work (tabulate shows both 1 and 0)
assert(~isempty(idx2))

% union of two idx for the new mask
assert(isequal(sz1, sz2))
idx = intersect(idx1,idx2);

% create new mask
newhdr = tpmhdr1;
newvol = zeros(sz1);
newhdr.fname = fullfile(ROIdir,'FromRSA',newname);
newvol(idx) = 1;
spm_write_vol(newhdr,newvol);

fprintf('Conjunction successfully! \n')

end

