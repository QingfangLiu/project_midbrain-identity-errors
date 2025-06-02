
% This script creates sphere amygdala ROIs using MNI coordinates

clear; close all; clc

%% define the work path
ROIpath = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs';
IndepROIpath = fullfile(ROIpath,'Indep');  % Independent ROI path
if ~exist(IndepROIpath,'dir')
   mkdir(IndepROIpath);
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROIs = {...
    'AMG_r',[16,-6,-14];...
    'AMG_l',[-20,-4,-20];...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% make spherical ROIs surrounding the coordinates, constrained by the tissue probability map
% these masks are in the normalized (MNI) space

rget = 5; 
% gray matter mask
tpmhdr = spm_vol_nifti(fullfile(ROIpath,'rTPM_point1_NoCerebellum.nii')); 

for i = 1:length(ROIs)
    
    filename = fullfile(IndepROIpath, sprintf('%s.nii',ROIs{i,1}));
    if ~exist(filename,'file')  %check if spherical ROIs have already been made
        [tpmvol,XYZ] = spm_read_vols(tpmhdr); 
        tpmidx = find(tpmvol>0); 
        sz = size(tpmvol);
        [MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3)); 
        ref_vox = round(sz/2);
        radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2); 
        radius_index = find(radii<rget) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3)); 
        coordidx = find(ismember(XYZ',ROIs{i,2},'rows')>0);
        idx = intersect(radius_index + coordidx, tpmidx);
        newhdr = tpmhdr;
        newvol = zeros(sz);
        newhdr.fname = filename;
        newvol(idx) = 1;
        spm_write_vol(newhdr,newvol); % write an image volume to disk
    else
        fprintf('This ROI has already been made\n');
    end
end


% combine left and right into one
name1 = fullfile(IndepROIpath,'AMG_r.nii');
name2 = fullfile(IndepROIpath,'AMG_l.nii');
newname = 'SeedRegion_OFC_bilateral.nii'; % define new name for the conjunctive mask

% read mask1
tpmhdr1 = spm_vol(name1);
[Y1,~] = spm_read_vols(tpmhdr1);
idx1 = find(Y1 == 1); % 1 for included 
sz = size(Y1);

% read mask2
tpmhdr2 = spm_vol(name2);
[Y2,~] = spm_read_vols(tpmhdr2);
idx2 = find(Y2 == 1); % 1 for included

% union of two idx for the new mask
idx = union(idx1,idx2);

% create new mask
newhdr = tpmhdr1;
newvol = zeros(sz);
newhdr.fname = fullfile(IndepROIpath,'AMG_b.nii');
newvol(idx) = 1;
spm_write_vol(newhdr,newvol);

