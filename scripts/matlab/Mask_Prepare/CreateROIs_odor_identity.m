
% This script creates ROIs using MNI coordinates
% using coordinates from odor identity sensitive regions from Howard 2015
% PNAS

clear; close all; clc

%% define the work path
ROIpath = '/Users/qingfangliu/Library/CloudStorage/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs';
IndepROIpath = fullfile(ROIpath,'Indep');  % Independent ROI path
if ~exist(IndepROIpath,'dir')
   mkdir(IndepROIpath);
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROIs = {...
    'lOFC',[42,36,-16];...
    'ACC',[6,38,16];...
    'HPC',[38,-16,-16];...
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


