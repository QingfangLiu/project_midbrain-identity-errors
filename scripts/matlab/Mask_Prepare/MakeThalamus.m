

% to create thalamus ROI from the functional ROI & anatomical ROI

ROIpath = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs';

% combine left and right into one

name1 = fullfile(ROIpath,'atlas','Left Thalamus Proper.nii');
name2 = fullfile(ROIpath,'atlas','Right Thalamus Proper.nii');
name3 = fullfile(ROIpath,'WB_corrected_masks','Thalamus_b.img');

%% reslice 
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref = {name3};
matlabbatch{1}.spm.spatial.coreg.write.source = {name1
                                                 name2
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


%%
name1 = fullfile(ROIpath,'atlas','rLeft Thalamus Proper.nii');
name2 = fullfile(ROIpath,'atlas','rRight Thalamus Proper.nii');

%%

% read mask1
tpmhdr1 = spm_vol(name1);
[Y1,~] = spm_read_vols(tpmhdr1);
idx1 = find(Y1 == 1); % 1 for included 
sz = size(Y1);

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
newhdr.fname = fullfile(ROIpath,'WB_corrected_masks','Conj_Thalamus.nii');
newvol(idx) = 1;
spm_write_vol(newhdr,newvol);



