

% Do reslice first on the original tpm and then select volume and do
% thresholding
% in this way we only slice once
% and we can use the 4th degree spline method for interpolation

% 10%-90% probabilistic gray matter mask
% - Threshold into binary values (0,1) from continuous values 
% - reslice them (when choosing interpolation method using nearest
% neighbor, not do spline)

% - Another way to do this is to reslice them first and threshold later (no
% need to worry about interpolation method anymore)

clear
sourceImg = 'TPM.nii';

% reslice
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref = {'C:\Data\SPEEDTMS\Piloting\Sub2\func\Run1\meanfSPEEDTMS_test-0004-00001-000001-01.nii,1'};
matlabbatch{1}.spm.spatial.coreg.write.source = {sourceImg};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4; % use 4th degree spline
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0]; % no wrap
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r1';

% % save job
%fname = 'ResliceJobAllInOne.mat';
%save(fname, 'matlabbatch');

% run job
spm_jobman('run', matlabbatch)


hdr = spm_vol('r1TPM.nii');
nhdr = hdr(1);

% gray matter
x = spm_read_vols(hdr(1)); 
y = zeros(size(x));
y(x > 0.1) = 1; % thereshold at 10%
nhdr.fname = 'New_gm_0.1.nii';
spm_write_vol(nhdr,y);


%% compare the gray matter mask from this script with what we got by reslicing first
v = spm_vol('New_gm_0.1.nii');
y = spm_read_vols(v);

old_gray_V = spm_vol('rgm_0.1.nii');
old_gray_y = spm_read_vols(old_gray_V);

tabulate(reshape(y,1,[]))
tabulate(reshape(old_gray_y,1,[]))


