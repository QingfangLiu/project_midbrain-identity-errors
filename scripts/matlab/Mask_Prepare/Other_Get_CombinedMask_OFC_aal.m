
% combines aal OFC masks

clear; clc
Dir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs/aal_OFC';

%%
% ant OFC
name1 = 'OFC_med_R.img';
name2 = 'OFC_med_L.img';

%%
% med OFC
%name1 = fullfile(parentDir,'ROIs','aal_OFC','OFC_med_L.hdr');
%name2 = fullfile(parentDir,'ROIs','aal_OFC','OFC_med_R.hdr');
%newname = 'OFC_med_b.nii'; 

%% reslice 
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref = {'/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs/LPFC_func_1e-7.img'};
matlabbatch{1}.spm.spatial.coreg.write.source = {fullfile(Dir,name1)
                                                 fullfile(Dir,name2)
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
