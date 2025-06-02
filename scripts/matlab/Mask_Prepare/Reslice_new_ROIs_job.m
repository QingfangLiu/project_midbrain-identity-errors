
% this script reslices 

clear; clc;

ResDir = '/Volumes/QF10TB/SPEEDTMS_results'; 
ROIpath = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs';

%%
refImg = fullfile(ROIpath,'GlobalConn','func_OFC.img');
%-----------------------------------------------------------------------
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref = {refImg};
matlabbatch{1}.spm.spatial.coreg.write.source = {
                                                 %'/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs/Indep/AMG_b.nii,1'
                                                 %'/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs/Indep/AMG_l.nii,1'
                                                 %'/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs/Indep/AMG_r.nii,1'
                                                 '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ROIs/WB_corrected_masks/Conj_Thalamus.nii,1'
                                                 };
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r3mm';

% run job
spm_jobman('run', matlabbatch)