
% This script will be run on the laptop to get bias corrected image before
% doing TMS motor threshold
% the data will be temporarily saved in the laptop download folder

clear
close all

tic

%% Convert DAY1 DICOMs to nifti (if it hasn't been done already)

studydir = '/Users/qingfangliu/Downloads/SPEEDTMS/Experiment';
sub = 46;
subno = sprintf('Sub%d',sub);
subdir = fullfile(studydir,subno);

dpath = dir(fullfile(subdir, '*MPRAGE*'));
dicompath = fullfile(subdir, dpath(1).name);

% where to save the nifti file
outdir = fullfile(subdir,'TMS','T1');

if ~isempty(dir(fullfile(outdir,'s*.nii')))
    
    fprintf('\nDICOM files already converted\n');
else
    
    if ~exist(outdir,'dir')
        mkdir(outdir);
    else
    end
     
    dicomnames = dir(fullfile(dicompath,'*'));
    sdicomname = strings(length(dicomnames),1);
     for i = 3:length(dicomnames)
        sname = fullfile(dicompath, dicomnames(i).name);
        sdicomname(i) = sname; 
     end
    fprintf('%d DICOM files found for Anatomy\n',length(dicomnames)-2);
    fprintf('Files converted: ');

    spm_dicom_convert(spm_dicom_headers(sdicomname),'all','flat','nii',outdir);
    fprintf('\n');
end

%% do bias correction with the Nifti file

% set filenames and paths
subfolder = outdir; % path with anatomical data
anat = dir(fullfile(subfolder,'s*'));

% change folder to anatomical folder
cd(subfolder)

fname = fullfile(subfolder, anat(1).name);
matlabbatch{1}.spm.spatial.preproc.channel.vols = {fname};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
ent = [1, 1, 2, 3, 4, 2];
for i = 1:6
    matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm = {fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%01d', i))};
    matlabbatch{1}.spm.spatial.preproc.tissue(i).ngaus = ent(i);
    matlabbatch{1}.spm.spatial.preproc.tissue(i).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(i).warped = [0 0];
end
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

spm_jobman('run', matlabbatch);

for i = 1:6
    delete(fullfile(subfolder, sprintf('c%01d%s', i, anat(1).name)))
end
clear matlabbatch

    fprintf('\nCOPY FILES TO DRIVE!! :) \n');

toc



