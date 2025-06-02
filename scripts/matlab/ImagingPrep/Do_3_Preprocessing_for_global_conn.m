
%% This script is to re-preprocess the task fMRI images
% This code was updated 9/2/2023
% purpose is to 
% (1) change the voxel size to (3*3*3)
% (2) change the smoothing size to (6*6*6)

% it starts from the normalization write step (not from the beginning realignment)
% no need to redo the realignment & coregistration
% as those have been written into the header files of the 'f' images
% Also no need to redo anatomical image normalization


tpmfile = '/Applications/spm12/tpm/TPM.nii';
ImagingDir = fullfile(External,'ImagingData');

%% To organize functional images 
func_ctr = 0;      % this counts functional images across runs and sessions
clear fimages      % this prepares to store all functional images 

% change where functional scans are saved
subdirfunc = fullfile(ImagingDir,subno,'func');
sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));
nSess = length(sessdirsfunc);

for ss = 1:nSess
    sess_name = sessdirsfunc(ss).name;
    for r = 1:nruns
        path = fullfile(subdirfunc, sess_name, sprintf('Run%d',r));
        n = dir(fullfile(path, sprintf('f*.nii'))); % get all files starting with 'f'
        for i = 1:length(n)
            func_ctr = func_ctr + 1; % count functional image numbers
            fname = fullfile(path, n(i).name); % current functional image
            fimages{func_ctr,1} = sprintf('%s,1', fname); % all functional images across sessions and runs
        end
    end
end


%% To organize anatomical image
anatpath = fullfile(ImagingDir, subno, 'anat');   
n = dir(fullfile(anatpath, sprintf('s*.nii')));
anatfile = n.name;                        
anat_name = fullfile(anatpath, anatfile);

%% Normalize
deformation_field = fullfile(anatpath, sprintf('y_%s', anatfile)); % get the deformation field

clear matlabbatch
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deformation_field};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = fimages; % all functional images
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3 3 3]; % change it to 3*3*3mm voxels
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w3'; % use a different prefix

%% Spatial smoothing
func_ctr = 0; 
clear wfilename

for ss = 1:nSess    
    sess_name = sessdirsfunc(ss).name;
    for r = 1:nruns
        path = fullfile(subdirfunc,sess_name,sprintf('Run%d', r) );
        n = dir(fullfile(path, sprintf('f*.nii')));
        for i = 1:length(n)
            func_ctr = func_ctr + 1;
            fname = fullfile(path, sprintf('w3%s',n(i).name)); % find all images starting with 'w3'
            wfilename{func_ctr,1} = sprintf('%s,1', fname);
        end
    end
end

matlabbatch{2}.spm.spatial.smooth.data = wfilename;
matlabbatch{2}.spm.spatial.smooth.fwhm = [6 6 6]; % smoothing size 6*6*6
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's6mm';


%% run the jobs
jobpath = fullfile(External, 'jobs', subno);
if ~exist(jobpath,'dir')
    mkdir(jobpath);
end
jobname = fullfile(jobpath, sprintf('Preprocessing_for_global_conn_%s.mat',subno));
save(jobname, 'matlabbatch');

%% run the jobs
spm_jobman('run', matlabbatch);    % SPM jobman will run the batch



