
% The purpose of this script is to convert functional and anatomical dicom files into
% NifTi images and save them under folders:
% 'func' for task-related functional images, with a prefix 'f' 
% 'WB' for whole-brain functional images, with a prefix 'f'
% 'anat' for structural image, with a prefix 's'

subdir = fullfile(studydir,subno,'OutDicom');
sessdirs = dir(fullfile(subdir,'Day*'));
nruns = 3; % number of runs
nSess = 3; % number of sessions


%% task runs - convert Dicom to NifTi images
% only for sessions day2 & day3 (also day4 for subj1)
for ss = 2:3
for r = 1:nruns

sess_name = sessdirs(ss).name;
runname = sprintf('*r%d_22_58_2000_2mmiso',r);
rundir = dir(fullfile(subdir,sess_name,runname));
dicomnames = dir(fullfile(rundir.folder,rundir.name));
dicomnames = dicomnames(~[dicomnames.isdir]); 

% folders to save functional images
func_dir = fullfile(studydir,subno,'func',sess_name,sprintf('Run%d',r));
if ~exist(func_dir,'dir')
    mkdir(func_dir)
end

if ~isempty(dir(fullfile(func_dir,'f*.nii')))
    fprintf('\nDICOM files already converted for run %d\n',r);
else
    fprintf('%d DICOM files found for run %d \n',length(dicomnames),r);
    fprintf('Files converted: \n');
    for i = 1:length(dicomnames)
        if ~mod(i,50)
            fprintf('DICOM to NIFTI %d%% done \n',round(i/length(dicomnames)*100));
        else
        end
        spm_dicom_convert(spm_dicom_headers(fullfile(dicomnames(i).folder,dicomnames(i).name)),'all','flat','nii',func_dir);
    end
end

end % end of runs
end % end of sessions


%% WB functional images
% each day has WB scan

for ss = 1:nSess
    fprintf('Session %d/%d\n',ss,nSess)
    sess_name = sessdirs(ss).name;
    WBdir = dir(fullfile(subdir,sess_name,'*WB*'));
    
    dicomnames = dir(fullfile(WBdir.folder,WBdir.name));
    dicomnames = dicomnames(~[dicomnames.isdir]); 

    % folders to save functional images
    func_dir = fullfile(studydir,subno,'WB',sess_name);
    if ~exist(func_dir,'dir')
        mkdir(func_dir)
    end

    if ~isempty(dir(fullfile(func_dir,'f*.nii')))
         fprintf('\nDICOM files already converted for WB\n');
    else

    fprintf('%d DICOM files found for WB \n',length(dicomnames));
    fprintf('Files converted: \n');
    for i = 1:length(dicomnames)
        spm_dicom_convert(spm_dicom_headers(fullfile(dicomnames(i).folder,dicomnames(i).name)),'all','flat','nii',func_dir);
    end
    end

end

%% anatomical image

% MPRAGE_SAG folder (T1-weighted image)
% magnetization-prepared rapid gradient-echo
% in the sagittal plane
% A potential issue here is that T1 DCM files might have been copied twice

dicompath = dir(fullfile(subdir,'Day1','MPR*'));  % ONLY for day1
adir = fullfile(dicompath.folder,dicompath.name);
T1dir = fullfile(studydir,subno,'anat');
if ~exist(T1dir,'dir')
    mkdir(T1dir)
end

% check if NifTi file has been converted
if isempty(dir(fullfile(T1dir,'s*.nii')))
    dicomnames = dir(adir); % all the Dicom files
    dicomnames = dicomnames(~[dicomnames.isdir]);
    sdicomname = strings(length(dicomnames),1);
    for i=1:length(dicomnames)
        sname = fullfile(adir, dicomnames(i).name);
        sdicomname(i) = sname;        
    end
    fprintf('%d DICOM files found for Anatomy\n',length(dicomnames));
    fprintf('Files converted: ');   
    spm_dicom_convert(spm_dicom_headers(sdicomname),'all','flat','nii',T1dir);
    fprintf('\n');
else
    fprintf('DICOM files already converted for T1\n');
end

