
%% This script is to preprocessing the task fMRI images to be used for decoding
% (1) motion correction (reslice all images): generating 'rf' images
% (2) coregistration (info passed to headers of rf images)
% (3) 2mm FWHM smoothing
% no normalization 

ImagingDir = fullfile(External,'ImagingData');

%% Functional scan realignment
rcntr = 0;                          % rnctr is run counter across runs and sessions
func_ctr = 0;                       % this counts functional images across runs and sessions
subdirfunc = fullfile(ImagingDir,subno,'func');
sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));
clear fimages

for ss = 1:nSess
    sess_name = sessdirsfunc(ss).name;
    for r = 1:nruns
        rcntr = rcntr + 1;   
        path = fullfile(subdirfunc, sess_name, sprintf('Run%d',r));
        n = dir(fullfile(path, sprintf('f*.nii')));
        for i = 1:length(n)
            func_ctr = func_ctr + 1;
            fname = fullfile(path, n(i).name);   
            fimages{func_ctr,1} = sprintf('%s,1', fname); 
        end
    end
end

clear matlabbatch
matlabbatch{1}.spm.spatial.realign.write.data = fimages;                    % images are organized as a big list, not by runs
matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 0];            % reslice all images (2), but not the mean image (already done, 0)
matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';


%% COREGISTRATION

% anatomical file
n = dir(fullfile(ImagingDir, subno, 'anat','s*.nii'));
anat_name = fullfile(n.folder, n.name);

% whole-brain file
n = dir(fullfile(ImagingDir,subno,'WB','Day*','mean*'));     
wb_name = fullfile(n.folder,n.name);

% function file
n = dir(fullfile(ImagingDir,subno,'func','Day*','Run*','mean*'));     
func_name = fullfile(n.folder,n.name);

% organize 'rf' func images for coregistration
% rf_files = dir(fullfile(ImagingDir,subno,'func','Day*','Run*','rf*'));
% rfimages = cell(length(rf_files),1);
% for i = 1:length(rf_files)
%     rfimages{i,1} = fullfile(rf_files(i).folder,rf_files(i).name);
% end


func_ctr = 0; 
clear filename
subdirfunc = fullfile(ImagingDir,subno,'func');
sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));

for ss = 1:nSess    
    sess_name = sessdirsfunc(ss).name;
    for r = 1:nruns
        path = fullfile(subdirfunc,sess_name,sprintf('Run%d',r));
        n = dir(fullfile(path, sprintf('f*.nii')));
        for i = 1:length(n)
            func_ctr = func_ctr + 1;
            fname = fullfile(path, sprintf('r%s',n(i).name));
            filename{func_ctr,1} = sprintf('%s,1', fname);
        end
    end
end

% Anatomical and whole brain coreg
matlabbatch{2}.spm.spatial.coreg.estimate.ref = {anat_name};
matlabbatch{2}.spm.spatial.coreg.estimate.source = {wb_name};
matlabbatch{2}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

% Whole brain and functional scan coreg
matlabbatch{3}.spm.spatial.coreg.estimate.ref = {wb_name};
matlabbatch{3}.spm.spatial.coreg.estimate.source = {func_name};
matlabbatch{3}.spm.spatial.coreg.estimate.other = filename; % func images in all runs of all sessions
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];


%% Spatial smoothing
matlabbatch{4}.spm.spatial.smooth.data = filename;
matlabbatch{4}.spm.spatial.smooth.fwhm = [2 2 2];
matlabbatch{4}.spm.spatial.smooth.dtype = 0;
matlabbatch{4}.spm.spatial.smooth.im = 0;
matlabbatch{4}.spm.spatial.smooth.prefix = 's2';

%% run the jobs
jobpath = fullfile(External, 'jobs', subno);
if ~exist(jobpath,'dir')
    mkdir(jobpath);
end
jobname = fullfile(jobpath, sprintf('Preprocessing_for_decoding_native_%s.mat',subno));
save(jobname, 'matlabbatch');


%%
spm_jobman('run', matlabbatch);

