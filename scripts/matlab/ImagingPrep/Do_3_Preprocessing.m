
%% This script is to preprocessing the task fMRI images

%cd(fullfile(studydir,subno)) % change working directory so ps file can be saved properly
nruns = 3; % number of runs
tpmfile = '/Applications/spm12/tpm/TPM.nii';
ImagingDir = fullfile(External,'ImagingData');

%% REALIGNMENT OF WHOLE BRAIN

subdirWB = fullfile(ImagingDir,subno,'WB');
sessdirsWB = dir(fullfile(subdirWB,'Day*')); % only use folders starting from 'Day'
nSess = length(sessdirsWB); % number of sessions
scntr = 0; % count session numbers

WBfilename = [];
for ss = 1:nSess    
    scntr = scntr + 1;
    sess_name = sessdirsWB(ss).name;
    path = fullfile(subdirWB, sess_name); 
    n = dir(fullfile(path, 'f*.nii'));  
    if ss == 1
       wbpath = path;  % used in coregistration
       wbfile = n(1).name;   
    end

    for i = 1:length(n)
        fname = fullfile(path, n(i).name);
        WBfilename{scntr}{i,1} = sprintf('%s,1', fname);        
    end
end

matlabbatch = [];
matlabbatch{1}.spm.spatial.realign.estwrite.data = WBfilename;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1]; % only to reslice the mean image, not all images
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';


%% Functional scan realignment
rcntr = 0;         % rnctr is run counter across runs and sessions
func_ctr = 0;      % this counts functional images across runs and sessions
clear filename fimages

subdirfunc = fullfile(ImagingDir,subno,'func');
sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));
nSess = length(sessdirsfunc);

for ss = 1:nSess
    sess_name = sessdirsfunc(ss).name;

    for r = 1:nruns
        rcntr = rcntr + 1;   % For next run. After 3rd run the next run (session-2 1st run) will be 4th run, total of 6 runs for each subject
        path = fullfile(subdirfunc, sess_name, sprintf('Run%d',r));
        n = dir(fullfile(path, sprintf('f*.nii')));
        if rcntr == 1
           funcpath = path; % used in coregistration
           funcfile = n(1).name;
        end

        for i = 1:length(n)
            func_ctr = func_ctr + 1;
            fname = fullfile(path, n(i).name);
            func_filename{rcntr}{i,1} = sprintf('%s,1', fname);     % because each run is in its own curly bracket, first defined the rcntr to move through the runs and {i,1} define sthe strcture of the file witha 1 after comma
            fimages{func_ctr,1} = sprintf('%s,1', fname);            % for coregistration we dot need brakes in between different runs so we created ii and made a list of images without rcntr
        end
    end
end

matlabbatch{2}.spm.spatial.realign.estwrite.data = func_filename;           % do not need curly brackets here because already defined curly brackets
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [0 1];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';


%% COREGISTRATION

anatpath = fullfile(ImagingDir, subno, 'anat');   % define path for anatomical image
n = dir( fullfile(anatpath, sprintf('s*.nii')));
anatfile = n.name;                        

anat_name = fullfile(anatpath, anatfile);
wb_name = fullfile(wbpath, sprintf('mean%s', wbfile));          % since mean file does not exist, we are creating a name for the mean file
func_name = fullfile(funcpath, sprintf('mean%s', funcfile));

% Anatomical and whole brain coreg
matlabbatch{3}.spm.spatial.coreg.estimate.ref = {anat_name};
matlabbatch{3}.spm.spatial.coreg.estimate.source = {wb_name};
matlabbatch{3}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

% Whole brain and functional scan coreg
matlabbatch{4}.spm.spatial.coreg.estimate.ref = {wb_name};
matlabbatch{4}.spm.spatial.coreg.estimate.source = {func_name};
matlabbatch{4}.spm.spatial.coreg.estimate.other = fimages; % func images in all runs of all sessions
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];


%% Normalize
matlabbatch{5}.spm.spatial.normalise.estwrite.subj.vol = {anat_name};
matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample =  {anat_name};
matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.tpm = {tpmfile};
matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{5}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{5}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
matlabbatch{5}.spm.spatial.normalise.estwrite.woptions.interp = 4;

% Normalise write
fname = fullfile(anatpath, sprintf('y_%s', anatfile));
matlabbatch{6}.spm.spatial.normalise.write.subj.def = {fname};
matlabbatch{6}.spm.spatial.normalise.write.subj.resample = fimages;
matlabbatch{6}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{6}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{6}.spm.spatial.normalise.write.woptions.interp = 4;


%% Spatial smoothing

func_ctr = 0; 
clear filename
subdirfunc = fullfile(ImagingDir,subno,'func');
sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));
nSess = length(sessdirsfunc);

for ss = 1:nSess    
    sess_name = sessdirsfunc(ss).name;
    for r = 1:nruns
        path = fullfile(subdirfunc,sess_name,sprintf('Run%d', r) );
        n = dir(fullfile(path, sprintf('f*.nii'))); % look for f* not wf* because these files are not there yet
        for i = 1:length(n)
            func_ctr = func_ctr + 1;
            fname = fullfile(path, sprintf('w%s',n(i).name));
            filename{func_ctr,1} = sprintf('%s,1', fname);
        end
    end
end

matlabbatch{7}.spm.spatial.smooth.data = filename;
matlabbatch{7}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{7}.spm.spatial.smooth.dtype = 0;
matlabbatch{7}.spm.spatial.smooth.im = 0;
matlabbatch{7}.spm.spatial.smooth.prefix = 's6';

matlabbatch{8}.spm.spatial.smooth.data = filename;
matlabbatch{8}.spm.spatial.smooth.fwhm = [2 2 2];
matlabbatch{8}.spm.spatial.smooth.dtype = 0;
matlabbatch{8}.spm.spatial.smooth.im = 0;
matlabbatch{8}.spm.spatial.smooth.prefix = 's2mm';

%% run the jobs
% jobpath = fullfile(studydir, subno, 'jobs');
% if ~exist(jobpath,'dir')
%     mkdir(jobpath);
% end
% jobname = fullfile(jobpath, sprintf('Preprocessing_%s.mat',subno));
% save(jobname, 'matlabbatch');

%spm fmri
spm_jobman('run', matlabbatch);    % SPM jobman will run the batch

