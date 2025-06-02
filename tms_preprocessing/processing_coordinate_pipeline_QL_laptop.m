
% This script should be run on Qingfang's laptop

clear; close all; clc

%% define the work path

homedir = '/Users/qingfangliu/Downloads/SPEEDTMS/Experiment';
ScriptFolder = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/TMS_preproc';

sub = 46;
subno = sprintf('Sub%d',sub);
subdir = fullfile(homedir,subno); % to save anat and func outputs

dpath = dir(fullfile(subdir, '*MPRAGE*')); % find the T1 scans
studyadir = fullfile(subdir,dpath(1).name); % T1 folder

dpath = dir(fullfile(subdir, '*MB*')); % find the rs scans
studyfdir = fullfile(subdir,dpath(1).name); % rs folder

subadir = fullfile(subdir,'TMS','T1'); % to save anat outputs
subfdir = fullfile(subdir,'TMS','rsfunc'); % to save func outputs
subJobdir = fullfile(subdir,'TMS','jobs'); % to save matlab batch jobs
if ~exist(subJobdir,'dir')
   mkdir(subJobdir);
else
end

maskdir = fullfile(subfdir,'ConnectivityMasks'); % path to masks for study
if ~exist(maskdir,'dir')
    mkdir(maskdir);
end
unnormdir = fullfile(subfdir, 'UnnormedMasks');  %path to native space masks for each subejct
if ~exist(unnormdir,'dir')
    mkdir(unnormdir);
end
fdir = fullfile(subfdir,'RestingState'); % path to converted nifti images and processed images
% Create a folder to store transformed Nifti files
if ~exist(fdir,'dir')
    mkdir(fdir);
else
end

% scanning parameter
TR = 2;
nvols = 250; % number of volumes acquired in each resting state run
% we use the same rs-fMRI protocol for all studies

%% Define seed and target coordinates and other variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%
SeedMNI_R = [28 38 -16]; %Right OFC
SeedName_R = 'rOFC';

SeedMNI_L = [-28 38 -16]; %left OFC
SeedName_L = 'lOFC';

TargetMNI_R = [48 38 20]; %Right LPFC
TargetName_R = 'rLPFC';

TargetMNI_L = [-48 38 20]; %Left LPFC
TargetName_L = 'lLPFC';
%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% make spherical ROIs surrounding the coordinates, constrained by the tissue probability map
% these masks are in the normalized (MNI) space

seedmask_R = fullfile(maskdir, sprintf('SeedRegion_%s.nii',SeedName_R));
seedmask_L = fullfile(maskdir, sprintf('SeedRegion_%s.nii',SeedName_L));

targetmask_R = fullfile(maskdir, sprintf('TargetRegion_%s.nii',TargetName_R));
targetmask_L = fullfile(maskdir, sprintf('TargetRegion_%s.nii',TargetName_L));


%check if spherical ROIs have already been made
if ~exist(seedmask_R,'file') && ~exist(targetmask_R,'file')  && ...
        ~exist(targetmask_L,'file')  && ~exist(targetmask_L,'file')
    % only do if all files don't exist
    
    % tpmthresh = 1; not useful 
    rget = 4; 
    %spm fmri
    
    % gray matter mask
    tpmhdr = spm_vol_nifti(fullfile(ScriptFolder,'rTPM_point1_NoCerebellum.nii')); 
        
    [tpmvol,XYZ] = spm_read_vols(tpmhdr); 
    tpmidx = find(tpmvol>0); % tpmvol only contains 0 and 255 - linear indices of 255 location
    sz = size(tpmvol);
    [MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3)); % MX,MY,MZ have the size dim and indicate x,y,z coordinates of each point
    ref_vox = round(sz/2);
    radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2); % find the distance of each point to the center
    radius_index = find(radii<rget) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3)); % sub2ind: convert subscripts to linear indices
    % why do we need this shift? before this radius_index is relative to the
    % center of the 3D, and after the center doesn't matter - radius_index
    % can be relative to any location in order to make a sphere
    
    %seed region Right
    coordidx = find(ismember(XYZ',SeedMNI_R,'rows')>0); % find the location of the seed coordinate in the 3D array - a scalar
    idx = intersect(radius_index + coordidx, tpmidx);
    newhdr = tpmhdr;
    newvol = zeros(sz);
    newhdr.fname = seedmask_R;
    newvol(idx) = 1;
    spm_write_vol(newhdr,newvol); % write an image volume to disk
    
    %seed region Left
    coordidx = find(ismember(XYZ',SeedMNI_L,'rows')>0);
    idx = intersect(radius_index + coordidx, tpmidx);
    newhdr = tpmhdr;
    newvol = zeros(sz);
    newhdr.fname = seedmask_L;
    newvol(idx) = 1;
    spm_write_vol(newhdr,newvol);
    
    %target region Right
    coordidx = find(ismember(XYZ',TargetMNI_R,'rows')>0);
    idx = intersect(radius_index + coordidx, tpmidx);
    newhdr = tpmhdr;
    newvol = zeros(sz);
    newhdr.fname = targetmask_R;
    newvol(idx) = 1;
    spm_write_vol(newhdr,newvol);
    
    % target region Left
    coordidx = find(ismember(XYZ',TargetMNI_L,'rows')>0);
    idx = intersect(radius_index + coordidx, tpmidx);
    newhdr = tpmhdr;
    newvol = zeros(sz);
    newhdr.fname = targetmask_L;
    newvol(idx) = 1;
    spm_write_vol(newhdr,newvol);
    
    clear MX MY MZ r2vox ref_vox rget sz tpmhdr tpmvol tpmidx XYZ radii radius_index
else
    fprintf('\nSeed and Target masks already made\n');
end

%% Convert DAY1 DICOMs to nifti (if it hasn't been done already)

dicompath = studyfdir;

if ~isempty(dir(fullfile(fdir,'f*.nii')))
     fprintf('\nDICOM files already converted\n');
else
     dicomnames = dir(fullfile(dicompath));
     % sdicomname = strings(length(dicomnames),1);
     
     fprintf('%d DICOM files found for DAY1\n',length(dicomnames));
     fprintf('Files converted: ');
    for i = 3:length(dicomnames)
        if ~mod(i,25)
            fprintf('DICOM to NIFTI %d%% done \n',i/length(dicomnames)*100);
        else
        end
        spm_dicom_convert(spm_dicom_headers(fullfile(dicompath,dicomnames(i).name)),'all','flat','nii',fdir);
        % sdicomname(i) = fullfile(dicompath,dicomnames(i).name); 
    end
   
    % spm_dicom_convert(spm_dicom_headers(sdicomname),'all','flat','img',maskdir);
    % note: this gives the same output (i.e. individual Nifti images
    % instead of a 4D iamge, see help spm_dicom_convert
    fprintf('\n');
    % clear dicomnames i
end


%% Do functional image preprocessing

% get anatomical image (this has been bias corrected)
aname = dir(fullfile(subadir,'m*.nii'));
afile = fullfile(subadir,aname(1).name); % the corrected T1 image

cd(fdir);

% Check if the data has already been preprocessed
if length(dir(fullfile(fdir,'s6rf*.nii'))) ~= 250
    
    % get functional images
    fnames = dir(fullfile(fdir,'f*.nii'));
    
    ffiles = [];
    for i = 1:length(fnames)
        ffiles = [ffiles; {fullfile(fdir,fnames(i).name)}];
    end
    
    spm fmri
    
    % Realignment (estimate & reslice)
    matlabbatch{1}.spm.spatial.realign.estwrite.data = {ffiles};
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % Num Passes: register to first
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1]; % default: All images + Mean image
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1; % Mask images
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    % Coregister mean realigned EPI to T1, and bring along the realigned
    % EPIs (coregister: estimate)
    matlabbatch{2}.spm.spatial.coreg.estimate.ref = {afile}; % reference image
    matlabbatch{2}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{2}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','sess', '()',{1}, '.','rfiles'));
    matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    % Normalize T1 to generate warping parameters (estimate warps & write)
    matlabbatch{3}.spm.spatial.normalise.estwrite.subj.vol = {afile};
    matlabbatch{3}.spm.spatial.normalise.estwrite.subj.resample = {afile};
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.tpm =  {fullfile(spm('dir'),'tpm','TPM.nii')};
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
    matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.prefix = 'w'; 
    
    % Smooth EPI
    matlabbatch{4}.spm.spatial.smooth.data(1) = cfg_dep('Coregister: Estimate: Coregistered Images', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{4}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{4}.spm.spatial.smooth.dtype = 0;
    matlabbatch{4}.spm.spatial.smooth.im = 0;
    matlabbatch{4}.spm.spatial.smooth.prefix = 's6'; 
    
    % make inverse deformation field
    matlabbatch{5}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(subadir,sprintf('y_%s',aname(1).name))};
    matlabbatch{5}.spm.util.defs.comp{1}.inv.space = {fullfile(fdir,sprintf('s6r%s',fnames(1).name))};
    matlabbatch{5}.spm.util.defs.out{1}.savedef.ofname = 'inverse_deformation';
    matlabbatch{5}.spm.util.defs.out{1}.savedef.savedir.saveusr = {subadir};
    
    % save batch job
    jobname = fullfile(subJobdir,'Job_Preprocessing.mat');
    save(jobname,'matlabbatch')
    
    % Run preprocessing
    spm_jobman('run', matlabbatch)
    clear matlabbatch
    
else
    fprintf('\nData already preprocessed\n');
end


%% un-norm and reslice spherical seed and target ROIs
fnames = dir(fullfile(fdir,'s6rf*.nii'));% get functional images
unnormfile = fullfile(subadir,'y_inverse_deformation.nii');% inverse deformation map

spm fmri

matlabbatch = [];
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {unnormfile};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {...
    strcat(seedmask_R,',1');...
    strcat(targetmask_R,',1');...
    strcat(seedmask_L,',1');...
    strcat(targetmask_L,',1');...
    };
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

matlabbatch{2}.spm.spatial.coreg.write.ref = {fullfile(fdir,fnames(1).name)};
matlabbatch{2}.spm.spatial.coreg.write.source(1) = ...
cfg_dep('Normalise: Write: Normalised Images (Subj 1)',... 
substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';
spm_jobman('run', matlabbatch)


% move un-norm files to the subject mask folder

movefile(fullfile(maskdir, sprintf('rwSeedRegion_%s.nii',SeedName_R)), unnormdir);
movefile(fullfile(maskdir, sprintf('wSeedRegion_%s.nii', SeedName_R)),unnormdir);

movefile(fullfile(maskdir, sprintf('rwTargetRegion_%s.nii',TargetName_R)),unnormdir);
movefile(fullfile(maskdir, sprintf('wTargetRegion_%s.nii',TargetName_R)),unnormdir);

movefile(fullfile(maskdir, sprintf('rwSeedRegion_%s.nii',SeedName_L)),unnormdir);
movefile(fullfile(maskdir, sprintf('wSeedRegion_%s.nii',SeedName_L)),unnormdir);

movefile(fullfile(maskdir, sprintf('rwTargetRegion_%s.nii',TargetName_L)),unnormdir);
movefile(fullfile(maskdir, sprintf('wTargetRegion_%s.nii',TargetName_L)),unnormdir);




%% set up GLM
[seedvol, seedXYZ] = spm_read_vols(spm_vol_nifti(fullfile(unnormdir,sprintf('rwSeedRegion_%s.nii',SeedName_R))));
R_seedidx = find(seedvol==1);

[seedvol_L, seedXYZ_L] = spm_read_vols(spm_vol_nifti(fullfile(unnormdir,sprintf('rwSeedRegion_%s.nii',SeedName_L))));
L_seedidx = find(seedvol_L==1);

% make GLM folder
modeldir = fullfile(subfdir,sprintf('GLM_%s_seed',SeedName_R));
modeldir2 = fullfile(subfdir,sprintf('GLM_%s_seed',SeedName_L));

if ~exist(modeldir,'dir')
    mkdir(modeldir);
end
if ~exist(modeldir2,'dir')
    mkdir(modeldir2);
end

% get resting state files
fnames = dir(fullfile(fdir,'s6rf*.nii'));
ffiles = [];
for i = 1:length(fnames)
    ffiles = [ffiles; {fullfile(fdir,fnames(i).name)}];
end

%make vectors of mean seed region activity
r = zeros(nvols,1);
l = zeros(nvols,1);
for i = 1:length(ffiles)
    tmpvol = spm_read_vols(spm_vol_nifti(ffiles{i}));
    r(i) = nanmean(tmpvol(R_seedidx));
    l(i) = nanmean(tmpvol(L_seedidx));
end

nname = dir(fullfile(fdir,'rp*.txt'));
nfile = fullfile(fdir, nname(1).name); % realignment parameters

%specify the models
%%%% R %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = {modeldir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = ffiles;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress.name = SeedName_R;
matlabbatch{1}.spm.stats.fmri_spec.sess.regress.val = r; %  timeseries of mean seed region activity
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {nfile};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

%estimate  model
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

%do +1 contrast (t-test) on the seed region regressor
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = SeedName_R;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{3}.spm.stats.con.delete = 0;
spm_jobman('run', matlabbatch)

%%%%%%% L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = {modeldir2};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = ffiles;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress.name = SeedName_L;
matlabbatch{1}.spm.stats.fmri_spec.sess.regress.val = l; %  timeseries of mean seed region activity
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {nfile};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

%estimate the model
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

%do a +1 contrast (t-test) on the seed region regressor
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = SeedName_L;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{3}.spm.stats.con.delete = 0;
spm_jobman('run', matlabbatch)
clear matlabbatch % need to clear it before going to the next loop


%% find coordinates of maximum connectivity in target regions
% R
[targetvol_R, targetXYZ_R] = spm_read_vols(spm_vol_nifti(fullfile(unnormdir,sprintf('rwTargetRegion_%s.nii',TargetName_R))));
targetidx_R = find(targetvol_R==1);

Rconvol = spm_read_vols(spm_vol_nifti(fullfile(modeldir,'spmT_0001.nii')));
Rconvals = Rconvol(targetidx_R);

RStimCoord = round(targetXYZ_R(:,targetidx_R(Rconvals==max(Rconvals)))');
sprintf('right coordinate: x=%01d, y=%01d, z=%01d', RStimCoord(1), RStimCoord(2), RStimCoord(3))

%L
[targetvol_L, targetXYZ_L] = spm_read_vols(spm_vol_nifti(fullfile(unnormdir,sprintf('rwTargetRegion_%s.nii',TargetName_L))));
targetidx_L = find(targetvol_L==1);

Lconvol = spm_read_vols(spm_vol_nifti(fullfile(modeldir2,'spmT_0001.nii')));
Lconvals = Lconvol(targetidx_L);

LStimCoord = round(targetXYZ_L(:,targetidx_L(Lconvals==max(Lconvals)))');
sprintf('left coordinate: x=%01d, y=%01d, z=%01d', LStimCoord(1), LStimCoord(2), LStimCoord(3))

% save coordinates
save(fullfile(subdir, 'stimulation_coordinates.mat'), 'RStimCoord', 'LStimCoord')

fileID = fopen(fullfile(subdir,'Stimulation_coordinates.txt'),'wt');
fprintf(fileID,'%s %d %d %d\n','Right',RStimCoord);
fprintf('\n');
fprintf(fileID,'%s %d %d %d\n','Left',LStimCoord);
fclose(fileID);


