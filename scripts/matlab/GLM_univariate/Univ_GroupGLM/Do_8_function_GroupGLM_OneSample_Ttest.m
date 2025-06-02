
function Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)
      
% function inputs
% modelname: a character of model name, e.g. iPE_fourpt
% connames: contrast names and index from 1st-level analysis, consistent
% with contrast index

% (you can get this directly from {SPM.xCon(:).name}
% with SPM from any subject)

% function has no output variable

% This script runs across multiple contrasts
% For each contrast specified in the first level, a pair of contrasts
% (positive & negative) is created at group level.

% This was updated to be done in external drive

%% set up path
External = '/Volumes/QF10TB/SPEEDTMS_results'; 
%studydir = '/Users/qingfangliu/Experiment';
Groupdir = fullfile(External,'UnivariateGLM','GroupRes'); % create a folder to store group analysis result
if ~exist(Groupdir,'dir')
    mkdir(Groupdir)
end

%% subject info
ResDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ExptRes';
parentDir = fileparts(ResDir);
SubInfoFile = fullfile(parentDir,'SubjectConds.XLSX');
SubInfo = readtable(SubInfoFile) ;
Subs = SubInfo.Sub;  
Excluded = SubInfo.Excluded; % whether excluding from analysis
Subs = Subs(Excluded==0); % excluding subs from analysis
nSubs = length(Subs); % number of subjects

%% Specify contrast images
ModelDir = fullfile(Groupdir,modelname); % create a model directory
Jobdir = fullfile(ModelDir, 'jobs'); % create a job folder under model directory
if ~exist(Jobdir,'dir')
    mkdir(Jobdir);
end

maskmap = fullfile(parentDir,'ROIs','SPEEDTMS_external_mask.nii'); % use explicit mask

%% Specify contrast images

for c = 1:length(connames)

    ConImages = cell(nSubs,1);
    for j = 1:nSubs
        tempConImage = fullfile(External,'UnivariateGLM',sprintf('Sub%d',Subs(j)),'fxUnivariate',modelname,sprintf('con_%04d.nii', c));
        ConImages{j} = tempConImage;
    end

    %% Define batch scripts
    clear matlabbatch
    
    pathname = fullfile(ModelDir, connames{c});
    if ~exist(pathname, 'dir')
        mkdir(pathname);
    end
         
    matlabbatch{1}.spm.stats.factorial_design.dir = {pathname}; 
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = ConImages;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {maskmap};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    % model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    % contrast
    % one sample t-test at group level, so weight is either 1 or -1
    % two T results for each contrast 
    matlabbatch{3}.spm.stats.con.spmmat(1) = {fullfile(pathname, 'SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = sprintf('positive_%s', connames{c});
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = sprintf('negative_%s', connames{c});
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 1; % delete existing contrasts
 
    %% save the job
    fname = fullfile(Jobdir, sprintf('Group_Specify_Univ_%s_%s.mat',modelname,connames{c}));
    save(fname, 'matlabbatch');

    %% run the job
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    
end % end of loop of contrasts

end

