
clc; clear;

External = '/Volumes/QF10TB/SPEEDTMS_results';
studydir = '/Users/qingfangliu/Experiment';
ResDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ExptRes';
parentDir = fileparts(ResDir); % the entire Analysis folder
maskpath = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/Global_connectedness';



SubInfoFile = fullfile(parentDir,'SubjectConds.XLSX');
SubInfo = readtable(SubInfoFile) ;
Subs = SubInfo.Sub; % subject index
Conds = SubInfo.Cond; % TMS condition
Excluded = SubInfo.Excluded; % whether excluding from analysis

Subs = Subs(Excluded==0); % excluding subs from analysis
Conds = Conds(Excluded==0);

nSess = 2; % number of sessions
nruns = 3; % number of runs
nscans = 430; % number of scans in each run

%% load masks & find indexes
% mask should be resliced to match the dimension of scan data (i.e. normalized space)
% created using SPM's tissue probability map

gm_nii = fullfile(maskpath,'gm_0.1_3mm.nii'); % gray matter mask
gm_dat = spm_read_vols(spm_vol(gm_nii));
gm_idx = find(gm_dat > 0);
nvox = length(gm_idx); % number of voxels

wm_nii = fullfile(maskpath,'wm_0.9_3mm.nii'); % white matter mask
wm_dat = spm_read_vols(spm_vol(wm_nii));
wm_idx = find(wm_dat > 0);

csf_nii = fullfile(maskpath,'csf_0.9_3mm.nii'); % CSF mask
csf_dat = spm_read_vols(spm_vol(csf_nii));
csf_idx = find(csf_dat > 0);

%%

for subj =  Subs'
    subno = sprintf('Sub%d',subj);
    subdir = fullfile(studydir,subno);    
    subdirfunc = fullfile(External,'TempFuncCopied',subno,'func');
    sessdirsfunc = dir(fullfile(subdirfunc,'Day*'));
    
    for ss = 1:nSess
        
        sess_name = sessdirsfunc(ss).name;
        
    for r = 1:nruns
 
        fprintf('Extraxt data from %s Sess %d Run %d\n', subno, ss, r);
        dat = zeros(nscans,nvox); % initiliaze to get data
        wm_mean = zeros(nscans,1);
        csf_mean = zeros(nscans,1);
        
        % load data
        path = fullfile (subdirfunc, sess_name, sprintf('Run%d',r));
        n = dir(fullfile(path, 's6mmw3f*.nii')); % data with 3mm voxel size
               
        for i = 1:nscans
            tmp = spm_read_vols(spm_vol(fullfile(path, n(i).name)));
            dat(i,:) = tmp(gm_idx); % gray matter
            
            wm_mean(i) = mean(tmp(wm_idx)); % white matter: take the mean
            csf_mean(i) = mean(tmp(csf_idx)); % CSF: take the mean
        end
                
        % filter data
        mreg = load(fullfile(subdir, 'NuisanceRegressor', sprintf('nuisance_regressors_%s_%s_Run%d.txt', subno, sess_name ,r)));
        mreg = [zscore([mreg, mean(dat,2), wm_mean, csf_mean, [1:nscans]']), ones(nscans,1)]; % add mean(gm), mean(wm), mean(csf), drift, and constant
        b = inv(mreg'*mreg)*(mreg'*dat);
        dat = dat - mreg*b;
        
        % save data
        WhereToSave = fullfile(External, 'GlobalConn', subno, sess_name, sprintf('Run%d',r));
        if ~exist(WhereToSave,'dir')
            mkdir(WhereToSave)
        end
        
        save(fullfile(WhereToSave,'tc_filtered_3mm.mat'), 'dat');       
       
    end % end of run loop
    end % end of session loop
end

