
% This script is used for making nuisance regressors to be used for GLM
% analysis.
% The nuisance regressors include
% 24 motion parameters (6 motion parameters * 4)
% 4 processed breathing signals (4)
% 6 slice differ and slice var signals (6)
% bad volumes (if applicable)
% so the final col # is equal or above 34
% regressors from 35th indicate bad volumes

%% setup
fprintf('Make nuisance regressors for %s\n',subno)

% load preprocessed labchart breathing data
LabChartdir = fullfile(studydir,subno,'ProcessedLC');
LabChartMat = dir(fullfile(LabChartdir,'*.mat'));
load(fullfile(LabChartMat.folder,LabChartMat.name)) 

nSess = 2; % 2 sessions 
nruns = 3; % 3 runs each session
sessdirs = dir(fullfile(studydir,subno,'func','Day*'));
NROutdir = fullfile(studydir, subno, 'NuisanceRegressor');
if ~exist(NROutdir,'dir') % create a NR folder
    mkdir(NROutdir)
end

TR = 2; % TR of the scan in sec
nTR = zeros(nSess,nruns);  % number of volumes
NR = cell(nSess,nruns);  %NR Nuisanace Regressors
slicediff = cell(nSess,nruns);
slicevar = cell(nSess,nruns);

%% loop over sessions and runs
for ss = 1:nSess % sessions
    sess_name = sessdirs(ss).name; % current session name
    
for r = 1:nruns % runs
    
    run_name = sprintf('Run%d',r);
    fprintf('Working on %s %s\n',sess_name,run_name)
    
    % find the motion parameters
    funcRunDir = fullfile(studydir,subno,'func',sess_name,run_name);
    mpfilename = dir(fullfile(funcRunDir, 'rp_*.txt'));    %rp file contains 6 columns represnting positions  
    mp = load(fullfile(funcRunDir, mpfilename(1).name));          % mp motion parameters
    nTR(ss,r) = size(mp,1); % number of TR

    % add diff, squared mp, and squared diff
    NR{ss,r} = [mp, [zeros(1,6); diff(mp)], mp.^2, [zeros(1,6); diff(mp).^2]];

    figure;
    subplot(2,4,1); plot(NR{ss,r}(:,1:3)); xlabel('scans'); ylabel('mm'); 
    subplot(2,4,2); plot(NR{ss,r}(:,4:6)); xlabel('scans'); ylabel('deg');
    subplot(2,4,3); plot(NR{ss,r}(:,7:9)); xlabel('scans'); ylabel('diff(mm)'); 
    subplot(2,4,4); plot(NR{ss,r}(:,10:12)); xlabel('scans'); ylabel('diff(deg)');
    subplot(2,4,5); plot(NR{ss,r}(:,13:15)); xlabel('scans'); ylabel('mm^2');
    subplot(2,4,6); plot(NR{ss,r}(:,16:18)); xlabel('scans'); ylabel('deg^2');
    subplot(2,4,7); plot(NR{ss,r}(:,19:21)); xlabel('scans'); ylabel('diff(mm)^2'); 
    subplot(2,4,8); plot(NR{ss,r}(:,22:24)); xlabel('scans'); ylabel('diff(deg)^2');
    set(gcf, 'PaperPosition', [2 1 18 8]);
    tmpfigname = fullfile(NROutdir,sprintf('MotionParameters_%s_%s_%s.bmp',subno,sess_name,run_name));
    saveas(gcf,tmpfigname)
    SubResDir = fullfile(ResDir,subno); % copy mp figures to dropbox
    copyfile(tmpfigname,SubResDir)

    % add four breathing traces to nuisance regressor
    NR{ss,r} = [AllssR{ss,r}, NR{ss,r}];

    %% analyze slice differences and variances
    idx1 = 1:2:58;   % odd slices
    idx2 = 2:2:58;   % even slices

    funcfilenames = dir(fullfile(funcRunDir,'f*.nii')); % using images before preprocessing

    for vol = 1:length(funcfilenames)
        v = spm_read_vols(spm_vol(fullfile(funcRunDir,funcfilenames(vol).name)));   %dimensions of v show that there are 58 slices
        odd_slice = nanmean(nanmean(nanmean(v(:,:,idx1),3)));
        even_slice = nanmean(nanmean(nanmean(v(:,:,idx2),3)));
        slicediff{ss, r}(vol) = (odd_slice - even_slice); % diff b/t odd and even slices
        slicevar{ss, r}(vol) = var(nanmean(squeeze(nanmean(v)))); % var across all slices
    end

           
end % end of run
end % end of session


%% get mean and std across all runs and all sessions

all_slicediff = [slicediff{:}]; % aggregate all slice differences across sessions and runs
all_slicevar = [slicevar{:}]; % aggregate all slice variances across sessions and runs
aztmp = [max(abs(zscore(all_slicediff))),max(abs(zscore(all_slicevar)))];

adj_slicediff = cell(nSess,nruns); % to store standardized slice differences
adj_slicevar = cell(nSess,nruns); % to store standardlized slice variances
badvol = cell(nSess,nruns); % to store bad volumes in each run of each session
ctr_badvol = zeros(nSess,nruns); % to store number of bad volumes in each run of each session
ctr_NR = zeros(nSess,nruns); % to store number of nuisance regressors

% set criterion for above mean
crit = [4,4];        % criterion for slicediff and slicevar, eyeball it looking at the peaks

for ss = 1:nSess
    sess_name = sessdirs(ss).name;
    figure
    for r = 1:nruns
        run_name = sprintf('Run%d',r);
        
        adj_slicediff{ss, r} = abs((slicediff{ss, r} - mean(all_slicediff))./std(all_slicediff));
        adj_slicevar{ss, r} = abs((slicevar{ss, r} - mean(all_slicevar))./std(all_slicevar));

        sh = subplot(nruns, 2, r*2-1);    
        plot(adj_slicediff{ss, r})
        set(sh, 'YLim', [0,max(aztmp)+1]);
        title(sprintf('slice diff, %s %s',sess_name,run_name));
        hold on
        plot(find(adj_slicediff{ss, r}> crit(1)), adj_slicediff{ss,r}(adj_slicediff{ss, r}> crit(1)), '*r');
        
        sh = subplot(nruns, 2, r*2);
        plot(adj_slicevar{ss, r})
        set(sh, 'YLim', [0,max(aztmp)+1]);
        title(sprintf('slice var, %s %s',sess_name,run_name));
        hold on
        plot(find(adj_slicevar{ss, r}> crit(2)), adj_slicevar{ss,r}(adj_slicevar{ss, r}> crit(2)), '*r');

        % include regressors for bad volumes
        tmpbadvol = unique([find(adj_slicediff{ss, r}> crit(1)), find(adj_slicevar{ss, r}> crit(2))]);
        badvol{ss,r} = tmpbadvol;
        ctr_badvol(ss,r) = length(tmpbadvol);         % count bad volumes
        bv = zeros(nTR(ss,r),ctr_badvol(ss,r));
        for i = 1:ctr_badvol(ss,r)
            bv(tmpbadvol(i),i) = 1;
        end
            % the way it works is that if no vol is bad then bv is an empty
            % matrix - will have no effect on NR; if anyone vol is bad then
            % bv is not empty and will be added to NR
            % After finding bad volumes of each run, a nuisance regressor
            % will be created with 1 in that volume and 0 everywhere else.
        
        % add slicediff and diff of slicediff and bad volume regressors
        NR{ss,r} = [NR{ss,r},adj_slicediff{ss, r}', adj_slicevar{ss, r}',...
            [0; diff(adj_slicediff{ss, r})'], [0; diff(adj_slicevar{ss, r})'],...
            [0; diff(adj_slicediff{ss, r})'].^2, [0; diff(adj_slicevar{ss, r})'].^2,bv];   
        ctr_NR(ss,r) = size(NR{ss,r},2);
        % save nuisance regressors as txt
        NRfilename = fullfile(NROutdir,sprintf('nuisance_regressors_%s_%s_%s.txt', subno, sess_name, run_name));
        dlmwrite(NRfilename, zscore(NR{ss,r}), 'delimiter', ' ', 'precision', '%.24f');   %z scoring
    end % end of runs
    
    % Save Slice diff and var figures and copy them to dropbox
    tmpfigname = fullfile(NROutdir,sprintf('SliceDiffVar_%s_%s.bmp',subno,sess_name));
    saveas(gcf,tmpfigname)
    SubResDir = fullfile(ResDir,subno); 
    copyfile(tmpfigname,SubResDir)
    
end % end of sessions

% finally save bad vols of each run and each session
filename = fullfile(NROutdir,sprintf('BadVolumes_%s.mat',subno));
save(filename,'badvol','ctr_badvol','ctr_NR')
copyfile(filename,SubResDir)

%close all % close all figures

