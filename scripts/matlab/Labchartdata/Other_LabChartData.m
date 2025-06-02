
% This script is used to analyze all we need from the LabChart data:
    % Compare breathing trace across conditions;
    % Make nuisance regressors for GLM analysis;
        % (saving under NuisanceRegressor folder)
    % Determine t0 adjustment for obtaining condition onsets in GLM
        % (saving under subject folder)
    
%%
addpath('/Users/qingfangliu/Downloads/cprintf') % package for printing in color
fprintf('Working on %s Lab Chart data\n',subno)

% load the saved lab chart mat file
load(fullfile(LabChartPath,sprintf('LabChart_sub%d',subj)))

% to save processed lab chart output
LabChartdir = fullfile(studydir,subno,'ProcessedLC');

if ~exist(LabChartdir,'dir')
    mkdir(LabChartdir)
end

%% set up

AllssR = cell(nSess,nruns);
alltimecorr = zeros(nSess,nruns);
allstarts = zeros(nSess,nruns); 

if UseDay4(subj==Subs)==1 
    cd(fullfile(studydir,subno,'DAY4')) % go to folder day4
else
    cd(fullfile(studydir,subno,'DAY3')) % go to folder day3
end

dat = dir('comp*');

% load data
load(dat.name)

% whether spirometer sides were reversed
SpirometerCorrect = [SubInfo.Breathing1(Subs==subj),SubInfo.Breathing2(Subs==subj)];

%% loop over sessions and runs

for ss = 1:nSess % loop over sessions
    fprintf('Session %d/%d\n',ss,nSess)
    
for r = 1:nruns % loop over runs
    fprintf('Run %d/%d\n',r,nruns)
         
    % load labchart data for this run
    resp_data = resp_dat{ss,r};
    trig_data = trig_dat{ss,r};
    olf_data = olf_dat{ss,r};
    
    if SpirometerCorrect(ss)==0
        resp_data = (-1) * resp_data; % reverse breathing signals if needed
    end
        
    %% find time adjustment
    % determine trigger locations
    trs = find(round(trig_data./max(trig_data))==1);
    
    % Sub25-Day2-Run2: restarted scan so don't use the first few triggers before 2min
    if subj==25 && ss==1 && r==2
        trs = trs(trs>120*1000);
    end
    
    % Find the first trigger 
    tr = trs(1); 
    if subj==9 && ss==1 && r==2 % For Sub9-Day2-Run2, missed 4 triggers
        tr = tr - 4*2000;
    end
    
    tr_in_s = floor(tr/1000);
    tr_in_min = floor(tr_in_s/60);
    tr_in_s_rem = rem(tr_in_s,60);
    fprintf('First trigger location at %d min %d s\n',tr_in_min,tr_in_s_rem) % print out to visually compare with labchart data

    % determine t0 - the time point when Enter was hit
    if (subj==8 && ss==2) || (subj==4 && ss==1) % thresholding values are different for sessions with Tantor
        start = find(round(olf_data * 10)<=5 & round(olf_data * 10)>=0); % according to olf_data
        start = start(1);
    else % for all other sessions with Benjamin
        start = find(round(olf_data * 10)<=25 & round(olf_data * 10)>=20); % according to olf_data
        start = start(1);
    end

    % handle sub1-day4-run1 differently
    if subj==1 && ss==2 && r==1   % Scan missed the first four trials (included in labchart recording though)
        odor_starts = find(diff(olf_data)<-1.1); % find the odor starts of all 64 trials (arbitrary thresh -1.1 to have exactly 64 trials)
        trial_intervals = diff(odor_starts); % should be around 12-14s
        start = odor_starts(1) - 12000; % Infer "real" start based on the first odor_starts, should be negative value
    end
    
    start_in_s = floor(start/1000);
    start_in_min = floor(start_in_s/60);
    start_in_s_rem = rem(start_in_s,60);
    fprintf('Task start location: %d min %d s\n',start_in_min,start_in_s_rem) % print out to visually compare with labchart data
    
    adjust = start - tr;  % delay of hitting Enter compared with the 1st volume (in ms)
     
    % check if all adjusted times look ok (if hitting enter after 4 volumes)
    if adjust<=6000 || adjust>=8000
        cprintf('red','Adjusted time does not look correct.\n'); % print in color 
    end
    
    alltimecorr(ss,r) = adjust; % add adjust to all time corrections
    allstarts(ss,r) = start; 

    %% process breathing traces
    Use_start = tr; % start of scan
    Use_end = tr + (nTR(ss,r)-1) * TR * SampleRate; % end of scan (update to minus 1 to account for 1st vol)

    if subj==9 && ss==1 && r==2 % can only start from the 5th trigger (first 4 not recorded)
        Use_start = tr + 4*2000;
    end
    
    R = resp_data(Use_start:Use_end); % truncate the breathing data  

    % smooth breathing data
    R = smooth(R,250);

    % downsample breathing data
    R = R(1:100:end);
    sr = SampleRate/100;

    % high pass filter
    K.RT = 1/sr;
    K.row = ones(length(R),1);
    K.HParam = 50;
    R = spm_filter(K,R);

    R = zscore(R); % Zscore R before integration - important!
    IntR = cumtrapz(R); % integration to get breathing volume
    IntR = zscore(IntR); % Zscore IntR

    % down sample to match scan resolution
    ssR = R(1:20:end);
    IntssR = IntR(1:20:end);
    
    % keep flow rate, flow volume, and their squared values
    AllssR{ss,r} = [ssR, ssR.^2, IntssR, IntssR.^2];
    
    if subj==9 && ss==1 && r==2 % missing the first 4 vols, so filling with mean values
        temp_mean = mean(AllssR{ss,r},1);
        AllssR{ss,r} = [repmat(temp_mean,4,1);AllssR{ss,r}];
    end
    
    %% load and organize behavioral data
    
    if ss == 1 && UseDay4(subj==Subs)==0
       agg_rundat{ss,r} = res.reversal_learning_task_DAY2{1,r};
    elseif ss == 2 && UseDay4(subj==Subs)==0
       agg_rundat{ss,r} = res.reversal_learning_task_DAY3{1,r};
    elseif ss == 1 && UseDay4(subj==Subs)==1 
       agg_rundat{ss,r} = res.reversal_learning_task_DAY3{1,r};
    elseif ss == 2 && UseDay4(subj==Subs)==1 
       agg_rundat{ss,r} = res.reversal_learning_task_DAY4{1,r};
    end
    
    % Get cue onsets directly from data recording, with the time adjustment
    rundat = agg_rundat{ss,r}; 
    cue_onsets = start + rundat(:,11);  % col 11 is the cue onset
    time_window = -3000:4000;  % time window of breathing data signals each trial
    seg_resp_data = nan(size(rundat,1),size(time_window,2)); % nrow: # trials; ncol: # samples
        
    % Sub1-Day4-Run1: skip trial 1 of breathing trace
    if subj==1 && ss==2 && r==1
        for z = 2:length(cue_onsets)
            seg_resp_data(z,:) = resp_data((cue_onsets(z)-3000):(cue_onsets(z)+4000));
        end
    else
        for z = 1:length(cue_onsets)
            seg_resp_data(z,:) = resp_data((cue_onsets(z)-3000):(cue_onsets(z)+4000));
        end
    end
    agg_resp_data{ss,r} = seg_resp_data;
    
end % end of loop over runs

fprintf('\n')
end % end of loop over sessions

%% save all time corrections
%save(fullfile(studydir,subno,sprintf('AllTimeCorr_%s.mat',subno)),'alltimecorr');
%fprintf('Time Correction saved.\n')

%% save all time corrections
save(fullfile(studydir,subno,sprintf('AllStarts_%s.mat',subno)),'allstarts');
fprintf('Starting Time saved.\n')

%% save trial-wise respiratory data for each run
%save(fullfile(LabChartdir,sprintf('Agg_resp_data_%s.mat',subno)),'agg_resp_data');
%fprintf('agg_resp_data saved.\n')

%% Plot aggregated breathing traces for each run

max_breath = cellfun(@(x) max(nanmean(x,1)),agg_resp_data);
max_breath = max(max(max_breath))+5;
min_breath = cellfun(@(x) min(nanmean(x,1)),agg_resp_data);
min_breath = min(min(min_breath))-5;

figure
for ss = 1:nSess
    subplot(1,2,ss)
    plot(-3000:4000,nanmean(agg_resp_data{ss,1},1),'LineWidth',3)
    hold on
    plot(-3000:4000,nanmean(agg_resp_data{ss,2},1),'LineWidth',3)
    plot(-3000:4000,nanmean(agg_resp_data{ss,3},1),'LineWidth',3)
    hold off
    xlabel('Time')
    ylabel('Breathing')
    title(sprintf('Day%d',ss+1))
    legend('Run 1','Run 2','Run 3')
    ylim([min_breath,max_breath])
end
set(gcf,'position',[100,200,1000,500])
%saveas(gcf,fullfile(LabChartdir,sprintf('BreathingTraces_%s.bmp',subno)))


%% Compare breathing traces across cues and odors
all_agg_rundat = vertcat(agg_rundat{:});
all_agg_resp_data = vertcat(agg_resp_data{:});

idx_cue1 = all_agg_rundat(:,1)==1;
idx_cue2 = all_agg_rundat(:,1)==2;

idx_odor1 = all_agg_rundat(:,2)==1;
idx_odor2 = all_agg_rundat(:,2)==2;
idx_odor3 = all_agg_rundat(:,2)==3;

figure
subplot(1,2,1) % Compare breathing traces for two cues
plot(-3000:4000,nanmean(all_agg_resp_data(idx_cue1,:),1),'LineWidth',3)
hold on
plot(-3000:4000,nanmean(all_agg_resp_data(idx_cue2,:),1),'LineWidth',3)
xlabel('Time')
ylabel('Breathing')
legend('Cue 1','Cue 2')
ylim([min_breath,max_breath])

subplot(1,2,2) % Compare breathing traces for three odors
plot(-3000:4000,nanmean(all_agg_resp_data(idx_odor1,:),1),'LineWidth',3)
hold on
plot(-3000:4000,nanmean(all_agg_resp_data(idx_odor2,:),1),'LineWidth',3)
plot(-3000:4000,nanmean(all_agg_resp_data(idx_odor3,:),1),'LineWidth',3)
xlabel('Time')
ylabel('Breathing')
legend('Odor 1','Odor 2','Odor 3')
ylim([min_breath,max_breath])

set(gcf,'position',[100,200,1000,500])
%saveas(gcf,fullfile(LabChartdir,sprintf('CompareBreathing_%s.bmp',subno)))

%% Plot processed breathing signals 
figure
for ss = 1:nSess
    for r = 1:nruns
        subplot(2,3,r+(ss-1)*3)
        plot(AllssR{ss,r}(:,1)); 
        xlabel('Scans')
        ylabel('Breathing')
        title(sprintf('Day%d-Run%d',ss+1,r))
    end
end
set(gcf,'position',[100,200,1200,700])
%saveas(gcf,fullfile(LabChartdir,sprintf('BreathingFlowRate_%s.bmp',subno)))

figure
for ss = 1:nSess
    for r = 1:nruns
        subplot(2,3,r+(ss-1)*3)
        plot(AllssR{ss,r}(:,3)); 
        xlabel('Scans')
        ylabel('Breathing')
        title(sprintf('Day%d-Run%d',ss+1,r))
    end
end
set(gcf,'position',[100,200,1200,700])
%saveas(gcf,fullfile(LabChartdir,sprintf('BreathingFlowVolume_%s.bmp',subno)))

close all % close all open figures

%% save processed breathing signals
filename = fullfile(LabChartdir, sprintf('ProcessedLabChart_%s.mat',subno));
%save(filename, 'AllssR');
%fprintf('All processed LabChart mat saved - to be used as nuisance regressors.\n')
fprintf('\n')

