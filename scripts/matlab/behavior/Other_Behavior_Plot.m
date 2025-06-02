
% clear matlab
clear; clc; close all

ResDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ExptRes';
ScriptDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ForExperiment';

% load variables from behavioral analysis
load(fullfile(ResDir,'BehAllSubjects.mat'))

AvgIDAcc = [mean(AggID_acc.A,2),mean(AggID_acc.P,2)];
horzcat(Subs,AvgIDAcc)

% Identify three subjects with lower than 0.8 ID task accuracy on one of
% two sessions

%% Exclude subjects
% Subjects to exclude (only run this part once!)
% Find index location of subjects to be excluded
BadSubIdx = ismember(Subs,[5,9,35,42]); % which subjects to exclude

AggRevAcc.P = AggRevAcc.P(:,:,~BadSubIdx);
AggRevAcc.A = AggRevAcc.A(:,:,~BadSubIdx);

AggRevRT.P = AggRevRT.P(:,:,~BadSubIdx);
AggRevRT.A = AggRevRT.A(:,:,~BadSubIdx);

AggAcc = AggAcc(~BadSubIdx);
AggmeanRT = AggmeanRT(~BadSubIdx);
AggRTs = AggRTs(~BadSubIdx);

AggID_acc.P = AggID_acc.P(~BadSubIdx,:);
AggID_acc.A = AggID_acc.A(~BadSubIdx,:);
AggID_RT.P = AggID_RT.P(~BadSubIdx,:);
AggID_RT.A = AggID_RT.A(~BadSubIdx,:);

Subs = Subs(~BadSubIdx,:);
nsubs = length(Subs);

% save updated variables
save(fullfile(ResDir,'BehAllSubjects_cleaned.mat'),'AggAcc','AggmeanRT','AggRTs','AggRevAcc','AggRevRT','AggID_acc','AggID_RT','Subs')

%% Plot summary of odor prediction task

SummaryOdorPredRev = arrayfun(@(x) mean(x.Rev,2),AggAcc,'UniformOutput',false);
SummaryOdorPredRev = cell2mat(SummaryOdorPredRev);

SummaryOdorPredNonrev = arrayfun(@(x) mean(x.Nonrev,2),AggAcc,'UniformOutput',false);
SummaryOdorPredNonrev = cell2mat(SummaryOdorPredNonrev);

% to use violin plots
addpath('/Users/qingfangliu/Downloads/Violinplot-Matlab-master')

%%
figure
subplot(2,2,1)
violinplot(SummaryOdorPredRev');
xticklabels({'P','A'})
title('Reversal trials')
ylabel('Accuracy')

subplot(2,2,2)
violinplot(SummaryOdorPredNonrev');
xticklabels({'P','A'})
title('Nonreversal trials')
ylabel('Accuracy')

subplot(2,2,3)
violinplot(diff(SummaryOdorPredRev)');
title('Reversal trials')
ylabel('Accuracy diff (A-P)')

subplot(2,2,4)
violinplot(diff(SummaryOdorPredNonrev)');
title('Nonreversal trials')
ylabel('Accuracy diff (A-P)')

set(gcf,'position',[100,100,800,800])
saveas(gcf,fullfile(ResDir,'OdorPred_AccRevNonrev.bmp'))

%% Plot summary of odor prediction task for each run and each session
OdorPredNonrev.P = arrayfun(@(x) x.Nonrev(1,:),AggAcc,'UniformOutput',false);
OdorPredNonrev.P = cell2mat(OdorPredNonrev.P(:));
OdorPredNonrev.A = arrayfun(@(x) x.Nonrev(2,:),AggAcc,'UniformOutput',false);
OdorPredNonrev.A = cell2mat(OdorPredNonrev.A(:));

%% Compare nonrev trial acc across runs
figure
subplot(1,3,1)
violinplot(OdorPredNonrev.P);
ylabel('Accuracy')
title('Placebo')
xlabel('Runs')
ylim([0.2,1])

subplot(1,3,2)
violinplot(OdorPredNonrev.A);
ylabel('Accuracy')
title('Active')
xlabel('Runs')
ylim([0.2,1])

subplot(1,3,3)
violinplot(OdorPredNonrev.A - OdorPredNonrev.P);
ylabel('Accuracy Diff (A-P)')
title('Active - Placebo')
xlabel('Runs')
ylim([-0.3,0.3])

set(gcf,'position',[100,100,1200,400])
saveas(gcf,fullfile(ResDir,'OdorPred_AccRevCompRuns.bmp'))


%% Plot RT comparison

% for all trials (both correct and wrong)
OdorPredRT.P = arrayfun(@(x) x.All(1,:),AggmeanRT,'UniformOutput',false);
OdorPredRT.P = cell2mat(OdorPredRT.P(:));
OdorPredRT.A = arrayfun(@(x) x.All(2,:),AggmeanRT,'UniformOutput',false);
OdorPredRT.A = cell2mat(OdorPredRT.A(:));

RTmax = max([OdorPredRT.P;OdorPredRT.A],[],'all');
RTmin = min([OdorPredRT.P;OdorPredRT.A],[],'all');

figure
subplot(1,3,1)
violinplot(OdorPredRT.P);
ylim([RTmin,RTmax])
xlabel('Runs')
title('Placebo')

subplot(1,3,2)
violinplot(OdorPredRT.A);
ylim([RTmin,RTmax])
xlabel('Runs')
title('Active')

subplot(1,3,3)
violinplot(OdorPredRT.A - OdorPredRT.P);
xlabel('Runs')
title('Active - Placebo')

set(gcf,'position',[100,100,1200,400])
saveas(gcf,fullfile(ResDir,'OdorPred_RTCompRuns.bmp'))

% for correct trials only
OdorPredRTC.P = arrayfun(@(x) x.Corr(1,:),AggmeanRT,'UniformOutput',false);
OdorPredRTC.P = cell2mat(OdorPredRTC.P(:));
OdorPredRTC.A = arrayfun(@(x) x.Corr(2,:),AggmeanRT,'UniformOutput',false);
OdorPredRTC.A = cell2mat(OdorPredRTC.A(:));

RTCmax = max([OdorPredRTC.P;OdorPredRTC.A],[],'all');
RTCmin = min([OdorPredRTC.P;OdorPredRTC.A],[],'all');

figure
subplot(1,3,1)
violinplot(OdorPredRTC.P);
ylim([RTCmin,RTCmax])
xlabel('Runs')
title('Placebo')

subplot(1,3,2)
violinplot(OdorPredRTC.A);
ylim([RTCmin,RTCmax])
xlabel('Runs')
title('Active')

subplot(1,3,3)
violinplot(OdorPredRTC.A - OdorPredRTC.P);
xlabel('Runs')
title('Active - Placebo')

set(gcf,'position',[100,100,1200,400])
saveas(gcf,fullfile(ResDir,'OdorPred_RTCompRuns_CorrOnly.bmp'))


%% Plot RT distributions

for i = 1:length(AggRTs)
    tmpRTs = AggRTs{i};
    tmpRTsP{i} = vertcat(tmpRTs{1,:});
    tmpRTsA{i} = vertcat(tmpRTs{2,:});
    
    tmpRTsPR1{i} = tmpRTs{1,1};
    tmpRTsPR2{i} = tmpRTs{1,2};
    tmpRTsPR3{i} = tmpRTs{1,3};
    tmpRTsAR1{i} = tmpRTs{2,1};
    tmpRTsAR2{i} = tmpRTs{2,2};
    tmpRTsAR3{i} = tmpRTs{2,3};
end

AggRTsP = cat(1,tmpRTsP{:});
AggRTsA = cat(1,tmpRTsA{:});

figure
subplot(2,2,1)
histogram(AggRTsP)
hold on
histogram(AggRTsA)
hold off
legend('Placebo','Active')

AggRTsPR1 = cat(1,tmpRTsPR1{:});
AggRTsPR2 = cat(1,tmpRTsPR2{:});
AggRTsPR3 = cat(1,tmpRTsPR3{:});
AggRTsAR1 = cat(1,tmpRTsAR1{:});
AggRTsAR2 = cat(1,tmpRTsAR2{:});
AggRTsAR3 = cat(1,tmpRTsAR3{:});

subplot(2,2,2)
histogram(AggRTsPR1)
hold on
histogram(AggRTsAR1)
hold off
legend('Placebo','Active')
title('Run1')

subplot(2,2,3)
histogram(AggRTsPR2)
hold on
histogram(AggRTsAR2)
hold off
legend('Placebo','Active')
title('Run2')

subplot(2,2,4)
histogram(AggRTsPR3)
hold on
histogram(AggRTsAR3)
hold off
legend('Placebo','Active')
title('Run3')

set(gcf,'position',[100,100,800,600])
saveas(gcf,fullfile(ResDir,'OdorPred_RTDist_CorrOnly.bmp'))


%% Plot odor prediction task performance using 4-trial window

% averaged accuracy across all subjects
figure
subplot(1,2,1)
plot(-1:2,squeeze(mean(AggRevAcc.P,[2,3])),'-go','LineWidth',2)
hold on
plot(-1:2,squeeze(mean(AggRevAcc.A,[2,3])),'-ro','LineWidth',2)
hold off
xlim([-1 2])
ylim([0 1])
xlabel('Trials')
ylabel('Accuracy')
title('All three runs')
legend('Placebo','Active','Location','southeast')

% averaged accuracy for each subject
PlotSubRevAccP = squeeze(mean(AggRevAcc.P,2));
PlotSubRevAccA = squeeze(mean(AggRevAcc.A,2));

jitter = 0.2;
subplot(1,2,2)
scatter(ones(1,nsubs)*(-1),PlotSubRevAccP(1,:),'go','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(zeros(1,nsubs),PlotSubRevAccP(2,:),'go','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(ones(1,nsubs),PlotSubRevAccP(3,:),'go','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(ones(1,nsubs)*2,PlotSubRevAccP(4,:),'go','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(ones(1,nsubs)*(-1),PlotSubRevAccA(1,:),'ro','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(zeros(1,nsubs),PlotSubRevAccA(2,:),'ro','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(ones(1,nsubs),PlotSubRevAccA(3,:),'ro','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(ones(1,nsubs)*2,PlotSubRevAccA(4,:),'ro','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold off
xlim([-1 2])
ylim([0 1])
xlabel('Trials')
ylabel('Accuracy')
title('All three runs')

set(gcf,'position',[100,100,800,400])
saveas(gcf,fullfile(ResDir,'OdorPred_Acc.bmp'))


%% Plot acc 1st trial after reversal
figure
subplot(2,3,1)
violinplot(PlotSubRevAccP(3,:));
ylim([0.4 1])
title('Placebo TMS')
ylabel('1st trial')

subplot(2,3,2)
violinplot(PlotSubRevAccA(3,:));
ylim([0.4 1])
title('Active TMS')

subplot(2,3,3)
violinplot(PlotSubRevAccA(3,:)-PlotSubRevAccP(3,:));
ylim([-0.4 0.4])
title('Active - Placebo')

% Plot acc 2nd trial after reversal
subplot(2,3,4)
violinplot(PlotSubRevAccP(4,:));
ylim([0.4 1])
title('Placebo TMS')
ylabel('2nd trial')

subplot(2,3,5)
violinplot(PlotSubRevAccA(4,:));
ylim([0.4 1])
title('Active TMS')

subplot(2,3,6)
violinplot(PlotSubRevAccA(4,:)-PlotSubRevAccP(4,:));
ylim([-0.4 0.4])
title('Active - Placebo')

set(gcf,'position',[100,100,800,600])
saveas(gcf,fullfile(ResDir,'OdorPred_1st2ndPostRev.bmp'))


%% plot acc with RT

RTCmax = max([OdorPredRTC.P;OdorPredRTC.A],[],'all');
RTCmin = min([OdorPredRTC.P;OdorPredRTC.A],[],'all');

figure
subplot(1,2,1)
plot(OdorPredNonrev.A,OdorPredRTC.A,'.','MarkerSize',12)
xlabel('Accuracy (Nonrev)')
ylabel('Mean RT (ms)')
title('Active TMS')
ylim([RTCmin,RTCmax])

subplot(1,2,2)
plot(OdorPredNonrev.P,OdorPredRTC.P,'.','MarkerSize',12)
xlabel('Accuracy (Nonrev)')
ylabel('Mean RT (ms)')
title('Placebo TMS')
ylim([RTCmin,RTCmax])
legend('Run1','Run2','Run3','Location','southwest')

set(gcf,'position',[100,100,600,300])
saveas(gcf,fullfile(ResDir,'OdorPred_AccRT.bmp'))

corr(OdorPredNonrev.A(:),OdorPredRTC.A(:))
corr(OdorPredNonrev.P(:),OdorPredRTC.P(:))

corr(OdorPredNonrev.A(:,1),OdorPredRTC.A(:,1))
corr(OdorPredNonrev.A(:,2),OdorPredRTC.A(:,2))
corr(OdorPredNonrev.A(:,3),OdorPredRTC.A(:,3))


%% Plot odor prediction task performance using 4-trial window (RT)

% averaged accuracy for each subject
PlotSubRevRTP = squeeze(nanmean(AggRevRT.P,2));
PlotSubRevRTA = squeeze(nanmean(AggRevRT.A,2));

figure
jitter = 0.2;
subplot(2,3,1)
scatter(ones(1,nsubs)*(-1),PlotSubRevRTP(1,:),'go','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(zeros(1,nsubs),PlotSubRevRTP(2,:),'go','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(ones(1,nsubs),PlotSubRevRTP(3,:),'go','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(ones(1,nsubs)*2,PlotSubRevRTP(4,:),'go','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(ones(1,nsubs)*(-1),PlotSubRevRTA(1,:),'ro','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(zeros(1,nsubs),PlotSubRevRTA(2,:),'ro','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(ones(1,nsubs),PlotSubRevRTA(3,:),'ro','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
scatter(ones(1,nsubs)*2,PlotSubRevRTA(4,:),'ro','LineWidth',2,'jitter', 'on', 'jitterAmount', jitter)
hold on
plot(-1:2,squeeze(median(PlotSubRevRTP,2,'omitnan')),'-go','LineWidth',2)
hold on
plot(-1:2,squeeze(median(PlotSubRevRTA,2,'omitnan')),'-ro','LineWidth',2)
hold off

xlim([-1 2])
ylim([500 1500])
xlabel('Trials')
ylabel('RT (ms)')
title('All three runs')


% averaged RT across all subjects
subplot(2,3,2)
plot(-1:2,squeeze(nanmean(PlotSubRevRTP,2)),'-go','LineWidth',2)
hold on
plot(-1:2,squeeze(nanmean(PlotSubRevRTA,2)),'-ro','LineWidth',2)
hold off
xlim([-1 2])
ylim([800 1100])
xlabel('Trials')
ylabel('Mean RT (ms)')
title('All three runs')
legend('Placebo','Active','Location','northeast')


% plot RT diff (A-P) at each of the four point window 

PlotSubRevRT_diff = PlotSubRevRTA - PlotSubRevRTP;
yPlotSubRevRT_diff = reshape(PlotSubRevRT_diff,[],1);
xxWindow = reshape(repmat(-1:2,1,nsubs),[],1);

subplot(2,3,4)
violinplot(yPlotSubRevRT_diff,xxWindow); % both inputs are vectors of the same length (2nd being group variable)
title('Active - Placebo')
ylabel('RT difference (ms)')
xlabel('Trials')

yPlotSubRevRTA = reshape(PlotSubRevRTA,[],1);
yPlotSubRevRTP = reshape(PlotSubRevRTP,[],1);

subplot(2,3,3)
violinplot(yPlotSubRevRTA,xxWindow,'ViolinColor',[1,0,0]);
hold on
violinplot(yPlotSubRevRTP,xxWindow,'ViolinColor',[0,1,0]);
ylabel('RT (ms)')
xlabel('Trials')

subplot(2,3,5)
plot(-1:2,PlotSubRevRT_diff)
hold on
plot(-1:2,mean(PlotSubRevRT_diff,2),'LineWidth',5,'Color','r')
hold off
title('Active - Placebo')
ylabel('RT difference (ms)')

subplot(2,3,6)
plot([0,1],[zeros(1,nsubs);diff(PlotSubRevRT_diff([2,3],:),1)]);
title('Active - Placebo')
ylabel('Shifted RT difference (ms)')

set(gcf,'position',[100,100,900,600])
saveas(gcf,fullfile(ResDir,'OdorPred_RT4TrialWin.bmp'))


%% repeated-anova 

% test 2*4
addpath('/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/MatlabFuns')
datSubRevRT = [PlotSubRevRTP',PlotSubRevRTA'];
repanova(datSubRevRT,[2,4]);

% test 2*2 (on and after reversal)
datSubRevRT = [PlotSubRevRTP',PlotSubRevRTA'];
datSubRevRT_revOnly = datSubRevRT(:,[2,3,6,7]);
repanova(datSubRevRT_revOnly,[2,2]);



%% plot acc diff b/t 1st trial and reversal trial

Rev1P = squeeze(AggRevAcc.P(3,:,:));
Rev0P = squeeze(AggRevAcc.P(2,:,:));
MeanRev1P = mean(Rev1P,1);
MeanRev0P = mean(Rev0P,1);

Rev1A = squeeze(AggRevAcc.A(3,:,:));
Rev0A = squeeze(AggRevAcc.A(2,:,:));
MeanRev1A = mean(Rev1A,1);
MeanRev0A = mean(Rev0A,1);

figure
subplot(2,2,1)
scatter(MeanRev0P,MeanRev1P,'filled')
corr(MeanRev0P',MeanRev1P')
xlabel('Accuracy on reversal trial (P)')
ylabel('Accuracy on 1st trial after reversal (P)')
subplot(2,2,2)
scatter(MeanRev0A,MeanRev1A,'filled')
corr(MeanRev0A',MeanRev1A')
xlabel('Accuracy on reversal trial (A)')
ylabel('Accuracy on 1st trial after reversal (A)')

DiffRevP = MeanRev1P - MeanRev0P;
DiffRevA = MeanRev1A - MeanRev0A;
subplot(2,2,3)
scatter(DiffRevP,DiffRevA,'filled')
xlabel('Accuracy change (P)')
ylabel('Accuracy change (A)')
xlim([0.3,1])
ylim([0.3,1])
h = refline(1,0);
h.LineStyle = ':';

subplot(2,2,4)
violinplot(DiffRevA - DiffRevP);
title('Accuracy change (A - P)')

set(gcf,'position',[100,100,700,700])
saveas(gcf,fullfile(ResDir,'AccChange.bmp'))


%% plot ID task performance 

jitter = 0.2;

figure
subplot(2,2,1)
bar(mean(AggID_acc.A,1),'FaceAlpha',.3,'FaceColor','r')
hold on
scatter(ones(1,nsubs),AggID_acc.A(:,1), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nsubs)*2,AggID_acc.A(:,2), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nsubs)*3,AggID_acc.A(:,3), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
title('Active TMS')
ylabel('Odor identification accuracy')
xlabel('Runs')

subplot(2,2,2)
bar(mean(AggID_acc.P,1),'FaceAlpha',.3,'FaceColor','g')
hold on
scatter(ones(1,nsubs),AggID_acc.P(:,1), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nsubs)*2,AggID_acc.P(:,2), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nsubs)*3,AggID_acc.P(:,3), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
title('Placebo TMS')
xlabel('Runs')

RTmax = max(max(max(AggID_RT.A),max(AggID_RT.P)));
subplot(2,2,3)
bar(mean(AggID_RT.A,1),'FaceAlpha',.3,'FaceColor','r')
hold on
scatter(ones(1,nsubs),AggID_RT.A(:,1), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nsubs)*2,AggID_RT.A(:,2), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nsubs)*3,AggID_RT.A(:,3), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
ylabel('Odor identification RT (ms)')
xlabel('Runs')
ylim([0,RTmax])

subplot(2,2,4)
bar(mean(AggID_RT.P,1),'FaceAlpha',.3,'FaceColor','g')
hold on
scatter(ones(1,nsubs),AggID_RT.P(:,1), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nsubs)*2,AggID_RT.P(:,2), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nsubs)*3,AggID_RT.P(:,3), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
xlabel('Runs')
ylim([0,RTmax])

set(gcf,'position',[100,100,700,700])
saveas(gcf,fullfile(ResDir,'ID_AccRT.bmp'))

%% close after printing
close all
cd(ScriptDir)

