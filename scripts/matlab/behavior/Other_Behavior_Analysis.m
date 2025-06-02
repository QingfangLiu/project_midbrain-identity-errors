
% This script is to analyze all behavioral choice and RT data, of odor
% prediction task, and odor identification task;
% and plot individual behavioral results

% clear matlab
clear; clc; close all

% set up path
ResDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ExptRes';
ScriptDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ForExperiment';
studydir = '/Users/qingfangliu/Experiment';

SubCtr = 0;
RunCtr = 0;
nruns = 3;
nSess = 2;

parentDir = fileparts(ResDir);
SubInfoFile = fullfile(parentDir,'SubjectConds.XLSX');
SubInfo = readtable(SubInfoFile) ;
Subs = SubInfo.Sub; % subject index
Conds = SubInfo.Cond; % TMS condition
OdorDay1 = SubInfo.OdorDay1; % whether use day1 odor ratings
UseDay4 = SubInfo.UseDay4; % if use data from day3 and day4

%% loop over subject
for subj = Subs'
    fprintf('Subject %d\n',subj)
    SubCtr = SubCtr + 1;
    subfolder = sprintf('Sub%d',subj);
    SubResDir = fullfile(ResDir,subfolder);
    if ~exist(SubResDir,'dir')
        mkdir(SubResDir)
    end
    
    if UseDay4(SubCtr)==1     
        cd(fullfile(studydir,subfolder,'DAY4')) % only load from day4 for subj1 and subj35
    else
        cd(fullfile(studydir,subfolder,'DAY3')) % load data from day3 (should contain all data from all sessions)
    end
    dat = dir('comp*');
    
    % load data
    load(dat.name)
    SubRevAcc = cell(nSess,nruns);
    SubRevRT = cell(nSess,nruns);
    IDacc = cell(nSess,nruns);
    IDRT = cell(nSess,nruns);
    
%%    
for ss = 1:nSess % two sessions
    
    % code if this sess is active or placebo
    TMS = 1; % TMS: 1 for placebo, 2 for active
    if strcmp(Conds{SubCtr},'PA') && ss==2 
        TMS = 2;
    end
    if strcmp(Conds{SubCtr},'AP') && ss==1
        TMS = 2;
    end
    
%%
for r = 1:nruns % 3 runs
    
    % load data
    fprintf('Session %d Run%d\n',ss,r)
    RunCtr = RunCtr + 1;
    if ss == 1 && UseDay4(SubCtr)==0
       rundat = res.reversal_learning_task_DAY2{1,r};
       IDdat = res.identification_task_DAY2{1,r};
    elseif ss == 2 && UseDay4(SubCtr)==0 
       rundat = res.reversal_learning_task_DAY3{1,r};
       IDdat = res.identification_task_DAY3{1,r};
    elseif ss == 1 && UseDay4(SubCtr)==1  
        rundat = res.reversal_learning_task_DAY3{1,r};
        IDdat = res.identification_task_DAY3{1,r};
    elseif ss == 2 && UseDay4(SubCtr)==1  
        rundat = res.reversal_learning_task_DAY4{1,r};
        IDdat = res.identification_task_DAY4{1,r};
    end
    
    rundat(:,14) = double(rundat(:,2) == rundat(:,13)); % add correct or wrong 
    Reversal = rundat(:,3);
    Correct = rundat(:,14); % note: non-resp is coded as 0
    a = sum(Reversal==1 & Correct==1);
    b = sum(Reversal==0 & Correct==1);
    c = sum(Reversal==1 & Correct==0);
    d = sum(Reversal==0 & Correct==0);
    tmpRTA = nanmean(rundat(:,6)); % mean RT of all trials, excluding nan 
    tmpRTC = nanmean(rundat(Correct==1,6)); % mean RT of all correct trials, excluding nan
    tmpRTs = rundat(Correct==1,6); % all RTs of all correct trials
    
    % save current rundat into spreadsheet
    ntrials = size(rundat,1);
    datmat{RunCtr} = [repmat(subj,ntrials,1),repmat(ss,ntrials,1),repmat(TMS,ntrials,1),...
        repmat(r,ntrials,1),(1:ntrials)',rundat];
   
    % separating the data matrix into two cues
    for cueID = 1:2 % loop cue 1 and 2
        locCue = rundat(:,1) == cueID; % trial location of a certain cue
        cuedat = rundat(locCue,:); % data with this cue
        reversal_loc = find(cuedat(:,3) == 1)';
        for k = 1:6 % 6 times of reversals for each cue
            rows = (reversal_loc(k)-1):(reversal_loc(k)+2);
            table = cuedat(rows,[1:3,13,14,6]); % cue,odor,rev,resp,acc,RT
            Resps{cueID,k} = table(:,5);
            RTs{cueID,k} = table(:,6);
            %tmpcodedAcc = table(:,5);
            %tmpcodedAcc(2) = 1 - tmpcodedAcc(2); % revert acc for the reversal location
            %RTs{cueID,k} = table(:,6) .* tmpcodedAcc; % code 'wrong' trial RT as zero
            %RTs{cueID,k}(RTs{cueID,k}==0) = nan; % code 'wrong' trial RT as nan
        end
    end

    SubRevAcc{TMS,r} = cell2mat(reshape(Resps,1,[])); % keep all sessions and all runs 
    SubRevRT{TMS,r} = cell2mat(reshape(RTs,1,[])); % keep all sessions and all runs 
    SubRTs{TMS,r} = tmpRTs;
    
    Acc.Rev(TMS,r) = a/(a+c); % accuracy for reversal trials
    Acc.Nonrev(TMS,r) = b/(b+d); % accuracy for nonreversal trials
    meanRT.All(TMS,r) = tmpRTA; % mean RT across all trials
    meanRT.Corr(TMS,r) = tmpRTC; % mean RT across all correct trials
    
    IDacc{TMS,r} = IDdat(:,1)==IDdat(:,4);
    IDRT{TMS,r} = IDdat(:,5);
    % placebo is on row 1 and active is on row 2
    
end % end of runs
end % end of sessions

AllRevAcc.P = horzcat(SubRevAcc{1,:}); % placebo acc (combining runs)
AllRevAcc.A = horzcat(SubRevAcc{2,:}); % active acc (combining runs)
AllRevRT.P = horzcat(SubRevRT{1,:}); % placebo RT (combining runs)
AllRevRT.A = horzcat(SubRevRT{2,:}); % active RT (combining runs)

% save result of current subject
AggRevAcc.P(:,:,SubCtr) = AllRevAcc.P;
AggRevAcc.A(:,:,SubCtr) = AllRevAcc.A;
AggRevRT.P(:,:,SubCtr) = AllRevRT.P;
AggRevRT.A(:,:,SubCtr) = AllRevRT.A;
AggRTs{SubCtr} = SubRTs;

AggAcc(SubCtr) = Acc;
AggmeanRT(SubCtr) = meanRT;

AllID.acc.P = horzcat(IDacc{1,:}); % placebo ID acc
AllID.acc.A = horzcat(IDacc{2,:}); % active ID acc
AllID.RT.P = horzcat(IDRT{1,:}); % placebo ID RT
AllID.RT.A = horzcat(IDRT{2,:}); % active ID RT

% save result of current subject
AggID_acc.P(SubCtr,:,:) = AllID.acc.P;
AggID_acc.A(SubCtr,:,:) = AllID.acc.A;
AggID_RT.P(SubCtr,:,:) = AllID.RT.P;
AggID_RT.A(SubCtr,:,:) = AllID.RT.A;

%% Plot accuracy
% use green for placebo, use red for active
% figure
% % plot for each session and each run
% for ss = 1:nSess
%     for r = 1:nruns
%         subplot(2,6,r+(ss-1)*6)
%         if ss==1
%             plot(-1:2,mean(SubRevAcc{1,r},2),'--go','LineWidth',2)
%         else
%             plot(-1:2,mean(SubRevAcc{2,r},2),'--ro','LineWidth',2)
%         end
%         xlim([-1 2])
%         ylim([0 1])
%         xlabel('Trials')
%         ylabel('Accuracy')
%         title(sprintf('Run %d',r))
%     end
% end
% 
% % plot across runs
% subplot(2,6,[4 5 10 11])
% plot(-1:2,mean(AllRevAcc.P,2),'-go','LineWidth',2)
% hold on
% plot(-1:2,mean(AllRevAcc.A,2),'-ro','LineWidth',2)
% hold off
% xlim([-1 2])
% ylim([0 1])
% xlabel('Trials')
% ylabel('Accuracy')
% title(sprintf('Subject %d - All three runs',subj))
% legend('Placebo','Active','Location','southeast')
% 
% subplot(2,6,6)
% h = bar(Acc.Rev');
% ylabel('Accuracy')
% xlabel('Runs')
% title('On Reversal Trials')
% ylim([0 0.5])
% h(1).FaceColor = [0,1,0]; % placebo
% h(2).FaceColor = [1,0,0]; % active
% 
% subplot(2,6,12)
% h = bar(Acc.Nonrev');
% ylabel('Accuracy')
% xlabel('Runs')
% title('On Nonreversal Trials')
% ylim([0 1])
% h(1).FaceColor = [0,1,0]; % placebo
% h(2).FaceColor = [1,0,0]; % active
% set(gcf,'position',[100,300,1600,500])
% 
% saveas(gcf,fullfile(SubResDir,sprintf('BehOdorPredAcc_%s.bmp',subfolder)))
% 
% %% plot RT
% 
% figure
% % plot RT for each session and each run
% for ss = 1:nSess
%     for r = 1:nruns
%         subplot(2,5,r+(ss-1)*5)
%         if ss==1
%             plot(-1:2,nanmean(SubRevRT{1,r},2),'--go','LineWidth',2)
%         else
%             plot(-1:2,nanmean(SubRevRT{2,r},2),'--ro','LineWidth',2)
%         end
%         xlim([-1 2])
%         %ylim([500 1100])
%         xlabel('Trials')
%         ylabel('RT')
%         title(sprintf('Run %d',r))
%     end
% end
%     
% % plot RT across runs
% subplot(2,5,[4 5 9 10])
% plot(-1:2,nanmean(AllRevRT.P,2),'-go','LineWidth',2)
% hold on
% plot(-1:2,nanmean(AllRevRT.A,2),'-ro','LineWidth',2)
% hold off
% xlim([-1 2])
% xlabel('Trials')
% ylabel('RT')
% title(sprintf('Subject %d - All three runs',subj))
% legend('Placebo','Active','Location','southeast')
% 
% set(gcf,'position',[100,300,1600,500])
% saveas(gcf,fullfile(SubResDir,sprintf('BehOdorPredRT_%s.bmp',subfolder)))
% 
% %% plot RT distributions
% 
% figure
% SubRTsP = SubRTs(1,:);
% SubRTsP = cat(1,SubRTsP{:});
% 
% SubRTsA = SubRTs(2,:);
% SubRTsA = cat(1,SubRTsA{:});
% 
% histogram(SubRTsP)
% hold on
% histogram(SubRTsA)
% hold off
% legend('Placebo','Active')
% xlabel('Response Time (ms)')
% 
% set(gcf,'position',[100,300,700,500])
% saveas(gcf,fullfile(SubResDir,sprintf('BehOdorPredRTDist_%s.bmp',subfolder)))
% 
% 
% %% Plot Odor ID task acc and RT
% figure
% subplot(2,1,1)
% h = bar([AllID.acc.P;AllID.acc.A]');
% ylabel('Accuracy')
% xlabel('Runs')
% h(1).FaceColor = [0,1,0]; % placebo
% h(2).FaceColor = [1,0,0]; % active
% title(sprintf('Odor ID task Sub%d',subj))
% 
% subplot(2,1,2)
% h = bar([AllID.RT.P;AllID.RT.A]');
% ylabel('RT')
% xlabel('Runs')
% h(1).FaceColor = [0,1,0]; % placebo
% h(2).FaceColor = [1,0,0]; % active
% 
% set(gcf,'position',[100,300,700,500])
% saveas(gcf,fullfile(SubResDir,sprintf('BehOdorID_%s.bmp',subfolder)))
% 
% close all
fprintf('\n')

end % end of subs

%% Aggregate all datmat across subjects
All_datmat = vertcat(datmat{:});
All_datmat = All_datmat(:,[1:11,18:19]); % remove those timestamps
All_datmat = array2table(All_datmat);
All_datmat.Properties.VariableNames = {'Sub','Sess','TMS','Run','Trial','Cue','Odor','Rev','CorButton','ActualButton','RT','Resp','Acc'};
%writetable(All_datmat,fullfile(ResDir,'OdorPredRes.xlsx'))

%% Save variables 
save(fullfile(ResDir,'BehAllSubjects.mat'),'AggAcc','AggmeanRT','AggRTs','AggRevAcc','AggRevRT','AggID_acc','AggID_RT','Subs')



