
% This script is used to summarize data from screening session (or combined
% with day1) across all subjects
% updated 2/21/23

clc; clear; close all

studydir = '/Users/qingfangliu/Experiment';
parentDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis';

SubInfoFile = fullfile(parentDir,'SubjectConds.XLSX');
SubInfo = readtable(SubInfoFile) ;
Subs = SubInfo.Sub;             % subject index
OdorDay1 = SubInfo.OdorDay1;    % whether use day1 odor ratings
Excluded = SubInfo.Excluded;    % whether excluding from analysis

Subs = Subs(Excluded==0);       % excluding subs from analysis
OdorDay1 = OdorDay1(Excluded==0);
nSubs = size(Subs,1);

SubCtr = 0;
saveDir = fullfile(parentDir,'ScreeningRes');

%%
for subj = Subs'
    
    SubCtr = SubCtr + 1;
    subno = sprintf('Sub%d',subj);
    fprintf('Subject %d\n',subj)
    subdir = fullfile(studydir,subno);
    
    if OdorDay1(ismember(Subs,subj))==1
        file = dir(fullfile(studydir,subno,'DAY1','complete*'));
    else
        file = dir(fullfile(studydir,subno,'DAY0','complete*'));
    end
    
    load(fullfile(file.folder,file.name))

    for i = 1:3
        odor_name{SubCtr,i} = res.odornames{res.odors(i),1};
    end

    % follow-up pleasantness ratings
    for i = 1:3 % loop three odors
        tmploc = res.pleasantness_ratings(:,1)==i;
        plea_ratings(i,:,SubCtr) = res.pleasantness_ratings(tmploc,5);
    end
    
    % intensity ratings
    for i = 1:3 % loop three odors
        tmploc = res.intensity_ratings(:,1)==i;
        inten_ratings(i,:,SubCtr) = res.intensity_ratings(tmploc,5);
    end
    
    % discrimination 
    discrimination(SubCtr,1) = mean(res.discrimination(:,13)); % mean discrimination score
        
    % similarity ratings
    similarity(SubCtr,1) = mean(res.similarity_ratings(:,9)); % mean similarity score
end

Avgplea_ratings = squeeze(median(plea_ratings,2)); % median over three times of ratings
Avginten_ratings = squeeze(median(inten_ratings,2)); % median over three times of ratings

OutMat = [Avgplea_ratings',Avginten_ratings',discrimination,similarity];
OutMat = num2cell(OutMat);
OutMat = [odor_name,OutMat];
OutMat = [num2cell(Subs),OutMat];

OutMat = cell2table(OutMat,...
    "VariableNames",["Sub" "Odor1" "Odor2" "Odor3"...
    "plea1" "plea2" "plea3" "inten1" "inten2" "inten3"...
    "disc" "simi"]);

writetable(OutMat,fullfile(saveDir,'ScreeningSubOdors.csv'))

%% Plotting

%% Pleasantness ratings

Avgplea_ratings = sort(Avgplea_ratings,1); % sort three odors based on ratings (be careful when using three odors elsewhere)
h = figure('Position', [10 10 1200 1000]);
subplot(3,1,2)
bar(Avgplea_ratings')
title('Individual Pleasantness Ratings')
xlabel('Subjects')

jitter = 0.2;
subplot(3,2,1)
boxplot(Avgplea_ratings','Labels',{'Lowest','Middle','Highest'})
hold on
scatter(ones(1,nSubs),Avgplea_ratings(1,:), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nSubs)*2,Avgplea_ratings(2,:), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nSubs)*3,Avgplea_ratings(3,:), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
title('Pleasantness Ratings')
xlabel('Odors')

%% Intensity ratings

Avginten_ratings = sort(Avginten_ratings,1); % sort three odors based on ratings (be careful when using three odors elsewhere)
subplot(3,1,3)
bar(Avginten_ratings')
title('Individual Intensity Ratings')
xlabel('Subjects')

jitter = 0.2;
subplot(3,2,2)
boxplot(Avginten_ratings','Labels',{'Lowest','Middle','Highest'})
hold on
scatter(ones(1,nSubs),Avginten_ratings(1,:), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nSubs)*2,Avginten_ratings(2,:), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
scatter(ones(1,nSubs)*3,Avginten_ratings(3,:), 'filled', 'jitter', 'on', 'jitterAmount', jitter);
title('Intensity Ratings')
xlabel('Odors')

saveas(h,fullfile(saveDir,'ScreeningOdorRatings_Result.png'))


