
%%
clc
close all
clear 

%% subject info
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'; 
SubInfoDir = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/SubInfo';

%%
SubInfoFile = fullfile(SubInfoDir,'SubjectConds.XLSX');
SubInfo = readtable(SubInfoFile) ;
Subs = SubInfo.Sub;  
Conds = SubInfo.Cond; % TMS condition
OdorDay1 = SubInfo.OdorDay1; % whether use day1 odor ratings
UseDay4 = SubInfo.UseDay4; % whether use data from day4
Excluded = SubInfo.Excluded; % whether excluding from analysis

Subs = Subs(Excluded==0); % excluding subs from analysis
nSubs = length(Subs); % number of subjects
Conds = Conds(Excluded==0); 

%% Specify contrast images
modelname = 'iPE_fourpt_clean_semiconcat_noResp_Nz';


%% Analyze extracted beta values

clear aggMeanBOLD aggMidBetas

for j = 1:nSubs
    tempMidBetas = fullfile(ResDir,'Neural',sprintf('Sub%d',Subs(j)),sprintf('Betas_%s_others.mat',modelname));
    load(tempMidBetas)
    aggMeanBOLD(:,:,:,j) = MeanBOLD;    
    aggMidBetas{j} = beta_vals;
end


%% organize into a dataframe to do mixed-effect modeling

beta_file_name = sprintf('Betas_%s_others.mat',modelname);
save(fullfile(ResDir,'Neural',beta_file_name),'aggMeanBOLD')



