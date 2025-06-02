
% This script does group analysis for all univariate GLM models

addpath('/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/ForExperiment')

%% iPE models

modelname = 'iPE';
connames = {'Reversal>Non-reversal','Cue','Response','Rev-Sess-Interaction','Rev-TMS-Interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)
connames = {'Rev_Nonrev_Real_vs_Sham','Rev_Nonrev_1st_vs_2nd'};
IndContrastIdx = [6,7;8,9]; 
Do_8_function_GroupGLM_Paired_Ttest(modelname,connames,IndContrastIdx)


%%
modelname = 'iPE_clean';
connames = {'Reversal>Non-reversal','Cue','Response','Rev-Sess-Interaction','Rev-TMS-Interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)
connames = {'Rev_Nonrev_Real_vs_Sham','Rev_Nonrev_1st_vs_2nd'};
IndContrastIdx = [6,7;8,9]; 
Do_8_function_GroupGLM_Paired_Ttest(modelname,connames,IndContrastIdx)


%% 
modelname = 'iPE_noResp';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)
connames = {'Rev_Nonrev_Real_vs_Sham','Rev_Nonrev_1st_vs_2nd'};
IndContrastIdx = [5,6;7,8]; 
Do_8_function_GroupGLM_Paired_Ttest(modelname,connames,IndContrastIdx)


%% iPE concat models
modelname = 'iPE_concat';
connames = {'Rev_vs_Nonrev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%% iPE pseudoconcat models
modelname = 'iPE_pseudoconcat';
connames = {'Reversal>Non-reversal','Rev-Sess-Interaction','Rev-TMS-Interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

modelname = 'iPE_pseudoconcat_noResp';
connames = {'Reversal>Non-reversal','Rev-Sess-Interaction','Rev-TMS-Interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

modelname = 'iPE_clean_pseudoconcat';
connames = {'Reversal>Non-reversal','Rev-Sess-Interaction','Rev-TMS-Interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)


%% iPE semi-concat model

modelname = 'iPE_semiconcat_noResp';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'iPE_clean_semiconcat_noResp';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'iPE_clean_semiconcat_noResp_Nz_rename';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)


%% iPE_fourpt models

modelname = 'iPE_fourpt';
connames = {'Rev_vs_Nonrev','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

modelname = 'iPE_fourpt_clean';
connames = {'Rev_vs_Nonrev','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'iPE_fourpt_concat';
connames = {'Rev_vs_Nonrev','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames);

%% iPE_fourpt pseudoconcat models

modelname = 'iPE_fourpt_pseudoconcat';
connames = {'Reversal>Non-reversal','Cue','Response','Rev-Sess-Interaction','Rev-TMS-Interaction',...
    'PostRev-Sess-Interaction','PostRev-TMS-Interaction','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames);

%%
modelname = 'iPE_fourpt_pseudoconcat_noResp';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction',...
    'PostRev-Sess-Interaction','PostRev-TMS-Interaction','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames);

%%
modelname = 'iPE_fourpt_clean_pseudoconcat_noResp';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction',...
    'PostRev-Sess-Interaction','PostRev-TMS-Interaction','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames);

%%
modelname = 'iPE_fourpt_semiconcat_noResp';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction',...
    'PostRev-Sess-Interaction','PostRev-TMS-Interaction','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)


%%
modelname = 'iPE_fourpt_clean_semiconcat_noResp';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction',...
    'PostRev-Sess-Interaction','PostRev-TMS-Interaction','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)


%%
modelname = 'iPE_fourpt_clean_semiconcat_noResp_2sec';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction',...
    'PostRev-Sess-Interaction','PostRev-TMS-Interaction','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'iPE_fourpt_clean_semiconcat_noResp_Nz';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction',...
    'PostRev-Sess-Interaction','PostRev-TMS-Interaction','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'iPE_fourpt_semiconcat_noResp_Nz';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction',...
    'PostRev-Sess-Interaction','PostRev-TMS-Interaction','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)


%%
modelname = 'Expect_fourpt_clean_semiconcat_noResp_Nz';
connames = {'Expect: Pre-rev > Reversal','Expect: Post-rev > Reversal','PreRev-Sess-Interaction',...
    'PreRev-TMS-Interaction','PostRev-Sess-Interaction','PostRev-TMS-Interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'iPE_2pt_clean_semiconcat_noResp_Nz';
connames = {'Reversal > Post-rev','Cue','PostRev-Sess-Interaction','PostRev-TMS-Interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)


%%
modelname = 'iPE_fourpt_clean_semiconcat_noResp_Nz_peak';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction',...
    'PostRev-Sess-Interaction','PostRev-TMS-Interaction','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'iPE_fourpt_clean_semiconcat_noResp_Nz_start';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction',...
    'PostRev-Sess-Interaction','PostRev-TMS-Interaction','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'Run1';
connames = {'Reversal>Non-reversal','Cue','Rev-Sess-Interaction','Rev-TMS-Interaction',...
    'PostRev-Sess-Interaction','PostRev-TMS-Interaction','Reversal_vs_Post-rev'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'Run1_iPE_combSess';
connames = {'iPE_para'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'Run1_iPE_combSess_addTMS';
connames = {'iPE_para','iPE_para_TMS_interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'AllRun_iPE_combSess_addTMS';
connames = {'iPE_para','iPE_para_TMS_interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'Run1_iPE_combSess_addTMS_v4a';
connames = {'iPE_para','iPE_para_TMS_interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

%%
modelname = 'Run1_iPE_combSess_addTMS_v7a';
connames = {'iPE_para','iPE_para_TMS_interaction'};
Do_8_function_GroupGLM_OneSample_Ttest(modelname,connames)

