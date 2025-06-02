
% The purpose of this script is to add contrast to the GLM estimation
% result, based on two conditions: Outcomes_iPE and Outcomes_noniPE

% Note that this is a pseudoconcat model, so SPM considers only one
% session. I need to manually specify multiple contrasts.

%% Setup
modelname = 'iPE_fourpt_clean_semiconcat_noResp_Nz';
modeldir = fullfile(studydir,subno,'fxUnivariate',modelname); 

%% define contrast names and weights
% compare all reversal trials with nonreversal trials, across all runs
contr_name{1} = 'Reversal>Non-reversal'; 
contr{1} = repmat([3,-1,-1,-1],1,2);
sessrep{1} = 'none';  

% find cue-related neural activation
contr_name{2} = 'Cue'; 
contr{2} = [zeros(1,8),1];
sessrep{2} = 'none';

% interaction between session-order and reversal types
% this contrast indicates how much more difference between reversal and
% non-reversal trials in the first session than in the second session
contr_name{3} = 'Rev-Sess-Interaction'; 
contr{3} = [3,-1,-1,-1,-3,1,1,1];
sessrep{3} = 'none';

% interaction between TMS condition and reversal types
% this contrast indicates how much more difference between reversal and
% non-reversal trials in the placebo condition than in the active condition
contr_name{4} = 'Rev-TMS-Interaction'; 
tmpweights = [3,-1,-1,-1,-3,1,1,1];

if strcmp(Conds{Subs==subj},'AP')
    tmpweights = tmpweights*(-1);
end

contr{4} = tmpweights;
sessrep{4} = 'none';

% interaction between session-order and post-reversal change
contr_name{5} = 'PostRev-Sess-Interaction'; 
contr{5} = [[1,0,-1,0],[-1,0,1,0]];
sessrep{5} = 'none';

% interaction between TMS condition and post-reversal change
contr_name{6} = 'PostRev-TMS-Interaction'; 
tmpweights = [[1,0,-1,0],[-1,0,1,0]];

if strcmp(Conds{Subs==subj},'AP')
    tmpweights = tmpweights*(-1);
end

contr{6} = tmpweights;
sessrep{6} = 'none';

% compare all reversal trials with their subsequent trials, across all runs
contr_name{7} = 'Reversal > Post-rev'; 
contr{7} = repmat([1,0,-1,0],1,2);
sessrep{7} = 'none';


%% SPM contrast
matlabbatch = cell(1,length(contr_name));
for c = 1:length(contr_name)
    matlabbatch{c}.spm.stats.con.spmmat = {fullfile(modeldir,'SPM.mat')};
    matlabbatch{c}.spm.stats.con.consess{1}.tcon.name = contr_name{c};
    matlabbatch{c}.spm.stats.con.consess{1}.tcon.weights = contr{c};
    matlabbatch{c}.spm.stats.con.consess{1}.tcon.sessrep = sessrep{c};
    matlabbatch{c}.spm.stats.con.delete = 0; 
    if c == 1
        matlabbatch{c}.spm.stats.con.delete = 1; % delete existing contrasts if rerun this script
    end
end

%% run the jobs
fname = fullfile(Jobdir, sprintf('AddContrasts_%s_%s.mat',modelname,subno));
save(fname, 'matlabbatch');
spm_jobman('run', matlabbatch);  
clear matlabbatch
