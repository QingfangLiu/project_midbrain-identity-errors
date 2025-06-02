
# This script conducts testing & plotting on ROI-based MVPA 
# also saving the summarized data into xlsx

rm(list = ls())

# define directory paths
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
MVPADir = file.path(ResDir,'MVPA_ROI')
AnaDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis'
SubInfoDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/SubInfo'

# ggplot setup
SetupFile = file.path(AnaDir,'Scripts_R','Setup.R')
source(SetupFile)

# source data 
sourceDir = '~/Desktop/SourceData'

##################################################################

SubInfo = read.xlsx(file.path(SubInfoDir,'SubjectConds.xlsx'))
UseTMSConds = SubInfo$Cond[SubInfo$Excluded==0]
Subs = SubInfo$Sub[SubInfo$Excluded==0]
nsubs = length(Subs)

ROIs = c('Func_OFC','Func_LPFC')    # ROI names in savd results
#ROI.labs = c('Lateral OFC','LPFC')  # labels to use in plotting
nROIs = length(ROIs)

All_corr_vals_df = NULL
sub_corr_rev_array = array(NA,c(nsubs,nROIs,6,2,12))

ctr = 0
for(j in 1:nsubs){
subno = paste0('Sub',Subs[j])
print(subno)
for(r in 1:nROIs){
ctr = ctr + 1
roi = ROIs[r]
tmp_dir = file.path(MVPADir,subno,'STCues_LSS',roi)
sub_corr_vars = R.matlab::readMat(file.path(tmp_dir,'corr_vars.mat'))
sub_corr = sub_corr_vars$corr.res
sub_corr_rev = sub_corr_vars$corr.res.rev   # also for each rev
sub_corr_rev_array[j,r,,,] = sub_corr_rev
sub_corr = as.numeric(sub_corr)
sub_corr_vals_df = data.frame('Sess' = rep(c(1,1,1,2,2,2),2),
                              'Run' = rep(c(1,2,3,1,2,3),2),
                              'Type' = rep(c('pre-rev','post-rev'),each=6),
                              'corr' = sub_corr)
sub_corr_vals_df$Sub = Subs[j]
sub_corr_vals_df$ROI = roi
All_corr_vals_df[[ctr]] = sub_corr_vals_df

}
}

Corr_vals_df <- do.call(rbind, All_corr_vals_df)

# convert Sess to TMS
# be very very careful at this step !!!
Corr_vals_df$TMS = Corr_vals_df$Sess
for(j in 1:nsubs){
  if(UseTMSConds[j]=='AP'){
    Corr_vals_df$TMS[Corr_vals_df$Sub==Subs[j]] = 3 - 
      Corr_vals_df$TMS[Corr_vals_df$Sub==Subs[j]]
  }
}

Corr_vals_df$TMS = plyr::mapvalues(Corr_vals_df$TMS,from = c(1,2),to = c('sham','cTBS'))
Corr_vals_df$TMS = factor(Corr_vals_df$TMS)
Corr_vals_df$ROI = factor(Corr_vals_df$ROI,levels = ROIs)

# reduce the matrix by pre-rev - post-rev (identity expectation index)
Corr_diff_vals = subset(Corr_vals_df,Type == 'pre-rev')$corr -
  subset(Corr_vals_df,Type == 'post-rev')$corr
Corr_diff_vals_other = subset(Corr_vals_df,Type == 'pre-rev')
Corr_diff_vals_df = Corr_diff_vals_other
Corr_diff_vals_df$corr = Corr_diff_vals
Corr_diff_vals_df$Type = 'pre-post'

# summarize across runs
Corr_diff_vals_df_summary = Corr_diff_vals_df %>%
  group_by(TMS,Sub,ROI,Sess) %>%
  dplyr::summarise(Representation = mean(corr)) %>% 
  arrange(Sub,ROI) %>%
  ungroup() %>%
  dplyr::mutate(paired = rep(1:(n()/2),each=2))


Corr_diff_vals_df_summary = Corr_diff_vals_df_summary %>%
  subset(ROI == 'Func_OFC')

## save source data
write.xlsx(Corr_diff_vals_df_summary,
           file = file.path(sourceDir,'MVPA.xlsx'))

################ plots ###################

p = NULL
for(i in 1:nROIs){
  p[[i]] = ggplot(data = subset(Corr_diff_vals_df_summary,ROI == ROIs[i]),
                mapping = aes(x = TMS, y = Representation, fill = TMS)) + 
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=TMS,group=paired),size=2,shape=21,alpha=1,
             position = position_dodge(0.2)) +
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  labs(y = "", x = "", title = '') + 
  common +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.text.x = element_blank(),legend.position = "none" ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted")
}

high_y2 = 0.23
low_y2 = 0.22
p[[1]] = p[[1]] +
  geom_segment(aes(x=1,y=high_y2,xend=2,yend=high_y2),linewidth = 0.2) +
  geom_segment(aes(x=1,y=low_y2,xend=1,yend=high_y2),linewidth = 0.2) +
  geom_segment(aes(x=2,y=low_y2,xend=2,yend=high_y2),linewidth = 0.2) +
  annotate(geom = "text", x = c(1.5), y = c(0.235), label = '*', size = 6)

pdf(file.path(MVPADir,'MVPA_OFC.pdf'),3.5,4)
par(mar=c(1,1,1,1))
print(p[[1]])
dev.off()

high_y2 = 0.18
low_y2 = 0.17
p[[2]] = p[[2]] +
  geom_segment(aes(x=1,y=high_y2,xend=2,yend=high_y2),linewidth = 0.2) +
  geom_segment(aes(x=1,y=low_y2,xend=1,yend=high_y2),linewidth = 0.2) +
  geom_segment(aes(x=2,y=low_y2,xend=2,yend=high_y2),linewidth = 0.2) +
  annotate(geom = "text", x = c(1.5), y = c(0.2), label = 'n.s.', size = 6)

pdf(file.path(MVPADir,'MVPA_LPFC.pdf'),3.5,4)
par(mar=c(1,1,1,1))
print(p[[2]])
dev.off()

for(roi in ROIs){
  print(roi)
  x_sham = subset(Corr_diff_vals_df_summary,
                  ROI == roi & TMS == 'sham')$Representation
  x_cTBS = subset(Corr_diff_vals_df_summary,
                  ROI == roi & TMS == 'cTBS')$Representation
  print(wilcox.test(x_sham)) # representation > 0 from sham trials?
  print(wilcox.test(x_cTBS)) # representation > 0 from cTBS trials?
  print(wilcox.test(x_sham,x_cTBS,paired = T,alternative = 'g'))
}


#####################################################
## so far is what I reported in the paper ##
## below is more analysis 
#####################################################

# test TMS and Sess effects jointly
# lmer1: including both TMS and Sess
# lmer2: including only TMS
# lmer3: including only Sess

# focus on OFC
lmer1_OFC <- lmer(Representation ~ Sub + TMS + Sess + (1 | Sub), 
                 data = subset(Corr_diff_vals_df_summary,ROI=='Func_OFC'))
lmer2_OFC <- lmer(Representation ~ Sub + TMS + (1 | Sub), 
                  data = subset(Corr_diff_vals_df_summary,ROI=='Func_OFC'))
lmer3_OFC <- lmer(Representation ~ Sub + Sess + (1 | Sub), 
                  data = subset(Corr_diff_vals_df_summary,ROI=='Func_OFC'))
anova(lmer1_OFC,lmer2_OFC) # Sess useful? (with TMS) p = 0.38
anova(lmer1_OFC,lmer3_OFC) # TMS useful? (with Sess) p = 0.058

# focus on LPFC
lmer1_LPFC <- lmer(Representation ~ Sub + TMS + Sess + (1 | Sub), 
                  data = subset(Corr_diff_vals_df_summary,ROI=='Func_LPFC'))
lmer2_LPFC <- lmer(Representation ~ Sub + TMS + (1 | Sub), 
                  data = subset(Corr_diff_vals_df_summary,ROI=='Func_LPFC'))
lmer3_LPFC <- lmer(Representation ~ Sub + Sess + (1 | Sub), 
                  data = subset(Corr_diff_vals_df_summary,ROI=='Func_LPFC'))
anova(lmer1_LPFC,lmer2_LPFC) # Sess useful? (with TMS) p = 0.0095
anova(lmer1_LPFC,lmer3_LPFC) # TMS useful? (with Sess) p = 0.46


# save summarized result into xlsx
Corr_diff_vals_df_summary %>%
  select(Sub,ROI,TMS,Sess,Representation) %>%
  write.xlsx(file = file.path(MVPADir,'IdentityRep_SubsTMS.xlsx'),
             rowNames = F)


## results below have not been updated from Oct 2023 

######### test if cTBS effect on identity expectation is correlated with
######### other between-subject variables 

roi = 'lOFC_SimiComp'
x_sham = subset(Corr_diff_vals_df_summary,ROI == roi & TMS == 'sham')$Representation
x_cTBS = subset(Corr_diff_vals_df_summary,ROI == roi & TMS == 'cTBS')$Representation
rep_effect = x_cTBS - x_sham

# compare TMS effect on different orders
group1 = SubInfo$Cond[SubInfo$Excluded==0]=='PA'
group2 = SubInfo$Cond[SubInfo$Excluded==0]=='AP'
t.test(rep_effect[group1],rep_effect[group2])

# correlate age with TMS effect
cor.test(SubInfo$Age[SubInfo$Excluded==0],rep_effect)

# correlate TMS intensity with TMS effect
cor.test(SubInfo$TMS[SubInfo$Excluded==0],rep_effect)

# compare TMS effect on different sex
group1 = SubInfo$Sex[SubInfo$Excluded==0]=='F'
group2 = SubInfo$Sex[SubInfo$Excluded==0]=='M'
t.test(rep_effect[group1],rep_effect[group2])





####### Check identity representation averaged across TMS conditions

Corr_diff_vals_df_summary_both = Corr_diff_vals_df %>%
  group_by(Sub,ROI) %>%
  dplyr::summarise(Representation = mean(corr)) %>% 
  arrange(Sub,ROI) %>%
  ungroup() 

ggplot(data = Corr_diff_vals_df_summary_both,
       mapping = aes(ROI,Representation)) + 
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  labs(y = 'Representation', x = '') + 
  common +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted")

# test against zero
for(roi in ROIs){
  print(roi)
  x = subset(Corr_diff_vals_df_summary,ROI == roi)$Representation
  print(wilcox.test(x))
  print(t.test(x))
}

####### another way to look at the results is to look at pre-rev
####### and post-rev similarity separately

Corr_vals_df_summary = Corr_vals_df %>%
  group_by(TMS,Sub,ROI,Type) %>%
  dplyr::summarise(Representation = mean(corr)) %>% 
  arrange(Sub,ROI,Type) %>%
  ungroup() %>%
  dplyr::mutate(paired = rep(1:(n()/2),each=2))

Corr_vals_df_summary$Type = factor(Corr_vals_df_summary$Type,
                                   levels = c('pre-rev','post-rev'))

p = NULL
for(i in 1:nROIs){
  roi = ROIs[i]
  print(roi)
  p[[i]] = ggplot(data = subset(Corr_vals_df_summary, ROI == ROIs[i]),
                  mapping = aes(x = TMS, y = Representation, fill = TMS)) + 
    geom_line(aes(group=paired), position = position_dodge(0.2)) +
    geom_point(aes(fill=TMS,group=paired),size=2,shape=21,alpha=1,
               position = position_dodge(0.2)) +
    facet_wrap(~Type) +
    geom_boxplot(alpha=0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
    labs(y = 'Pattern Similarity', x = '', title = ROIs[i]) + 
    common +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none" ) 
  
  x_sham = subset(Corr_vals_df_summary,
                  ROI == roi & TMS == 'sham' & Type == 'pre-rev')$Representation
  x_cTBS = subset(Corr_vals_df_summary,
                  ROI == roi & TMS == 'cTBS' & Type == 'pre-rev')$Representation
  print(wilcox.test(x_sham,x_cTBS,paired = T))
  
  x_sham = subset(Corr_vals_df_summary,
                  ROI == roi & TMS == 'sham' & Type == 'post-rev')$Representation
  x_cTBS = subset(Corr_vals_df_summary,
                  ROI == roi & TMS == 'cTBS' & Type == 'post-rev')$Representation
  print(wilcox.test(x_sham,x_cTBS,paired = T))
}

#pdf('figs/MVPSA_Comp_pre_post_rev.pdf',6,8)
ggpubr::ggarrange(p[[1]], p[[2]], align = 'h', nrow = 2, ncol = 1)
#dev.off()

#############
# go back to work on the array that contains rev dim
#############

sub_corr_rev_df = structure(sub_corr_rev_array,
                       .Dim = c(nsubs,nROIs,6,2,12), 
                       .Dimnames = structure(list(Subs = Subs,
                                                  ROI = ROIs, 
                                                  Run = 1:6,
                                                  Type = c('pre-rev','post-rev'),
                                                  Reversal = 1:12), 
                      .Names = c("Subs", "ROI" ,"Run", "Type", "Reversal")))
sub_corr_rev_df = plyr::adply(sub_corr_rev_df, c(1,2,3,4,5))
names(sub_corr_rev_df) = c(names(sub_corr_rev_df)[1:5],'corr')

sub_corr_rev_df <- sub_corr_rev_df %>%
  mutate(Sess = recode(Run, `1` = 1, `2` = 1, `3` = 1, 
                            `4` = 2, `5` = 2, `6` = 2)) %>%
  mutate(Run = recode(Run, `1` = 1, `2` = 2, `3` = 3, 
                      `4` = 1, `5` = 2, `6` = 3)) %>%
  mutate(Run = factor(Run)) %>%
  mutate(TMS = Sess)

for(j in 1:length(Subs)){
  if(UseTMSConds[j]=='AP'){
    sub_corr_rev_df$TMS[sub_corr_rev_df$Subs==Subs[j]] = 3 - 
      sub_corr_rev_df$TMS[sub_corr_rev_df$Subs==Subs[j]]
  }
}

sub_corr_rev_df$TMS = plyr::mapvalues(sub_corr_rev_df$TMS,
                                      from = c(1,2),to = c('sham','cTBS'))
sub_corr_rev_df$TMS = factor(sub_corr_rev_df$TMS)

sub_corr_rev_diff = subset(sub_corr_rev_df,Type == 'pre-rev')$corr -
  subset(sub_corr_rev_df,Type == 'post-rev')$corr
sub_corr_rev_df_other = subset(sub_corr_rev_df,Type == 'pre-rev')
sub_corr_rev_diff_df = sub_corr_rev_df_other
sub_corr_rev_diff_df$corr = sub_corr_rev_diff
sub_corr_rev_diff_df$Type = 'pre-post'

# save as csv file
write_csv(sub_corr_rev_diff_df,'IdentityRepresent_per_rev.csv')

sub_corr_rev_df_summary = sub_corr_rev_diff_df %>%
  subset(ROI == ROIs[1]) %>%
  group_by(Subs,TMS,Run) %>%
  dplyr::summarise(Representation = mean(corr)) %>% 
  ungroup() 

p = ggplot(data = sub_corr_rev_df_summary,
       mapping = aes(TMS,Representation)) + 
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  facet_wrap(~ Run) +
  labs(y = 'Representation', x = '') + 
  common

# compare representation across runs
pdf('figs/MVPSA_Comp_runs.pdf',8,6)
print(p)
dev.off()

# load ST iPE estimates from GLM
# beta extracted using matlab code
wd = dirname(getwd())
all_betas_iPE = array(NA,c(nsubs,2,2,2,3,12))

for(j in 1:nsubs){
  subno = paste0('Sub',Subs[j])
  print(subno)
  tmp_dir = file.path(wd,'NeuralAnalyzeRes',subno)
  sub_beta_mat = R.matlab::readMat(file.path(tmp_dir,'betas_per_rev_from_ST.mat'))
  sub_betas = sub_beta_mat$all.betas
  sub_betas_iPE = sub_betas[,,,,,2] - sub_betas[,,,,,3] # here using rev - (rev+1)
  all_betas_iPE[j,,,,,] = sub_betas_iPE
  
}

betas_iPE_df = structure(all_betas_iPE,
                            .Dim = c(nsubs,2,2,2,3,12), 
                            .Dimnames = structure(list(Subs = Subs,
                                                       Model = c('ST_Outcome_LSS','STOdorD_pseudoconcat_Nz'),
                                                       ROI = c('Midbrain','LPFC'), 
                                                       Sess = 1:2,
                                                       Run = 1:3,
                                                       Reversal = 1:12), 
                                                  .Names = c("Subs", "Model", "ROI" ,"Sess", "Run", "Reversal")))
betas_iPE_df = plyr::adply(betas_iPE_df, c(1,2,3,4,5,6))
names(betas_iPE_df) = c(names(betas_iPE_df)[1:6],'iPE')
betas_iPE_df$TMS = unfactor(betas_iPE_df$Sess)

for(j in 1:length(Subs)){
  if(UseTMSConds[j]=='AP'){
    betas_iPE_df$TMS[betas_iPE_df$Subs==Subs[j]] = 3 - 
      betas_iPE_df$TMS[betas_iPE_df$Subs==Subs[j]]
  }
}

betas_iPE_df$TMS = plyr::mapvalues(betas_iPE_df$TMS,
                                      from = c(1,2),to = c('sham','cTBS'))
betas_iPE_df$TMS = factor(betas_iPE_df$TMS)

# combine the three dfs and do some correlations

Use_identity_df = sub_corr_rev_diff_df %>%
  subset(ROI == 'lOFC_SimiComp') %>%
  select(Subs,TMS,Run,Reversal,corr)

Use_beta_df_Midbrain = betas_iPE_df %>%
  subset(Model == 'ST_Outcome_LSS' & ROI == 'Midbrain') %>%
  select(Subs,TMS,Run,Reversal,iPE)

Use_beta_df_LPFC = betas_iPE_df %>%
  subset(Model == 'ST_Outcome_LSS' & ROI == 'LPFC') %>%
  select(Subs,TMS,Run,Reversal,iPE)

# 3 dfs should have been organized in an identical way
Use_combined = cbind(Use_identity_df,
                     Use_beta_df_Midbrain$iPE,
                     Use_beta_df_LPFC$iPE)  
names(Use_combined) = c(names(Use_combined)[1:4],'identity','Midbrain','LPFC')

df_corr <- Use_combined %>% select(identity, Midbrain, LPFC)
cor(df_corr)
pairs(df_corr)

# correlation within each subject
cors = numeric(nsubs)
for(j in Subs){
  use_df = subset(Use_combined,Subs == j)
  df_corr <- use_df %>% select(identity, Midbrain, LPFC)
  cor_mat <- cor(df_corr)
  cors[j] = cor_mat[1,2] # focus on b/t identity & Midbrain
}

hist(cors)
summary(cors)



