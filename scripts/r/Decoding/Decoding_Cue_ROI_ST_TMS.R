
rm(list = ls())
library(varhandle)
library(tidyverse)
library(plyr)

source("~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/figs/ggplotSetup.R")

External = '/Volumes/QF10TB/SPEEDTMS_results'
subdirs = dir(file.path(External,'Decoding_ROI'))
nsubs = length(subdirs)

modelname_list = c('STCues_pseudoconcat_Nz','STCues_LSS')
nmodels = length(modelname_list)

mask_labels = c('OFC_cTBS','lOFC_cTBS','rOFC_cTBS',
  'LPFC_cTBS','lLPFC_cTBS','rLPFC_cTBS',
  'OFC_func','lOFC_med','rOFC_med','OFC_med')
n_masks = length(mask_labels)

dat_use_strings = c('beta') 
dat_scale_strings = c('raw','Mean_across_conds','Mean_across_all_voxels')
TMSconds = c('A','P')

dat_table = expand.grid(subdirs,modelname_list,mask_labels,
                        dat_scale_strings,TMSconds)
names(dat_table) = c('Sub','Model','Mask','Scale','TMS')
dat_table$acc = NA
dat_table$acc_clean = NA
dat_table$acc_trialn1 = NA
dat_table$acc_trial0 = NA
dat_table$acc_trial1 = NA
dat_table$acc_trial2 = NA

for(i in 1:nsubs){
  print(subdirs[i])
  for(j in 1:nmodels){
    for(k in 1:n_masks){
      for(m in 1:length(dat_scale_strings)){
        for(t in TMSconds){
          tmpdir = file.path(External,'Decoding_ROI',subdirs[i],modelname_list[j],
                             mask_labels[k],dat_use_strings,
                             dat_scale_strings[m],t)
          acc = as.numeric(read.csv(file.path(tmpdir,'Acc_LOOCV.csv'),header = F))
          acc_clean = as.numeric(read.csv(file.path(tmpdir,'Acc_LOOCV_clean.csv'),header = F))
          acc_trialn1 = as.numeric(read.table(file.path(tmpdir,'Acc_LOOCV_rev_trial-1.txt'),header = F))
          acc_trial0 = as.numeric(read.table(file.path(tmpdir,'Acc_LOOCV_rev_trial0.txt'),header = F))
          acc_trial1 = as.numeric(read.table(file.path(tmpdir,'Acc_LOOCV_rev_trial1.txt'),header = F))
          acc_trial2 = as.numeric(read.table(file.path(tmpdir,'Acc_LOOCV_rev_trial2.txt'),header = F))
          find_row = dat_table$Sub==subdirs[i] & 
            dat_table$Model==modelname_list[j] &
            dat_table$Mask==mask_labels[k] &
            dat_table$Scale==dat_scale_strings[m] &
            dat_table$TMS==t
          dat_table$acc[find_row] = acc
          dat_table$acc_clean[find_row] = acc_clean
          dat_table$acc_trialn1[find_row] = acc_trialn1
          dat_table$acc_trial0[find_row] = acc_trial0
          dat_table$acc_trial1[find_row] = acc_trial1
          dat_table$acc_trial2[find_row] = acc_trial2
}}}}}

dat_table$TMS = plyr::mapvalues(dat_table$TMS,from = c('P','A'),to = c('sham','cTBS'))

# save this dat_table
write.csv(dat_table, 'ROI-TMS-ST/Decoding_Cue_ROI_ST_TMS_Acc.csv', row.names=FALSE)

###########################################################

dat_table = read.csv('ROI-TMS-ST/Decoding_Cue_ROI_ST_TMS_Acc.csv')

options = c('acc','acc_clean','acc_trialn1','acc_trial0','acc_trial1','acc_trial2')

# plotting
for(i in options){
for (m in modelname_list){
  fig = dat_table %>%
    subset(Model == m) %>%
    ggplot(aes(y= get(i), x=TMS)) + 
    geom_boxplot(alpha=0.5,outlier.alpha = 0.3) +
    labs(y = "Decoding Accuracy (%)", x = '',title = '') + 
    #ylim(0, 100) + 
    common +
    facet_wrap(~Scale + Mask, nrow = 3) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0)) +
    theme(legend.position = "none" ) +
    geom_hline(yintercept=100/3, linetype="dotted", color = "black", size=0.5)
  pdf(paste0('ROI-TMS-ST/DecodingROIsST_TMSCue_',m,'_',i,'.pdf'),20,14)
  print(fig)
  dev.off()
}}

#################

# quick t-tests
# all failed !!!
use.dat = subset(dat_table, 
                 Model == 'STCues_LSS' & 
                   Mask == 'lOFC_cTBS' &
                   Scale %in% c('Mean_across_conds'))

dat_to_test_sham = subset(use.dat, TMS == 'sham')$acc_trial0
dat_to_test_cTBS = subset(use.dat, TMS == 'cTBS')$acc_trial0
wilcox.test(dat_to_test_sham,dat_to_test_cTBS,paired = T,alternative = 'g')
wilcox.test(dat_to_test_cTBS,mu = 1/3*100)
wilcox.test(dat_to_test_sham,mu = 1/3*100)


