
# load decoding accuracy & permutation test accuracy
# then do boostrap across subjects

rm(list = ls())
library(varhandle)
library(tidyverse)
library(plyr)

source("~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/figs/ggplotSetup.R")
External = '/Volumes/QF10TB/SPEEDTMS_results'
subdirs = dir(file.path(External,'Decoding_ROI'))
nsubs = length(subdirs)
modelname_list = c('Expect_pseudoconcat_Nz',
                   'Expect_pseudoconcat_Nz_1s',
                   'Expect_pseudoconcat_Nz_4s',
  'Expect_pseudoconcat_Nz_clean',
  'Expect_pseudoconcat_Nz_clean_addtd',
  'Expect_pseudoconcat_Nz_rev_trial-1',
  'Expect_pseudoconcat_Nz_rev_trial0',
  'Expect_pseudoconcat_Nz_rev_trial1',
  'Expect_pseudoconcat_Nz_rev_trial2')
nmodels = length(modelname_list)

mask_labels = c('OFC_cTBS','lOFC_cTBS','rOFC_cTBS',
  'LPFC_cTBS','lLPFC_cTBS','rLPFC_cTBS',
  'OFC_func','lOFC_med','rOFC_med','OFC_med')
n_masks = length(mask_labels)

dat_use_strings = c('beta') 
dat_scale_strings = c('raw','Mean_across_conds','Mean_across_all_voxels')

dat_table = expand.grid(subdirs,modelname_list,mask_labels,
                        dat_scale_strings)
names(dat_table) = c('Sub','Model','Mask','Scale')
dat_table$Acc = NA

for(i in 1:nsubs){
  print(subdirs[i])
  for(j in 1:nmodels){
    for(k in 1:n_masks){
      for(m in 1:length(dat_scale_strings)){
          tmpdir = file.path(External,'Decoding_ROI',subdirs[i],modelname_list[j],
                             mask_labels[k],dat_use_strings,
                             dat_scale_strings[m])
          csvfile = file.path(tmpdir,'Acc_LOOCV.csv')
          permutation_file = file.path(tmpdir,'PermutationAcc.csv')
          acc = as.numeric(read.csv(csvfile,header = F))
          dat_table$Acc[dat_table$Sub==subdirs[i] & 
                dat_table$Model==modelname_list[j] &
                dat_table$Mask==mask_labels[k] &
                dat_table$Scale==dat_scale_strings[m]] = acc
      }}}}


# save this dat_table
write.csv(dat_table, 'ROI/Decoding_ROI_Expect_Acc.csv', row.names=FALSE)

for (m in modelname_list){
  use.dat = subset(dat_table,Model == m)
  fig = ggplot(data = use.dat, aes(y=Acc, x=Mask)) + 
    geom_boxplot(alpha=0.5,outlier.alpha = 0) +
    labs(y = "Decoding Accuracy (%)", x = '',title = '') + 
    ylim(0, 100) + common +
    facet_wrap(~Scale) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90)) +
    theme(legend.position = "none" ) +
    geom_hline(yintercept=100/3, linetype="dotted", color = "black", size=0.5)

  pdf(paste0('ROI/DecodingROIsCue_',m,'.pdf'),15,5)
  print(fig)
  dev.off()

}


