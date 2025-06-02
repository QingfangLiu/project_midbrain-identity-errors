
rm(list = ls())
library(varhandle)
library(tidyverse)
library(plyr)

source("~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/figs/ggplotSetup.R")
External = '/Volumes/QF10TB/SPEEDTMS_results'
subdirs = dir(file.path(External,'Decoding_native_ROI'))
nsubs = length(subdirs)
modelname_list = c('Expect_pseudoconcat_Nz',
                  'Expect_pseudoconcat_Nz_clean',
                  'Expect_pseudoconcat_Nz_clean_addtd')
nmodels = length(modelname_list)

mask_labels = c('TMS_lOFC','TMS_rOFC')
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
          tmpdir = file.path(External,'Decoding_native_ROI',subdirs[i],modelname_list[j],
                             mask_labels[k],dat_use_strings,
                             dat_scale_strings[m])
          csvfile = file.path(tmpdir,'Acc_LOOCV.csv')
          #permutation_file = file.path(tmpdir,'PermutationAcc.csv')
          acc = as.numeric(read.csv(csvfile,header = F))
          dat_table$Acc[dat_table$Sub==subdirs[i] & 
                dat_table$Model==modelname_list[j] &
                dat_table$Mask==mask_labels[k] &
                dat_table$Scale==dat_scale_strings[m]] = acc
      }}}}


# save this dat_table
write.csv(dat_table, 'ROI/Decoding_native_ROI_Expect_Acc.csv', row.names=FALSE)

