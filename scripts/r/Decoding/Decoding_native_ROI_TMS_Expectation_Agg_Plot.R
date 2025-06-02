
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
TMSconds = c('A','P')

dat_table = expand.grid(subdirs,modelname_list,mask_labels,
                        dat_scale_strings,TMSconds)
names(dat_table) = c('Sub','Model','Mask','Scale','TMS')
dat_table$Acc = NA

for(i in 1:nsubs){
  print(subdirs[i])
  for(j in 1:nmodels){
    for(k in 1:n_masks){
      for(m in 1:length(dat_scale_strings)){
        for(t in TMSconds){
          tmpdir = file.path(External,'Decoding_native_ROI',subdirs[i],modelname_list[j],
                             mask_labels[k],dat_use_strings,
                             dat_scale_strings[m],t)
          csvfile = file.path(tmpdir,'Acc_LOOCV.csv')
          #permutation_file = file.path(tmpdir,'PermutationAcc.csv')
          acc = as.numeric(read.csv(csvfile,header = F))
          dat_table$Acc[dat_table$Sub==subdirs[i] & 
                dat_table$Model==modelname_list[j] &
                dat_table$Mask==mask_labels[k] &
                dat_table$Scale==dat_scale_strings[m] &
            dat_table$TMS==t] = acc
      }}}}}


# save this dat_table
write.csv(dat_table, 'ROI/Decoding_native_ROI_TMS_Expect_Acc.csv', row.names=FALSE)


#dat_table = read.csv('ROI_TMS/Decoding_ROI_TMS_Expect_Acc.csv')

# ignore these ROIs from plotting
dat_table = subset(dat_table, !Mask %in% c('lOFC_med','rOFC_med',
                                           'OFC_med','OFC_cTBS'))
# plotting
for (m in modelname_list){
  use.dat = subset(dat_table,Model == m)
  fig = ggplot(data = use.dat, aes(y=Acc, x=TMS)) + 
    geom_boxplot(alpha=0.5,outlier.alpha = 0) +
    labs(y = "Decoding Accuracy (%)", x = '',title = '') + 
    ylim(0, 100) + common +
    facet_wrap(~Scale + Mask) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90)) +
    theme(legend.position = "none" ) +
    geom_hline(yintercept=100/3, linetype="dotted", color = "black", size=0.5)
  pdf(paste0('ROI_TMS/DecodingNativeROIsTMSCue_',m,'.pdf'),15,14)
  print(fig)
  dev.off()
}


#######
rois = c('lOFC_cTBS','rOFC_cTBS','lLPFC_cTBS','rLPFC_cTBS')
roi_label = c('left OFC','right OFC','left LPFC','right LPFC')

use.dat = subset(dat_table, Model == 'Expect_pseudoconcat_Nz_clean_addtd' & 
                   Mask %in% rois &
                   Scale %in% c('Mean_across_conds'))
use.dat = use.dat %>% 
  arrange(Mask,Sub) %>%
  dplyr::mutate(paired = rep(1:(n()/2),each=2))

p = NULL
for(i in 1:length(rois)){
  p[[i]] = ggplot(data = subset(use.dat, Mask == rois[i]), 
                  aes(x = TMS, y = Acc, fill = TMS)) + 
    geom_boxplot(alpha=0.5,outlier.alpha = 0) +
    geom_line(aes(group=paired), position = position_dodge(0.2)) +
    geom_point(aes(fill=TMS,group=paired),size=2,shape=21,alpha=1,
               position = position_dodge(0.2)) +
    scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
    labs(y = "Decoding Accuracy (%)", x = '',title = roi_label[i]) + 
    ylim(0, 100) + common +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none" ) +
    geom_hline(yintercept=100/3, linetype="dotted", color = "black", size=0.5)
  #geom_segment(aes(x=1,y=95,xend=1,yend=100),size = 0.2)
}

pdf('ROI_TMS/CompareDecoding_allTrials.pdf',12,4)
fig = ggpubr::ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], 
                        ncol = 4, nrow = 1, align = "hv")
print(fig)
dev.off()

# test sig

roi = rois[1]
dat_sham = subset(use.dat, Mask == roi & TMS == 'sham')$Acc
dat_cTBS = subset(use.dat, Mask == roi & TMS == 'cTBS')$Acc
t.test(dat_sham,mu = 1/3*100,alternative = 'g')
t.test(dat_cTBS,mu = 1/3*100,alternative = 'l')
t.test(dat_sham,dat_cTBS,alternative = 'g',paired = T)



