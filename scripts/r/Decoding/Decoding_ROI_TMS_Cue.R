
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
                   'Expect_pseudoconcat_Nz_clean_addtd')

nmodels = length(modelname_list)

mask_labels = c('OFC_cTBS','lOFC_cTBS','rOFC_cTBS',
  'LPFC_cTBS','lLPFC_cTBS','rLPFC_cTBS',
  'OFC_func','lOFC_med','rOFC_med','OFC_med')
n_masks = length(mask_labels)

dat_use_strings = c('beta') 
dat_scale_strings = c('raw','Mean_across_conds','Mean_across_all_voxels')
TMSconds = c('A','P')

###########################################################
# skip this block when coming back for plotting

dat_table = expand.grid(subdirs,modelname_list,mask_labels,
                        dat_scale_strings,TMSconds)
names(dat_table) = c('Sub','Model','Mask','Scale','TMS')
dat_table$Acc = NA
dat_table[1,7:1006] = NA

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
          permutation_acc = read.csv(file.path(tmpdir,'PermutationAcc.csv'),header = F)
          find_row = dat_table$Sub==subdirs[i] & 
            dat_table$Model==modelname_list[j] &
            dat_table$Mask==mask_labels[k] &
            dat_table$Scale==dat_scale_strings[m] &
            dat_table$TMS==t
          dat_table$Acc[find_row] = as.numeric(acc)
          dat_table[find_row,7:1006] = as.numeric(permutation_acc$V1)
}}}}}

dat_table$TMS = plyr::mapvalues(dat_table$TMS,from = c('P','A'),to = c('sham','cTBS'))
write.csv(dat_table, 'ROI-TMS/Decoding_ROI_TMS_Expect_Acc.csv', row.names=FALSE)

###########################################################

dat_table = read.csv('ROI-TMS/Decoding_ROI_TMS_Expect_Acc.csv')
dat_table$Mask = factor(dat_table$Mask, levels = mask_labels) # specify mask orders
dat_table$Scale = factor(dat_table$Scale, levels = dat_scale_strings)

# plotting
for (m in modelname_list){
  fig = dat_table %>%
    subset(Model == m) %>%
    ggplot(aes(y=Acc, x=TMS)) +
    facet_wrap(~Scale + Mask,nrow = 3) +
    geom_boxplot(alpha=0.5,outlier.alpha = 0.3) +
    labs(y = "Decoding Accuracy (%)", x = '',title = '') + 
    ylim(0, 100) + common +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90)) +
    theme(legend.position = "none" ) +
    geom_hline(yintercept=100/3, linetype="dotted", color = "black", size=0.5)
  
  pdf(paste0('ROI-TMS/DecodingROIsTMSCue_',m,'.pdf'),20,10)
  print(fig)
  dev.off()
}


####### plot for paper

use.model = 'Expect_pseudoconcat_Nz_1s'
use.scale = 'Mean_across_conds'#'raw'
use.mask = 'lOFC_cTBS'

rois = c('lOFC_cTBS','rOFC_cTBS','lLPFC_cTBS','rLPFC_cTBS','OFC_func')
roi_label = c('left OFC','right OFC','left LPFC','right LPFC','OFC_func')

use.dat = dat_table %>% 
  subset(Model == use.model & 
           Mask %in% rois &
           Scale %in% c(use.scale)) %>%
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

pdf('ROI-TMS/CompareDecoding_allTrials.pdf',12,4)
fig = ggpubr::ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], 
                  ncol = 4, nrow = 1, align = "hv")
print(fig)
dev.off()

# test sig roughly using t-tests here
dat_sham = subset(use.dat, Mask == use.mask & TMS == 'sham')$Acc
dat_cTBS = subset(use.dat, Mask == use.mask & TMS == 'cTBS')$Acc
t.test(dat_sham,mu = 1/3*100,alternative = 'g')
t.test(dat_cTBS,mu = 1/3*100,alternative = 'g')
wilcox.test(dat_sham,dat_cTBS,alternative = 'g',paired = T)

###################################################
# test sig using permutation + bootstrap approach
###################################################

use_dat_table = subset(dat_table, 
                       Model == use.model & 
                       Mask == use.mask &
                       Scale == use.scale)
  
# reoranize permutated acc
perm_acc_df = data.frame(
  'Sub' = rep(use_dat_table$Sub,each=1000),
  'TMS' = rep(use_dat_table$TMS,each=1000),
  'Acc' = as.numeric(t(as.matrix(use_dat_table[,7:1006])))
)

pdf('ROI-TMS/Dist_perm_acc_TMS.pdf',15,12)
ggplot(perm_acc_df) +
  geom_histogram(aes(Acc, fill = TMS),alpha = 0.5, position = "identity", bins = 20) +
  facet_wrap(~Sub)
dev.off()

# do bootstrap for group-level null distribution
n_repetitions <- 1e4 # (this high number is needed; 1e3 gave unstable dists)
matrix_sham = use_dat_table[use_dat_table$TMS=='sham',7:1006]
matrix_cTBS = use_dat_table[use_dat_table$TMS=='cTBS',7:1006]

averages <- matrix(0,n_repetitions,2)
for (i in 1:n_repetitions) {
  if (i %% 100 == 0) cat("Iteration", i, "\n")
  averages[i,1] <- mean(apply(matrix_sham, 1, sample, size = 1))
  averages[i,2] <- mean(apply(matrix_cTBS, 1, sample, size = 1))
}

Obs_Acc_P = dat_table$Acc[dat_table$TMS=='sham']
Obs_Acc_A = dat_table$Acc[dat_table$TMS=='cTBS']

sum(averages[,1] >= mean(Obs_Acc_P)) / n_repetitions
sum(averages[,2] >= mean(Obs_Acc_A)) / n_repetitions

hist(averages[,1],freq = F,col = alpha(UseColor[2],0.5),
     xlab = '',main = 'Decoding Acc from permutation & bootstrap')
hist(averages[,2],add=T,freq = F,col = alpha(UseColor[1],0.5))

abline(v=mean(Obs_Acc_A),col=UseColor[1],lwd=2)
abline(v=mean(Obs_Acc_P),col=UseColor[2],lwd=2)

##################################
# check TMS effect across rev trials
##################################

modelname_list = c('Expect_pseudoconcat_Nz_rev_trial-1',
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
TMSconds = c('A','P')

#########
## this block can be skipped when coming back

dat_table = expand.grid(subdirs,modelname_list,mask_labels,
                        dat_scale_strings,TMSconds)
names(dat_table) = c('Sub','Model','Mask','Scale','TMS')
dat_table$Acc = NA
dat_table[1,7:1006] = NA

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
          permutation_acc = read.csv(file.path(tmpdir,'PermutationAcc.csv'),header = F)
          find_row = dat_table$Sub==subdirs[i] & 
            dat_table$Model==modelname_list[j] &
            dat_table$Mask==mask_labels[k] &
            dat_table$Scale==dat_scale_strings[m] &
            dat_table$TMS==t
          dat_table$Acc[find_row] = as.numeric(acc)
          dat_table[find_row,7:1006] = as.numeric(permutation_acc$V1)
        }}}}}

dat_table$TMS = plyr::mapvalues(dat_table$TMS,from = c('P','A'),to = c('sham','cTBS'))

# save this dat_table
write.csv(dat_table, 'ROI-TMS/Decoding_ROI_TMS_Expect_Acc_rev_trials.csv', row.names=FALSE)

#########
## plotting can start from here
dat_table = read.csv('ROI-TMS/Decoding_ROI_TMS_Expect_Acc_rev_trials.csv')
dat_table$Mask = factor(dat_table$Mask, levels = mask_labels) # specify mask orders
dat_table$Scale = factor(dat_table$Scale, levels = dat_scale_strings)

pdf('ROI-TMS/CompareDecoding_acrossTrials.pdf',12,10)
fig = dat_table %>%
      subset(Model %in% modelname_list & 
               #Mask == 'lOFC_cTBS' &
               Scale %in% c('Mean_across_conds','raw')) %>%
  ggplot(aes(x = Model, y = Acc, color = TMS)) + 
  geom_boxplot(alpha=0.5,outlier.alpha = 0.2,size=1.1) +
  facet_wrap(~Scale + Mask) +
  scale_colour_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  labs(y = "Decoding Accuracy (%)", x = '',title = '') + 
  common + theme(legend.position = "none" ) +
  scale_x_discrete(labels=c(expression(rev[-1]),"rev",expression(rev[+1]),
                            expression(rev[+2]))) +
  geom_hline(yintercept=100/3, linetype="dotted", color = "black", size=0.5)
print(fig)
dev.off()
  
  
