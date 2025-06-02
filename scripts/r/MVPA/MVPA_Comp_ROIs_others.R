
# This script conducts testing & plotting on ROI-based MVPA 
# also saving the summarized data into xlsx
# test other ROIs

rm(list = ls())

# define directory paths
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
AnaDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis'

MVPADir = file.path(ResDir,'MVPA_ROI')
SubInfoDir = file.path(AnaDir,'SubInfo')

# ggplot setup
SetupFile = file.path(AnaDir,'Scripts_R','Setup.R')
source(SetupFile)

##################################################################

SubInfo = read.xlsx(file.path(SubInfoDir,'SubjectConds.xlsx'))
UseTMSConds = SubInfo$Cond[SubInfo$Excluded==0]
Subs = SubInfo$Sub[SubInfo$Excluded==0]
nsubs = length(Subs)

ROIs = c('Func_DS','Func_VS',
         'Func_MB','Func_mPFC',
         'Func_Insula','Func_HPC',
         'Func_Thalamus','Func_Acc')    
nROIs = length(ROIs)

ROIs.labs = c('Dorsal Striatum','Ventral Striatum','Midbrain',
         'mPFC','Insula','Hippocampus','Thalamus','Acc')    

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
  labs(y = "", title = "", x = '') + 
  common +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.text.x = element_blank(),legend.position = "none" ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted")
}

pdf(file.path(MVPADir,'MVPA_additionalROIs.pdf'),14,8)
fig = ggpubr::ggarrange(p[[1]], p[[2]], 
                        p[[3]], p[[4]], 
                        p[[5]], p[[6]],
                        p[[7]], p[[8]],
                        ncol = 4, nrow = 2, align = "hv")
print(fig)
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
  print(wilcox.test(x_sham,x_cTBS,paired = T,alternative = 't'))
  print(t.test(x_sham,x_cTBS,paired = T,alternative = 't'))
}



