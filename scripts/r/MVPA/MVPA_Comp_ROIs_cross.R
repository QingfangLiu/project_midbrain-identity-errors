
# This script conducts testing & plotting on ROI-based MVPA 
# also saving the summarized data into xlsx

rm(list = ls())

# define directory paths
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
MVPADir = file.path(ResDir,'MVPA_ROI_cross')
SubInfoDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/SubInfo'

# ggplot setup
source("../Setup.R")

##################################################################

SubInfo = read.xlsx(file.path(SubInfoDir,'SubjectConds.xlsx'))
UseTMSConds = SubInfo$Cond[SubInfo$Excluded==0]
Subs = SubInfo$Sub[SubInfo$Excluded==0]
nsubs = length(Subs)

ROIs = c('Func_OFC','Func_LPFC')    # ROI names in savd results
nROIs = length(ROIs)

All_corr_vals_df = NULL

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

#pdf(file.path(MVPADir,'MVPSA_OFC.pdf'),3.5,4)
par(mar=c(1,1,1,1))
print(p[[1]])
#dev.off()

high_y2 = 0.18
low_y2 = 0.17
p[[2]] = p[[2]] +
  geom_segment(aes(x=1,y=high_y2,xend=2,yend=high_y2),linewidth = 0.2) +
  geom_segment(aes(x=1,y=low_y2,xend=1,yend=high_y2),linewidth = 0.2) +
  geom_segment(aes(x=2,y=low_y2,xend=2,yend=high_y2),linewidth = 0.2) +
  annotate(geom = "text", x = c(1.5), y = c(0.2), label = 'n.s.', size = 6)

#pdf(file.path(MVPADir,'MVPSA_LPFC.pdf'),3.5,4)
par(mar=c(1,1,1,1))
print(p[[2]])
#dev.off()

for(roi in ROIs){
  print(roi)
  x_sham = subset(Corr_diff_vals_df_summary,
                  ROI == roi & TMS == 'sham')$Representation
  x_cTBS = subset(Corr_diff_vals_df_summary,
                  ROI == roi & TMS == 'cTBS')$Representation
  
  x = subset(Corr_diff_vals_df_summary,ROI == roi)$Representation
  print(wilcox.test(x)) # representation > 0 from cTBS trials?
  print(wilcox.test(x_sham,x_cTBS,paired = T,alternative = 't'))
}


