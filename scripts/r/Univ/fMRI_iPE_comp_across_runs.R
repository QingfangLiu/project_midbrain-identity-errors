
rm(list = ls())
library(varhandle)
library(readxl)
library(plyr)

# ggplot setup
source("~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/figs/ggplotSetup.R")

# this is after excluding a few subjects
SubInfo = read_excel('SubjectConds.XLSX')
UseTMSConds = SubInfo$Cond[SubInfo$Excluded==0]
Subs = SubInfo$Sub[SubInfo$Excluded==0]
nsubs = length(Subs)

use_dir = "~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/NeuralAnalyzeRes"

All_beta_vals_df = NULL
for(j in 1:nsubs){

subno = paste0('Sub',Subs[j])
tmp_dir = file.path(use_dir,subno)
matdat = R.matlab::readMat(file.path(tmp_dir,'betas.mat'))
ROIs = unlist(matdat$roi.labels)
nROIs = length(ROIs)
sub_beta_vals = matdat$beta.vals
sub_beta_vals_df = structure(sub_beta_vals,
                       .Dim = c(nROIs, 2, 3, 4), 
                       .Dimnames = structure(list(ROI = ROIs, 
                                                  Sess = 1:2, 
                                                  Run = 1:3,
                                                  bidx = 1:4), 
                       .Names = c("ROI" ,"Sess", "Run", "bidx")))
sub_beta_vals_df = plyr::adply(sub_beta_vals_df, c(1,2,3,4))
names(sub_beta_vals_df) = c(names(sub_beta_vals_df)[1:4],'vals')
sub_beta_vals_df$Sub = subno
All_beta_vals_df[[j]] = sub_beta_vals_df

}

beta_vals_df <- do.call(rbind, All_beta_vals_df)
beta_vals_df$TMS = unfactor(beta_vals_df$Sess)
for(j in 1:nsubs){
  if(UseTMSConds[j]=='AP'){
    beta_vals_df$TMS[beta_vals_df$Sub==Subs[j]] = 3 - 
      beta_vals_df$TMS[beta_vals_df$Sub==Subs[j]]
  }
}

beta_vals_df$TMS = plyr::mapvalues(beta_vals_df$TMS,from = c(1,2),to = c('sham','cTBS'))
beta_vals_df$TMS = factor(beta_vals_df$TMS)

##
# reduce the matrix by bt12 - bt23
beta_diff_vals = subset(beta_vals_df,bidx == 3)$vals -
  subset(beta_vals_df,bidx == 1)$vals

beta_diff_vals_df = subset(beta_vals_df,bidx == 3)
beta_diff_vals_df$vals = beta_diff_vals
beta_diff_vals_df$bidx = 'post-rev'

roi = ROIs
ggplot(data = subset(beta_diff_vals_df,ROI %in% roi),
       mapping = aes(x = TMS, y = vals, fill = TMS)) + 
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  facet_wrap(~ Run)
  
beta_diff_vals_df %>%
    group_by(TMS,Run) %>%
    dplyr::summarise(
      iPE = mean(vals)
    )

beta_diff_vals_df_summary = beta_diff_vals_df %>%
  group_by(TMS,Sub,Run,ROI) %>%
  dplyr::summarise(
    iPE = mean(vals)
  ) %>% 
  arrange(Sub) %>%
  ungroup() %>%
  dplyr::mutate(paired = rep(1:(n()/12),each=12))
  

################ plots ###################

ggplot(data = beta_diff_vals_df_summary,
              aes(x = TMS, y = iPE, fill = TMS)) + 
  facet_wrap(~ ROI + Run) +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=TMS,group=paired),size=2,shape=21,alpha=1,
             position = position_dodge(0.2)) +
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  labs(y = 'iPE change', x = '') + 
  common +
  theme(plot.title = element_text(hjust = 0.5)) +
  #theme(strip.text.x = element_blank()) +
  theme(legend.position = "none" ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted")

  write.csv(beta_diff_vals_df_summary,
              file = './NeuralAnalyzeRes/beta_diff_vals_df_summary.csv')
  
  # t-test or wilcoxon test says only the effect of MB at run2 was sig
  # so don't think it's very helpful to include this information in the paper
  test_dat_sham = subset(beta_diff_vals_df_summary, ROI == 'MB' & Run == 3 & TMS == 'sham')$iPE
  test_dat_cTBS = subset(beta_diff_vals_df_summary, ROI == 'MB' & Run == 3 & TMS == 'cTBS')$iPE
  wilcox.test(test_dat_cTBS,test_dat_sham,paired = T,alternative = 'g')
  
  

