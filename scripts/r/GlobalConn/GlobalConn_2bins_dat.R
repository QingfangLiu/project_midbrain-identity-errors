
rm(list = ls())

# define directory paths
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
AnaDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis'

# dependent directory paths
GCDir = file.path(ResDir,'GlobalConn_ROI')
SubInfoDir = file.path(AnaDir,'SubInfo')

# ggplot setup
SetupFile = file.path(AnaDir,'Scripts_R','Setup.R')
source(SetupFile)

# source data 
sourceDir = '~/Desktop/SourceData'

##################################################################

nSess = 2
nruns = 3
ntimes = nruns * 2

SubInfo = read.xlsx(file.path(SubInfoDir,'SubjectConds.xlsx'))
UseTMSConds = SubInfo$Cond[SubInfo$Excluded==0]
Subs = SubInfo$Sub[SubInfo$Excluded==0]
nSubs = length(Subs)

# load GC ROI-based data
dat_file = file.path(GCDir,'GC_2bins_dat.mat') 
dat = readMat(dat_file)$dat
ROIs = readMat(dat_file)$ROIlabels 
ROIs = unlist(ROIs)
nROIs = length(ROIs)

# convert dat array to df
dfGC = structure(dat,.Dim = c(nROIs, nSubs, nSess, ntimes), 
                  .Dimnames = structure(list(ROI = ROIs, 
                                             Sub = Subs,
                                             Sess = c(1,2),
                                             Time = 1:6), 
                  .Names = c("ROI", "Sub", "Sess", "Time")))
dfGC = adply(dfGC, c(1,2,3,4))
names(dfGC) = c(names(dfGC)[1:4],'GC')

# add TMS column
dfGC$Time = unfactor(dfGC$Time) # unfactor time
dfGC$TMS = unfactor(dfGC$Sess)
for(j in 1:nSubs){
  if(UseTMSConds[j]=='AP'){
    dfGC$TMS[dfGC$Sub==Subs[j]] = 3 - dfGC$TMS[dfGC$Sub==Subs[j]]
  }
}

dfGC$TMS = mapvalues(dfGC$TMS,from = c(1,2),to = c('sham','cTBS'))
dfGC$TMS = factor(dfGC$TMS)


## save df into excel file
dfGC = dfGC %>%
  subset(ROI %in% c('OFC','LPFC','func_OFC','func_LPFC',
                    'Insula_b','MB_b','mPFC','Striatum_b',
                    'AMG_b','Thalamus'))

write.xlsx(dfGC,file = file.path(sourceDir,'GlobalConnectedness.xlsx'))

## plotting

p = NULL # for line, point, ribbon plots
for(i in 1:2){
  tmp_dfGC = subset(dfGC, ROI == ROIs[i]) 
  avg_df = tmp_dfGC %>%
           group_by(Time, TMS) %>% 
           dplyr::summarise(mean = mean(GC), sd = sd(GC),
                   n = n(), sde = sd/sqrt(n))
  
p[[i]] = ggplot(avg_df,aes(Time)) + # notice here 'Time' can't be a factor
  geom_point(aes(y = mean, color = TMS), size = 3) +
  geom_line(aes(y = mean, color = TMS), linewidth = 1.5) +
  geom_ribbon(aes(ymin = mean - sde, ymax = mean + sde, fill = TMS), alpha = 0.3) +
  scale_color_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  common +
  labs(x ="", y = "") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  geom_vline(xintercept=c(2.5,4.5), linetype="dotted", color = "black", size=0.5)

}

p[[1]] = p[[1]] +
  labs(title = 'Targeted OFC seed') +
  coord_cartesian(ylim = c(0.095, 0.115), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") +
  annotate(geom = "text", x = c(0.9,1.1,1.9,2.1), y = c(0.103,0.103,0.104,0.104), 
           label = '*', size = 6)

p[[2]] = p[[2]] +
  labs(title = 'LPFC stimulation site') +
  coord_cartesian(ylim = c(0.109, 0.128), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") +
  annotate(geom = "text", x = c(1.9,2.1,3,4), y = c(0.117,0.117,0.119,0.12), 
           label = '*', size = 6)


pdf(file.path(GCDir,'GC_2bins.pdf'),5,8)
fig = ggarrange(p[[1]], p[[2]], ncol = 1, nrow = 2)
print(fig)
dev.off()

###### stat testing on OFC & LPFC ROIs

for(r in c('OFC','LPFC')){
  for(t in c('sham','cTBS')){ # test linear trend for each TMS condition
    print(t)
    print('Simple linear regression')
    model <- lm(GC ~ Time, data = subset(dfGC, ROI == r & TMS == t)) 
    print(summary(model))
    print('linear mixed with time')
    lme1 = lmer(GC ~ Time + Sub + (1 | Sub), data = subset(dfGC, ROI == r & TMS == t))
    print(summary(lme1))
    print('linear mixed without time')
    lme2 = lmer(GC ~ Sub + (1 | Sub), data = subset(dfGC, ROI == r & TMS == t))
    print(summary(lme2))
    print('Compare LME models')
    print(anova(lme1,lme2))
  }
  
  # test the difference of TMS conditions
  model1 <- lmer(GC ~ Time * TMS + Sub + (1 | Sub), data = subset(dfGC, ROI == r))
  model2 <- lmer(GC ~ Time + TMS + Sub + (1 | Sub), data = subset(dfGC, ROI == r))
  summary(model)
  anova(model1,model2)
  
}

# test TMS effect at each time point for each ROI

for(r in c('OFC','LPFC')){
for(t in 1:ntimes){
  print(paste('Time',t,'for',r))
  test_data = subset(dfGC, ROI == r & Time == t)
  test_data = test_data[order(test_data$Sub),] # make sure the order is correct!
  
  test_sham = subset(test_data,TMS == 'sham')$GC
  test_cTBS = subset(test_data,TMS == 'cTBS')$GC
  print(wilcox.test(test_cTBS,test_sham,paired = T,alternative = 'l'))
  
  lme1 = lmer(GC ~ TMS + (1 | Sub),data = test_data)
  lme2 = lmer(GC ~ Sess + (1 | Sub),data = test_data)
  lme3 = lmer(GC ~ TMS + Sess + (1 | Sub),data = test_data)
  anova(lme1,lme2)
  anova(lme1,lme3) # Sess is useful?
  anova(lme2,lme3) # TMS is useful?
              
}
}

# conclusion here from LME model comparison
# for OFC
# run1,2: effect of TMS after accounting for Sess
# run3,4: effect of Sess after accounting for TMS
# run5,6: none

# for LPFC
# run1,5: none
# run2: effect of TMS after accounting for Sess
# run3,4: effect of both TMS and Sess
# run6: effect of Sess after accounting for TMS


##################################################################
# reorganize by putting two ROIs in separate columns
# convenient to be correlated with other experiment variables

vecs = matrix(NA,372,2)
for(i in 1:2){
  reduced = dfGC %>%
    select(ROI,Sub,TMS,Time,GC) %>%
    subset(ROI == ifelse(i==1,'OFC','LPFC')) %>%
    arrange(Sub,TMS,Time) 
  vecs[,i] = reduced$GC
}

re_df = data.frame(vecs)
names(re_df) = c('OFC','LPFC')
re_df$Sub = reduced$Sub
re_df$TMS = reduced$TMS
re_df$Time = reduced$Time
write.xlsx(re_df,file = file.path(GCDir,'GC_SubsTMS.xlsx'),rowNames = F)


##################################################################
# look at two functional ROIs
##################################################################

for(i in 3:4){
  for(t in 1:ntimes){
    print(paste('Time',t,'for',ROIs[i]))
    test_data = subset(dfGC, ROI == ROIs[i] & Time == t)
    test_data = test_data[order(test_data$Sub),] # make sure the order is correct!
    
    test_sham = subset(test_data,TMS == 'sham')$GC
    test_cTBS = subset(test_data,TMS == 'cTBS')$GC
    print(wilcox.test(test_cTBS,test_sham,paired = T,alternative = 'l'))
  }}

# for functional OFC
# time 1: p = 9.77e-5
# time 2: p = 0.00012

# for functional LPFC
# time 1: p = 0.0027
# time 2: p = 3.89e-5
# time 4: p = 0.045


### cb plots

cb = NULL
for(i in 3:4){
  tmp_dfGC = subset(dfGC, ROI == ROIs[i]) 
  avg_df = tmp_dfGC %>%
    group_by(Time, TMS) %>% 
    dplyr::summarise(mean = mean(GC), sd = sd(GC),
                     n = n(), sde = sd/sqrt(n))
  
  cb[[i]] = ggplot(tmp_dfGC,aes(Time)) + # notice here 'Time' can't be a factor
    geom_line(aes(y = GC, group = interaction(TMS, Sub), color = TMS), 
              linewidth = 0.3, linetype = 'solid',show.legend = F) +
    geom_point(data = avg_df,aes(y = mean, color = TMS), size = 3) +
    geom_line(data = avg_df,aes(y = mean, color = TMS), linewidth = 1.5,show.legend = F) +
    geom_ribbon(data = avg_df,aes(ymin = mean - sde, ymax = mean + sde, fill = TMS), alpha = 0.3) +
    scale_color_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
    scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
    common +  
    labs(x ="", y = "Global Connectedness \n (z-scored)") +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.ticks.x=element_blank(),
          axis.text.x = element_blank()) +
    geom_vline(xintercept=c(2.5,4.5), linetype="dotted", color = "black", size=0.5)
  
}


cb[[3]] = cb[[3]] +
  coord_cartesian(ylim = c(0.068, 0.185), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") +
  labs(title = 'OFC (functional ROI)') +
  annotate(geom = "text", x = c(1.6,3.5,5.4), y = 0.061, 
           label = c("Run 1", "Run 2", "Run 3"), size = 6) +
  annotate(geom = "text", x = c(0.9,1.1,1.9,2.1), y = c(0.11,0.11,0.11,0.11), 
           label = '*', size = 6)

cb[[4]] = cb[[4]] +
  coord_cartesian(ylim = c(0.07, 0.225), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") +
  labs(title = 'LPFC (functional ROI)') +
  annotate(geom = "text", x = c(1.6,3.5,5.4), y = 0.06, 
           label = c("Run 1", "Run 2", "Run 3"), size = 6) +
  annotate(geom = "text", x = c(0.9,1.1,1.9,2.1,4), y = c(0.105,0.105,0.108,0.108,0.111), 
           label = '*', size = 6)

fig = ggarrange(cb[[3]], cb[[4]], 
                common.legend = T,
                legend = 'right',
                ncol = 2, nrow = 1, align = "hv")
pdf(NULL)
pdf(file.path(GCDir,'GC_2bins_funcROIs_overlay.pdf'),10,4)
print(fig)
dev.off()


##################################################################
### another plot to look at difference value distributions

bplt = NULL
for(i in 3:4){
  tmp_dfGC = subset(dfGC, ROI == ROIs[i]) 
  tmp_dfGC_cTBS = tmp_dfGC %>%
                  subset(TMS == 'cTBS') %>%
                  arrange(Sub,Time)
  tmp_dfGC_sham = tmp_dfGC %>%
                  subset(TMS == 'sham') %>%
                  arrange(Sub,Time)
  GC_diff = tmp_dfGC_cTBS$GC - tmp_dfGC_sham$GC
  tmp_dfGC_diff = cbind(tmp_dfGC_cTBS,GC_diff)
  tmp_dfGC_diff$GC = NULL
  avg_df = tmp_dfGC_diff %>%
    group_by(Time) %>% 
    dplyr::summarise(mean = mean(GC_diff), sd = sd(GC_diff),
                     n = n(), sde = sd/sqrt(n))
  bplt[[i]] = ggplot(tmp_dfGC_diff,aes(x=factor(Time),y=GC_diff)) + # notice here 'Time' can't be a factor
    geom_line(aes(y = GC_diff, group = Sub), 
              linewidth = 0.3, linetype = 'solid',
              show.legend = F, color = 'gray') +
    geom_point(data = avg_df,aes(y = mean), 
               size = 3,color = 'black') +
    geom_line(data = avg_df,aes(x = Time,y = mean), linewidth = 1.5,
              show.legend = F,color = 'black') +
    geom_ribbon(data = avg_df,aes(x = Time, y = mean, ymin = mean - sde, 
                                  ymax = mean + sde), alpha = 0.3) +
    common +
    labs(x ="", y = "") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.ticks.x=element_blank(),
          axis.text.x = element_blank()) +
    geom_vline(xintercept=c(2.5,4.5), linetype="dotted", color = "black", size=0.5) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black", size=1)
  
}

bplt[[3]] = bplt[[3]] +
  labs(title = 'OFC (functional ROI)') +
  coord_cartesian(ylim = c(-0.085, 0.06), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off")
bplt[[4]] = bplt[[4]] +
  labs(title = 'LPFC (functional ROI)') +
  coord_cartesian(ylim = c(-0.14, 0.07), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off")
fig = ggarrange(bplt[[3]], bplt[[4]], ncol = 1, nrow = 2)
pdf(NULL)
pdf(file.path(GCDir,'GC_2bins_funcROIs_diff.pdf'),5,8)
print(fig)
dev.off()


##################################################################
# look at other ROIs
# insula, midbrain, mPFC, striatum, thalamus, amygdala
##################################################################

library(rstatix)
ridx = c(7,10,11,14,18,17)
ROIs[ridx]
# perform linear mixed effect models for each ROI

for(i in ridx){
  print(ROIs[i])
  test_data = subset(dfGC, ROI == ROIs[i])
  res.aov <- anova_test(
    data = test_data, dv = GC, wid = Sub,
    within = c(TMS, Time)
  )
  print(get_anova_table(res.aov))
  
}
  
  # # anova (not working well)
  # anova_result = aov(GC ~ TMS * Time + Error(Sub/(TMS * Time)), data = test_data)
  # print(summary(anova_result))
  
  # test the effect of time
  lme1 = lmer(GC ~ TMS * Time + Sub + (1 | Sub), 
              data = test_data, control = lmerControl(optimizer ="Nelder_Mead"))
  summary(lme1)
  lme2 = lmer(GC ~ TMS + Time + Sub + (1 | Sub), 
              data = test_data, control = lmerControl(optimizer ="Nelder_Mead"))
  lme3 = lmer(GC ~ TMS + Sub + (1 | Sub), 
              data = test_data, control = lmerControl(optimizer ="Nelder_Mead"))
  lme4 = lmer(GC ~ Time + Sub + (1 | Sub), 
              data = test_data, control = lmerControl(optimizer ="Nelder_Mead"))
  lme5 = lmer(GC ~ Sub + (1 | Sub), 
              data = test_data, control = lmerControl(optimizer ="Nelder_Mead"))
  
  print('Compare LME models: interaction effect?')
  print(anova(lme1,lme2))
  
  print('Compare LME models: effects of time & TMS?')
  print(anova(lme2,lme5))
  
  print('Compare LME models: effect of time beyond TMS?')
  print(anova(lme2,lme3))
  
  print('Compare LME models: effect of TMS beyond time?')
  print(anova(lme2,lme4))
  



for(i in ridx){
  for(t in 1:ntimes){
    print(paste('Time',t,'for',ROIs[i]))
    test_data = subset(dfGC, ROI == ROIs[i] & Time == t)
    test_data = test_data[order(test_data$Sub),] # make sure the order is correct!
    
    test_sham = subset(test_data,TMS == 'sham')$GC
    test_cTBS = subset(test_data,TMS == 'cTBS')$GC
    print(wilcox.test(test_cTBS,test_sham,paired = T,alternative = 'l'))
  }}

# for insula: ns for all times

# for midbrain:
# time 2: p = 0.0043
# time 6: p = 0.045

# for mPFC: ns for all times

# for striatum: 
# time 1: p = 0.022

# for thalamus: ns for all times

# for AMG:
# time 1: p = 0.014
# time 2: p = 0.021



for(i in ridx){
  tmp_dfGC = subset(dfGC, ROI == ROIs[i]) 
  avg_df = tmp_dfGC %>%
    group_by(Time, TMS) %>% 
    dplyr::summarise(mean = mean(GC), sd = sd(GC),
                     n = n(), sde = sd/sqrt(n))
  
  p[[i]] = ggplot(avg_df,aes(Time)) + # notice here 'Time' can't be a factor
    geom_point(aes(y = mean, color = TMS), size = 3) +
    geom_line(aes(y = mean, color = TMS), linewidth = 1.5) +
    geom_ribbon(aes(ymin = mean - sde, ymax = mean + sde, fill = TMS), alpha = 0.3) +
    scale_color_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
    scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
    common +
    labs(x ="", y = "") +
    #labs(x ="", y = "Global Connectedness \n (z-scored)") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.ticks.x=element_blank(),
          axis.text.x = element_blank()) +
    geom_vline(xintercept=c(2.5,4.5), linetype="dotted", color = "black", size=0.5)
}

p[[ridx[1]]] = p[[ridx[1]]] +
  coord_cartesian(ylim = c(0.105, 0.12), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") 
  #labs(title = 'Insula') +
  #annotate(geom = "text", x = c(1.6,3.5,5.4), y = 0.104, 
  #         label = c("Run 1", "Run 2", "Run 3"), size = 6)

p[[ridx[2]]] = p[[ridx[2]]] +
  coord_cartesian(ylim = c(0.075, 0.095), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") 
  #labs(title = 'Midbrain') +
  #annotate(geom = "text", x = c(1.6,3.5,5.4), y = 0.0738, 
  #         label = c("Run 1", "Run 2", "Run 3"), size = 6)

p[[ridx[3]]] = p[[ridx[3]]] +
  coord_cartesian(ylim = c(0.113, 0.130), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") 
  #labs(title = 'mPFC') +
  #annotate(geom = "text", x = c(1.6,3.5,5.4), y = 0.112, 
  #         label = c("Run 1", "Run 2", "Run 3"), size = 6)

p[[ridx[4]]] = p[[ridx[4]]] +
  coord_cartesian(ylim = c(0.082, 0.105), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") 
  # labs(title = 'Striatum') +
  # annotate(geom = "text", x = c(1.6,3.5,5.4), y = 0.0805, 
  #          label = c("Run 1", "Run 2", "Run 3"), size = 6)

p[[ridx[5]]] = p[[ridx[5]]] +
  coord_cartesian(ylim = c(0.081, 0.102), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off")
  #labs(title = 'Thalamus') +
  # annotate(geom = "text", x = c(1.6,3.5,5.4), y = 0.0795, 
  #          label = c("Run 1", "Run 2", "Run 3"), size = 6)

p[[ridx[6]]] = p[[ridx[6]]] +
  coord_cartesian(ylim = c(0.0995, 0.111), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off")
  #labs(title = 'Amygdala') +
  # annotate(geom = "text", x = c(1.6,3.5,5.4), y = 0.0988, 
  #          label = c("Run 1", "Run 2", "Run 3"), size = 6)

pdf(file.path(GCDir,'GC_2bins_otherROIs.pdf'),15,8)
fig = ggpubr::ggarrange(p[[ridx[1]]], p[[ridx[2]]], 
                        p[[ridx[3]]], p[[ridx[4]]], 
                        p[[ridx[5]]], p[[ridx[6]]],
                        ncol = 3, nrow = 2, align = "hv")
print(fig)
dev.off()


