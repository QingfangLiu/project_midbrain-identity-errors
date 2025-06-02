
library(tidyverse)
library(readxl)
library(lme4)
rm(list = ls())

source("~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/figs/ggplotSetup.R")

dat = R.matlab::readMat('GC_2bins_dat_partial.mat')$dat
ROIs = c('OFC','LPFC')
nROIs = length(ROIs)
nSess = 2
nruns = 3
ntimes = nruns * 2

# this is after excluding a few subjects
SubInfo = read_excel('../../SubjectConds.XLSX')
UseTMSConds = SubInfo$Cond[SubInfo$Excluded==0]
Subs = SubInfo$Sub[SubInfo$Excluded==0]
nSubs = length(Subs)

# convert this array to df
dfGC = structure(dat,.Dim = c(nROIs, nSubs, nSess, ntimes), 
                  .Dimnames = structure(list(ROI = ROIs, 
                                             Sub = Subs,
                                             Sess = c(1,2),
                                             Time = 1:6), 
                  .Names = c("ROI", "Sub", "Sess", "Time")))
dfGC = plyr::adply(dfGC, c(1,2,3,4))
names(dfGC) = c(names(dfGC)[1:4],'GC')

# add TMS column
dfGC$Time = varhandle::unfactor(dfGC$Time) # unfactor time
dfGC$TMS = varhandle::unfactor(dfGC$Sess)
for(j in 1:nSubs){
  if(UseTMSConds[j]=='AP'){
    dfGC$TMS[dfGC$Sub==Subs[j]] = 3 - dfGC$TMS[dfGC$Sub==Subs[j]]
  }
}

dfGC$TMS = plyr::mapvalues(dfGC$TMS,from = c(1,2),to = c('sham','cTBS'))
dfGC$TMS = factor(dfGC$TMS)

## plotting

p = NULL
for(i in 1:length(ROIs)){
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
  common +  theme(plot.title = element_text(hjust = 0.5)) +
  #labs(x ="", y = "Global Connectedness") +
  labs(x ="", y = "") +
  theme(legend.position = "none",
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank()) +
  geom_vline(xintercept=c(2.5,4.5), linetype="dotted", color = "black", size=0.5)
}

p[[1]] = p[[1]] +
  coord_cartesian(ylim = c(0.095, 0.115), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") +
  #labs(title = 'Targeted OFC seed') +
  labs(title = '') +
  #annotate(geom = "text", x = c(1.6,3.5,5.4), y = 0.094, 
  #         label = c("Run 1", "Run 2", "Run 3"), size = 6) +
  annotate(geom = "text", x = c(0.9,1.1,1.9,2.1), y = c(0.103,0.103,0.104,0.104), 
           label = '*', size = 6)

p[[2]] = p[[2]] +
  coord_cartesian(ylim = c(0.109, 0.128), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") +
  #labs(title = 'LPFC stimulation site') +
  labs(title = '') +
  #annotate(geom = "text", x = c(1.6,3.5,5.4), y = 0.108, 
  #         label = c("Run 1", "Run 2", "Run 3"), size = 6) +
  annotate(geom = "text", x = c(1.9,2.1,3,4), y = c(0.117,0.117,0.119,0.12), 
           label = '*', size = 6)

pdf('GC_2bins_partial.pdf',5,8)
fig = ggpubr::ggarrange(p[[1]], p[[2]], ncol = 1, nrow = 2, align = "hv")
print(fig)
dev.off()

###### stat testing

for(i in 1:length(ROIs)){
  #####################
  # test linear trend for each condition
  for(t in c('sham','cTBS')){
    print(t)
    print('Simple linear regression')
    model <- lm(GC ~ Time, data = subset(dfGC, ROI == ROIs[i] & TMS == t)) 
    print(summary(model))
    print('linear mixed with time')
    lme1 = lmer(GC ~ Time + Sub + (1 | Sub), data = subset(dfGC, ROI == ROIs[i] & TMS == t))
    print(summary(lme1))
    print('linear mixed without time')
    lme2 = lmer(GC ~ Sub + (1 | Sub), data = subset(dfGC, ROI == ROIs[i] & TMS == t))
    print(summary(lme2))
    print('Compare LME models')
    print(anova(lme1,lme2))
  }
  
  # test the difference of TMS conditions
  model1 <- lmer(GC ~ Time * TMS + Sub + (1 | Sub), data = subset(dfGC, ROI == ROIs[i]))
  model2 <- lmer(GC ~ Time + TMS + Sub + (1 | Sub), data = subset(dfGC, ROI == ROIs[i]))
  summary(model)
  anova(model1,model2)
  
}

# test TMS effect at each time point for each ROI

for(i in 1:2){
for(t in 1:6){
  print(paste('Time',t,'for',ROIs[i]))
  test_data = subset(dfGC, ROI == ROIs[i] & Time == t)
  test_data = test_data[order(test_data$Sub),] # make sure the order is correct!
  
  test_sham = subset(test_data,TMS == 'sham')$GC
  test_cTBS = subset(test_data,TMS == 'cTBS')$GC
  
  print(t.test(test_cTBS,test_sham,paired = T,alternative = 'l'))
  print(wilcox.test(test_cTBS,test_sham,paired = T,alternative = 'l'))
  
  test_Sess1 = subset(test_data,Sess == 1)$GC
  test_Sess2 = subset(test_data,Sess == 2)$GC
  print(wilcox.test(test_Sess1,test_Sess2,paired = T))
  
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

