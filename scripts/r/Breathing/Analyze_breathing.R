

rm(list = ls())

library(plyr)
library(tidyverse)
library(R.matlab)
library(lme4)
library(readxl)
# here to load plyr first because we want to use dplyr::summarise 

# load subject info
datinfo = read_excel("SubjectConds.XLSX")
subs = subset(datinfo,Excluded==0)$Sub # exclude 4 subs based on datinfo
Conds = subset(datinfo,Excluded==0)$Cond
nsubs = length(subs)

# read breathing data from each subject
studydir = '/Users/qingfangliu/Experiment'
AllBreathing = array(NA,c(nsubs,6,64,7001))
ATMS_Breathing = array(NA,c(nsubs,3,64,7001))
PTMS_Breathing = array(NA,c(nsubs,3,64,7001))

# 6 runs in time order

for(i in 1:nsubs){
  
  subno = paste('Sub',subs[i],sep = '')
  filename = paste('Agg_resp_data_',subno,'.mat',sep = '')
  tmp = readMat(paste(studydir,subno,'ProcessedLC',filename,sep = '/'))
  # put the lists into arrays
  AllBreathing[i,1,,] = tmp$agg.resp.data[[1]][[1]]
  AllBreathing[i,2,,] = tmp$agg.resp.data[[3]][[1]]
  AllBreathing[i,3,,] = tmp$agg.resp.data[[5]][[1]]
  AllBreathing[i,4,,] = tmp$agg.resp.data[[2]][[1]]
  AllBreathing[i,5,,] = tmp$agg.resp.data[[4]][[1]]
  AllBreathing[i,6,,] = tmp$agg.resp.data[[6]][[1]]
  
  if (Conds[i]=='PA')  ATMS_Breathing[i,,,] = AllBreathing[i,4:6,,]
  if (Conds[i]=='PA')  PTMS_Breathing[i,,,] = AllBreathing[i,1:3,,]
  if (Conds[i]=='AP')  ATMS_Breathing[i,,,] = AllBreathing[i,1:3,,]
  if (Conds[i]=='AP')  PTMS_Breathing[i,,,] = AllBreathing[i,4:6,,]
}

# get colors to use for two TMS conditions
library(scales)
UseColor = hue_pal()(2)

# trace for each subject
SubATMS_Breathing = apply(ATMS_Breathing,c(1,4),mean,na.rm=T)
SubPTMS_Breathing = apply(PTMS_Breathing,c(1,4),mean,na.rm=T)
Sub_Breathing = apply(AllBreathing,c(1,4),mean,na.rm=T)

# baseline for each subject
# mean of 1s before cue 
ABase = apply(SubATMS_Breathing[,2001:3001],1,mean)
ABasemat = matrix(ABase,31,7001,byrow = F)
SubATMS_Breathing = SubATMS_Breathing - ABasemat

PBase = apply(SubPTMS_Breathing[,2001:3001],1,mean)
PBasemat = matrix(PBase,31,7001,byrow = F)
SubPTMS_Breathing = SubPTMS_Breathing - PBasemat

pdf('BreathingRes/Comp_breath_TMS.pdf',8,6)
par(mfrow=c(1,1),cex=1.5,cex.lab=1.5,mar=c(5,5,2,2))
plot(seq(-3,4,by=0.001),apply(SubATMS_Breathing,2,mean,na.rm=T),col= UseColor[1],
     xlab = 'Time from sniff cue (second)', ylab = 'Breathing Signals (a.u.)',
     xlim = c(-1,4),type = 'l',lwd = 3)
lines(seq(-3,4,by=0.001),apply(SubPTMS_Breathing,2,mean,na.rm=T),
      col= UseColor[2],lwd = 3)
abline(h=0, lty=3)
dev.off()

pdf('BreathingRes/Subs_Comp_breath_TMS.pdf',8,6)
par(mfrow=c(3,3))
for(i in 1:nsubs){
plot(seq(-3,4,by=0.001),SubATMS_Breathing[i,],col= UseColor[1],
     xlab = 'Time from sniff cue (second)', ylab = 'Breathing Signals (a.u.)',
     xlim = c(-1,4),type = 'l',lwd = 2,main = paste('Sub',subs[i],sep = ''))
lines(seq(-3,4,by=0.001),SubPTMS_Breathing[i,],
      col= UseColor[2],lwd = 2)
}
dev.off()

##### find individual peak amplitude 
PeakA = sapply(1:31,function(x) max(SubATMS_Breathing[x,]))
PeakP = sapply(1:31,function(x) max(SubPTMS_Breathing[x,]))

library(ggplot2)
data <- data.frame(
  TMS=c('Stim','Sham'),
  Amplitude=c(mean(PeakA),mean(PeakP)),
  sde=c(sd(PeakA),sd(PeakP))/sqrt(31)
)

p = NULL
p[[1]] = ggplot(data) +
  geom_bar(aes(x=TMS, y=Amplitude,fill=TMS), stat="identity") +
  geom_errorbar( aes(x=TMS, ymin=Amplitude-sde, ymax=Amplitude+sde), width=0.1, 
                 colour="black", alpha=0.9, size=1.3) +
  theme_classic(base_size = 22) +
  theme(legend.position = 'none')

FindDuration = function(x){
  out = which(x<1) # use 1 instead of 0 as threshold to accomodate sub30
  (out[out>4000][1] - 3000)/1000
}

DurationA = sapply(1:31,function(x) FindDuration(SubATMS_Breathing[x,]))
DurationP = sapply(1:31,function(x) FindDuration(SubPTMS_Breathing[x,]))

data <- data.frame(
  TMS=c('Stim','Sham'),
  Duration=c(mean(DurationA),mean(DurationP)),
  sde=c(sd(DurationA),sd(DurationP))/sqrt(31)
)

p[[2]] = ggplot(data) +
  geom_bar(aes(x=TMS, y=Duration,fill=TMS), stat="identity") +
  geom_errorbar( aes(x=TMS, ymin=Duration-sde, ymax=Duration+sde), width=0.1, 
                 colour="black", alpha=0.9, size=1.3) +
  theme_classic(base_size = 22) +
  theme(legend.position = 'none')

pdf('BreathingRes/Comp_breath_peak_duration.pdf',8,6)
library(gridExtra)
grid.arrange(p[[1]], p[[2]], nrow = 1)
dev.off()


