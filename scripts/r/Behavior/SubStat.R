
rm(list = ls())
library(openxlsx)
library(tidyverse)
library(lubridate)

# folders
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
behdir = file.path(ResDir,'Behavior')

SubInfoDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/SubInfo'

SubInfo = read.xlsx(file.path(SubInfoDir,'SubjectConds.xlsx'))
ExpInfo = read.xlsx(file.path(SubInfoDir,'SubMaskNumbers.xlsx'),sheet = 'Sheet3')

UseTMSConds = SubInfo$Cond[SubInfo$Excluded==0]
Subs = SubInfo$Sub[SubInfo$Excluded==0]
nSubs = length(Subs)

summary(SubInfo$Age)
sd(SubInfo$Age)
table(SubInfo$Sex)

dat = subset(SubInfo,Excluded==0)
summary(dat$Age)
sd(dat$Age)
table(dat$Sex)

# which finger movement was observed for rTMS (left or right)
# then you get the opposite side of brain for stimulation
table(dat$Finger_observed)

# age by TMS cond
mean(dat$Age[dat$Cond=='PA'])
mean(dat$Age[dat$Cond=='AP'])

sd(dat$Age[dat$Cond=='PA'])
sd(dat$Age[dat$Cond=='AP'])


par(mfrow=c(1,2))
hist(dat$MT,breaks = 8,main = '',xlab = 'Resting motor threshold (%)')
summary(dat$MT)
sd(dat$MT)

hist(dat$TMS,breaks = 8,main = '',xlab = 'TMS amplitude (%)')
summary(dat$TMS)
sd(dat$TMS)


hist(dat$SessGap,breaks = 20,main = '',xlab = 'Days between two TMS sessions')
summary(dat$SessGap)
sd(dat$SessGap)

summary(dat$Screen2MTGap)
sd(dat$Screen2MTGap)

summary(dat$MT2Sess1Gap)
sd(dat$MT2Sess1Gap)


ExpInfo = ExpInfo[SubInfo$Excluded==0,]
nsubs = nrow(ExpInfo)
TimeElapsed = array(NA,dim = c(nsubs,6))
TimeElapsed[,1] = minute(ExpInfo$TS1R1)
TimeElapsed[,2] = minute(ExpInfo$TS1R2)
TimeElapsed[,3] = minute(ExpInfo$TS1R3)
TimeElapsed[,4] = minute(ExpInfo$TS2R1)
TimeElapsed[,5] = minute(ExpInfo$TS2R2)
TimeElapsed[,6] = minute(ExpInfo$TS2R3)

# average time elapsed in the first run
FR_time = as.numeric(TimeElapsed[,c(1,4)])
summary(FR_time)
sd(FR_time)

TMSmat = matrix(NA,nsubs,6)
PASub = which(dat$Cond=='PA')
APSub = which(dat$Cond=='AP')

TMSmat[PASub,] = matrix(c(0,0,0,1,1,1),length(PASub),6,byrow = T)
TMSmat[APSub,] = matrix(c(1,1,1,0,0,0),length(APSub),6,byrow = T)

DatTimeElapsed = cbind(sort(rep(dat$Sub,6)),
      rep(c(1,1,1,2,2,2),nsubs),
      as.numeric(t(TMSmat)),
      rep(c(1,2,3,1,2,3),nsubs),
      as.numeric(t(TimeElapsed)))

DatTimeElapsed = as.data.frame(DatTimeElapsed)
names(DatTimeElapsed) = c('Sub','Session','TMS','Run','Time')
DatTimeElapsed$TMS[DatTimeElapsed$TMS==0] = 'Sham'
DatTimeElapsed$TMS[DatTimeElapsed$TMS==1] = 'cTBS'
DatTimeElapsed$Run = factor(DatTimeElapsed$Run)
DatTimeElapsed$TMS = factor(DatTimeElapsed$TMS)
DatTimeElapsed$Session = factor(DatTimeElapsed$Session)

# comparing first run between cTBs and sham
FR_time_cTBS = DatTimeElapsed$Time[DatTimeElapsed$TMS=='cTBS' &
                      DatTimeElapsed$Run==1]
FR_time_Sham = DatTimeElapsed$Time[DatTimeElapsed$TMS=='Sham' &
                      DatTimeElapsed$Run==1]
t.test(FR_time_cTBS,FR_time_Sham,paired = T)

pdf('./figs/TimeAftercTBS.pdf',12,5)
library(scales)
UseColor = hue_pal()(2)
p1 = ggplot(DatTimeElapsed,aes(x = Run, y = Time, fill = TMS)) +
  geom_violin(alpha = 0.8,position = position_dodge(width = 0.9)) +
  geom_point(size = 0.5,position = position_jitterdodge(seed = 1, dodge.width = 0.9)) +
  theme_classic(base_size = 22) +
  scale_fill_manual(values = c("Sham" = UseColor[2],"cTBS" = UseColor[1]))+
  labs(y='Time (minutes)',title = 'Time Elapsed after cTBS')

UseColor = hue_pal()(4)
p2 = ggplot(DatTimeElapsed,aes(x = Run, y = Time, fill = Session)) +
  geom_violin(alpha = 0.8,position = position_dodge(width = 0.9)) +
  geom_point(size = 0.5,position = position_jitterdodge(seed = 1, dodge.width = 0.9)) +
  theme_classic(base_size = 22) +
  scale_fill_manual(values = c("1" = UseColor[3],"2" = UseColor[4]))+
  labs(y='Time (minutes)',title = 'Time Elapsed after cTBS')

ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1)
dev.off()

p1
Use.time.sham = subset(DatTimeElapsed,Run == 1 & TMS == 'Sham')$Time
Use.time.cTBS = subset(DatTimeElapsed,Run == 1 & TMS == 'cTBS')$Time
summary(Use.time.sham)
sd(Use.time.sham)
summary(Use.time.cTBS)
sd(Use.time.cTBS)

Use.time.both = subset(DatTimeElapsed,Run == 1)$Time
summary(Use.time.both)
sd(Use.time.both)

#############################
# look at stimulation coordinates
Coordinates_df_L = dat[,c('Sub','TargetLX','TargetLY','TargetLZ')]
Coordinates_df_L$Label = 'OFC.L'
names(Coordinates_df_L) = c('Sub','X','Y','Z','Label')

Coordinates_df_R = dat[,c('Sub','TargetRX','TargetRY','TargetRZ')]
Coordinates_df_R$Label = 'OFC.R'
names(Coordinates_df_R) = c('Sub','X','Y','Z','Label')

Coordinates_df = rbind(Coordinates_df_L,Coordinates_df_R)


###########################
## look at data from survey

data = cbind(dat$Sub,dat[,c(32:39)])
TMSCond = dat$Cond

# reorganize uncomfort scores to A and P sessions
uncomfort = matrix(NA,31,2)
# two columns: P and A
for(i in 1:31){
  uncomfort[i,1] = ifelse(TMSCond[i]=='PA',dat$uncomfort1[i],dat$uncomfort2[i])
  uncomfort[i,2] = ifelse(TMSCond[i]=='PA',dat$uncomfort2[i],dat$uncomfort1[i])
}

boxplot(uncomfort)
summary(uncomfort)
wilcox.test(uncomfort[,1],uncomfort[,2],paired = T)
t.test(uncomfort[,1],uncomfort[,2],paired = T)

strong = matrix(NA,31,2)
for(i in 1:31){
  strong[i,1] = ifelse(TMSCond[i]=='PA',dat$strong1[i],dat$strong2[i])
  strong[i,2] = ifelse(TMSCond[i]=='PA',dat$strong2[i],dat$strong1[i])
}

boxplot(strong)
summary(strong)
wilcox.test(strong[,1],strong[,2],paired = T)
t.test(strong[,1],strong[,2],paired = T)

# write uncomfort and strong into excel
df_ratings = data.frame(cbind(Subs,uncomfort,strong))
names(df_ratings) = c('Subs','uncomfort_sham','uncomfort_cTBS',
                        'strong_sham','strong_cTBS')
write.xlsx(df_ratings,
          rowNames = F,
          file = file.path(ResDir,'Behavior','TMSratings.xlsx'))

# if noticing TMS difference from two TMS sessions
dat$notice_diff
table(dat$notice_diff)

dat$guess1
dat$guess2
guess_res = matrix(cbind(dat$guess1,dat$guess2),ncol = 2)
guess_res = apply(guess_res,1,paste,collapse = "")
cbind(guess_res,TMSCond)
TMSCond

table(guess_res == TMSCond)

table(dat$considered_TMS)


