
# process behavior data
# starting with OdorPredRes.xlsx data
# end with saving OdorPredRes_processed.xlsx

rm(list = ls())

library(plyr)
library(tidyverse)
library(openxlsx)

# define directory paths
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
AnaDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis'

Behdir = file.path(ResDir,'Behavior')
dat = read.xlsx(file.path(Behdir,'OdorPredRes.xlsx'))
datinfo = read.xlsx(file.path(Behdir,'SubjectConds.xlsx'))

# exclude 4 subs: only use 31 subjects' data
subs = subset(datinfo,Excluded==0)$Sub 
nsubs = length(subs)
dat = subset(dat,Sub %in% subs)

# add Age, Sex, SideFirst, MT into dataset
dat['Age'] <- NA
dat['Gender'] <- NA
dat['SideFirst'] <- NA
dat['MT'] <- NA

for(j in 1:nsubs){
  dat$Age[dat$Sub==subs[j]] = datinfo$Age[datinfo$Sub==subs[j]]
  dat$Gender[dat$Sub==subs[j]] = datinfo$Sex[datinfo$Sub==subs[j]]
  dat$SideFirst[dat$Sub==subs[j]] = datinfo$SideFirst[datinfo$Sub==subs[j]]
  dat$MT[dat$Sub==subs[j]] = datinfo$MT[datinfo$Sub==subs[j]]
}

#################################

Add_processed_vars = function(x){
  
  # a function that adds variables to run-wise data under the same cue
  # reversal locations: -1,0,1,2
  # valid: 0 or 1
  # block: 1,2,3,4,5,6; add var block 1-6 to each run & each cue
  
  # to test
  # x=subset(dat1,Sub==1 & Sess==1 & Run==1 & Cue==1)
  
  LocRev = which(x$Rev==1)
  RevLoc = Valid = Block = rep(NA,nrow(x))
  RevLoc[LocRev-1] = -1
  RevLoc[LocRev] = 0
  RevLoc[LocRev+1] = 1
  RevLoc[LocRev+2] = 2
  
  Valid[x$Rev==0] = x$Acc[x$Rev==0]
  Valid[LocRev] = as.numeric(x$Resp[LocRev] == x$Odor[LocRev-1])
  
  for(z in c(-1,0,1,2)){
    Block[!is.na(RevLoc) & RevLoc==z] = c(1:6)
  }
  
  x$RevLoc = RevLoc
  x$Valid = Valid
  x$Block = Block
  return(x)
}

# for each subject, remove all trials with RT beyond 3 std
# then normalize RT by substracting mean RTs
ProcessRT = function(x){ # x: sub-wise data
  meanRT = mean(x$RT)
  sdRT = sd(x$RT)
  GoodTrials = (x$RT > meanRT - 3*sdRT) & (x$RT < meanRT + 3*sdRT)
  x = x[GoodTrials,]  
  meanRT = mean(x$RT)
  sdRT = sd(x$RT)
  x$NormRT = x$RT - meanRT
  x$zRT = x$NormRT/sdRT
  return(x)
}


AllRevdf = dat %>% 
  ############# add columns RevLoc, Valid, Block
  group_by(Sub,Sess,Run,Cue) %>%
  do(Add_processed_vars(.)) %>% 
  ungroup() %>% 
  ############
  filter(complete.cases(RT))  %>%  # remove all NA RT values
  ############# process RT for each subject
  group_by(Sub) %>%
  do(ProcessRT(.)) %>% 
  ungroup() %>% 
  ############## only including 4-trial window trials
  filter(complete.cases(RevLoc)) %>% 
  mutate(TMS=recode(factor(TMS),'1'='sham','2'='cTBS'),
         RevLoc=factor(RevLoc),
         logRT=log(RT)) %>%
  ############# reorder rows and select columns
  arrange(Sub,Sess,Run,Trial) %>%
  select(Sub,Sess,TMS,Run,Trial, 
       Cue,Odor,Rev,RevLoc,Block,
       RT,NormRT,zRT,logRT,
       Resp,Acc,Valid,
       Age,Gender,SideFirst,MT) 
  
# write this df to disk
write.xlsx(AllRevdf,file.path(Behdir,'OdorPredRes_processed.xlsx'))

