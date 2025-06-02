
# test fMRI univariate effect using functional OFC ROI
rm(list = ls())
library(varhandle)
library(openxlsx) # read and write xlsx files
library(plyr)
library(ggpubr)

# ggplot setup
source("../Setup.R")

ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
NeuralDir = file.path(ResDir,'Neural')
SubInfoDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/SubInfo'

# source data 
sourceDir = '~/Desktop/SourceData'

# read in our data
filename = file.path(ResDir,'Neural','Betas_iPE_fourpt_clean_semiconcat_noResp_Nz.mat')
matdat = R.matlab::readMat(filename)
aggMeanBOLD = matdat$aggMeanBOLD

ROIs = c('OFC','lOFC','rOFC') # check matlab consistency
nROIs = length(ROIs)

SubInfo = read.xlsx(file.path(SubInfoDir,'SubjectConds.xlsx'))
UseTMSConds = SubInfo$Cond[SubInfo$Excluded==0]
Subs = SubInfo$Sub[SubInfo$Excluded==0]

# convert this array to df
UseAggBOLD = structure(aggMeanBOLD,.Dim = c(4L, nROIs, 2L, length(Subs)), 
                       .Dimnames = structure(list(RevLoc = c(0,-1,1,2), 
                       ROI = ROIs, Sess = 1:2, Subs = as.vector(Subs)), 
                       .Names = c("RevLoc", "ROI" ,"Sess", "Subs")))
UseAggBOLD = plyr::adply(UseAggBOLD, c(1,2,3,4))
names(UseAggBOLD) = c(names(UseAggBOLD)[1:4],'Beta')

# reorder revloc factors for the convenience of plotting
UseAggBOLD$RevLoc <- factor(UseAggBOLD$RevLoc, levels = c("-1", "0", "1","2"))
UseAggBOLD$TMS = unfactor(UseAggBOLD$Sess)
for(j in 1:length(Subs)){
  if(UseTMSConds[j]=='AP'){
    UseAggBOLD$TMS[UseAggBOLD$Subs==Subs[j]] = 3 - 
      UseAggBOLD$TMS[UseAggBOLD$Subs==Subs[j]]
  }
}

UseAggBOLD$TMS = plyr::mapvalues(UseAggBOLD$TMS,from = c(1,2),to = c('sham','cTBS'))
UseAggBOLD$TMS = factor(UseAggBOLD$TMS)

## save source data
write.xlsx(subset(UseAggBOLD,ROI=='lOFC'),
           file = file.path(sourceDir,'UseAggBOLD_lOFC.xlsx'))

################ summarize agg BOLD by TMS ###################

SummaryAggBOLD = UseAggBOLD %>%
  group_by(TMS, RevLoc, ROI) %>%
  dplyr::summarise(
    mean = mean(Beta),
    sd = sd(Beta),
    n = n(), 
    sde = sd/sqrt(n)
    )

################## Check out beta change 

DiffAgg = aggMeanBOLD[3,,,] - aggMeanBOLD[1,,,] # (rev+1) - rev
UseDiffAggBOLD = structure(DiffAgg,.Dim = c(nROIs, 2L, length(Subs)), 
                       .Dimnames = structure(list(ROI = ROIs, 
                                                  Sess = 1:2, 
                                                  Subs = as.vector(Subs)), 
                                             .Names = c("ROI", "Sess", "Subs")))
UseDiffAggBOLD = adply(UseDiffAggBOLD, c(1,2,3))
names(UseDiffAggBOLD) = c(names(UseDiffAggBOLD)[1:3],'Beta')
UseDiffAggBOLD$TMS = unfactor(UseDiffAggBOLD$Sess)

for(j in 1:length(Subs)){
  if(UseTMSConds[j]=='AP'){
    UseDiffAggBOLD$TMS[UseDiffAggBOLD$Subs==Subs[j]] = 3 - 
      UseDiffAggBOLD$TMS[UseDiffAggBOLD$Subs==Subs[j]]
  }
}

UseDiffAggBOLD$TMS = mapvalues(UseDiffAggBOLD$TMS,from = c(1,2),to = c('sham','cTBS'))
UseDiffAggBOLD = UseDiffAggBOLD %>% mutate(TMS = factor(TMS)) %>%
  arrange(ROI,Subs) %>%
  dplyr::mutate(paired = rep(1:(n()/2),each=2))

write.xlsx(subset(UseDiffAggBOLD,ROI=='lOFC'),
           file = file.path(sourceDir,'UseDiffAggBOLD_lOFC.xlsx'))

p = NULL
p[[1]] = ggplot(data = subset(SummaryAggBOLD,ROI == 'lOFC'),
       aes(x = RevLoc, y = mean, group = TMS)) + 
  geom_errorbar(aes(ymin=mean-sde, ymax=mean+sde,color=TMS), width=.2, size=1.1) +
  geom_line(aes(color=TMS),size=1.1)+
  geom_point(aes(color=TMS),show.legend=FALSE) +
  scale_colour_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  labs(x="", y = "", title = '') + 
  scale_x_discrete(labels= c("-1"=expression(rev[-1]), "0"="rev",
                         "1"=expression(rev[+1]),"2"=expression(rev[+2]))) +
  common + theme(legend.position = 'none', strip.text.x = element_blank())

high_y1 = 1.8
low_y1 = 1.6

# dat_text <- data.frame(
#   label = c("**", "n.s."),
#   ROI   = c('lOFC', 'rOFC'),
#   x     = c(1.5, 1.5),
#   y     = c(1.9, 2)
# )

dat_text <- data.frame(
  label = c("**"),
  ROI   = c('lOFC'),
  x     = c(1.5),
  y     = c(1.9)
)

p[[2]] = ggplot(data = subset(UseDiffAggBOLD,ROI == 'lOFC'),
                mapping = aes(x = TMS, y = Beta)) + 
  geom_boxplot(aes(fill=TMS),alpha=0.5, outlier.shape = NA) +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=TMS,group=paired),size=2,shape=21,alpha=1,
             position = position_dodge(0.2)) +
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  ylim(-1.8, 2) +
  labs(x="", y = "", title = '') + 
  common + 
  theme(strip.text.x = element_blank(),legend.position = "none" ) +
  geom_segment(aes(x=1,y=high_y1,xend=2,yend=high_y1),size = 0.2) +
  geom_segment(aes(x=1,y=low_y1,xend=1,yend=high_y1),size = 0.2) +
  geom_segment(aes(x=2,y=low_y1,xend=2,yend=high_y1),size = 0.2) +
  geom_text(
    data    = dat_text,
    mapping = aes(x = x, y = y, label = label),size=4.5)

pdf(file.path(NeuralDir,'fMRIUni_OFC.pdf'),8,4)
fig = ggarrange(p[[1]], p[[2]], widths = c(1.3,1), align = "hv")
print(fig)
dev.off()


# t tests and save post-rev beta decrease
attach(UseDiffAggBOLD)
t.test(Beta[ROI=='OFC' & TMS=='cTBS'], Beta[ROI=='OFC' & TMS=='sham'], 
       paired = T, alternative = 'g')
t.test(Beta[ROI=='lOFC' & TMS=='cTBS'], Beta[ROI=='lOFC' & TMS=='sham'],
       paired = T, alternative = 'g')
t.test(Beta[ROI=='rOFC' & TMS=='cTBS'], Beta[ROI=='rOFC' & TMS=='sham'],
       paired = T, alternative = 'g')

BetaChangeSave = data.frame(
  'Subs' = Subs,
  'OFC-cTBS' = Beta[ROI=='OFC' & TMS=='cTBS'], 
  'OFC-sham' = Beta[ROI=='OFC' & TMS=='sham'],
  'lOFC-cTBS' = Beta[ROI=='lOFC' & TMS=='cTBS'], 
  'lOFC-sham' = Beta[ROI=='lOFC' & TMS=='sham'],
  'rOFC-cTBS' = Beta[ROI=='rOFC' & TMS=='cTBS'], 
  'rOFC-sham' = Beta[ROI=='rOFC' & TMS=='sham'])

write.xlsx(BetaChangeSave,
           file = file.path(NeuralDir,'BetaDecreasePostRev_SubsTMS_OFC.xlsx'),
           rowNames = F)
detach(UseDiffAggBOLD)

##################################################################
# compare post-reversal decrease from different sessions

attach(UseDiffAggBOLD)
t.test(Beta[ROI=='OFC' & Sess==1], Beta[ROI=='OFC' & Sess==2], paired = T)
t.test(Beta[ROI=='lOFC' & Sess==1], Beta[ROI=='lOFC' & Sess==2], paired = T)
t.test(Beta[ROI=='rOFC' & Sess==1], Beta[ROI=='rOFC' & Sess==2], paired = T)
detach(UseDiffAggBOLD)


