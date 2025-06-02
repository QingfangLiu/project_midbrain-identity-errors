
# test fMRI univariate effect using other ROIs

rm(list = ls())

# define directory paths
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
AnaDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis'

# dependent directory paths
NeuralDir = file.path(ResDir,'Neural')
BehDir = file.path(ResDir,'Behavior')
MVPADir = file.path(ResDir,'MVPA_ROI')
SubInfoDir = file.path(AnaDir,'SubInfo')

# ggplot setup
SetupFile = file.path(AnaDir,'Scripts_R','Setup.R')
source(SetupFile)

# source data 
sourceDir = '~/Desktop/SourceData'

##################################################################

# read in our data
filename = file.path(ResDir,'Neural','Betas_iPE_fourpt_clean_semiconcat_noResp_Nz_others.mat')
matdat = readMat(filename)

ROIs = c('Insula','mPFC','Striatum','Thalamus','Amygdala') # check matlab consistency
nROIs = length(ROIs)
aggMeanBOLD = matdat$aggMeanBOLD

SubInfo = read.xlsx(file.path(SubInfoDir,'SubjectConds.xlsx'))
UseTMSConds = SubInfo$Cond[SubInfo$Excluded==0]
Subs = SubInfo$Sub[SubInfo$Excluded==0]

# convert this array to df
UseAggBOLD = structure(aggMeanBOLD,
                       .Dim = c(4L, nROIs, 2L, length(Subs)), 
                       .Dimnames = structure(list(RevLoc = c(0,-1,1,2), 
                                                  ROI = ROIs, 
                                                  Sess = 1:2, 
                                                  Subs = as.vector(Subs)), 
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
write.xlsx(UseAggBOLD,file = file.path(sourceDir,'UseAggBOLD_others.xlsx'))

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

DiffAgg = aggMeanBOLD[3,,,] - aggMeanBOLD[1,,,]  # (rev+1) - rev
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
UseDiffAggBOLD = UseDiffAggBOLD %>% 
  mutate(TMS = factor(TMS)) %>%
  arrange(ROI,Subs) %>%
  mutate(paired = rep(1:(n()/2),each=2))

write.xlsx(UseDiffAggBOLD,file = file.path(sourceDir,'UseDiffAggBOLD_others.xlsx'))

################ plots ###################

rctr = 0
p = NULL
for(roi in ROIs){
p[[2*rctr+1]] = ggplot(data = subset(SummaryAggBOLD,ROI == roi),
       aes(x = RevLoc, y = mean, group = TMS)) + 
  geom_errorbar(aes(ymin=mean-sde, ymax=mean+sde,color=TMS), width=.2, size=1.1) +
  geom_line(aes(color=TMS),size=1.1)+
  geom_point(aes(color=TMS),show.legend=FALSE) +
  scale_colour_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  labs(x="", y = "", title = "") + 
  scale_x_discrete(labels=c("-1"=expression(rev[-1]), "0"="rev",
                         "1"=expression(rev[+1]),"2"=expression(rev[+2]))) + common +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none", strip.text.x = element_blank())
p[[2*rctr+2]] = ggplot(data = subset(UseDiffAggBOLD,ROI == roi),
                mapping = aes(x = TMS, y = Beta, fill = TMS)) + 
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=TMS,group=paired),size=2,shape=21,alpha=1,
             position = position_dodge(0.2)) +
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  labs(x="", y="", title = '') + 
  common +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.text.x = element_blank(),legend.position = "none" )
rctr = rctr + 1
}


# range of y limits
ylim_low = c(-2,-3,-3,-3,-2)
ylim_high = c(0.9,1.4,2.2,1.4,2.2)

# line segments
high_y_line = c(0.8,0.9,2.1,0.9,2.1)
low_y_line = c(0.7,0.8,2,0.8,2)
  
# somehow it doesn't work well so I gave up & drew them in AI
# for(k in 1:5){
#   id = 2*k
#   p[[id]] = p[[id]] +
#     ylim(ylim_low[k],ylim_high[k]) +
#     geom_segment(aes(x=1,y=high_y_line[k],xend=2,yend=high_y_line[k]),size = 0.2) +
#     geom_segment(aes(x=1,y=low_y_line[k],xend=1,yend=high_y_line[k]),size = 0.2) +
#     geom_segment(aes(x=2,y=low_y_line[k],xend=2,yend=high_y_line[k]),size = 0.2)
#     #annotate(geom = "text", x = c(1.5), y = c(0.6), label = 'ns.', size = 6)
#   dev.off()
#   
# }

pdf(file.path(NeuralDir,'fMRIUni_others.pdf'),14,6)
fig = ggarrange(p[[1]], p[[3]], p[[5]], p[[7]], p[[9]], 
                p[[2]], p[[4]], p[[6]], p[[8]], p[[10]], 
                ncol = 5, nrow = 2, align = "hv")
print(fig)
dev.off()



################ stat tests ###################
# t tests and save (rev+1) - rev

for(i in 1:nROIs){
  print(ROIs[i])
  dat_cTBS = subset(UseDiffAggBOLD,ROI==ROIs[i] & TMS =='cTBS')$Beta
  dat_sham = subset(UseDiffAggBOLD,ROI==ROIs[i] & TMS =='sham')$Beta
  print(t.test(dat_cTBS, dat_sham, paired = T, alternative = 'g'))
}

# testing result:
# only mPFC shows t(30)=1.94, p=0.031

##################################################################
## post-hoc analysis on rev and rev+1 locations separately
# all ns

for(i in 1:nROIs){
  print(ROIs[i])
  
# test on rev trials
beta_rev_cTBS = subset(UseAggBOLD,RevLoc == 0 & ROI == ROIs[i] & TMS == 'cTBS')$Beta
beta_rev_sham = subset(UseAggBOLD,RevLoc == 0 & ROI == ROIs[i] & TMS == 'sham')$Beta
print(t.test(beta_rev_cTBS,beta_rev_sham,paired = T,alternative = 'l'))

# test on rev+1 trials
beta_rev1_cTBS = subset(UseAggBOLD,RevLoc == 1 & ROI == ROIs[i] & TMS == 'cTBS')$Beta
beta_rev1_sham = subset(UseAggBOLD,RevLoc == 1 & ROI == ROIs[i] & TMS == 'sham')$Beta
print(t.test(beta_rev1_cTBS,beta_rev1_sham,paired = T,alternative = 'g'))

}


