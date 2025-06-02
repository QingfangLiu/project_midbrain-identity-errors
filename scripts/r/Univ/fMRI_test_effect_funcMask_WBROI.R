
# test fMRI univariate effect using WB ROIs (Midbrain, LPFC)

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
filename = file.path(ResDir,'Neural','Betas_iPE_fourpt_clean_semiconcat_noResp_Nz_WB.mat')
matdat = readMat(filename)
ROIs = c('MB','LPFC') # check matlab consistency
roi_label = c('Midbrain','LPFC')
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
UseAggBOLD = adply(UseAggBOLD, c(1,2,3,4))
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

UseAggBOLD$TMS = mapvalues(UseAggBOLD$TMS,from = c(1,2),to = c('sham','cTBS'))
UseAggBOLD$TMS = factor(UseAggBOLD$TMS)

## save source data
write.xlsx(UseAggBOLD,
           file = file.path(sourceDir,'UseAggBOLD.xlsx'))

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
  dplyr::mutate(paired = rep(1:(n()/2),each=2))


write.xlsx(UseDiffAggBOLD,
           file = file.path(sourceDir,'UseDiffAggBOLD.xlsx'))

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
  labs(x="", y = "", title = '') + 
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
  labs(x="", y = "", title = '') + 
  common +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.text.x = element_blank(),legend.position = "none" )
rctr = rctr + 1
}

high_y1 = 1.5
low_y1 = 1.3
p[[2]] = p[[2]] + ylim(-2.5,1.6) +
   geom_segment(aes(x=1,y=high_y1,xend=2,yend=high_y1),size = 0.2) +
   geom_segment(aes(x=1,y=low_y1,xend=1,yend=high_y1),size = 0.2) +
   geom_segment(aes(x=2,y=low_y1,xend=2,yend=high_y1),size = 0.2) +
   annotate(geom = "text", x = c(1.5), y = c(1.59), label = '*', size = 6)

pdf(file.path(NeuralDir,'fMRIUni.pdf'),8,4)
fig = ggpubr::ggarrange(p[[1]], p[[2]], widths = c(1.3,1),
                ncol = 2, nrow = 1, align = "hv")
print(fig)
dev.off()

high_y2 = 1
low_y2 = 0.8
p[[4]] = p[[4]] +
  geom_segment(aes(x=1,y=high_y2,xend=2,yend=high_y2),size = 0.2) +
  geom_segment(aes(x=1,y=low_y2,xend=1,yend=high_y2),size = 0.2) +
  geom_segment(aes(x=2,y=low_y2,xend=2,yend=high_y2),size = 0.2) +
  annotate(geom = "text", x = c(1.5), y = c(1.1), label = '**', size = 6)

pdf(file.path(NeuralDir,'fMRIUni_LPFC.pdf'),8,4)
fig = ggpubr::ggarrange(p[[3]], p[[4]], widths = c(1.3,1), ncol = 2, nrow = 1, align = "hv")
print(fig)
dev.off()

################ stat tests ###################
# t tests and save (rev+1) - rev

attach(UseDiffAggBOLD)
BetaChangeSave = data.frame(
                      'Subs' = Subs,
                      'MB-cTBS' = Beta[ROI=='MB' & TMS=='cTBS'], 
                      'MB-sham' = Beta[ROI=='MB' & TMS=='sham'],
                      'LPFC-cTBS' = Beta[ROI=='LPFC' & TMS=='cTBS'], 
                      'LPFC-sham' = Beta[ROI=='LPFC' & TMS=='sham'])
write.xlsx(BetaChangeSave,
           file = file.path(NeuralDir,'BetaDecreasePostRev_SubsTMS.xlsx'),
           rowNames = F)
detach(UseDiffAggBOLD)

# paired t-test result reported in the paper
t.test(BetaChangeSave$MB.cTBS, BetaChangeSave$MB.sham, 
       paired = T, alternative = 'g')
t.test(BetaChangeSave$LPFC.cTBS, BetaChangeSave$LPFC.sham,
       paired = T, alternative = 'g')

##################################################################
## post-hoc analysis on rev and rev+1 locations separately

# Midbrain
# test on rev trials
beta_MB_rev_cTBS = subset(UseAggBOLD,RevLoc == 0 & ROI == 'MB' & TMS == 'cTBS')$Beta
beta_MB_rev_sham = subset(UseAggBOLD,RevLoc == 0 & ROI == 'MB' & TMS == 'sham')$Beta
t.test(beta_MB_rev_cTBS,beta_MB_rev_sham,paired = T,alternative = 'l')

# test on rev+1 trials
beta_MB_rev1_cTBS = subset(UseAggBOLD,RevLoc == 1 & ROI == 'MB' & TMS == 'cTBS')$Beta
beta_MB_rev1_sham = subset(UseAggBOLD,RevLoc == 1 & ROI == 'MB' & TMS == 'sham')$Beta
t.test(beta_MB_rev1_cTBS,beta_MB_rev1_sham,paired = T,alternative = 'g')

# also on rev-1 and rev+2 trials
beta_MB_revn1_cTBS = subset(UseAggBOLD,RevLoc == -1 & ROI == 'MB' & TMS == 'cTBS')$Beta
beta_MB_revn1_sham = subset(UseAggBOLD,RevLoc == -1 & ROI == 'MB' & TMS == 'sham')$Beta
t.test(beta_MB_revn1_cTBS,beta_MB_revn1_sham,paired = T,alternative = 'g')

beta_MB_rev2_cTBS = subset(UseAggBOLD,RevLoc == 2 & ROI == 'MB' & TMS == 'cTBS')$Beta
beta_MB_rev2_sham = subset(UseAggBOLD,RevLoc == 2 & ROI == 'MB' & TMS == 'sham')$Beta
t.test(beta_MB_rev2_cTBS,beta_MB_rev2_sham,paired = T,alternative = 'g')

# LPFC
# test on rev trials
beta_LPFC_rev_cTBS = subset(UseAggBOLD,RevLoc == 0 & ROI == 'LPFC' & TMS == 'cTBS')$Beta
beta_LPFC_rev_sham = subset(UseAggBOLD,RevLoc == 0 & ROI == 'LPFC' & TMS == 'sham')$Beta
t.test(beta_LPFC_rev_cTBS,beta_LPFC_rev_sham,paired = T,alternative = 'l')

# test on rev+1 trials
beta_LPFC_rev1_cTBS = subset(UseAggBOLD,RevLoc == 1 & ROI == 'LPFC' & TMS == 'cTBS')$Beta
beta_LPFC_rev1_sham = subset(UseAggBOLD,RevLoc == 1 & ROI == 'LPFC' & TMS == 'sham')$Beta
t.test(beta_LPFC_rev1_cTBS,beta_LPFC_rev1_sham,paired = T,alternative = 'g')

# adding a test
beta_LPFC_rev_cTBS = subset(UseAggBOLD,RevLoc == 0 & ROI == 'LPFC' & TMS == 'cTBS')$Beta
beta_LPFC_nonrev_cTBS = UseAggBOLD %>%
  subset(RevLoc!=0 & ROI == 'LPFC' & TMS == 'cTBS') %>%
  group_by(Subs) %>%
  dplyr::summarise(Beta = mean(Beta))
beta_LPFC_diff_cTBS = beta_LPFC_rev_cTBS - beta_LPFC_nonrev_cTBS$Beta

beta_LPFC_rev_sham = subset(UseAggBOLD,RevLoc == 0 & ROI == 'LPFC' & TMS == 'sham')$Beta
beta_LPFC_nonrev_sham = UseAggBOLD %>%
  subset(RevLoc!=0 & ROI == 'LPFC' & TMS == 'sham') %>%
  group_by(Subs) %>%
  dplyr::summarise(Beta = mean(Beta))
beta_LPFC_diff_sham = beta_LPFC_rev_sham - beta_LPFC_nonrev_sham$Beta
t.test(beta_LPFC_diff_cTBS,beta_LPFC_diff_sham,paired = T,alternative = 'l')


##################################################################
# compare post-reversal decrease from different sessions

SummaryAggBOLD = UseAggBOLD %>%
  group_by(Sess, RevLoc, ROI) %>%
  dplyr::summarise(
    mean = mean(Beta),
    sd = sd(Beta),
    n = n(), 
    sde = sd/sqrt(n)
  )

q = NULL
q[[1]] = ggplot(data = SummaryAggBOLD,
                aes(x = RevLoc, y = mean, group = Sess)) + 
  geom_errorbar(aes(ymin=mean-sde, ymax=mean+sde,color=Sess), width=.2, size=1.1) +
  geom_line(aes(color=Sess),size=1.1)+
  geom_point(aes(color=Sess),show.legend=FALSE) +
  facet_wrap(~ ROI) + labs(x="", y = "") + 
  scale_x_discrete(labels=
                     c("-1"=expression(rev[-1]), "0"="rev",
                       "1"=expression(rev[+1]),"2"=expression(rev[+2]))) +
  common + theme(legend.position = c(0.15, 0.95), 
        legend.title=element_blank(),
        strip.text.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=3))) # only change legend size

q[[2]] = ggplot(data = UseDiffAggBOLD,mapping = aes(x = Sess, y = Beta, fill = Sess)) + 
  geom_boxplot(alpha=0.5) + facet_wrap(~ROI) +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Sess,group=paired),size=2,shape=21,alpha=1,
             position = position_dodge(0.2)) +
  labs(y = "", x = '') + ylim(-1.8, 4.5) +
  common + theme(strip.text.x = element_blank(),legend.position = "none" )

pdf(file.path(NeuralDir,'fMRIUni_Sess.pdf'),9,4)
gridExtra::grid.arrange(q[[1]], q[[2]], nrow = 1)
dev.off()

# test sess effect 
attach(UseDiffAggBOLD)
t.test(Beta[ROI=='MB' & Sess==1], Beta[ROI=='MB' & Sess==2], paired = T)
t.test(Beta[ROI=='LPFC' & Sess==1], Beta[ROI=='LPFC' & Sess==2], paired = T)
detach(UseDiffAggBOLD)


# test TMS and Sess effects jointly
# lmer1: including both TMS and Sess
# lmer2: including only TMS
# lmer3: including only Sess

# focus on MB
lmer1_MB <- lmer(Beta ~ Subs + TMS + Sess + (1 | Subs), 
               data = subset(UseDiffAggBOLD,ROI=='MB'))
lmer2_MB <- lmer(Beta ~ Subs + TMS + (1 | Subs), 
               data = subset(UseDiffAggBOLD,ROI=='MB'))
lmer3_MB <- lmer(Beta ~ Subs + Sess + (1 | Subs), 
              data = subset(UseDiffAggBOLD,ROI=='MB'))
anova(lmer1_MB,lmer2_MB) # Sess useful? (when already having TMS)
anova(lmer1_MB,lmer3_MB) # TMS useful? (when already having Sess)

# focus on LPFC
lmer1_LPFC <- lmer(Beta ~ Subs + TMS + Sess + (1 | Subs), 
              data = subset(UseDiffAggBOLD,ROI=='LPFC'))
lmer2_LPFC <- lmer(Beta ~ Subs + TMS + (1 | Subs), 
              data = subset(UseDiffAggBOLD,ROI=='LPFC'))
lmer3_LPFC <- lmer(Beta ~ Subs + Sess + (1 | Subs), 
              data = subset(UseDiffAggBOLD,ROI=='LPFC'))
anova(lmer1_LPFC,lmer2_LPFC) # Sess useful? (when already having TMS)
anova(lmer1_LPFC,lmer3_LPFC) # TMS useful? (when already having Sess)

