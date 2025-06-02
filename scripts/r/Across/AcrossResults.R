

# this script conducts analyses across results from diff sources
# TMS ratings & behavior
# TMS ratings & ROI-based MVPA
# TMS ratings & ROI-based iPE
# TMS ratings & Global Conn

rm(list = ls())

# define directory paths
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
AnaDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis'

# dependent directory paths
BehDir = file.path(ResDir,'Behavior')
NeuralDir = file.path(ResDir,'Neural')
MVPADir = file.path(ResDir,'MVPA_ROI')
GCDir = file.path(ResDir,'GlobalConn_ROI')

SubInfoDir = file.path(AnaDir,'SubInfo')

# ggplot setup
SetupFile = file.path(AnaDir,'Scripts_R','Setup.R')
source(SetupFile)


##################################################################
# test if TMS ratings affect behavioral effect
##################################################################

df_ratings = read.xlsx(file.path(BehDir,'TMSratings.xlsx'))
df_beh = read.xlsx(file.path(BehDir,'beh_signature_df_SubsTMS.xlsx'))

use.df = df_beh %>%
  arrange(TMS) %>%
  mutate(uncomfort = c(df_ratings$uncomfort_sham,df_ratings$uncomfort_cTBS),
         strong = c(df_ratings$strong_sham,df_ratings$strong_cTBS),
         TMS = factor(TMS),
         Sub = factor(Sub)
  )

# lmer1: baseline model with subs and ratings
# lmer2: adding TMS to lmer1
lmer1 <- lmer(beh_signature ~ Sub + uncomfort + strong + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
lmer2 <- lmer(beh_signature ~ Sub + TMS + uncomfort + strong + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
anova(lmer1,lmer2) 
# is TMS useful after accounting for two ratings?
# p = 0.011

##################################################################
# test if TMS ratings affect ROI-based MVPA result
##################################################################

df_ratings = read.xlsx(file.path(BehDir,'TMSratings.xlsx'))
df_IdRep = read.xlsx(file.path(MVPADir,'IdentityRep_SubsTMS.xlsx'))

# what ROIs are in the data?
table(df_IdRep$ROI)

use.df = df_IdRep %>%
  subset(ROI == 'Func_OFC') %>%
  arrange(TMS) %>%
  mutate(uncomfort = c(df_ratings$uncomfort_cTBS,df_ratings$uncomfort_sham),
         strong = c(df_ratings$strong_cTBS,df_ratings$strong_sham),
         TMS = factor(TMS),
         Sub = factor(Sub)
         )

# lmer1: baseline model with subs and ratings
# lmer2: adding TMS to lmer1
lmer1 <- lmer(Representation ~ Sub + uncomfort + strong + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
lmer2 <- lmer(Representation ~ Sub + TMS + uncomfort + strong + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
anova(lmer1,lmer2) # is TMS useful after accounting for two ratings?
# Func_OFC: p = 0.043
# Func_LPFC: p = 0.52


##################################################################
# test if TMS ratings affect ROI-based iPE result
##################################################################

df_ratings = read.xlsx(file.path(BehDir,'TMSratings.xlsx'))
df_iPE = read.xlsx(file.path(NeuralDir,'BetaDecreasePostRev_SubsTMS_comb.xlsx'))

use.df = UseDiffAggBOLD %>%
  subset(ROI == 'LPFC') %>%
  arrange(TMS) %>%
  mutate(uncomfort = c(df_ratings$uncomfort_cTBS,df_ratings$uncomfort_sham)) %>%
  mutate(strong = c(df_ratings$strong_cTBS,df_ratings$strong_sham))

# lmer1: baseline model with subs and ratings
# lmer2: adding TMS to lmer1
# lmer3: adding Sess to lmer1
# lmer4: adding both TMS and Sess to lmber1
lmer1 <- lmer(Beta ~ Subs + uncomfort + strong + (1 | Subs), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
lmer2 <- lmer(Beta ~ Subs + TMS + uncomfort + strong + (1 | Subs), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
lmer3 <- lmer(Beta ~ Subs + Sess + uncomfort + strong + (1 | Subs), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
lmer4 <- lmer(Beta ~ Subs + TMS + Sess + uncomfort + strong + (1 | Subs), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))

anova(lmer1,lmer2)
anova(lmer3,lmer4)


##################################################################
# test if TMS ratings affect GC ROI results
##################################################################

df_ratings = read.xlsx(file.path(BehDir,'TMSratings.xlsx'))
df_GC = read.xlsx(file.path(GCDir,'GC_SubsTMS.xlsx'))

# test goal: OFC, Time = 1,2
use.df = df_GC %>%
  mutate(Sub = as.double(Sub)) %>% # Sub can't be char if want to match order
  subset(Time %in% c(2)) %>%
  group_by(Sub,TMS) %>%
  dplyr::summarise(GC = mean(OFC)) %>%
  ungroup() %>%
  arrange(TMS,Sub) %>%
  mutate(uncomfort = c(df_ratings$uncomfort_cTBS,df_ratings$uncomfort_sham),
         strong = c(df_ratings$strong_cTBS,df_ratings$strong_sham)
         ) 

# lmer1: baseline model with subs and ratings
# lmer2: adding TMS to lmer1

lmer1_OFC_GC <- lmer(GC ~ Sub + uncomfort + strong + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
lmer2_OFC_GC <- lmer(GC ~ Sub + TMS + uncomfort + strong + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
anova(lmer1_OFC_GC,lmer2_OFC_GC)

# test goal: LPFC, Time = 2,3,4
use.df = df_GC %>%
  mutate(Sub = as.double(Sub)) %>%
  subset(Time %in% c(3)) %>%
  group_by(Sub,TMS) %>%
  dplyr::summarise(GC = mean(LPFC)) %>%
  ungroup() %>%
  arrange(TMS,Sub) %>%
  mutate(uncomfort = c(df_ratings$uncomfort_cTBS,df_ratings$uncomfort_sham),
         strong = c(df_ratings$strong_cTBS,df_ratings$strong_sham)
  ) 

# lmer1: baseline model with subs and ratings
# lmer2: adding TMS to lmer1
lmer1_LPFC_GC <- lmer(GC ~ Sub + uncomfort + strong + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
lmer2_LPFC_GC <- lmer(GC ~ Sub + TMS + uncomfort + strong + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
anova(lmer1_LPFC_GC,lmer2_LPFC_GC)


##################################################################
# test if subject-level variables are correlated with ROI-based iPE result
##################################################################

##################################################################
############ check across-subjects correlation between ROIs
# MB and LPFC
##################################################################


# post-reversal fMRI response decrease
attach(UseDiffAggBOLD)
RespMB = Beta[ROI=='MB']
RespLPFC = Beta[ROI=='LPFC']
TMSconds = unfactor(TMS[ROI=='LPFC'])
TMScolors = TMSconds
detach(UseDiffAggBOLD)

cor_df = data.frame('RespMB' = RespMB,
                    'RespLPFC' = RespLPFC,
                    'TMSconds' = TMSconds,
                    'Sex'=rep(SubInfo$Sex[SubInfo$Excluded==0],each=2))

p = NULL
p[[1]] = ggplot(cor_df,aes(x = RespMB, y = RespLPFC, color = TMSconds)) +
  labs(y = "", x = '') +
  geom_smooth(aes(group = TMSconds),show.legend=FALSE,
              method=lm, se=F, linetype="solid",size = 1.5) +
  geom_point(show.legend=F, size = 1.5, aes(shape = Sex)) +
  scale_color_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  common 

cor.test(RespMB[TMSconds=='cTBS'],RespLPFC[TMSconds=='cTBS'])
cor.test(RespMB[TMSconds=='sham'],RespLPFC[TMSconds=='sham'])

## then compare TMS-effects on post-rev decrease between ROIs
ChangeMB = RespMB[TMSconds=='sham'] - RespMB[TMSconds=='cTBS']
ChangeLPFC = RespLPFC[TMSconds=='sham'] - RespLPFC[TMSconds=='cTBS']
change_cor_df = data.frame('ChangeMB' = ChangeMB,
                           'ChangeLPFC' = ChangeLPFC,
                           'Sex'=SubInfo$Sex[SubInfo$Excluded==0])
change_cor_df$Sex = mapvalues(change_cor_df$Sex,c('F','M'),c('females','males'))

p[[2]] = ggplot(change_cor_df,aes(x = ChangeMB, y = ChangeLPFC)) +
  labs(y = "", x = '') +
  geom_smooth(method=lm, se=F, linetype="solid",size = 1.5, color = 'gray') +
  stat_cor(label.sep='\n',
           label.x.npc = 0.6,label.y.npc = 0.2,
           r.accuracy = 0.01, size = 5) +
  geom_point(aes(shape = Sex),size = 2.5) +
  common +
  theme(legend.position = c(0.2,0.9),
        legend.text = element_text(size=16),
        legend.title= element_blank())
cor.test(ChangeMB,ChangeLPFC)

pdf(file.path(NeuralDir,'fMRI_corr_ROIs.pdf'),4.5,4)
ggpubr::ggarrange(p[[2]], align = 'h', nrow = 1, ncol = 1)
dev.off()


# compare TMS effect on different orders
group1 = SubInfo$Cond[SubInfo$Excluded==0]=='PA'
group2 = SubInfo$Cond[SubInfo$Excluded==0]=='AP'
t.test(ChangeMB[group1],ChangeMB[group2])
t.test(ChangeLPFC[group1],ChangeLPFC[group2])

# correlate age with TMS effect
cor.test(SubInfo$Age[SubInfo$Excluded==0],ChangeMB)
cor.test(SubInfo$Age[SubInfo$Excluded==0],ChangeLPFC)

# correlate TMS intensity with TMS effect
cor.test(SubInfo$TMS[SubInfo$Excluded==0],ChangeMB)
cor.test(SubInfo$TMS[SubInfo$Excluded==0],ChangeLPFC)

# compare TMS effect on different sex
group1 = SubInfo$Sex[SubInfo$Excluded==0]=='F'
group2 = SubInfo$Sex[SubInfo$Excluded==0]=='M'
t.test(ChangeMB[group1],ChangeMB[group2])
#wilcox.test(ChangeMB[group1],ChangeMB[group2])
t.test(ChangeLPFC[group1],ChangeLPFC[group2])
#wilcox.test(ChangeLPFC[group1],ChangeLPFC[group2])



##################################################################
## check across-subjects correlation b/t beh effect & neural effect
##################################################################

beh_effect = df_beh$beh_signature[df_beh$TMS=='cTBS'] - 
  df_beh$beh_signature[df_beh$TMS=='sham']

plot(beh_effect,df_iPE$MB.comp)
cor.test(beh_effect,df_iPE$MB.comp)

plot(beh_effect,df_iPE$LPFC.comp)
cor.test(beh_effect,df_iPE$LPFC.comp)

plot(beh_effect,df_iPE$lOFC.comp)
cor.test(beh_effect,df_iPE$lOFC.comp)

df_IdRep = subset(df_IdRep,ROI == 'Func_OFC')
IdRep_effect = df_IdRep$Representation[df_IdRep$TMS=='cTBS'] - 
  df_IdRep$Representation[df_IdRep$TMS=='sham']
plot(beh_effect,IdRep_effect)
cor.test(beh_effect,IdRep_effect)



