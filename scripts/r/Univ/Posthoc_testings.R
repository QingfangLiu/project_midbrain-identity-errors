
rm(list = ls())
library(varhandle)
library(openxlsx) # read and write xlsx files
library(plyr)
library(lme4)

# define directory paths
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
NeuralDir = file.path(ResDir,'Neural')
BehDir = file.path(ResDir,'Behavior')
SubInfoDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/SubInfo'

# ggplot setup
SetupFile = file.path(AnaDir,'Scripts_R','ggplotSetup.R')
source(SetupFile)

##################################################################

# correlation b/t left OFC with other ROIs (MB, LPFC)
df1 = read.xlsx(file.path(NeuralDir,'BetaDecreasePostRev_SubsTMS.xlsx'))
df2 = read.xlsx(file.path(NeuralDir,'BetaDecreasePostRev_SubsTMS_OFC.xlsx'))
df = cbind(df1,df2[,c(4,5)]) # only including lOFC columns

pairs(df[,-1])

# add columns of cTBS - sham comparison
df = df %>%
  dplyr::mutate(MB.comp = MB.cTBS - MB.sham,
                LPFC.comp = LPFC.cTBS - LPFC.sham,
                lOFC.comp = lOFC.cTBS - lOFC.sham)
  
pairs(df[,c(8,9,10)])

# write combined MB, LPFC, left OFC results
write.xlsx(df,file = file.path(NeuralDir,'BetaDecreasePostRev_SubsTMS_comb.xlsx'),
           rowNames = F)

#########################################
# load TMS rating data
df_ratings = read.xlsx(file.path(BehDir,'TMSratings.xlsx'))

# create a df 
use.df = data.frame('Subs'=rep(df$Subs,2),
                    'TMS'=c(rep('cTBS',31),rep('sham',31)),
                    'MB'=c(df$MB.cTBS,df$MB.sham),
                    'LPFC'=c(df$LPFC.cTBS,df$LPFC.sham),
                    'lOFC'=c(df$lOFC.cTBS,df$lOFC.sham),
                    'uncomfort'=c(df_ratings$uncomfort_cTBS,df_ratings$uncomfort_sham),
                    'strong'=c(df_ratings$strong_cTBS,df_ratings$strong_sham)
                    )
# pay attention to this! They should be factors!!!
use.df = use.df %>%
  mutate(Subs = factor(Subs),
         TMS = factor(TMS)) 

# test if TMS effect still exists after accounting for uncomfort and strong
# on Midbrain
lmer1_MB <- lmer(MB ~ Subs + uncomfort + strong + (1 | Subs), 
              data = use.df)
lmer2_MB <- lmer(MB ~ Subs + TMS + uncomfort + strong + (1 | Subs), 
              data = use.df)
anova(lmer1_MB,lmer2_MB)

# on LPFC
lmer1_LPFC <- lmer(LPFC ~ Subs + uncomfort + strong + (1 | Subs), 
              data = use.df)
lmer2_LPFC <- lmer(LPFC ~ Subs + TMS + uncomfort + strong + (1 | Subs), 
              data = use.df)
anova(lmer1_LPFC,lmer2_LPFC)

# on lOFC
lmer1_lOFC <- lmer(lOFC ~ Subs + uncomfort + strong + (1 | Subs), 
                   data = use.df)
lmer2_lOFC <- lmer(lOFC ~ Subs + TMS + uncomfort + strong + (1 | Subs), 
                   data = use.df)
anova(lmer1_lOFC,lmer2_lOFC)

# plot pairwise correlations
pairs(use.df[,c(3,4,5,6,7)])

## another way to look at it is to use cTBS - sham
# draw pairwise plot showing the relationship between effects of TMS on
# ratings and on different ROIs iPE effect

df_ratings = df_ratings %>%
  dplyr::mutate(uncomfort.comp = uncomfort_cTBS - uncomfort_sham,
                strong.comp = strong_cTBS - strong_sham)
aggdf = cbind(df,df_ratings)
pairs(aggdf[,c('MB.comp','LPFC.comp','lOFC.comp',
               'uncomfort.comp','strong.comp')])


