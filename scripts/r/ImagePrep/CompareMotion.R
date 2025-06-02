
rm(list = ls())

# folders
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
AnaDir = '~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis'

Behdir = file.path(ResDir,'Behavior')
MRIPreDir = file.path(ResDir,'MRIPreprocessOut')
GCDir = file.path(ResDir,'GlobalConn_ROI')

# ggplot setup
SetupFile = file.path(AnaDir,'Scripts_R','Setup.R')
source(SetupFile)

# source data 
#sourceDir = '~/Desktop/SourceData'
sourceDir = '~/Documents/Projects/SPEEDTMS_GitHub'

##################################################################

# read in our data
matdat = readMat(file.path(MRIPreDir,'all_motions.mat'))
all_motions = matdat$all.motions

datinfo = read.xlsx(file.path(Behdir,'SubjectConds.xlsx'))
Subs = subset(datinfo,Excluded==0)$Sub 
nSubs = length(Subs)
UseTMSConds = subset(datinfo,Excluded==0)$Cond

nSess = 2
nruns = 3
ntimes = nruns * 2

# Para.names = c('x','y','z','pitch','roll','yaw')

# average translation and rotation parameters
avg_motions_translation = apply(all_motions[1:3,, , ], c(2, 3, 4), mean)
avg_motions_rotation = apply(all_motions[4:6,, , ], c(2, 3, 4), mean)
# note that those rotation parameters read from saved files
# are on radians (not degrees)
# so to transform rotation parameters to arc length: 
# simply multiplying radians with assumed brain radius (50mm)

# Combine the results into a new array
avg_motions = array(NA,dim = c(3,2,6,31))
avg_motions[1,,,] = avg_motions_translation
avg_motions[2,,,] = avg_motions_rotation

# reorganize all motions into mm unit
avg_motions_mm = array(NA,dim = c(2,2,6,31))
avg_motions_mm[1,,,] = avg_motions_translation
avg_motions_mm[2,,,] = avg_motions_rotation * 50
avg_motions_mm = apply(avg_motions_mm,c(2,3,4),mean)

# then add this to the big array
avg_motions[3,,,] = avg_motions_mm

Para.names = c('translation','rotation','mm')
nParas = length(Para.names)

# convert this array to df
df_motions = structure(avg_motions,.Dim = c(3, nSess, ntimes, nSubs), 
                       .Dimnames = structure(list(
                         Paras = Para.names,
                         Sess = 1:2, 
                         Time = 1:6, 
                         Subs = as.vector(Subs)), 
                       .Names = c("Paras", "Sess" ,"Time", "Subs")))
df_motions = plyr::adply(df_motions, c(1,2,3,4))
names(df_motions) = c(names(df_motions)[1:4],'value')
df_motions$TMS = unfactor(df_motions$Sess)
for(j in 1:nSubs){
  if(UseTMSConds[j]=='AP'){
    df_motions$TMS[df_motions$Subs==Subs[j]] = 3 - 
      df_motions$TMS[df_motions$Subs==Subs[j]]
  }
}

df_motions$TMS = plyr::mapvalues(df_motions$TMS,from = c(1,2),to = c('sham','cTBS'))
df_motions$TMS = factor(df_motions$TMS)

write.xlsx(df_motions,
           file = file.path(sourceDir,'df_motions.xlsx'))

##################################################################
# summarize across subjects

p = NULL
for(i in 1:nParas){
Avgdf = df_motions %>% 
  subset(Paras == Para.names[i]) %>% 
  group_by(TMS,Time,Paras) %>%
  dplyr::summarise(mean = mean(value), sd = sd(value),
                   n = n(), sde = sd/sqrt(n))

p[[i]] = ggplot(Avgdf,aes(unfactor(Time))) + # notice here 'Time' can't be a factor
  geom_point(aes(y = mean, color = TMS), size = 3) +
  geom_line(aes(y = mean, color = TMS), linewidth = 1.5, show.legend = F) +
  geom_ribbon(aes(ymin = mean - sde, ymax = mean + sde, fill = TMS), alpha = 0.3) +
  scale_colour_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  labs(x="", y = "", title = '') +
  common + 
  geom_vline(xintercept=c(2.5,4.5), linetype="dotted", color = "black", size=0.5) +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank())
}

p[[1]] = p[[1]] +
  coord_cartesian(ylim = c(9.6, 12.2), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") 
  #labs(y = 'mm', title = 'Translation parameters') +
  #annotate(geom = "text", x = c(1.6,3.5,5.4), y = 9.45, 
  #         label = c("Run 1", "Run 2", "Run 3"), size = 6)

p[[2]] = p[[2]] +
  coord_cartesian(ylim = c(0.075, 0.13), xlim = c(0.6,6.4),
                  expand = FALSE, clip = "off") 
  #labs(y = 'degrees', title = 'Rotation parameters') +
  #annotate(geom = "text", x = c(1.6,3.5,5.4), y = 0.072, 
  #         label = c("Run 1", "Run 2", "Run 3"), size = 6)

fig = ggpubr::ggarrange(p[[1]], p[[2]], 
                        common.legend = T,
                        legend = 'none',
                        ncol = 2, nrow = 1, align = "hv")

pdf(NULL)
pdf(file.path(MRIPreDir,'Comp2Motion_runs.pdf'), 9, 4)
print(fig)
dev.off()



library(stats)
# perform a 2-way repeated-measure ANOVA
for(i in 1:nParas){
  print(Para.names[i])
  use.df = subset(df_motions,Paras == Para.names[i])
  anova_result = aov(value ~ TMS * Time + Error(Subs/(TMS * Time)), data = use.df)
  print(summary(anova_result))
  
  lme1 = lmer(value ~ unfactor(Time) + Subs + (1 | Subs), data = use.df)
  print(summary(lme1))
  anova(lme1)
  
  lme2 = lmer(value ~ Subs + (1 | Subs), data = use.df)
  print(summary(lme2))
  
  print('Compare LME models')
  print(anova(lme1,lme2))
  
}


##################################################################
# re-organize df to look at covariation of motion parameters 
# with global connectedness across subjects
##################################################################

vecs = matrix(NA,372,nParas)
for(i in 1:nParas){
reduced = df_motions %>%
  subset(Paras == Para.names[i]) %>%
  select(Subs,TMS,Time,value) %>%
  arrange(Subs,TMS,Time) 
  vecs[,i] = reduced$value
}

re_df_motion = data.frame(vecs)
names(re_df_motion) = Para.names
re_df_motion$Subs = reduced$Subs
re_df_motion$TMS = reduced$TMS
re_df_motion$Time = reduced$Time


##################################################################
# load similarly reorganized data from global connectedness
GC_SubsTMS = read.xlsx(file.path(GCDir,'GC_SubsTMS.xlsx'))

# check the rows are aligned correctly !!!
# combine motion & GC results
df = cbind(GC_SubsTMS,re_df_motion[,Para.names])

use.df = df %>%
  subset(Time %in% c(2)) %>%
  mutate(TMS = factor(TMS),
         Sub = factor(Sub))

# lmer1: baseline model with Sub, Time, and motion parameters
# lmer2: adding TMS to lmer1
lmer1_OFC <- lmer(OFC ~ Sub + Time + translation + rotation + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
lmer2_OFC <- lmer(OFC ~ Sub + Time + TMS + translation + rotation + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
anova(lmer1_OFC,lmer2_OFC) # is TMS useful after accounting for motions?
# p = 0.0001055

use.df = df %>%
  subset(Time %in% c(4)) %>%
  mutate(TMS = factor(TMS),
         Sub = factor(Sub))

# lmer1: baseline model with Sub and motion parameters
# lmer2: adding TMS to lmer1
lmer1_LPFC <- lmer(LPFC ~ Sub + Time + translation + rotation + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
lmer2_LPFC <- lmer(LPFC ~ Sub + Time + TMS + translation + rotation + (1 | Sub), 
              data = use.df, control = lmerControl(optimizer ="Nelder_Mead"))
anova(lmer1_LPFC,lmer2_LPFC) # is TMS useful after accounting for motions?
# p = 1.697e-6



