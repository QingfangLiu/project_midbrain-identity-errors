
rm(list = ls())

source("~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/figs/ggplotSetup.R")

conn_name = 'ConnectednessMap_3mm'
Option.list = c('Signed','Unsigned','Squared')
Option = 'Unsigned'

#for(Option in Option.list){
  print(Option)

# where to save the result
saveDir = file.path(getwd(),'Summary',Option)
if (!dir.exists(saveDir)) {dir.create(saveDir,recursive = T)}

matname = paste0('Extracted_GC_values_',conn_name,'_',Option,'.mat')
matdat <- readMat(matname)
names(matdat)

MeanGC = matdat$MeanGC
ROIs = unlist(matdat$roi.labels) 
nROIs = length(ROIs)
nSess = 2
nruns = 3

# this is after excluding a few subjects
SubInfo = read_excel('../../SubjectConds.XLSX')
UseTMSConds = SubInfo$Cond[SubInfo$Excluded==0]
Subs = SubInfo$Sub[SubInfo$Excluded==0]
nSubs = length(Subs)

# convert this array to df
dfGC = structure(MeanGC,
                  .Dim = c(nROIs, nSess, nruns, nSubs), 
                  .Dimnames = structure(list(ROI = ROIs, 
                                             Sess = c(1,2),
                                             Run = c(1,2,3), 
                                             Subs = Subs), 
                  .Names = c("ROI", "Sess" ,"Run", "Subs")))
dfGC = adply(dfGC, c(1,2,3,4))
names(dfGC) = c(names(dfGC)[1:4],'GC')

# add TMS column
dfGC$TMS = varhandle::unfactor(dfGC$Sess)
for(j in 1:length(Subs)){
  if(UseTMSConds[j]=='AP'){
    dfGC$TMS[dfGC$Subs==Subs[j]] = 3 - dfGC$TMS[dfGC$Subs==Subs[j]]
  }
}

dfGC$TMS = plyr::mapvalues(dfGC$TMS,from = c(1,2),to = c('sham','cTBS'))
dfGC$TMS = factor(dfGC$TMS)

p = ggplot(data = dfGC,mapping = aes(x = Run, y = GC, fill = TMS)) + 
  geom_boxplot() +
  facet_wrap(~ ROI) +
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  labs(y = "Global Connectedess") + 
  common

pdf(file.path(saveDir,'GC_byROIs.pdf'),9,6)
print(p)
dev.off()

#########################################
# only look at bilateral ROIs
#########################################
ROIs = c('OFC','LPFC') # replace the ROIs variable here
dfGC = subset(dfGC, ROI %in% ROIs)

p = ggplot(data = dfGC, aes(x = Run, y = GC, fill = TMS)) + 
  geom_boxplot() +
  facet_wrap(~ ROI) +
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  labs(y = "Global Connectedess") + 
  common +
  theme(legend.position="none") 
pdf(file.path(saveDir,'GC_byROIs_bilateral.pdf'),12,6)
print(p)
dev.off()

p = NULL
p[[1]] = ggplot(data = subset(dfGC, Run == 1), aes(x = ROI, y = GC, fill = TMS)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  labs(y = "Global Connectedess") + 
  common + 
  ggtitle("Run 1 only") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") 
p[[2]] = ggplot(data = dfGC, mapping = aes(fill = TMS, y = GC, x = ROI)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  labs(y = "Global Connectedess") + 
  common +
  ggtitle("Across 3 runs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") 

pdf(file.path(saveDir,'GC_byROIs_bilateral_summary.pdf'),12,6)
par(mfrow=c(1,2),cex.lab=1.8,mar=c(5,5,3,2),cex.axis=1.5,cex.main=2)
ggpubr::ggarrange(p[[1]], p[[2]], ncol = 2, nrow = 1, widths = c(1,1))
dev.off()

# do some simple t-test
print('Compare GC between TMS conditions on run 1 only')
for (ROI in ROIs){
  print(ROI)
  x = dfGC$GC[dfGC$ROI==ROI & dfGC$Run == 1 & dfGC$TMS == 'cTBS']
  y = dfGC$GC[dfGC$ROI==ROI & dfGC$Run == 1 & dfGC$TMS == 'sham']
  print(t.test(x,y,paired = T, alternative = 'l'))
}

print('Compare GC between TMS conditions across 3 runs')
for (ROI in ROIs){
  print(ROI)
  x = dfGC$GC[dfGC$ROI==ROI & dfGC$TMS == 'cTBS']
  y = dfGC$GC[dfGC$ROI==ROI & dfGC$TMS == 'sham']
  print(t.test(x,y,paired = T, alternative = 'l'))
}

#########################################
# add linear mixed effect models to test TMS effect
#########################################

library(lme4)
for (r in ROIs){
  for(i in 1:3){
    print(paste(r,'Run',i))
    df = subset(dfGC, ROI==r & Run==i)
    lmer1 <- lmer(GC ~ TMS + Subs + Sess + (1 | Subs), data = df)
    lmer2 <- lmer(GC ~ Subs + Sess + (1 | Subs), data = df)
    print(anova(lmer1,lmer2))
}
}
#########################################
# simply another way to look at the data
# calculate difference between sham and cTBS
#########################################

dfGCdiff = dfGC$GC[dfGC$TMS == 'sham'] - dfGC$GC[dfGC$TMS == 'cTBS']
df_unique = dfGC[dfGC$TMS == 'sham',]
dfGCdiff = cbind(df_unique,dfGCdiff)
dfGCdiff <- subset(dfGCdiff, select = -c(Sess, TMS, GC))

p = NULL
for(i in 1:length(ROIs)){
  p[[i]] = ggplot(subset(dfGCdiff,ROI == ROIs[i]), aes(x = Run, y = dfGCdiff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.4, width = 0.2) + 
  labs(y = '', x = '') + common +
  scale_x_discrete(labels=c("1"="Run 1","2"="Run 2","3"="Run 3")) +
  geom_hline(yintercept=0, linetype="dotted", color = "black", size=0.5) +
  theme(strip.text.x = element_blank())
}

pdf(file.path(saveDir,'GC_byROIs_bilateral_Change.pdf'),4,9)
par(mfrow=c(1,2),cex.lab=1.8,mar=c(5,5,3,2),cex.axis=1.5,cex.main=2)
fig = ggpubr::ggarrange(p[[1]], p[[2]], ncol = 1, nrow = 2, align = "hv")
print(fig)
dev.off()

# compare TMS effect on different orders
group1 = SubInfo$Cond[SubInfo$Excluded==0]=='PA'
group2 = SubInfo$Cond[SubInfo$Excluded==0]=='AP'

group1 = SubInfo$Sex[SubInfo$Excluded==0]=='F'
group2 = SubInfo$Sex[SubInfo$Excluded==0]=='M'

## stat test
print('Compare GC difference at each run')
for (r in ROIs){
  for(i in 1:3){
    print(paste(r,'Run',i))
    x = dfGCdiff$dfGCdiff[dfGCdiff$ROI == r & dfGCdiff$Run == i]
    #print(mean(x))
    print(t.test(x,alternative = 'g'))
    print(wilcox.test(x,mu=0,alternative = 'g'))
    #print(t.test(x[group1],x[group2]))
    #print(wilcox.test(x[group1],x[group2]))
  }
}


# convert run from a factor to numeric variable
dfGCdiff$Run = varhandle::unfactor(dfGCdiff$Run)
model <- lm(dfGCdiff ~ Run, data = subset(dfGCdiff,ROI == 'OFC'))
summary(model)

model <- lm(dfGCdiff ~ Run, data = subset(dfGCdiff,ROI == 'LPFC'))
summary(model)

#########################################
# correlate with behavioral and neural effects of TMS
#########################################

beh_post_rev_acc = read_csv('../../BehAnalyzeRes/SubPostAcc.csv')
comp_post_rev_acc = beh_post_rev_acc$Acc[beh_post_rev_acc$TMS=='Sham'] -
  beh_post_rev_acc$Acc[beh_post_rev_acc$TMS=='Stim']

alpha_beh = read.csv('../../BehModeling/PPD/v4az/alpha_map.csv')
logit = function(p) log(p/(1-p))
comp_alpha_beh = logit(alpha_beh$V1) - logit(alpha_beh$V2)

# compare global conn diff with behavioral diff
for(r in ROIs){
  print(r)
  SubGCdiff = filter(dfGCdiff, Run==1 & ROI == r)
  
  pdf(file.path(saveDir,paste0('GC_beh_corr_',r,'.pdf')),10,5)
  par(mfrow=c(1,2))
  plot(comp_alpha_beh,SubGCdiff$dfGCdiff,
       xlab = 'Learning rate (sham - cTBS)',
       ylab = 'Global Connect (sham - cTBS)')
  print(cor.test(comp_alpha_beh,SubGCdiff$dfGCdiff))
  
  plot(comp_post_rev_acc,SubGCdiff$dfGCdiff,
       xlab = 'rev+1 acc (sham - cTBS)',
       ylab = 'Global Connect (sham - cTBS)')
  print(cor.test(comp_post_rev_acc,SubGCdiff$dfGCdiff))
  dev.off()
}

# neural
SubGCdiff = dfGCdiff %>% group_by(ROI, Subs) %>%
  summarise(mean = mean(dfGCdiff))
colnames(SubGCdiff)[3] = 'GCdiff'

# save this as a subject-level of GC change (sham - cTBS)
# of each ROI
ToSave = subset(SubGCdiff,ROI == 'OFC')
ToSave = ToSave[,-1]
write.table(ToSave, 
            file.path(saveDir,"SubGCdiff_Sham_minus_cTBS_OFC.txt"),
            sep="\t",row.names=FALSE)

ToSave = subset(SubGCdiff,ROI == 'LPFC')
ToSave = ToSave[,-1]
write.table(ToSave, 
            file.path(saveDir,"SubGCdiff_Sham_minus_cTBS_LPFC.txt"),
            sep="\t",row.names=FALSE)



neural_post_rev_beta = read.table('../../NeuralAnalyzeRes/BetaDecreasePostRev_SubsTMS.txt',
                                  header = T)

for(r in ROIs){
  print(r)
  SubGCdiffROI = filter(SubGCdiff, ROI == r)$GCdiff
  pdf(file.path(saveDir,paste0('GC_neural_corr_',r,'.pdf')),10,5)
  par(mfrow=c(1,2))
  TMS_LPFC = neural_post_rev_beta$LPFC.Sham - neural_post_rev_beta$LPFC.cTBS
  plot(TMS_LPFC,SubGCdiffROI,
       xlab = 'TMS effect on LPFC (sham - cTBS)',
       ylab = 'Global Connect (sham - cTBS)')
  cor.test(TMS_LPFC,SubGCdiffROI)
  
  TMS_MB = neural_post_rev_beta$MB.Sham - neural_post_rev_beta$MB.cTBS
  plot(TMS_MB,SubGCdiffROI,
       xlab = 'TMS effect on MB (sham - cTBS)',
       ylab = 'Global Connect (sham - cTBS)')
  cor.test(TMS_MB,SubGCdiffROI)
  dev.off()
}

#}


