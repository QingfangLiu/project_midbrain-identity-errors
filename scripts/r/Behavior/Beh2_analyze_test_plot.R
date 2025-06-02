
rm(list = ls())

# dir of beh data
ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
Behdir = file.path(ResDir,'Behavior')
modeldir = file.path(ResDir,'BehavioralModeling/JagsPosteriors/diff')
  
source("../Setup.R")

datinfo = read.xlsx(file.path(Behdir,'SubjectConds.xlsx'))
subs = subset(datinfo,Excluded==0)$Sub 
nsubs = length(subs)

AllRevdf = read.xlsx(file.path(Behdir,'OdorPredRes_processed.xlsx'))
SubRevAcc = AllRevdf %>%
  subset(Run == 1) %>%
  group_by(Sub, RevLoc, TMS) %>%
  dplyr::summarize(Acc = mean(Acc)) %>%
  ungroup()

# average across subjects 
AvgRevAcc = SubRevAcc %>% 
  group_by(RevLoc, TMS) %>%
  mutate(RevLoc = factor(RevLoc),TMS = factor(TMS)) %>%
  dplyr::summarize(mean = mean(Acc),sde = sd(Acc)/sqrt(n())) %>% 
  ungroup()

################################################
## start plotting
p = NULL

################################################
# panel a: plotting error bar plot of acc across trials

p[[1]] = AvgRevAcc %>%
  ggplot(aes(x = RevLoc, y = mean, group = TMS)) + 
  geom_errorbar(aes(ymin=mean-sde, ymax=mean+sde,color=TMS), 
                width=.2, linewidth=1.1) +
  geom_line(aes(color=TMS),linewidth=1.1)+
  geom_point(aes(color=TMS)) +
  scale_color_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) + 
  scale_x_discrete(labels=c("-1"=expression(rev[-1]), "0"="rev",
                            "1"=expression(rev[+1]),"2"=expression(rev[+2]))) + 
  common + 
  theme(legend.position = "none") +
  labs(x="",y="") + ylim(0,1)


###########################################
# change of accuracy from pre-reversal to post-reversal trial
###########################################

Acc_post_sham = subset(SubRevAcc, TMS == 'sham' & RevLoc == 1)$Acc
Acc_post_cTBS = subset(SubRevAcc, TMS == 'cTBS' & RevLoc == 1)$Acc
Acc_pre_sham = subset(SubRevAcc, TMS == 'sham' & RevLoc == -1)$Acc
Acc_pre_cTBS = subset(SubRevAcc, TMS == 'cTBS' & RevLoc == -1)$Acc

t.test(Acc_post_cTBS,Acc_post_sham,alternative = 't',paired = T)
# t = -1.6067, df = 30, p-value = 0.1186
t.test(Acc_pre_cTBS,Acc_pre_sham,alternative = 't',paired = T)
# t = 1.4756, df = 30, p-value = 0.1505

# acc decreased from pre to post rev in cTBS
t.test(Acc_pre_cTBS,Acc_post_cTBS,paired = T, alternative = 'g')
# t = 3.0134, df = 30, p-value = 0.002605

# acc did not decrease from pre to post rev in sham
t.test(Acc_pre_sham,Acc_post_sham,paired = T, alternative = 'g')
# t = -0.3189, df = 30, p-value = 0.624

# change of acc from pre to post reversal trials (interaction test)
change_sham = Acc_post_sham - Acc_pre_sham
change_cTBS = Acc_post_cTBS - Acc_pre_cTBS

t.test(change_sham,change_cTBS,alternative = 'g',paired = T)
# t = 2.3442, df = 30, p-value = 0.01294

beh_signature_df_SubsTMS = data.frame('Sub'=rep(subs,2),
                   'beh_signature'=c(change_sham,change_cTBS),
                   'TMS'=rep(c('sham','cTBS'),each=nsubs))


df2 = data.frame(Sub = c(subs,subs),
                 TMS = c(rep('sham',nsubs),rep('cTBS',nsubs)),
                 Change = c(change_sham,change_cTBS))
df2 = df2 %>% 
  mutate(Sub = factor(Sub), TMS = factor(TMS))

# plot the change of accuracy
p[[2]] = ggplot(df2, aes(x = TMS, y = Change, group = TMS)) +
  geom_boxplot(aes(fill = TMS),alpha = 0.5,outlier.alpha = 0) +
  geom_line(aes(group= Sub), position = position_dodge(0.2), size = 0.2) +
  geom_point(aes(fill=TMS,group=Sub),size=2,shape=21,alpha=1,
             position = position_dodge(0.2)) +
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  common +
  theme(legend.position="none") +
  labs(x="",y="") 

# add line segments and *
high_y1 = 0.57
low_y1 = 0.52
p[[2]] = p[[2]] +
  geom_segment(aes(x=1,y=high_y1,xend=2,yend=high_y1),size = 0.2) +
  geom_segment(aes(x=1,y=low_y1,xend=1,yend=high_y1),size = 0.2) +
  geom_segment(aes(x=2,y=low_y1,xend=2,yend=high_y1),size = 0.2) +
  annotate(geom = "text", x = c(1.5), y = c(0.59), label = '*', size = 6)

################################################

AggSimudf = readRDS(file.path(modeldir,'SimuDat_MAP.rds'))
samples = readRDS(file.path(modeldir,'samples.rds'))
alpha_beh = read.xlsx(file.path(modeldir,'alpha_map.xlsx'))

################################################
# panel c: plotting observed acc against model simulated acc

p[[3]] = AggSimudf %>%
  mutate(TMS = factor(TMS,levels = c('cTBS','sham'))) %>%
  ggplot(aes(x=RevLoc, y=Acc, fill=TMS)) + 
  geom_boxplot(alpha=0.5, outlier.shape = NA, show.legend=F) + 
  scale_fill_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  geom_point(data = AvgRevAcc, 
             position = position_dodge(width=0.75),
             aes(x = RevLoc, y = mean, group = TMS, color = TMS),
             show.legend = F) +
  geom_errorbar(data=AvgRevAcc, 
                position=position_dodge(width=0.75),
                aes(x = RevLoc, ymin=mean-sde, ymax=mean+sde, y = mean,
                    color = TMS), width=.2,show.legend = F) +
  scale_x_discrete(labels=c("-1"=expression(rev[-1]), "0"="rev",
                            "1"=expression(rev[+1]),"2"=expression(rev[+2]))) +
  scale_color_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) + 
  common +
  labs(x="",y="") 

################################################
# panel d: posterior estimation of hyper learning rates
################################################

mu <- samples$BUGSoutput$sims.list$mu
nsamples = nrow(mu)
x = as.numeric(mu)
group = c(rep("sham", nsamples), rep("cTBS", nsamples))
df_mu <- data.frame(x, group)

prior_x = seq(0.001,0.999,0.001)
prior_y = dbeta(prior_x,8,2)
df_prior = data.frame(prior_x,prior_y)

p[[4]] = ggplot(df_prior)+
  geom_line(aes(prior_x,prior_y),colour = 'gray',linetype = "dotted", size = 1.5)+
  geom_density(data = df_mu, aes(x = x, colour = group),size = 1.5) +
  scale_color_manual(values = c("sham" = UseColor[2],"cTBS" = UseColor[1])) +
  common + 
  lims(x = c(0.6,1), y = c(0,9.5)) +
  labs(x="",y="") +
  theme(legend.position="none") 

################################################
# panel e: correlation between alpha difference and post-acc diff
# load alpha estimates from behavioral data

logit = function(p) log(p/(1-p))
comp_alpha_beh = logit(alpha_beh$sham) - logit(alpha_beh$cTBS)

# calculate acc change (rev+1 - rev-1)
# using sham - cTBS on the acc change as a measure of behavioral TMS effect

tmp_acc_change_sham = subset(SubRevAcc,RevLoc==1 & TMS=='sham')$Acc - 
  subset(SubRevAcc,RevLoc==-1 & TMS=='sham')$Acc
tmp_acc_change_cTBS = subset(SubRevAcc,RevLoc==1 & TMS=='cTBS')$Acc - 
  subset(SubRevAcc,RevLoc==-1 & TMS=='cTBS')$Acc
tmp_acc_change_comp = tmp_acc_change_sham - tmp_acc_change_cTBS
tmpdf = data.frame('Sub'=subs,
                   'alpha_diff'=comp_alpha_beh,
                   'acc_diff'=tmp_acc_change_comp)

p[[5]] = ggplot(tmpdf,aes(alpha_diff,acc_diff)) +
  geom_point(size=5, shape=16, alpha = 0.6, color = 'black') +
  geom_smooth(method=lm, se=FALSE, linetype="dotted",color = "gray",size = 1.5) +
  ggpubr::stat_cor(label.sep='\n',label.x.npc = 0.6,label.y.npc = 0.2,
           r.accuracy = 0.01, method = 'spearman') +
  common +
  labs(x="",y="")


## plotting!
pdf(file.path(Behdir,'Behavior.pdf'), 8, 12)
fig = ggpubr::ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], NULL,
                  ncol = 2, nrow = 3,
                  heights = c(1, 1), 
                  widths = c(1,1,1))
print(fig)
dev.off()

pdf(file.path(Behdir,'Behavior_reduced.pdf'), 8, 8)
fig = ggpubr::ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], 
                        ncol = 2, nrow = 2,
                        heights = c(1,1), 
                        widths = c(1,1))
print(fig)
dev.off()


# save individual rev performance for future use
write.xlsx(SubRevAcc, file.path(Behdir,'SubRevAcc.xlsx'))
write.xlsx(tmpdf,file.path(Behdir,'BehavioralTMSeffect.xlsx'))

write.xlsx(beh_signature_df_SubsTMS,
           file.path(Behdir,'beh_signature_df_SubsTMS.xlsx'))

