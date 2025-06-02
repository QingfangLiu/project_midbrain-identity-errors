
rm(list = ls())
SCdat = read.csv('ScreeningSubOdors.csv')

# prepare for plotting
source("~/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/figs/ggplotSetup.R")

# 3 plea ratings are corresponding to 3 odors
# same for the intensity ratings

# what are the scales of the ratings?
# pleasant: -10 ~ 10 (10 is the most pleasant)
# intensity: -10 ~ 10 (10 is the most strong)
# discrimination: accuracy (1 is correct)
# similarity: -10 ~ 10 (10 is identical)

# what to report?
nsubs = nrow(SCdat)
plea_mean = sapply(1:nsubs,function(x) mean(as.numeric(SCdat[x,5:7])))
plea_sd = sapply(1:nsubs,function(x) sd(as.numeric(SCdat[x,5:7])))

inten_mean = sapply(1:nsubs,function(x) mean(as.numeric(SCdat[x,8:10])))
inten_sd = sapply(1:nsubs,function(x) sd(as.numeric(SCdat[x,8:10])))

df = data.frame(plea_mean = plea_mean,
                plea_sd = plea_sd,
                inten_mean = inten_mean,
                inten_sd = inten_sd,
                disc = SCdat$disc,
                simi = SCdat$simi)
#tt = geom_histogram(aes(y=..density..),bins = 8)
#tt = geom_histogram(binwidth = 3)
  
p = NULL
library(ggplot2)
p[[1]] = ggplot(df,aes(x=plea_mean)) + labs(x ="Average pleasantness") + 
         common + geom_histogram(binwidth = 1)
p[[2]] = ggplot(df,aes(x=plea_sd)) + labs(x ="SD. (pleasantness)") + 
         common + geom_histogram(binwidth = 0.5)
p[[3]] = ggplot(df,aes(x=inten_mean)) + labs(x ="Average intensity") + 
         common + geom_histogram(binwidth = 1)
p[[4]] = ggplot(df,aes(x=inten_sd)) + labs(x ="SD. (intensity)") + 
         common + geom_histogram(binwidth = 0.5)
p[[5]] = ggplot(df,aes(x=disc)) + labs(x ="Discrimination score") + 
         common + geom_histogram(binwidth = 0.05)
p[[6]] = ggplot(df,aes(x=simi)) + labs(x ="Similarity") + 
         common + geom_histogram(binwidth = 1)

#pdf('../figs/finals/FigS1_Odors_screening.pdf', 9, 12)
pdf('FigS1_Odors_screening.pdf', 9, 12)
ggpubr::ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]],
                  ncol = 2, nrow = 3, 
                  #labels = c('a','','b','','c','d'),
                  align = "v", heights = c(1, 1, 1), widths = c(1,1),
                  font.label = list(size = 24))
dev.off()

