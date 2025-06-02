
rm(list = ls())
library(openxlsx)
library(stats)
library(dplyr)

ResDir = '/Volumes/QF10TB/SPEEDTMS_results'
MVPADir = file.path(ResDir,'MVPA_Searchlight')

df = read.xlsx(file.path(MVPADir,'TrialOnsetDelay.xlsx'))
df$TMS = factor(df$TMS)
df$PrePost = factor(df$PrePost)

df %>%
  group_by(TMS,PrePost) %>%
  dplyr::summarise(mean = mean(Delay),
                   sd = sd(Delay))

# perform a 2-way repeated-measure ANOVA
anova_result <- aov(Delay ~ TMS * PrePost + 
                      Error(Subs/(TMS * PrePost)), data = df)
summary(anova_result)

# see blog post
# http://agroninfotech.blogspot.com/2020/06/two-way-repeated-measures-analysis-in-r.html