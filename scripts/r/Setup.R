
print('Loading R packages ...')
library(scales)    # to use 'show_col' function
library(ggpubr)    # to use 'ggarrange' function
                   # also has a mutate function so load it early
library(plyr)      # load plyr before dplyr (within tidyverse)
                   # to use 'adply', mapvalues function
library(tidyverse) # load ggplot2, dplyr packages
library(openxlsx)  # open and save xlsx files
library(lme4)      # conduct linear-mixed effect models
library(stats)     # to use 'aov' function for Anova
library(R.matlab)  # load matlab mat files
library(varhandle) # to use 'unfactor' function
library(GGally)    # pair-wise corr plot working with ggplot

# some functions in multiple packages
# to use with package name to be safe
# dplyr::summarise()
# dplyr::mutate()

print('Setting up ggplot basics ...')
UseColor = c('#D35FB7','#159191') # colors for cTBS (magenta) and sham (green) respectively
#show_col(UseColor)
common = theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) 