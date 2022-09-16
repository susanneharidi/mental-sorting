# Analysis of the RTs
# This script contains all Models described in the Behavioral results section

#########################################################################
# load packages
#########################################################################

library(readr)
library(tidyverse)
library(ggplot2)
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)
library(plyr)
library(dplyr)
library(data.table)
library(gridExtra)
library(ggthemes)
library(ggsignif)
library(brms)
library(lmerTest)
library(sjPlot)
library(insight)
library(httr)
library(brms)
library(here)

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project

#########################################################################
# functions
#########################################################################

#standard error function
se <- function(x){sd(x)/sqrt(length(x))}


#########################################################################
# load data 
#########################################################################

d <- read.csv(here("Data", "Experiment", "data_to_work_with2.csv"))

#if you are not using rproject or here packages, than use the lines below and 
#d<-read.csv("Data/Experiment/data_to_work_with.csv")

################################################################################
# Set some parameters
################################################################################

# set my labels and colors
myColors <- c("#a2f9ff","#00AFBB", "#ffe682", "#E7B800", "#fd8453", "#FC4E07")
mylabels <- c("None Memory", "None Sort", "Query Memory", "Query Sort", "Sequence Memory", "Sequence Sort")
alpha <- 0.7
pd <- position_dodge(0.1)
textsize <- 16

#color pallette for three different groups
cbbPalette <- c("grey50", "#0072B2", "#D55E00")

#other parameters in s
rt_cutoff <- 10
seed <- 2022

############################################################################
# determining the number of incorret trials and trials above cutoff
############################################################################

incorrect <- subset(d, correct == 0)
# 670 trials
proportionIncorrect = length(incorrect$rt)/length(d$rt)
# 0.04370515, so 4.37% of the data
incorrectAndLong <- subset(subset(d, correct == 1), rt > rt_cutoff | queryRT > rt_cutoff)
# 1216 trials
pCorrectLong = length(incorrectAndLong$rt)/length(subset(d, correct == 1)$rt)
#0.08328786


#########################################################################
# Modelling
#########################################################################


############################################################################
# rt analysis, looking at the effects of structure and the sequence length
############################################################################
# important!!!!! This is where all the models are saved, so we don't have to recalculate them every time
# so RUN THIS!
setwd(here("analyses", "brmModel"))
# setwd("./analyses/brmModel")
###################################################################
# Full RT Model with RT Type as interaction
###################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & otherRT <= rt_cutoff)

# The full Encoding RT model 
mrtallStructConditionRecallInteraction <- brm(rt ~ (Number_of_Bars + Structure + Condition) * Stimulus_type + ((Number_of_Bars + Structure + Condition + Block)*Stimulus_type|subject_id), 
                             data = dp,
                             iter = 10000,
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             file = "mrtallStructConditionNotLogRecallInteractionFinal")
# Model summary
summary(mrtallStructConditionRecallInteraction)
# Check for convergence
mcmc_plot(mrtallStructConditionRecallInteraction, type = "trace")
# Plot effects 
mcmc_plot(mrtallStructConditionRecallInteraction, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(mrtallStructConditionRecallInteraction, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", 
                              "Sequence Structure", "Sort Task", "Recall", 
                              "Seq Len*Recall", "QueryStruc*Recall", "SeqStruc*Recall",
                              "SortTask*Recall"))+# 
  ggtitle("Ecoding vs. Recall RT analysis")
pp_check(mrtallStructConditionRecallInteraction, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mrtallStructConditionRecallInteraction), points = FALSE)
tab_model(mrtallStructConditionRecallInteraction)

# Results 16.09.22
#                                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                               -0.21      0.15    -0.49     0.08 1.00     1396     3087
# Number_of_Bars                           0.89      0.05     0.80     0.99 1.00      799     1592
# StructureQuery                          -0.31      0.09    -0.48    -0.14 1.00     2161     4373
# StructureSequence                       -0.22      0.05    -0.31    -0.12 1.00     3210     7911
# ConditionSort                            0.40      0.06     0.28     0.53 1.00     2461     5063
# Stimulus_typequery                       0.18      0.16    -0.13     0.48 1.00     1572     3219
# Number_of_Bars:Stimulus_typequery       -0.51      0.05    -0.61    -0.40 1.00     1185     2370
# StructureQuery:Stimulus_typequery        0.18      0.06     0.05     0.31 1.00     5291     9608
# StructureSequence:Stimulus_typequery     0.21      0.05     0.12     0.30 1.00     8012    12061
# ConditionSort:Stimulus_typequery        -0.15      0.07    -0.29    -0.01 1.00     2793     6114

