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

d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))

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
# Full encoding RT model
###################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# The full Encoding RT model 
mrtallStructCondition <- brm(rt ~ Number_of_Bars + Structure + Condition + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                             data = dp,
                             iter = 10000,
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             family = exgaussian(),
                             file = "mrtallStructConditionNotLogFinalEx")
# Model summary
summary(mrtallStructCondition)
# Check for convergence
mcmc_plot(mrtallStructCondition, type = "trace")
# Plot effects 
mcmc_plot(mrtallStructCondition, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(mrtallStructCondition, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Raw RT analysis")
pp_check(mrtallStructCondition, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mrtallStructCondition), points = FALSE)
tab_model(mrtallStructCondition, transform = NULL)


# results exgaussian 23.09.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.60      0.11     0.38     0.82 1.00     2362     4347
# Number_of_Bars        0.66      0.04     0.57     0.75 1.00     1410     2866
# StructureQuery       -0.26      0.07    -0.41    -0.12 1.00     3283     6181
# StructureSequence    -0.17      0.04    -0.25    -0.09 1.00     3355     8744
# ConditionSort         0.26      0.05     0.17     0.36 1.00     4700     8819

##############
# Plot the Residuals
###############
residulas <- residuals(mrtallStructCondition)
dp$Residuals <- residulas[,1]

#get summarize data
dRTSummary <- ddply(dp, ~Number_of_Bars + Condition + Structure + Stimulus_type, summarize, mu = mean(Residuals), se = se(Residuals))

#start plotting the residuals
pResiduals <- ggplot(dRTSummary, aes(x = Number_of_Bars, y = mu, col = Structure)) +
  #points
  geom_point(position = pd)+
  #error bars
  geom_errorbar(aes(ymin = mu-se, ymax = mu+se), width = 0, size = 1, position = pd) +
  #lines
  geom_line(position = pd, size = 1.2) +
  scale_color_manual(values = cbbPalette)+
  #legend on top
  theme(strip.background = element_blank())+
  #change x ticks
  scale_x_continuous(breaks = round(seq(min(0), max(10), by = 1),1)) +
  #change y ticks
  scale_y_continuous(breaks = round(seq(min(0), max(10), by = 1),1)) +
  #minimal theme
  theme_minimal()+
  #add xlab
  xlab("Sequence Length")+
  #add ylab
  ylab('Mean Residualsin s (±SE)')+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"))+
  #legend on top
  theme(legend.position = "top")+
  facet_grid(cols = vars(Condition))

#show!
pResiduals


###################################################################
# Full encoding RT model With Interactions 
###################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# The full Encoding RT model 
mrtallStructConditionInterExGaussian <- brm(rt ~ Number_of_Bars*(Structure + Condition) + (Number_of_Bars*(Structure + Condition) + Block|subject_id), 
                                  data = dp,
                                  iter = 10000,
                                  cores = 4,
                                  save_pars = save_pars(all = TRUE),
                                  seed = seed,
                                  family = exgaussian(),
                                  file = "mrtallStructConditionInterExGaussianFinal")
# Model summary
summary(mrtallStructConditionInterExGaussian)
# Check for convergence
mcmc_plot(mrtallStructConditionInterExGaussian, type = "trace")
# Plot effects 
mcmc_plot(mrtallStructConditionInterExGaussian, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(mrtallStructConditionInterExGaussian, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task", "Seq length*Query Struc", "Seq length*Seq Struc", "Seq length*Sort Task"))+# 
  ggtitle("Raw RT analysis")
pp_check(mrtallStructConditionInterExGaussian, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mrtallStructConditionInterExGaussian), points = FALSE)
tab_model(mrtallStructConditionInterExGaussian, transform = NULL)

# results 20.09.2022
# Intercept                            0.63      0.11     0.42     0.83 1.00     2292     4677
# Number_of_Bars                       0.66      0.05     0.56     0.75 1.00     1051     2574
# StructureQuery                       0.31      0.07     0.16     0.45 1.00     2826     6701
# StructureSequence                   -0.06      0.05    -0.16     0.04 1.00     6627    12078
# ConditionSort                       -0.19      0.05    -0.29    -0.09 1.00     5454    11549
# Number_of_Bars:StructureQuery       -0.16      0.04    -0.24    -0.08 1.00     2809     6619
# Number_of_Bars:StructureSequence    -0.04      0.02    -0.08    -0.00 1.00     4952     9750
# Number_of_Bars:ConditionSort         0.14      0.02     0.09     0.18 1.00     5160    10137


##############
# Plot the Residuals of Model with Intercations
###############
residulas <- residuals(mrtallStructConditionInterExGaussian)
dp$Residuals <- residulas[,1]

#get summarize data
dRTSummary <- ddply(dp, ~Number_of_Bars + Condition + Structure + Stimulus_type, summarize, mu = mean(Residuals), se = se(Residuals))

#start plotting the residuals
pResiduals <- ggplot(dRTSummary, aes(x = Number_of_Bars, y = mu, col = Structure)) +
  #points
  geom_point(position = pd)+
  #error bars
  geom_errorbar(aes(ymin = mu-se, ymax = mu+se), width = 0, size = 1, position = pd) +
  #lines
  geom_line(position = pd, size = 1.2) +
  scale_color_manual(values = cbbPalette)+
  #legend on top
  theme(strip.background = element_blank())+
  #change x ticks
  scale_x_continuous(breaks = round(seq(min(0), max(10), by = 1),1)) +
  #change y ticks
  scale_y_continuous(breaks = round(seq(min(0), max(10), by = 1),1)) +
  #minimal theme
  theme_minimal()+
  #add xlab
  xlab("Sequence Length")+
  #add ylab
  ylab('Mean Residualsin s (±SE)')+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"))+
  #legend on top
  theme(legend.position = "top")+
  facet_grid(cols = vars(Condition))

#show!
pResiduals


###################################################################
# Checking the effect of Block Number for encoding RT
###################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# The full Encoding RT model plus block effects
mrtallStructConditionBlock <- brm(rt ~ Number_of_Bars + Structure + Condition + Block + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                                  data = dp,
                                  iter = 10000,
                                  cores = 4,
                                  save_pars = save_pars(all = TRUE),
                                  seed = seed,
                                  family = exgaussian(),
                                  file = "mrtallStructConditionNotLogFinalBlockEx")
# Model summary
summary(mrtallStructConditionBlock)
# Check for convergence
mcmc_plot(mrtallStructConditionBlock, type = "trace")
# Plot effects 
mcmc_plot(mrtallStructConditionBlock, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(mrtallStructConditionBlock, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Raw RT analysis")
pp_check(mrtallStructConditionBlock, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mrtallStructConditionBlock), points = FALSE)


# results Exgaussian 23.09.2022
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             1.04      0.10     0.84     1.25 1.00     2491     4598
# Number_of_Bars        0.66      0.03     0.60     0.73 1.00     1647     3318
# StructureQuery       -0.33      0.05    -0.44    -0.22 1.00     4282     7589
# StructureSequence    -0.16      0.03    -0.22    -0.11 1.00     8921    12676
# ConditionSort         0.25      0.04     0.18     0.32 1.00     6035    10085
# Block                -0.09      0.01    -0.12    -0.07 1.00     4416     7858


###########################################################################
# Encoding RT testing for the importance of Sequence length, Structure and Memory vs. Sort
###########################################################################

#################################################################################
# testing for Sequence length vs no Sequence length
#################################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# regression
mrtnoNoB <- brm(rt ~ Structure + Condition + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                data = dp,
                iter = 10000,
                cores = 4,
                save_pars = save_pars(all = TRUE),
                control = list(max_treedepth = 15, adapt_delta = 0.99),
                seed = seed,
                family = exgaussian(),
                file = "mrtnoNoBNotLogFinalEx")
summary(mrtnoNoB)

# exgaussian 23.09.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             1.45      0.18     1.11     1.82 1.01      563     1245
# StructureQuery       -0.42      0.12    -0.67    -0.18 1.00     3810     6190
# StructureSequence    -0.09      0.06    -0.20     0.05 1.00     5072     9363
# ConditionSort         0.17      0.08     0.01     0.33 1.00     2884     7765


bayes_factor(mrtallStructCondition, mrtnoNoB)
# exgaussian 23.09.2022
# Warning messages:
#   1: logml could not be estimated within maxiter, rerunning with adjusted starting value. 
# Estimate might be more variable than usual. 
# 2: logml could not be estimated within maxiter, rerunning with adjusted starting value. 
# Estimate might be more variable than usual. 
# BF = 6191194495895217918082062666.00000

#################################################################################
# testing for Structure vs no Structure
#################################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# regression
mrtnoStructure <- brm(rt ~ Number_of_Bars + Condition + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                      data = dp,
                      iter = 10000,
                      cores = 4,
                      save_pars = save_pars(all = TRUE),
                      control = list(max_treedepth = 15, adapt_delta = 0.99),
                      seed = seed,
                      family = exgaussian(),
                      file = "mrtnoStructureNotLogFinalEx")

summary(mrtnoStructure)
# exgaussian 23.09.2022
#                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept          0.48      0.12     0.25     0.71 1.00     2171     3598
# Number_of_Bars     0.66      0.05     0.57     0.75 1.00     1323     2514
# ConditionSort      0.26      0.05     0.16     0.37 1.00     3965     7872

bayes_factor(mrtallStructCondition, mrtnoStructure)
# exgaussian 23.09.2022
# BF =  395.17883
# Warning messages:
#   1: logml could not be estimated within maxiter, rerunning with adjusted starting value. 
# Estimate might be more variable than usual. 
# 2: logml could not be estimated within maxiter, rerunning with adjusted starting value. 
# Estimate might be more variable than usual. 

#################################################################################
# testing for task Condition vs no task Condition
#################################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# regression
mrtnoCondition <- brm(rt ~ Number_of_Bars + Structure + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                      data = dp,
                      iter = 10000,
                      cores = 4,
                      save_pars = save_pars(all = TRUE),
                      control = list(max_treedepth = 15, adapt_delta = 0.99),
                      seed = seed,
                      family = exgaussian(),
                      file = "mrtnoConditionNotLogFinalEx")

summary(mrtnoCondition)

# exgaussian 23.09.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.90      0.12     0.67     1.13 1.00     1216     2667
# Number_of_Bars        0.62      0.05     0.52     0.72 1.01      659     1664
# StructureQuery       -0.31      0.08    -0.47    -0.15 1.00     1428     3290
# StructureSequence    -0.16      0.05    -0.25    -0.08 1.00     1189     3552

bayes_factor(mrtallStructCondition, mrtnoCondition)
# exgaussian 23.09.2022
# bf = 222100.09266
# Warning messages:
#  1: logml could not be estimated within maxiter, rerunning with adjusted starting value. 
# Estimate might be more variable than usual. 
# 2: logml could not be estimated within maxiter, rerunning with adjusted starting value. 
# Estimate might be more variable than usual.

###################################################################
# Full RT Model with RT Type as interaction
###################################################################
d2 <- read.csv(here("Data", "Experiment", "data_to_work_with2.csv"))

#create data set as before
dp <- subset(d2, rt <= rt_cutoff & correct == 1 & otherRT <= rt_cutoff)

# The full Encoding RT model 
mrtallStructConditionRecallInteraction <- brm(rt ~ (Number_of_Bars + Structure + Condition) * Stimulus_type + ((Number_of_Bars + Structure + Condition + Block)*Stimulus_type|subject_id), 
                                              data = dp,
                                              iter = 10000,
                                              cores = 4,
                                              save_pars = save_pars(all = TRUE),
                                              seed = seed,
                                              family = exgaussian(),
                                              file = "mrtallStructConditionNotLogRecallInteractionFinalEx")
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

# Results Exgaussian 23.09.2022
#                                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                                0.60      0.09     0.42     0.77 1.00     2836     5329
# Number_of_Bars                           0.60      0.04     0.52     0.67 1.00     1351     3047
# StructureQuery                          -0.26      0.06    -0.39    -0.14 1.00     3605     6720
# StructureSequence                       -0.16      0.03    -0.22    -0.11 1.00     5492     9719
# ConditionSort                            0.23      0.04     0.15     0.31 1.00     5517    10222
# Stimulus_typequery                       0.26      0.10     0.07     0.45 1.00     3340     6273
# Number_of_Bars:Stimulus_typequery       -0.38      0.04    -0.45    -0.31 1.00     1783     3859
# StructureQuery:Stimulus_typequery        0.17      0.06     0.06     0.28 1.00     6899    11675
# StructureSequence:Stimulus_typequery     0.14      0.03     0.09     0.20 1.00    10989    12719
# ConditionSort:Stimulus_typequery        -0.11      0.04    -0.19    -0.03 1.00     6221    10252


###################################################################
# Summed RT Model
###################################################################

#create data set as before
dp <- subset(d2, rt <= rt_cutoff & correct == 1 & otherRT <= rt_cutoff & Stimulus_type == "bars")
dp$summedRT <- dp$rt + dp$otherRT

# The full Encoding RT model 
mrtSummedallStructConditionNotLogFinalEx <- brm(summedRT ~ Number_of_Bars + Structure + Condition + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                                              data = dp,
                                              iter = 10000,
                                              cores = 4,
                                              save_pars = save_pars(all = TRUE),
                                              seed = seed,
                                              control = list(max_treedepth = 15, adapt_delta = 0.99),
                                              family = exgaussian(),
                                              file = "mrtSummedallStructConditionNotLogFinalExControl")
# Model summary
summary(mrtSummedallStructConditionNotLogFinalEx)
# Check for convergence
mcmc_plot(mrtSummedallStructConditionNotLogFinalEx, type = "trace")
# Plot effects 
mcmc_plot(mrtSummedallStructConditionNotLogFinalEx, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(mrtSummedallStructConditionNotLogFinalEx, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Ecoding vs. Recall RT analysis")
pp_check(mrtSummedallStructConditionNotLogFinalEx, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mrtSummedallStructConditionNotLogFinalEx), points = FALSE)
tab_model(mrtSummedallStructConditionNotLogFinalEx)


# Results Exgaussian 16.03.2023 
#                      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.75      0.13     0.51     1.01 1.00     2649     4793
# Number_of_Bars        0.99      0.05     0.88     1.10 1.00     1043     2243
# StructureQuery       -0.42      0.12    -0.66    -0.18 1.00     2474     4728
# StructureSequence    -0.22      0.05    -0.33    -0.12 1.00     4120     9611
# ConditionSort         0.44      0.08     0.29     0.59 1.00     3860     7749



##################################################################################
# the following model does not converge!!!
###################################################################
# Full RT Model with RT Type and seq length as interaction
###################################################################

#create data set as before
dp <- subset(d2, rt <= rt_cutoff & correct == 1 & otherRT <= rt_cutoff)

# The full Encoding RT model 
mrtallRecallSeqLenInteraction <- brm(rt ~ Number_of_Bars * Stimulus_type * (Structure + Condition) +
                                       ((Number_of_Bars + Stimulus_type + Structure + Condition + Block)|subject_id), 
                                     data = dp,
                                     iter = 10000,
                                     cores = 4,
                                     save_pars = save_pars(all = TRUE),
                                     seed = seed,
                                     family = exgaussian(),
                                     file = "mrtallRecallSeqLenInteractionFinal")
# Model summary
summary(mrtallRecallSeqLenInteraction)
# Check for convergence
mcmc_plot(mrtallRecallSeqLenInteraction, type = "trace")
# Plot effects 
mcmc_plot(mrtallRecallSeqLenInteraction, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(mrtallRecallSeqLenInteraction, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  # scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", 
  #                             "Sequence Structure", "Sort Task", "Recall", 
  #                             "Seq Len*Recall", "QueryStruc*Recall", "SeqStruc*Recall",
  #                             "SortTask*Recall"))+# 
  ggtitle("Ecoding vs. Recall RT analysis")
pp_check(mrtallRecallSeqLenInteraction, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mrtallRecallSeqLenInteraction), points = FALSE)
tab_model(mrtallRecallSeqLenInteraction)

# Results 26.09.2022
#                                                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                                               0.89      0.37     0.66     1.17 1.50        8       20
# Number_of_Bars                                          0.53      0.09     0.42     0.59 1.45        9       11
# Stimulus_typequery                                      0.16      0.18    -0.05     0.74 1.43        8       11
# StructureQuery                                          0.32      0.14     0.16     0.78 1.36        9       11
# StructureSequence                                      -0.20      0.33    -1.59    -0.01 1.55        7       11
# ConditionSort                                          -0.05      0.08    -0.26     0.03 1.19       14       13
# Number_of_Bars:Stimulus_typequery                      -0.34      0.17    -0.48    -0.25 1.38        9       11
# Number_of_Bars:StructureQuery                          -0.16      0.15    -0.27    -0.14 1.57       10       13
# Number_of_Bars:StructureSequence                       -0.02      0.19    -0.03     0.19 1.24       15       14
# Number_of_Bars:ConditionSort                            0.05      0.20     0.04     0.10 1.35       15       25
# Stimulus_typequery:StructureQuery                      -0.32      0.37    -1.77    -0.12 1.45        8       11
# Stimulus_typequery:StructureSequence                    0.04      0.08    -0.14     0.17 1.28       19       26
# Stimulus_typequery:ConditionSort                        0.02      0.05    -0.07     0.10 1.37        9       26
# Number_of_Bars:Stimulus_typequery:StructureQuery        0.14      0.14     0.08     0.45 1.32       10       11
# Number_of_Bars:Stimulus_typequery:StructureSequence     0.02      0.04    -0.01     0.06 1.18       22       27
# Number_of_Bars:Stimulus_typequery:ConditionSort        -0.04      0.04    -0.06    -0.01 1.21       14       29

# Warning messages:
#   1: There were 268 divergent transitions after warmup. See
# https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 2: There were 4731 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
# https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
# https://mc-stan.org/misc/warnings.html#bfmi-low 
# 4: Examine the pairs() plot to diagnose sampling problems
# 
# 5: The largest R-hat is NA, indicating chains have not mixed.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#r-hat 
# 6: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#bulk-ess 
# 7: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#tail-ess 

###################################################################
# Full encoding RT model with recall RT as a predictor
###################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# The full Encoding RT model 
mrtallStructConditionPlusRecallRT <- brm(rt ~ Number_of_Bars + Structure + Condition + queryRT+ (Number_of_Bars + Structure + Condition + Block + queryRT|subject_id), 
                             data = dp,
                             iter = 10000,
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             family = exgaussian(),
                             file = "mrtallStructConditionPlusRecallRTNotLogFinalEx")
# Model summary
summary(mrtallStructConditionPlusRecallRT)
# Check for convergence
mcmc_plot(mrtallStructConditionPlusRecallRT, type = "trace")
# Plot effects 
mcmc_plot(mrtallStructConditionPlusRecallRT, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(mrtallStructConditionPlusRecallRT, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Raw RT analysis")
pp_check(mrtallStructConditionPlusRecallRT, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mrtallStructConditionPlusRecallRT), points = FALSE)
tab_model(mrtallStructConditionPlusRecallRT, transform = NULL)

# Results 23.03.2023
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.55      0.11     0.34     0.77 1.00     2168     4154
# Number_of_Bars        0.58      0.04     0.50     0.66 1.00     1383     3134
# StructureQuery       -0.19      0.05    -0.29    -0.08 1.00     2985     6580
# StructureSequence    -0.17      0.04    -0.24    -0.09 1.00     2831     8094
# ConditionSort         0.22      0.04     0.14     0.30 1.00     4125     8110
# queryRT               0.29      0.04     0.21     0.36 1.00     5462     9239



###################################################################
# Full encoding RT with exponential term model reply to reviewer 1
###################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)
dp$NoBsquared <- dp$Number_of_Bars**2

# The full Encoding RT model 
mrtallStructConditionWithQuadraticTerm <- brm(rt ~ Number_of_Bars + NoBsquared + Structure + Condition + (Number_of_Bars + NoBsquared + Structure + Condition + Block|subject_id), 
                             data = dp,
                             iter = 10000,
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             family = exgaussian(),
                             file = "mrtallStructConditionNotLogFinalExWithQuadraticTerm")
# Model summary
summary(mrtallStructConditionWithQuadraticTerm)
# Check for convergence
mcmc_plot(mrtallStructConditionWithQuadraticTerm, type = "trace")
# Plot effects 
mcmc_plot(mrtallStructConditionWithQuadraticTerm, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(mrtallStructConditionWithQuadraticTerm, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Sequence Length^2","Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Raw RT analysis")
pp_check(mrtallStructConditionWithQuadraticTerm, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mrtallStructConditionWithQuadraticTerm), points = FALSE)
tab_model(mrtallStructConditionWithQuadraticTerm, transform = NULL)


# results exgaussian  29.03.2023
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.81      0.08     0.65     0.98 1.00     2418     4762
# Number_of_Bars        0.47      0.09     0.29     0.65 1.00     1331     2754
# NoBsquared            0.03      0.01     0.00     0.06 1.00     1309     2526
# StructureQuery       -0.26      0.08    -0.41    -0.11 1.00     2782     5510
# StructureSequence    -0.18      0.04    -0.26    -0.09 1.00     3826     8019
# ConditionSort         0.27      0.05     0.17     0.36 1.00     4524     8188

# BF with not quadratic model
bayes_factor(mrtallStructConditionWithQuadraticTerm, mrtallStructCondition)
