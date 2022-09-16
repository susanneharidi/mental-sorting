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
                             file = "mrtallStructConditionNotLogFinal")
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
tab_model(mrtallStructCondition)


# results form 26.7.2022
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            -0.26      0.17    -0.60     0.09 1.00     2504     4533
# Number_of_Bars        0.89      0.05     0.79     1.00 1.00     1560     3300
# StructureQuery       -0.29      0.09    -0.47    -0.11 1.00     3373     5909
# StructureSequence    -0.24      0.06    -0.35    -0.12 1.00     3857     8564
# ConditionSort         0.42      0.07     0.29     0.56 1.00     4945     9155

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
mrtallStructConditionInter <- brm(rt ~ Number_of_Bars*(Structure + Condition) + (Number_of_Bars*(Structure + Condition) + Block|subject_id), 
                             data = dp,
                             iter = 10000,
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             file = "mrtallStructConditionNotLogInteractionFinal")
# Model summary
summary(mrtallStructConditionInter)
# Check for convergence
mcmc_plot(mrtallStructConditionInter, type = "trace")
# Plot effects 
mcmc_plot(mrtallStructConditionInter, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(mrtallStructConditionInter, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task", "Seq length*Query Struc", "Seq length*Seq Struc", "Seq length*Sort Task"))+# 
  ggtitle("Raw RT analysis")
pp_check(mrtallStructConditionInter, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mrtallStructConditionInter), points = FALSE)
tab_model(mrtallStructConditionInter)


# results form 14.09.2022
#                                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                           -0.17      0.16    -0.48     0.13 1.00     2239     3879
# Number_of_Bars                       0.88      0.06     0.76     0.99 1.00     1225     2188
# StructureQuery                       0.23      0.08     0.08     0.38 1.00     4904    10394
# StructureSequence                   -0.10      0.06    -0.23     0.03 1.00    12236    14677
# ConditionSort                       -0.13      0.07    -0.27     0.00 1.00     5923    12727
# Number_of_Bars:StructureQuery       -0.14      0.04    -0.21    -0.07 1.00     3740     8916
# Number_of_Bars:StructureSequence    -0.03      0.02    -0.08     0.01 1.00     9875    12994
# Number_of_Bars:ConditionSort         0.16      0.03     0.11     0.21 1.00     7294    12127

##############
# Plot the Residuals of Model with Intercations
###############
residulas <- residuals(mrtallStructConditionInter)
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
                                  file = "mrtallStructConditionNotLogFinalBlock")
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


# results form 26.07.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.34      0.15     0.05     0.63 1.00     2233     3963
# Number_of_Bars        0.89      0.04     0.81     0.96 1.00     1371     2828
# StructureQuery       -0.37      0.06    -0.49    -0.24 1.00     4370     8299
# StructureSequence    -0.20      0.04    -0.28    -0.11 1.00     9517    13123
# ConditionSort         0.38      0.05     0.28     0.48 1.00     5576     9152
# Block                -0.13      0.02    -0.16    -0.10 1.00     4349     9398

###################################################################
# Full recall rt model 
###################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# The full Recall RT model
mRecallRTallStructCondition <- brm(queryRT ~ Number_of_Bars + Structure + Condition + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                                   data = dp,
                                   iter = 10000,
                                   cores = 4,
                                   save_pars = save_pars(all = TRUE),
                                   seed = seed,
                                   file = "mmRecallRTallStructConditionNotLogFinal")
# Model summary
summary(mRecallRTallStructCondition)
# Check for convergence
mcmc_plot(mRecallRTallStructCondition, type = "trace")
# Plot effects 
mcmc_plot(mRecallRTallStructCondition, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(mRecallRTallStructCondition, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_classic()+
  theme(text = element_text(size = 15, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Number of Bars", "Query Structure", "Sequence Structure", "Sort Task"))+
  ggtitle("Recall RT")
pp_check(mRecallRTallStructCondition, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mRecallRTallStructCondition), points = FALSE)
tab_model(mRecallRTallStructCondition)

# results form 26.07.2022
# Intercept             0.04      0.06    -0.08     0.17 1.00     2328     5941
# Number_of_Bars        0.35      0.03     0.30     0.41 1.00     1411     2931
# StructureQuery       -0.15      0.06    -0.26    -0.04 1.00     3877     7702
# StructureSequence    -0.01      0.02    -0.06     0.03 1.00    10692    13260
# ConditionSort         0.18      0.05     0.07     0.28 1.00     2616     5700


###############################################################################
###################################################################
# include recall rt as a factor in encoding RT model
###################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# regression
mrtallincludingRecallRT <- brm(rt ~ Number_of_Bars + Structure + Condition + queryRT + (Number_of_Bars + Structure + Condition + queryRT + Block|subject_id), 
                               data = dp,
                               iter = 10000,
                               cores = 4,
                               save_pars = save_pars(all = TRUE),
                               seed = seed,
                               file = "mrtallincludingRecallRTNotLogFinal")
summary(mrtallincludingRecallRT)
# results 27.07.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            -0.25      0.17    -0.58     0.07 1.00     2523     5367
# Number_of_Bars        0.81      0.05     0.71     0.90 1.00     2217     4489
# StructureQuery       -0.25      0.07    -0.39    -0.12 1.00     4047     8765
# StructureSequence    -0.23      0.05    -0.33    -0.12 1.00     5152     9022
# ConditionSort         0.36      0.06     0.24     0.48 1.00     5906    10569
# queryRT               0.29      0.05     0.19     0.38 1.00     7826    12171

# only plot main population level effects
mcmc_plot(mrtallincludingRecallRT, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme(text = element_text(size = 15, family = "sans"))+
  theme_classic()+
  scale_y_discrete(labels = c("Intercept",
                              "Sequence Length",
                              "Query Structure",
                              "Sequence Structure",
                              "Sort Task",
                              "Recall RT"))+
  ggtitle("Raw RT analysis including Recall RT as predictor")

# Model summary
summary(mrtallincludingRecallRT)
# Check for convergence
mcmc_plot(mrtallincludingRecallRT, type = "trace")
# Plot effects 
mcmc_plot(mrtallincludingRecallRT, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
pp_check(mrtallincludingRecallRT, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mrtallincludingRecallRT), points = FALSE)
tab_model(mrtallincludingRecallRT)


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
                file = "mrtnoNoBNotLogFinal")
summary(mrtnoNoB)
# results 07.08.2022
# Intercept             1.06      0.30     0.46     1.65 1.01      337      719
# StructureQuery       -0.29      0.15    -0.58     0.01 1.00     4602     8909
# StructureSequence    -0.17      0.09    -0.34     0.03 1.00     7777     9366
# ConditionSort         0.31      0.12     0.06     0.54 1.00     4220     8578

bayes_factor(mrtallStructCondition, mrtnoNoB)
# results 07.08.2022
#Warning messages:
#  1: logml could not be estimated within maxiter, rerunning with adjusted starting value. 
#Estimate might be more variable than usual. 
#2: logml could not be estimated within maxiter, rerunning with adjusted starting value. 
#Estimate might be more variable than usual. 
# BF =  1726925481417150973224220426.00000

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
                      file = "mrtnoStructureNotLogFinal")

summary(mrtnoStructure)
# results 07.08.2022
# Intercept         -0.40      0.18    -0.75    -0.05 1.00     1945     3994
# Number_of_Bars     0.88      0.06     0.77     1.00 1.00     1150     2529
# ConditionSort      0.45      0.07     0.31     0.59 1.00     4712     9510

bayes_factor(mrtallStructCondition, mrtnoStructure)
# results 07.08.2022
#Warning message:
#  logml could not be estimated within maxiter, rerunning with adjusted starting value. 
#Estimate might be more variable than usual. 
# BF = 9.28885

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
                      file = "mrtnoConditionNotLogFinal")

summary(mrtnoCondition)
# results 07.08.2022
# Intercept             0.11      0.19    -0.27     0.48 1.00     1313     3152
# Number_of_Bars        0.85      0.07     0.72     0.98 1.00      786     1403
# StructureQuery       -0.37      0.10    -0.58    -0.17 1.00     1529     2877
# StructureSequence    -0.23      0.06    -0.36    -0.10 1.01     1540     4601


bayes_factor(mrtallStructCondition, mrtnoCondition)
# results 29.06.2022
# BF = 275273816.09508
#Warning messages:
#  1: logml could not be estimated within maxiter, rerunning with adjusted starting value. 
#Estimate might be more variable than usual. 
#2: logml could not be estimated within maxiter, rerunning with adjusted starting value. 
#Estimate might be more variable than usual. 


############################################################################
# analyse the sum of recall Rt and Encoding Rt together
############################################################################
d$fullRT <- d$rt + d$queryRT

###################################################################
# model rt all structures + memory condition not log and scaled
###################################################################

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# regression
mFullRTallStructCondition <- brm(fullRT ~ Number_of_Bars + Structure + Condition + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                                 data = dp,
                                 iter = 10000,
                                 cores = 4,
                                 save_pars = save_pars(all = TRUE),
                                 seed = seed,
                                 file = "mFullRTallStructConditionNotLogFinal")

# Model summary
summary(mFullRTallStructCondition)
# results 07.08.2022
# Intercept            -0.36      0.19    -0.72     0.03 1.00     3326     5759
# Number_of_Bars        1.27      0.06     1.15     1.40 1.00     1355     3063
# StructureQuery       -0.43      0.15    -0.72    -0.14 1.00     3216     5698
# StructureSequence    -0.23      0.07    -0.37    -0.10 1.00     6255    10813
# ConditionSort         0.62      0.11     0.41     0.83 1.00     5386    10352


# Check for convergence
mcmc_plot(mFullRTallStructCondition, type = "trace")
# Plot effects 
mcmc_plot(mFullRTallStructCondition, type = "intervals", prob = 0.95)
# only plot main population level effects
mcmc_plot(mFullRTallStructCondition, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_classic()+
  theme(text = element_text(size = 15, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Number of Bars", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Raw Total (recall+encoding) RT analysis ")
pp_check(mFullRTallStructCondition, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mFullRTallStructCondition), points = FALSE)
tab_model(mFullRTallStructCondition)


