# Looking at the effect of Structure
# final analysis for publication
# corresponds to the results in the "Structure helps, but is not used to its full extent" section of the paper


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
library(here)

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project


#########################################################################
# functions
#########################################################################

#standard error function
se<-function(x){sd(x)/sqrt(length(x))}


#########################################################################
# load data and exclude not independent data and people who do not meet the performance threshold
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

#color pallette for three different groups
cbbPalette <- c("grey50", "#0072B2", "#D55E00")

#other parameters in s
rt_cutoff <- 10
seed <- 2022


##############################################################################################################################################################
# Modelling
#########################################################################

# important!!!!! This is where all the models are saved, so we don't have to recalculate them every time
# so RUN THIS!
setwd(here("analyses", "brmModel"))
# setwd("./analyses/brmModel")

##########################################
#prepare data
#############################################
# Only use the correct and the data that meet the rt cutoff criteria
dStrucDiff <- subset(d, rt <= rt_cutoff & correct == 1 & Condition == "Sort" & queryRT <= rt_cutoff)

# calculate the differences of the structure conditions
dNone_Query <- ddply(dStrucDiff, ~ TrialID + subject_id + Number_of_Bars + Stimulus_type, 
                     summarize, m = rt[Conditions == "noStructureExp"] - rt[Conditions == "queryStructureExp"])
dNone_Sequence <- ddply(dStrucDiff, ~ TrialID + subject_id + Number_of_Bars + Stimulus_type, 
                        summarize, m = rt[Conditions == "noStructureExp"] - rt[Conditions == "sortStructureExp"])

##############################################################
# difference between no structure and query structure
################################################################
# set up dataframe
dp <- dNone_Query

# regression for the constant model
mrtStrucDiffQueryNull <- brm(m ~ 1 + (Number_of_Bars |subject_id), 
                         data = dp,
                         iter = 10000,
                         cores = 4,
                         save_pars = save_pars(all = TRUE),
                         control = list(max_treedepth = 15, adapt_delta = 0.99),
                         seed = seed,
                         file = "mrtStrucDiffQueryNullFinal")

# Model summary
summary(mrtStrucDiffQueryNull)
#results 08.08.2022
# Intercept     0.02      0.07    -0.12     0.15 1.00    14663    14004

############################################################################
# regression for the linear model
mrtStrucDiffQuery <- brm(m ~ Number_of_Bars + (Number_of_Bars |subject_id), 
                             data = dp,
                             iter = 10000,
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             control = list(max_treedepth = 15, adapt_delta = 0.99),
                             seed = seed,
                             file = "mrtStrucDiffQueryFinal")

# Model summary
summary(mrtStrucDiffQuery)
#results 08.08.2022
# Intercept         -0.16      0.08    -0.32    -0.01 1.00    22336    17015
# Number_of_Bars     0.14      0.04     0.06     0.21 1.00    15928    14696


###############################################################################
# 1233333 hypothesis
dp$Number_of_Bars <- ifelse(dp$Number_of_Bars > 3, 3, dp$Number_of_Bars)

# regression for the 1233333 model
mrtStrucDiffQuery2 <- brm(m ~ Number_of_Bars + (Number_of_Bars |subject_id), 
                         data = dp,
                         iter = 10000,
                         cores = 4,
                         save_pars = save_pars(all = TRUE),
                         control = list(max_treedepth = 15, adapt_delta = 0.99),
                         seed = seed,
                         file = "mrtStrucDiffQuery2Final")

# Model summary
summary(mrtStrucDiffQuery2)
#results 08.08.2022
# Intercept         -0.27      0.12    -0.51    -0.02 1.00    17266    14779
# Number_of_Bars     0.23      0.07     0.09     0.38 1.00    10015    13249

##################################
# BF for query structure
###################################

# Linear vs Constant
bayes_factor(mrtStrucDiffQuery, mrtStrucDiffQueryNull)
# results 08.08.2022
# BF = 50.84994
# Linear Wins

# Linear vs 1233333
bayes_factor(mrtStrucDiffQuery, mrtStrucDiffQuery2)
# results 08.08.2022
# BF = 1414209381796690.25000
# Linear Wins

# 1233333 vs constant
bayes_factor(mrtStrucDiffQuery2, mrtStrucDiffQueryNull)
# results 08.08.2022
# BF = 0.00000
# constant wins

# constant vs 1233333
bayes_factor(mrtStrucDiffQueryNull, mrtStrucDiffQuery2)
# results 08.08.2022
# BF = 35397537893822.34375
# constant wins

# 1233333 vs Linear
bayes_factor(mrtStrucDiffQuery2, mrtStrucDiffQuery)
# results 08.08.2022
# BF = 0.00000
# Linear wins




####################################################################
# difference between no structure and sequence structure
###################################################################

#create data set as before
dp <- dNone_Sequence


# regression for the constant model
mrtStrucDiffSequenceNull <- brm(m ~ 1 + (Number_of_Bars |subject_id), 
                            data = dp,
                            iter = 10000,
                            cores = 4,
                            save_pars = save_pars(all = TRUE),
                            control = list(max_treedepth = 15, adapt_delta = 0.99),
                            seed = seed,
                            file = "mrtStrucDiffSequenceNullFinal")

# Model summary
summary(mrtStrucDiffSequenceNull)
#results 08.08.2022
#Warning message:
#Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#Running the chains for more iterations may help. See
#https://mc-stan.org/misc/warnings.html#bulk-ess 
# Intercept     0.11      0.06    -0.01     0.22 1.00    18974    15156

################################################################
# regression for the Linear model
mrtStrucDiffSequence <- brm(m ~ Number_of_Bars + (Number_of_Bars |subject_id), 
                         data = dp,
                         iter = 10000,
                         cores = 4,
                         save_pars = save_pars(all = TRUE),
                         control = list(max_treedepth = 15, adapt_delta = 0.99),
                         seed = seed,
                         file = "mrtStrucDiffSequenceFinal")

# Model summary
summary(mrtStrucDiffSequence)
#results 21.02.2022
#Warning message:
#Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#Running the chains for more iterations may help. See
#https://mc-stan.org/misc/warnings.html#bulk-ess 
# Intercept          0.03      0.07    -0.11     0.16 1.00    26836    15793
# Number_of_Bars     0.05      0.03    -0.00     0.11 1.00    12325    13524


#########################################################################
# 1233333 hypothesis
dp$Number_of_Bars <- ifelse(dp$Number_of_Bars > 3, 3, dp$Number_of_Bars)


# regression for the 1233333 model
mrtStrucDiffSequence2 <- brm(m ~ Number_of_Bars + (Number_of_Bars |subject_id), 
                            data = dp,
                            iter = 10000,
                            cores = 4,
                            save_pars = save_pars(all = TRUE),
                            control = list(max_treedepth = 15, adapt_delta = 0.99),
                            seed = seed,
                            file = "mrtStrucDiffSequence2Final")

# Model summary
summary(mrtStrucDiffSequence2)
#results 08.08.2022
# Intercept         -0.14      0.11    -0.36     0.07 1.00    25182    15442
# Number_of_Bars     0.14      0.06     0.03     0.25 1.00    17486    16245



##################################
# BF for sequence structure
##################################

# Linear vs Constant
bayes_factor(mrtStrucDiffSequence, mrtStrucDiffSequenceNull)
#results 08.08.2022
# BF = 0.18812

bayes_factor(mrtStrucDiffSequenceNull, mrtStrucDiffSequence)
# BF = 2.32601

# Constant wins

# 1233333 vs Constant
bayes_factor(mrtStrucDiffSequence2, mrtStrucDiffSequenceNull)
#results 08.08.2022
#BF = 0.80448

bayes_factor(mrtStrucDiffSequenceNull, mrtStrucDiffSequence2)
#BF = 1.25963

# results unclear

# 1233333 vs Linear
bayes_factor(mrtStrucDiffSequence2, mrtStrucDiffSequence)
# results 08.08.2022
# BF = 1.97287

bayes_factor(mrtStrucDiffSequence, mrtStrucDiffSequence2)
# BF = 0.91989

# 1233333 Wins