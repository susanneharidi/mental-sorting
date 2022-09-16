
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

#########################################################################
# Looking at the residulas of the no structure model to see, 
# weather structure is learned over trials
##############################################################################
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

##############
# Plot the Residuals according to trial number
###############
dp$TrialNr <- dp$Trial_index-(dp$Block-1)*35
residulas <- residuals(mrtnoStructure)
dp$Residuals <- residulas[,1]

#get summarize data
dRTSummary <- ddply(dp, ~TrialNr + Condition + Structure, summarize, mu = mean(Residuals), se = se(Residuals))

#start plotting the residuals
pResiduals <- ggplot(dRTSummary, aes(x = TrialNr, y = mu, col = Structure)) +
  #points
  geom_point(position = pd)+
  #error bars
  #geom_errorbar(aes(ymin = mu-se, ymax = mu+se), width = 0, size = 1, position = pd) +
  #lines
  geom_line(position = pd, size = 1.2) +
  scale_color_manual(values = cbbPalette)+
  geom_smooth(method = lm)+
  #legend on top
  theme(strip.background = element_blank())+
  #change x ticks
  scale_x_continuous(breaks = round(seq(min(0), max(35), by = 5),1)) +
  #minimal theme
  theme_minimal()+
  #add xlab
  xlab("Trial")+
  #add ylab
  ylab('Mean Residualsin s (±SE)')+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"))+
  #legend on top
  theme(legend.position = "top")+
  facet_grid(cols = vars(Condition))

#show!
pResiduals
