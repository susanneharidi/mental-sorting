# This script contains the accuracy time trade offs described in the corresponding section in the supplements

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
# load data and exclude not independent data and people who do not meet the performance threshold
#########################################################################


d <- read.csv(here("Data","Experiment","data_to_work_with.csv"))


# if you are not using rproject or here packages, than use the lines below and 
# d <- read.csv("Data/Experiment/data_to_work_with.csv")

#########################################################################
# parameters an functions
########################################################################
# rt cutoff in s
cutoff = 10

#dodged position
pd <- position_dodge(.2)
#color pallette for three different groups
cbbPalette <- c("grey50","#0072B2", "#D55E00")

#standard error
se<-function(x){sd(x)/sqrt(length(x))}


#########################################################################
#Calculate participant specific effect of errors and RT
#########################################################################
#Only use the data that meet the rt cutoff criteria, of course here we also need the incorrect trials
valid_dataSort <- subset(d, rt <= cutoff & Condition == "Sort" & queryRT <= cutoff)

setwd(here("analyses","brmModel"))

# Acurracy Model
mAccuracyIndividual <- brm(correct ~ Number_of_Bars + (Number_of_Bars|subject_id),
                           data = valid_dataSort,
                           save_pars = save_pars(all = TRUE),
                           iter = 10000,
                           cores = 4,
                           chains = 4,
                           control = list(max_treedepth = 15, adapt_delta = 0.99), # the tree depth is an efficiency thing, while the delta is supposed to help with the validity issue due to divergent transitions after warmup
                           family = "bernoulli",
                           file = "mAccuracyIndividualSort")
print(summary(mAccuracyIndividual))
# 27.07.2022
#                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept          5.51      0.30     4.97     6.12 1.00     6560    10972
# Number_of_Bars    -0.49      0.05    -0.60    -0.39 1.00     7349    11136
RandomEffectsAccuracyIndividual = data.frame(ranef(mAccuracyIndividual))
RandomEffectsAccuracyIndividual$subject_id = sort(unique(valid_dataSort$subject_id))
RandomEffectsAccuracyIndividual$AccuracyEstimate = RandomEffectsAccuracyIndividual$subject_id.Estimate.Number_of_Bars
RandomEffectsAccuracyIndividual = RandomEffectsAccuracyIndividual[, c("subject_id", "AccuracyEstimate")]

# individual rt estimates
mrtScalingIndividual <- brm(rt ~ Number_of_Bars + (Number_of_Bars|subject_id),
                            data = valid_dataSort,
                            iter = 10000,
                            chains = 4,
                            cores = 4,
                            save_pars = save_pars(all = TRUE),
                            control = list(max_treedepth = 15, adapt_delta = 0.99),
                            file = "mrtScalingIndividual_with_incorrectSort")
print(summary(mrtScalingIndividual))
# 27.07.2022
#                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept         -0.13      0.13    -0.39     0.13 1.01     1722     3432
# Number_of_Bars     0.94      0.04     0.86     1.03 1.01      901     1922
RandomEffectsRtIndividual = data.frame(ranef(mrtScalingIndividual))
RandomEffectsRtIndividual$subject_id = sort(unique(valid_dataSort$subject_id))
RandomEffectsRtIndividual$RTEstimate = RandomEffectsRtIndividual$subject_id.Estimate.Number_of_Bars
RandomEffectsRtIndividual = RandomEffectsRtIndividual[, c("subject_id", "RTEstimate")]

AccuracyRTRandomEffects = merge(RandomEffectsAccuracyIndividual, RandomEffectsRtIndividual, by = "subject_id")

################################################################################
# calculate the correlation
################################################################################
cor(AccuracyRTRandomEffects$RTEstimate, AccuracyRTRandomEffects$AccuracyEstimate, 
    method = "pearson")
# results 27.07.2022
# persons R  = -0.04511333
cor.test(AccuracyRTRandomEffects$RTEstimate, AccuracyRTRandomEffects$AccuracyEstimate,
         method = "pearson")
# results 27.07.2022
# data:  AccuracyRTRandomEffects$RTEstimate and AccuracyRTRandomEffects$AccuracyEstimate
# t = -0.38052, df = 71, p-value = 0.7047
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2723538  0.1868937
# sample estimates:
#   cor 
# -0.04511333  
#############################
# plot just for some visuals
#############################

plotAccuracyVStime <- ggplot(AccuracyRTRandomEffects, aes(x = RTEstimate, y = AccuracyEstimate)) +
  #points
  geom_point()+
  geom_smooth(method="lm")+ 
  #minimal theme
  theme_minimal()+
  #add xlab
  xlab("Effect of Sequence Length on RT")+
  #add ylab
  ylab('Effect of Sequence Length on Accuracy')+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"))+
  #add title
  theme(legend.position="none")

plotAccuracyVStime

ggsave(here("Figures", "AccurcayTimeTradeOffFigure.pdf"), plotAccuracyVStime,  
       width = 10, 
       height = 10)

###########################################################################
# model logistic regression for the accuracy time trade of
#########################################################################
data <- subset(d, rt <= cutoff & Condition == "Sort" & queryRT <= cutoff)
#accuracy ~ speed + ... + (1|participant ID)
accuracyLogistic <- glm(correct ~ Number_of_Bars + rt + Structure + (Number_of_Bars + rt + Structure|subject_id),
                        family = binomial(link='logit'),
                        data = data)
