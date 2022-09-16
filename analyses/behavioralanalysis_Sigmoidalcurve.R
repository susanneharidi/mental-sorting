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
library(stats)

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project

#########################################################################
# functions
#########################################################################

#standard error function
se <- function(x){sd(x)/sqrt(length(x))}

# function needed for visualization purposes
sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3])))
}

###############################################################
# variables 
###############################################################

rt_cutoff <- 10
seed <- 2022
myColors <- c("#a2f9ff","#00AFBB", "#ffe682", "#E7B800", "#fd8453", "#FC4E07")
mylabels <- c("None Memory", "None Sort", "Query Memory", "Query Sort", "Sequence Memory", "Sequence Sort")
alpha <- 0.7
pd <- position_dodge(0.1)
textsize <- 16

#color pallette for three different groups
cbbPalette <- c("grey50", "#0072B2", "#D55E00")


#########################################################################
# load data 
#########################################################################

d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))

######################################################################
# fit a sigmoid for each participant
####################################################################

dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

ids <- unique(dp$subject_id)
asyms <- c()
xmids <- c()
scales <- c()

Params <- data.frame(Asym = numeric(),
                    xmid = numeric(),
                    scales = numeric(),
                    id = character())
Preds <- data.frame(id = character(),
                          Number_of_Bars = numeric(),
                          predictions = numeric(),
                          rt = numeric())

outs = c("rbrzxgywqun5ned", "k45h5lxk13xjkxs", "mqng9gal5nqwgz9")

for (id in ids)
{ if (id %in% outs){
    print("there where errors with the following id")
    print(id)
  } else {
    subject = subset(dp, subject_id == id)
    x = subject$Number_of_Bars
    y = subject$mu
    print(id)
    
    # fitting code
    fitmodel <- nls(rt ~ SSlogis(Number_of_Bars, Asym, xmid, scal), data = subject)
    
    # visualization code
    # get the coefficients using the coef function
    
    params = coef(fitmodel)
    currentParams = data.frame(Asym = params[1],
                              xmid = params[2],
                              scales = params[3],
                              id = id
    )
    Params = rbind(Params, currentParams)
    currentPreds = data.frame(id = id,
                              Number_of_Bars = subject$Number_of_Bars,
                              predictions = predict(fitmodel),
                              rt = subject$rt)
    
    Preds = rbind(Preds, currentPreds)
  }
}

# calculate mean Fit
meanRT = ddply(dp, ~Number_of_Bars+subject_id, summarize, mean = mean(rt))
fitmodel <- nls(mean ~ SSlogis(Number_of_Bars, Asym, xmid, scal), data = meanRT)
meanFit = data.frame(Number_of_Bars = meanRT$Number_of_Bars,
                     predictions = predict(fitmodel))
meanFit = ddply(meanFit, ~Number_of_Bars, summarize, mean = mean(predictions))

summarizedRT = ddply(dp, ~Number_of_Bars, summarize, mean = mean(rt))
ggplot(Preds, aes(x = Number_of_Bars, y = predictions, col = id))+
  geom_line(alpha = 0.3)+
  geom_line(aes(x = Number_of_Bars, y = mean), data = summarizedRT, color = "Black", size = 1)+
  geom_line(aes(x = Number_of_Bars, y = mean), data = meanFit, color = "red", size = 1)+
  theme(legend.position = "none")


#Arguments
#input - a numeric vector of values at which to evaluate the model.
#Asym - a numeric parameter representing the asymptote. 
#xmid - a numeric parameter representing the x value at the inflection point of the curve. The value of SSlogis will be Asym/2 at xmid.
#scal - a numeric scale parameter on the input axis.

##################################################
# make a new Number of Bars Variable
#################################################
dp$SigmoidNoB = dp$Number_of_Bars
for (i in 1:7){
  dp["SigmoidNoB"][dp["SigmoidNoB"] == i] <- meanFit$mean[i]
}

# plot the mean of that
SigmoidRT = ddply(dp, ~SigmoidNoB, summarize, mean = mean(rt))
SigmoidRTSub = ddply(dp, ~SigmoidNoB + subject_id, summarize, mean = mean(rt))
ggplot(SigmoidRTSub, aes(x = SigmoidNoB, y = mean, col = subject_id))+
  geom_point()+
  geom_line(alpha = 0.3)+
  geom_line(aes(x = SigmoidNoB, y = mean), data = SigmoidRT, color = "Black", size = 1)+
  theme(legend.position = "none")

##############################################################################
# fit the model with the sigmoid Sequence length
##############################################################################

setwd(here("analyses", "brmModel"))
# The full Encoding RT model 
mrtallStructConditionInterSigmoid <- brm(rt ~ SigmoidNoB*(Structure + Condition) + (SigmoidNoB*(Structure + Condition) + Block|subject_id), 
                                  data = dp,
                                  iter = 10000,
                                  cores = 4,
                                  save_pars = save_pars(all = TRUE),
                                  seed = seed,
                                  file = "mrtallStructConditionNotLogInteractionSigmoidFinal")
  
# Model summary
summary(mrtallStructConditionInterSigmoid)
# Check for convergence
mcmc_plot(mrtallStructConditionInterSigmoid, type = "trace")
# Plot effects 
mcmc_plot(mrtallStructConditionInterSigmoid, type = "intervals", prob = 0.95)
# only plot main populationlevel effects
mcmc_plot(mrtallStructConditionInterSigmoid, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("RT in s")+
  theme_minimal()+
  theme(text = element_text(size = 20, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task", "Seq length*Query Struc", "Seq length*Seq Struc", "Seq length*Sort Task"))+# 
  ggtitle("Raw RT analysis")
pp_check(mrtallStructConditionInterSigmoid, type = "ecdf_overlay")
#visualize all effects
plot(conditional_effects(mrtallStructConditionInterSigmoid), points = FALSE)
tab_model(mrtallStructConditionInterSigmoid)

# results 15.09
#                                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                       -0.13      0.15    -0.43     0.17 1.00     2197     3570
# SigmoidNoB                       1.03      0.07     0.90     1.16 1.00     1081     2933
# StructureQuery                   0.21      0.07     0.07     0.36 1.00     5778    10876
# StructureSequence               -0.11      0.06    -0.23     0.02 1.00    11371    14288
# ConditionSort                   -0.11      0.07    -0.24     0.02 1.00     5685    11592
# SigmoidNoB:StructureQuery       -0.16      0.04    -0.24    -0.08 1.00     4257     8994
# SigmoidNoB:StructureSequence    -0.04      0.02    -0.09     0.01 1.00     8995    13574
# SigmoidNoB:ConditionSort         0.18      0.03     0.12     0.24 1.00     7405    12420

##############
# Plot the Residuals of Model with Intercations
###############
residulas <- residuals(mrtallStructConditionInterSigmoid)
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
