# Scaling analysis
# what is analysed here is the scaling of the  Structure condition difference between
# the RT in the sort and the memory task for each prematched trial (matched for length and Structure)
# the results here correspond to the results reported in the "Scaling anlysis" section
# in the paper

##############################################################################
# load libraries
################################################################################

library(ggsignif)
library(gridExtra)
library(ggthemes)
library(ggplot2)
library(plyr)
library(ggpubr)
library(here)
library(brms)
library(ggpubr)
library(MASS)
library(plot.matrix)
library(RColorBrewer)
library(corrplot)

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project

#################################################################################
# load data
###########################################################################

d <- read.csv(here("Data","Experiment","data_to_work_with.csv"))

#if you are not using rproject or here packages, than use the lines below and 
#data <- read.csv("Data/Experiment/data_to_work_with.csv")

########################################################################
# parameters an functions
########################################################################
# set the maximum rt we still include for each trial in s:
cutoff <- 10
# set the seed for the models
seed <- 2022

#dodged position
pd <- position_dodge(.2)
#color pallette for three different groups
cbbPalette <- c("grey50","#0072B2", "#D55E00")

########################################################################
# functions
########################################################################

#standard error
se<-function(x){sd(x)/sqrt(length(x))}

# this is a very specific function to do a round robin comparison
# of the BFs of models where the predictor "sequence length" has been 
# transformed to describe different scaling factors
# I used that here, because I needed to make the same comparison for all
# three structure conditions
BFCalulcation <- function(ConstantModel, LogModel, LinearModel, NLogNModel, Poli2Model, ExpoModel){
  # Constant
  BFConstLog = bayes_factor(ConstantModel, LogModel)
  BFConstLinear = bayes_factor(ConstantModel, LinearModel)
  BFConstNLogN = bayes_factor(ConstantModel, NLogNModel)
  BFConstPoli2 = bayes_factor(ConstantModel, Poli2Model)
  BFConstExpo = bayes_factor(ConstantModel, ExpoModel)
  BFConst = c(NaN, BFConstLog$bf, BFConstLinear$bf, 
              BFConstNLogN$bf,  BFConstPoli2$bf, BFConstExpo$bf)
  # Log
  BFLogConst = bayes_factor(LogModel, ConstantModel)
  BFLogLinear = bayes_factor(LogModel, LinearModel)
  BFLogNLogN = bayes_factor(LogModel, NLogNModel)
  BFLogPoli2 = bayes_factor(LogModel, Poli2Model)
  BFLogExpo = bayes_factor(LogModel, ExpoModel)
  BFLog = c(BFLogConst$bf, NaN, BFLogLinear$bf, 
            BFLogNLogN$bf,  BFLogPoli2$bf, BFLogExpo$bf)
  # Linear
  BFLinearConst = bayes_factor(LinearModel, ConstantModel)
  BFLinearLog = bayes_factor(LinearModel, LogModel)
  BFLinearNLogN = bayes_factor(LinearModel, NLogNModel)
  BFLinearPoli2 = bayes_factor(LinearModel, Poli2Model)
  BFLinearExpo = bayes_factor(LinearModel, ExpoModel)
  BFLinear = c(BFLinearConst$bf, BFLinearLog$bf, NaN, 
               BFLinearNLogN$bf,  BFLinearPoli2$bf, BFLinearExpo$bf)
  #NLogN
  BFNLogNConst = bayes_factor(NLogNModel, ConstantModel)
  BFNLogNLog = bayes_factor(NLogNModel, LogModel)
  BFNLogNLinear = bayes_factor(NLogNModel, LinearModel) 
  BFNLogNPoli2 = bayes_factor(NLogNModel, Poli2Model)
  BFNLogNExpo = bayes_factor(NLogNModel, ExpoModel)
  BFNLogN = c(BFNLogNConst$bf, BFNLogNLog$bf, BFNLogNLinear$bf, 
              NaN, BFNLogNPoli2$bf, BFNLogNExpo$bf)
  # Poli2
  BFPoli2Const = bayes_factor(Poli2Model, ConstantModel)
  BFPoli2Log = bayes_factor(Poli2Model, LogModel)
  BFPoli2Linear = bayes_factor(Poli2Model, LinearModel)
  BFPoli2NLogN = bayes_factor(Poli2Model, NLogNModel)
  BFPoli2Expo = bayes_factor(Poli2Model, ExpoModel)
  BFPoli2 = c(BFPoli2Const$bf, BFPoli2Log$bf, BFPoli2Linear$bf, 
              BFPoli2NLogN$bf,  NaN, BFPoli2Expo$bf)
  # Expo
  BFExpoConst = bayes_factor(ExpoModel, ConstantModel)
  BFExpoLog = bayes_factor(ExpoModel, LogModel)
  BFExpoLinear = bayes_factor(ExpoModel, LinearModel)
  BFExpoNLogN = bayes_factor(ExpoModel, NLogNModel)
  BFExpoPoli2 = bayes_factor(ExpoModel, Poli2Model)
  BFExpo = c(BFExpoConst$bf, BFExpoLog$bf, BFExpoLinear$bf, 
             BFExpoNLogN$bf,  BFExpoPoli2$bf, NaN)
  #define the row and column names
  rown = c("Constant", "Log", "Linear", "NLogN", "Poli2", "Expo")
  coln = c("Constant", "Log", "Linear", "NLogN", "Poli2", "Expo")
  
  # create a matrix with all BF values
  mBF <- matrix(c(BFConst, BFLog, BFLinear, BFNLogN, BFPoli2, BFExpo), 
                nrow = 6, byrow = TRUE, 
                dimnames = list(rown, coln))
  return(mBF)
}


#########################################################
# prepare data
############################################################

# Only use the correct and the data that meet the rt cutoff criteria
valid_data <- subset(d, rt <= cutoff & MatchRTMemory <= cutoff & correct == 1 & MatchCorrectMemory == 1 & Condition == "Sort" & queryRT <= cutoff)
# the match Correct and matchRt are needed, because otherwise we allow high rts for the memory condition to influence our results

# get summarized data
dScalingSummary <- ddply(valid_data, ~ Number_of_Bars + Condition + Structure + Stimulus_type, 
                               summarize, mu = mean(RTSort_Memory), se=se(RTSort_Memory))


#########################################################################

#########################################################################


setwd(here("analyses", "brmModel"))

###################
# the following code is still buggy, because the diff of diffs have issues, 
# because some participant had to long rts for the higher sequence lengths
# for now, to solve this, I only took participants that had data for all values
# that is             n
# no structure:       57
# query structure:    59
# sequence structure  53
# see: 
hist(subset(valid_data, Stimulus_type == "bars")$Number_of_Bars)


######################################################################################
# prepare data for the diff of diff analysis and plot
##############################################################################

###############################
#  structure (valid_data)
###############################

dScaling <- ddply(valid_data, 
                        ~TrialID + subject_id + Number_of_Bars + Structure, summarize, mu = mean(RTSort_Memory))


#get the mean over n per subject
dScaling <- ddply(dScaling, ~ Number_of_Bars + subject_id + Structure, summarize, mu = mean(mu))
#initialize diff frame
ddiff <- data.frame(id = numeric(), subject_id = numeric(), name = character(), Structure = character())
#go through all subjects
for (i in 1:length(unique(dScaling$subject_id))){
  
  for (Struc in unique(dScaling$Structure)){
    # set the current participant
    current_participant = subset(dScaling, subject_id == unique(dScaling$subject_id)[i] & Structure == Struc)
    
    # for now to get rid of the weird bug, lets exclude all participants, that don't have means for all seven sequence lengths
    if (length(current_participant$mu) == 7){
      #assign id and get difference (diff is a function that returns lagged differences of a vector)
      ddiff<-rbind(ddiff, data.frame(id = current_participant$subject_id[1], 
                                                 diff = diff(current_participant$mu), 
                                                 name = c("2-1","3-2","4-3","5-4","6-5","7-6"),
                                                 Structure = Struc))
      # the result is for each participant: c(meanrt_2_bars -meanmeanrt_1_bars, meanrt_3_bars-meanrt_2_bars,...)
    }
  }
}

#new n per participant
ddiff$n <- ave(ddiff$diff, ddiff$id, ddiff$Structure, FUN = seq_along)


##########################################################################################################################################################
# start BRM scaling analysis
##########################################################################################################################################

########################################################################
# No structure condition all scaling models
#########################################################################
setwd(here("analyses", "brmModel"))

dScalingNone = subset(dScaling, Structure == "None")

###################################################
# Looking at scaling directly via log log analysis
###################################################

#scaling log log model
mScalingNoneLogLog <- brm(log10(mu+10) ~ log(Number_of_Bars) + (log(Number_of_Bars)|subject_id), 
                          data = dScalingNone, 
                          iter = 20000, # for 10000 the ESS was to low
                          cores = 4,
                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                          save_pars = save_pars(all = TRUE),
                          seed = seed,
                          file = "mScalingNoneLogLog10")
looR2LogLog <- loo_R2(mScalingNoneLogLog)
# results 27.06.2022
#    Estimate  Est.Error      Q2.5     Q97.5
# R2 0.2316741 0.04602355 0.1385132 0.3183788

summary(mScalingNoneLogLog)
# results 27.06.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             1.00      0.00     0.99     1.00 1.00    62019    29421
# logNumber_of_Bars     0.02      0.00     0.01     0.03 1.00    26273    29855

##################################################################
# compare different scaling factors directly
##################################################################

ScalingNoneModelResults = data.frame(ScalingTime = dScalingNone$mu,
                                           Model = "human Data",
                                           NoB = dScalingNone$Number_of_Bars)

# Constant Model
mScalingNoneConstant <- brm(mu ~ 1 + (Number_of_Bars|subject_id), 
                              data = dScalingNone,
                              iter = 10000,
                              cores = 4,
                              control = list(max_treedepth = 15, adapt_delta = 0.99),
                              save_pars = save_pars(all = TRUE),
                              seed = seed,
                              file = "mScalingNoneConstant")
looR2Constant <- loo_R2(mScalingNoneConstant)
# results 27.06.2022
#     Estimate  Est.Error     Q2.5     Q97.5
# R2 0.222554 0.04919648 0.122872 0.3161753

summary(mScalingNoneConstant)
# results 27.06.2022
# Intercept     0.21      0.07     0.06     0.35 1.00    15045    14421

ConstantModelEstimatesNone = fitted(mScalingNoneConstant)
ModelRTConstant = data.frame(ScalingTime = ConstantModelEstimatesNone[,1],
                             Model = "Constant Model",
                             NoB = dScalingNone$Number_of_Bars)

# Linear scaling
mScalingNoneLinear <- brm(mu ~ Number_of_Bars + (Number_of_Bars|subject_id), 
                          data = dScalingNone, 
                          iter = 10000,
                          cores = 4,
                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                          save_pars = save_pars(all = TRUE),
                          seed = seed,
                          file = "mScalingNoneLinear")
looR2Linear <- loo_R2(mScalingNoneLinear)
# results 27.06.2022
#     Estimate  Est.Error      Q2.5     Q97.5
# R2 0.2370997 0.04997404 0.1348706 0.3301106

summary(mScalingNoneLinear)
# results 27.06.2022
# Intercept         -0.12      0.09    -0.30     0.05 1.00    25714    15775
# Number_of_Bars     0.16      0.03     0.10     0.23 1.00    15589    14861

LinearModelEstimatesNone = fitted(mScalingNoneLinear)
ModelRTLinear = data.frame(ScalingTime = LinearModelEstimatesNone[,1],
                             Model = "linear Model",
                           NoB = dScalingNone$Number_of_Bars)


#scaling log model
dScalingNone$NoBLog <- log10(dScalingNone$Number_of_Bars)
mScalingNoneLog <- brm(mu ~ NoBLog + (NoBLog|subject_id), 
                                data = dScalingNone, 
                                iter = 20000, # for 10000 the ESS was to low
                                cores = 4,
                                control = list(max_treedepth = 15, adapt_delta = 0.99),
                                save_pars = save_pars(all = TRUE),
                                seed = seed,
                                file = "mScalingNoneLog10")
looR2Log <- loo_R2(mScalingNoneLog)
# results 27.06.2022
#    Estimate  Est.Error      Q2.5     Q97.5
# R2 0.2499941 0.04499141 0.1580782 0.3350307

summary(mScalingNoneLog)
# results 27.06.2022
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -0.09      0.09    -0.26     0.08 1.00    65619    30110
# NoBLog        1.17      0.21     0.75     1.59 1.00    27052    29731

LogModelEstimatesNone = fitted(mScalingNoneLog)
ModelRTLog = data.frame(ScalingTime = LogModelEstimatesNone[,1],
                             Model = "Log Model",
                        NoB = dScalingNone$Number_of_Bars)

#scaling log model
dScalingNone$NoBNLogN <- log10(dScalingNone$Number_of_Bars)*dScalingNone$Number_of_Bars
mScalingNoneNLogN <- brm(mu ~ NoBNLogN + (NoBNLogN|subject_id), 
                             data = dScalingNone, 
                             iter = 10000,
                             cores = 4,
                             control = list(max_treedepth = 15, adapt_delta = 0.99),
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             file = "mScalingNoneNLogN10")
looR2NLogN <- loo_R2(mScalingNoneNLogN)
# results 27.06.2022
#     Estimate  Est.Error      Q2.5     Q97.5
# R2 0.2209895 0.04997552 0.1193938 0.3144741

summary(mScalingNoneNLogN)
# results 27.06.2022
#            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.11      0.07    -0.03     0.24 1.00    25040    15748
# NoBNLogN      0.16      0.03     0.09     0.23 1.00    12377    15109


NLogNModelEstimatesNone = fitted(mScalingNoneNLogN)
ModelRTNLogN = data.frame(ScalingTime = NLogNModelEstimatesNone[,1],
                             Model = "NLogN Model",
                          NoB = dScalingNone$Number_of_Bars)

dScalingNone$NoBPoli2 <- dScalingNone$Number_of_Bars*dScalingNone$Number_of_Bars
#scaling Plinomiyal model
mScalingNonePoli2 <- brm(mu ~ NoBPoli2 + (NoBPoli2|subject_id), 
                                data = dScalingNone, 
                                iter = 10000,
                                cores = 4,
                                control = list(max_treedepth = 15, adapt_delta = 0.99),
                                save_pars = save_pars(all = TRUE),
                                seed = seed,
                                file = "mScalingNonePoli2")
looR2Poli2 <- loo_R2(mScalingNonePoli2)
# results 27.06.2022
#      Estimate  Est.Error      Q2.5     Q97.5
# R2 0.2037354 0.05179009 0.0972519 0.3014786

summary(mScalingNonePoli2)
# results 27.06.2022
#                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept          0.15      0.07     0.01     0.28 1.00    22762    15383
# NoBPoli2           0.02      0.00     0.01     0.03 1.00    10643    13732

Poli2ModelEstimatesNone = fitted(mScalingNonePoli2)
ModelRTPoli2 = data.frame(ScalingTime = Poli2ModelEstimatesNone[,1],
                          Model = "Poli2 Model",
                          NoB = dScalingNone$Number_of_Bars)


dScalingNone$NoBexponential <- exp(dScalingNone$Number_of_Bars)
#scaling Exponential model
mScalingNoneExponential <- brm(mu ~ NoBexponential + (NoBexponential|subject_id), 
                                     data = dScalingNone, 
                                     iter = 10000,
                                     cores = 4,
                                     control = list(max_treedepth = 15, adapt_delta = 0.99),
                                     save_pars = save_pars(all = TRUE),
                                     seed = seed,
                                     file = "mScalingNoneExponential")
looR2Exponential <- loo_R2(mScalingNoneExponential)
# results 27.06.2022
#    Estimate  Est.Error     Q2.5     Q97.5
# R2 0.1868809 0.04303111 0.100193 0.2678115

summary(mScalingNoneExponential)
# results 27.06.2022
#                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept          0.35      0.07     0.21     0.50 1.00    12447    15103
# NoBexponential     0.00      0.00     0.00     0.00 1.00    14082    14931

ExponentialModelEstimatesNone = fitted(mScalingNoneExponential)
ModelRTExponential = data.frame(ScalingTime = ExponentialModelEstimatesNone[,1],
                          Model = "Exponential Model",
                          NoB = dScalingNone$Number_of_Bars)

ScalingNoneModelResults = rbind(ScalingNoneModelResults,
                                      ModelRTConstant,
                                      ModelRTLog,
                                      ModelRTLinear,
                                      ModelRTNLogN,
                                      ModelRTPoli2,
                                      ModelRTExponential)
ScalingNoneModelResults$Structure = "None"
ModelSummary = ddply(ScalingNoneModelResults, ~ NoB + Model, summarize, MeanScalingTime = mean(ScalingTime))

# plot to compare behavior to model predictions
ggplot(ModelSummary , aes(x = NoB, y = MeanScalingTime, color = Model))+
  geom_point()+
  geom_line()

LooValues <- data.frame(Scaling = c("1Constant", "2Log", "3Linear", "4NLogN", "5Poli2", "6Exponential"),
                        loo_R2 = c(looR2Constant[1], looR2Log[1], looR2Linear[1], 
                                   looR2NLogN[1], looR2Poli2[1], looR2Exponential[1]),
                        Structure = "None")
# results 27.06.2022
#        Scaling    loo_R2 Structure
# 1    1Constant 0.2222912      None
# 2         2Log 0.2501398      None
# 3      3Linear 0.2366419      None
# 4       4NLogN 0.2200753      None
# 5       5Poli2 0.2037007      None
# 6 6Exponential 0.1868994      None

################################
# comparing the different scaling factors
###################################

mNoneBF <- BFCalulcation(mScalingNoneConstant, mScalingNoneLog, mScalingNoneLinear, 
                          mScalingNoneNLogN, mScalingNonePoli2, mScalingNoneExponential)
# results 27.06.2022
#               Constant          Log       Linear        NLogN        Poli2         Expo
# Constant          NaN 2.064729e-05 8.378649e-05 1.893300e-04 1.950936e-01 1.476547e+10
# Log      3.908458e+04          NaN 3.375056e+00 7.344518e+00 6.495618e+03 6.639973e+14
# Linear   9.618832e+03 3.053886e-01          NaN 2.113267e+00 1.808159e+03 1.690864e+14
# NLogN    5.229147e+03 1.228315e-01 5.120023e-01          NaN 9.911016e+02 9.293710e+13
# Poli2    5.975385e+00 1.701025e-04 4.554102e-04 9.112205e-04          NaN 7.154843e+10
# Expo     6.870846e-11 2.216095e-15 5.935479e-15 1.294346e-14 1.236235e-11          NaN

plot(log10(mNoneBF), key=list(tick=FALSE), 
     digits = 2, breaks = range(-21,21), cex = 0.9,
     col = brewer.pal(name = "RdBu", n = 10),
     border = NA)

BFNonePlot = corrplot(log10(mNoneBF), type="full", tl.col="black", tl.srt=45, is.corr = FALSE, 
         col = brewer.pal(name = "RdBu", n = 10),
         col.lim = c(-21, 21),
         number.cex = 0.8,
         na.label = "-")

#################################

######################################
# diff of diff
#########################################
dScalingNoneDiff = subset(ddiff, Structure == "None")

#let's assess sub, super, or linear?
#constant diff
mScalingNoneConstantDiff <- brm(diff ~ 1 + (n|id), 
                                  data = dScalingNoneDiff,
                                  iter = 10000,
                                  cores = 4,
                                  control = list(max_treedepth = 15, adapt_delta = 0.99),
                                  save_pars = save_pars(all = TRUE),
                                  seed = seed,
                                  file = "mScalingNoneNullDiff")

summary(mScalingNoneConstantDiff)
# results 09.08.2022
# Intercept     0.16      0.07     0.03     0.30 1.00    27753    14318


#sub or super
mScalingNoneDiff <- brm(diff ~ n + (n|id), 
                              data = dScalingNoneDiff, 
                              iter = 10000,
                              cores = 4,
                              control = list(max_treedepth = 15, adapt_delta = 0.99), # the delta is supposed to help with the validity issue due to the divergant transitions after warmup
                              save_pars = save_pars(all = TRUE),
                              seed = seed,
                              file = "mScalingNoneDiff")

summary(mScalingNoneDiff)
# results 09.08.2022
# Intercept     0.27      0.16    -0.04     0.57 1.00    27981    15405
# n            -0.03      0.04    -0.11     0.05 1.00    26993    15265



BFDiffDiffLinearConstantNone <- bayes_factor(mScalingNoneDiff, mScalingNoneConstantDiff)
# # results 29.08.2022 0.13174
# evidence for lineare scaling since adding the number of bars does not improve the diff of diff model

BFDiffDiffConstantLinearNone <- bayes_factor(mScalingNoneConstantDiff, mScalingNoneDiff)
# results 29.08.2022 7.53848
# evidence for lineare scaling since adding the number of bars does not improve the diff of diff model





################################################################
# Query condition
####################################################################
setwd(here("analyses", "brmModel"))

dScalingQuery = subset(dScaling, Structure == "Query")

##################################################
# directly looking at the scaling factor (log Log)
##################################################

#scaling log log model
mScalingQueryLogLog <- brm(log10(mu+10) ~ log10(Number_of_Bars) + (log10(Number_of_Bars)|subject_id), 
                           data = dScalingQuery, 
                           iter = 10000,
                           cores = 4,
                           control = list(max_treedepth = 15, adapt_delta = 0.99),
                           save_pars = save_pars(all = TRUE),
                           seed = seed,
                           file = "mScalingQueryLogLog10")
looR2LogLog <- loo_R2(mScalingQueryLogLog)
# results 27.06.2022
#    Estimate  Est.Error      Q2.5     Q97.5
# R2 0.3933627 0.05094206 0.2860824 0.4863637

summary(mScalingQueryLogLog)
# results 27.06.2022
#                      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept               1.00      0.00     0.99     1.01 1.00    22333    14962
# log10Number_of_Bars     0.03      0.01     0.01     0.06 1.00    15779    15211



###################################################################
# comparing different scaling factors
###################################################################

ScalingQueryModelResults = data.frame(ScalingTime = dScalingQuery$mu,
                                           Model = "human Data",
                                           NoB = dScalingQuery$Number_of_Bars)
# constant Model
mScalingQueryConstant <- brm(mu ~ 1 + (Number_of_Bars|subject_id), 
                               data = dScalingQuery, 
                               iter = 10000,
                               cores = 4,
                               control = list(max_treedepth = 15, adapt_delta = 0.99),
                               save_pars = save_pars(all = TRUE),
                               seed = seed,
                               file = "mScalingQueryConstant")
looR2Constant <- loo_R2(mScalingQueryConstant)
# results 27.06.2022
#    Estimate  Est.Error      Q2.5     Q97.5
# R2 0.4066043 0.05467989 0.2924912 0.5054329

summary(mScalingQueryConstant)
# results 27.06.2022
# Intercept     0.17      0.08     0.01     0.32 1.00    24920    14378

ConstantModelEstimatesQuery = fitted(mScalingQueryConstant)
ModelRTConstantQuery = data.frame(ScalingTime = ConstantModelEstimatesQuery[,1],
                             Model = "Constant Model",
                             NoB = dScalingQuery$Number_of_Bars)

# linear model
mScalingQueryLinear <- brm(mu ~ Number_of_Bars + (Number_of_Bars|subject_id), 
                           data = dScalingQuery,
                           iter = 10000,
                           cores = 4,
                           control = list(max_treedepth = 15, adapt_delta = 0.99),
                           save_pars = save_pars(all = TRUE),
                           seed = seed,
                           file = "mScalingQueryLinear")
looR2Linear <- loo_R2(mScalingQueryLinear)
# results 27.06.2022
#    Estimate Est.Error      Q2.5   Q97.5
# R2 0.4135506 0.0549111 0.2965552 0.51417

summary(mScalingQueryLinear)
# results 27.06.2022
# Intercept         -0.09      0.11    -0.31     0.13 1.00    18303    14569
# Number_of_Bars     0.14      0.04     0.06     0.22 1.00    12410    14019


LinearModelEstimatesQuery = fitted(mScalingQueryLinear)
ModelRTLinearQuery = data.frame(ScalingTime = LinearModelEstimatesQuery[,1],
                           Model = "linear Model",
                           NoB = dScalingQuery$Number_of_Bars)



#scaling log model
dScalingQuery$NoBLog <- log10(dScalingQuery$Number_of_Bars)
mScalingQueryLog <- brm(mu ~ NoBLog + (NoBLog|subject_id), 
                             data = dScalingQuery, 
                             iter = 10000,
                             cores = 4,
                             control = list(max_treedepth = 15, adapt_delta = 0.99),
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             file = "mScalingQueryLog10")
looR2Log <- loo_R2(mScalingQueryLog)
# results 27.06.2022
#    Estimate  Est.Error      Q2.5     Q97.5
# R2 0.4031513 0.04829689 0.3019266 0.4899983

summary(mScalingQueryLog)
# results 27.06.2022
#             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -0.06      0.10    -0.25     0.14 1.00    24872    15586
# NoBLog        1.00      0.28     0.44     1.55 1.00    14119    15518


LogModelEstimatesQuery = fitted(mScalingQueryLog)
ModelRTLogQuery = data.frame(ScalingTime = LogModelEstimatesQuery[,1],
                        Model = "Log Model",
                        NoB = dScalingQuery$Number_of_Bars)


#scaling log model
dScalingQuery$NoBNLogN <- log10(dScalingQuery$Number_of_Bars)*dScalingQuery$Number_of_Bars
mScalingQueryNLogN <- brm(mu ~ NoBNLogN + (NoBNLogN|subject_id), 
                               data = dScalingQuery, 
                               iter = 10000,
                               cores = 4,
                               control = list(max_treedepth = 15, adapt_delta = 0.99),
                               save_pars = save_pars(all = TRUE),
                               seed = seed,
                               file = "mScalingQueryNLogN10")
looR2NLogN <- loo_R2(mScalingQueryNLogN)
# results 27.06.2022
#     Estimate  Est.Error      Q2.5     Q97.5
# R2 0.4096006 0.05567487 0.2933641 0.5120582

summary(mScalingQueryNLogN)
# results 27.06.2022
#            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.10      0.08    -0.05     0.25 1.00    21199    14675
# NoBNLogN      0.14      0.05     0.05     0.23 1.00     9183    12188


NLogNModelEstimatesQuery = fitted(mScalingQueryNLogN)
ModelRTNLogNQuery = data.frame(ScalingTime = NLogNModelEstimatesQuery[,1],
                          Model = "NLogN Model",
                          NoB = dScalingQuery$Number_of_Bars)

#scaling Polinomiyal model
dScalingQuery$NoBPoli2 <- dScalingQuery$Number_of_Bars*dScalingQuery$Number_of_Bars
mScalingQueryPoli2 <- brm(mu ~ NoBPoli2 + (NoBPoli2|subject_id), 
                               data = dScalingQuery, 
                               iter = 10000,
                               cores = 4,
                               control = list(max_treedepth = 15, adapt_delta = 0.99),
                               save_pars = save_pars(all = TRUE),
                               seed = seed,
                               file = "mScalingQueryPoli2")
looR2Poli2 <- loo_R2(mScalingQueryPoli2)
# results 27.06.2022
#     Estimate  Est.Error      Q2.5     Q97.5
# R2 0.4015291 0.05757447 0.2809284 0.5069607

summary(mScalingQueryPoli2)
# results 27.06.2022
#             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.14      0.08    -0.02     0.29 1.00    21199    15688
# NoBPoli2      0.02      0.01     0.01     0.03 1.00     9152    12520

Poli2ModelEstimatesQuery = fitted(mScalingQueryPoli2)
ModelRTPoli2Query = data.frame(ScalingTime = Poli2ModelEstimatesQuery[,1],
                          Model = "Poli2 Model",
                          NoB = dScalingQuery$Number_of_Bars)

#scaling Exponential model
dScalingQuery$NoBexponential <- exp(dScalingQuery$Number_of_Bars)
mScalingQueryExponential <- brm(mu ~ NoBexponential + (NoBexponential|subject_id), 
                                     data = dScalingQuery, 
                                     iter = 10000,
                                     cores = 4,
                                     control = list(max_treedepth = 15, adapt_delta = 0.99),
                                     save_pars = save_pars(all = TRUE),
                                     seed = seed,
                                     file = "mScalingQueryExponential")
looR2Exponential <- loo_R2(mScalingQueryExponential)
# results 27.06.2022
#     Estimate  Est.Error      Q2.5     Q97.5
# R2 0.3727436 0.04746024 0.2776012 0.4648028

summary(mScalingQueryExponential)
# results 27.06.2022
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept          0.31      0.09     0.13     0.49 1.00     8866    13162
# NoBexponential     0.00      0.00     0.00     0.00 1.00    13018    14860



ExponentialModelEstimatesQuery = fitted(mScalingQueryExponential)
ModelRTExponentialQuery = data.frame(ScalingTime = ExponentialModelEstimatesQuery[,1],
                                Model = "Exponential Model",
                                NoB = dScalingQuery$Number_of_Bars)


ScalingQueryModelResults = rbind(ScalingQueryModelResults,
                                      ModelRTConstantQuery,
                                      ModelRTLogQuery,
                                      ModelRTLinearQuery,
                                      ModelRTNLogNQuery,
                                      ModelRTPoli2Query,
                                      ModelRTExponentialQuery)
ScalingQueryModelResults$Structure = "Query"
ModelSummary = ddply(ScalingQueryModelResults, ~ NoB + Model, summarize, MeanScalingTime = mean(ScalingTime))
ggplot(ModelSummary , aes(x = NoB, y = MeanScalingTime, color = Model))+
  geom_point()+
  geom_line()


LooValuesQuery <- data.frame(Scaling = c("1Constant", "2Log", "3Linear", "4NLogN", "5Poli2", "6Exponential"),
                        loo_R2 = c(looR2Constant[1], looR2Log[1], looR2Linear[1], 
                                   looR2NLogN[1], looR2Poli2[1], looR2Exponential[1]),
                        Structure = "Query")

LooValues <- rbind(LooValues, LooValuesQuery)

################################
# comparing the different scaling factors
###################################


mQueryBF <- BFCalulcation(mScalingQueryConstant, mScalingQueryLog, mScalingQueryLinear, 
                        mScalingQueryNLogN, mScalingQueryPoli2, mScalingQueryExponential)
# results 28.06.2022
#              Constant          Log       Linear        NLogN        Poli2         Expo
# Constant          NaN 2.865688e+01 4.754937e-02 1.425725e-02 2.512656e+01 3.249771e+18
# Log      4.425142e-02          NaN 8.239216e-04 1.356278e-03 8.221015e-01 4.026003e+16
# Linear   2.320359e+01 9.496584e+02          NaN 2.351841e-01 7.111223e+02 9.550608e+19
# NLogN    1.119982e+02 2.097981e+03 3.039238e+00          NaN 2.137500e+03 2.997682e+20
# Poli2    3.394097e-02 1.058622e+00 1.521573e-03 4.481065e-04          NaN 1.208803e+17
# Expo     3.118529e-19 7.775989e-18 1.158252e-20 3.284781e-21 7.123575e-18          NaN

plot(log10(mQueryBF), key=list(tick=FALSE), 
     digits = 2, breaks = range(-21,21), cex = 0.9,
     col = brewer.pal(name = "RdBu", n = 10),
     border = NA)

corrplot(log10(mQueryBF), type = "full", tl.col = "black", tl.srt = 45, is.corr = FALSE, 
         col = brewer.pal(name = "RdBu", n = 10),
         col.lim = c(-21, 21),
         number.cex = 0.8,
         na.label = "-")

#################################
# diff of diff
####################################
dScalingQueryDiff = subset(ddiff, Structure == "Query")

#test for constant scaling
mScalingQueryDiffNull <- brm(diff ~ 1 + (n|id), 
                                   data = dScalingQueryDiff, 
                                   iter = 10000,
                                   cores = 4,
                                   control = list(max_treedepth = 15, adapt_delta = 0.99),
                                   save_pars = save_pars(all = TRUE),
                                   seed = seed,
                                   file = "mScalingQueryDiffNull")

summary(mScalingQueryDiffNull)
# results 09.08.2022
# Intercept     0.16      0.07     0.02     0.30 1.00    31464    15497


mScalingQueryDiff <- brm(diff ~ n + (n|id), 
                               data = dScalingQueryDiff, 
                               iter = 10000,
                               cores = 4,
                               control = list(max_treedepth = 15, adapt_delta = 0.99), # the delta is supposed to help with the validity issue due to the divergant transitions after warmup
                               save_pars = save_pars(all = TRUE),
                               seed = seed,
                               file = "ScalingQueryDiff")

summary(mScalingQueryDiff)
# results 09.08.2022
# Intercept     0.14      0.16    -0.18     0.46 1.00    30984    14106
# n             0.01      0.04    -0.08     0.09 1.00    31154    14309

BFDiffDiffLinearConstantQuery <- bayes_factor(mScalingQueryDiff, mScalingQueryDiffNull)
# results 29.08.2022: 0.10646

BFDiffDiffConstantLinearQuery <- bayes_factor(mScalingQueryDiffNull, mScalingQueryDiff)
# results 29.08.2022: 9.67007

######################################################################
# Sequence structure condition
####################################################################
setwd(here("analyses","brmModel"))

dScalingSequence = subset(dScaling, Structure == "Sequence")

####################################################################
# calculating scaling factors directly
####################################################################

#scaling log log model
mScalingSequenceLogLog <- brm(log10(mu+10) ~ log10(Number_of_Bars) + (log10(Number_of_Bars)|subject_id), 
                              data = dScalingSequence, 
                              iter = 10000,
                              cores = 4,
                              control = list(max_treedepth = 15, adapt_delta = 0.99),
                              save_pars = save_pars(all = TRUE),
                              seed = seed,
                              file = "mScalingSequenceLogLog10")
looR2LogLog <- loo_R2(mScalingSequenceLogLog)
# results 27.06.2022
#    Estimate  Est.Error       Q2.5     Q97.5
# R2 0.1521627 0.05223227 0.04541176 0.2492254
 
summary(mScalingSequenceLogLog)
# results 27.06.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept               1.00      0.00     0.99     1.01 1.00    32262    14154
# log10Number_of_Bars     0.04      0.01     0.02     0.06 1.00    15231    15656

############################################################
# comparing different scaling models
############################################################

ScalingSequenceModelResults = data.frame(ScalingTime = dScalingSequence$mu,
                                            Model = "human Data",
                                            NoB = dScalingSequence$Number_of_Bars)

# constant model
mScalingSequenceConstant <- brm(mu ~ 1 + (Number_of_Bars|subject_id), 
                                  data = dScalingSequence, 
                                  iter = 10000,
                                  cores = 4,
                                  control = list(max_treedepth = 15, adapt_delta = 0.99),
                                  save_pars = save_pars(all = TRUE),
                                  seed = seed,
                                  file = "mScalingSequenceConstant")
looR2Constant <- loo_R2(mScalingSequenceConstant)
# results 27.06.2022
#     Estimate  Est.Error       Q2.5     Q97.5
# R2 0.1456547 0.05623426 0.03302302 0.2528533

summary(mScalingSequenceConstant)
# results 27.06.2022
# Intercept     0.26      0.07     0.11     0.40 1.00    17648    14850

ConstantModelEstimatesSequence = fitted(mScalingSequenceConstant)
ModelRTConstantSequence = data.frame(ScalingTime = ConstantModelEstimatesSequence[,1],
                                  Model = "Constant Model",
                                  NoB = dScalingSequence$Number_of_Bars)

# linear Model
mScalingSequenceLinear <- brm(mu ~ Number_of_Bars + (Number_of_Bars|subject_id), 
                              data = dScalingSequence,
                              iter = 10000,
                              cores = 4,
                              control = list(max_treedepth = 15, adapt_delta = 0.99), # the delta is supposed to help with the validity issue due to the divergant transitions after warmup
                              save_pars = save_pars(all = TRUE),
                              seed = seed,
                              file = "mScalingSequenceLinear")
looR2Linear <- loo_R2(mScalingSequenceLinear)
# results 27.06.2022
#    Estimate  Est.Error       Q2.5     Q97.5
# R2 0.168446 0.05482803 0.05769054 0.2708901

summary(mScalingSequenceLinear)
# results 27.06.2022
# Intercept         -0.08      0.10    -0.28     0.12 1.00    24372    15898
# Number_of_Bars     0.15      0.03     0.09     0.22 1.00    17392    16376


LinearModelEstimatesSequence = fitted(mScalingSequenceLinear)
ModelRTLinearSequence = data.frame(ScalingTime = LinearModelEstimatesSequence[,1],
                                Model = "linear Model",
                                NoB = dScalingSequence$Number_of_Bars)


#scaling log model
dScalingSequence$NoBLog <- log10(dScalingSequence$Number_of_Bars)
mScalingSequenceLog <- brm(mu ~ NoBLog + (NoBLog|subject_id), 
                              data = dScalingSequence, 
                              iter = 10000,
                              cores = 4,
                              control = list(max_treedepth = 15, adapt_delta = 0.99),
                              save_pars = save_pars(all = TRUE),
                              seed = seed,
                              file = "mScalingSequenceLog10")
looR2Log <- loo_R2(mScalingSequenceLog)
# results 27.06.2022
#     Estimate  Est.Error      Q2.5     Q97.5
# R2 0.4031513 0.04829689 0.3019266 0.4899983

summary(mScalingSequenceLog)
# results 27.06.2022
#             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -0.04      0.09    -0.22     0.15 1.00    31635    14526
# NoBLog        1.05      0.22     0.61     1.48 1.00    13852    14311


LogModelEstimatesSequence = fitted(mScalingSequenceLog)
ModelRTLogSequence = data.frame(ScalingTime = LogModelEstimatesSequence[,1],
                             Model = "Log Model",
                             NoB = dScalingSequence$Number_of_Bars)


#scaling NlogN model
dScalingSequence$NoBNLogN <- log10(dScalingSequence$Number_of_Bars)*dScalingSequence$Number_of_Bars
mScalingSequenceNLogN <- brm(mu ~ NoBNLogN + (NoBNLogN|subject_id), 
                                data = dScalingSequence, 
                                iter = 10000,
                                cores = 4,
                                control = list(max_treedepth = 15, adapt_delta = 0.99),
                                save_pars = save_pars(all = TRUE),
                                seed = seed,
                                file = "mScalingSequenceNLogN10")
looR2NLogN <- loo_R2(mScalingSequenceNLogN)
# results 27.06.2022
#     Estimate  Est.Error       Q2.5     Q97.5
# R2 0.1729944 0.05803515 0.05788498 0.2856702

summary(mScalingSequenceNLogN)
# results 27.06.2022
#             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.13      0.07    -0.01     0.27 1.00    28893    15374
# NoBNLogN      0.15      0.04     0.08     0.23 1.00    10488    13357


NLogNModelEstimatesSequence = fitted(mScalingSequenceNLogN)
ModelRTNLogNSequence = data.frame(ScalingTime = NLogNModelEstimatesSequence[,1],
                               Model = "NLogN Model",
                               NoB = dScalingSequence$Number_of_Bars)


#scaling polinomiyal model
dScalingSequence$NoBPoli2 <- dScalingSequence$Number_of_Bars*dScalingSequence$Number_of_Bars
mScalingSequencePoli2 <- brm(mu ~ NoBPoli2 + (NoBPoli2|subject_id), 
                                data = dScalingSequence, 
                                iter = 10000,
                                cores = 4,
                                control = list(max_treedepth = 15, adapt_delta = 0.99),
                                save_pars = save_pars(all = TRUE),
                                seed = seed,
                                file = "mScalingSequencePoli2")
looR2Poli2 <- loo_R2(mScalingSequencePoli2)
# results 27.06.2022

summary(mScalingSequencePoli2)
# results 27.06.2022
#             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.16      0.07     0.02     0.29 1.00    27491    14683
# NoBPoli2      0.02      0.00     0.01     0.03 1.00    12102    15137


Poli2ModelEstimatesSequence = fitted(mScalingSequencePoli2)
ModelRTPoli2Sequence = data.frame(ScalingTime = Poli2ModelEstimatesSequence[,1],
                               Model = "Poli2 Model",
                               NoB = dScalingSequence$Number_of_Bars)

#scaling log log model
dScalingSequence$NoBexponential <- exp(dScalingSequence$Number_of_Bars)
mScalingSequenceExponential <- brm(mu ~ NoBexponential + (NoBexponential|subject_id), 
                                      data = dScalingSequence, 
                                      iter = 10000,
                                      cores = 4,
                                      control = list(max_treedepth = 15, adapt_delta = 0.99),
                                      save_pars = save_pars(all = TRUE),
                                      seed = seed,
                                      file = "mScalingSequenceExponential")
looR2Exponential <- loo_R2(mScalingSequenceExponential)
# results 27.06.2022
#     Estimate  Est.Error       Q2.5     Q97.5
# R2 0.1760787 0.05661102 0.06820906 0.2902152

summary(mScalingSequenceExponential)
# results 27.06.2022
#                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept          0.33      0.07     0.20     0.46 1.00    16770    17063
# NoBexponential     0.00      0.00     0.00     0.00 1.00    16905    15985


ExponentialModelEstimatesSequence = fitted(mScalingSequenceExponential)
ModelRTExponentialSequence = data.frame(ScalingTime = ExponentialModelEstimatesSequence[,1],
                                     Model = "Exponential Model",
                                     NoB = dScalingSequence$Number_of_Bars)


ScalingSequenceModelResults = rbind(ScalingSequenceModelResults,
                                       ModelRTConstantSequence,
                                       ModelRTLogSequence,
                                       ModelRTLinearSequence,
                                       ModelRTNLogNSequence,
                                       ModelRTPoli2Sequence,
                                       ModelRTExponentialSequence)
ScalingSequenceModelResults$Structure = "Sequence"

ModelSummary = ddply(ScalingSequenceModelResults, ~ NoB + Model, summarize, MeanScalingTime = mean(ScalingTime))
ggplot(ModelSummary , aes(x = NoB, y = MeanScalingTime, color = Model))+
  geom_point()+
  geom_line()

###############################
# plot all three conditions together
##############################
ScalingModelResults = rbind(ScalingNoneModelResults, 
                                  ScalingQueryModelResults, 
                                  ScalingSequenceModelResults)
ModelSummary = ddply(ScalingModelResults , ~ NoB + Model + Structure, summarize, MeanScalingTime = mean(ScalingTime))
ggplot(ModelSummary , aes(x = NoB, y = MeanScalingTime, color = Model))+
  geom_point()+
  geom_line(size = 1, alpha = 0.5)+
  geom_line(data = subset(ModelSummary, Model == "human Data"),  
            aes(x = NoB, y = MeanScalingTime, color = Model), size = 1.5, alpha = 0.8)+
  theme_classic()+
  theme(legend.position = "top")+
  facet_grid(cols = vars(Structure))
####################################

LooValuesSequence <- data.frame(Scaling = c("1Constant", "2Log", "3Linear", "4NLogN", "5Poli2", "6Exponential"),
                             loo_R2 = c(looR2Constant[1], looR2Log[1], looR2Linear[1], 
                                        looR2NLogN[1], looR2Poli2[1], looR2Exponential[1]),
                             Structure = "Sequence")

LooValues <- rbind(LooValues, LooValuesSequence)

ggplot(LooValues, aes(x = Scaling, y = loo_R2, color = Structure, fill = Structure))+
  geom_bar(stat = "identity", position = position_dodge())+
  theme_minimal()


################################
# comparing the different scaling factors
###################################

mSeqBF <- BFCalulcation(mScalingSequenceConstant, mScalingSequenceLog, mScalingSequenceLinear, 
                        mScalingSequenceNLogN, mScalingSequencePoli2, mScalingSequenceExponential)
# results 28.06.2022
#              Constant          Log       Linear        NLogN        Poli2         Expo
# Constant          NaN 2.283473e-03 7.202449e-04 7.729877e-05 5.158672e-03     210162.5
# Log      4.255170e+02          NaN 2.672584e-01 4.138816e-02 2.328313e+00   76653948.1
# Linear   1.473105e+03 2.859592e+00          NaN 1.232740e-01 5.551769e+00  287311123.0
# NLogN    1.328612e+04 2.625342e+01 1.117454e+01          NaN 6.230335e+01 3134758267.8
# Poli2    1.819364e+02 5.231047e-01 1.297588e-01 2.393221e-02          NaN   37672324.4
# Expo     4.464677e-06 1.099361e-08 3.031845e-09 3.111157e-10 1.705994e-08          NaN

plot(log10(mSeqBF), key=list(tick=FALSE), 
     digits = 2, breaks = range(-21, 21), cex = 0.9,
     col = brewer.pal(name = "RdBu", n = 10),
     border = NA)

BFSeqPlot = corrplot(log10(mSeqBF), type="full", tl.col="black", tl.srt = 45, is.corr = FALSE, 
         col = brewer.pal(name = "RdBu", n = 10),
         col.lim = c(-21, 21),
         number.cex = 0.8,
         na.label = "-")

#############################
# diff of diff
##############################
dScalingSequenceDiff = subset(ddiff, Structure == "Sequence")

mScalingSequenceDiffNull <- brm(diff ~ 1 + (n|id), 
                                      data = dScalingSequenceDiff, 
                                      iter = 10000,
                                      cores = 4,
                                      control = list(max_treedepth = 15, adapt_delta = 0.99),
                                      save_pars = save_pars(all = TRUE),
                                      seed = seed,
                                      file = "mScalingSequenceDiffNull")
summary(mScalingSequenceDiffNull)
# results 09.08.2022
# Intercept     0.17      0.08     0.02     0.32 1.00    33389    14850

mScalingSequenceDiff <- brm(diff ~ n + (n|id), 
                                  data = dScalingSequenceDiff, 
                                  iter = 10000,
                                  cores = 4,
                                  control = list(max_treedepth = 15, adapt_delta = 0.99), # the delta is supposed to help with the validity issue due to the divergant transitions after warmup
                                  save_pars = save_pars(all = TRUE),
                                  seed = seed,
                                  file = "mScalingSequenceDiff")
summary(mScalingSequenceDiff)
# results 09.08.2022
# Intercept     0.07      0.17    -0.26     0.42 1.00    38355    14578
# n             0.03      0.04    -0.06     0.11 1.00    35368    14212

BFDiffDiffLinearConstantSeq <- bayes_factor(mScalingSequenceDiff, mScalingSequenceDiffNull)
# results 29.08.2022 BF = 0.12985

BFDiffDiffConstantLinearSeq <- bayes_factor(mScalingSequenceDiffNull, mScalingSequenceDiff)
# results 29.08.2022 BF = 7.30704

###########################################################
# Plotting all Diff od Diff bayes factors
################################################################

#define the row and column names
rown = c("Constant vs Linear", "Linear vs Constant")
coln = c("None", "Query", "Sequence")

# create a matrix with all BF values
mDiffDiffBF <- matrix(c(BFDiffDiffConstantLinearNone$bf, BFDiffDiffConstantLinearQuery$bf, BFDiffDiffConstantLinearSeq$bf, 
                BFDiffDiffLinearConstantNone$bf, BFDiffDiffLinearConstantQuery$bf, BFDiffDiffLinearConstantSeq$bf), 
              nrow = 2, byrow = TRUE, 
              dimnames = list(rown, coln))
corrplot(log10(mDiffDiffBF), 
         type = "full", tl.col = "black", tl.srt = 0, is.corr = FALSE, method="shade",
         col = brewer.pal(name = "RdBu", n = 10),
         col.lim = c(-17, 17),
         number.cex = 0.8,
         na.label = "-")
dfDiffDiffBF <- data.frame(BF = c(BFDiffDiffConstantLinearNone$bf, 
                                  BFDiffDiffConstantLinearQuery$bf, 
                                  BFDiffDiffConstantLinearSeq$bf),
                           Structure = c("None", "Query", "sequence"))
#color pallette for three different groups
cbbPalette <- c("grey50", "#0072B2", "#D55E00")
ggplot(dfDiffDiffBF, aes(x = Structure, y = BF, fill = Structure))+
  geom_bar(stat="identity")+
  #change font and text type
  theme(text = element_text(size = 20, family = "sans"))+
  theme_minimal()+
  theme(legend.position = "none")+
  geom_hline(yintercept = 1, linetype = "longdash", color = "gray")+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)

#define the row and column names
coln = c("Constant vs Linear", "Linear vs Constant")
rown = c("None", "Query", "Sequence")

# create a matrix with all BF values
mDiffDiffBF <- matrix(c(BFDiffDiffConstantLinearNone$bf, BFDiffDiffLinearConstantNone$bf,
                        BFDiffDiffConstantLinearQuery$bf, BFDiffDiffLinearConstantQuery$bf, 
                        BFDiffDiffConstantLinearSeq$bf,BFDiffDiffLinearConstantSeq$bf), 
                      nrow = 3, byrow = TRUE, 
                      dimnames = list(rown, coln))
corrplot(log10(mDiffDiffBF), 
         type = "full", tl.col = "black", tl.srt = 45, is.corr = FALSE, 
         col = brewer.pal(name = "RdBu", n = 10),
         col.lim = c(-17, 17),
         number.cex = 0.8,
         na.label = "-")

######################################################################################
# all Structures Together
#################################################################################
#scaling log log model
mScalingLogLog10 <- brm(log10(mu+10) ~ log10(Number_of_Bars) + Structure + (log10(Number_of_Bars) + Structure|subject_id), 
                            data = dScaling, 
                            iter = 10000,
                            cores = 4,
                            control = list(max_treedepth = 15, adapt_delta = 0.99),
                            save_pars = save_pars(all = TRUE),
                            seed = seed,
                            file = "mScalingLogLog10")
looR2AllLogLog10 <- loo_R2(mScalingLogLog10)
# results 27.06.2022
#    Estimate  Est.Error      Q2.5     Q97.5
# R2 0.237759 0.02510598 0.1880203 0.2861738

summary(mScalingLogLog10)
# results 27.06.2022
#                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept               1.00      0.00     0.99     1.01 1.00    11725    14593
# log10Number_of_Bars     0.04      0.01     0.02     0.05 1.00     7072    11122
# StructureQuery         -0.00      0.01    -0.02     0.01 1.00     8508    11993
# StructureSequence      -0.00      0.00    -0.01     0.01 1.00    11349    13945


##############################################################################
# Doing the whole scaling analysis for all Structure together
############################################################################

# Constant
mScalingConstant <- brm(mu ~ 1 + Structure + (Number_of_Bars + Structure|subject_id), 
                      data = dScaling, 
                      iter = 10000,
                      cores = 4,
                      control = list(max_treedepth = 15, adapt_delta = 0.99),
                      save_pars = save_pars(all = TRUE),
                      seed = seed,
                      file = "mScalingConstant")
summary(mScalingConstant)
# results 27.06.2022
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.29      0.09     0.11     0.48 1.00     5292     9618
# StructureQuery       -0.19      0.16    -0.51     0.13 1.00     5517     9557
# StructureSequence    -0.02      0.12    -0.27     0.22 1.00     6467    10953

# get the posterior Rts from the model
ConstantModelEstimates = fitted(mScalingConstant)
ModelRTConstant = data.frame(ScalingTime = ConstantModelEstimates[,1],
                        Model = "Constant Model",
                        NoB = dScaling$Number_of_Bars)

# Log
dScaling$NoBLog <- log10(dScaling$Number_of_Bars)
mScalingLog <- brm(mu ~ NoBLog + Structure + (NoBLog + Structure|subject_id), 
                     data = dScaling, 
                     iter = 10000,
                     cores = 4,
                     control = list(max_treedepth = 15, adapt_delta = 0.99),
                     save_pars = save_pars(all = TRUE),
                     seed = seed,
                     file = "mScalingLog10")
summary(mScalingLog)
# results 27.06.2022
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            -0.04      0.09    -0.21     0.14 1.00    15851    15669
# NoBLog                1.07      0.16     0.77     1.38 1.00     9450    12870
# StructureQuery       -0.06      0.13    -0.33     0.20 1.00    13022    14214
# StructureSequence    -0.02      0.10    -0.21     0.18 1.00    15643    15202

# get the posterior Rts from the model
LogModelEstimates = fitted(mScalingLog)
ModelRTLog = data.frame(ScalingTime = LogModelEstimates[,1],
                             Model = "Log Model",
                             NoB = dScaling$Number_of_Bars)

# linear
mScalingLinear <- brm(mu ~ Number_of_Bars + Structure + (Number_of_Bars + Structure|subject_id), 
                      data = dScaling, 
                      iter = 10000,
                      cores = 4,
                      control = list(max_treedepth = 15, adapt_delta = 0.99),
                      save_pars = save_pars(all = TRUE),
                      seed = seed,
                      file = "mScalingLinear")
summary(mScalingLinear)
# results 27.06.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            -0.07      0.10    -0.26     0.11 1.00    14041    16261
# Number_of_Bars        0.15      0.02     0.10     0.20 1.00     9459    12923
# StructureQuery       -0.06      0.13    -0.33     0.20 1.00    11501    12577
# StructureSequence    -0.02      0.10    -0.22     0.19 1.00    14362    14966

# get the posterior Rts from the model
LinearModelEstimates = fitted(mScalingLinear)
ModelRTLinear = data.frame(ScalingTime = LinearModelEstimates[,1],
                           Model = "Linear Model",
                           NoB = dScaling$Number_of_Bars)

# NLogN
dScaling$NoBNLogN <- log10(dScaling$Number_of_Bars)*dScaling$Number_of_Bars
mScalingNLogN <- brm(mu ~ NoBNLogN + Structure + (NoBNLogN + Structure|subject_id), 
                      data = dScaling, 
                      iter = 10000,
                      cores = 4,
                      control = list(max_treedepth = 15, adapt_delta = 0.99),
                      save_pars = save_pars(all = TRUE),
                      seed = seed,
                      file = "mScalingNLogN10")
summary(mScalingNLogN)
# results 27.06.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.14      0.08    -0.02     0.30 1.00    11935    14171
# NoBNLogN              0.15      0.03     0.10     0.20 1.00     6354    11199
# StructureQuery       -0.06      0.14    -0.32     0.21 1.00     9001    12308
# StructureSequence    -0.02      0.10    -0.22     0.18 1.00    11394    14106

# get the posterior Rts from the model
NLogNModelEstimates = fitted(mScalingNLogN)
ModelRTNLogN = data.frame(ScalingTime = NLogNModelEstimates[,1],
                           Model = "NLogN Model",
                           NoB = dScaling$Number_of_Bars)

# Poli2
dScaling$NoBPoli2 <- dScaling$Number_of_Bars*dScaling$Number_of_Bars
mScalingPoli2 <- brm(mu ~ NoBPoli2 + Structure + (NoBPoli2 + Structure|subject_id), 
                   data = dScaling, 
                   iter = 10000,
                   cores = 4,
                   control = list(max_treedepth = 15, adapt_delta = 0.99),
                   save_pars = save_pars(all = TRUE),
                   seed = seed,
                   file = "mScalingPoli2")
summary(mScalingPoli2)
# results 27.06.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.18      0.08     0.01     0.34 1.00    12611    14654
# NoBPoli2              0.02      0.00     0.01     0.02 1.00     7867    10765
# StructureQuery       -0.06      0.13    -0.32     0.20 1.00    10306    12286
# StructureSequence    -0.02      0.10    -0.22     0.18 1.00    12828    14129

# get the posterior Rts from the model
Poli2ModelEstimates = fitted(mScalingPoli2)
ModelRTPoli2 = data.frame(ScalingTime = Poli2ModelEstimates[,1],
                              Model = "Poli2 Model",
                              NoB = dScaling$Number_of_Bars)


# Expo
dScaling$NoBExpo <- exp(dScaling$Number_of_Bars)
mScalingExpo <- brm(mu ~ NoBExpo + Structure + (NoBExpo + Structure|subject_id), 
                     data = dScaling, 
                     iter = 10000,
                     cores = 4,
                     control = list(max_treedepth = 15, adapt_delta = 0.99),
                     save_pars = save_pars(all = TRUE),
                     seed = seed,
                     file = "mScalingExpo")
summary(mScalingExpo)
# results 27.06.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.35      0.08     0.20     0.50 1.00    11226    14956
# NoBExpo               0.00      0.00     0.00     0.00 1.00     8763    12658
# StructureQuery       -0.05      0.14    -0.31     0.22 1.00     7355    11806
# StructureSequence    -0.01      0.10    -0.21     0.18 1.00    10217    12287

# get the posterior Rts from the model
ExpoModelEstimates = fitted(mScalingExpo)
ModelRTExpo = data.frame(ScalingTime = ExpoModelEstimates[,1],
                             Model = "Expo Model",
                             NoB = dScaling$Number_of_Bars)
Human = data.frame(ScalingTime = dScalingNone$mu,
                                           Model = "Human Data",
                                           NoB = dScalingNone$Number_of_Bars)
ScalingModelResults = rbind(Human,
                                  ModelRTConstant,
                                  ModelRTLog,
                                  ModelRTLinear,
                                  ModelRTNLogN,
                                  ModelRTPoli2,
                                  ModelRTExpo)

ModelSummary = ddply(ScalingModelResults, ~ NoB + Model, summarize, MeanScalingTime = mean(ScalingTime))
ggplot(ModelSummary , aes(x = NoB, y = MeanScalingTime, color = Model))+
  geom_point()+
  geom_line()+
  theme_classic()+
  geom_line(data = subset(ModelSummary, Model == "Human Data"),  
                        aes(x = NoB, y = MeanScalingTime, color = Model), size = 1.5, alpha = 0.8)

#####################################################
# BF Calculation
#####################################################
mBFAll <- BFCalulcation(mScalingConstant, mScalingLog, mScalingLinear, 
                         mScalingNLogN, mScalingPoli2, mScalingExpo)

# results 29.08.2022
# Constant          Log       Linear        NLogN        Poli2         Expo
# Constant          NaN 1.362715e-04 6.222138e-07 3.429190e-07 1.861972e-04 7.820836e+10
# Log      8.032897e+03          NaN 4.959072e-03 3.067387e-03 1.089450e+00 6.105754e+14
# Linear   1.415471e+06 2.199709e+02          NaN 5.594098e-01 2.494943e+02 1.206497e+17
# NLogN    3.035440e+06 3.989868e+02 2.154860e+00          NaN 4.537120e+02 2.065289e+17
# Poli2    5.835510e+03 8.085687e-01 4.116539e-03 1.836241e-03          NaN 4.456402e+14
# Expo     1.412471e-11 2.010889e-15 7.546046e-18 4.916933e-18 2.082991e-15          NaN

# find the maximum bayes factor
mBFAllvalid <- mBFAll[which(!is.na(mBFAll))]
mBFAllvalid <- mBFAllvalid[which(mBFAllvalid < Inf)]
maxValue <- max(mBFAllvalid)
minValue <- min(mBFAllvalid)
# set all infs to the max value
mBFAll[mBFAll == Inf] <- maxValue
# set all 0 to the min value
mBFAll[mBFAll == 0] <- minValue

logmBFAll = log10(mBFAll)
absLogMax = ceiling(max(abs(logmBFAll[which(!is.na(logmBFAll))])))
plot(logmBFAll, key=list(tick=FALSE), 
     digits = 2, breaks = range(-absLogMax, absLogMax), cex = 0.9,
     col = brewer.pal(name = "RdBu", n = 10),
     border = NA)

BFPlot = corrplot(log10(mBFAll), type = "full", tl.col = "black", tl.srt = 45, is.corr = FALSE, 
                      col = brewer.pal(name = "RdBu", n = 10),
                      col.lim = c(-absLogMax, absLogMax),
                      number.cex = 0.8,
                      na.label = "-")
# plot only nLog N
NLogN <- data.frame("M2" = c("Constant","Log","Linear","NLogN","Poli2","Expo"),
                    "BF" = log10(mBFAll["NLogN",1:6]))
NLogNFigure <- ggplot(NLogN, aes(x = M2, y = BF, color = M2, fill = M2))+
  geom_bar(stat="identity")+
  #minimal theme
  theme_minimal()+
  theme(legend.position = "None")+
  #add xlab
  xlab("Models")+
  #add ylab
  ylab("BF")+
  #change font and text type
  theme(text = element_text(size = 20, family = "sans"))
NLogNFigure
#############################################################################
# comparing a model without structure to a model that contains structure for the sorting times
#############################################################################
# linear
mScalingLinear <- brm(mu ~ Number_of_Bars + Structure + (Number_of_Bars + Structure|subject_id), 
                      data = dScaling, 
                      iter = 10000,
                      cores = 4,
                      control = list(max_treedepth = 15, adapt_delta = 0.99),
                      save_pars = save_pars(all = TRUE),
                      seed = seed,
                      file = "mScaling")
summary(mScalingLinear)
# results 28.06.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            -0.07      0.10    -0.26     0.11 1.00    14041    16261
# Number_of_Bars        0.15      0.02     0.10     0.20 1.00     9459    12923
# StructureQuery       -0.06      0.13    -0.33     0.20 1.00    11501    12577
# StructureSequence    -0.02      0.10    -0.22     0.19 1.00    14362    14966

mScalingLinearNoStructure <- brm(mu ~ Number_of_Bars + (Number_of_Bars + Structure|subject_id), 
                      data = dScaling, 
                      iter = 10000,
                      cores = 4,
                      control = list(max_treedepth = 15, adapt_delta = 0.99),
                      save_pars = save_pars(all = TRUE),
                      seed = seed,
                      file = "mScalingNoStructure")
summary(mScalingLinearNoStructure)
# results 27.06.2022
# Intercept         -0.10      0.07    -0.23     0.04 1.00    17316    15483
# Number_of_Bars     0.15      0.02     0.11     0.20 1.00     7964    12144

bayes_factor(mScalingLinear, mScalingLinearNoStructure)
# results 28.06.2022 BF = 0.09289
bayes_factor(mScalingLinearNoStructure, mScalingLinear)
# results 20.05.2022 BF = 12.47824

# no evicence for effects of Structure on the sorting times

#############################################################################
# Plot how the different tranformations would scale
##############################################################################
SL = seq(1,10)
# make a dataframe with all Sl tranformations
constant = data.frame(scaling = "constant",
                      values = SL*0+1,
                      SL = SL)
Log = data.frame(scaling = "Log",
                      values = log10(SL),
                 SL = SL)
Linear = data.frame(scaling = "Linear",
                 values = SL-1,
                 SL = SL)
NLogN = data.frame(scaling = "NLogN",
                 values = SL*log10(SL),
                 SL = SL)
Poli = data.frame(scaling = "Polynomial (2)",
                 values = SL^2,
                 SL = SL)
Expo = data.frame(scaling = "Exponential",
                 values = exp(SL),
                 SL = SL)
SLdf = rbind(Log, Linear, NLogN)
# plot the whole thing
SlScalingPlot = ggplot(SLdf, aes(x = SL, y = values, color = scaling))+
  geom_line()+
  geom_vline(xintercept = 7, linetype="dashed")+
  theme_minimal()+
  xlab("Sequence Lengths")+
  #add ylab
  ylab('Scaled Time')+
  #change fonts
  theme(text = element_text(size = 25, family = "sans"), legend.position = "top")#

ggsave(here("Figures", "SlScaling.png"), SlScalingPlot,  
       width = 6, 
       height = 6)
