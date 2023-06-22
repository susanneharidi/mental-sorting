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
library(ggbreak) 
library(viridis)

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
# No structure condition second analysis derivative of scaling
#########################################################################
setwd(here("analyses", "brmModel"))

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
# Query condition second analysis derivative of scaling
####################################################################
setwd(here("analyses", "brmModel"))


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
# Sequence structure condition second analysis derivative of scaling
####################################################################
setwd(here("analyses","brmModel"))


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


##############################################################################
# Doing the whole scaling analysis with different f(s) for all Structure together
############################################################################
# prepare vector to gather R2 values
loo_R2s <- c()

# Constant
mScalingConstant <- brm(mu ~ 1 + Structure + (1 + Structure|subject_id), 
                      data = dScaling, 
                      iter = 10000,
                      cores = 4,
                      control = list(max_treedepth = 15, adapt_delta = 0.99),
                      save_pars = save_pars(all = TRUE),
                      seed = seed,
                      file = "mScalingConstant1")
summary(mScalingConstant)

# results 16.06.2023 after correcting the constant model("mScalingConstant1")
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.49      0.08     0.34     0.65 1.00    14224    15552
# StructureQuery       -0.03      0.14    -0.31     0.24 1.00     7498    12859
# StructureSequence    -0.02      0.10    -0.22     0.19 1.00    14776    15347


loo_R2_constant <- loo_R2(mScalingConstant)
loo_R2s <- append(loo_R2s, loo_R2_constant[1])
# results 16.06.2023 after correcting the constant model
#      Estimate  Est.Error     Q2.5     Q97.5
# R2 0.1614404 0.01999896 0.121902 0.2001055

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

loo_R2_Log <- loo_R2(mScalingLog)
loo_R2s <- append(loo_R2s, loo_R2_Log[1])
#    Estimate  Est.Error      Q2.5     Q97.5
#R2 0.2574409 0.02323915 0.2105828 0.3020583

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

loo_R2_Linear <- loo_R2(mScalingLinear)
loo_R2s <- append(loo_R2s, loo_R2_Linear[1])
#    Estimate  Est.Error      Q2.5     Q97.5
#R2 0.2591055 0.02560232 0.2075066 0.3079029

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

loo_R2_NLogN <- loo_R2(mScalingNLogN)
loo_R2s <- append(loo_R2s, loo_R2_NLogN[1])
#    Estimate  Est.Error      Q2.5     Q97.5
#0.259061 0.02655287 0.2055914 0.3098298

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

loo_R2_Poli2 <- loo_R2(mScalingPoli2)
loo_R2s <- append(loo_R2s, loo_R2_Poli2[1])
#    Estimate  Est.Error      Q2.5     Q97.5
#R2 0.2552156 0.02755365 0.1993971 0.3078159

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

loo_R2_Expo <- loo_R2(mScalingExpo)
loo_R2s <- append(loo_R2s, loo_R2_Expo[1])
#    Estimate  Est.Error      Q2.5     Q97.5
#R2 0.2156919 0.02782176 0.1592456 0.2685296

# get the posterior Rts from the model
ExpoModelEstimates = fitted(mScalingExpo)
ModelRTExpo = data.frame(ScalingTime = ExpoModelEstimates[,1],
                             Model = "Expo Model",
                             NoB = dScaling$Number_of_Bars)
Human = data.frame(ScalingTime = dScaling$mu,
                                           Model = "Human Data",
                                           NoB = dScaling$Number_of_Bars)
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

#######################################
# Loo R2 Comparison
Loo_R2s_df = data.frame(Loo_R2 = loo_R2s,
                        Loo_R2_diff = loo_R2s-rep(max(loo_R2s), length(loo_R2s)),
                        Scaling = c("constant","log", "linear", 
                                    "n log n", "polinomiyal", "exponential"))
ggplot(Loo_R2s_df, aes(x = Scaling, y = Loo_R2_diff, fill = Scaling))+
  scale_x_discrete(limits=c("constant", "log", "linear", 
                            "n log n", "polinomiyal", "exponential"))+
  theme_minimal()+
  geom_bar(stat="identity")+
  #add xlab
  xlab("Scaling factor")+
  #add ylab
  ylab("Loo R2 difference")+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "none")
  

#####################################################
# BF Calculation
#####################################################
mBFAll <- BFCalulcation(mScalingConstant, mScalingLog, mScalingLinear, 
                         mScalingNLogN, mScalingPoli2, mScalingExpo)

# results 16.06.23 after fixing the error in the constant model
#              Constant          Log       Linear        NLogN        Poli2         Expo
# Constant          NaN 1.575721e-36 8.221548e-39 3.849751e-39 2.032698e-36 9.771046e-22
# Log      5.392020e+35          NaN 5.855371e-03 2.541319e-03 1.356692e+00 4.767445e+14
# Linear   1.211177e+38 1.603187e+02          NaN 6.475510e-01 2.279002e+02 1.391250e+17
# NLogN    2.302288e+38 4.141968e+02 1.921638e+00          NaN 3.942044e+02 2.358027e+17
# Poli2    5.951126e+35 7.365663e-01 4.015539e-03 2.267465e-03          NaN 5.146436e+14
# Expo     9.333458e+20 1.751376e-15 7.143167e-18 4.143137e-18 2.176766e-15          NaN


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
# Plot how the different tranformations would scale
##############################################################################

SL = seq(1,10)
# make a dataframe with all Sl tranformations
constant = data.frame(scaling = "Constant",
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

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
SlScalingPlot = ggplot(SLdf, aes(x = SL, y = values, color = scaling))+
  geom_line(size = 1.5)+
  scale_color_manual(values = cbPalette)+
  geom_vline(xintercept = 7, linetype="dashed")+
  theme_minimal()+
  xlab("Sequence Lengths")+
  #add ylab
  ylab('Time')+
  labs(color = "Complexity") +
  #change fonts
  theme(text = element_text(size = 25, family = "sans"), legend.position = "top")#

SlScalingPlot

ggsave(here("Figures", "SlScaling.pdf"), SlScalingPlot,  
       width = 10, 
       height = 6)
