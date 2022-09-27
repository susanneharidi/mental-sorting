# Script to generate the figures for the publication

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
library(patchwork)  #https://github.com/thomasp85/patchwork
library(scales)
library(plot.matrix)
library(RColorBrewer)
library(corrplot)
library(ggpubr)
library(png)

#########################################################################
# functions
#########################################################################

#standard error function
se <- function(x){sd(x)/sqrt(length(x))}

# this is a very specific function to do a round robin comparison
# of the BFs of models where the predictor "sequence length" has been 
# transformed to describe different scaling factors
# I used that here, because I needed to make the same comparison for all
# three structure conditions
BFCalulcation <- function(ConstantModel, LogModel, LinearModel, Poli2Model, ExpoModel){
  # Constant
  BFConstLog = bayes_factor(ConstantModel, LogModel)
  BFConstLinear = bayes_factor(ConstantModel, LinearModel)
  BFConstPoli2 = bayes_factor(ConstantModel, Poli2Model)
  BFConstExpo = bayes_factor(ConstantModel, ExpoModel)
  BFConst = c(NaN, BFConstLog$bf, BFConstLinear$bf, 
               BFConstPoli2$bf, BFConstExpo$bf)
  # Log
  BFLogConst = bayes_factor(LogModel, ConstantModel)
  BFLogLinear = bayes_factor(LogModel, LinearModel)
  BFLogPoli2 = bayes_factor(LogModel, Poli2Model)
  BFLogExpo = bayes_factor(LogModel, ExpoModel)
  BFLog = c(BFLogConst$bf, NaN, BFLogLinear$bf, 
            BFLogPoli2$bf, BFLogExpo$bf)
  # Linear
  BFLinearConst = bayes_factor(LinearModel, ConstantModel)
  BFLinearLog = bayes_factor(LinearModel, LogModel)
  BFLinearPoli2 = bayes_factor(LinearModel, Poli2Model)
  BFLinearExpo = bayes_factor(LinearModel, ExpoModel)
  BFLinear = c(BFLinearConst$bf, BFLinearLog$bf, NaN, 
               BFLinearPoli2$bf, BFLinearExpo$bf)

  # Poli2
  BFPoli2Const = bayes_factor(Poli2Model, ConstantModel)
  BFPoli2Log = bayes_factor(Poli2Model, LogModel)
  BFPoli2Linear = bayes_factor(Poli2Model, LinearModel)
  BFPoli2Expo = bayes_factor(Poli2Model, ExpoModel)
  BFPoli2 = c(BFPoli2Const$bf, BFPoli2Log$bf, BFPoli2Linear$bf, 
              NaN, BFPoli2Expo$bf)
  # Expo
  BFExpoConst = bayes_factor(ExpoModel, ConstantModel)
  BFExpoLog = bayes_factor(ExpoModel, LogModel)
  BFExpoLinear = bayes_factor(ExpoModel, LinearModel)
  BFExpoPoli2 = bayes_factor(ExpoModel, Poli2Model)
  BFExpo = c(BFExpoConst$bf, BFExpoLog$bf, BFExpoLinear$bf, 
             BFExpoPoli2$bf, NaN)
  #define the row and column names
  rown = c("Constant", "Log", "Linear", "Poli2", "Expo")
  coln = c("Constant", "Log", "Linear", "Poli2", "Expo")
  
  # create a matrix with all BF values
  mBF <- matrix(c(BFConst, BFLog, BFLinear, BFPoli2, BFExpo), 
                nrow = 5, byrow = TRUE, 
                dimnames = list(rown, coln))
  return(mBF)
}

#########################################################################
# load data and exclude not independent data and people who do not meet the performance threshold
#########################################################################

d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))

################################################################################
# Set some parameters
################################################################################

pd <- position_dodge(0.1)
textsize <- 12

#color pallette for three different groups
cbbPalette <- c("grey50", "#0072B2", "#D55E00")

#other parameters in s
rt_cutoff <- 10
seed <- 2022
alpha <- 0.7


#############################################################################
# Behavioural Results Figure
#############################################################################
textsize <- 16
# raw Rt plot
#get time frame for correct only again
dRT <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

#get summarize data
dRTSummary <- ddply(dRT, ~Number_of_Bars + Condition + Structure + Stimulus_type, summarize, mu = mean(rt), se = se(rt))

#dodged position
pd <- position_dodge(.2)

#start plotting
pRawRT<-ggplot(dRTSummary, aes(x = Number_of_Bars, y = mu, col = Structure)) +
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
  ylab('Mean RT in s (±SE)')+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"))+
  #legend on top
  theme(legend.position = "top")+
  facet_grid(cols = vars(Condition))

#show!
pRawRT

##################################################
# behavioural model results plot
setwd(here("analyses", "brmModel"))
#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

#second regression contains everything
mrtallStructCondition <- brm(rt ~ Number_of_Bars + Structure + Condition + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                             data = dp,
                             iter = 10000,
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             family = exgaussian(),
                             file = "mrtallStructConditionNotLogFinalEx")

Model_Plot <-mcmc_plot(mrtallStructCondition, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
        geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
        xlab("Posterior Weights")+
        theme_minimal()+
        theme(text = element_text(size = textsize, family = "sans"))+
        scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))
Model_Plot
##############################################################
# RT Gain through Structure plot

#color pallette for three different groups
cbbPalette2 <- c("#0072B2", "#D55E00")

#Only use the correct and the data that meet the rt cutoff criteria
dStrucDiff <- subset(d, rt <= rt_cutoff & correct == 1 & MatchRTMemory <= rt_cutoff & MatchCorrectMemory == 1 & Condition == "Sort" & queryRT <= rt_cutoff)

# calculate the differences of the structure conditions
dNone_Query <- ddply(dStrucDiff, ~ TrialID + subject_id + Number_of_Bars + Stimulus_type, 
                     summarize, m = rt[Conditions == "noStructureExp"] - rt[Conditions == "queryStructureExp"])
dNone_Sequence <- ddply(dStrucDiff, ~ TrialID + subject_id + Number_of_Bars + Stimulus_type, 
                        summarize, m = rt[Conditions == "noStructureExp"] - rt[Conditions == "sortStructureExp"])

#means over n +se
dNone_Query2 <- ddply(dNone_Query, ~ Number_of_Bars + Stimulus_type, summarize, mu = mean(m), se = se(m))
dNone_Sequence2 <- ddply(dNone_Sequence, ~ Number_of_Bars + Stimulus_type, summarize, mu = mean(m), se = se(m))
#bind them
dStrucDiffSummary <- rbind(dNone_Query2, dNone_Sequence2)
dStrucDiffSummary$Structure <- rep(c('None - Query', 'None - Sequence'), each = 7)


#start plotting
plotStrucDiffRT <- ggplot(dStrucDiffSummary, aes(x = Number_of_Bars, y = mu, col = Structure)) +
  #points
  geom_point(position = pd)+
  #error bars
  geom_errorbar(aes(ymin = mu-se, ymax = mu+se), width = 0, size = 1, position = pd) +
  #lines
  geom_line(position = pd, size = 1.2) +
  scale_color_manual(values = cbbPalette2)+
  #legend on top
  theme(strip.background = element_blank())+
  #change x ticks
  scale_x_continuous(breaks = round(seq(min(0), max(10), by = 1),1)) +
  #change y ticks
  scale_y_continuous(breaks = round(seq(min(0), max(10), by = 1),1)) +
  #minimal theme
  theme_minimal()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  #add xlab
  xlab("Sequence Length")+
  #add ylab
  ylab('Mean RT Gain in s')+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"))+
  #legend on top
  theme(legend.position = "None")+
  ylim(-0.2,1.2)

plotStrucDiffRT

simulatedData <- data.frame(x = seq(1,7),
                            y = c(seq(1,7)*0.75, seq(1,3)*0.75, rep(3*0.75, 4)),
                            Structure = c(rep("None-Query", 7), rep("None-Seq", 7)))
#start plotting
plotStrucDiffSim <- ggplot(simulatedData, aes(x = x, y = y-0.5, col = Structure)) +
  #lines
  geom_line(position = pd, size = 1.2) +
  scale_color_manual(values = cbbPalette2)+
  #legend on top
  theme(strip.background = element_blank())+
  #change x ticks
  scale_x_continuous(breaks = round(seq(min(0), max(10), by = 1),1)) +
  #change y ticks
  scale_y_continuous(breaks = round(seq(min(0), max(10), by = 1),1)) +
  #minimal theme
  theme_minimal()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  #add xlab
  xlab("Sequence Length")+
  #add ylab
  ylab('RT Gain')+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"))+
  #legend on top
  ylim(-1, 6)+
  theme(legend.position = "top",
        axis.text.y = element_blank())

plotStrucDiffSim

PlotRTGain <- ggarrange(plotStrucDiffSim + rremove("xlab"), plotStrucDiffRT, 
                        ncol = 1, 
                        nrow = 2, 
                        common.legend = TRUE, 
                        legend = "top")
PlotRTGain

BehaviouralFigure <- pRawRT + Model_Plot + (plotStrucDiffSim / plotStrucDiffRT) + 
  plot_annotation(tag_levels = list(c("A", "B", "C")))
BehaviouralFigure


ggsave(here("Figures", "BehaviouralFigure.pdf"), BehaviouralFigure,  
       width = 14, 
       height = 5)

############################################################################
# Scaling figure
################################################################################

#########################################################
# prepare data
#########################################################
textsize <- 16
#Only use the correct and the data that meet the rt cutoff criteria
valid_data <- subset(d, rt <= rt_cutoff & MatchRTMemory <= rt_cutoff & correct == 1 & MatchCorrectMemory == 1 & Condition == "Sort" & queryRT <= rt_cutoff)
# the match Correct and matchRt are needed, because otherwise we allow high rts for the memory condition to influence our results

#get summarize data
dScalingSummaryWithin <- ddply(valid_data, ~ Number_of_Bars + Condition + Structure + Stimulus_type, 
                               summarize, mu = mean(RTSort_Memory), se=se(RTSort_Memory))


######################################################################################
# prepare data for the diff of diff analysis and plot
##############################################################################

###############################
# within structure (valid_data)
###############################

dScalingWithin <- ddply(valid_data, 
                        ~TrialID + subject_id + Number_of_Bars + Structure, summarize, mu = mean(RTSort_Memory))


#get the mean over n per subject
dScalingWithin <- ddply(dScalingWithin, ~ Number_of_Bars + subject_id + Structure, summarize, mu = mean(mu))
#initialize diff frame
ddiffWithin <- data.frame(id = numeric(), subject_id = numeric(), name = character(), Structure = character())
#go through all subjects
for (i in 1:length(unique(dScalingWithin$subject_id))){
  
  for (Struc in unique(dScalingWithin$Structure)){
    # set the current participant
    current_participant = subset(dScalingWithin, subject_id == unique(dScalingWithin$subject_id)[i] & Structure == Struc)
    
    # for now to get rid of the weird bug, lets exclude all participants, that don't have means for all seven sequence lengths
    if (length(current_participant$mu) == 7){
      #assign id and get difference (diff is a function that returns lagged differences of a vector)
      ddiffWithin<-rbind(ddiffWithin, data.frame(id = current_participant$subject_id[1], 
                                                 diff = diff(current_participant$mu), 
                                                 name = c("2-1","3-2","4-3","5-4","6-5","7-6"),
                                                 Structure = Struc))
      # the result is for each participant: c(meanrt_2_bars -meanmeanrt_1_bars, meanrt_3_bars-meanrt_2_bars,...)
    }
  }
}

#new n per participant
ddiffWithin$n <- ave(ddiffWithin$diff, ddiffWithin$id, ddiffWithin$Structure, FUN = seq_along)

#rank for sub or super analysis
ddiffWithin$rank <- ave(ddiffWithin$diff, ddiffWithin$id, ddiffWithin$Structure, FUN = rank)

#####################################################
#Plot Rt difference within
plot_scalingWithin <- ggplot(dScalingSummaryWithin, aes(x=Number_of_Bars, y=mu, col=Structure)) +
  #points
  geom_point(position = pd)+
  #error bars
  geom_errorbar(aes(ymin=mu-se, ymax=mu+se), width=0, size=1, position=pd) +
  #lines
  geom_line(position=pd, size=1.2) +
  scale_color_manual(values=cbbPalette)+
  #legend on top
  theme(strip.background=element_blank())+
  #change x ticks
  scale_x_continuous(breaks = round(seq(min(0), max(10), by = 1),1)) +
  #change y ticks
  scale_y_continuous(breaks = round(seq(min(0), max(10), by = 1),1)) +
  #minimal theme
  theme_minimal()+
  #add xlab
  xlab("Sequence Length")+
  #add ylab
  ylab('Mean RT Difference in s')+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"))+
  #legend on top
  theme(legend.position = "top")+
  ylim(-0.2,1.2)

plot_scalingWithin
###############################################################
# make the diff of diff scaling plot within

ddiffSummaryWithin<-ddply(ddiffWithin, ~n+Structure, summarize, mu=mean(diff), se = se(diff))

# plot this:
#start plotting
plot_scalingDiffWithin<-ggplot(ddiffSummaryWithin, aes(x=n, y=mu, col = Structure)) +
  #points
  geom_point(position =pd)+
  # add line at 0
  geom_hline(yintercept = 0, linetype="dashed", color = "black")+
  #lines
  geom_line(position=pd, size=1.2) +
  scale_color_manual(values=cbbPalette)+
  #error bars
  geom_errorbar(aes(ymin=mu-se, ymax=mu+se), width=0, size=1, position=pd) +
  #legend on top
  theme(strip.background=element_blank())+
  #change x ticks
  scale_x_discrete(name ="Difference of Sequence Lengths", 
                   limits=c("2-1","3-2","4-3","5-4","6-5","7-6"))+
  #change y ticks
  scale_y_continuous(breaks = round(seq(min(0), max(10), by = 1),1)) +
  #minimal theme
  theme_minimal()+
  #add ylab
  ylab('Mean RT Diff of Diff')+
  #change fonts
  theme(text = element_text(size = textsize, family="sans"))+
  #legend on top
  theme(legend.position = "None")+
  ylim(-0.75,0.75)

plot_scalingDiffWithin

#############################################################################
# Bf of scaling factors plot

##############################################################################
# Doing the whole scaling analysis for all Structure together
############################################################################
setwd(here("analyses", "brmModel"))
# Constant
mScalingConstant <- brm(mu ~ 1 + Structure + (Number_of_Bars + Structure|subject_id), 
                        data = dScalingWithin,
                        file = "mScalingConstant")

# Log
mScalingLog <- brm(mu ~ NoBLog + Structure + (NoBLog + Structure|subject_id), 
                   data = dScalingWithin, 
                   file = "mScalingLog10")

# linear
mScalingLinear <- brm(mu ~ Number_of_Bars + Structure + (Number_of_Bars + Structure|subject_id), 
                      data = dScalingWithin, 
                      file = "mScalingLinear")

# Poli2
mScalingPoli2 <- brm(mu ~ NoBPoli2 + Structure + (NoBPoli2 + Structure|subject_id), 
                     data = dScalingWithin, 
                     file = "mScalingPoli2")

# Expo
mScalingExpo <- brm(mu ~ NoBExpo + Structure + (NoBExpo + Structure|subject_id), 
                    data = dScalingWithin, 
                    file = "mScalingExpo")

#####################################################
# BF Calculation
#####################################################
mBFAll <- BFCalulcation(mScalingConstant, mScalingLog, mScalingLinear, 
                        mScalingPoli2, mScalingExpo)
mBFAll
# 29.08.2022
#               Constant          Log       Linear        Poli2         Expo
# Constant          NaN 1.314472e-04 6.760500e-07 1.350176e-04 8.809369e+10
# Log      7.395254e+03          NaN 4.988331e-03 1.121718e+00 6.096429e+14
# Linear   1.636488e+06 1.991883e+02          NaN 2.005533e+02 1.436129e+17
# Poli2    6.216751e+03 8.678909e-01 3.874278e-03          NaN 4.404615e+14
# Expo     1.195163e-11 1.631470e-15 8.122222e-18 1.531652e-15          NaN

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
a <- plot(logmBFAll, key=list(tick=FALSE), 
     digits = 2, breaks = range(-absLogMax, absLogMax), cex = 0.9,
     col = brewer.pal(name = "RdBu", n = 10),
     border = NA)


BFPlot = corrplot(log10(mBFAll), type = "full", tl.col = "black", tl.srt = 45, is.corr = FALSE, 
                  cex.col = 3,
                  col = brewer.pal(name = "RdBu", n = 10),
                  col.lim = c(-absLogMax, absLogMax),
                  number.cex = 2,
                  tl.cex = 1,
                  na.label = "-")


################################################################
# BF for diff of diff analysis
dScalingNoneWithinDiff = subset(ddiffWithin, Structure == "None")

#constant diff
mScalingNoneNullDiff <- brm(diff ~ 1 + (n|id), 
                                  data = dScalingNoneWithinDiff,
                                  file = "mScalingNoneNullDiff")

#sub or super
mScalingNoneDiff <- brm(diff ~ n + (n|id), 
                              data = dScalingNoneWithinDiff, 
                              file = "mScalingNoneDiff")

BFDiffDiffLinearConstantNone <- bayes_factor(mScalingNoneDiff, mScalingNoneNullDiff)
BFDiffDiffConstantLinearNone <- bayes_factor(mScalingNoneNullDiff, mScalingNoneDiff)

# query ####################################

dScalingQueryWithinDiff = subset(ddiffWithin, Structure == "Query")

#test for constant scaling
mScalingQueryDiffNull <- brm(diff ~ 1 + (n|id), 
                                   data = dScalingQueryWithinDiff, 
                                   file = "mScalingQueryDiffNull")

mScalingQueryDiff <- brm(diff ~ n + (n|id), 
                               data = dScalingQueryWithinDiff, 
                               file = "ScalingQueryDiff")


BFDiffDiffLinearConstantQuery <- bayes_factor(mScalingQueryDiff, mScalingQueryDiffNull)
BFDiffDiffConstantLinearQuery <- bayes_factor(mScalingQueryDiffNull, mScalingQueryDiff)

# sequence ################################

dScalingSequenceWithinDiff = subset(ddiffWithin, Structure == "Sequence")

mScalingSequenceDiffNull <- brm(diff ~ 1 + (n|id), 
                                      data = dScalingSequenceWithinDiff, 
                                      file = "mScalingSequenceDiffNull")

mScalingSequenceDiff <- brm(diff ~ n + (n|id), 
                                  data = dScalingSequenceWithinDiff, 
                                  file = "mScalingSequenceDiff")

BFDiffDiffLinearConstantSeq <- bayes_factor(mScalingSequenceDiff, mScalingSequenceDiffNull)
BFDiffDiffConstantLinearSeq <-bayes_factor(mScalingSequenceDiffNull, mScalingSequenceDiff)

##########

dfDiffDiffBF <- data.frame(BF = c(BFDiffDiffConstantLinearNone$bf, 
                                  BFDiffDiffConstantLinearQuery$bf, 
                                  BFDiffDiffConstantLinearSeq$bf),
                           Structure = c("None", "Query", "Sequence"))
#color pallette for three different groups
cbbPalette <- c("grey50", "#0072B2", "#D55E00")
EvidenceForLinearPlot <- ggplot(dfDiffDiffBF, aes(x = Structure, y = BF, fill = Structure))+
            geom_bar(stat = "identity")+
            #change font and text type
            theme(text = element_text(size = textsize, family = "sans"))+
            theme_minimal()+
            ylab('Bayes Factor')+
            theme(legend.position = "none",
                  text = element_text(size = textsize, family="sans"))+
            geom_hline(yintercept = 1, linetype = "longdash", color = "black")+
            #change fill
            scale_fill_manual(values = cbbPalette)+
            #change color
            scale_color_manual(values = cbbPalette)
EvidenceForLinearPlot

# the BF plot needs to be added outside of R, since it is not a ggplot.
# for now I use the plot_scalingWithin figure as a plce holder
ScalingFigure <- (plot_scalingWithin + plot_scalingWithin) / (plot_scalingDiffWithin + EvidenceForLinearPlot)+ 
  plot_annotation(tag_levels = "A")
ScalingFigure 

ggsave(here("Figures", "ScalingFigure.png"), ScalingFigure,  
       width = 9, 
       height = 7)

######################################################################################
# The Model comparison figure is generated in the Model_Comparison script
######################################################################################

######################################################################################
# Model Output Figure
######################################################################################

source(here("Model", "hypoGeneratorFunctions.R")) 

###############################################################################################################
#######################################################################################
# run hypo learner on data
#######################################################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)

learningrate = 0.09
rt_cutoff = 10
#######################################################################################
# set up hypothesis space
#######################################################################################

start_particle = generateStartParticles(nOfBars = 7, nOfColors = 10)

#######################################################################################
# run sorter on all data
#######################################################################################
ModelToData <- runHypoLearnerOnAllData(start_particle,
                                       AllData = d,
                                       learningrate)

#######################################################################################
# plot all thresholds
cbbPalette <- c("grey50","#0072B2", "#D55E00")
pd <- position_dodge(.1)
allfinalHypo <- ModelToData
textsize <- 18

plotFinalThresholds <- ggplot(subset(allfinalHypo, trial == 35), aes(x = threshold, fill = Structure, color = Structure)) +
  #points
  geom_histogram(alpha = .4, position="dodge")+
  #minimal theme
  theme_minimal()+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #add xlab
  xlab("Threshold")+
  #add ylab
  ylab('Runs')+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"), legend.position = "None")+
  xlim(1,7)
  #add title
  #ggtitle("Final Thresholds")

plotFinalThresholds

#######################################################################################
# final connections

# labe the correct connection for the structure condition:
allfinalHypo$correctConnection = "dother 3"
for (i in 1:length(allfinalHypo$connection)){
  if (allfinalHypo$connection[i] == allfinalHypo$trueConnection[i]){
    allfinalHypo$correctConnection[i] = "a3correct"
  }else if(allfinalHypo$connection[i] == substr(allfinalHypo$trueConnection[i], 1, 2)){
    allfinalHypo$correctConnection[i] = "b2correct"
  }else if(allfinalHypo$connection[i] == substr(allfinalHypo$trueConnection[i], 2, 3)){
    allfinalHypo$correctConnection[i] = "b2correct"
  }else if(length(string_to_vec(allfinalHypo$connection[i])) <= 1){
    allfinalHypo$correctConnection[i] = "cno"
  }else if(length(string_to_vec(allfinalHypo$connection[i])) == 2){
    allfinalHypo$correctConnection[i] = "b2other"
  }
}

plotFinalConnections <- ggplot(subset(allfinalHypo, trial == 35) , aes(x = correctConnection, fill = Structure, color = Structure)) +
  #points
  geom_bar(alpha = .4, position = position_dodge2(preserve = "single"))+
  #minimal theme
  theme_minimal()+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels = c("3 correct", " 2 correct", "no Connections"))+
  #add xlab
  xlab("Connections")+
  #add ylab
  ylab('Runs')+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"), legend.position = "None")#+
  #add title
  #ggtitle("Final Connections")
#facet_grid(rows = vars(startAtSmall))

plotFinalConnections

#######################################################################################
# Time over NoB

rtNoB = ddply(allfinalHypo, ~NoB + Structure , summarize, meanTime = mean(time), se = se(time))
plotNoB <- ggplot(rtNoB, aes(x = NoB, y = meanTime, fill = Structure, color = Structure)) +
  geom_line(position = pd) +
  #points
  geom_point(position = pd)+
  #error bars
  geom_errorbar(aes(ymin = meanTime-se, ymax = meanTime+se), width = 0, size = 1, position = pd)+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  theme_minimal()+
  xlab("Sequence Length")+
  #add ylab
  ylab('Mean Model Time')+
  #change color
  scale_color_manual(values = cbbPalette)+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"), legend.position = "top")#+
  #add title
  #ggtitle("Model Time")

plotNoB

#######################################################################################
# time over trial index

rtTrials = ddply(allfinalHypo, ~trial + Structure , summarize, meanTime = mean(time), se = se(time))

plotrtTrials <- ggplot(rtTrials, aes(x = trial, y = meanTime, fill = Structure, color = Structure)) +
  geom_line(position = pd, size = 0.5) +
  #points
  geom_point(position = pd)+
  #error bars
  geom_errorbar(aes(ymin = meanTime-se, ymax = meanTime+se), width = 0, size = 1, position = pd)+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  theme_minimal()+
  xlab("Trial")+
  #add ylab
  ylab('Mean Model Time')+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"), legend.position = "top")#+
  #add title
  #ggtitle("Model Time over Trials")

plotrtTrials

#######################################################################################
# show how well the real rts correspond to the models RTs

# this part is necesarry, so I can add the correct columns to the model dataframe
relevantD <- subset(d, Condition == "Sort")
relevantD1 <- subset(relevantD, Structure == "Query")
relevantD2 <- subset(relevantD, Structure == "Sequence")
relevantD3 <- subset(relevantD, Structure == "None")
relevantD <- rbind(relevantD1, relevantD2, relevantD3)

ModelToData$rawRt <- relevantD$rt
ModelToData$queryRT <- relevantD$queryRT
modelData <- (subset(ModelToData, realAccuracy == 1 & rawRt <= rt_cutoff & queryRT <= rt_cutoff))

ModelvsDataPlot <- ggplot(modelData, aes(x = time, y = rawRt, color = Structure))+
  geom_point(alpha = 0.2)+
  theme(legend.position = "None")+
  geom_smooth(data = modelData, aes(x = time, y = rawRt), method = lm)+
  theme_minimal()+
  scale_x_continuous(breaks = c(0, 100, 200))+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  theme(text = element_text(size = textsize, family = "sans"), legend.position = "none")+
  xlab("Model Time")+
  #add ylab
  ylab('Participant RT in s')+
  #change fonts
  facet_grid(cols = vars(Structure))

ModelvsDataPlot
#######################################################################################
# plotting the model that tries to emulate the effects of participant RT

setwd(here("Model", "hypoFittingBrmModels"))
EffectsModel <- brm(time ~ Structure + NoB + (Structure + NoB|participant),
                    data = modelData,
                    chains = 2,
                    save_pars = save_pars(all = TRUE),
                    control = list(max_treedepth = 15, adapt_delta = 0.99),
                    file = "hypoGeneratorEmulatingRTresultsBucketSort")

ModelofModelPlot <- mcmc_plot(EffectsModel, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("Model Time")+
  theme_minimal()+
  theme(text = element_text(size = textsize, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Query Structure", "Sequence Structure", "Sequence Length"))# 

ModelofModelPlot

ModelOutputPlot <- (plotNoB + ModelofModelPlot + plotrtTrials)/(plotFinalThresholds + plotFinalConnections + ModelvsDataPlot) +
  plot_annotation(tag_levels = "A")+ 
  plot_layout(guides = 'collect')
ModelOutputPlot

ggsave(here("Figures", "ModelOutputFigure.pdf"), ModelOutputPlot,  
       width = 14, 
       height = 7)
