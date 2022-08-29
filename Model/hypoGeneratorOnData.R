# Run the hypothesis generator on the experimental data (i.e. the trials the participants saw)
# the output of this script is used for the model selection

################################################################################
# libraries
###############################################################################
library(ggplot2)
library(cowplot)
library(plyr)
library(ggthemes)
library(ggsignif)
library(gridExtra)
library(ggpubr)
library(here)
library(brms)
library(ggpubr)
library(stringr)
library(readr)

source(here("Model", "hypoGeneratorFunctions.R")) 

###############################################################################################################
###############################################################################################################
# run hypo learner on data
####################################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)

learningrate = 0.09

##################################################################################
# set up hypothesis space
###########################################################################

start_particle = generateStartParticles(nOfBars = 7, nOfColors = 10)

###########################
# run sorter on all data
###############################
ModelToData <- runHypoLearnerOnAllData(start_particle,
                                       AllData = d,
                                       learningrate)

#################################
# save model to data dataframe for further analysis
###################################
write_csv(ModelToData, here("Model", "FittedModelData", paste(as.character(learningrate), "hypoGeneratorfitted_BucketSort.csv", sep = "_")))


##################################################################
# plots of model run on real data
###################################################################

################################################################
# plot all thresholds
cbbPalette <- c("grey50","#0072B2", "#D55E00")
pd <- position_dodge(.2)
allfinalHypo = ModelToData

plotFinalThresholds <- ggplot(subset(allfinalHypo, trial == 35) , aes(x = threshold, fill = Structure, color = Structure)) +
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
  theme(text = element_text(size = 20, family = "sans"))+
  xlim(1,7)+
  #add title
  ggtitle("Final Thresholds")

print(plotFinalThresholds)

thresholdsTime = ddply(allfinalHypo, ~trial + Structure , summarize, meanthreshold = mean(threshold), se = se(threshold))
plotThresholdsoverTime <- ggplot(thresholdsTime, aes(y = meanthreshold, x = trial , fill = Structure, color = Structure)) +
  #points
  geom_line(position = pd, size = 0.5) +
  #points
  geom_point(position = pd)+
  #error bars
  geom_errorbar(aes(ymin = meanthreshold-se, ymax = meanthreshold+se), width = 0, size = 1, position = pd)+
  #minimal theme
  theme_minimal()+
  #add xlab
  xlab("Trial")+
  #add ylab
  ylab('Mean Threshold')+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"))+
  #add title
  ggtitle("Thresholds over Trials")

print(plotThresholdsoverTime)

plotFinalDirections <- ggplot(subset(allfinalHypo, trial == 35) , aes(x = factor(startAtSmall), fill = Structure, color = Structure)) +
  #points
  geom_bar(alpha = .4, position="dodge")+
  #minimal theme
  theme_minimal()+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #add xlab
  xlab("Start from the smallest Bar")+
  #add ylab
  ylab('Runs')+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"))+
  #add title
  ggtitle("Final Directions")

print(plotFinalDirections)


# labe the correct connection for the structure condition:
allfinalHypo$correctConnection = "other 3"
for (i in 1:length(allfinalHypo$connection)){
  if (allfinalHypo$connection[i] == allfinalHypo$trueConnection[i]){
    allfinalHypo$correctConnection[i] = "3correct"
  }else if(allfinalHypo$connection[i] == substr(allfinalHypo$trueConnection[i], 1, 2)){
    allfinalHypo$correctConnection[i] = "2correct"
  }else if(allfinalHypo$connection[i] == substr(allfinalHypo$trueConnection[i], 2, 3)){
    allfinalHypo$correctConnection[i] = "2correct"
  }else if(length(string_to_vec(allfinalHypo$connection[i])) <= 1){
    allfinalHypo$correctConnection[i] = "no"
  }else if(length(string_to_vec(allfinalHypo$connection[i])) == 2){
    allfinalHypo$correctConnection[i] = "2other"
  }
}

plotFinalConnections <- ggplot(subset(allfinalHypo, trial == 35) , aes(x = correctConnection, fill = Structure, color = Structure)) +
  #points
  geom_bar(alpha = .4, position="dodge")+
  #minimal theme
  theme_minimal()+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #add xlab
  xlab("Connections")+
  #add ylab
  ylab('Runs')+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"))+
  #add title
  ggtitle("Final Connections")
#facet_grid(rows = vars(startAtSmall))

print(plotFinalConnections)


rtNoB = ddply(allfinalHypo, ~NoB + Structure , summarize, meanTime = mean(time), se = se(time))
pd <- position_dodge(.1)
plotNoB <- ggplot(rtNoB, aes(x = NoB, y = meanTime, fill = Structure, color = Structure)) +
  geom_line(position = pd, size = 0.5) +
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
  theme(text = element_text(size = 20, family = "sans"), legend.position = "top")+
  #add title
  ggtitle("Model Time")

print(plotNoB)

# plot rt gain through Structure
rtNoBGain = subset(rtNoB, Structure == "Query")
rtNoBGain2 = subset(rtNoB, Structure == "Sequence")
rtNoBGain$RTgain = subset(rtNoB, Structure == "None")$meanTime - rtNoBGain$meanTime 
rtNoBGain2$RTgain = subset(rtNoB, Structure == "None")$meanTime - rtNoBGain2$meanTime 
rtNoBGain <- rbind(rtNoBGain, rtNoBGain2)
pd <- position_dodge(.1)
plotNoBRTGain <- ggplot(rtNoBGain, aes(x = NoB, y = RTgain, fill = Structure, color = Structure)) +
  geom_line(position = pd, size = 0.5) +
  #points
  geom_point(position = pd)+
  theme_minimal()+
  xlab("Sequence Length")+
  #add ylab
  ylab('Model Time Gain')+
  #change fill
  scale_fill_manual(values = cbbPalette[-1])+
  #change color
  scale_color_manual(values = cbbPalette[-1])+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"), legend.position = "top")+
  #add title
  ggtitle("Model Time Gain through Structure")

print(plotNoBRTGain)


rtTrialsNoB = ddply(allfinalHypo, ~NoB + trial + Structure , summarize, meanTime = mean(time), se = se(time))

plotrtTrialsNoB <- ggplot(rtTrialsNoB, aes(x = trial, y = meanTime, fill = factor(NoB), color = factor(NoB))) +
  geom_line(position = pd, size = 0.5) +
  #points
  geom_point(position = pd)+
  #error bars
  geom_errorbar(aes(ymin = meanTime-se, ymax = meanTime+se), width = 0, size = 1, position = pd)+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"), legend.position = "top")+
  #add title
  ggtitle("Model time over Trials")+
  facet_grid(rows = vars(Structure))

print(plotrtTrialsNoB)


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
  theme(text = element_text(size = 20, family = "sans"), legend.position = "top")+
  #add title
  ggtitle("Model Time over Trials")

print(plotrtTrials)


##############################################################
# show how well the real rts correspond to the models RTs
##########################################################
ggplot(ModelToData, aes(x = time, y = realRT))+
  geom_point()+
  ylim(-10,10)+
  geom_smooth()


###########
# raw RT
##############
# this part is necesarry, so I can add the correct columns to the model dataframe
relevantD <- subset(d, Condition == "Sort" & Stimulus_type == "bars")
relevantD1 <- subset(relevantD, Structure == "Query")
relevantD2 <- subset(relevantD, Structure == "Sequence")
relevantD3 <- subset(relevantD, Structure == "None")
relevantD <- rbind(relevantD1, relevantD2, relevantD3)

ModelToData$rawRt <- relevantD$rt
modelData <- (subset(ModelToData, realAccuracy == 1 & rawRt <= 10 ))

ggplot(modelData, aes(x = time, y = rawRt, color = Structure))+
  geom_point()+
  theme(legend.position = "None")+
  geom_smooth(data = modelData, aes(x = time, y = rawRt))+
  theme_minimal()+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  theme(text = element_text(size = 20, family = "sans"), legend.position = "none")+
  xlab("Model Time")+
  #add ylab
  ylab('Participant RT')+
  #change fonts
  facet_grid(cols = vars(Structure))

########################################################################
# plotting the model that tries to emulate the effects of participant RT
##########################################################################
setwd(here("Model", "hypoFittingBrmModels"))
EffectsModel <- brm(time ~ Structure + NoB + (Structure + NoB|participant),
              data = modelData,
              chains = 2,
              save_pars = save_pars(all = TRUE),
              control = list(max_treedepth = 15, adapt_delta = 0.99),
              file = "hypoGeneratorEmulatingRTresultsBucketSort")
summary(EffectsModel)
# results bucket sort 08.06.2022
#                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            -0.18      0.03    -0.23    -0.13 1.00     2739     1692
# StructureQuery       -0.28      0.02    -0.33    -0.24 1.00     2187     1611
# StructureSequence    -0.15      0.03    -0.20    -0.09 1.00     1363     1437
# NoB                   1.68      0.01     1.66     1.69 1.00     2362     1460

# Results for no Structure Learning for Bucket sort 08.06.22
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            -1.00      0.00    -1.00    -1.00 1.91        3       24
# StructureQuery        0.00      0.00     0.00     0.00 1.84        3       33
# StructureSequence     0.00      0.00     0.00     0.00 1.33        5       50
# NoB                   2.00      0.00     2.00     2.00 1.87        3       35

mcmc_plot(EffectsModel, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("Model Time")+
  theme_minimal()+
  theme(text = element_text(size = 15, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Query Structure", "Sequence Structure", "Sequence Length"))# 





