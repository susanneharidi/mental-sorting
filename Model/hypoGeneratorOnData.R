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
library(stringr)
library(readr)
library(patchwork)  #https://github.com/thomasp85/patchwork
library(scales)
library(plot.matrix)
library(RColorBrewer)
library(corrplot)
library(png)

source(here("Model", "hypoGeneratorFunctions.R")) 

###############################################################################################################
###############################################################################################################
# run hypo learner on data
####################################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)



##################################################################################
# set up hypothesis space
###########################################################################

start_particle = generateStartParticles(nOfBars = 7, nOfColors = 10)

###############################
# run sorter on all data with the best LRs
###############################
bestLRs <- read.csv(here("Model", "FittedModelData", "hypoGeneratorBestHyperPsIndLL_afterBugFix05_23.csv"))


source(here("Model", "hypoGeneratorFunctionsIndLR.R")) 
ModelToData <- runHypoLearnerOnAllDataIndLR(start_particle = generateStartParticles(nOfBars = 7, nOfColors = 10),
                                            d,
                                            bestLRs)

## add the LR to the dataframe
participants = unique(ModelToData$participant)
for (i in 1:length(participants)){
  currentData = subset(ModelToData, participant == participants[i])
  currentData$LR = subset(bestLRs, subject == participants[i])$LR
  if (i == 1){
    ModelToDataWithLR = currentData
  }else{
    ModelToDataWithLR = rbind(ModelToDataWithLR, currentData)
  }
}

#################################
# save model to data dataframe for further analysis
###################################

write_csv(ModelToData, here("Model", "FittedModelData", "hypoGeneratorfitted_BucketSortIndHypoPsLL_afterBugFix05_23.csv"))
                    # orignal data for first submission "hypoGeneratorfitted_BucketSortIndHypoPsLL0_1.csv
ModelToData <- read.csv(here("Model", "FittedModelData", "hypoGeneratorfitted_BucketSortIndHypoPsLL_afterBugFix05_23.csv"))

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
RTCutoff <- 10
textsize <- 16

modelData <- (subset(ModelToData, realAccuracy == 1 & realRT <= RTCutoff))

ModelvsDataPlot <- ggplot(modelData, aes(x = time, y = realRT, color = Structure))+
  geom_point(alpha = 0.2)+
  theme(legend.position = "None")+
  geom_smooth(data = modelData, aes(x = time, y = realRT), method = lm)+
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
#############################################################################

setwd(here("Model", "hypoFittingBrmModels"))
EffectsModel <- brm(time ~ Structure + NoB + (Structure + NoB|participant),
                    data = modelData,
                    chains = 2,
                    save_pars = save_pars(all = TRUE),
                    control = list(max_treedepth = 15, adapt_delta = 0.99),
                    file = "hypoGeneratorEmulatingRTresultsBucketSortIndLRLL_afterBugFix05_23")
                    # old: hypoGeneratorEmulatingRTresultsBucketSortIndLRLL0_1
summary(EffectsModel)

# results 14.10.2022 (LL)
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.15      0.03     0.08     0.21 1.00     2084     1551
# StructureQuery       -0.51      0.03    -0.56    -0.46 1.00     2861     1365
# StructureSequence    -0.27      0.05    -0.37    -0.17 1.00     1129     1315
# NoB                   1.52      0.01     1.50     1.54 1.00     1136     1380

# results 11.05.2023 (LL) after bug fix
#                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            -0.12      0.06    -0.24    -0.00 1.00      172      592
# StructureQuery       -0.35      0.03    -0.41    -0.28 1.00      394      918
# StructureSequence    -0.32      0.06    -0.44    -0.20 1.00      273      658
# NoB                   1.65      0.02     1.60     1.69 1.00      155      349



ModelofModelPlot <- mcmc_plot(EffectsModel, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("Posterior estimates")+
  theme_minimal()+
  theme(text = element_text(size = 18, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Query Structure", "Sequence Structure", "Sequence Length"))# 

ModelofModelPlot

ModelOutputPlot <- (plotNoB + ModelofModelPlot + plotrtTrials)/(plotFinalThresholds + plotFinalConnections + ModelvsDataPlot) +
  plot_annotation(tag_levels = "A")+ 
  plot_layout(guides = 'collect')
ModelOutputPlot


