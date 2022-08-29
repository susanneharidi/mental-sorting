# run hypothesis mutator on the real trials
# the output of this script is used for the model comparison


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
library(readr)

source(here("Model", "hypoMutatorFunctions.R")) 

##########################################################################
#load the data and set some parameters
##############################################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)



nStartParticles = 15
particles_to_mutate = 1

particles = generateRandomParticles(nOfParticles = nStartParticles, 
                                    nOfBars = 7, 
                                    nOfMaxConnections = 3, 
                                    threshold = TRUE, 
                                    connectedness = TRUE,
                                    startAtSmall = TRUE)

##############################################################
# run hypo mutator on all data
################################################################

output <- runHypoMutateOnAllData(particles,
                                 d,
                                 particles_to_mutate)

ModelToData <- output[[1]]
allfinalHypo <- output[[2]]

#################################
# save model to data dataframe for further analysis
###################################
write_csv(ModelToData, here("Model", "FittedModelData", paste(as.character(nStartParticles), as.character(particles_to_mutate), "hypoMutatorfitted_BucketSort.csv", sep = "_")))
# load the data if it has already been saved before (this will not allow all plotting below)
ModelToData <- read.csv(here("Model", "FittedModelData",paste(as.character(nStartParticles), as.character(particles_to_mutate), "hypoMutatorfitted_BucketSort.csv", sep = "_")))

############################################
# plot the output of the model
############################################
# plot all thresholds
cbbPalette <- c("grey50","#0072B2", "#D55E00")
pd <- position_dodge(.2)


plotFinalThresholds <- ggplot(allfinalHypo , aes(x = threshold, fill = Structure, color = Structure)) +
  #points
  geom_histogram(alpha = .4, position="dodge")+
  #minimal theme
  theme_minimal()+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #add xlab
  xlab("threshold")+
  #add ylab
  ylab('count')+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"))+
  #add title
  ggtitle("Thresholds")

print(plotFinalThresholds)

plotFinalDirections <- ggplot(allfinalHypo , aes(x = factor(startAtSmall), fill = Structure, color = Structure)) +
  #points
  geom_bar(alpha = .4, position="dodge")+
  #minimal theme
  theme_minimal()+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #add xlab
  xlab("direction")+
  #add ylab
  ylab('count')+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"))+
  #add title
  ggtitle("Directions")

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
  }else if(length(string_to_vec(as.character(allfinalHypo$connection[i]))) <= 1){
    allfinalHypo$correctConnection[i] = "no"
  }else if(length(string_to_vec(as.character(allfinalHypo$connection[i]))) == 2){
    allfinalHypo$correctConnection[i] = "2other"
  }
}

plotFinalConnections <- ggplot(allfinalHypo , aes(x = correctConnection, fill = Structure, color = Structure)) +
  #points
  geom_bar(alpha = .4, position="dodge")+
  #minimal theme
  theme_minimal()+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #add xlab
  xlab("connections")+
  #add ylab
  ylab('count')+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"))+
  #add title
  ggtitle("final connections")
#facet_grid(rows = vars(startAtSmall))

print(plotFinalConnections)


rtNoB = ddply(ModelToData , ~NoB + Structure , summarize, meanTime = mean(time), se = se(time))
pd <- position_dodge(.2)
plotNoB <- ggplot(rtNoB, aes(x = NoB, y = meanTime, fill = Structure, color = Structure)) +
  geom_line(position = pd, size = 0.5) +
  #points
  geom_point(position = pd)+
  #minimal theme
  theme_minimal()+
  #error bars
  geom_errorbar(aes(ymin = meanTime-se, ymax = meanTime+se), width = 0, size = 1, position = pd)+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"), legend.position = "top")+
  #add title
  ggtitle("Model time per NoB")

print(plotNoB)


rtTrialsNoB = ddply(ModelToData , ~NoB + trial + Structure , summarize, meanTime = mean(time), se = se(time))

plotrtTrialsNoB <- ggplot(rtTrialsNoB, aes(x = trial, y = meanTime, fill = factor(NoB), color = factor(NoB))) +
  geom_line(position = pd, size = 0.5) +
  #points
  geom_point(position = pd)+
  #minimal theme
  theme_minimal()+
  #error bars
  geom_errorbar(aes(ymin = meanTime-se, ymax = meanTime+se), width = 0, size = 1, position = pd)+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"), legend.position = "top")+
  #add title
  ggtitle("Model time over Trials")+
  facet_grid(rows = vars(Structure))

print(plotrtTrialsNoB)


rtTrials = ddply(ModelToData , ~trial + Structure , summarize, meanTime = mean(time), se = se(time))

plotrtTrials <- ggplot(rtTrials, aes(x = trial, y = meanTime, fill = Structure, color = Structure)) +
  geom_line(position = pd, size = 0.5) +
  #points
  #minimal theme
  theme_minimal()+
  geom_point(position = pd)+
  #error bars
  geom_errorbar(aes(ymin = meanTime-se, ymax = meanTime+se), width = 0, size = 1, position = pd)+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"), legend.position = "top")+
  #add title
  ggtitle("Model time over trials")

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
  geom_smooth(data = modelData, aes(x = time, y = rawRt), method = lm)+
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
                    file = "hypoMutatorEmulatingRTresultsBucketSort_15_1")
summary(EffectsModel)
# results 15.06.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.13      0.04     0.05     0.21 1.00     1176     1509
# StructureQuery       -1.14      0.05    -1.26    -1.04 1.00      896     1407
# StructureSequence     0.01      0.05    -0.09     0.09 1.00     1156     1601
# NoB                   1.59      0.01     1.57     1.62 1.00      829     1254

mcmc_plot(EffectsModel, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("Model Time")+
  theme_minimal()+
  theme(text = element_text(size = 15, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Query Structure", "Sequence Structure", "Sequence Length"))
