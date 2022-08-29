# Finding good hyperparameters (learningrate) for the
# hypothesis generator
# (via LooR^2 values)

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
library(lme4)
library(stats)

source(here("Model", "hypoGeneratorFunctions.R")) 

###############################################################################################################
###############################################################################################################
# load data and set some parameters
####################################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)
RTCutoff = 10
learningrate = seq(from = 0.01, to = 0.1, by = 0.01)

# this part is necessary, so I can add the correct columns to the model dataframe
relevantD <- subset(d, Condition == "Sort" & Stimulus_type == "bars")
relevantD1 <- subset(relevantD, Structure == "Query")
relevantD2 <- subset(relevantD, Structure == "Sequence")
relevantD3 <- subset(relevantD, Structure == "None")
relevantD <- rbind(relevantD1, relevantD2, relevantD3)

##################################################################################
# set up hypothesis space
###########################################################################

start_particle = generateStartParticles(nOfBars = 7, nOfColors = 10)

###########################
# find the best hyperparameter (e.g. the learningrate)
###############################
setwd(here("Model", "hypoFittingBrmModels"))

# if the models have already been calculated set the following parameter to False
ModelNOTCalculatedYet == FALSE

loo_R2 <- c(1:length(learningrate))
for (i in 1: length(learningrate)){
  if (ModelNOTCalculatedYet == TRUE){
    ModelToData <- runHypoLearnerOnAllData(start_particle,
                                           AllData = d,
                                           learningrate[i])
    # add the raw RT to the model output (the real RT in the model is the RT difference between memory and sort, which is the wrong metric here)
    ModelToData$rawRt <- relevantD$rt
    modelData <- (subset(ModelToData, realAccuracy == 1 & rawRt <= RTCutoff ))
  }
  filename = paste("hypeGeneratorMaxTHProbRawRTLRBucketSort", as.character(learningrate[i]), sep = "_")
  model <- brm(rawRt ~ time + (time|participant), 
               data = modelData,
               chains = 2,
               save_pars = save_pars(all = TRUE),
               file = filename)
  goodness <- loo_R2(model)
  loo_R2[i] <- goodness
  
  #this just serves as a signal that a model is finished
  plot(i,goodness[1])
} 

###############################################################################
# finding good hyperparameters by inspecting the final LooR2
###############################################################################

# set up a dataframe for the looR2 values and the corresponding hyperparameters
Hyperparameters <- data.frame(learningrate = learningrate,
                              loo_R2 = loo_R2)

# plot the LooR2 for the hyperparameters
ggplot(Hyperparameters, aes(x = learningrate, y = loo_R2))+
  geom_point()+
  ylim(0.54, 0.7)+
  geom_line()



