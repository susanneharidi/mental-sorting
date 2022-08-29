# Finding good hyperparameters for the
# hypothesis evaluator
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

source(here("Model", "hypoEvaluatorFunctions.R"))   


###########################################################
# get data and set parameters
###########################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)

particles_to_eliminate = c(0.01, 0.015, 0.02,  0.025, 0.03, 0.035, 0.04, 0.05, 0.075, 0.1, 0.18, 0.25 ) # 
nofColors = 10    # 
RTCutoff = 10


# this part is necesarry, so I can add the correct columns to the model dataframe
relevantD <- subset(d, Condition == "Sort" & Stimulus_type == "bars")
relevantD1 <- subset(relevantD, Structure == "Query")
relevantD2 <- subset(relevantD, Structure == "Sequence")
relevantD3 <- subset(relevantD, Structure == "None")
relevantD <- rbind(relevantD1, relevantD2, relevantD3)
##################################################################################
# set up hypothesis space
###########################################################################

#all hypotheses of size 2
two <- t(combn(letters[1:nofColors], 2)) # it is important that the nofColors here is set to 10 as in the experiment

#all hypotheses of size 3
three <- t(combn(letters[1:nofColors], 3))

# all Start at small variations
startAtSmall <- c(0, 1)

#create all possible hypotheses
connections <- c('a', paste0(two[,1], two[,2]), paste0(three[,1], three[,2], three[,3]))

#generate particle hypotheses
particles <- expand.grid(threshold = 1:7, connections = connections, startAtSmall = startAtSmall)


################################################################
# run the model on all Data for different hyperparameters
##############################################################
setwd(here("Model", "hypoFittingBrmModels"))

# if the models have already been calculated set the following parameter to False
ModelNOTCalculatedYet == FALSE

loo_R2 <- c(1:length(particles_to_eliminate))
for (i in 1:length(particles_to_eliminate)){
  if (ModelNOTCalculatedYet == TRUE){
    output <- runHypoSortOnAllData(start_particles = particles,
                                    AllData = d,
                                    particles_to_eliminate[i])
    ModelToData <- output[[1]]
    allfinalHypo <- output[[2]]
    # add the raw RT to the model output (the real RT in the model is the RT difference between memory and sort, which is the wrong metric here)
    ModelToData$rawRt <- relevantD$rt
    # only chose the trials, in which participants met the inclusion criteria
    modelData <- (subset(ModelToData, realAccuracy == 1 & rawRt <= RTCutoff ))
  }
  
  filename = paste("hypoEvaluatorRawRTPtoEBucketSort", as.character(particles_to_eliminate[i]), sep = "_")
  model <- brm(rawRt ~ time + (time|participant), 
               data = modelData,
               chains = 2,
               save_pars = save_pars(all = TRUE),
               file = filename)
  goodness <- loo_R2(model)
  loo_R2[i] <- goodness[1]
  # this just serves as a signal that a model is finished
  plot(i, goodness[1])
}

###############################################################################
# finding good hyperparameters by inspecting the final LooR2
################################################################################

# set up a dataframe for the looR2 values and the corresponding hyperparameters
Hyperparameters <- data.frame(particles_to_eliminate = particles_to_eliminate,
                            loo_R2 = loo_R2)

# plot the LooR2 for the hyperparaemeters
ggplot(Hyperparameters, aes(x = particles_to_eliminate, y = loo_R2))+
  geom_point()+
  ylim(0.54, 0.7)+
  geom_line()


