# Finding good hyperparameters (learningrate) for the
# hypothesis mutator
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

source(here("Model", "hypoMutatorFunctions.R")) 

##########################################################################
# load data and set some parameters
##############################################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)

RTCutoff = 10

nStartParticles = c(5, 10, 15, 20) # 5, 10,
particles_to_mutate = c(1, 2, 3, 4, 5) # 



# this part is necesarry, so I can add the correct columns to the model dataframe
relevantD <- subset(d, Condition == "Sort" & Stimulus_type == "bars")
relevantD1 <- subset(relevantD, Structure == "Query")
relevantD2 <- subset(relevantD, Structure == "Sequence")
relevantD3 <- subset(relevantD, Structure == "None")
relevantD <- rbind(relevantD1, relevantD2, relevantD3)

##############################################################
# run hypo mutator on all data
################################################################
setwd(here("Model", "hypoFittingBrmModels"))

# if the models have already been calculated set the following parameter to False
ModelNOTCalculatedYet == FALSE

theLength <- (length(particles_to_mutate)*length(nStartParticles))
loo_R2 <- c(1:theLength)
nParticles <- c(1:theLength)
mutatingParticles <- c(1:theLength)
for (i in 1:length(nStartParticles)){
  particles = generateRandomParticles(nOfParticles = nStartParticles[i], 
                                      nOfBars = 7, 
                                      nOfMaxConnections = 3, 
                                      threshold = TRUE, 
                                      connectedness = TRUE,
                                      startAtSmall = TRUE)
  for (l in 1:length(particles_to_mutate)){
    if (ModelNOTCalculatedYet == TRUE){
      output <- runHypoMutateOnAllData(particles,
                                       d,
                                       particles_to_mutate[l])
  
      ModelToData <- output[[1]]
      allfinalHypo <- output[[2]]
      # add the raw RT to the model output (the real RT in the model is the RT difference between memory and sort, which is the wrong metric here)
      ModelToData$rawRt <- relevantD$rt
      modelData <- (subset(ModelToData, realAccuracy == 1 & rawRt <= RTCutoff ))
    }
    filename = paste("hypoMutatorRawRTBucketSortNoP_MP", as.character(nStartParticles[i]),
                     as.character(particles_to_mutate[l]), sep = "_")
    model <- brm(rawRt ~ time + Structure + (time + Structure|participant), 
                 data = modelData,
                 chains = 2,
                 save_pars = save_pars(all = TRUE),
                 file = filename)
    goodness <- loo_R2(model)
    position = (i-1)*length(particles_to_mutate)+l
    loo_R2[position] <- goodness
    nParticles[position] <- nStartParticles[i]
    mutatingParticles[position] <- particles_to_mutate[l]
    plot(nStartParticles[i],particles_to_mutate[l])
  }
}
###############################################################################
# finding good hyperparameters by inspecting the final LooR2
###############################################################################

# set up a dataframe for the looR2 values and the corresponding hyperparameters
Hyperparameters <- data.frame(mutatingParticles = mutatingParticles,
                              nParticles = nParticles,
                              loo_R2 = loo_R2)

# plot the LooR2 for the hyperparameters
ggplot(subset(Hyperparameters, mutatingParticles < 6 & loo_R2 < 1), aes(x = factor(nParticles), y = loo_R2, group = factor(mutatingParticles)))+
  geom_point(aes(color = factor(mutatingParticles)))+
  ylim(0.54, 0.7)+ 
  geom_line(aes(color = factor(mutatingParticles)))

  