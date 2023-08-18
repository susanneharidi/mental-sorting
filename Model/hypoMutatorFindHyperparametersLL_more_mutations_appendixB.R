# Finding good hyperparameters for the
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
library(readr)
library(patchwork)  #https://github.com/thomasp85/patchwork
library(scales)
library(plot.matrix)
library(RColorBrewer)
library(corrplot)
library(png)

source(here("Model", "hypoMutatorFunctions.R"))

set.seed(0)
#################################################################
# Functions
##################################################################
ModelLLrealRT <- function(AllData, particles_to_mutate, particles){
  output <- runHypoMutateOnAllData(particles,
                                   AllData,
                                   particles_to_mutate)
  #####################################
  # exlude all data where the participant was incorrect or th real RT is larger than 10/-10
  ############################################
  ModelToData <- output[[1]]
  ModelToData <- subset(ModelToData, realAccuracy == 1 & realRT <= RTCutoff)
  output <- lm(realRT ~ time, data = ModelToData)
  intercept <- as.numeric(output$coefficients[1])
  beta <- as.numeric(output$coefficients[2])
  ModelToData$timeasRT = intercept + ModelToData$time*beta
  #calculate the loglikelihood
  sd <- sd(ModelToData$realRT, na.rm = FALSE)
  LL <- 0
  for (i in 1:length(ModelToData$timeasRT)){
    LL <- LL + dnorm(x = ModelToData$realRT[i], mean = ModelToData$timeasRT[i], sd = sd, log = TRUE)
  }

  return(-LL)
}


##########################################################################
# load data and set some parameters
##############################################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)

RTCutoff = 10

# extended hyper parameters with far larger number of start Particles for model in appendix
nStartParticles = c(5, 10, 15, 20, 100, 1000, 10000)
particles_to_mutate = c(1, 2, 3, 4, 5, 10, 50, 500, 750, 1000, 5000, 7500) #


##############################################################
# run hypo mutator on all data
################################################################


# either load the dataframe with the LogLikelihood values
LLdf <- read.csv(here("Model", "FittedModelData", "hypoMutator_hyperparameterfitting_LL_with_more_mutations.csv"))

# or calculate the dataframe anew with the code below
LLdf <- data.frame(subject  = character(),
                    nStartParticles = numeric(),
                    particles_to_mutate = numeric(),
                    LL = numeric())


for (subject in subjects){
  
  print("subject:")
  print(subject)
  
  for (i in 1:length(nStartParticles)){
    particles = generateRandomParticles(nOfParticles = nStartParticles[i], 
                                        nOfBars = 7, 
                                        nOfMaxConnections = 3, 
                                        threshold = TRUE, 
                                        connectedness = TRUE,
                                        startAtSmall = TRUE)
    
    current_particles_to_mutate <- particles_to_mutate[particles_to_mutate < nStartParticles[i]]
    current_particles_to_mutate <- current_particles_to_mutate[current_particles_to_mutate > nStartParticles[i]/21]
    LL <- c(1:length(current_particles_to_mutate))
    for (l in 1:length(current_particles_to_mutate)){
        allData <- subset(d, subject_id == subject)
        LL[l] <- ModelLLrealRT(allData, min(current_particles_to_mutate[l], nStartParticles[i]), particles)
        print(paste(LL[l], current_particles_to_mutate[l], nStartParticles[i]))
    }
    LLdfTemp <- data.frame(subject  = subject,
                            nStartParticles = nStartParticles[i],
                            particles_to_mutate = current_particles_to_mutate,
                            LL = LL)

    LLdf <- rbind(LLdf, LLdfTemp)
  }
  # plot partricipantwise results
  ggplot(subset(LLdf, subject == subject), aes(x = particles_to_mutate, y = LL, color = nStartParticles))+
    geom_point()+
    geom_line()+
    ggtitle(subject)+
    theme(legend.position = "top")
}


ggplot(LLdf, aes(x = particles_to_mutate, y = LL, color = subject))+
  geom_point()+
  geom_line()+
  facet_grid(rows = vars(nStartParticles))+
  theme(legend.position = "none")

#write_csv(LLdf, here("Model", "FittedModelData", "hypoMutator_hyperparameterfitting_LL_with_more_mutations.csv"))

#########################################################
# identify the hyperparameters with the minimal LL
#-> I minimize here, because my LL function returns the negative LL
######################################################
bestHyperPs <- data.frame(nStartParticles = numeric(),
                          particles_to_mutate = numeric(),
                          subject = character(),
                          minLL = numeric())

for(currentsubject in subjects){
  data <- subset(LLdf, subject == currentsubject)
  minIndex <- which.min(data$LL)
  nStartParticles <- data$nStartParticles[minIndex]
  particles_to_mutate <- data$particles_to_mutate[minIndex]
  minLL <- data$LL[minIndex]
  bestHyperP <- data.frame(nStartParticles = nStartParticles,
                            particles_to_mutate = particles_to_mutate,
                            subject = currentsubject,
                            minLL = minLL)
  bestHyperPs <- rbind(bestHyperPs, bestHyperP)
}


ggplot(bestHyperPs, aes(x = log(nStartParticles), fill = as.factor(particles_to_mutate))) + 
  geom_bar(position = "dodge2")+
  theme(legend.position = "top")


write_csv(bestHyperPs, here("Model", "FittedModelData", "hypoMutatorBestHyperPsReview_with_more_mutations.csv")) 
