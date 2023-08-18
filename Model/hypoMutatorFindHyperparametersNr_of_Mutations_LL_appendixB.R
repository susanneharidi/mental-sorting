# Finding good hyperparameters for the
# hypothesis mutator chnaged version with only one hypothesis Appendix B

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

source(here("Model", "hypoMutatorFunctionsReview_several_Mutations_per_Trial.R"))

set.seed(0)
#################################################################
# Functions
##################################################################
ModelLLrealRT <- function(AllData, Nr_of_mutations, particle){
  output <- runHypoMutateOnAllData(particle,
                                   AllData,
                                   Nr_of_mutations)
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

nStartParticles = 1
possible_Nr_of_mutations = c(1, 5, 10, 100, 1000, 10000) 


##############################################################
# run hypo mutator on all data
################################################################

# load dataframe of LogLokelihoods 
LLdf <- read.csv(here("Model", "FittedModelData", "hypoMutator_hyperparameterfitting_LL_one_hypo_several_muts.csv"))

# or caculate it anew with the code below

subjects <- unique(d$subject_id)

LLdf <- data.frame(subject  = character(),
                    nStartParticles = numeric(),
                    Nr_of_mutations = numeric(),
                    LL = numeric())


for (subject in subjects){
  
  print("subject:")
  print(subject)
  

  particle = generateRandomParticles(nOfParticles = nStartParticles, 
                                      nOfBars = 7, 
                                      nOfMaxConnections = 3, 
                                      threshold = TRUE, 
                                      connectedness = TRUE,
                                      startAtSmall = TRUE)
  
  LL <- c(1:length(possible_Nr_of_mutations))
  for (l in 1:length(possible_Nr_of_mutations)){
      allData <- subset(d, subject_id == subject)
      LL[l] <- ModelLLrealRT(allData, possible_Nr_of_mutations[l], particle)
      print("current NR of Mutations")
      print(possible_Nr_of_mutations[l])
  }
  LLdfTemp <- data.frame(subject  = subject,
                          nStartParticles = nStartParticles,
                          Nr_of_mutations = possible_Nr_of_mutations,
                          LL = LL)

  LLdf <- rbind(LLdf, LLdfTemp)

  # plot partricipantwise results
  
  ggplot(LLdf, aes(x = Nr_of_mutations, y = LL, color = subject))+
    geom_point()+
    geom_line()+
    ggtitle(subject)+
    theme(legend.position = "top")
}


ggplot(LLdf, aes(x = log(Nr_of_mutations), y = LL, color = subject))+
  geom_point()+
  geom_line()+
  theme_minimal()+
  theme(legend.position = "none")

# write_csv(LLdf, here("Model", "FittedModelData", "hypoMutator_hyperparameterfitting_LL_one_hypo_several_muts.csv"))

#########################################################
# identify the hyperparameters with the minimal LL
#-> I minimize here, because my LL function returns the negative LL
######################################################
bestHyperPs <- data.frame(Nr_of_mutations = numeric(),
                          subject = character(),
                          minLL = numeric())

for(currentsubject in subjects){
  data <- subset(LLdf, subject == currentsubject)
  minIndex <- which.min(data$LL)
  nNr_of_mutations <- data$Nr_of_mutations[minIndex]
  minLL <- data$LL[minIndex]
  bestHyperP <- data.frame(Nr_of_mutations = nNr_of_mutations,
                            subject = currentsubject,
                            minLL = minLL)
  bestHyperPs <- rbind(bestHyperPs, bestHyperP)
}


ggplot(bestHyperPs, aes(x = log(Nr_of_mutations))) + 
  geom_histogram()+
  theme_minimal()


write_csv(bestHyperPs, here("Model", "FittedModelData", "hypoMutatorBestHyperPs_one_hypo_several_muts.csv"))

