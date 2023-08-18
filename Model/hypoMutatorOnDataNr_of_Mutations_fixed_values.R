# run hypothesis mutator version with one hypothesis (Appendix B) on the real trials
# the output from this script is used fro the model comparion for Appendix B



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

##########################################################################
#load the data and set some parameters
##############################################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)

RTCutoff <- 10
seed <- 2022

##############################################################
# set hyperparameter I want:
################################################################

HyperPs <-  c(1, 10, 100, 1000, 10000)

###########################################################
# run this on one participant at a time
############################################################

subjects <- unique(d$subject_id)

# load the data
ModelToData <- read.csv(here("Model", "FittedModelData", "hypoMutatorfitted_BucketSort_Nr_of_mutations_fixed_values.csv"))

# or calculate it anew with the code below

ModelToData <- data.frame(NoB = numeric(),
                          time = numeric(),
                          realRT = numeric(),
                          accuracy = numeric(),
                          realAccuracy = numeric(),
                          trueConnection = character(),
                          participant = character(),
                          Structure = character(),
                          trial = numeric(),
                          Nr_of_mutations = numeric())
allfinalHypo <- data.frame(threshold = numeric(),
                           connections = character(),
                           startAtSmall = numeric(),
                           time = numeric(),
                           correct = numeric(),
                           run = numeric(),
                           Structure = character(),
                           trueConnection = character())

for (nNr_of_mutations in c(1000, 10000)){
  print("Number of Mutations")
  print(nNr_of_mutations)
  for (currentsubject in subjects){
    if (currentsubject == "bvoobt12gur23nj"){
      continue = TRUE
    }
    print(currentsubject)
    
    if (continue){
      
      nStartParticles <- 1
      currtentData <- subset(d, subject_id == currentsubject)
      
      
      particles = generateRandomParticles(nOfParticles = nStartParticles, 
                                          nOfBars = 7, 
                                          nOfMaxConnections = 3, 
                                          threshold = TRUE, 
                                          connectedness = TRUE,
                                          startAtSmall = TRUE)
      
      
      output <- runHypoMutateOnAllData(particles,
                                       currtentData,
                                       nNr_of_mutations)
      
      output[[1]]$Nr_of_mutations = nNr_of_mutations
      ModelToData <- rbind(ModelToData, output[[1]])
      allfinalHypo <- rbind(allfinalHypo, output[[2]])
    }
  }
}
#################################
# save model to data dataframe for further analysis
###################################

#write_csv(ModelToData, here("Model", "FittedModelData", "hypoMutatorfitted_BucketSortReview_Nr_of_mutations_fixed_values.csv"))
ModelToData <- read.csv(here("Model", "FittedModelData", "hypoMutatorfitted_BucketSort_Nr_of_mutations_fixed_values.csv"))




