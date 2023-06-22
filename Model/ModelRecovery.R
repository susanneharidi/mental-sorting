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
for (j in 1:20){

  start_particle = generateStartParticles(nOfBars = 7, nOfColors = 10)
  
  ###############################
  # run sorter on all data with the best LRs
  ###############################
  bestLRs <- read.csv(here("Model", "FittedModelData", "hypoGeneratorBestHyperPsIndLL_afterBugFix05_23.csv"))
  
  source(here("Model", "hypoGeneratorFunctionsIndLR.R")) 
  
  
  
  for (i in 1:2){
    ModelToData <- runHypoLearnerOnAllDataIndLR(start_particle = generateStartParticles(nOfBars = 7, nOfColors = 10),
                                                d,
                                                bestLRs)
    modelData = data.frame(time = ModelToData$time,
                           NoB = ModelToData$NoB,
                           model = "generator",
                           FullRun = i,
                           trial = ModelToData$trial,
                           structure = ModelToData$Structure,
                           participant = ModelToData$participant)
    if (i == 1){
      allModelData = modelData
    }else{
      allModelData = rbind(allModelData, modelData)
    }
      
  }
  
  
  ####################################################################
  # hypothesis mutator
  #####################################################################
  
  source(here("Model", "hypoMutatorFunctions.R")) 
  
  ##########################################################################
  #load the data and set some parameters
  ##############################################################################
  d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
  d$subject_id <- factor(d$subject_id)
  
  
  ##############################################################
  # run hypo mutator on all data
  ################################################################
  
  bestHyperPs <- read.csv(here("Model", "FittedModelData", "hypoMutatorBestHyperPsIndLL_afterBugFix05_23.csv"))
  
  ###########################################################
  # run this on one participant at a time
  ############################################################
  
  subjects <- unique(d$subject_id)
  for (i in 1:2){
  
    
    ModelToData <- data.frame(NoB = numeric(),
                              time = numeric(),
                              realRT = numeric(),
                              accuracy = numeric(),
                              realAccuracy = numeric(),
                              trueConnection = character(),
                              participant = character(),
                              Structure = character(),
                              trial = numeric())
    allfinalHypo <- data.frame(threshold = numeric(),
                               connections = character(),
                               startAtSmall = numeric(),
                               time = numeric(),
                               correct = numeric(),
                               run = numeric(),
                               Structure = character(),
                               trueConnection = character())
    
    for (currentsubject in subjects){
      currentHyperPs <- subset(bestHyperPs, subject == currentsubject)
      nStartParticles <- currentHyperPs$nStartParticles
      particles_to_mutate = currentHyperPs$particles_to_mutate
      currtentData <- subset(d, subject_id == currentsubject)
      
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
                                       currtentData,
                                       particles_to_mutate)
      
      ModelToData <- rbind(ModelToData, output[[1]])
      allfinalHypo <- rbind(allfinalHypo, output[[2]])
    }
    
    modelData = data.frame(time = ModelToData$time,
                           NoB = ModelToData$NoB,
                           model = "mutator",
                           FullRun = i,
                           trial = ModelToData$trial,
                           structure = ModelToData$Structure,
                           participant = ModelToData$participant)
    allModelData = rbind(allModelData, modelData)
  
  }
  
  # calculate the correlations
  CorrMutatorMutator = cor(subset(allModelData, model == "mutator" & FullRun == 1)$time,
                           subset(allModelData, model == "mutator" & FullRun == 2)$time,
                           method = "pearson")
  
  CorrMutatorGenerator = cor(subset(allModelData, model == "mutator" & FullRun == 1)$time,
                           subset(allModelData, model == "generator" & FullRun == 1)$time,
                           method = "pearson")
  
  CorrGeneratorGenerator = cor(subset(allModelData, model == "generator" & FullRun == 1)$time,
                           subset(allModelData, model == "generator" & FullRun == 2)$time,
                           method = "pearson")
  
  CorrGeneratorMutator = cor(subset(allModelData, model == "generator" & FullRun == 2)$time,
                           subset(allModelData, model == "mutator" & FullRun == 2)$time,
                           method = "pearson")
  
  # Plot the matrix
  library('plot.matrix')
  # numeric matrix
  x <- matrix(1:4, nrow = 2, ncol = 2) # create a numeric matrix object
  x[1,1] = CorrMutatorMutator
  x[1,2] = CorrMutatorGenerator
  x[2,1] = CorrGeneratorMutator
  x[2,2] = CorrGeneratorGenerator
  class(x)
  plot(x)
  
  Correlations = data.frame(MutatorMutator = CorrMutatorMutator,
                            MutatorGenerator = CorrMutatorGenerator,
                            GeneratorMutator = CorrGeneratorMutator,
                            GeneratorGenerator = CorrGeneratorGenerator)
  if (j == 1){
    AllCorrelations = Correlations
  }else{
    AllCorrelations = rbind(AllCorrelations, Correlations)
  } 
  
}

# reorder data for easier plotting

AllCorrelationsdf1 = data.frame(
  correlation = AllCorrelations$MutatorMutator,
  model1 = "mutator",
  model2 = "mutator"
)

AllCorrelationsdf2 = data.frame(
  correlation = AllCorrelations$MutatorGenerator,
  model1 = "mutator",
  model2 = "generator"
)

AllCorrelationsdf3 = data.frame(
  correlation = AllCorrelations$GeneratorMutator,
  model1 = "generator",
  model2 = "mutator"
)


AllCorrelationsdf4 = data.frame(
  correlation = AllCorrelations$GeneratorGenerator,
  model1 = "generator",
  model2 = "generator"
)

AllCorrelationsdf = rbind(AllCorrelationsdf1, AllCorrelationsdf2,
                          AllCorrelationsdf3, AllCorrelationsdf4)

# save the dataframe
write.csv(AllCorrelationsdf, here("Data", "ModelRecoveryCorrelations_afterBugFix05_23.csv"), row.names=FALSE, quote=FALSE) 
# and load it
AllCorrelationsdf = read.csv(here("Data", "ModelRecoveryCorrelations_afterBugFix05_23.csv"))

textsize = 23 
modelRecoveryPlot <- ggplot(AllCorrelationsdf, aes(x=model1, y=model2, fill=correlation))+
      geom_tile()+
      xlab("Model 2")+
      #add ylab
      ylab('Model 1')+
      #change fonts
      #minimal theme
      theme_minimal()+
      theme(text = element_text(size = textsize, family = "sans"))
      
modelRecoveryPlot

meanCorrelations = ddply(AllCorrelationsdf, ~ model1 + model2, summarize, meanCorrelation = mean(correlation))
meanCorrelations
# results 10.05.2023
# model1    model2 meanCorrelation
# generator generator      1.00000000
# generator   mutator      0.04610458
#   mutator generator      0.04746784
#   mutator   mutator      0.93347756

# results 12.05.2023 after bigfix
#    model1    model2 meanCorrelation
# 1 generator generator      1.00000000
# 2 generator   mutator      0.03821748
# 3   mutator generator      0.04101716
# 4   mutator   mutator      0.92696374


ggsave(here("Figures", "ModelRecovery_Reviewed.pdf"), modelRecoveryPlot,  
       width = 10, 
       height = 7)
