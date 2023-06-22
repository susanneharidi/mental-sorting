#
#
#checking different leanringrates for the hypo generator model
#



################################################################################
# libraries
###############################################################################
library(ggplot2)
library(cowplot)
packages <- c('plyr', 'ggplot2', 'ggthemes', 'ggsignif', 'gridExtra')
#load
lapply(packages, library, character.only = TRUE)
library(ggpubr)
library(here)
library(brms)
library(ggpubr)
library(stringr)
library(lme4)
library(stats)
library(patchwork)  #https://github.com/thomasp85/patchwork
library(scales)
library(plot.matrix)
library(RColorBrewer)
library(corrplot)
library(ggpubr)
library(png)
library(readr)

source(here("Model", "hypoGeneratorFunctions.R")) 

##########################################################################
# Functions
############################################################################

ModelLLrealRT <- function(AllData, learningrate){
  ModelToData <- runHypoLearnerOnAllData(start_particle = generateStartParticles(nOfBars = 7, nOfColors = 10),
                                         AllData,
                                         learningrate)
  #####################################
  # exlude all data where the participant was incorrect or th real RT is larger than 10/-10
  ############################################
  ModelToData <- (subset(ModelToData, realAccuracy == 1 & realRT <= RTCutoff))
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
  # plot = ggplot(ModelToData, aes(x = realRT, y = timeasRT, color = Structure))+
  #   geom_point()+
  #   ggtitle(paste(LL, learningrate))
  # print(plot)
  return(-LL)
}


###############################################################################################################
###############################################################################################################
# run hypo learner on data
####################################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)
RTCutoff = 10
learningrate = seq(from = 0.0, to = 1, by = 0.01) #c(0.001, seq(from = 0.01, to = 1, by = 0.01))


###########################
# find the best hyperparameter (e.g. the learningrate)
###############################
subjects <- unique(d$subject_id)

LLdf <- data.frame(subject  = character(),
                  learningrate = numeric(),
                  LL = numeric()
)

for (subject in subjects){
  
  LL <- c(1:length(learningrate))
  print("subject:")
  print(subject)
  
  for (i in 1:length(learningrate)){
    print(i)
    allData <- subset(d, subject_id == subject)
    LL[i] <- ModelLLrealRT(allData, learningrate[i])
  }
  
  LLdfTemp <- data.frame(subject  = subject,
                       learningrate = learningrate,
                       LL = LL)
  LLdf <- rbind(LLdf, LLdfTemp)
  
}

ggplot(LLdf, aes(x = learningrate, y = LL, color = subject))+
  geom_point()+
  geom_line()+
  theme(legend.position = "none")


#########################################################
# identify the learningrate with the minimal LL
######################################################
bestLRs <- data.frame(LR = numeric(),
                     subject = character(),
                     minLL = numeric())

for(currentsubject in subjects){
  data <- subset(LLdf, subject == currentsubject)
  minIndex <- which.min(data$LL)
  LR <- data$learningrate[minIndex]
  minLL <- data$LL[minIndex]
  bestLR <- data.frame(LR = LR,
                      subject = currentsubject,
                      minLL = minLL)
  bestLRs <- rbind(bestLRs, bestLR)
}

ggplot(bestLRs, aes(x = LR)) + 
  geom_histogram()
ggplot(bestLRs, aes(x = minLL, y = LR)) + 
  geom_point()

write_csv(bestLRs, here("Model", "FittedModelData", "hypoGeneratorBestHyperPsIndLL_afterBugFix05_23.csv"))
bestLRs <- read.csv(here("Model", "FittedModelData", "hypoGeneratorBestHyperPsIndLL_afterBugFix05_23.csv"))

