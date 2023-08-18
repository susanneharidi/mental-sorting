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
library(readr)
library(patchwork)  #https://github.com/thomasp85/patchwork
library(scales)
library(plot.matrix)
library(RColorBrewer)
library(corrplot)
library(png)

source(here("Model", "hypoMutatorFunctions.R")) 

##########################################################################
#load the data and set some parameters
##############################################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)


##############################################################
# run hypo mutator on all data
################################################################

# hyperparameters from 
bestHyperPs <- read.csv(here("Model", "FittedModelData", "hypoMutatorBestHyperPsReview_with_more_mutations.csv"))

###########################################################
# run this on one participant at a time
############################################################

#subjects <- unique(d$subject_id)
subjects <- unique(bestHyperPs$subject)

# load model to dataframe
ModelToData <- read.csv(here("Model", "FittedModelData", "hypoMutatorfitted_BucketSortIndHypoMoreMutations.csv"))

# or recreate it with the code from below

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

number = 1
for (currentsubject in subject){
  print(currentsubject)
  print(number)
  number = number + 1
  currentHyperPs <- subset(bestHyperPs, subject == currentsubject)
  nStartParticles <- currentHyperPs$nStartParticles
  particles_to_mutate = currentHyperPs$particles_to_mutate
  currtentData <- subset(d, subject_id == currentsubject)
  
  print(currentHyperPs$nStartParticles)
  
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






#################################
# save model to data dataframe for further analysis
###################################

# write_csv(ModelToData, here("Model", "FittedModelData", "hypoMutatorfitted_BucketSortIndHypoMoreMutations.csv"))

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
##############
# real RT
##############
modelData <- (subset(ModelToData, realAccuracy == 1 & realRT <= 10 ))

ModelvsDataPlot <- ggplot(modelData, aes(x = time, y = realRT, color = Structure))+
  geom_point(alpha = 0.2)+
  theme(legend.position = "None")+
  geom_smooth(data = modelData, aes(x = time, y = realRT), method = lm)+
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

ModelvsDataPlot

################### Plotting it all together

ModelOutputPlot <- (plotNoB  + plotrtTrials)/(plotFinalThresholds + plotFinalConnections + ModelvsDataPlot) +
  plot_annotation(tag_levels = "A")+ 
  plot_layout(guides = 'collect')
ModelOutputPlot
