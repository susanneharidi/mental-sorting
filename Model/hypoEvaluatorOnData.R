# Run the hypothesis evaluator on the actual data
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

source(here("Model", "hypoEvaluatorFunctions.R"))   

#####################################################################################
# run sorter on real data
#####################################################################################
###########################################################
# get data and set parameters
###########################################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)
particles_to_eliminate =  0.01 
nofColors = 10

##################################################################################
# set up hypothesis space
###########################################################################

#all hypotheses of size 2
two <- t(combn(letters[1:nofColors], 2))

#all hypotheses of size 3
three <- t(combn(letters[1:nofColors], 3))

# all Start at small
startAtSmall <- c(0, 1)

#create all possible hypotheses
connections <- c('a', paste0(two[,1], two[,2]), paste0(three[,1], three[,2], three[,3]))

#generate particle hypotheses
particles <- expand.grid(threshold = 1:7, connections = connections, startAtSmall = startAtSmall)


################################################################
# run the model on all Data
##############################################################
output <- runHypoSortOnAllData(start_particles = particles,
                                AllData = d,
                                particles_to_eliminate)
ModelToData <- output[[1]]
allfinalHypo <- output[[2]]

#################################
# save model to data dataframe for further analysis
###################################
write_csv(ModelToData, here("Model", "FittedModelData", paste(as.character(particles_to_eliminate), "hypoEvaluatorfitted_BucketSort.csv", sep = "_")))
# load the data if its already been saved before
ModelToData <- read.csv(here("Model", "FittedModelData", paste(as.character(particles_to_eliminate), "hypoEvaluatorfitted_BucketSort.csv", sep = "_")))

############################################
# plot the model output
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
  xlab("threshold")+
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

# time over sequence length
rtNoB = ddply(ModelToData , ~NoB + Structure , summarize, meanTime = mean(time), se = se(time))
pd <- position_dodge(.2)
plotNoB <- ggplot(rtNoB, aes(x = NoB, y = meanTime, fill = Structure, color = Structure)) +
  geom_line(position = pd, size = 0.5) +
  #points
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
  ggtitle("Model time per NoB")

print(plotNoB)


rtTrialsNoB = ddply(ModelToData , ~NoB + trial + Structure , summarize, meanTime = mean(time), se = se(time))
# time over trials
plotrtTrialsNoB <- ggplot(rtTrialsNoB, aes(x = trial, y = meanTime, fill = factor(NoB), color = factor(NoB))) +
  geom_line(position = pd, size = 0.5) +
  #points
  geom_point(position = pd)+
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
