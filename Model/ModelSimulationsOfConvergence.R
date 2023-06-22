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

##################################################################################
# set up hypothesis space
###########################################################################

start_particle = generateStartParticles(nOfBars = 7, nOfColors = 10)

total_runs = 10
trials = 500

# Since the average best learning rate for participants was 0.12, thats what I will use here
LR = 0.12
# same for particles 
nStartParticles = 12
particles_to_mutate = 0.2


###############################
# Simulate the hypothesis Generator
#################################

simDataGeneratorQuery <- runSimulation(start_particle, # the current particle
              total_runs, # int
              trials, # int
              learningrate = LR,
              structure_query = TRUE, 
              structure_task = FALSE)

timesGeneratorQuery <- simDataGeneratorQuery[[1]]
accuracyGeneratorQuery <- simDataGeneratorQuery[[2]]
directionsGeneratorQuery <- simDataGeneratorQuery[[3]]
NoBGeneratorQuery <- simDataGeneratorQuery[[4]]
thresholdsGeneratorQuery <- simDataGeneratorQuery[[5]]
connectionsGeneratorQuery <- simDataGeneratorQuery[[6]]

# make a dataframe of the final hypotheses
finalHypothesesGeneratorQuery = data.frame(
  threshold = thresholdsGeneratorQuery[trials,],
  connections =connectionsGeneratorQuery[trials,],
  startAtSmall = directionsGeneratorQuery[trials,],
  time = timesGeneratorQuery[trials,],
  correct = accuracyGeneratorQuery[trials,],
  run = seq(total_runs),
  Structure = "Query",
  Model = "Generator")

# in the sequence structure e-f-g is always connected

simDataGeneratorSequence <- runSimulation(start_particle, # the current particle
                              total_runs, # int
                              trials, # int
                              learningrate = LR,
                              structure_query = FALSE, 
                              structure_task = TRUE)

timesGeneratorSequence <- simDataGeneratorSequence[[1]]
accuracyGeneratorSequence <- simDataGeneratorSequence[[2]]
directionsGeneratorSequence <- simDataGeneratorSequence[[3]]
NoBGeneratorSequence <- simDataGeneratorSequence[[4]]
thresholdsGeneratorSequence <- simDataGeneratorSequence[[5]]
connectionsGeneratorSequence <- simDataGeneratorSequence[[6]]

for (i in 1:total_runs){
  print(i)
  temp <- data.frame(
    Connection = connectionsGeneratorSequence[,i],
    trial = seq(trials),
    thresholds = thresholdsGeneratorSequence[,i],
    correct = accuracyGeneratorSequence[,i],
    run = i)
  if (i == 1){
    HypothesisGeneratorSequence = temp
  }else{
    HypothesisGeneratorSequence = rbind(HypothesisGeneratorSequence, temp)
  }
}

# make plot of development over time

ggplot(HypothesisGeneratorSequence, aes(x = trial, y = run, color = Connection))+
  geom_point()+
  facet_grid(rows = vars(as.factor(correct)))


finalHypothesesGeneratorSequence = data.frame(
  threshold = thresholdsGeneratorSequence[trials,],
  connections =connectionsGeneratorSequence[trials,],
  startAtSmall = directionsGeneratorSequence[trials,],
  time = timesGeneratorSequence[trials,],
  correct = accuracyGeneratorSequence[trials,],
  run = seq(total_runs),
  Structure = "Sequence",
  Model = "Generator")

finalHypothesesGenerator = rbind(finalHypothesesGeneratorSequence, finalHypothesesGeneratorQuery)

############################################
# Simulate hypothesis Mutator
#############################################
source(here("Model", "hypoMutatorFunctions.R")) 

simDataMuatorQuery <- runSimulation(nStartParticles, # generated particles
                total_runs, # int
                trials, # int
                particles_to_mutate, # int, number of particles that get mutated each trial
                structure_query = TRUE, 
                structure_task = FALSE)

timesMuatorQuery <- simDataMuatorQuery[[1]]
accuracyMuatorQuery <- simDataMuatorQuery[[2]]
finalHypothesesMuatorQuery <- simDataMuatorQuery[[3]]
thresholdsMutatorQuery <- finalHypothesesMuatorQuery$threshold
connectionsMutatorQuery <- finalHypothesesMuatorQuery$connections
directionsMutatorQuery <- finalHypothesesMuatorQuery$startAtSmall
NoBMuatorQuery <- simDataMuatorQuery[[4]]


simDataMuatorSequence <- runSimulation(nStartParticles, # generated particles
                                    total_runs, # int
                                    trials, # int
                                    particles_to_mutate, # int, number of particles that get mutated each trial
                                    structure_query = FALSE, 
                                    structure_task = TRUE)

timesMuatorSequence <- simDataMuatorSequence[[1]]
accuracyMuatorSequence <- simDataMuatorSequence[[2]]
finalHypothesesMuatorSequence <- simDataMuatorSequence[[3]]
thresholdsMutatorSequence <- finalHypothesesMuatorSequence$threshold
connectionsMutatorSequence <- finalHypothesesMuatorSequence$connections
directionsMutatorSequence <- finalHypothesesMuatorSequence$startAtSmall
NoBMuatorSequence <- simDataMuatorSequence[[4]]

finalHypothesesMuatorSequence$Structure = "Sequence"
finalHypothesesMuatorQuery$Structure = "Query"
finalHypothesesMuator = rbind(finalHypothesesMuatorSequence, finalHypothesesMuatorQuery)
finalHypothesesMuator$Model = "Mutator"

# Put the final hypotheses of both models in one dataframe
finalHypotheses = rbind(finalHypothesesMuator, finalHypothesesGenerator)

############################################
# plot the output of the mutator model 
############################################
# plot all thresholds
cbbPalette <- c("#0072B2", "#D55E00")
pd <- position_dodge(.2)


plotFinalThresholdsMutator <- ggplot(finalHypotheses, aes(x = threshold, fill = Structure, color = Structure)) +
  #points
  geom_histogram(alpha = .4, position="dodge")+
  #minimal theme
  theme_minimal()+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #add xlab
  xlab("Final Threshold")+
  #add ylab
  ylab('Count')+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"), legend.position = "top")+
  facet_grid(rows = vars(Model), scales = "free_y")

print(plotFinalThresholdsMutator)

plotFinalDirectionsMutator <- ggplot(finalHypotheses , aes(x = factor(startAtSmall), fill = Structure, color = Structure)) +
  #points
  geom_bar(alpha = .4, position = position_dodge2(preserve = "single"))+
  #minimal theme
  theme_minimal()+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #add xlab
  xlab("Final Direction")+
  #add ylab
  ylab('Count')+
  #change fonts
  theme(text = element_text(size = 20, family = "sans"), legend.position = "top")+
  facet_grid(rows = vars(Model), scales = "free_y")

print(plotFinalDirectionsMutator)


# labe the correct connection for the structure condition:
finalHypotheses$correctConnection = "d3 other"
for (i in 1:length(finalHypotheses$connection)){
  if (finalHypotheses$connection[i] == "efg"){
    finalHypotheses$correctConnection[i] = "a3correct"
  }else if(finalHypotheses$connection[i] == "fg"){
    finalHypotheses$correctConnection[i] = "b2correct"
  }else if(finalHypotheses$connection[i] == "ef"){
    finalHypotheses$correctConnection[i] = "b2correct"
  }else if(length(string_to_vec(finalHypotheses$connection[i])) <= 1){
    finalHypotheses$correctConnection[i] = "cno"
  }else if(length(string_to_vec(finalHypotheses$connection[i])) == 2){
    finalHypotheses$correctConnection[i] = "d2other"
  }
}

plotFinalConnectionsMutator <- ggplot(finalHypotheses, aes(x = correctConnection, fill = Structure, color = Structure)) +
  #points
  geom_bar(alpha = .4, position = position_dodge2(preserve = "single"))+
  #minimal theme
  theme_minimal()+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #add xlab
  xlab("Final Connections")+
  #add ylab
  ylab('Count')+
  scale_x_discrete(labels = c("3 correct", "2 correct", "none", "2 other", "3 other"))+ 
  #change fonts
  theme(text = element_text(size = 20, family = "sans"), legend.position = "top")+
  facet_grid(rows = vars(Model), scales = "free_y")

print(plotFinalConnectionsMutator)


######
# get model recovery plot
######

AllCorrelationsdf = read.csv(here("Data", "ModelRecoveryCorrelations_afterBugFix05_23.csv"))

textsize = 20 
modelRecoveryPlot <- ggplot(AllCorrelationsdf, aes(x=model1, y=model2, fill=correlation))+
  geom_tile()+
  xlab("Model 2")+
  #add ylab
  ylab('Model 1')+
  #change fonts
  #minimal theme
  theme_minimal()+
  labs(fill = "Correlation")+
  theme(text = element_text(size = textsize, family = "sans"))+
  theme(legend.position = "top")+
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1))

modelRecoveryPlot


ModelConvergencePlot <- ((plotFinalThresholdsMutator + plotFinalDirectionsMutator)/ 
                                         (plotFinalConnectionsMutator + (modelRecoveryPlot +  plot_layout(guides = 'keep'))))+
  plot_annotation(tag_levels = "A")+ 
  plot_layout(guides = 'collect') &  theme(legend.position = 'top')
ModelConvergencePlot

ggsave(here("Figures", "FigureB2_Simulation_model_Convergence.pdf"), ModelConvergencePlot,  
       width = 14, 
       height = 10)

# just some script for debugging

trial = generatetrial(structure_task = TRUE)

noisy_sort_backwards(failTrial, threshold = 7, connections = 'efg')
