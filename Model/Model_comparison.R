
#####################################################
# load libraries
#####################################################

library(ggplot2)
library(dplyr)
library(here)
library(brms)

#####################################
# load data and set some parameters
#####################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)
RTCutoff <- 10 # in s
learningrate = seq(from = 0.01, to = 0.1, by = 0.01)

# this part is necessary, so I can add the correct columns to the model dataframe
relevantD <- subset(d, Condition == "Sort" & Stimulus_type == "bars")
relevantD1 <- subset(relevantD, Structure == "Query")
relevantD2 <- subset(relevantD, Structure == "Sequence")
relevantD3 <- subset(relevantD, Structure == "None")
relevantD <- rbind(relevantD1, relevantD2, relevantD3)
#############################################
# Load model data
#############################################

dGenerator <- read.csv(here("Model", "FittedModelData", "0.09_hypoGeneratorfitted_BucketSort.csv"))
dMutator <- read.csv(here("Model", "FittedModelData", "15_1_hypoMutatorfitted_BucketSort.csv"))
dEvaluator <- read.csv(here("Model", "FittedModelData", "0.01_hypoEvaluatorfitted_BucketSort.csv"))

# and add the real RTs
dGenerator$rawRT <- relevantD$rt
dMutator$rawRT <- relevantD$rt
dEvaluator$rawRT <- relevantD$rt

dGenerator$queryRT <- relevantD$queryRT
dMutator$queryRT <- relevantD$queryRT
dEvaluator$queryRT <- relevantD$queryRT

#################################################
# fit the models
##################################################
setwd(here("Model", "hypoFittingBrmModels"))
dGeneratorModel <- (subset(dGenerator, realAccuracy == 1 & rawRT <= RTCutoff & queryRT <= RTCutoff))
modelGenerator <- brm(rawRT ~ time + (time |participant), 
                      data = dGeneratorModel,
                      iter = 6000,
                      save_pars = save_pars(all = TRUE),
                      file = "hypoGeneratorFinalrawRT_BucketSort")

dMutatorModel <- (subset(dMutator, realAccuracy == 1 & rawRT <= RTCutoff & queryRT <= RTCutoff))
modelMutator <- brm(rawRT ~ time + (time |participant), 
                    data = dMutatorModel,
                    iter = 6000,
                    save_pars = save_pars(all = TRUE),
                    file = "hypoMutatorFinalrawRT_BucketSort")

dEvaluatorModel <- (subset(dEvaluator, realAccuracy == 1 & rawRT <= RTCutoff & queryRT <= RTCutoff))
modelEvaluator <- brm(rawRT ~ time + (time |participant), 
                    data = dEvaluatorModel,
                    iter = 6000,
                    save_pars = save_pars(all = TRUE),
                    file = "hypoEvaluatorFinalrawRT_BucketSort")

lGenerator <- loo_R2(modelGenerator)[1]
lMutator <- loo_R2(modelMutator)[1]
lEvaluator <- loo_R2(modelEvaluator)[1]

dl <- data.frame(r2 = c(lEvaluator, lMutator, lGenerator), 
               model = c('1Hypothesis Evaluator', '2Hypothesis Mutator', '3Hypothesis Generator'))

dl
# results 14.08.06.2022
#       r2            model
# 1 0.6158147 1Hypothesis Evaluator
# 2 0.6399365   2Hypothesis Mutator
# 3 0.6637963 3Hypothesis Generator

############################################################
# Calulate BFs for the models
############################################################
BFGeneratorMutator <- bayes_factor(modelGenerator, modelMutator)
# results 14.08.2022
# BF = 383460497425364809408824046840644622820842466468246640264082260884806840460262662886660288864202488.00000
BFGeneratorEvaluator <- bayes_factor(modelGenerator, modelEvaluator)
# results 14.08.2022
# BF = 53672460394143118196466842680806646448806884684222480422640200002266820442004844004206088000668226864680442484240464640222820808022468642444066246200266802426248620646826204680808288802660662.00000
BFGeneratorEvaluator <- bayes_factor(modelMutator, modelEvaluator)
# results 14.08.2022
# BF = 141816691697024135468462824088648404684040268482428204462006082208462404664620462686844460202.00000

############################################################
# plot the results
############################################################
cbbPalette <- c("#FF5050", "#00C85A", "#049CCC")
ModelComparisonFigure <- ggplot(dl, aes(x = model, y = r2, color = model, fill = model))+
        geom_bar(stat="identity")+
        #minimal theme
        theme_minimal()+
        theme(legend.position = "None")+
        #change fill
        scale_fill_manual(values = cbbPalette)+
        #change color
        scale_color_manual(values = cbbPalette)+
        scale_x_discrete(guide = guide_axis(n.dodge = 2), 
                   labels = c("Evaluator", "Mutator", "Generator"))+
        #add xlab
        xlab("Models")+
        #add ylab
        ylab(expression(paste("Loo R"^2*" value")))+
        #change font and text type
        theme(text = element_text(size = 20, family = "sans"))
ModelComparisonFigure
ggsave(here("Figures", "ModelComparisonFigure.png"), ModelComparisonFigure,  
       width = 3, 
       height = 5)
# last saved figure 06.07.2022



##############################################################################
# model how well the mutator and the generator explain each other
##############################################################################
setwd(here("Model", "hypoFittingBrmModels"))
dAll <- data.frame(participant = dGenerator$participant,
                   timeGenerator = dGenerator$time,
                   timeMutator = dMutator$time)
modelGeneratorOnMutator <- brm(timeMutator ~ timeGenerator + (timeGenerator|participant), 
                      data = dAll,
                      iter = 6000,
                      save_pars = save_pars(all = TRUE),
                      file = "modelGeneratorOnMutator_BucketSort")
summary(modelGeneratorOnMutator)
# results 14.06.2022
#                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept         0.08      0.03     0.02     0.13 1.00    11382    10034
# timeGenerator     1.01      0.01     0.99     1.02 1.00     8271     8026
loo_R2(modelGeneratorOnMutator)[1]
# results 14.06.2022
# 0.928479

modelMutatorOnGenerator <- brm(timeGenerator ~ timeMutator + (timeMutator|participant), 
                               data = dAll,
                               iter = 6000,
                               save_pars = save_pars(all = TRUE),
                               file = "modelMutatorOnGenerator_BucketSort")
summary(modelMutatorOnGenerator)
# results 14.06.2022
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept       0.38      0.02     0.33     0.42 1.00    13894     9180
# timeMutator     0.92      0.00     0.91     0.93 1.00     9223     8947

loo_R2(modelMutatorOnGenerator)[1]
# results 14.06.2022
# 0.9278038
