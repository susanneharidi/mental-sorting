
#####################################################
# load libraries
#####################################################

library(ggplot2)
library(dplyr)
library(here)
library(brms)

#############################################
# functions
###############################################
# calculate the BIC
# BIC = -2 * loglikelihood + d * log(N)
BIC <- function(ModelToData, d){ # d is the number of free parameters in the model
  N <- length(ModelToData$realRT) # sample size
  subjects <- unique(ModelToData$participant)
  LL <- 0
  for (subject in subjects){
    print("subject:")
    print(subject)
    
    Datasubset <- subset(ModelToData, subject == subject)
    sd <- sd(Datasubset$realRT, na.rm = FALSE)
    
    output <- lm(realRT ~ time, data = Datasubset)
    intercept <- as.numeric(output$coefficients[1])
    beta <- as.numeric(output$coefficients[2])
    Datasubset$timeasRT <- intercept + Datasubset$time*beta
    
    for (i in 1:length(Datasubset$timeasRT)){
      LL <- LL + dnorm(x = Datasubset$realRT[i], mean = Datasubset$timeasRT[i], sd = sd, log = TRUE)
    }
    print("LL:")
    print(LL)
  }
  bic = -2 * LL + d * log(N)
  return(c(bic, LL))
}


#####################################
# load data and set some parameters
#####################################
d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))
d$subject_id <- factor(d$subject_id)
RTCutoff <- 10 # in s
seed <- 2022
#############################################
# Load model data
#############################################

dGenerator <- read.csv(here("Model", "FittedModelData", "hypoGeneratorfitted_BucketSortIndHypoPsLL_afterBugFix05_23.csv"))
dMutator <- read.csv(here("Model", "FittedModelData", "hypoMutatorfitted_BucketSortIndHypoPsLL_afterBugFix05_23.csv"))


#################################################
# fit the models
##################################################
setwd(here("Model", "hypoFittingBrmModels"))
dGeneratorModel <- (subset(dGenerator, realAccuracy == 1 & realRT <= RTCutoff))
modelGenerator <- brm(realRT ~ time + (time |participant), 
                      data = dGeneratorModel,
                      iter = 6000,
                      seed = seed, # added 12.05.2023
                      family = exgaussian(), # added 12.05.2023
                      save_pars = save_pars(all = TRUE),
                      file = "hypoGeneratorFinalrawRT_BucketSortIndHypoLL_afterBugFix05_23")
summary(modelGenerator)
# results 28.10.2022 LL 0-1
#            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.20      0.11    -0.01     0.41 1.00     1822     3319
# time          0.57      0.02     0.52     0.61 1.01      778     1551

# results 12.05.23 LL 0-1 after bug fix
#            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     1.10      0.07     0.95     1.24 1.00     2634     4551
# time          0.38      0.02     0.34     0.42 1.01      872     1908


dMutatorModel <- (subset(dMutator, realAccuracy == 1 & realRT <= RTCutoff))
modelMutator <- brm(realRT ~ time + (time |participant), 
                    data = dMutatorModel,
                    iter = 6000,
                    seed = seed, # added 12.05.2023
                    family = exgaussian(), # added 12.05.2023
                    save_pars = save_pars(all = TRUE),
                    file = "hypoMutatorFinalrawRT_BucketSortIndHypoLL_afterBugFix05_23")
summary(modelMutator)
# results 13.10.2022
#            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.26      0.11     0.06     0.47 1.00     2067     3920
# time          0.56      0.02     0.51     0.61 1.00      914     1546
# results 18.10.2022 LL
#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.27      0.10     0.06     0.48 1.00     2028     3543
# time          0.56      0.02     0.51     0.60 1.00      977     2194

# results 12.05.23 LL 0-1 after bug fix
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     1.13      0.07     1.00     1.26 1.00     4132     7220
# time          0.40      0.02     0.36     0.43 1.01     1030     2494


lGenerator <- loo_R2(modelGenerator)[1]
lMutator <- loo_R2(modelMutator)[1]



dl <- data.frame(r2 = c(lMutator, lGenerator), 
                 model = c('1Hypothesis Mutator', '2Hypothesis Generator'))

dl
# 28.10.2022
#       r2                model
# 1 0.6365322 Hypothesis Mutator
# 2 0.6480957 Hypothesis Generator

# results 12.05.23 LL 0-1 after bug fix

# r2                 model
# 1 0.5659963   1Hypothesis Mutator
# 2 0.5906816 2Hypothesis Generator

############################################################
# Calulate BFs for the models
############################################################
BFGeneratorMutator <- bayes_factor(modelGenerator, modelMutator)
BFGeneratorMutator
# Results 31.10.2022
# BF = 1540457581781216155846866286620866082048484666.00000

# results 12.05.23 LL 0-1 after bug fix
# BF = 10589538106590291565604666866826620682288262248864040806.00000


###############################################################
# Calculate the BIC (the smaller the better)
#############################################################
BICMuator <- BIC(dMutatorModel, d = 2)
BICMuator
# Bic = 1969099.4
# LL = -984540.9

# results 12.05.23 LL 0-1 after bug fix
# Bic = 1972804.6 
# LL = -986393.5

BICGenerator <- BIC(dGeneratorModel, d = 1)
BICGenerator
# Bic = 1959558.5 
# LL = -979774.8

# results 12.05.23 LL 0-1 after bug fix
# Bic = 1957380.4 
# LL = -978685.8

BICDifference <- BICMuator[1] - BICGenerator[1]
# 9540.963

# results 12.05.23 LL 0-1 after bug fix
# 15424.2

#################################################################
# add model predictions as regressor
################################################################
setwd(here("Model", "hypoFittingBrmModels"))
# The full Encoding RT model + model time 
mrtallGenerator <- brm(realRT ~ NoB + Structure + time +
                       (NoB + Structure + time |participant), 
                             data = dGeneratorModel,
                             iter = 10000,
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             seed = seed,
                             family = exgaussian(),
                             file = "mrtallGenerator")
summary(mrtallGenerator)
# results 02.11.2022
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.78      0.09     0.59     0.96 1.00     2352     4956
# NoB                   0.54      0.07     0.39     0.68 1.00     2288     4831
# StructureQuery       -0.24      0.06    -0.36    -0.11 1.00     3854     7835
# StructureSequence    -0.15      0.05    -0.25    -0.06 1.00     5581    10117
# time                  0.13      0.04     0.06     0.20 1.00     2849     5652

mrtallMutator <- brm(realRT ~ NoB + Structure + time +
                         (NoB + Structure + time |participant), 
                         data = dMutatorModel,
                         iter = 10000,
                         cores = 4,
                         save_pars = save_pars(all = TRUE),
                         seed = seed,
                         family = exgaussian(),
                         file = "mrtallMutator")
summary(mrtallMutator)
# results 02.11.2022
#                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.80      0.09     0.62     0.98 1.00     2917     5844
# NoB                   0.57      0.07     0.44     0.70 1.00     2928     5824
# StructureQuery       -0.22      0.06    -0.34    -0.11 1.00     5637    10572
# StructureSequence    -0.17      0.04    -0.25    -0.08 1.00     8331    12527
# time                  0.10      0.03     0.04     0.16 1.00     3764     6630

bayes_factor(mrtallGenerator, mrtallMutator)
# results 02.11.2022
# BF = 2033919.94967

mrtall <- brm(realRT ~ NoB + Structure + 
                    (NoB + Structure|participant), 
                     data = dMutatorModel,
                     iter = 10000,
                     cores = 4,
                     save_pars = save_pars(all = TRUE),
                     seed = seed,
                     family = exgaussian(),
                     file = "mrtall")
summary(mrtall)
# results 02.11.2022
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.83      0.09     0.64     1.01 1.00     2986     5846
# NoB                   0.72      0.04     0.64     0.79 1.00     1751     4042
# StructureQuery       -0.32      0.07    -0.47    -0.18 1.00     4166     7648
# StructureSequence    -0.16      0.04    -0.24    -0.07 1.00     8497    12546

looGenerator <- loo_R2(mrtallGenerator)[1]
looMutator <- loo_R2(mrtallMutator)[1]
looNoModel <- loo_R2(mrtall)[1]
# results 02.11.2022
# looGenerator = 0.6669775
# looMutator = 0.6592212
# looNoModel = 0.6428746

############################################################
# plot the results
############################################################
cbbPalette <- c("#00C85A", "#049CCC")
ModelComparisonFigure <- ggplot(dl, aes(x = model, y = r2, color = model, fill = model))+
  geom_bar(stat="identity")+
  #minimal theme
  theme_minimal()+
  theme(legend.position = "None")+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  scale_x_discrete(#guide = guide_axis(n.dodge = 2), 
                   labels = c("Mutator", "Generator"))+
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
# last saved figure 31.10.2022


