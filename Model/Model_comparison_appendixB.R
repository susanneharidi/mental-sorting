
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
dMutatorOriginal <- read.csv(here("Model", "FittedModelData", "hypoMutatorfitted_BucketSortIndHypoPsLL.csv"))
dMutatorMoreMutations <- read.csv(here("Model", "FittedModelData", "hypoMutatorfitted_BucketSortIndHypoMoreMutations.csv"))
dMutatorNrOfMutationsFixedValues <- read.csv(here("Model", "FittedModelData", "hypoMutatorfitted_BucketSort_Nr_of_mutations_fixed_values.csv"))
dMutatorNrOfMutations <- read.csv(here("Model", "FittedModelData", "hypoMutatorfitted_BucketSortIndHypoPsLL_Nr_of_mutations.csv"))


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
# results 12.05.23 LL 0-1 after bug fix
#            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     1.10      0.07     0.95     1.24 1.00     2634     4551
# time          0.38      0.02     0.34     0.42 1.01      872     1908


dMutatorModelMoreMutations <- (subset(dMutatorMoreMutations, realAccuracy == 1 & realRT <= RTCutoff))
modelMutatorMoreMutations <- brm(realRT ~ time + (time |participant), 
                    data = dMutatorModelMoreMutations,
                    iter = 6000,
                    seed = seed, # added 12.05.2023
                    family = exgaussian(), # added 12.05.2023
                    save_pars = save_pars(all = TRUE),
                    file = "hypoMutatorFinalrawRT_BucketSortIndHyperPs_more_mutations")
summary(modelMutatorMoreMutations)


# original hypothesis mutator model
dMutatorModelOriginal <- (subset(dMutatorOriginal, realAccuracy == 1 & realRT <= RTCutoff))
modelMutatorOriginal <- brm(realRT ~ time + (time |participant), 
                            data = dMutatorModelOriginal,
                            iter = 6000,
                            seed = seed, # added 12.05.2023
                            family = exgaussian(), # added 12.05.2023
                            save_pars = save_pars(all = TRUE),
                            file = "hypoMutatorFinalrawRT_BucketSortIndHypoLL_afterBugFix05_23")
summary(modelMutatorOriginal)

# 1 Mutation
dMutatorModel_1 <- (subset(dMutatorNrOfMutationsFixedValues, realAccuracy == 1 & realRT <= RTCutoff & Nr_of_mutations == 1))
modelMutator_1 <- brm(realRT ~ time + (time |participant), 
                      data = dMutatorModel_1,
                      iter = 6000,
                      seed = seed, # added 12.05.2023
                      family = exgaussian(), # added 12.05.2023
                      save_pars = save_pars(all = TRUE),
                      file = "hypoMutatorFinalrawRT_BucketSorthypoBucketSortIndHypoPsLLReview_Nr_of_mutation_1")
summary(modelMutator_1)

# 10 Mutations
dMutatorModel_10 <- (subset(dMutatorNrOfMutationsFixedValues, realAccuracy == 1 & realRT <= RTCutoff & Nr_of_mutations == 10))
modelMutator_10 <- brm(realRT ~ time + (time |participant), 
                       data = dMutatorModel_10,
                       iter = 6000,
                       seed = seed, # added 12.05.2023
                       family = exgaussian(), # added 12.05.2023
                       save_pars = save_pars(all = TRUE),
                       file = "hypoMutatorFinalrawRT_BucketSorthypoBucketSortIndHypoPsLLReview_Nr_of_mutation_10")
summary(modelMutator_10)

# 100 Mutations
dMutatorModel_100 <- (subset(dMutatorNrOfMutationsFixedValues, realAccuracy == 1 & realRT <= RTCutoff & Nr_of_mutations == 100))
modelMutator_100 <- brm(realRT ~ time + (time |participant), 
                        data = dMutatorModel_100,
                        iter = 6000,
                        seed = seed, # added 12.05.2023
                        family = exgaussian(), # added 12.05.2023
                        save_pars = save_pars(all = TRUE),
                        file = "hypoMutatorFinalrawRT_BucketSorthypoBucketSortIndHypoPsLLReview_Nr_of_mutation_100")
summary(modelMutator_100)

# 1000 Mutations
dMutatorModel_1000 <- (subset(dMutatorNrOfMutationsFixedValues, realAccuracy == 1 & realRT <= RTCutoff & Nr_of_mutations == 1000))
modelMutator_1000 <- brm(realRT ~ time + (time |participant), 
                         data = dMutatorModel_1000,
                         iter = 6000,
                         seed = seed, # added 12.05.2023
                         family = exgaussian(), # added 12.05.2023
                         save_pars = save_pars(all = TRUE),
                         file = "hypoMutatorFinalrawRT_BucketSorthypoBucketSortIndHypoPsLLReview_Nr_of_mutation_1000")
summary(modelMutator_1000)


# 1000 Mutations
dMutatorModel_10000 <- (subset(dMutatorNrOfMutationsFixedValues, realAccuracy == 1 & realRT <= RTCutoff & Nr_of_mutations == 10000))
modelMutator_10000 <- brm(realRT ~ time + (time |participant), 
                          data = dMutatorModel_10000,
                          iter = 6000,
                          seed = seed, # added 12.05.2023
                          family = exgaussian(), # added 12.05.2023
                          save_pars = save_pars(all = TRUE),
                          file = "hypoMutatorFinalrawRT_BucketSorthypoBucketSortIndHypoPsLLReview_Nr_of_mutation_10000")
summary(modelMutator_10000)

# fitted number of mutations 
dMutatorModelfitted <- (subset(dMutatorNrOfMutations, realAccuracy == 1 & realRT <= RTCutoff))
modelMutatorfitted <- brm(realRT ~ time + (time |participant), 
                          data = dMutatorModelfitted,
                          iter = 6000,
                          seed = seed, # added 12.05.2023
                          family = exgaussian(), # added 12.05.2023
                          save_pars = save_pars(all = TRUE),
                          file = "hypoMutatorFinalrawRT_BucketSortIndHyperPs_mut_per_trial_one_hypo")
summary(modelMutatorfitted)


lGenerator <- loo_R2(modelGenerator)[1]
lMutator_1 <- loo_R2(modelMutator_1)[1]
lMutator_10 <- loo_R2(modelMutator_10)[1]
lMutator_100 <- loo_R2(modelMutator_100)[1]
lMutator_1000 <- loo_R2(modelMutator_1000)[1]
lMutator_10000 <- loo_R2(modelMutator_10000)[1]
lMutator_fitted <- loo_R2(modelMutatorfitted)[1]
lMutator_original <- loo_R2(modelMutatorOriginal)[1]
lMutatorMoreMutations <- loo_R2(modelMutatorMoreMutations)[1]


dl <- data.frame(r2 = c(lMutator_1,
                        lMutator_10,
                        lMutator_100,
                        lMutator_1000,
                        lMutator_10000,
                        lMutator_fitted,
                        lMutatorMoreMutations,
                        lMutator_original,
                        lGenerator), 
                 model = c('1Hypothesis Mutator',
                           '10Hypothesis Mutator',
                           '100Hypothesis Mutator',
                           '1000Hypothesis Mutator',
                           '10000Hypothesis Mutator',
                           '2Hypothesis Mutator fitted',
                           '3Hypothesis Mutator original more Mutations',
                           '4Hypothesis Mutator original',
                           '5Hypothesis Generator'))


dl
# results 16.08.23 
# r2                                       model
# 1 0.5043392                         1Hypothesis Mutator
# 2 0.4998010                        10Hypothesis Mutator
# 3 0.4527709                       100Hypothesis Mutator
# 4 0.4459111                      1000Hypothesis Mutator
# 5 0.4515842                     10000Hypothesis Mutator
# 6 0.5281273                  2Hypothesis Mutator fitted
# 7 0.5759633 3Hypothesis Mutator original more Mutations
# 8 0.5659011                4Hypothesis Mutator original
# 9 0.5907812                       5Hypothesis Generator

############################################################
# Calulate BFs for the models
############################################################
BFGeneratorMutator <- bayes_factor(modelGenerator, modelMutatorMoreMutations)
BFGeneratorMutator
# 16.08.2023
# BF = 1641611268309706994344468280240602802880868840.00000

############################################################
# plot the results
############################################################
cbbPalette <- c("#FFFFE0", "#FFFF99", "#FFFF00", "#FFEF00", "#FFDF00", "#FFD700", "#7CFC00", "#00C85A", "#049CCC")
# cbbPalette <- c("#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#31A354", "#006D2C", "#00C85A", "#049CCC") # green tones
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
    labels = c("Mutator B 10000", "Mutator B 1000", "Mutator B 100", "Mutator B 10", "Mutator B 1", "Mutator B fitted", "Mutator A fitted", "Mutator A", "Generator"))+
  #add xlab
  xlab("Models")+
  #add ylab
  ylab(expression(paste("Loo R"^2*" value")))+
  #change font and text type
  theme(text = element_text(size = 20, family = "sans"))
ModelComparisonFigure

ggsave(here("Figures", "ModelComparisonFigure_appendixB.pdf"), ModelComparisonFigure,  
       width = 10, 
       height = 6)
# last saved figure 31.10.2022


