# Figures for the supplemets

#########################################################################
# load packages
#########################################################################

library(readr)
library(tidyverse)
library(ggplot2)
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)
library(plyr)
library(dplyr)
library(data.table)
library(gridExtra)
library(ggthemes)
library(ggsignif)
library(brms)
library(lmerTest)
library(here)
library(patchwork)  #https://github.com/thomasp85/patchwork
library(scales)
library(plot.matrix)
library(RColorBrewer)
library(corrplot)
library(ggpubr)
library(png)
library(tidyr)
library(distributional)
library(ggdist)
library(sjPlot)

#########################################################################
# functions
#########################################################################

#standard error function
se <- function(x){sd(x)/sqrt(length(x))}

#########################################################################
# load data and exclude not independent data and people who do not meet the performance threshold
#########################################################################

d <- read.csv(here("Data", "Experiment", "data_to_work_with.csv"))

################################################################################
# Set some parameters
################################################################################

pd <- position_dodge(0.1)
textsize <- 12

#color pallette for three different groups
cbbPalette <- c("grey50", "#0072B2", "#D55E00")

#other parameters in s
rt_cutoff <- 10
seed <- 2022
alpha <- 0.7

############################################################################
# Accuracy Figure
#########################################################################
dcorrect <- data.frame(
  Response = c("correct", "incorrect"),
  value = c(length(subset(d, correct == 1)$correct), length(subset(d, correct == 0)$correct)))

# Compute the position of labels
dcorrect <- dcorrect %>% 
  arrange(desc(Response)) %>%
  mutate(prop = value / sum(dcorrect$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )


# Incorrect Trials piechart
IncorrectPlot <-ggplot(dcorrect, aes(x = Response, y = prop, fill = Response)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_minimal() + 
  ylab('Proportion of Trials')+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  xlab('Response')+
  theme(text = element_text(size = textsize), legend.position = "None") +
  geom_text(aes(y = ypos, label = percent(prop/100)), color = "black", size = 2.5) +
  scale_fill_manual(values = c("darkolivegreen3", "brown2" ))
IncorrectPlot


#########################################################################

dcutoff <- data.frame(
  RT = c("< 10s", "> 10s"),
  value = c(length(subset(d, correct == 1 & rt <= rt_cutoff & queryRT <= rt_cutoff)$correct), length(subset(subset(d, correct == 1), rt > rt_cutoff | queryRT > rt_cutoff)$correct)))

# Compute the position of labels
dcutoff <- dcutoff %>% 
  arrange(desc(RT)) %>%
  mutate(prop = value / sum(dcutoff$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )


RTExclusionPlot <- ggplot(dcutoff, aes(x = RT, y = prop, fill = RT)) +
  geom_bar(stat = "identity", width = 0.5)+
  theme_minimal() + 
  ylab('Proportion of Trials')+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  xlab('RT')+
  theme(text = element_text(size = textsize), legend.position = "None") +
  geom_text(aes(y = ypos, label = percent(prop/100)), color = "black", size = 2.5) +
  scale_fill_manual(values = c( "darkolivegreen3", "brown2"))
RTExclusionPlot

##############################################################################

#for accuracy we need all points, not only correct ones
dp<-subset(d, rt <= rt_cutoff & queryRT <= rt_cutoff)

#failure rate per subject and condition
dp<-ddply(dp, ~Structure + subject_id + Condition, summarize, mistakes = length(correct)-sum(correct))

# calculate the mean data
vline.data <- ddply(dp , ~Structure + Condition, summarize, mean = mean(mistakes))

#plot errors
plotMistakes <- ggplot(dp, aes(x = mistakes, fill = Structure, col = Structure), dotsize = 0) +
  #violin plot
  geom_bar(position = "dodge", alpha = 0.6)+
  geom_vline(aes(xintercept = mean, color = Structure), vline.data, linetype = "longdash")+
  #change fill
  scale_fill_manual(values = cbbPalette)+
  #change color
  scale_color_manual(values = cbbPalette)+
  #change ylab
  ylab('Number of Participants')+
  xlab('Number of Mistakes')+
  #minimal theme
  theme_minimal()+
  #change fonts
  theme(text = element_text(size = textsize, family = "sans"))+
  #no legend
  theme(legend.position = "None")+
  #add title
  #ggtitle("Mean Errors of the Participants")+
  facet_grid(cols = vars(Condition))

plotMistakes

##############################################################################################
# Accuracy over number of Bars
mylabels = c("None", "Query", "Sequence ")
barwidth = 0.7
pErrors = ggplot(subset(d, correct == 0), aes( x= Number_of_Bars, fill=Structure)) + 
  geom_bar(position = "dodge")+
  theme_minimal()+
  #change fonts
  theme(text = element_text(size = textsize, family="sans"))+
  scale_fill_manual(values = cbbPalette,
                    name = "Structure",
                    labels = mylabels)+
  labs(x = "Sequence Length",
       y = "Number of Mistakes")+
  facet_grid(cols = vars(Condition))+
  scale_y_continuous(limits = c(0,75), expand = c(0, 0))

pErrors

# Accuracy over number of Bars
pErrors2 = ggplot(subset(d, correct == 0), aes( x= Correct_Position, fill=Structure)) + 
  geom_bar(position = "dodge")+
  theme_minimal()+
  #change fonts
  theme(text = element_text(size=textsize, family="sans"))+
  scale_fill_manual(values = cbbPalette,
                    name="Structure",
                    labels=mylabels)+
  labs(x = "Queried Position",
       y = "Number of Mistakes")+
  facet_grid(cols = vars(Condition))+
  scale_y_continuous(limits = c(0,75), expand = c(0, 0))
pErrors2

# Accuracy over number of Bars
pCorrect_Position = ggplot(d, aes( x= Correct_Position, fill=Structure)) + 
  geom_bar(position = "dodge")+
  theme_minimal()+
  #change fonts
  theme(text = element_text(size = textsize, family="sans"))+
  scale_fill_manual(values = cbbPalette,
                    name="Structure",
                    labels=mylabels)+
  labs(x = "Queried Position",
       y = "Number of Queries")+
  facet_grid(cols = vars(Condition))+
  scale_y_continuous(expand = c(0, 0))

pCorrect_Position

figure_mistakes <- ggarrange(pCorrect_Position, pErrors2, pErrors,
                             ncol = 3, nrow = 1,
                             common.legend = TRUE, legend = "top")

figure_mistakes


ExlusionPlot <- IncorrectPlot + RTExclusionPlot + plotMistakes +
  plot_layout(ncol = 3, width = c(1, 1, 4))

ExlusionPlot

AccuracyFigure <- ExlusionPlot / figure_mistakes +
  plot_annotation(tag_levels = "A")  +
  plot_layout(guides = 'collect') +
  plot_layout(ncol = 1, heights = c(1, 2))
AccuracyFigure 

ggsave(here("Figures", "FigureC1_Accuracy.pdf"), AccuracyFigure,  
       width = 7, 
       height = 8)

# last saved 22.06.2023 just before revised submission

###################################################################
# Recall RT Figure
###################################################################
#####################################################################
# look at trade-off between recall RT and encoding RT
######################################################################

dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)
dpSummary <- ddply(dp, ~ subject_id + Structure + Condition + Number_of_Bars, summarize, rt = mean(rt), queryRT = mean(queryRT))
RTTradeoff_plot <- ggplot(dpSummary, aes(x = rt, y = queryRT, color = Structure, fill = Structure))+
  geom_point(alpha = 0.2)+
  scale_color_manual(values = cbbPalette)+
  geom_smooth(method = lm, se = FALSE)+
  ylab("Recall RT in s")+
  xlab("Encoding RT in s")+
  scale_x_continuous(breaks = c(5, 10))+
  #scale_y_continuous(breaks = c(5, 10))+
  theme_minimal()+
  theme(text = element_text(size = textsize, family = "sans"),legend.position = "top")+
  facet_grid(rows = vars(Condition), cols = vars(Number_of_Bars))
RTTradeoff_plot 


#######################################################################
# Recall and Encoding RT Model
######################################################################
setwd(here("analyses", "brmModel"))
dp <- subset(d, rt <= rt_cutoff & correct == 1 & otherRT <= rt_cutoff)
m <- brm(rt ~ (Number_of_Bars + Structure + Condition) * Stimulus_type + ((Number_of_Bars + Structure + Condition + Block)*Stimulus_type|subject_id), 
              data = dp,
              iter = 10000,
              cores = 4,
              save_pars = save_pars(all = TRUE),
              seed = seed,
              family = exgaussian(),
              file = "mrtallStructConditionNotLogRecallInteractionFinalEx")

tab_model(m, transform = NULL)
# plotting the effects seperatly
drawsNob <- posterior_samples(m)$b_Number_of_Bars #effect of number of bars on encoding RT
drawsNobRecall <- drawsNob + 
  posterior_samples(m)$`b_Number_of_Bars:Stimulus_typequery`

df1 = data.frame(draws = drawsNob,
                 predictor = "SeqLen",
                 RTType = "Encoding")
df2 = data.frame(draws = drawsNobRecall,
                 predictor = "SeqLen",
                 RTType = "Recall")

drawsStructureQuery <- posterior_samples(m)$b_StructureQuery 
drawsStructureQueryRecall <- drawsStructureQuery + 
  posterior_samples(m)$`b_StructureQuery:Stimulus_typequery`

df3 = data.frame(draws = drawsStructureQuery,
                 predictor = "StructureQuery",
                 RTType = "Encoding")
df4 = data.frame(draws = drawsStructureQueryRecall,
                 predictor = "StructureQuery",
                 RTType = "Recall")

drawsStructureSequence <- posterior_samples(m)$b_StructureSequence #effect of number of bars on encoding RT
drawsStructureSequenceRecall <- drawsStructureSequence + 
  posterior_samples(m)$`b_StructureSequence:Stimulus_typequery`

df5 = data.frame(draws = drawsStructureSequence,
                 predictor = "StructureSequence",
                 RTType = "Encoding")
df6 = data.frame(draws = drawsStructureSequenceRecall,
                 predictor = "StructureSequence",
                 RTType = "Recall")

drawsConditionSort <- posterior_samples(m)$b_ConditionSort #effect of number of bars on encoding RT
drawsConditionSortRecall <- drawsConditionSort + 
  posterior_samples(m)$`b_ConditionSort:Stimulus_typequery`

df7 = data.frame(draws = drawsConditionSort,
                 predictor = "ConditionSort",
                 RTType = "Encoding")
df8 = data.frame(draws = drawsConditionSortRecall,
                 predictor = "ConditionSort",
                 RTType = "Recall")

# put all dataframes together for later plotting
allEffects <- rbind(df1, df2, df3, df4, df5, df6, df7, df8)
meanEffect <- ddply(allEffects, ~ predictor + RTType, summarize, meanEffect = mean(draws))
# get confidence intervals
for (predict in unique(allEffects$predictor)){
  for (RTtype in unique(allEffects$RTType)){
    a = subset(allEffects, RTType == RTtype & predictor == predict)
    print(RTtype)
    print(predict)
    print("lower bound:")
    print(round(sort(a$draws)[500],2))
    print("upper bound:")
    print(round(sort(a$draws)[19500],2))
  }
}

# let the plotting begin
EncodingRecallplot <- ggplot(allEffects, aes(x = draws, y = predictor, fill = RTType))+
  stat_halfeye()+
  labs(x = "Posterior Estimate", y = NULL,
       fill = "RT type") +
  geom_vline(xintercept = 0)+
  theme_minimal()+
  theme(text = element_text(size = textsize, family = "sans"),legend.position = "bottom")+
  scale_y_discrete(labels = c("Sort Task", "Sequence Length", "Query Structure", "Sequence Structure" ))+
  facet_grid(cols = vars(RTType))+
  theme(legend.position = "None")
EncodingRecallplot

#######################################################################
# Summed RT Model
######################################################################
setwd(here("analyses", "brmModel"))
d2 <- read.csv(here("Data", "Experiment", "data_to_work_with2.csv"))
#create data set as before
dp <- subset(d2, rt <= rt_cutoff & correct == 1 & otherRT <= rt_cutoff & Stimulus_type == "bars")
dp$summedRT <- dp$rt + dp$otherRT

# The full Encoding RT model 
mrtSummedallStructConditionNotLogFinalEx <- brm(summedRT ~ Number_of_Bars + Structure + Condition + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                                                data = dp,
                                                iter = 10000,
                                                cores = 4,
                                                save_pars = save_pars(all = TRUE),
                                                seed = seed,
                                                control = list(max_treedepth = 15, adapt_delta = 0.99),
                                                family = exgaussian(),
                                                file = "mrtSummedallStructConditionNotLogFinalEx")

summedRTPlot <- mcmc_plot(mrtSummedallStructConditionNotLogFinalEx, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("Posterior Estimate")+
  theme_minimal()+
  theme(text = element_text(size = textsize, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))# 
summedRTPlot

##################################################################
# Final Plot Recall RT sub figure
##################################################################

# adjusted figure after review
RecallRTFigure <- (EncodingRecallplot / 
                     ((RTTradeoff_plot + theme(axis.title.y = element_text(margin = margin(r = -75, unit = "pt")))) + summedRTPlot)) +
  plot_annotation(tag_levels = list(c("A", "B", "C")))
  #plot_layout(widths = c(3, 2, 2))
RecallRTFigure

# figure before review
RecallRTFigure <- (EncodingRecallplot + RTTradeoff_plot) +
  plot_annotation(tag_levels = list(c("A", "B")))+
  plot_layout(widths = c(3, 2))
RecallRTFigure


ggsave(here("Figures", "FigureA1_RecallRT.pdf"), RecallRTFigure,  
       width = 10, 
       height = 8)

#################################################################################
# Model behavioural output evaluator and mutator
#################################################################################

setwd(here("Model", "hypoFittingBrmModels"))
EffectsModelMutator <- brm(time ~ Structure + NoB + (Structure + NoB|participant),
                    data = modelData,
                    chains = 2,
                    save_pars = save_pars(all = TRUE),
                    control = list(max_treedepth = 15, adapt_delta = 0.99),
                    file = "hypoMutatorEmulatingRTresultsBucketSortIndHypoLL")

ModelofModelPlotMutator <- mcmc_plot(EffectsModelMutator, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("Posterior Estimates")+
  ggtitle("Hypotheses Mutator")+
  theme_minimal()+
  xlim(-1.5, 2)+
  theme(text = element_text(size = 15, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Query Structure", "Sequence Structure", "Sequence Length"))
ModelofModelPlotMutator


ggsave(here("Figures", "FigureB1_hypothesisMutator.pdf"), ModelofModelPlotMutator,  
       width = 6, 
       height = 4)

##############################################################
# the accuracy time trade off figure can be found in the accuracy_time_trade_offs.R script
###############################################################


