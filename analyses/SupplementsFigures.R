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

ggsave(here("Figures", "AccuracyFigure.pdf"), AccuracyFigure,  
       width = 7, 
       height = 8)

# last saved 16.09.2022 after exclusion correction

###################################################################
# Recall RT Figure
###################################################################
#####################################################
# Full Recall RT model
textsize <- 16
setwd(here("analyses", "brmModel"))
#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# The full Recall RT model
mRecallRTallStructCondition <- brm(queryRT ~ Number_of_Bars + Structure + Condition + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                                   data = dp,
                                   iter = 10000,
                                   cores = 4,
                                   save_pars = save_pars(all = TRUE),
                                   seed = seed,
                                   file = "mmRecallRTallStructConditionNotLogFinal")
# Model summary
summary(mRecallRTallStructCondition)
# results form 26.07.2022
# Intercept             0.04      0.06    -0.08     0.17 1.00     2328     5941
# Number_of_Bars        0.35      0.03     0.30     0.41 1.00     1411     2931
# StructureQuery       -0.15      0.06    -0.26    -0.04 1.00     3877     7702
# StructureSequence    -0.01      0.02    -0.06     0.03 1.00    10692    13260
# ConditionSort         0.18      0.05     0.07     0.28 1.00     2616     5700

# only plot main populationlevel effects
RecallRTmodel_plot <- mcmc_plot(mRecallRTallStructCondition, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("Posterior Weights")+
  theme_minimal()+
  theme(text = element_text(size = textsize, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+
  ggtitle("Recall RT")
RecallRTmodel_plot

###################################################################
# include recall rt as a factor in encoding RT model

#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# regression
mrtallincludingRecallRT <- brm(rt ~ Number_of_Bars + Structure + Condition + queryRT + (Number_of_Bars + Structure + Condition + queryRT + Block|subject_id), 
                               data = dp,
                               iter = 10000,
                               cores = 4,
                               save_pars = save_pars(all = TRUE),
                               seed = seed,
                               file = "mrtallincludingRecallRTNotLogFinal")
summary(mrtallincludingRecallRT)
# results 28.06.2022
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            -0.25      0.16    -0.58     0.07 1.00     2064     3572
# Number_of_Bars        0.82      0.05     0.72     0.91 1.00     1674     3394
# StructureQuery       -0.26      0.07    -0.40    -0.12 1.00     3521     7460
# StructureSequence    -0.23      0.05    -0.34    -0.12 1.00     3753     9952
# ConditionSort         0.37      0.06     0.25     0.50 1.00     4160     8932
# queryRT               0.26      0.05     0.17     0.35 1.00     5490    10482

# only plot main population level effects
model_plotwithRecallRT <- mcmc_plot(mrtallincludingRecallRT, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("Posterior Weights")+
  theme_minimal()+
  theme(text = element_text(size = textsize, family = "sans"))+
  scale_y_discrete(labels = c("Intercept",
                              "Sequence Length",
                              "Query Structure",
                              "Sequence Structure",
                              "Sort Task",
                              "Recall RT"))+
  ggtitle("Encoding RT")
model_plotwithRecallRT

#############################################################################################################################################################
# sum of recall Rt and Encoding Rt together
d$fullRT <- d$rt + d$queryRT


#create data set as before
dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)

# regression
mFullRTallStructCondition <- brm(fullRT ~ Number_of_Bars + Structure + Condition + (Number_of_Bars + Structure + Condition + Block|subject_id), 
                                 data = dp,
                                 iter = 10000,
                                 cores = 4,
                                 save_pars = save_pars(all = TRUE),
                                 seed = seed,
                                 file = "mFullRTallStructConditionNotLogFinal")

# Model summary
summary(mFullRTallStructCondition)
# results 29.06.2022
# Intercept            -0.71      0.23    -1.17    -0.25 1.00     2533     4894
# Number_of_Bars        1.37      0.08     1.21     1.54 1.00     1333     2699
# StructureQuery       -0.35      0.17    -0.69    -0.02 1.00     3022     6378
# StructureSequence    -0.23      0.09    -0.40    -0.06 1.00     5172    11796
# ConditionSort         0.85      0.14     0.58     1.12 1.00     5137     9157


SumRTModel_plot <- mcmc_plot(mFullRTallStructCondition, type = "areas", prob = 0.95, variable = "^b_", regex = TRUE)+
  geom_vline(xintercept = 0, linetype = "longdash", color = "gray")+
  xlab("Posterior Weights")+
  theme_minimal()+
  theme(text = element_text(size = textsize, family = "sans"))+
  scale_y_discrete(labels = c("Intercept", "Sequence Length", "Query Structure", "Sequence Structure", "Sort Task"))+# 
  ggtitle("Recall + Encoding RT")
SumRTModel_plot
#####################################################################
# look at trade-off between recall RT and encoding RT

dp <- subset(d, rt <= rt_cutoff & correct == 1 & queryRT <= rt_cutoff)
dpSummary <- ddply(dp, ~ subject_id + Structure + Condition + Number_of_Bars, summarize, rt = mean(rt), queryRT = mean(queryRT))
RTTradeoff_plot <- ggplot(dpSummary, aes(x = rt, y = queryRT, color = Structure, fill = Structure))+
  geom_point(alpha = 0.2)+
  scale_color_manual(values = cbbPalette)+
  geom_smooth(method = lm, se = FALSE)+
  ylab("Recall RT")+
  xlab("Encoding RT")+
  scale_x_continuous(breaks = c(5, 10))+
  scale_y_continuous(breaks = c(5, 10))+
  theme_minimal()+
  theme(text = element_text(size = textsize, family = "sans"),legend.position = "bottom")+
  facet_grid(rows = vars(Condition), cols = vars(Number_of_Bars))
RTTradeoff_plot 


RecallRTFigure <- (RecallRTmodel_plot + RTTradeoff_plot)/
  (model_plotwithRecallRT + SumRTModel_plot) +  
  plot_annotation(tag_levels = list(c("A", "B", "C")))
RecallRTFigure


ggsave(here("Figures", "RecallRTFigure.pdf"), RecallRTFigure,  
       width = 10, 
       height = 10)


##############################################################
# the accuracy time trade off figure can be found in the accuracy_time_trade_offs.R script
###############################################################
