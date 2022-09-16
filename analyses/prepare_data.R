# This script cleans all the data and prepares it for further analysis
# by creating a cleaned data file

# House keeping
rm(list = ls())

#########################################################################
# load packages
#########################################################################

library(readr)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)
library(gridExtra)
library(ggthemes)
library(ggsignif)
library(brms)
library(here)

# SET WD:
# If you are not using the here package and R projects, this is the place to 
# set the working directory to the inside of the github project
# otherwise:
setwd(here("Data","Experiment"))

###############################################################################
# set important variables
###########################################################################
# minimum mean performance required to be included in the study
cutoff <- 0.75

#########################################################################
# load data and exclude not independent data and people who do not meet the performance threshold
#########################################################################


#read in data
d1<-read.csv("batch_1_and_2.csv", sep=';')
d2<-read.csv("batch_3.csv", sep=';')
d3<-read.csv("batch_4.csv", sep=';')
d4<-read.csv("batch_5.csv", sep=';')
d5<-read.csv("batch_6.csv", sep=';')
d6<-read.csv("batch_7.csv", sep=';')
d7<-read.csv("batch_8_and_9.csv", sep=';')
d8<-read.csv("batch_10.csv", sep=';')
d9<-read.csv("batch_11_12.csv", sep=';')
d10<-read.csv("batch_13.csv", sep=';')

#bind together
d<-rbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10)

# exlusion list
d <-subset(d, subject_id != "gjfcv4xzcw54pq6") # participated in pilot and batch 2 and 4
d <-subset(d, subject_id != "q0g32jefcy9uusg") # same as above
d <-subset(d, subject_id != "o6w9fl1y5jgmx4o") # participated twice
d <-subset(d, subject_id != "ddd9z40nfc89n9e") # participated twice

# getting the trials, which contain the length of the experiment
times <- subset(d,Stimulus_type == "end_of_experiment")
times$timeMs <- as.numeric(as.character(times$timeMs))

# calculating the mean time in mins
mean(times$timeMs)/60000

# histogram of the time in mins
hist(times$timeMs/60000)

# only select the experimental data (i.e. throw out the pilot data and between trials)
d <- d[ d$TrialType == "experimental_data", ]
d <- d[ d$trial_type == "psychophysics", ]

# number of participants = 93

##################################################################################
# do the exclusion based on the cutoff
##############################################################################

#recode the correct label for easier processing
d$correct <- ifelse(d$correct == "true", 1, 0)
#calculate the average accuracy per block and exclude everyone who has a block that is below the cutoff
did <- ddply(subset(d, Block != 'NULL' & Block != "0"),  ~ subject_id + Block, summarize, m = mean(correct))
#find the partcicipants below the cutoff
did <- did[did$m < cutoff,]

#get their id (its 20 participants)
throwout <- paste(unique(did$subject_id))

#throw them out (remaining: 73)
d <- subset(d, !(subject_id %in% throwout))

# number of participants after exlusion = 73



################################################################################
# cleaning: recode some of the variables for easy processing
################################################################################

#mark control condition
d$Condition <- ifelse(grepl('Cont', d$Conditions, fixed = TRUE), "Memory", "Sort")

#Add column for the structure manipulation: None vs Query vs Sequence
d$Structure <- 'None' #defaulting to None
sequenceStructureConditions <- c('sortStructureExp', 'sortStructureCont') 
d[d$Conditions %in% sequenceStructureConditions, 'Structure'] <- 'Sequence' 
queryStructureConditions <- c('queryStructureExp', 'queryStructureCont') 
d[d$Conditions %in% queryStructureConditions, 'Structure'] <- 'Query' 
d$Structure <- factor(d$Structure) # transforms it into a level

# make factors where necessary
d$Stimulus_type <- factor(d$Stimulus_type)
d$Query_color <- factor(d$Query_color)
d$subject_id <- factor(d$subject_id)

# unlevel
d$TrialID <- as.numeric(as.character(d$TrialID))
d$Correct_Position <- as.numeric(as.character(d$Correct_Position))
d$response_position <- as.numeric(as.character(d$response_position))
d$Number_of_Bars <- as.numeric(as.character(d$Number_of_Bars))
d$Trial_index <- as.numeric(as.character(d$Trial_index))
d$Block <- as.numeric(as.character(d$Block))

# add column that shows condition plus order
d$CondOrder <- 'none'
for (i in 1:nrow(d)){
  d[i,]$CondOrder = paste(d[i,]$Conditions, "-", as.character(d[i,]$Block), sep = "")
}

################################################################################
# change RT to seconds
#########################################################################
d$rt <- d$rt/1000 

#################################################################################
# calculate the scaling times (sort vs. memory)
#########################################################################

n_participants <- length(unique(d$subject_id))
blocks <- 6
rts <- 35*2

# Add column for the scaling times
#subtracting the rt of the matching condition
d$RTSort_Memory <- 'None' 
# subtracting the rt of the none condition
d$RTSort_MemoryNone <- 'None' 
# add a column that codes if the matching trial was solved correctly
d$MatchCorrectMemory <- 0
d$MatchRTMemory <- 0

for (participant in 1:n_participants){
  for (block in 1:blocks){
    for (rt in 1:rts ){
      current_position = (participant-1)*(blocks*rts) + (block-1)*rts + rt
      trial_id = d$TrialID[current_position]
      structure = d$Structure[current_position]
      condition = d$Condition[current_position]
      stimulus = d$Stimulus_type[current_position]
      for (block_2 in 1:blocks){
        for (rt_2 in 1:rts ){
          current_position_2 = (participant-1)*(blocks*rts) + (block_2-1)*rts + rt_2
          trial_id_2 = d$TrialID[current_position_2]
          structure_2 = d$Structure[current_position_2]
          condition_2 = d$Condition[current_position_2]
          stimulus_2 = d$Stimulus_type[current_position_2]
          correct_2 = d$correct[current_position_2]
          if (current_position == current_position_2){
            # do nothing of significance, since all you want is to find the correct matchuing trial and this would be the trial itself
          }else{
            if (trial_id == trial_id_2){
              if(structure == structure_2){
                if(condition != condition_2){
                  if(stimulus == stimulus_2){
                    if (condition == "Sort"){
                      d$RTSort_Memory[current_position] = d$rt[current_position] - d$rt[current_position_2]
                    }else{
                      d$RTSort_Memory[current_position] = d$rt[current_position_2] - d$rt[current_position]
                    }
                    # save the rt of the match, this is needed so we can di the correct exclusions later
                    d$MatchRTMemory[current_position] = d$rt[current_position_2]
                    if(correct_2 == 1){
                      d$MatchCorrectMemory[current_position] = 1
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

# add a column that codes if the matching trial was solved correctly
d$MatchCorrectMemoryNone <- 0
d$MatchRTMemoryNone <- 0

for (participant in 1:n_participants){
  for (block in 1:blocks){
    for (rt in 1:rts ){
      current_position = (participant-1)*(blocks*rts) + (block-1)*rts + rt
      trial_id = d$TrialID[current_position]
      structure = d$Structure[current_position]
      condition = d$Condition[current_position]
      stimulus = d$Stimulus_type[current_position]
      for (block_2 in 1:blocks){
        for (rt_2 in 1:rts ){
          current_position_2 = (participant-1)*(blocks*rts) + (block_2-1)*rts + rt_2
          trial_id_2 = d$TrialID[current_position_2]
          structure_2 = d$Structure[current_position_2]
          condition_2 = d$Condition[current_position_2]
          stimulus_2 = d$Stimulus_type[current_position_2]
          correct_2 = d$correct[current_position_2]
          if (current_position == current_position_2){
            # do nothing of significance, since all you want is to find the correct matchuing trial and this would be the trial itself
            #print("we are at trial:")
            #print(current_position)
          }else{
            if (trial_id == trial_id_2){
              if(structure_2 == "None"){    # this is what differ to the code above
                if(condition != condition_2){
                  if(stimulus == stimulus_2){
                    #print("success")
                    if (condition == "Sort"){
                      d$RTSort_MemoryNone[current_position] = d$rt[current_position] - d$rt[current_position_2]
                    }else{
                      d$RTSort_MemoryNone[current_position] = d$rt[current_position_2] - d$rt[current_position]
                    }
                    # save the rt of the match, this is needed so we can di the correct exclusions later
                    d$MatchRTMemoryNone[current_position] = d$rt[current_position_2]
                    if(correct_2 == 1){
                      d$MatchCorrectMemoryNone[current_position] = 1
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

d$RTSort_Memory <- as.numeric(as.character(d$RTSort_Memory))
d$RTSort_MemoryNone <- as.numeric(as.character(d$RTSort_MemoryNone))


#################################################################################
# add the query RT as a variable to the bars trial and only save those trials
##############################################################################
dbars <- subset(d, Stimulus_type == "bars")
dquery <- subset(d, Stimulus_type == "query")
dbars$queryRT <- dquery$rt
d1 <- dbars

# drop irrelevant columns
drops <- c("timeMs", "timeMin", "score", "reward", "center_x", "center_y", "avg_frame_time", "trial_type", "TrialType")
d1 <- d1[ , !(names(d1) %in% drops)]

################################################################
# save the data frame with all the new variables
##############################################################

write.csv(x = d1, file = "data_to_work_with.csv")

################################################################################
# make a second frame for interaction analysis
###############################################################################

dbars <- subset(d, Stimulus_type == "bars")
dquery <- subset(d, Stimulus_type == "query")
dbars$otherRT <- dquery$rt
dquery$otherRT <- dbars$rt
d2 <- rbind(dbars,dquery)

# drop irrelevant columns
d2 <- d2[ , !(names(d2) %in% drops)]

################################################################
# save the data frame with all the new variables
##############################################################

write.csv(x = d2, file = "data_to_work_with2.csv")
