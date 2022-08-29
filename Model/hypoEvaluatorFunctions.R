# hypothesis evaluator: all functions
# this is basically the implementation, without running the model yet


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

#################################################################################
# functions
##############################################################################
# function to generate a trial
# conditions can either manipulate task structure or query structure
# the length of each trial is randomly determined from 1-7
# only necessary for simulations, later the model is run on the actual trial data
generatetrial <- function(structure_task=FALSE, structure_query=FALSE){
  sequence_length = sample(1:7, 1)
  print(sequence_length)
  #normal condition
  if (structure_task==FALSE &  structure_query==FALSE){
    #letter are assigned to sizes randomly
    #target letter is assigned at random
    d <- data.frame(bar=sample(letters[1:sequence_length]), size=1:sequence_length, target=sample(letters[1:sequence_length], 1))
  }
  #structure in the task
  if (structure_task==TRUE &  structure_query==FALSE){
    #initialize unpermuted data frame
    d <- data.frame(bar=letters[1:sequence_length], size=1:sequence_length, target=sample(letters[1:sequence_length], 1))
    #where to insert connected items 
    insert <- sample(1:max(sequence_length-3, 1), 1)
    #scramble the other letters
    l <- sample(letters[1:4], 4)
    #put them together, e-f-g is always connected
    if (insert == 1){
      l <- c(letters[5:7], l)
    }else{
      l <- c(l[1:insert-1], letters[5:7], l[(insert):7])
    }
    #only need the first how ever long the sequence is
    d$bar <- l[1:sequence_length]
    d$target = sample(d$bar, 1)
  }
  #structure in the query
  if (structure_task==FALSE &  structure_query==TRUE){
    #initialize unpermuted frame
    d <- data.frame(bar=sample(letters[1:7]), size=sample(1:7), target=sample(letters[1:7], 1))
    #sample a target that is within the first three sorted bars
    #d$target <- sample(subset(d, size<=3)$bar, 1) this is the old version, where the first three things got sampled
    d = d[1:sequence_length,]
    d$size = sample(1:sequence_length)
    d$target <- sample(subset(d, size == max(d$size))$bar, 1) # new version where the tallest bar always gets samples
    
  }
  #randomly permute for giggles
  d <- d[sample(1:sequence_length),]
  #return task data frame
  return(d)
}


#string to vectore: "ab"->[a, b], etc. 
#this will make it easier for the sort later on
string_to_vec <- function(string){
  #initialize empty vector
  out <- numeric()
  #loop over string length
  for (i in 1:nchar(string)){
    #concatenate letters to vector
    out <- c(out, substr(string, i, i))
  }
  #return vector
  return(out)
}

# function for execution of one trial
# trial: generate trial or use trial from actual data
# threshold: when to stop sorting
# connect string indicating hypotheses of connected letters
execute <- function(trial, threshold, connect, startAtSmall){
  if(startAtSmall==1){#sort given a trial, threshold, and hypotheses
    out <- noisy_sort(trial, threshold, connect)
    #check if sorted response is correct
    correct <- ifelse(which(paste(trial$target[1])==out$collect)==subset(trial, paste(bar)==target)$size, 1, 0)
    #if not present, then it means it stopped too early, which is also incorrect
    correct <- ifelse(length(correct)==0, 0, correct)
    #return time and whether or not correct
    return(data.frame(time = out$time[1], correct = correct))
  }else{
    out <- noisy_sort_backwards(trial, threshold, connect)
    #check if sorted response is correct
    correct <- ifelse(which(paste(trial$target[1])==out$collect)==subset(trial, paste(bar)==target)$size, 1, 0)
    #if not present, then it means it stopped too early, which is also incorrect
    correct <- ifelse(length(correct)==0, 0, correct)
    #return time and whether or not correct
    return(data.frame(time = out$time[1], correct = correct))
  }
  
}

#function that replaces bottom 5% with top 5%
die_and_replace <- function(data, fitness){
  #top 5% fittest programs
  top <- (fitness %in% sort(fitness)[1:20])
  #bottom 5% fittest programs
  bottom <- (fitness %in% sort(fitness)[(length(fitness)-19):length(fitness)])
  #replace
  data[bottom==1,] <- data[top==1,]
  #return new data
  return(data)
}

########################################
# bucket sort functions for the sorter
########################################

# bucket sort algorithm should scale 2n if no structure is learned
# trial: trial data frame, generated from the trials the participants saw
# threshold: until what point to look at before stopping
# connect: hypothesis about which bars are connected
noisy_sort <- function(trial, threshold = 7, connections = 'a'){
  # create a vector from the connect string
  connections <- string_to_vec(connections)
  connectNotFound = TRUE
  tookConnectAway = FALSE
  # initialize a vector to collect the outcomes of sorting algorithm
  collect <- seq(1, length(trial$bar))
  collect[1:length(collect)] = NaN
  # length of the vector has to be bigger than 1 to make a time-saving connections
  if (length(connections)>1 & connections[1] %in% trial$bar){
    # remove the connected items that comes after the first connect from the to-be-sorted data 
    # that should reduce the data to be sorted and therefore speed up the process
    trial <- subset(trial, !(bar %in% connections[-1]))
    tookConnectAway = TRUE
  }
  # count and time are 0 at start
  count <- 0
  time <- 0
  # number of to-be-sorted rectangles
  nr <- nrow(trial)
  # go through the trials to determine the shortest bar
  MinSize = trial$size[1]
  smallest = 1
  for (i in 1:nrow(trial)){
    # check if current rectangle is smaller than the current smallest rectangle
    if (trial$size[i] < MinSize){
      MinSize = trial$size[i]
      smallest = i
    }
    # increment time
    time <- time + 1
  }
  # add the smallest bar to the list
  collect[1] <- trial$bar[smallest]
  # bar gets removed from consideration set
  trial <- trial[-smallest,]
  # count increments
  count <- count + 1
  # if the hypothesis contained at least two elements
  if (length(connections)>1 & connectNotFound){
    # if the first part was collected
    if(connections[1] %in% collect){
      connectNotFound = FALSE
      # increase the count accordingly
      count <- count + length(connections)-1
    }
  }
  # the count has to be smaller than the threshold or total number of trials
  while (count < pmin(threshold, nr) & length(trial$size) > 0){
    position  <- trial$size[1] # this is the relative size of the bar which allows it to be put in the correct bucket
    collect[position] <- trial$bar[1]
    # bar gets removed from consideration set
    trial <- trial[-1,]
    # count increments
    count <- count + 1
    time = time + 1
    if (length(connections)>1 & connectNotFound){
      # if the first part was collected
      if(connections[1] %in% collect){
        connectNotFound = FALSE
        # increase the count accordingly
        count <- count + length(connections)-1
      }
    }
  }
  
  
  # add the connection if it was found
  if (connectNotFound == FALSE){
    # if the first part was collected
    if(connections[1] %in% collect){
      # mark where it appeared  
      mark <- which(collect == connections[1])
      # bind it together
      if (mark == length(collect)){
        collect <- c(collect, connections[-1])
      }else{
        for (connect in range(1, (length(connections)-1))){
          if (is.nan(as.numeric(collect[(mark+connect)]))){
            collect[(mark+connect)] = connections[(connect+1)]
          } else {
            collect[(mark+connect):(length(collect)+1)] <- c(connections[(connect+1)], collect[(mark+connect):length(collect)])
          }
        }
      }
    }
  }
  
  # attach the rest of the trial in the unsorted order:
  if (connectNotFound & tookConnectAway){
    toInsert <- c(trial$bar, connections[-1])
  } else {
    toInsert <- c(trial$bar)
  }
  foundBar = 1
  for (i in 1:length(collect)){
    if (foundBar == length(toInsert)+1){
      break
    }
    if(is.nan(as.numeric(collect[i]))){
      collect[i] = toInsert[foundBar]
      foundBar = foundBar + 1
    }
  }
  
  
  collect <- collect[!is.nan(as.numeric(collect))]

  
  # collect into data frame
  dout <- data.frame(collect, time)
  # return data frame
  return(dout)
}


noisy_sort_backwards <- function(trial, threshold = 7, connections = 'a'){
  # create a vector from the connect string
  connections <- string_to_vec(connections)
  connectNotFound = TRUE
  tookConnectAway = FALSE
  # initialize a vector to collect the outcomes of sorting algorithm
  collect <- seq(1, length(trial$bar))
  collect[1:length(collect)] = NaN
  # length of the vector has to be bigger than 1 to make a time-saving connections
  if (length(connections)>1 & connections[1] %in% trial$bar){
    # remove the connected items that comes after the first connect from the to-be-sorted data 
    # that should reduce the data to be sorted and therefore speed up the process
    trial <- subset(trial, !(bar %in% connections[-1]))
    tookConnectAway = TRUE
  }
  # count and time are 0 at start
  count <- 0
  time <- 0
  # number of to-be-sorted rectangles
  nr <- nrow(trial)
  # go through the trials to determine the shortest bar
  MaxSize = trial$size[1]
  tallest = 1
  for (i in 1:nrow(trial)){
    # check if current rectangle is smaller than the current tallest rectangle
    if (trial$size[i] > MaxSize){
      MaxSize = trial$size[i]
      tallest = i
    }
    # increment time
    time <- time + 1
  }
  # add the tallest bar to the list
  collect[length(collect)] <- trial$bar[tallest]
  # bar gets removed from consideration set
  trial <- trial[-tallest,]
  # count increments
  count <- count + 1
  # if the hypothesis contained at least two elements
  if (length(connections)>1 & connectNotFound){
    # if the first part was collected
    if(connections[1] %in% collect){
      connectNotFound = FALSE
      # increase the count accordingly
      count <- count + length(connections)-1
    }
  }
  # the count has to be smaller than the threshold or total number of trials
  while (count < pmin(threshold, nr) & length(trial$size) > 0){
    position  <- trial$size[1] # this is the relative size of the bar which allows it to be put in the correct bucket
    collect[position] <- trial$bar[1]
    # bar gets removed from consideration set
    trial <- trial[-1,]
    # count increments
    count <- count + 1
    time = time + 1
    if (length(connections)>1 & connectNotFound){
      # if the first part was collected
      if(connections[1] %in% collect){
        connectNotFound = FALSE
        # increase the count accordingly
        count <- count + length(connections)-1
      }
    }
  }
  
  
  # add the connection if it was found
  if (connectNotFound == FALSE){
    # if the first part was collected
    if(connections[1] %in% collect){
      # mark where it appeared  
      mark <- which(collect == connections[1])
      # bind it together
      if (mark == length(collect)){
        collect <- c(collect, connections[-1])
      }else{
        for (connect in range(1, (length(connections)-1))){
          if (is.nan(as.numeric(collect[(mark+connect)]))){
            collect[(mark+connect)] = connections[(connect+1)]
          } else {
            collect[(mark+connect):(length(collect)+1)] <- c(connections[(connect+1)], collect[(mark+connect):length(collect)])
          }
        }
      }
    }
  }
  
  # attach the rest of the trial in the unsorted order:
  if (connectNotFound & tookConnectAway){
    toInsert <- c(trial$bar, connections[-1])
  } else {
    toInsert <- c(trial$bar)
  }
  foundBar = 1
  for (i in 1:length(collect)){
    if (foundBar == length(toInsert)+1){
      break
    }
    if(is.nan(as.numeric(collect[i]))){
      collect[i] = toInsert[foundBar]
      foundBar = foundBar + 1
    }
  }
  
  
  collect <- collect[!is.nan(as.numeric(collect))]

  
  
  # collect into data frame
  dout <- data.frame(collect, time)
  # return data frame
  return(dout)
}




###########################
# function to run the simulation
############################
runSimulation = function(start_particles, # generated particles
                         total_runs, # int
                         trials, # int
                         particles_to_eliminate, # float, percentage of particles that get eliminated each trial
                         structure_query = FALSE, 
                         structure_task = FALSE){
  # matrix to collect run time
  times = matrix(0, trials, total_runs)
  accuracy = matrix(0, trials, total_runs)
  NoB = matrix(0, trials, total_runs) # Number of Bars
  finalHypotheses = data.frame(threshold = as.numeric(),
                               conncetions = as.character(),
                               startAtSmall = as.integer(),
                               time = as.numeric(),
                               run = as.numeric())
  # determine the number of hypotheses 
  nhypotheses = length(start_particles$threshold)
  print("number of hypotheses")
  print(nhypotheses)
  # for 20 simulations in total
  for (n_runs in 1:total_runs){
    print("started run:")
    print(n_runs)
    # generate a trial
    if(structure_query){
      trial = generatetrial(structure_query = TRUE)
    }else if(structure_task){
      trial = generatetrial(structure_task = TRUE)
    }else{
      trial = generatetrial()
    }
    
    # create an out file
    # bind the output lists into data frame
    dout <- do.call(rbind, 
                    #apply trial to each program defined by the particles
                    apply(start_particles, 1, function(x) {execute(trial = trial, threshold = x[1], connect=x[2], startAtSmall = x[3])})
    )
    # bind them togethers
    hypotheses <- cbind(start_particles, dout)
    
    # go through 50 trials
    for (i in 1:trials){
      print("trial no:")
      print(i)
      print("n of hypotheses before selection")
      print(length(hypotheses$threshold))
      nhypotheses = length(hypotheses$threshold)
      # generate new trial
      # generate a trial
      if(structure_query){
        trial = generatetrial(structure_query = TRUE)
      }else if(structure_task){
        trial = generatetrial(structure_task = TRUE)
      }else{
        trial = generatetrial()
      }
      
      NoB[i,n_runs] <- length(trial$bar)
      
      print("current trial")
      print(trial)
      
      # collect results of all programs
      dout <- do.call(rbind, apply(hypotheses[,1:3], 1, function(x) {execute(trial = trial, threshold = x[1], connect = x[2], startAtSmall = x[3])}))
      print("dout")
      print(dout)
      # update time so that each hypotheses time reflects the average time
      hypotheses$time <- hypotheses$time*((i-1)/i) + (1/i)*dout$time
      # update accuracy
      hypotheses$correct <- dout$correct
      # only keep correct hypotheses
      incorrect = length(subset(hypotheses, correct == 0)$correct)
      hypotheses = subset(hypotheses, correct == 1)
      neliminate = ceiling(particles_to_eliminate*nhypotheses)
      if (incorrect < neliminate & length(hypotheses$correct) > 1){
        remaining = min(neliminate - incorrect, length(hypotheses$correct) - 1)
        count = 0
        while (count < remaining){
          print("slowest")
          print(subset(hypotheses, time == max(hypotheses$time)))
          count = count + length(subset(hypotheses, time == max(hypotheses$time)))
          hypotheses = subset(hypotheses, time != max(hypotheses$time))
        }
      }else{
        print("Not eliminating more, because I had")
        print(neliminate)
        print(incorrect)
      }
      
      
      # track max time
      times[i, n_runs] <- max(dout$time)
      # track mean accuracy
      accuracy[i,n_runs] <- mean(dout$correct)
      
    }
    print("run:")
    print(n_runs)
    print("final hypothese")
    print(hypotheses)
    dummy = hypotheses
    dummy$run = n_runs
    finalHypotheses = rbind(finalHypotheses, dummy)
  }
  output <- list()
  output[[1]] <- times
  output[[2]] <- accuracy
  output[[3]] <- finalHypotheses
  output[[4]] <- NoB
  return(output)
}

###################################################################
# function to run model on participant data. 
# should take dataframe as input and should output:
# 1 the estimated rt
# the real rt
# accuracy???? 
# final hypotheses?
################################################################
runHypoSortOnData = function(start_particles, # generated particles
                             data,
                             particles_to_eliminate # float, percentage of particles that get eliminated each trial
){
  trials = 35
  n_participants = length(unique(data$subject_id))
  # matrix to collect run time
  times = matrix(0, trials, n_participants)
  accuracy = matrix(0, trials, n_participants)
  # participant data
  realRTs = matrix(0, trials, n_participants)
  realAccuracy = matrix(0, trials, n_participants)
  NoB = matrix(0, trials, n_participants) # Number of Bars
  correctConnection = matrix(0, trials, n_participants)
  # save the correct connection
  if (data$Structure[1] == "Sequence"){
    for (i in 1:length(unique(data$subject_id))){
      theConnection = unlist(lapply(strsplit(subset(data, Number_of_Bars == 3 & subject_id == unique(data$subject_id)[i])$color_sequence[1], ",")[[1]], convert_colors))
      correctConnection[,i] = paste(theConnection, collapse = "")
    }
  }
  finalHypotheses = data.frame(threshold = as.numeric(),
                               conncetions = as.character(),
                               startAtSmall = as.integer(),
                               time = as.numeric(),
                               run = as.numeric())
  # determine the number of hypotheses 
  nhypotheses = length(start_particles$threshold)
  print("number of hypotheses")
  print(nhypotheses)
  # for 20 simulations in total
  for (participant in 1:n_participants){
    print("started run:")
    print(participant)
    # get current trial
    current_trial = subset(data, subject_id == unique(data$subject_id)[participant])[1,]
    trial = datacolum_to_trial(current_trial)
    
    # create an out file
    # bind the output lists into data frame
    dout <- do.call(rbind, 
                    # apply trial to each program defined by the particles
                    apply(start_particles, 1, function(x) {execute(trial = trial, threshold = x[1], connect=x[2], startAtSmall = x[3])})
    )
    # bind them togethers
    hypotheses <- cbind(start_particles, dout)
    
    # go through 50 trials
    for (i in 1:trials){
      print("trial no:")
      print(i)
      print("n of hypotheses before selection")
      print(length(hypotheses$threshold))
      nhypotheses = length(hypotheses$threshold)
      # get current trial
      current_trial = subset(data, subject_id == unique(data$subject_id)[participant])[i,]
      trial = datacolum_to_trial(current_trial)
      
      NoB[i,participant] <- length(trial$bar)
      
      print("current trial")
      print(trial)
      
      # collect results of all programs
      dout <- do.call(rbind, apply(hypotheses[,1:3], 1, function(x) {execute(trial = trial, threshold = x[1], connect = x[2], startAtSmall = x[3])}))
      print("dout")
      print(dout)
      # update time so that each hypotheses time reflects the average time
      hypotheses$time <- hypotheses$time*((i-1)/i) + (1/i)*dout$time
      # update accuracy
      hypotheses$correct <- dout$correct
      # only keep correct hypotheses
      incorrect = length(subset(hypotheses, correct == 0)$correct)
      if(incorrect < length(hypotheses$time)){
        hypotheses = subset(hypotheses, correct == 1)
      }else{
        incorrect = 0
      }
      neliminate = ceiling(particles_to_eliminate*nhypotheses)
      if (incorrect < neliminate & length(hypotheses$correct) > 1){
        remaining = min(neliminate - incorrect, length(hypotheses$correct) - 1)
        count = 0
        while (count < remaining & length(subset(hypotheses, time == max(hypotheses$time))$time) < length(hypotheses$time)){
          print("slowest")
          print(subset(hypotheses, time == max(hypotheses$time)))
          count = count + length(subset(hypotheses, time == max(hypotheses$time))$time)
          hypotheses = subset(hypotheses, time != max(hypotheses$time))
        }
      }else{
        print("Not eliminating more, because I had")
        print(neliminate)
        print(incorrect)
      }
      
      
      # track max time
      times[i, participant] <- max(dout$time)
      # track mean accuracy
      accuracy[i, participant] <- mean(dout$correct)
      # track real RT
      realRTs[i, participant] <- current_trial$RTSort_Memory
      # track real accuracy
      realAccuracy[i, participant] <- current_trial$correct
      
      
    }
    print("run:")
    print(participant)
    print("final hypothese")
    print(hypotheses)
    dummy = hypotheses
    dummy$run = participant
    finalHypotheses = rbind(finalHypotheses, dummy)
  }
  output <- list()
  output[[1]] <- times
  output[[2]] <- realRTs
  output[[3]] <- accuracy
  output[[4]] <- realAccuracy
  output[[5]] <- finalHypotheses
  output[[6]] <- NoB
  output[[7]] <- correctConnection
  return(output)
}

####################################################################
# function that runs the simulation on all the data of all conditions
####################################################################
runHypoSortOnAllData = function(start_particles,
                                AllData,
                                particles_to_eliminate){
  ###########################################
  # query Structure
  ############################################
  dQuery <- subset(AllData, Structure == "Query" & Condition == "Sort" & Stimulus_type == "bars")
  
  QueryDataOutput <- runHypoSortOnData(start_particles,
                                       data = dQuery,
                                       particles_to_eliminate)
  
  Querytimes <- QueryDataOutput[[1]]
  QueryrealRTs <- QueryDataOutput[[2]] 
  Queryaccuracy <- QueryDataOutput[[3]] 
  QueryrealAccuracy <- QueryDataOutput[[4]] 
  QueryfinalHypotheses <- QueryDataOutput[[5]] 
  QueryNoB <- QueryDataOutput[[6]]
  QuerytrueConnection <- QueryDataOutput[[7]]
  
  QueryfinalHypotheses$Structure = "Query"
  QueryfinalHypotheses$trueConnection = 0
  
  ModelToData = data.frame(NoB = as.numeric(),
                           time = as.numeric(),
                           realRT = as.numeric(),
                           accurcay = as.integer(),
                           realAccuracy = as.integer(),
                           trueConnection = as.character(),
                           participant = as.numeric(),
                           Structure = as.character(), 
                           trial = as.numeric())
  
  for (i in 1:length(unique(dQuery$subject_id))){
    ModelToData = rbind(ModelToData, data.frame(NoB = QueryNoB[,i],
                                                time = Querytimes[,i],
                                                realRT = QueryrealRTs[,i],
                                                accuracy = Queryaccuracy[,i],
                                                realAccuracy = QueryrealAccuracy[,i],
                                                trueConnection = "a",
                                                participant = unique(dQuery$subject_id)[i],
                                                Structure = "Query", 
                                                trial = 1:35))
  }
  ###########################################
  # Sequence Structure
  ############################################
  dSequence <- subset(AllData, Structure == "Sequence" & Condition == "Sort" & Stimulus_type == "bars")
  
  SequenceDataOutput <- runHypoSortOnData(start_particles,
                                          data = dSequence,
                                          particles_to_eliminate)
  
  Sequencetimes <- SequenceDataOutput[[1]]
  SequencerealRTs <- SequenceDataOutput[[2]] 
  Sequenceaccuracy <- SequenceDataOutput[[3]] 
  SequencerealAccuracy <- SequenceDataOutput[[4]] 
  SequencefinalHypotheses <- SequenceDataOutput[[5]] 
  SequenceNoB <- SequenceDataOutput[[6]]
  SequencetrueConnection <- SequenceDataOutput[[7]]
  
  SequencefinalHypotheses$Structure = "Sequence"
  SequencefinalHypotheses$trueConnection = "a"#SequencetrueConnection[1,]
  
  counter = 1
  for (i in 1: length(SequencetrueConnection[1,])){
    while(SequencefinalHypotheses$run[counter] == i & counter <= length(SequencefinalHypotheses$run)){
      SequencefinalHypotheses$trueConnection[counter] = SequencetrueConnection[1,i]
      counter = counter + 1
    }
  }
  
  for (i in 1:length(unique(dSequence$subject_id))){
    ModelToData = rbind(ModelToData, data.frame(NoB = SequenceNoB[,i],
                                                time = Sequencetimes[,i],
                                                realRT = SequencerealRTs[,i],
                                                accuracy = Sequenceaccuracy[,i],
                                                realAccuracy = SequencerealAccuracy[,i],
                                                trueConnection = SequencetrueConnection[1,i],
                                                participant = unique(dSequence$subject_id)[i],
                                                Structure = "Sequence", 
                                                trial = 1:35))
  }
  ###########################################
  # None Structure
  ############################################
  dNone <- subset(AllData, Structure == "None" & Condition == "Sort" & Stimulus_type == "bars")
  
  NoneDataOutput <- runHypoSortOnData(start_particles,
                                      data = dNone,
                                      particles_to_eliminate)
  
  Nonetimes <- NoneDataOutput[[1]]
  NonerealRTs <- NoneDataOutput[[2]] 
  Noneaccuracy <- NoneDataOutput[[3]] 
  NonerealAccuracy <- NoneDataOutput[[4]] 
  NonefinalHypotheses <- NoneDataOutput[[5]] 
  NoneNoB <- NoneDataOutput[[6]]
  NonetrueConnection <- NoneDataOutput[[7]]
  
  NonefinalHypotheses$Structure = "None"
  NonefinalHypotheses$trueConnection = 0
  
  for (i in 1:length(unique(dNone$subject_id))){
    ModelToData = rbind(ModelToData, data.frame(NoB = NoneNoB[,i],
                                                time = Nonetimes[,i],
                                                realRT = NonerealRTs[,i],
                                                accuracy = Noneaccuracy[,i],
                                                realAccuracy = NonerealAccuracy[,i],
                                                trueConnection = "a",
                                                participant = unique(dNone$subject_id)[i],
                                                Structure = "None", 
                                                trial = 1:35))
  }
  allfinalHypo = rbind(NonefinalHypotheses, QueryfinalHypotheses, SequencefinalHypotheses)
  output <- list()
  output[[1]] <- ModelToData
  output[[2]] <- allfinalHypo
  return(output)
}

#################################################
# necessary extra functions for working with the actual data
###############################################
convert_colors <- function(color){
  colors = c("green", "blue", "red", "purple", "saddlebrown", "white", "orange", "black", "yellow", "pink")
  myletters = c(letters[1:length(colors)])
  for (i in 1:length(colors)){
    if (colors[i] == color){
      return(myletters[i])
    }
  }
  print("no match found")
}

convert_postion_seq_to_height_order <- function(position_sequence){
  position_sequence <- lapply(strsplit(position_sequence, ","), as.numeric)[[1]]
  # the following is only necessary, because apperently somtiems the positions where wrongly saved (e.g. 390 as 39)
  for (i in 1:length(position_sequence)){
    while(position_sequence[i] < 100){
      position_sequence[i] <- position_sequence[i]*10
    } 
  }
  postions_in_order <- sort.int(position_sequence)
  positions <- 1:length(position_sequence)
  for (i in 1:length(position_sequence)){
    positions[i] = which(position_sequence == postions_in_order[i])
  }
  return(positions)
}

datacolum_to_trial <- function(trialData){  # d[i,]
  heights <- convert_postion_seq_to_height_order(trialData$position_sequence)
  target <- convert_colors(trialData$Query_color)
  bars <- unlist(lapply(strsplit(trialData$color_sequence, ",")[[1]][heights], convert_colors))
  trial <- data.frame(bar = bars, size = heights, target = target)
  return(trial)
}

####################################################################

#standard error
se <- function(x){sd(x)/sqrt(length(x))}