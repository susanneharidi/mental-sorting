# Hypothesis Generatir: all functions
# this is basically the implementation of the model without running it on the data yet
# for bucket sort
# both for simulations and for run on actual data


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
library(stringr)

#################################################################################
# functions
##############################################################################

# function that gets all transitions of a sequence
get_transitions = function(sortedSequence){ # the sorted sequence needs to be a list of characters
  transitions = list()
  NofTransitions = 0
  while (length(sortedSequence) > 1){
    NofTransitions = NofTransitions +1
    transitions[NofTransitions] = list(sortedSequence[1:2])
    sortedSequence = sortedSequence[-1]
  }
  return(transitions)
}

# function that converts letters into the number that corresponds to their position in the alphabet
convert_letter_to_num = function(letter){
  myLetters = letters[1:26]
  number = match(letter, myLetters)
  return(number)
}

# finally the transition matrix update
transitionmatrix_update = function(transitionmatrix, 
                                   sortedSequence,
                                   learningrate = 1,
                                   correct = 1){
  # make sure the sorted sequence is actually long enough to provide transitions
  if (length(sortedSequence) <= 1){
    print(" the sorted sequence is to short to update the transition matrix")
    return(transitionmatrix)
  }
  NumberOfRows = dim(transitionmatrix)[1]
  NumberOfColumns = dim(transitionmatrix)[2]
  transitions = get_transitions(sortedSequence)
  # if the sorter provided an incorrect response, increase uncertainty about the transitions
  # by adding the learning-rate to all transitions and normalizing (if this is done infinitly many times the 
  # matrix should return to the uniform prior)
  if (correct == 0){
    transitionmatrix = transitionmatrix + learningrate
    for (i in 1:NumberOfRows){
      transitionmatrix[i,] = transitionmatrix[i,]/sum(transitionmatrix[i,]) # normalize row, so the sum of all probabilities is 1
    }
    return(transitionmatrix)
  }else{
    # if the sorter provided the correct response, update the transition matrix based on the observed and sorted sequence
    for(i in 1:length(transitions)){
      # update correct cell
      transitionmatrix[convert_letter_to_num(transitions[[i]][1]), convert_letter_to_num(transitions[[i]][2])] = 
        learningrate + transitionmatrix[convert_letter_to_num(transitions[[i]][1]), convert_letter_to_num(transitions[[i]][2])]
      # normalize row
      transitionmatrix[convert_letter_to_num(transitions[[i]][1]),] = 
        transitionmatrix[convert_letter_to_num(transitions[[i]][1]),]/sum(transitionmatrix[convert_letter_to_num(transitions[[i]][1]),])
    }
    return(transitionmatrix)
  }
}

hypothesis_transitions = function(indeces_list){
  # check if there are any potentially longer lists
  dublicates = anyDuplicated(unlist(indeces_list))
  connections = c()
  if (length(indeces_list) > 0){
    for ( i in 1:length(indeces_list)){
      connections = c(connections, str_flatten(letters[indeces_list[[i]]]))
    }
  }
  if (dublicates > 0){
    for (i in 1:length(indeces_list)){
      check = 1:length(indeces_list)
      ckeck = check[-1]
      for (tocheck in 1:length(check)){
        # this only find connections of length three at most
        if (indeces_list[[i]][length(indeces_list[[i]])] == indeces_list[[tocheck]][1]){
          # great, a longer connection was found. remove the shorter connections and
          connections = connections[connections != str_flatten(letters[indeces_list[[i]]])]
          connections = connections[connections != str_flatten(letters[indeces_list[[tocheck]]])]
          # add the longer one
          connections = c(connections, paste(str_flatten(letters[indeces_list[[i]]]), 
                                             str_flatten(letters[indeces_list[[tocheck]][2:length(indeces_list[[tocheck]])]]), 
                                             sep = ""))
        }
      }
    }
  }
  
  return(connections)
}

# function that checks the transition matrix for any values above the having it as a hypothesis point and 
# forms the corresponding connections
# with lists
get_hypothesis_transitions = function(transitionmatrix){
  indeces =  which(transitionmatrix > 0.8, arr.ind = TRUE)
  # if there are no transition probabilities that exceed the threshold, return "a" as default transition
  if (length(indeces) == 0){
    return(c("a")) 
  }
  colnames(indeces) <- NULL
  # turn indeces into a list
  indeces_list = list()
  indeces_list = append(indeces_list, 1:dim(indeces)[1])
  for (i in 1:dim(indeces)[1]){
    indeces_list[[i]] = indeces[i,]
  }
  #check if there are any potentially longer lists
  hypothesis_transitions(indeces_list)
}

# function to generate a trial (for simulations)
# conditions can either manipulate task structure or query structure
# the length of each trial is randomly determined from 1-7
generatetrial <- function(structure_task=FALSE, structure_query=FALSE){
  sequence_length = sample(1:7, 1)
  print(sequence_length)
  # normal condition
  if (structure_task==FALSE &  structure_query==FALSE){
    # letter are assigned to sizes randomly
    # target letter is assigned at random
    d <- data.frame(bar=sample(letters[1:sequence_length]), size=1:sequence_length, target=sample(letters[1:sequence_length], 1))
  }
  # structure in the task
  if (structure_task==TRUE &  structure_query==FALSE){
    # initialize unpermuted data frame
    d <- data.frame(bar=letters[1:sequence_length], size=1:sequence_length, target=sample(letters[1:sequence_length], 1))
    # where to insert connected items 
    insert <- sample(1:max(sequence_length-3, 1), 1)
    # scramble the other letters
    l <- sample(letters[1:4], 4)
    # put them together, e-f-g is always connected
    if (insert == 1){
      l <- c(letters[5:7], l)
    }else{
      l <- c(l[1:insert-1], letters[5:7], l[(insert):7])
    }
    # only need the first how ever long the sequence is
    d$bar <- l[1:sequence_length]
    d$target = sample(d$bar, 1)
  }
  # structure in the query
  if (structure_task==FALSE &  structure_query==TRUE){
    # initialize unpermuted frame
    d <- data.frame(bar=sample(letters[1:7]), size=sample(1:7), target=sample(letters[1:7], 1))
    # sample a target that is within the first three sorted bars
    d = d[1:sequence_length,]
    d$size = sample(1:sequence_length)
    d$target <- sample(subset(d, size == max(d$size))$bar, 1) # new version where the tallest bar always gets samples
    
  }
  # randomly permute for giggles
  d <- d[sample(1:sequence_length),]
  # return task data frame
  return(d)
}


# string to vectore: "ab"->[a, b], etc. 
# this will make it easier for the sort later on
string_to_vec <- function(string){
  # initialize empty vector
  out <- numeric()
  # loop over string length
  for (i in 1:nchar(string)){
    # concatenate letters to vector
    out <- c(out, substr(string, i, i))
  }
  # return vector
  return(out)
}

# function for execution of one trial
# trial: generate trial
# threshold: when to stop sorting
# connect string indicating hypotheses of connected letters
execute <- function(trial, threshold, connections, startAtSmall){
  if(startAtSmall==1){#sort given a trial, threshold, and hypotheses
    out <- noisy_sort(trial, threshold, connections)
    # check if sorted response is correct
    correct <- ifelse(which(paste(trial$target[1]) == out$collect) == subset(trial, paste(bar)==target)$size, 1, 0)
    target_position <- which(paste(trial$target[1]) == out$collect)
    # if not present, then it means it stopped too early, which is also incorrect
    correct <- ifelse(length(correct)==0, 0, correct)
    # return time and whether or not correct
    return(list(time = out$time[1], correct = correct, sortedSequence = out$collect, target_position = target_position))
  }else{
    out <- noisy_sort_backwards(trial, threshold, connections)
    # check if sorted response is correct
    correct <- ifelse(which(paste(trial$target[1]) == out$collect)==subset(trial, paste(bar)==target)$size, 1, 0)
    target_position <- which(paste(trial$target[1]) == out$collect)
    # if not present, then it means it stopped too early, which is also incorrect
    correct <- ifelse(length(correct)==0, 0, correct)
    # return time and whether or not correct
    return(list(time = out$time[1], correct = correct, sortedSequence = out$collect, target_position = target_position))
  }
  
}

# function that replaces bottom 5% with top 5%
die_and_replace <- function(data, fitness){
  # top 5% fittest programs
  top <- (fitness %in% sort(fitness)[1:20])
  # bottom 5% fittest programs
  bottom <- (fitness %in% sort(fitness)[(length(fitness)-19):length(fitness)])
  # replace
  data[bottom==1,] <- data[top==1,]
  # return new data
  return(data)
}

########################################
# bucket sort functions for the sorter
########################################

# bucket sort algorithm should scale 2n if no structure is learned
# trial: trial data frame, generated by generate trial
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
    #increment time
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
  print("sorted Sequence")
  print(collect)
  
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
  print("sorted Sequence")
  print(collect)
  
  
  # collect into data frame
  dout <- data.frame(collect, time)
  # return data frame
  return(dout)
}



######################################
# functions for the initialization and mutation of particles
####################################

#function to generate most conservative initial hypothesis
generateStartParticles = function(nOfBars = 7,
                                  nOfColors = 10){
  particle = list(threshold = nOfBars,
                         connections = "a",
                         startAtSmall = 1,
                         transitionmatrix = matrix(1/nOfColors, nOfColors, nOfColors),
                         directionprobs = c(0.5, 0.5), #first is the probability for starting at small and the second teh prob for starting at large
                         thresholdprobs = c(integer(nOfBars-1), 1)) 
  return(particle)
}

update_directionprobs = function(relative_target_position, # target position/sequence length
                                directionprobs, # vector with two entries that sum to 1, first is prob for start at small
                                learningrate = 1,
                                correct = 1){
  if (correct == 0){
    # if the response was incorrect, increase uncertainty by adding the learningrate to all values
    # and normalizing the vector
    directionprobs = directionprobs + learningrate
    directionprobs = directionprobs/sum(directionprobs)
    return(directionprobs)
  }else{
    # if response was correct, update the directionbrobs, so that the prob for the 
    # direction that would have the shortest sort
    # gets updated based on how much shorter the sort would have been (relative to the sort length)
    directionprobs[1] = abs(directionprobs[1] + relative_target_position*learningrate) # start at small
    directionprobs[2] = abs(directionprobs[2] - relative_target_position*learningrate) # start at large
    directionprobs = directionprobs/sum(directionprobs)
    return(directionprobs)
  }
}

update_thresholdprobs = function(thresholdprobs,
                                 target_position,
                                 sequence_length,
                                 startAtSmall = 1,
                                 learning_rate = 1,
                                 correct = 1){
  if (correct == 0){
    # if the response was incorrect, move towards the default threshold of 7, which grantees better performance
    thresholdprobs[length(thresholdprobs)] = thresholdprobs[length(thresholdprobs)] + learning_rate
    thresholdprobs = thresholdprobs/sum(thresholdprobs)
    return(thresholdprobs)
  }else{
    # increse the probability of sorting just to the target position (based on the direction of the sort)
    if (startAtSmall == 1){
      thresholdprobs[target_position] = thresholdprobs[target_position] + learning_rate
      thresholdprobs = thresholdprobs/sum(thresholdprobs)
    }else{
      position = sequence_length - target_position + 1
      thresholdprobs[position] = thresholdprobs[position] + learning_rate
      thresholdprobs = thresholdprobs/sum(thresholdprobs)
    }
    return(thresholdprobs)
  }
}

updateParticle = function(particle,
                          sortedSequence,
                          target_position,
                          sequence_length,
                          learningrate,
                          correct,
                          nOfBars = 7, 
                          nOfMaxConnections = 3, 
                          threshold = TRUE, 
                          connectedness = TRUE,
                          startAtSmall = TRUE){

  if(threshold == FALSE & connectedness == FALSE & startAtSmall == FALSE){
    print("particle could not be updated, because both threshold and connectedness 
          mutations where set to false, change those settings if you want to get updates")
    return(particle)
  }
  if (connectedness == TRUE){
  # update connectedness
    transitionmatrix = transitionmatrix_update(particle$transitionmatrix, 
                                               sortedSequence,
                                               learningrate,
                                               correct)
    new_sequence_hypothesis = get_hypothesis_transitions(transitionmatrix)
    particle$transitionmatrix = transitionmatrix
    particle$connections = new_sequence_hypothesis[1]
  }
  if (startAtSmall == TRUE){
    # update direction
    new_directionprobs = update_directionprobs(relative_target_position = ((sequence_length + 1)/2 - target_position)/7, # This represent the deviation from the mean
                                               particle$directionprobs, 
                                               learningrate,
                                               correct)
    particle$directionprobs = new_directionprobs
    if (which.max(new_directionprobs) == 1){
      particle$startAtSmall = 1
    }else{
      particle$startAtSmall = 0
    }
    
  }
  if (threshold == TRUE){
    # update threshold
    new_thresholdprobs = update_thresholdprobs(particle$thresholdprobs,
                                               target_position,
                                               sequence_length,
                                               startAtSmall = particle$startAtSmall,
                                               learningrate/2,
                                               correct)
    particle$thresholdprobs = new_thresholdprobs
    # make the mean position the threshold
    particle$threshold = round(sum(c(1:length(particle$thresholdprobs))*new_thresholdprobs))

  }
  return(particle)
}

###########################
# function to run the simulation
############################
runSimulation = function(start_particle, # the current particle
                         total_runs, # int
                         trials, # int
                         learningrate = 1,
                         structure_query = FALSE, 
                         structure_task = FALSE){
  #matrix to collect run time
  times = matrix(0, trials, total_runs)
  thresholds = matrix(0, trials, total_runs)
  connections = matrix(0, trials, total_runs)
  accuracy = matrix(0, trials, total_runs)
  directions = matrix(0, trials, total_runs)
  NoB = matrix(0, trials, total_runs) # Number of Bars
  finalHypotheses = data.frame(threshold = as.numeric(),
                               conncetions = as.character(),
                               startAtSmall = as.integer(),
                               time = as.numeric(),
                               run = as.numeric())

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
    
    print("current trial")
    print(trial)
    
    # create an out file
    # bind the output (time and correct) lists into data frame
    dout <- execute(trial = trial, 
                    threshold = start_particle$threshold, 
                    connections = start_particle$connections, 
                    startAtSmall = start_particle$startAtSmall)
    
    # bind them together
    hypotheses = start_particle
    hypotheses$time = dout$time 
    hypotheses$correct = dout$correct
    
    # go through 50 trials
    for (i in 1:trials){
      # generate new trial
      if(structure_query){
        trial = generatetrial(structure_query = TRUE)
      }else if(structure_task){
        trial = generatetrial(structure_task = TRUE)
      }else{
        trial = generatetrial()
      }
      
      NoB[i,n_runs] <- length(trial$bar)
      
      
      # collect results of all programs
      dout <-  execute(trial = trial, 
                       threshold = hypotheses$threshold, 
                       connections = hypotheses$connections, 
                       startAtSmall = hypotheses$startAtSmall)
      # print("dout")
      # print(dout)
      # update time so that each hypotheses time reflects the average time
      hypotheses$time <- hypotheses$time*((i-1)/i) + (1/i)*dout$time
      # update accuracy
      hypotheses$correct <- dout$correct
      
      # update the particle
      if (length(trial$bar) > 1){
        hypotheses = updateParticle(hypotheses,
                                    sortedSequence = dout$sortedSequence, 
                                    target_position = dout$target_position,
                                    sequence_length = length(trial$bar),
                                    learningrate,
                                    correct = dout$correct,
                                    nOfBars = 7, 
                                    nOfMaxConnections = 3, 
                                    threshold = TRUE, 
                                    connectedness = TRUE,
                                    startAtSmall = TRUE)
      }
     
      print("sorted sequence:")
      print(dout$sortedSequence)
      print("updated hypotheses")
      print(hypotheses)
      # track max time
      times[i, n_runs] <- max(dout$time)
      # track mean accuracy
      accuracy[i, n_runs] <- mean(dout$correct)
      # track the thresholds and connections
      thresholds[i, n_runs] <- hypotheses$threshold
      connections[i, n_runs] <- hypotheses$connections
      directions[i, n_runs] <- hypotheses$startAtSmall
      
      
    }

  }
  output <- list()
  output[[1]] <- times
  output[[2]] <- accuracy
  output[[3]] <- directions
  output[[4]] <- NoB
  output[[5]] <- thresholds
  output[[6]] <- connections
  return(output)
}

#############################################################################
# function to run model on participant data. 
# should take dataframe as input and should output:
# 1 the estimated rt
# the real rt difference (the real rt is added later)
# the real accuracy
# final hypotheses
###############################################################################
runHypoLearnerOnData = function(start_particle, # the current particle
                         data,
                         learningrate = 1,
                         learnThresholds = TRUE,
                         learnConnectedness = TRUE,
                         learnDirection = TRUE){
  trials = 35
  n_participants = length(unique(data$subject_id))
  # matrix to collect run time
  times = matrix(0, trials, n_participants)
  thresholds = matrix(0, trials, n_participants)
  connections = matrix(0, trials, n_participants)
  accuracy = matrix(0, trials, n_participants)
  directions = matrix(0, trials, n_participants)
  NoB = matrix(0, trials, n_participants) # Number of Bars
  # participant data
  realRTs = matrix(0, trials, n_participants)
  realAccuracy = matrix(0, trials, n_participants)
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
  
  # for 20 simulations in total
  for (participant in 1:n_participants){
    print("started run:")
    print(participant)
    # get current trial
    current_trial = subset(data, subject_id == unique(data$subject_id)[participant])[1,]
    trial = datacolum_to_trial(current_trial)
    
    print("current trial")
    print(trial)
    
    # create an out file
    # bind the output (time and correct) lists into data frame
    dout <- execute(trial = trial, 
                    threshold = start_particle$threshold, 
                    connections = start_particle$connections, 
                    startAtSmall = start_particle$startAtSmall)
    
    # bind them togethers
    hypotheses = start_particle
    hypotheses$time = dout$time 
    hypotheses$correct = dout$correct
    
    # go through 50 trials
    for (i in 1:trials){
      #get current trial
      current_trial = subset(data, subject_id == unique(data$subject_id)[participant])[i,]
      trial = datacolum_to_trial(current_trial)
      
      NoB[i, participant] <- length(trial$bar)
      
      
      # collect results of all programs
      dout <-  execute(trial = trial, 
                       threshold = hypotheses$threshold, 
                       connections = hypotheses$connections, 
                       startAtSmall = hypotheses$startAtSmall)

      # update time so that each hypotheses time reflects the average time
      hypotheses$time <- hypotheses$time*((i-1)/i) + (1/i)*dout$time
      # update accuracy
      hypotheses$correct <- dout$correct
      
      # update the particle
      if (length(trial$bar) > 1){
        hypotheses = updateParticle(hypotheses,
                                    sortedSequence = dout$sortedSequence, 
                                    target_position = dout$target_position,
                                    sequence_length = length(trial$bar),
                                    learningrate,
                                    correct = dout$correct,
                                    nOfBars = 7, 
                                    nOfMaxConnections = 3, 
                                    threshold = learnThresholds, 
                                    connectedness = learnConnectedness,
                                    startAtSmall = learnDirection)
      }
      print("sorted sequence:")
      print(dout$sortedSequence)
      print("updated hypotheses")
      print(hypotheses)
      # track max time
      times[i, participant] <- max(dout$time)
      # track mean accuracy
      accuracy[i, participant] <- mean(dout$correct)
      # track the thresholds and connections
      thresholds[i, participant] <- hypotheses$threshold
      connections[i, participant] <- hypotheses$connections
      directions[i, participant] <- hypotheses$startAtSmall
      # track real RT
      realRTs[i, participant] <- current_trial$RTSort_Memory
      # track real accuracy
      realAccuracy[i, participant] <- current_trial$correct
      
      
    }

  }
  output <- list()
  output[[1]] <- times
  output[[2]] <- realRTs
  output[[3]] <- accuracy
  output[[4]] <- realAccuracy
  output[[5]] <- directions
  output[[6]] <- NoB
  output[[7]] <- thresholds
  output[[8]] <- connections
  output[[9]] <- correctConnection
  return(output)
}

#############################################################################
# run model on all data
##############################################################################

runHypoLearnerOnAllData = function(start_particle,
                                   AllData,
                                   learningrate,
                                   learnThresholds = TRUE,
                                   learnConnectedness = TRUE,
                                   learnDirection = TRUE){
  ######################################
  # query structure 
  #######################################
  dQuery <- subset(AllData, Structure == "Query" & Condition == "Sort" & Stimulus_type == "bars")
  
  QueryDataOutput <- runHypoLearnerOnData(start_particle, # the current particle
                                          data = dQuery,
                                          learningrate,
                                          learnThresholds,
                                          learnConnectedness,
                                          learnDirection)
  
  Querytimes <- QueryDataOutput[[1]]
  QueryrealRTs <- QueryDataOutput[[2]] 
  Queryaccuracy <- QueryDataOutput[[3]] 
  QueryrealAccuracy <- QueryDataOutput[[4]] 
  Querydirections <- QueryDataOutput[[5]] 
  QueryNoB <- QueryDataOutput[[6]]
  Querythresholds <- QueryDataOutput[[7]] 
  Queryconnections <- QueryDataOutput[[8]]
  QuerytrueConnection <- QueryDataOutput[[9]]
  
  ModelToData = data.frame(NoB = as.numeric(),
                           time = as.numeric(),
                           realRT = as.numeric(),
                           accurcay = as.integer(),
                           realAccuracy = as.integer(),
                           threshold = as.numeric(),
                           connection = as.character(),
                           startAtSmall = as.numeric(),
                           participant = as.numeric(),
                           Structure = as.character(), 
                           trial = as.numeric(),
                           trueConnection = as.numeric())
  
  for (i in 1:length(unique(dQuery$subject_id))){
    ModelToData = rbind(ModelToData, data.frame(NoB = QueryNoB[,i],
                                                time = Querytimes[,i],
                                                realRT = QueryrealRTs[,i],
                                                accuracy = Queryaccuracy[,i],
                                                realAccuracy = QueryrealAccuracy[,i],
                                                threshold = Querythresholds[,i],
                                                connection = Queryconnections[,i],
                                                startAtSmall = Querydirections[,i],
                                                participant = unique(dQuery$subject_id)[i],
                                                Structure = "Query", 
                                                trial = 1:35,
                                                trueConnection = QuerytrueConnection[,i]))
  }
  
  ######################################
  # sequence structure 
  #######################################
  dSequence <- subset(AllData, Structure == "Sequence" & Condition == "Sort" & Stimulus_type == "bars")
  
  SequenceDataOutput <- runHypoLearnerOnData(start_particle, # the current particle
                                             data = dSequence,
                                             learningrate,
                                             learnThresholds,
                                             learnConnectedness,
                                             learnDirection)
  
  Sequencetimes<- SequenceDataOutput[[1]]
  SequencerealRTs <- SequenceDataOutput[[2]] 
  Sequenceaccuracy <- SequenceDataOutput[[3]] 
  SequencerealAccuracy <- SequenceDataOutput[[4]] 
  Sequencedirections <- SequenceDataOutput[[5]] 
  SequenceNoB <- SequenceDataOutput[[6]]
  Sequencethresholds <- SequenceDataOutput[[7]] 
  Sequenceconnections <- SequenceDataOutput[[8]] 
  SequencetrueConnection <- SequenceDataOutput[[9]]
  
  for (i in 1:length(unique(dSequence$subject_id))){
    ModelToData = rbind(ModelToData, data.frame(NoB = SequenceNoB[,i],
                                                time = Sequencetimes[,i],
                                                realRT = SequencerealRTs[,i],
                                                accuracy = Sequenceaccuracy[,i],
                                                realAccuracy = SequencerealAccuracy[,i],
                                                threshold = Sequencethresholds[,i],
                                                connection = Sequenceconnections[,i],
                                                startAtSmall = Sequencedirections[,i],
                                                participant = unique(dSequence$subject_id)[i],
                                                Structure = "Sequence", 
                                                trial = 1:35,
                                                trueConnection = SequencetrueConnection[,i]))
  }
  
  ######################################
  # no structure 
  #######################################
  dNone <- subset(AllData, Structure == "None" & Condition == "Sort" & Stimulus_type == "bars")
  
  NoneDataOutput <- runHypoLearnerOnData(start_particle, # the current particle
                                         data = dNone,
                                         learningrate,
                                         learnThresholds,
                                         learnConnectedness,
                                         learnDirection)
  
  Nonetimes<- NoneDataOutput[[1]]
  NonerealRTs <- NoneDataOutput[[2]] 
  Noneaccuracy <- NoneDataOutput[[3]] 
  NonerealAccuracy <- NoneDataOutput[[4]] 
  Nonedirections <- NoneDataOutput[[5]] 
  NoneNoB <- NoneDataOutput[[6]]
  Nonethresholds <- NoneDataOutput[[7]] 
  Noneconnections <- NoneDataOutput[[8]] 
  NonetrueConnection <- NoneDataOutput[[9]]
  
  for (i in 1:length(unique(dNone$subject_id))){
    ModelToData = rbind(ModelToData, data.frame(NoB = NoneNoB[,i],
                                                time = Nonetimes[,i],
                                                realRT = NonerealRTs[,i],
                                                accuracy = Noneaccuracy[,i],
                                                realAccuracy = NonerealAccuracy[,i],
                                                threshold = Nonethresholds[,i],
                                                connection = Noneconnections[,i],
                                                startAtSmall = Nonedirections[,i],
                                                participant = unique(dNone$subject_id)[i],
                                                Structure = "None", 
                                                trial = 1:35,
                                                trueConnection = NonetrueConnection[,i]))
  }
  return(ModelToData)
}



################################################################################
# functions fo a potential MSE based hyperparameter selection (was not used)
################################################################################

calculateSquaredErrorMatrix <- function(output){
  SEMatrix = matrix(0, dim(output[[1]])[1], dim(output[[1]])[2])
  for (i in 1: dim(output[[1]])[2])
    SEMatrix[,i] = (scale(output[[2]][,i]) - scale(output[[1]][,i]))^2
  return(SEMatrix)
}

getMSE <- function(start_particle, # the current particle
                   data,
                   learningrates #vector of learn ingrates that sorter is run on
                   ){
  allMSE <- c(1:length(learningrates))
  for (i in 1:length(learningrates)){
    output = runHypoLearnerOnData(start_particle, # the current particle
                                   data,
                                   learningrate = learningrates[i])
    allMSE[i] = calculateSquaredErrorMatrix(output)
  }
  MSEandLR = data.frame(Learningrate = learningrates,
                        MSE = allMSE)
  return(MSEandLR)
}

#standard error
se <- function(x){sd(x)/sqrt(length(x))}

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

########################################################
# ll functions for hyper parameter selection (Not used)
########################################################
ModelLL <- function(AllData, learningrate){
  ModelToData <- runHypoLearnerOnAllData(start_particle = generateStartParticles(nOfBars = 7, nOfColors = 10),
                                         AllData,
                                         learningrate)
  #####################################
  # exlude all data where the participant was incorrect or th real RT is larger than 10/-10
  ############################################
  ModelToData <-subset(ModelToData, realAccuracy == 1 & realRT > -10 & realRT < 10)
  output <- lm(realRT ~ time, data = ModelToData)
  intercept <- as.numeric(output$coefficients[1])
  beta <- as.numeric(output$coefficients[2])
  ModelToData$timeasRT = intercept + ModelToData$time*beta
  #calculate the loglikelihood
  LL <- 0
  for (i in 1:length(ModelToData$timeasRT)){
    LL <- LL + dnorm(x = ModelToData$realRT[i], mean = ModelToData$timeasRT[i], sd = 10, log = TRUE)
  }
  
  return(LL)
}

ModelLLrawRT <- function(AllData, learningrate){
  ModelToData <- runHypoLearnerOnAllData(start_particle = generateStartParticles(nOfBars = 7, nOfColors = 10),
                                         AllData,
                                         learningrate)
  # this part is necesarry, so I can add the correct columns to the model dataframe
  relevantD <- subset(AllData, Condition == "Sort" & Stimulus_type == "bars")
  relevantD1 <- subset(relevantD, Structure == "Query")
  relevantD2 <- subset(relevantD, Structure == "Sequence")
  relevantD3 <- subset(relevantD, Structure == "None")
  relevantD <- rbind(relevantD1, relevantD2, relevantD3)
  #####################################
  # exlude all data where the participant was incorrect or th real RT is larger than 10/-10
  ############################################
  ModelToData$rawRt <- relevantD$rt
  ModelToData <- (subset(ModelToData, realAccuracy == 1 & rawRt <= RTCutoff ))
  output <- lm(rawRt ~ time, data = ModelToData)
  intercept <- as.numeric(output$coefficients[1])
  beta <- as.numeric(output$coefficients[2])
  ModelToData$timeasRT = intercept + ModelToData$time*beta
  # calculate the loglikelihood
  LL <- 0
  for (i in 1:length(ModelToData$timeasRT)){
    LL <- LL + dnorm(x = ModelToData$rawRt[i], mean = ModelToData$timeasRT[i], sd = 10, log = TRUE)
  }
  
  return(LL)
}

