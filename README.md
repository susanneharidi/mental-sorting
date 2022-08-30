# The Scaling of Mental Sorting
Susanne Haridi, Charley Wu, Ishita Dasgupta, and Eric Schulz

This is the Repository corresponding to the "Running head: MENTAL SORTING 1
The scaling of mental computation in a sorting task" manuscript

## This Repo contains three distinct parts

### 1. The Code for the behavioral experiment

This code has been altered slightly to protect the keys for the database (so data is no longer being saved). But otherwise all parts of the experiment are included and can be explored and reused freely. 
The experimental part of the data from this Experiment is contained in  batches the Data/Experiment folder

### 2. The analysis of the behavioral results 
(all scripts for this part are in the analyses folder)

The prepare_data.R script takes the raw data and transforms it into a more usable shape for all further analysis. The outputted data is in the data_to_work_with.csv file and only contains data of the 73 included participants. All other scripts use this file for further analysis.

Each row in the file is one trial and the variables of this dataframe are as follows:
"Trial_index" : The position of the trial in the experiment. Each participant has 210 trials
"TrialID" : For better comparison we matched trials with different structures (none, query, sequence) and tasks (memory and sort) for the length of the sequence, the heights of the rectangles and the order for the three structure conditions. Furthermore, within each structure condition (so across tasks) the same position was queried. The colors were not matched to prevent memory effects. Accordingly each participant should have 6 matched trials with the same ID.
"Block" : the current Block number. Each participant went through 6 Blocks, one for each condition.  
"Conditions" : unique identifiers for each condition
"Stimulus_type" : this was used to differentiate between the Recall and the query RT, in the current dataframe, this variable is no longer meaningful
"Number_of_Bars" : the sequence lengths of the current sequence
"rt" : the Encoding RT in seconds                 
"Query_color" : the color of the query
"color_sequence" : all colors of the sequence. The Order corresponds to the order of the sorted sequence
"position_sequence" : the x position corresponding to where on the screen each rectangle was presented. the first vale corresponds to the position of the smallest rectangle and so on.
"heights" : the heights of the rectangles.            
"key_press": This can be ignored
"Correct_Position" : the position of the queried rectangle in the sorted sequence
"response_position" : the position which the participant responded with
"subject_id" : a unique, randomly generated ID for each subject          
"time_elapsed" : the elapsed time since the beginning of the experiment
"date" :  no mystery there, its the date of the experiment                
"Condition" : the task of the current block (memory vs sort)
"Structure" : the structure of the current block (none, query, sequence)
"CondOrder" : a variable that codes for both the condition (task and structure) and the block
"RTSort_Memory" : the difference between the Encoding RT of the sort task and the matched memory task trial with the same structure in s        
"RTSort_MemoryNone" : the difference between the Encoding RT of the sort task and the matched memory task trial with no structure in s 
"MatchCorrectMemory" : whether the matched memory trial with the same structure was correct
"MatchRTMemory" : the encoding RT of the matched memory trial with the same structure
"MatchCorrectMemoryNone" : whether the matched memory trial with no structure was correct
"MatchRTMemoryNone" : the encoding RT of the matched memory trial with no structure
"queryRT" : the Recall RT in s

### 3. The Sorting Models


