##### Single Trial Analysis SCR in Behavioral Experiment on Affective Priming Effect 2018 ##### 
##### by Luisa Balzus
##### 09.11.2018


# clear environment
rm(list=ls())

####################     load     #######################################
####################   log data   #######################################

# Load logfiles
logfiles <- list.files("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_behav")       # lists files in folder
setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_behav")                        # path to folder containing the log files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 

for (subject in logfiles){                                                                              # loop reading txt-file by txt-file as table, ommit first 58 lines and added lines after trials; use first line as header
raw_log <- (read.table(subject, skip = 58, fill = TRUE, header = TRUE, nrows = 516))

if (subject=="ModERN_behav_01.txt")
{log <- raw_log}                                                                                        # at the first iteration, create a1
else {log <- rbind(log,raw_log)}                                                                        # rbind glues currently loaded file (tmp) to the end of dataframe a1
}


####################     load     #######################################
####################   scr data   #######################################

# Load SCR files (Ledalab output exorted as .csv)
scrfiles <- list.files("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/8_SCR_export/SCR_export_preprocessed") 
setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/8_SCR_export/SCR_export_preprocessed")    # path to folder containing the scr files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 

for (subject in scrfiles){                                                                              # loop reading csv-file
raw_scr <- read.csv2(subject)

if (subject=="ModERN_behav_01_SCR_Export_era_z.csv")
{scr <- raw_scr}                                                                                        # at the first iteration, create a1
else {scr <- rbind(scr,raw_scr)}                                                                        # rbind glues currently loaded file (tmp) to the end of dataframe a1
}

#### Add loop for each ID here



####To DO
# cbind relevant log rows to logs for all matching IDs (IDS, TrialID,TargetID,targetrespID, wordID and WordrespID have to be identical) -> find out what is missing for subj 5 
# make loop over log ID to do proprocessing for each subj in mastertable


