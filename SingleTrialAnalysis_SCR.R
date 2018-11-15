##### Single Trial Analysis SCR in Behavioral Experiment on Affective Priming Effect 2018 ##### 
##### by Luisa Balzus
##### 09.11.2018


library(dplyr)

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
scrfiles <- list.files("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/9_SCR_export_preprocessed") 
setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/9_SCR_export_preprocessed")               # path to folder containing the scr files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 

for (subject in scrfiles){                                                                              # loop reading csv-file
raw_scr <- read.csv(subject)

if (subject=="ModERN_behav_01_SCR_Export_era_z.csv")
{scr <- raw_scr}                                                                                        # at the first iteration, create a1
else {scr <- rbind(scr,raw_scr)}                                                                        # rbind glues currently loaded file (tmp) to the end of dataframe a1
}


####################     merge log and     #######################################
####################       scr data        #######################################

# return all rows of log, merged with matching columns of scr; rows in log with no match in scr will have NAs in the new columns; if there are multiple matches between log and scr, all combinations of the matches are returned
master <- left_join(log,scr,by = c("number"= "Participant.ID","count" = "Trial.ID","t_type"="Event.ID.Tar","t_resp"="Event.ID.TarResp","w_type"="Event.ID.Word","w_resp"="Event.ID.WordResp"))

# check whether there were rows in log with no match in scr
if (any(is.na(master$DDA.AreaSum.TarResp))){
  stop("There are rows in log with no match in SCR!")
}else{
    print("Matching of logfile and SCR-file is ok. Each row in log has a machtching row in SCR")
  }


# check whether there were there are multiple matches between log and scr
if (any(duplicated(master[,c('number','count')]))){
  stop("There are rows in log with multiple matches in SCR!")
}else{
  print("Matching of logfile and SCR-file is ok. No row in log has several matching row in SCR")
}








####################     loop over        #######################################
####################      subjects        #######################################


# create data frame for saving values for each subject in the for-loop
df4save = NULL

for (name in unique(master$name)){
  
 

  ####################       define      #####################################
  ####################  GNG conditions  #####################################
  

  GNG_FA <- subset(master, t_resp == 43 | t_resp == 44)
  GNG_FH <- subset(master, t_resp == 41)
  GNG_SH <- subset(master, t_resp == 42)

  
  ####################   calculate  #####################################
  ####################    GNG RT    ##################################### 
  
  GNG_mean_FA <- mean(GNG_FA$t_rt)
  GNG_mean_FH <- mean(GNG_FH$t_rt)
  GNG_mean_SH <- mean(GNG_SH$t_rt)
  
  

df4save <- rbind(df4save, data.frame(name,GNG_mean_FA,GNG_mean_FH,GNG_mean_SH))

}
















