##### Analysis Behavioral Experiment on Affective Priming Effect 2018 ##### 
##### by Luisa Balzus
##### 15.10.2018

# clear environment
rm(list=ls())

####################   load   #######################################
####################   data   #######################################

# load logfiles, create dataframe 

setwd("P:/LuisaBalzus/1_PhD_Project/1_Planning/Behavioral_Study_2/raw_data")                       #path to folder containing the log files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 

fl <- list.files("P:/LuisaBalzus/1_PhD_Project/1_Planning/Behavioral_Study_2/raw_data")            #lists files in folder

for (i in fl)                                                                                      # loop reading txt-file by txt-file as table, ommit first 58 lines and added lines after trials; use first line as header
{ tmp <- (read.table(i, skip = 58, fill = TRUE, header = TRUE, nrows = 516))
#tmp$name <- substr(tmp$name,8,9)                                                                  # change "names" column: have 02 insted of ModERN_02_Flanker
#tmp <- tmp[!duplicated(tmp$word),]                                                                # delete repeated words
#if (i=="Luisa_004.txt")
#{a1 <- tmp}                                                                                       # at the first iteration, create a1
#else {a1 <- rbind(a1,tmp)}                                                                        # rbind glues currently loaded file (tmp) to the end of dataframe a1
}


#setwd("P:/LuisaBalzus/1_PhD_Project/1_Planning/Behavioral_Study_2/analyses")  # setting a different folder as working directory to prevent saving stuff into the folder containing the logfiles
