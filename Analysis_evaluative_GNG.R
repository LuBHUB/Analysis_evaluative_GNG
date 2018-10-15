##### Analysis Behavioral Experiment on Affective Priming Effect 2018 ##### 
##### by Luisa Balzus
##### 15.10.2018

# clear environment
rm(list=ls())

####################   load   #######################################
####################   data   #######################################

# load logfiles, create dataframe 

setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/raw_data")                       #path to folder containing the log files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 

logfiles <- list.files("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/raw_data")            #lists files in folder

for (i in logfiles)                                                                                      # loop reading txt-file by txt-file as table, ommit first 58 lines and added lines after trials; use first line as header
{ tmp <- (read.table(i, skip = 58, fill = TRUE, header = TRUE, nrows = 516))
#tmp$name <- substr(tmp$name,8,9)                                                                  # change "names" column: have 02 insted of ModERN_02_Flanker
#tmp <- tmp[!duplicated(tmp$word),]                                                                # delete repeated words
#if (i=="Luisa_004.txt")
#{a1 <- tmp}                                                                                       # at the first iteration, create a1
#else {a1 <- rbind(a1,tmp)}                                                                        # rbind glues currently loaded file (tmp) to the end of dataframe a1



####################     define    #####################################
####################   conditions  #####################################

# FA = false alarm; FH = Fast Hit

neg_after_FA <- subset(tmp, w_type == 143 | w_type == 144 & w_resp == 51)
pos_after_FA <- subset(tmp, w_type == 243 | w_type == 244 & w_resp == 52)

neg_after_FH <- subset(tmp, w_type == 141 & w_resp == 51)
pos_after_FH <- subset(tmp, w_type == 241 & w_resp == 52)



####################     define    #####################################
####################    outliers   #####################################

# assign NA to rt_values deviating more than 3 standard deviations from median
neg_after_FA$w_rt[which((abs(neg_after_FA$w_rt - median(neg_after_FA$w_rt))/mad(neg_after_FA$w_rt))>3)] <- NA
pos_after_FA$w_rt[which((abs(pos_after_FA$w_rt - median(pos_after_FA$w_rt))/mad(pos_after_FA$w_rt))>3)] <- NA
neg_after_FH$w_rt[which((abs(neg_after_FH$w_rt - median(neg_after_FH$w_rt))/mad(neg_after_FH$w_rt))>3)] <- NA
pos_after_FH$w_rt[which((abs(pos_after_FH$w_rt - median(pos_after_FH$w_rt))/mad(pos_after_FH$w_rt))>3)] <- NA

# remove the rows with rt_values = NA
neg_after_FA <- neg_after_FA[!(is.na(neg_after_FA$w_rt)),]
pos_after_FA <- pos_after_FA[!(is.na(pos_after_FA$w_rt)),]
neg_after_FH <- neg_after_FH[!(is.na(neg_after_FH$w_rt)),]
pos_after_FH <- pos_after_FH[!(is.na(pos_after_FH$w_rt)),]




####################     calculate     ##################################
####################     mean, sd      ##################################

mean_neg_after_FA <- mean(neg_after_FA$w_rt)
mean_pos_after_FA <- mean(pos_after_FA$w_rt)
mean_neg_after_FH <- mean(neg_after_FH$w_rt)
mean_pos_after_FH <- mean(pos_after_FH$w_rt)

sd_neg_after_FA <- sd(neg_after_FA$w_rt)
sd_pos_after_FA <- sd(pos_after_FA$w_rt)
sd_neg_after_FH <- sd(neg_after_FH$w_rt)
sd_pos_after_FH <- sd(pos_after_FH$w_rt)



}



#setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/Analysis")  # setting a different folder as working directory to prevent saving stuff into the folder containing the logfiles
