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

neg_after_FA <- subset(tmp, (w_type == 143 | w_type == 144) & w_resp == 51)
pos_after_FA <- subset(tmp, (w_type == 243 | w_type == 244) & w_resp == 52)

neg_after_FH <- subset(tmp, w_type == 141 & w_resp == 51)
pos_after_FH <- subset(tmp, w_type == 241 & w_resp == 52)




##################   calculate priming   ##################################
##################  effect with median   ##################################

# calculate median before exclusion of outliers

median_neg_after_FA <- median(neg_after_FA$w_rt)
median_pos_after_FA <- median(pos_after_FA$w_rt)
median_neg_after_FH <- median(neg_after_FH$w_rt)
median_pos_after_FH <- median(pos_after_FH$w_rt)

median_priming_after_error <- median_pos_after_FA - median_neg_after_FA
median_priming_after_correct <- median_neg_after_FH - median_pos_after_FH
median_priming_overall <- (median_pos_after_FA + median_neg_after_FH) - (median_neg_after_FA + median_pos_after_FH) 




####################     exclude   #####################################
####################    outliers   #####################################

# assign NA to rt_values deviating more than 3 standard deviations from median
neg_after_FA$w_rt[which((abs(neg_after_FA$w_rt - median(neg_after_FA$w_rt))/mad(neg_after_FA$w_rt))>3)] <- NA
pos_after_FA$w_rt[which((abs(pos_after_FA$w_rt - median(pos_after_FA$w_rt))/mad(pos_after_FA$w_rt))>3)] <- NA
neg_after_FH$w_rt[which((abs(neg_after_FH$w_rt - median(neg_after_FH$w_rt))/mad(neg_after_FH$w_rt))>3)] <- NA
pos_after_FH$w_rt[which((abs(pos_after_FH$w_rt - median(pos_after_FH$w_rt))/mad(pos_after_FH$w_rt))>3)] <- NA

# count number of outliers
count_outlier <- length(which(is.na(neg_after_FA$w_rt))) + length(which(is.na(pos_after_FA$w_rt))) + length(which(is.na(neg_after_FH$w_rt))) + length(which(is.na(pos_after_FH$w_rt)))


# remove the rows with rt_values = NA
neg_after_FA <- neg_after_FA[!(is.na(neg_after_FA$w_rt)),]
pos_after_FA <- pos_after_FA[!(is.na(pos_after_FA$w_rt)),]
neg_after_FH <- neg_after_FH[!(is.na(neg_after_FH$w_rt)),]
pos_after_FH <- pos_after_FH[!(is.na(pos_after_FH$w_rt)),]




####################    calculat priming    ##################################
####################    effect with mean    ##################################

mean_neg_after_FA <- mean(neg_after_FA$w_rt)
mean_pos_after_FA <- mean(pos_after_FA$w_rt)
mean_neg_after_FH <- mean(neg_after_FH$w_rt)
mean_pos_after_FH <- mean(pos_after_FH$w_rt)

sd_neg_after_FA <- sd(neg_after_FA$w_rt)
sd_pos_after_FA <- sd(pos_after_FA$w_rt)
sd_neg_after_FH <- sd(neg_after_FH$w_rt)
sd_pos_after_FH <- sd(pos_after_FH$w_rt)

mean_priming_after_error <- mean_pos_after_FA - mean_neg_after_FA
mean_priming_after_correct <- mean_neg_after_FH - mean_pos_after_FH
mean_priming_overall <- (mean_pos_after_FA + mean_neg_after_FH) - (mean_neg_after_FA + mean_pos_after_FH) 




####################    count number    ##################################
####################     of events      ##################################

# counts number of events, after outlier are removed; if I want to work with median, shift this part before section "exclude outlier"

count_neg_after_FA <- length(neg_after_FA$w_rt)
count_pos_after_FA <- length(pos_after_FA$w_rt)
count_neg_after_FH <- length(neg_after_FH$w_rt)
count_pos_after_FH <- length(pos_after_FH$w_rt)

count_incorr_neg_after_FA <- length(which((tmp$w_type == 143 | tmp$w_type == 144) & tmp$w_resp == 53))
count_incorr_pos_after_FA <- length(which((tmp$w_type == 243 | tmp$w_type == 244) & tmp$w_resp == 54))
count_incorr_neg_after_FH <- length(which(tmp$w_type == 141 & tmp$w_resp == 53))
count_incorr_pos_after_FH <- length(which(tmp$w_type == 241 & tmp$w_resp == 54))

count_miss_neg_after_FA <- length(which((tmp$w_type == 143 | tmp$w_type == 144) & tmp$w_resp == 55))
count_miss_pos_after_FA <- length(which((tmp$w_type == 243 | tmp$w_type == 244) & tmp$w_resp == 56))
count_miss_neg_after_FH <- length(which(tmp$w_type == 141 & tmp$w_resp == 55))
count_miss_pos_after_FH <- length(which(tmp$w_type == 241 & tmp$w_resp == 56))

# misses are contained in overall number of events so that errors, correct and misses sum up to 100 %
count_all_neg_after_FA <- count_neg_after_FA + count_incorr_neg_after_FA + count_miss_neg_after_FA
count_all_pos_after_FA <- count_pos_after_FA + count_incorr_pos_after_FA + count_miss_pos_after_FA
count_all_neg_after_FH <- count_neg_after_FH + count_incorr_neg_after_FH + count_miss_neg_after_FH
count_all_pos_after_FH <- count_pos_after_FH + count_incorr_pos_after_FH + count_miss_pos_after_FH



####################    calculate    ##################################
####################    accuracy     ##################################

percent_correct_neg_after_FA <- 
accuracy_pos_after_FA <- 
accuracy_neg_after_FH <- 
accuracy_pos_after_FH <- 





}



#setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/Analysis")  # setting a different folder as working directory to prevent saving stuff into the folder containing the logfiles

# possibly use command "attach"; then I do not have to indicate columns with $