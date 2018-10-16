##### Analysis Behavioral Experiment on Affective Priming Effect 2018 ##### 
##### by Luisa Balzus
##### 15.10.2018


# install.packages('openxlsx')
library(openxlsx)
# install.packages("pastecs")
library(pastecs)
# install.packages("ggplot2")
library(ggplot2)


# clear environment
rm(list=ls())

# create data frame for saving values for each subject in the for-loop
df4save = NULL


####################   load   #######################################
####################   data   #######################################

# load logfiles, create dataframe 

setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/raw_data")                           #path to folder containing the log files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 

logfiles <- list.files("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/raw_data")          #lists files in folder

for (subject in logfiles)                                                                          # loop reading txt-file by txt-file as table, ommit first 58 lines and added lines after trials; use first line as header
{ tmp <- (read.table(subject, skip = 58, fill = TRUE, header = TRUE, nrows = 516))
#tmp$name <- substr(tmp$name,8,9)                                                                  # change "names" column: have 02 insted of ModERN_02_Flanker
#tmp <- tmp[!duplicated(tmp$word),]                                                                # delete repeated words
#if (subject=="Luisa_004.txt")
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

# assign NA to rt_values deviating more than 3 median absolute deviations

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

# counts number of events, after outliers are removed; if I want to work with median, shift this part before section "exclude outlier"

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

percent_correct_neg_after_FA <- count_neg_after_FA / count_all_neg_after_FA * 100 
percent_correct_pos_after_FA <- count_pos_after_FA / count_all_pos_after_FA * 100  
percent_correct_neg_after_FH <- count_neg_after_FH / count_all_neg_after_FH * 100  
percent_correct_pos_after_FH <- count_pos_after_FH / count_all_pos_after_FH * 100  

percent_incorrect_neg_after_FA <- count_incorr_neg_after_FA / count_all_neg_after_FA * 100 
percent_incorrect_pos_after_FA <- count_incorr_pos_after_FA / count_all_pos_after_FA * 100  
percent_incorrect_neg_after_FH <- count_incorr_neg_after_FH / count_all_neg_after_FH * 100  
percent_incorrect_pos_after_FH <- count_incorr_pos_after_FH / count_all_pos_after_FH * 100  

percent_miss_neg_after_FA <- count_miss_neg_after_FA / count_all_neg_after_FA * 100 
percent_miss_pos_after_FA <- count_miss_pos_after_FA / count_all_pos_after_FA * 100  
percent_miss_neg_after_FH <- count_miss_neg_after_FH / count_all_neg_after_FH * 100  
percent_miss_pos_after_FH <- count_miss_pos_after_FH / count_all_pos_after_FH * 100  




####################    save values   ##################################
####################   in dataframe   ##################################


df4save <- rbind(df4save, data.frame(subject,mean_neg_after_FA,mean_pos_after_FA,mean_pos_after_FH,mean_neg_after_FH,sd_neg_after_FA,sd_pos_after_FA,sd_pos_after_FH,sd_neg_after_FH,mean_priming_overall,mean_priming_after_error,mean_priming_after_correct,median_neg_after_FA,median_pos_after_FA,median_pos_after_FH,median_neg_after_FH,median_priming_overall,median_priming_after_error,median_priming_after_correct,count_neg_after_FA,count_pos_after_FA,count_pos_after_FH,count_miss_neg_after_FH,count_outlier,count_incorr_neg_after_FA,count_incorr_pos_after_FA,count_incorr_pos_after_FH,count_incorr_neg_after_FH,percent_correct_neg_after_FA,percent_correct_pos_after_FA,percent_correct_pos_after_FH,percent_correct_neg_after_FH,percent_incorrect_neg_after_FA,percent_incorrect_pos_after_FA,percent_incorrect_pos_after_FH,percent_incorrect_neg_after_FH,percent_miss_neg_after_FA,percent_miss_pos_after_FA,percent_miss_pos_after_FH,percent_miss_neg_after_FH))

}



####################    save dataframe   ##################################
####################    as excel file    ##################################


setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/Analysis")                                      # setting a different folder as working directory to prevent saving stuff into the folder containing the logfiles

date_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

filename <- paste("SummaryStatisticsEvaluativeGNG_For",length(logfiles),"_subjects_",date_time, ".xlsx", sep = "")

write.xlsx(df4save, file = filename)
#save(df4save, file = filename)

#current_file <- read.xlsx(filename)                                                                         # read in that file


###########################################################################
###########################################################################
###########################################################################
##############################Statistics###################################
###########################################################################
###########################################################################
###########################################################################


#####################    descriptive   ####################################
#####################    statistics    ####################################

# get descriptve statistics
descriptive_statistics <- stat.desc(df4save,basic=F)


#####################      bar     ####################################
#####################     plots    ####################################

# for RT 
df4plotRT <- data.frame(
  response = c("false alarm","false alarm","fast hit","fast hit"),
  valence = c("neg","pos","neg","pos"),
  conditions = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH"),
  mean = as.numeric(c(descriptive_statistics[2,2],descriptive_statistics[2,3],descriptive_statistics[2,5],descriptive_statistics[2,4])),
  se = c(70,100,90,80)
  #for now I had only one subject, so I made se up; later use: se = as.numeric(c(descriptive_statistics[3,2],descriptive_statistics[3,3],descriptive_statistics[3,5],descriptive_statistics[3,4]))
  )

# plot_rt <- 
  ggplot(df4plotRT,x = response,y = mean, aes(response, mean, fill = valence))+
  geom_bar(stat="identity", position=position_dodge()) +                                                       # add bars, based on stats values; dodge to avoid stacked bars
  geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.95), width=0.1) +    # add error bars
  ggtitle("Mean Response Time") +                                                                              # add title
  xlab("response") + ylab("Reaction Time in ms") +                                                             # label axes
  theme(plot.title = element_text(hjust = 0.5)) +                                                              # center title
  scale_fill_manual(values=c("mediumblue", "limegreen"))                                                       # change bar colors
  
  
# for ACCURACY
  df4plotACC <- data.frame(
    response = c("false alarm","false alarm","fast hit","fast hit"),
    valence = c("neg","pos","neg","pos"),
    conditions = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH"),
    mean = as.numeric(c(descriptive_statistics[2,29],descriptive_statistics[2,30],descriptive_statistics[2,32],descriptive_statistics[2,31])),
    se = c(5,4,3,2)
    #for now I had only one subject, so I made se up; later use: se = as.numeric(c(descriptive_statistics[3,29],descriptive_statistics[3,30],descriptive_statistics[3,32],descriptive_statistics[3,31]))
  )
  
  
  # plot_accuracy <- 
  ggplot(df4plotACC,x = response,y = mean, aes(response, mean, fill = valence))+
  geom_bar(stat="identity", position=position_dodge()) +                                                       # add bars, based on stats values; dodge to avoid stacked bars
  geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.95), width=0.1) +    # add error bars
  ggtitle("Mean Accuracy") +                                                                                   # add title
  xlab("response") + ylab("Accuracy in %") +                                                                   # label axes
  theme(plot.title = element_text(hjust = 0.5)) +                                                              # center title
  scale_fill_manual(values=c("mediumblue", "limegreen"))                                                       # change bar colors
  
  
