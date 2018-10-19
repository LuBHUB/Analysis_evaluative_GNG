##### Analysis Behavioral Experiment on Affective Priming Effect 2018 ##### 
##### by Luisa Balzus
##### 15.10.2018


# install.packages('openxlsx') for saving data in Excel sheet
library(openxlsx)
# install.packages("pastecs")
library(pastecs)
# install.packages("ggplot2") for plots
library(ggplot2)
# install.packages("ez") for ANOVA
library(ez)
# install.packages("nlme") for LMM
library(nlme)
# install.packages("tidyr")
library(tidyr)
# install.packages("lmerTest")
library(lmerTest)
# install.packages("psycho")
library(psycho)
# install.packages("multcomp")
library(multcomp)
# install.packages("emmeans")
library(emmeans)


# clear environment
rm(list=ls())

# create data frame for saving values for each subject in the for-loop
df4save = NULL

# force R to not use exponential notation
options(scipen = 999)

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



####################     exclude      #####################################
####################  GNG outliers   #####################################

# assign NA to GNG rt_values < 150 or > 500 ms (according to Pourtois)
tmp$t_rt[which((tmp$t_resp != 45 & tmp$t_resp != 46 & tmp$t_resp != 47) & tmp$t_rt < 150 | tmp$t_rt > 500)] <- NA

# count number of outliers
count_outlier_GNG <- length(which(is.na(tmp$t_rt))) 

# remove the rows with rt_values = NA
tmp_GNG <- tmp[!(is.na(tmp$t_rt)),]

# these outliers are excluded only for GNG trial, not for word categorization in corresponding trial (as in Pourtois)


####################       define      #####################################
####################  GNG conditions  #####################################

# FA = False Alarm; FH = Fast Hit; CI = Correctly Inhibited; SH = Slow Hit

GNG_FA <- subset(tmp_GNG, t_resp == 43 | t_resp == 44)
GNG_FH <- subset(tmp_GNG, t_resp == 41)
GNG_CI <- subset(tmp_GNG, t_resp == 45 | t_resp == 46)
GNG_SH <- subset(tmp_GNG, t_resp == 42)
GNG_miss <- subset(tmp_GNG, t_resp == 47)
GNG_wrong_key <- subset(tmp_GNG,  t_resp == 48 | t_resp == 49)


####################   calculate  #####################################
####################    GNG RT    ##################################### 

# not calculated for CI, because there is no rt
GNG_mean_FA <- mean(GNG_FA$t_rt)
GNG_mean_FH <- mean(GNG_FH$t_rt)
GNG_mean_SH <- mean(GNG_SH$t_rt)

GNG_sd_FA <- sd(GNG_FA$t_rt)
GNG_sd_FH <- sd(GNG_FH$t_rt)
GNG_sd_SH <- sd(GNG_SH$t_rt)

GNG_cov_FA <- sd(GNG_FA$t_rt)/mean(GNG_FA$t_rt) * 100
GNG_cov_FH <- sd(GNG_FH$t_rt)/mean(GNG_FH$t_rt) * 100
GNG_cov_SH <- sd(GNG_SH$t_rt)/mean(GNG_SH$t_rt) * 100


####################   calculate GNG     #####################################
####################  % response types   ##################################### 

GNG_count_FA <- length(GNG_FA$t_rt)
GNG_count_FH <- length(GNG_FH$t_rt)
GNG_count_CI <- length(GNG_CI$t_rt)
GNG_count_SH <- length(GNG_SH$t_rt)
GNG_count_miss <- length(GNG_miss$t_rt)

GNG_all_resp_without_outliers_wrong_keys <- GNG_count_FA + GNG_count_FH + GNG_count_CI + GNG_count_SH + GNG_count_miss 

GNG_percent_FA <- (GNG_count_FA/GNG_all_resp_without_outliers_wrong_keys)*100
GNG_percent_FH <- (GNG_count_FH/GNG_all_resp_without_outliers_wrong_keys)*100
GNG_percent_CI <- (GNG_count_CI/GNG_all_resp_without_outliers_wrong_keys)*100
GNG_percent_SH <- (GNG_count_SH/GNG_all_resp_without_outliers_wrong_keys)*100
GNG_percent_miss <- (GNG_count_miss/GNG_all_resp_without_outliers_wrong_keys)*100



####################       define      #####################################
####################  word conditions  #####################################

# FA = False Alarm; FH = Fast Hit; CI = Correctly Inhibited; SH = Slow Hit

neg_after_FA <- subset(tmp, (w_type == 143 | w_type == 144) & w_resp == 51)
pos_after_FA <- subset(tmp, (w_type == 243 | w_type == 244) & w_resp == 52)

neg_after_FH <- subset(tmp, w_type == 141 & w_resp == 51)
pos_after_FH <- subset(tmp, w_type == 241 & w_resp == 52)

neg_after_CI <- subset(tmp, (w_type == 145 | w_type == 146) & w_resp == 51)
pos_after_CI <- subset(tmp, (w_type == 245 | w_type == 246) & w_resp == 52)

neg_after_SH <- subset(tmp, w_type == 142 & w_resp == 51)
pos_after_SH <- subset(tmp, w_type == 242 & w_resp == 52)



##################   calculate priming   ##################################
##################  effect with median   ##################################

# calculate median before exclusion of word outliers

median_neg_after_FA <- median(neg_after_FA$w_rt)
median_pos_after_FA <- median(pos_after_FA$w_rt)
median_neg_after_FH <- median(neg_after_FH$w_rt)
median_pos_after_FH <- median(pos_after_FH$w_rt)
median_neg_after_CI <- median(neg_after_CI$w_rt)
median_pos_after_CI <- median(pos_after_CI$w_rt)
median_neg_after_SH <- median(neg_after_SH$w_rt)
median_pos_after_SH <- median(pos_after_SH$w_rt)

median_priming_after_FA <- median_pos_after_FA - median_neg_after_FA
median_priming_after_FH <- median_neg_after_FH - median_pos_after_FH
median_priming_overall <- (median_pos_after_FA + median_neg_after_FH) - (median_neg_after_FA + median_pos_after_FH) 




####################     exclude      #####################################
####################  word outliers   #####################################

# assign NA to rt_values deviating more than 3 median absolute deviations (MAD; better use 2.5?)
neg_after_FA$w_rt[which((abs(neg_after_FA$w_rt - median(neg_after_FA$w_rt))/mad(neg_after_FA$w_rt))>3)] <- NA
pos_after_FA$w_rt[which((abs(pos_after_FA$w_rt - median(pos_after_FA$w_rt))/mad(pos_after_FA$w_rt))>3)] <- NA
neg_after_FH$w_rt[which((abs(neg_after_FH$w_rt - median(neg_after_FH$w_rt))/mad(neg_after_FH$w_rt))>3)] <- NA
pos_after_FH$w_rt[which((abs(pos_after_FH$w_rt - median(pos_after_FH$w_rt))/mad(pos_after_FH$w_rt))>3)] <- NA
neg_after_CI$w_rt[which((abs(neg_after_CI$w_rt - median(neg_after_CI$w_rt))/mad(neg_after_CI$w_rt))>3)] <- NA
pos_after_CI$w_rt[which((abs(pos_after_CI$w_rt - median(pos_after_CI$w_rt))/mad(pos_after_CI$w_rt))>3)] <- NA
neg_after_SH$w_rt[which((abs(neg_after_SH$w_rt - median(neg_after_SH$w_rt))/mad(neg_after_SH$w_rt))>3)] <- NA
pos_after_SH$w_rt[which((abs(pos_after_SH$w_rt - median(pos_after_SH$w_rt))/mad(pos_after_SH$w_rt))>3)] <- NA

# Alternative: assign NA to rt_values deviating more than 2.5 standard deviations around the mean (this corresponds to Pourtois, but MAD is better practice)
# neg_after_FA$w_rt[which((abs(neg_after_FA$w_rt - mean(neg_after_FA$w_rt))/sd(neg_after_FA$w_rt)) > 2.5)] <- NA
# pos_after_FA$w_rt[which((abs(pos_after_FA$w_rt - mean(pos_after_FA$w_rt))/sd(pos_after_FA$w_rt)) > 2.5)] <- NA
# neg_after_FH$w_rt[which((abs(neg_after_FH$w_rt - mean(neg_after_FH$w_rt))/sd(neg_after_FH$w_rt)) > 2.5)] <- NA
# pos_after_FH$w_rt[which((abs(pos_after_FH$w_rt - mean(pos_after_FH$w_rt))/sd(pos_after_FH$w_rt)) > 2.5)] <- NA
# neg_after_CI$w_rt[which((abs(neg_after_CI$w_rt - mean(neg_after_CI$w_rt))/sd(neg_after_CI$w_rt)) > 2.5)] <- NA
# pos_after_CI$w_rt[which((abs(pos_after_CI$w_rt - mean(pos_after_CI$w_rt))/sd(pos_after_CI$w_rt)) > 2.5)] <- NA
# neg_after_SH$w_rt[which((abs(neg_after_SH$w_rt - mean(neg_after_SH$w_rt))/sd(neg_after_SH$w_rt)) > 2.5)] <- NA
# pos_after_SH$w_rt[which((abs(pos_after_SH$w_rt - mean(pos_after_SH$w_rt))/sd(pos_after_SH$w_rt)) > 2.5)] <- NA

# count number of outliers
count_outlier_words_FA_FH <- length(which(is.na(neg_after_FA$w_rt))) + length(which(is.na(pos_after_FA$w_rt))) + length(which(is.na(neg_after_FH$w_rt))) + length(which(is.na(pos_after_FH$w_rt)))

# remove the rows with rt_values = NA
neg_after_FA <- neg_after_FA[!(is.na(neg_after_FA$w_rt)),]
pos_after_FA <- pos_after_FA[!(is.na(pos_after_FA$w_rt)),]
neg_after_FH <- neg_after_FH[!(is.na(neg_after_FH$w_rt)),]
pos_after_FH <- pos_after_FH[!(is.na(pos_after_FH$w_rt)),]
neg_after_CI <- neg_after_CI[!(is.na(neg_after_CI$w_rt)),]
pos_after_CI <- pos_after_CI[!(is.na(pos_after_CI$w_rt)),]
neg_after_SH <- neg_after_SH[!(is.na(neg_after_SH$w_rt)),]
pos_after_SH <- pos_after_SH[!(is.na(pos_after_SH$w_rt)),]




####################    calculat priming    ##################################
####################    effect with mean    ##################################

mean_neg_after_FA <- mean(neg_after_FA$w_rt)
mean_pos_after_FA <- mean(pos_after_FA$w_rt)
mean_neg_after_FH <- mean(neg_after_FH$w_rt)
mean_pos_after_FH <- mean(pos_after_FH$w_rt)
mean_neg_after_CI <- mean(neg_after_CI$w_rt)
mean_pos_after_CI <- mean(pos_after_CI$w_rt)
mean_neg_after_SH <- mean(neg_after_SH$w_rt)
mean_pos_after_SH <- mean(pos_after_SH$w_rt)

mean_priming_after_FA <- mean_pos_after_FA - mean_neg_after_FA
mean_priming_after_FH <- mean_neg_after_FH - mean_pos_after_FH
mean_priming_overall <- (mean_pos_after_FA + mean_neg_after_FH) - (mean_neg_after_FA + mean_pos_after_FH) 




####################    count number    ##################################
####################     of events      ##################################

# counts number of events, after outliers are removed; if I want to work with median, shift this part before section "exclude outlier"

count_neg_after_FA <- length(neg_after_FA$w_rt)
count_pos_after_FA <- length(pos_after_FA$w_rt)
count_neg_after_FH <- length(neg_after_FH$w_rt)
count_pos_after_FH <- length(pos_after_FH$w_rt)
count_neg_after_CI <- length(neg_after_CI$w_rt)
count_pos_after_CI <- length(pos_after_CI$w_rt)
count_neg_after_SH <- length(neg_after_SH$w_rt)
count_pos_after_SH <- length(pos_after_SH$w_rt)

count_incorr_neg_after_FA <- length(which((tmp$w_type == 143 | tmp$w_type == 144) & tmp$w_resp == 53))
count_incorr_pos_after_FA <- length(which((tmp$w_type == 243 | tmp$w_type == 244) & tmp$w_resp == 54))
count_incorr_neg_after_FH <- length(which(tmp$w_type == 141 & tmp$w_resp == 53))
count_incorr_pos_after_FH <- length(which(tmp$w_type == 241 & tmp$w_resp == 54))
count_incorr_neg_after_CI <- length(which((tmp$w_type == 145 | tmp$w_type == 146) & tmp$w_resp == 53))
count_incorr_pos_after_CI <- length(which((tmp$w_type == 245 | tmp$w_type == 246) & tmp$w_resp == 54))
count_incorr_neg_after_SH <- length(which(tmp$w_type == 142 & tmp$w_resp == 53))
count_incorr_pos_after_SH <- length(which(tmp$w_type == 242 & tmp$w_resp == 54))

count_miss_neg_after_FA <- length(which((tmp$w_type == 143 | tmp$w_type == 144) & tmp$w_resp == 55))
count_miss_pos_after_FA <- length(which((tmp$w_type == 243 | tmp$w_type == 244) & tmp$w_resp == 56))
count_miss_neg_after_FH <- length(which(tmp$w_type == 141 & tmp$w_resp == 55))
count_miss_pos_after_FH <- length(which(tmp$w_type == 241 & tmp$w_resp == 56))
count_miss_neg_after_CI <- length(which((tmp$w_type == 145 | tmp$w_type == 146) & tmp$w_resp == 55))
count_miss_pos_after_CI <- length(which((tmp$w_type == 245 | tmp$w_type == 246) & tmp$w_resp == 56))
count_miss_neg_after_SH <- length(which(tmp$w_type == 142 & tmp$w_resp == 55))
count_miss_pos_after_SH <- length(which(tmp$w_type == 242 & tmp$w_resp == 56))

# misses are contained in overall number of events so that errors, correct and misses sum up to 100 %
count_all_neg_after_FA <- count_neg_after_FA + count_incorr_neg_after_FA + count_miss_neg_after_FA
count_all_pos_after_FA <- count_pos_after_FA + count_incorr_pos_after_FA + count_miss_pos_after_FA
count_all_neg_after_FH <- count_neg_after_FH + count_incorr_neg_after_FH + count_miss_neg_after_FH
count_all_pos_after_FH <- count_pos_after_FH + count_incorr_pos_after_FH + count_miss_pos_after_FH
count_all_neg_after_CI <- count_neg_after_CI + count_incorr_neg_after_CI + count_miss_neg_after_CI
count_all_pos_after_CI <- count_pos_after_CI + count_incorr_pos_after_CI + count_miss_pos_after_CI
count_all_neg_after_SH <- count_neg_after_SH + count_incorr_neg_after_SH + count_miss_neg_after_SH
count_all_pos_after_SH <- count_pos_after_SH + count_incorr_pos_after_SH + count_miss_pos_after_SH

pos_after_SH <- count_pos_after_SH + count_incorr_pos_after_SH + count_miss_pos_after_SH


####################    calculate    ##################################
####################    accuracy     ##################################

percent_correct_neg_after_FA <- count_neg_after_FA / count_all_neg_after_FA * 100 
percent_correct_pos_after_FA <- count_pos_after_FA / count_all_pos_after_FA * 100  
percent_correct_neg_after_FH <- count_neg_after_FH / count_all_neg_after_FH * 100  
percent_correct_pos_after_FH <- count_pos_after_FH / count_all_pos_after_FH * 100  
percent_correct_neg_after_CI <- count_neg_after_CI / count_all_neg_after_CI * 100 
percent_correct_pos_after_CI <- count_pos_after_CI / count_all_pos_after_CI * 100  
percent_correct_neg_after_SH <- count_neg_after_SH / count_all_neg_after_SH * 100  
percent_correct_pos_after_SH <- count_pos_after_SH / count_all_pos_after_SH * 100  

#percent_incorrect_neg_after_FA <- count_incorr_neg_after_FA / count_all_neg_after_FA * 100 
#percent_incorrect_pos_after_FA <- count_incorr_pos_after_FA / count_all_pos_after_FA * 100  
#percent_incorrect_neg_after_FH <- count_incorr_neg_after_FH / count_all_neg_after_FH * 100  
#percent_incorrect_pos_after_FH <- count_incorr_pos_after_FH / count_all_pos_after_FH * 100  

#percent_miss_neg_after_FA <- count_miss_neg_after_FA / count_all_neg_after_FA * 100 
#percent_miss_pos_after_FA <- count_miss_pos_after_FA / count_all_pos_after_FA * 100  
#percent_miss_neg_after_FH <- count_miss_neg_after_FH / count_all_neg_after_FH * 100  
#percent_miss_pos_after_FH <- count_miss_pos_after_FH / count_all_pos_after_FH * 100  




####################    save values   ##################################
####################   in dataframe   ##################################


df4save <- rbind(df4save, data.frame(subject,mean_neg_after_FA,mean_pos_after_FA,mean_neg_after_FH,mean_pos_after_FH,mean_neg_after_CI,mean_pos_after_CI,mean_neg_after_SH,mean_pos_after_SH,mean_priming_overall,mean_priming_after_FA,mean_priming_after_FH,median_neg_after_FA,median_pos_after_FA,median_neg_after_FH,median_pos_after_FH,median_neg_after_CI,median_pos_after_CI,median_neg_after_SH,median_pos_after_SH,median_priming_overall,median_priming_after_FA,median_priming_after_FH,count_neg_after_FA,count_pos_after_FA,count_neg_after_FH,count_pos_after_FH,count_outlier_words_FA_FH,percent_correct_neg_after_FA,percent_correct_pos_after_FA,percent_correct_neg_after_FH,percent_correct_pos_after_FH,percent_correct_neg_after_CI,percent_correct_pos_after_CI,percent_correct_neg_after_SH,percent_correct_pos_after_SH,GNG_mean_FA,GNG_mean_FH,GNG_mean_SH,GNG_sd_FA,GNG_sd_FH,GNG_sd_SH,GNG_cov_FA,GNG_cov_FH,GNG_cov_SH,GNG_percent_FA,GNG_percent_FH,GNG_percent_CI,GNG_percent_SH,GNG_percent_miss))
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

# plot rt 
  df4plotRT <- data.frame(
    response = c("false alarm","false alarm","fast hit","fast hit"),
    valence = c("neg","pos","neg","pos"),
    conditions = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH"),
    mean = as.numeric(c(descriptive_statistics[2,2],descriptive_statistics[2,3],descriptive_statistics[2,4],descriptive_statistics[2,5])),
    se = as.numeric(c(descriptive_statistics[3,2],descriptive_statistics[3,3],descriptive_statistics[3,4],descriptive_statistics[3,5])))

  ggplot(df4plotRT,x = response,y = mean, aes(response, mean, fill = valence))+
  geom_bar(stat="identity", position=position_dodge()) +                                                       # add bars, based on stats values; dodge to avoid stacked bars
  geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.95), width=0.1) +    # add error bars
  ggtitle("Mean Response Time") +                                                                              # add title
  xlab("response") + ylab("Reaction Time in ms") +                                                             # label axes
  theme(plot.title = element_text(hjust = 0.5)) +                                                              # center title
  scale_fill_manual(values=c("mediumblue", "limegreen"))                                                       # change bar colors
  
  
# plot accuracy
  df4plotACC <- data.frame(
    response = c("false alarm","false alarm","fast hit","fast hit"),
    valence = c("neg","pos","neg","pos"),
    conditions = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH"),
    mean = as.numeric(c(descriptive_statistics[2,29],descriptive_statistics[2,30],descriptive_statistics[2,31],descriptive_statistics[2,32])),
    se = as.numeric(c(descriptive_statistics[3,29],descriptive_statistics[3,30],descriptive_statistics[3,31],descriptive_statistics[3,32])))

  ggplot(df4plotACC,x = response,y = mean, aes(response, mean, fill = valence))+
  geom_bar(stat="identity", position=position_dodge()) +                                                       # add bars, based on stats values; dodge to avoid stacked bars
  geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.95), width=0.1) +    # add error bars
  ggtitle("Mean Accuracy") +                                                                                   # add title
  xlab("response") + ylab("Accuracy in %") +                                                                   # label axes
  theme(plot.title = element_text(hjust = 0.5)) +                                                              # center title
  scale_fill_manual(values=c("mediumblue", "limegreen"))                                                       # change bar colors
  
  
  #####################         ANOVA         ####################################
  #####################    words (mean rt)    #################################### 
  
  options(contrasts=c("contr.sum","contr.poly")) # to adhere to the sum-to-zero convention for effect weights, always do this before running ANOVAs in R. This matters sometimes (not always). If I don’t do it, the sum of squares calculations may not match what I get e.g. in SPSS
  
  
  # ANOVA requires several rows for each subject, each one row per factor: here I reorder data by creating new data frame
  # df4anova <- data.frame(
  #                       subject = df4save$subject,
  #                       response_type = rep(c("false alarm","false alarm","fast hit","fast hit","correctly inhibited","correctly inhibited", "slow hit", "slow hit"),length(df4save$subject)), # repeat response types as often as number of subjects
  #                       word_valence = rep(c("neg","pos","neg","pos","neg","pos","neg","pos"),length(df4save$subject)),                                                                        # repeat word valence as often as number of subjects
  #                       rt = as.vector(t(df4save[,2:9])),
  #                       accuracy = as.vector(t(df4save[,29:36]))
  #                       )   
  #
  # df4anova$subject <- sort(df4anova$subject)
  #
  # df4anova$subject <- factor(df4anova$subject)
  # df4anova$response_type <- factor(df4anova$response_type)
  # df4anova$word_valence <- factor(df4anova$word_valence)
  
  
  
  
  # ANOVA requires several rows for each subject, each one row per factor: here I reorder data by reshaping existing data frame
  df4anova <-   reshape(data = df4save[,c(1:9,29:36)], 
                direction = "long",
                varying = list(c("mean_neg_after_FA","mean_pos_after_FA","mean_neg_after_FH","mean_pos_after_FH","mean_neg_after_CI","mean_pos_after_CI" ,"mean_neg_after_SH","mean_pos_after_SH"),c("percent_correct_neg_after_FA", "percent_correct_pos_after_FA","percent_correct_neg_after_FH","percent_correct_pos_after_FH","percent_correct_neg_after_CI","percent_correct_pos_after_CI","percent_correct_neg_after_SH","percent_correct_pos_after_SH")),
                v.names = c("rt","accuracy"),
                idvar = "subject",
                timevar = "condition",
                times = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH","neg_after_CI","pos_after_CI","neg_after_SH","pos_after_SH")
                )
 
  row.names(df4anova) <- NULL
  df4anova <- df4anova[sort.list(df4anova$subject),]                        # sort df by subject
  
  df4anova$response_type <- substr(df4anova$condition, 11, 12)              # create columns needed as factors           
  df4anova$word_valence <- substr(df4anova$condition, 1, 3)

  df4anova$subject <- factor(df4anova$subject)                              # create factors response_type and word_valence and reorder factor levels! (this is important, because first level is used as reference factor in LMM)
  df4anova$response_type <- factor(df4anova$response_type, levels=c("FA","FH","CI","SH"))
  df4anova$word_valence <- factor(df4anova$word_valence)                    # automatic factor order is alphabetical (1 = neg, 2 = pos)
  
  

  # calculate ANOVA with ezANOVA for rt -> gives same result as SPSS :)
  anova_rt <- ezANOVA(data = df4anova, 
                    dv = .(rt), 
                    wid = .(subject), 
                    within = .(response_type, word_valence), 
                    detailed = TRUE)
  
  # calculate ANOVA with ezANOVA for accuracy -> gives same result as SPSS :)
  anova_accuracy <- ezANOVA(data = df4anova, 
                      dv = .(accuracy), 
                      wid = .(subject), 
                      within = .(response_type, word_valence), 
                      detailed = TRUE)
 


  # calculate ANOVA with aov for rt -> gives same result as SPSS :) (in this option, test for sphericity (Mauchly) AND post hoc tests are NOT possible) 
  anova_rt <- aov(rt ~ (response_type * word_valence) + 
                  Error(subject/(response_type * word_valence)), 
                  data = df4anova)
  summary(anova_rt)                                                                             # to get ANOVA output
  
  

  
  # calculate ANOVA with aov for accuracy -> gives same result as SPSS :) (in this option, test for sphericity (Mauchly) AND post hoc tests are NOT possible)
  anova_accuracy <- aov(accuracy ~ (response_type * word_valence) + 
                    Error(subject/(response_type * word_valence)), 
                    data = df4anova)
  summary(anova_accuracy)                                                                       # to get ANOVA output 
  
  
  
  
  # post hoc comparison: when calculating within-subject ANOVAs, no direct post hoc test is possible; I can only run paitwise t-tests with corresponding adjustment (e.g. Holm, which is better than Bonferroni)
  anova_rt_posthoc_response_type <- pairwise.t.test(df4anova$rt,df4anova$response_type,p.adjust.method="holm")
  anova_rt_posthoc_word_valence <- pairwise.t.test(df4anova$rt,df4anova$word_valence,p.adjust.method="holm") 
  anova_rt_posthoc_priming <- pairwise.t.test(df4anova$rt,df4anova$condition,p.adjust.method="holm") 
  
  anova_accuracy_posthoc_response_type <- pairwise.t.test(df4anova$accuracy,df4anova$response_type,p.adjust.method="holm")
  anova_accuracy_posthoc_word_valence <- pairwise.t.test(df4anova$accuracy,df4anova$word_valence,p.adjust.method="holm") 
  anova_accuracy_posthoc_priming <- pairwise.t.test(df4anova$accuracy,df4anova$condition,p.adjust.method="holm") 
  

  
  #####################    linear mixed models    ####################################
  #####################      words (mean RT)      #################################### 
  # advantage of LLMs over ANOVAs: ANOVA requires many strong assumptions, such as homoscedasticity, which is hard to verify; furthermore, LMMs are more flexible
  # an ANOVA is pretty much a condensed linear model where the predictors are factors. Therefore, I can run an ANOVA on a linear mixed model (which includes the “error” term, or random effect). The results are, for the important bits (the sum of squares, mean square and p value), very close to those of the traditional approach.

  
  # calculate LMM using lmer for rt and plot
  model_lmer_rt <- lmer(rt ~ response_type * word_valence + (1|subject), data=df4anova)        # possibly add REML = FALSE?
  anova(model_lmer_rt)                                                                         # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
  results_model_lmer_rt <- analyze(model_lmer_rt)                                              # print results
  print(results_model_lmer_rt)
  results_model_lmer_rt <- get_contrasts(model_lmer_rt, "response_type * word_valence")        # Provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
  print(results_model_lmer_rt$contrasts)                                                       # print contrasts
  print(results_model_lmer_rt$means)                                                           # investigate means

  ggplot(results_model_lmer_rt$means, aes(x=response_type, y=Mean, color=word_valence, group=word_valence)) +    # plot 
    geom_line(position = position_dodge(.3)) +
    geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher), 
                    position = position_dodge(.3)) +
    ylab("Response Time in ms") +
    xlab("Response Type") +
    theme_bw()
  
  
  
  # calculate LMM using lmer for accuracy and plot 
  model_lmer_accuracy <- lmer(accuracy ~ response_type * word_valence + (1|subject), data=df4anova)  # possibly add REML = FALSE?
  anova(model_lmer_accuracy)                                                                         # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
  results_model_lmer_accuracy <- analyze(model_lmer_accuracy)                                        # print results
  print(results_model_lmer_accuracy)
  results_model_lmer_accuracy <- get_contrasts(model_lmer_accuracy, "response_type * word_valence", adjust="holm")  # Provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
  print(results_model_lmer_accuracy$contrasts)                                                       # print contrasts
  print(results_model_lmer_accuracy$means)                                                           # investigate means
  
  ggplot(results_model_lmer_accuracy$means, aes(x=response_type, y=Mean, color=word_valence, group=word_valence)) +    # plot 
    geom_line(position = position_dodge(.3)) +
    geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher), 
                    position = position_dodge(.3)) +
    ylab("Accuracy in %") +
    xlab("Response Type") +
    theme_bw() 

  
  
  # calculate LMM using lme for rt
  model_lme_rt <- lme(rt ~ response_type * word_valence, random =~1|subject, data=df4anova)
  anova(model_lme_rt)
  summary(model_lme_rt)
  
  

  # calculate LMM using lme for accuracy
  model_lme_accuracy <- lme(accuracy ~ response_type * word_valence, random =~1|subject, data=df4anova)
  anova(model_lme_accuracy)
  summary(model_lme_accuracy)
  
  
  
  # post hoc tests for resolving main effects 
  summary(glht(model_lmer_rt, linfct=mcp(response_type = "Tukey")), test = adjusted("holm"))  
  summary(glht(model_lmer_rt, linfct=mcp(word_valence = "Tukey")), test = adjusted("holm"))   # same result can be obtained by lmer_rt_posthoc_response_type <- emmeans(model_lmer_rt, ~ word_valence);;pairs(lmer_rt_posthoc_response_type)
 
  summary(glht(model_lmer_accuracy, linfct=mcp(response_type = "Tukey")), test = adjusted("holm"))  
  summary(glht(model_lmer_accuracy, linfct=mcp(word_valence = "Tukey")), test = adjusted("holm"))
  
  # post hoc tests for resolving interaction effects 
  lmer_rt_posthoc_priming <- emmeans(model_lmer_rt, ~ word_valence|response_type, adjust="tukey")    # compare word valence within each level of response_type
  pairs(lmer_rt_posthoc_priming) 
   
  lmer_accuracy_posthoc_priming <- emmeans(model_lmer_accuracy, ~ word_valence|response_type, adjust="tukey")
  pairs(lmer_accuracy_posthoc_priming) 
  
   
  
  
  
  #####################         ANOVA         ####################################
  #####################     GNG (for rt)      #################################### 
  
  
  # ANOVA requires several rows for each subject, each one row per factor: here I reorder data by reshaping existing data frame
  df4anova_GNG_rt <-   reshape(data = df4save[,c(1,37:39)], 
                        direction = "long",
                        varying = list(c("GNG_mean_FA","GNG_mean_FH","GNG_mean_SH")),
                        v.names = c("rt"),
                        idvar = "subject",
                        timevar = "condition",
                        times = c("FA","FH","SH"))
  
  row.names(df4anova_GNG_rt) <- NULL
  df4anova_GNG_rt <- df4anova_GNG_rt[sort.list(df4anova_GNG_rt$subject),]                        # sort df by subject
  
  df4anova_GNG_rt$subject <- factor(df4anova_GNG_rt$subject)                                     # create factors
  df4anova_GNG_rt$condition <- factor(df4anova_GNG_rt$condition)

  
  # calculate ANOVA with ezANOVA for rt -> gives same result as SPSS :)
  anova_GNG_rt <- ezANOVA(data = df4anova_GNG_rt, 
                      dv = .(rt), 
                      wid = .(subject), 
                      within = .(condition), 
                      detailed = TRUE)
  
  
  # post hoc comparison: when calculating within-subject ANOVAs, no direct post hoc test is possible; I can only run paitwise t-tests with corresponding adjustment (e.g. Holm, which is better than Bonferroni)
  anova_GNG_rt_posthoc_response_type <- pairwise.t.test(df4anova_GNG_rt$rt,df4anova_GNG_rt$condition,p.adjust.method="holm")
  
  
  
  #####################            ANOVA                    ####################################
  #####################     GNG (for percent responses)     #################################### 
  
  
  # ANOVA requires several rows for each subject, each one row per factor: here I reorder data by reshaping existing data frame
  df4anova_GNG_percent <-   reshape(data = df4save[,c(1,46:49)], 
                               direction = "long",
                               varying = list(c("GNG_percent_FA","GNG_percent_FH","GNG_percent_CI","GNG_percent_SH")),
                               v.names = c("percent"),
                               idvar = "subject",
                               timevar = "condition",
                               times = c("FA","FH","CI","SH")
  )
  
  row.names(df4anova_GNG_percent) <- NULL
  df4anova_GNG_percent <- df4anova_GNG_percent[sort.list(df4anova_GNG_percent$subject),]                        # sort df by subject
  
  df4anova_GNG_percent$subject <- factor(df4anova_GNG_percent$subject)                                          # create factors
  df4anova_GNG_percent$condition <- factor(df4anova_GNG_percent$condition)
  
  
  

  # calculate ANOVA with ezANOVA for percent responses
  anova_GNG_percent <- ezANOVA(data = df4anova_GNG_percent, 
                            dv = .(percent), 
                            wid = .(subject), 
                            within = .(condition), 
                            detailed = TRUE)
  
  
  
  # post hoc comparison: when calculating within-subject ANOVAs, no direct post hoc test is possible; I can only run paitwise t-tests with corresponding adjustment (e.g. Holm, which is better than Bonferroni)
  anova_GNG_percent_posthoc_response_type <- pairwise.t.test(df4anova_GNG_percent$percent,df4anova_GNG_percent$condition,p.adjust.method="holm")
  

  
  #####################    linear mixed models    ####################################
  #####################          GNG (RT)         #################################### 
  
  # calculate LMM using lmer for rt and plot
  model_lmer_GNG_rt <- lmer(rt ~ condition + (1|subject), data=df4anova_GNG_rt)                    # possibly add REML = FALSE? 
  anova(model_lmer_GNG_rt)                                                                         # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
  results_model_lmer_GNG_rt <- analyze(model_lmer_GNG_rt)                                          # print results
  print(results_model_lmer_GNG_rt)
  results_model_lmer_GNG_rt <- get_contrasts(model_lmer_GNG_rt, "condition")                       # Provide the model and the factors to contrast;; by default, get_contrasts uses the Tukey method for p value adjustment; add ,adjust="none" to turn off automatic p value correction after Tucky
  print(results_model_lmer_GNG_rt$contrasts)                                                       # print contrasts
  print(results_model_lmer_GNG_rt$means)                                                           # investigate means
  
  ggplot(results_model_lmer_GNG_rt$means, aes(x=condition, y=Mean, color=condition)) +             # plot 
    geom_line(position = position_dodge(.3)) +
    geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher), 
                    position = position_dodge(.3)) +
    ylab("Response Time in ms") +
    xlab("Response Type") +
    theme_bw()
  
  
  
  # calculate LMM using lmer for percent and plot 
  model_lmer_GNG_percent <- lmer(percent ~ condition + (1|subject), data=df4anova_GNG_percent)          # possibly add REML = FALSE?
  anova(model_lmer_GNG_percent)                                                                         # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
  results_model_lmer_GNG_percent <- analyze(model_lmer_GNG_percent)                                     # print results
  print(results_model_lmer_GNG_percent)
  results_model_lmer_GNG_percent <- get_contrasts(model_lmer_GNG_percent, "condition")                  # Provide the model and the factors to contrast;; by default, get_contrasts uses the Tukey method for p value adjustment; add ,adjust="none" to turn off automatic p value correction after Tucky
  print(results_model_lmer_GNG_percent$contrasts)                                                       # print contrasts
  print(results_model_lmer_GNG_percent$means)                                                           # investigate means
  
  ggplot(results_model_lmer_GNG_percent$means, aes(x=condition, y=Mean, color=condition)) +             # plot 
    geom_line(position = position_dodge(.3)) +
    geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher), 
                    position = position_dodge(.3)) +
    ylab("Frequency in %") +
    xlab("Response Type") +
    theme_bw() 
  
  
  
  # calculate LMM using lme for rt
  model_lme_GNG_rt <- lme(rt ~ condition, random =~1|subject, data=df4anova_GNG_rt)
  anova(model_lme_GNG_rt)
  summary(model_lme_GNG_rt)
  
  
  
  # calculate LMM using lme for accuracy
  model_lme_GNG_percent <- lme(percent ~ condition, random =~1|subject, data=df4anova_GNG_percent)
  anova(model_lme_GNG_percent)
  summary(model_lme_GNG_percent)
  
  
  # post hoc tests for resolving main effects 
  summary(glht(model_lmer_GNG_rt, linfct=mcp(condition = "Tukey")), test = adjusted("holm"))  
  summary(glht(model_lmer_GNG_percent, linfct=mcp(condition = "Tukey")), test = adjusted("holm"))  