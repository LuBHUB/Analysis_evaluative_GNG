##### Single Trial Analysis SCR in Behavioral Experiment on Affective Priming Effect 2018 ##### 
##### by Luisa Balzus
##### 09.11.2018


library(dplyr)
library(openxlsx)
library(pastecs)
library(ggplot2)
library(lmerTest)
library(multcomp)
library(emmeans)
library(nlme)
library(tidyr)
library(psycho)
library(car)
library(foreign)
library(psych)


# clear environment
rm(list=ls())

# force R to not use exponential notation
options(scipen = 999)

####################     load     #######################################
####################   log data   #######################################

# Load logfiles
logfiles <- list.files("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_behav")       # lists files in folder
setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_behav")                        # path to folder containing the log files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 

for (subject in logfiles){                                                                              # loop reading txt-file by txt-file as table, ommit first 58 lines and added lines after trials; use first line as header
raw_log <- (read.table(subject, skip = 58, fill = TRUE, header = TRUE, nrows = 516))

if (subject=="ModERN_behav_01.txt")
{log <- raw_log}                                                                                        # at the first iteration, create log
else {log <- rbind(log,raw_log)}                                                                        # rbind glues currently loaded file (raw_log) to the end of dataframe log
}


####################     load     #######################################
####################   scr data   #######################################

# Load SCR files (Ledalab output exorted as .csv)
scrfiles <- list.files("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/9_SCR_export_preprocessed") 
setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/9_SCR_export_preprocessed")               # path to folder containing the scr files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 

for (subject in scrfiles){                                                                              # loop reading csv-file
raw_scr <- read.csv(subject)

if (subject=="ModERN_behav_01_SCR_Export_era_z.csv")
{scr <- raw_scr}                                                                                        # at the first iteration, create scr
else {scr <- rbind(scr,raw_scr)}                                                                        # rbind glues currently loaded file (raw_scr) to the end of dataframe scr
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
df4save <- data.frame()
master_GNG <- data.frame()
master_words <- data.frame()
master_scr <- data.frame()

for (number in unique(master$number)){
  
 singleID <- master[which(master$number==number), ]       # extract all rows for current subject from master data frame


 ####################       define      #####################################
 ####################      conditions   #####################################
 
 singleID <- mutate(singleID, 
                    valence = ifelse(w_type > 200, "pos", "neg"),                                                      
                    response_type = ifelse((t_resp == 43 | t_resp == 44), "FA", 
                                           ifelse(t_resp == 41, "FH", 
                                                  ifelse(t_resp == 45 | t_resp == 46, "CI", 
                                                         ifelse(t_resp == 42 , "SH",
                                                                "Miss_or_False_Key"))))
 )
 
 
 
 
 
 
 
 
  
  ####################     define      #####################################
  ####################  GNG  outliers   #####################################
  
  # assign TRUE to column GNG_outlier < 150 or > 500 ms (according to Pourtois)
  singleID <- mutate(singleID, outlier_GNG = (t_resp != 45 & t_resp != 46 & t_resp != 47) & t_rt < 150 | t_rt > 500) # create new column that has TRUE for each GNG outlier
  count_outlier_GNG <- length(singleID$outlier_GNG[singleID$outlier_GNG == TRUE])  # count number of outliers
 
  
  
  # here, GNG outliers are excluded only for GNG data frame, not for word categorization in corresponding trial (as in Pourtois)
  singleID_GNG <- singleID[(singleID$outlier_GNG == FALSE),]
  
  

  ####################   calculate  #####################################
  ####################    GNG RT    ##################################### 
  
  # FA = False Alarm; FH = Fast Hit; CI = Correctly Inhibited; SH = Slow Hit
  # RT not calculated for CI, because there is no rt
  GNG_mean_FA <- mean(singleID_GNG[singleID_GNG$response_type == "FA",]$t_rt)
  GNG_mean_FH <- mean(singleID_GNG[singleID_GNG$response_type == "FH",]$t_rt)
  GNG_mean_SH <- mean(singleID_GNG[singleID_GNG$response_type == "SH",]$t_rt)
  
  GNG_sd_FA <- sd(singleID_GNG[singleID_GNG$response_type == "FA",]$t_rt)
  GNG_sd_FH <- sd(singleID_GNG[singleID_GNG$response_type == "FH",]$t_rt)
  GNG_sd_SH <- sd(singleID_GNG[singleID_GNG$response_type == "SH",]$t_rt)
  
  GNG_cov_FA <- GNG_sd_FA/GNG_mean_FA * 100
  GNG_cov_FH <- GNG_sd_FH/GNG_mean_FH * 100
  GNG_cov_SH <- GNG_sd_SH/GNG_mean_SH * 100
  


  ####################   calculate GNG     #####################################
  ####################  % response types   ##################################### 
  
  GNG_count_FA <- length(singleID_GNG[singleID_GNG$response_type == "FA",]$t_rt)
  GNG_count_FH <- length(singleID_GNG[singleID_GNG$response_type == "FH",]$t_rt)
  GNG_count_CI <- length(singleID_GNG[singleID_GNG$response_type == "CI",]$t_rt)
  GNG_count_SH <- length(singleID_GNG[singleID_GNG$response_type == "SH",]$t_rt)
  GNG_count_miss <- length(singleID_GNG[singleID_GNG$t_resp == 47,]$t_rt)
  
  GNG_all_resp_without_outliers_wrong_keys <- GNG_count_FA + GNG_count_FH + GNG_count_CI + GNG_count_SH + GNG_count_miss 
  
  GNG_percent_FA <- (GNG_count_FA/GNG_all_resp_without_outliers_wrong_keys)*100
  GNG_percent_FH <- (GNG_count_FH/GNG_all_resp_without_outliers_wrong_keys)*100
  GNG_percent_CI <- (GNG_count_CI/GNG_all_resp_without_outliers_wrong_keys)*100
  GNG_percent_SH <- (GNG_count_SH/GNG_all_resp_without_outliers_wrong_keys)*100
  GNG_percent_miss <- (GNG_count_miss/GNG_all_resp_without_outliers_wrong_keys)*100
  
  

  
  ####################       define      #####################################
  ####################  word conditions  ##################################### 
  # FA = False Alarm; FH = Fast Hit; CI = Correctly Inhibited; SH = Slow Hit
  
  neg_after_FA <- subset(singleID, valence == "neg" & response_type == "FA" & w_resp == 51)
  pos_after_FA <- subset(singleID, valence == "pos" & response_type == "FA" & w_resp == 52)
  
  neg_after_FH <- subset(singleID, valence == "neg" & response_type == "FH" & w_resp == 51)
  pos_after_FH <- subset(singleID, valence == "pos" & response_type == "FH" & w_resp == 52)
  
  neg_after_CI <- subset(singleID, valence == "neg" & response_type == "CI" & w_resp == 51)
  pos_after_CI <- subset(singleID, valence == "pos" & response_type == "CI" & w_resp == 52)
  
  neg_after_SH <- subset(singleID, valence == "neg" & response_type == "SH" & w_resp == 51)
  pos_after_SH <- subset(singleID, valence == "pos" & response_type == "SH" & w_resp == 52)
  
  
 
   
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
  
  singleID <- mutate(singleID,
                     outlier_words = ifelse((valence == "neg" & w_resp == 51),
                                            ifelse(response_type == "FA", 
                                                   (abs(w_rt - median(neg_after_FA$w_rt))/mad(neg_after_FA$w_rt))>3,
                                            ifelse(response_type == "FH", 
                                                   (abs(w_rt - median(neg_after_FH$w_rt))/mad(neg_after_FH$w_rt))>3,
                                            ifelse(response_type == "CI", 
                                                   (abs(w_rt - median(neg_after_CI$w_rt))/mad(neg_after_CI$w_rt))>3,
                                                   (abs(w_rt - median(neg_after_SH$w_rt))/mad(neg_after_SH$w_rt))>3)
                                     )),
                                     # handle pos rows
                                     ifelse((valence == "pos" & w_resp == 52),
                                            ifelse(response_type == "FA", 
                                                   (abs(w_rt - median(pos_after_FA$w_rt))/mad(pos_after_FA$w_rt))>3,
                                            ifelse(response_type == "FH", 
                                                   (abs(w_rt - median(pos_after_FH$w_rt))/mad(pos_after_FH$w_rt))>3,
                                            ifelse(response_type == "CI", 
                                                   (abs(w_rt - median(pos_after_CI$w_rt))/mad(pos_after_CI$w_rt))>3,
                                            (abs(w_rt - median(pos_after_SH$w_rt))/mad(pos_after_SH$w_rt))>3))
                                    ), FALSE))
  )

  count_outlier_words_FA_FH <- length(singleID$outlier_words[singleID$outlier_words == TRUE & (singleID$response_type == "FA" | singleID$response_type == "FH" )])  # count number of outliers
  
  # exclude word outliers
  singleID_words <- singleID[(singleID$outlier_words == FALSE),]

  
  

  ####################    calculat priming    ##################################
  ####################    effect with mean    ##################################
  isNegative <- singleID_words$valence == "neg"
  isPositive <- singleID_words$valence == "pos"
  isFA <- singleID_words$response_type == "FA"
  isFH <- singleID_words$response_type == "FH"
  isCI <- singleID_words$response_type == "CI"
  isSH <- singleID_words$response_type == "SH"

  
  mean_neg_after_FA <- mean(singleID_words[isNegative & isFA & singleID_words$w_resp == 51,]$w_rt)
  mean_pos_after_FA <- mean(singleID_words[isPositive & isFA & singleID_words$w_resp == 52,]$w_rt)
  mean_neg_after_FH <- mean(singleID_words[isNegative & isFH & singleID_words$w_resp == 51,]$w_rt)
  mean_pos_after_FH <- mean(singleID_words[isPositive & isFH & singleID_words$w_resp == 52,]$w_rt)
  mean_neg_after_CI <- mean(singleID_words[isNegative & isCI & singleID_words$w_resp == 51,]$w_rt)
  mean_pos_after_CI <- mean(singleID_words[isPositive & isCI & singleID_words$w_resp == 52,]$w_rt)
  mean_neg_after_SH <- mean(singleID_words[isNegative & isSH & singleID_words$w_resp == 51,]$w_rt)
  mean_pos_after_SH <- mean(singleID_words[isPositive & isSH & singleID_words$w_resp == 52,]$w_rt)
  
  
  mean_priming_after_FA <- mean_pos_after_FA - mean_neg_after_FA
  mean_priming_after_FH <- mean_neg_after_FH - mean_pos_after_FH
  mean_priming_overall <- (mean_pos_after_FA + mean_neg_after_FH) - (mean_neg_after_FA + mean_pos_after_FH) 
  
  
  
  
  ####################    count number    ##################################
  ####################     of events      ##################################
  
  # counts number of events, after outliers are removed; if I want to work with median, shift this part before section "exclude outlier"
  
  count_neg_after_FA <- length(singleID_words[isNegative & isFA & singleID_words$w_resp == 51,]$w_rt)
  count_pos_after_FA <- length(singleID_words[isPositive & isFA & singleID_words$w_resp == 52,]$w_rt)
  count_neg_after_FH <- length(singleID_words[isNegative & isFH & singleID_words$w_resp == 51,]$w_rt)
  count_pos_after_FH <- length(singleID_words[isPositive & isFH & singleID_words$w_resp == 52,]$w_rt)
  count_neg_after_CI <- length(singleID_words[isNegative & isCI & singleID_words$w_resp == 51,]$w_rt)
  count_pos_after_CI <- length(singleID_words[isPositive & isCI & singleID_words$w_resp == 52,]$w_rt)
  count_neg_after_SH <- length(singleID_words[isNegative & isSH & singleID_words$w_resp == 51,]$w_rt)
  count_pos_after_SH <- length(singleID_words[isPositive & isSH & singleID_words$w_resp == 52,]$w_rt)
  
  count_incorr_neg_after_FA <- length(singleID_words[isNegative & isFA & singleID_words$w_resp == 53,]$w_rt)
  count_incorr_pos_after_FA <- length(singleID_words[isPositive & isFA & singleID_words$w_resp == 54,]$w_rt)
  count_incorr_neg_after_FH <- length(singleID_words[isNegative & isFH & singleID_words$w_resp == 53,]$w_rt)
  count_incorr_pos_after_FH <- length(singleID_words[isPositive & isFH & singleID_words$w_resp == 54,]$w_rt)
  count_incorr_neg_after_CI <- length(singleID_words[isNegative & isCI & singleID_words$w_resp == 53,]$w_rt)
  count_incorr_pos_after_CI <- length(singleID_words[isPositive & isCI & singleID_words$w_resp == 54,]$w_rt)
  count_incorr_neg_after_SH <- length(singleID_words[isNegative & isSH & singleID_words$w_resp == 53,]$w_rt)
  count_incorr_pos_after_SH <- length(singleID_words[isPositive & isSH & singleID_words$w_resp == 54,]$w_rt)
  
  count_miss_neg_after_FA <- length(singleID_words[isNegative & isFA & singleID_words$w_resp == 55,]$w_rt)
  count_miss_pos_after_FA <- length(singleID_words[isPositive & isFA & singleID_words$w_resp == 56,]$w_rt)
  count_miss_neg_after_FH <- length(singleID_words[isNegative & isFH & singleID_words$w_resp == 55,]$w_rt)
  count_miss_pos_after_FH <- length(singleID_words[isPositive & isFH & singleID_words$w_resp == 56,]$w_rt)
  count_miss_neg_after_CI <- length(singleID_words[isNegative & isCI & singleID_words$w_resp == 55,]$w_rt)
  count_miss_pos_after_CI <- length(singleID_words[isPositive & isCI & singleID_words$w_resp == 56,]$w_rt)
  count_miss_neg_after_SH <- length(singleID_words[isNegative & isSH & singleID_words$w_resp == 55,]$w_rt)
  count_miss_pos_after_SH <- length(singleID_words[isPositive & isSH & singleID_words$w_resp == 56,]$w_rt)
  
  
  # misses are contained in overall number of events so that errors, correct and misses sum up to 100 %
  count_all_neg_after_FA <- count_neg_after_FA + count_incorr_neg_after_FA + count_miss_neg_after_FA
  count_all_pos_after_FA <- count_pos_after_FA + count_incorr_pos_after_FA + count_miss_pos_after_FA
  count_all_neg_after_FH <- count_neg_after_FH + count_incorr_neg_after_FH + count_miss_neg_after_FH
  count_all_pos_after_FH <- count_pos_after_FH + count_incorr_pos_after_FH + count_miss_pos_after_FH
  count_all_neg_after_CI <- count_neg_after_CI + count_incorr_neg_after_CI + count_miss_neg_after_CI
  count_all_pos_after_CI <- count_pos_after_CI + count_incorr_pos_after_CI + count_miss_pos_after_CI
  count_all_neg_after_SH <- count_neg_after_SH + count_incorr_neg_after_SH + count_miss_neg_after_SH
  count_all_pos_after_SH <- count_pos_after_SH + count_incorr_pos_after_SH + count_miss_pos_after_SH
  
  
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
  
  
  
  ####################   calculate  #####################################
  ####################    SCR       #####################################   
  
 # Create master data frame for SCR analyses (for mean calculation excluding outlier in GNG only, exclude word outlier in single trial part, log-transform and standardize (z-score) SCR data)
 singleID_SCR <- singleID[(singleID$outlier_GNG == FALSE),]
 singleID_SCR$DDA.AMP.log <- log(singleID_SCR$DDA.AmpSum.TarResp+1)                              #log transform SCR 
 singleID_SCR$DDA.AMP.z_scorelog <- scale(singleID_SCR$DDA.AMP.log, center = TRUE, scale = TRUE) #log-transform SCR
 singleID_SCR$w_rt.z_score <- scale(singleID_SCR$w_rt, center = TRUE, scale = TRUE)              #z-transform word rt 

 
 # exclude trials followed by false alarm or wrong key in GNG or by incorrect word categorization or wrong key in word categorization
 singleID_SCR$exclude <- FALSE
 
 for (i in 1:(nrow(singleID_SCR)-1)) {
   current_row <- singleID_SCR[i,]
   next_row <- singleID_SCR[(i+1),]
   if ((next_row$t_resp %in% c(43,44,48,49) | next_row$w_resp %in% c(53,54,57,58))) {            
     singleID_SCR[i,]$exclude <- TRUE
   }     
 }
 singleID_SCR <- singleID_SCR[(singleID_SCR$exclude == FALSE),]
 
 
 # calculate mean SCR in conditions
 mean_SCR_after_FA <- mean(singleID_SCR[(singleID_SCR$t_resp == 43 | singleID_SCR$w_resp == 44) & (singleID_SCR$t_resp == 51 | singleID_SCR$w_resp == 52),]$DDA.AMP.z_scorelog)
 mean_SCR_after_FH <- mean(singleID_SCR[(singleID_SCR$t_resp == 41) & (singleID_SCR$t_resp == 51 | singleID_SCR$w_resp == 52),]$DDA.AMP.z_scorelog)
 mean_SCR_after_CI <- mean(singleID_SCR[(singleID_SCR$t_resp == 45 | singleID_SCR$w_resp == 46) & (singleID_SCR$t_resp == 51 | singleID_SCR$w_resp == 52),]$DDA.AMP.z_scorelog)
 mean_SCR_after_SH <- mean(singleID_SCR[(singleID_SCR$t_resp == 42) & (singleID_SCR$t_resp == 51 | singleID_SCR$w_resp == 52),]$DDA.AMP.z_scorelog)
 
 
 
 ####################     Add df in     #####################################
 ####################   in master df    #####################################
 
 df4save <- rbind(df4save, data.frame(number,mean_neg_after_FA,mean_pos_after_FA,mean_neg_after_FH,mean_pos_after_FH,mean_neg_after_CI,mean_pos_after_CI,mean_neg_after_SH,mean_pos_after_SH,mean_priming_overall,mean_priming_after_FA,mean_priming_after_FH,median_neg_after_FA,median_pos_after_FA,median_neg_after_FH,median_pos_after_FH,median_neg_after_CI,median_pos_after_CI,median_neg_after_SH,median_pos_after_SH,median_priming_overall,median_priming_after_FA,median_priming_after_FH,count_neg_after_FA,count_pos_after_FA,count_neg_after_FH,count_pos_after_FH,count_outlier_words_FA_FH,percent_correct_neg_after_FA,percent_correct_pos_after_FA,percent_correct_neg_after_FH,percent_correct_pos_after_FH,percent_correct_neg_after_CI,percent_correct_pos_after_CI,percent_correct_neg_after_SH,percent_correct_pos_after_SH,GNG_mean_FA,GNG_mean_FH,GNG_mean_SH,GNG_sd_FA,GNG_sd_FH,GNG_sd_SH,GNG_cov_FA,GNG_cov_FH,GNG_cov_SH,GNG_percent_FA,GNG_percent_FH,GNG_percent_CI,GNG_percent_SH,GNG_percent_miss, mean_SCR_after_FA,mean_SCR_after_FH,mean_SCR_after_CI,mean_SCR_after_SH))
 master_GNG <- rbind(master_GNG,singleID_GNG)         # excludes only GNG outlier
 master_words <- rbind(master_words,singleID_words)   # excludes only word outlier
 master_scr <- rbind(master_scr,singleID_SCR)         # exclude outlier GNG and trials with error in word or followed by error in GNG; exclude word outlier later, to keep these trials here for Analyses only on SCR, without Priming
}






 ####################     save df4save    ##################################
 ####################    as excel file    ##################################
 
 
 setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses")    # setting a different folder as working directory to prevent saving stuff into the folder containing the logfiles
 
 date_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
 
 filename <- paste("SummaryStatisticsSingleTrial_For_",length(logfiles),"_subjects_",date_time, ".xlsx", sep = "")
 
 write.xlsx(df4save, filename)
 #save(df4save, file = filename)
 
 
 
 
 
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
   ggtitle("Mean Response Time Word Categorization") +                                                          # add title
   xlab("Previous Response") + ylab("Reaction Time in ms") +                                                    # label axes
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
   ggtitle("Mean Accuracy Word Categorization") +                                                               # add title
   xlab("Previous Response") + ylab("Accuracy in %") +                                                          # label axes
   theme(plot.title = element_text(hjust = 0.5)) +                                                              # center title
   scale_fill_manual(values=c("mediumblue", "limegreen"))                                                       # change bar colors
 
 
 
 
 
 #### till here, results of both scripts (this and "Analysis_evaluative_GNG.R" are identical!!!!
 
 
 # plot rt all conditions
 df4plotRT <- data.frame(
   response = c("false alarm","false alarm","fast hit","fast hit","correctly inhibited","correctly inhibited","slow hit", "slow hit"),
   valence = c("neg","pos","neg","pos","neg","pos","neg","pos"),
   conditions = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH","neg_after_CI","pos_after_CI","neg_after_SH","pos_after_SH"),
   mean = as.numeric(c(descriptive_statistics[2,2],descriptive_statistics[2,3],descriptive_statistics[2,4],descriptive_statistics[2,5],descriptive_statistics[2,6],descriptive_statistics[2,7],descriptive_statistics[2,8],descriptive_statistics[2,9])),
   se = as.numeric(c(descriptive_statistics[3,2],descriptive_statistics[3,3],descriptive_statistics[3,4],descriptive_statistics[3,5],descriptive_statistics[3,6],descriptive_statistics[3,7],descriptive_statistics[3,8],descriptive_statistics[3,9])))
 
 df4plotRT$response <- factor(df4plotRT$response, levels = c("false alarm","fast hit","correctly inhibited","slow hit"))
 
 ggplot(df4plotRT,x = response,y = mean, aes(response, mean, fill = valence))+
   geom_bar(stat="identity", position=position_dodge()) +                                                       # add bars, based on stats values; dodge to avoid stacked bars
   geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.95), width=0.1) +    # add error bars
   ggtitle("Mean Response Time Word Categorization") +                                                          # add title
   xlab("Previous Response") + ylab("Reaction Time in ms") +                                                    # label axes
   theme(plot.title = element_text(hjust = 0.5)) +                                                              # center title
   scale_fill_manual(values=c("mediumblue", "limegreen"))                                                       # change bar colors
 
 
 # plot accuracy all conditions
 df4plotACC <- data.frame(
   response = c("false alarm","false alarm","fast hit","fast hit","correctly inhibited","correctly inhibited","slow hit", "slow hit"),
   valence = c("neg","pos","neg","pos","neg","pos","neg","pos"),
   conditions = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH","neg_after_CI","pos_after_CI","neg_after_SH","pos_after_SH"),
   mean = as.numeric(c(descriptive_statistics[2,29],descriptive_statistics[2,30],descriptive_statistics[2,31],descriptive_statistics[2,32],descriptive_statistics[2,33],descriptive_statistics[2,34],descriptive_statistics[2,35],descriptive_statistics[2,36])),
   se = as.numeric(c(descriptive_statistics[3,29],descriptive_statistics[3,30],descriptive_statistics[3,31],descriptive_statistics[3,32],descriptive_statistics[3,33],descriptive_statistics[3,34],descriptive_statistics[3,35],descriptive_statistics[3,36])))
 
 df4plotACC$response <- factor(df4plotACC$response, levels = c("false alarm","fast hit","correctly inhibited","slow hit"))
 
 ggplot(df4plotACC,x = response,y = mean, aes(response, mean, fill = valence))+
   geom_bar(stat="identity", position=position_dodge()) +                                                       # add bars, based on stats values; dodge to avoid stacked bars
   geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.95), width=0.1) +    # add error bars
   ggtitle("Mean Accuracy Word Categorization") +                                                               # add title
   xlab("Previous Response") + ylab("Accuracy in %") +                                                          # label axes
   theme(plot.title = element_text(hjust = 0.5)) +                                                              # center title
   scale_fill_manual(values=c("mediumblue", "limegreen"))                                                       # change bar colors
 
 
 # plot SCR all conditions
 df4plotSCR <- data.frame(
   response = c("false alarm","fast hit","correctly inhibited","slow hit"),
   valence = c("neg","pos","neg","pos"),
   conditions = c("mean_SCR_after_FA","mean_SCR_after_FH","mean_SCR_after_CI","mean_SCR_after_SH"),
   mean = as.numeric(c(descriptive_statistics[2,51],descriptive_statistics[2,52],descriptive_statistics[2,53],descriptive_statistics[2,54])),
   se = as.numeric(c(descriptive_statistics[3,51],descriptive_statistics[3,52],descriptive_statistics[3,53],descriptive_statistics[3,54])))
 
 df4plotSCR$response <- factor(df4plotSCR$response, levels = c("false alarm","fast hit","correctly inhibited","slow hit"))
 
 ggplot(df4plotSCR,x = response,y = mean, aes(response, mean))+
   geom_bar(stat="identity", position=position_dodge(), fill = "mediumblue") +                                  # add bars, based on stats values; dodge to avoid stacked bars; change bar color
   geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.95), width=0.1) +    # add error bars
   ggtitle("Mean SCR") +                                                                                        # add title
   xlab("Response type") + ylab("log SCR (z-score) in microsiemens") +                                          # label axes
   theme(plot.title = element_text(hjust = 0.5))                                                                # center title
 
 
 
 
 
 
 
 #######################################################################################
 #####################  Aggregate Data and calculate ANOVAs   ##########################
 #######################################################################################
 
 
 options(contrasts=c("contr.sum","contr.poly")) # to adhere to the sum-to-zero convention for effect weights, always do this before running ANOVAs in R. This matters sometimes (not always). If I don’t do it, the sum of squares calculations may not match what I get e.g. in SPSS
 
 
 
 #####################       ANOVA      ####################################
 #####################       Words      #################################### 
 
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
                     detailed = TRUE
 )
 
 # calculate ANOVA with ezANOVA for accuracy -> gives same result as SPSS :)
 anova_accuracy <- ezANOVA(data = df4anova, 
                           dv = .(accuracy), 
                           wid = .(subject), 
                           within = .(response_type, word_valence), 
                           detailed = TRUE)
 
 
 # post hoc comparison: when calculating within-subject ANOVAs, no direct post hoc test is possible; I can only run paitwise t-tests with corresponding adjustment (e.g. Holm, which is better than Bonferroni)
 anova_rt_posthoc_response_type <- pairwise.t.test(df4anova$rt,df4anova$response_type,p.adjust.method="holm")
 anova_rt_posthoc_word_valence <- pairwise.t.test(df4anova$rt,df4anova$word_valence,p.adjust.method="holm") 
 anova_rt_posthoc_priming <- pairwise.t.test(df4anova$rt,df4anova$condition,p.adjust.method="holm") 
 
 anova_accuracy_posthoc_response_type <- pairwise.t.test(df4anova$accuracy,df4anova$response_type,p.adjust.method="holm")
 anova_accuracy_posthoc_word_valence <- pairwise.t.test(df4anova$accuracy,df4anova$word_valence,p.adjust.method="holm") 
 anova_accuracy_posthoc_priming <- pairwise.t.test(df4anova$accuracy,df4anova$condition,p.adjust.method="holm") 
 
 
 
 
 
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
 
 
 
 #####################    ANOVA SCR and   ####################################
 #####################     GNG percent    #################################### 
 
 
 # ANOVA requires several rows for each subject, each one row per factor: here I reorder data by reshaping existing data frame
 df4anova_GNG_scr <-   reshape(data = df4save[,c(1,46:49,51:54)], 
                               direction = "long",
                               varying = list(c("GNG_percent_FA","GNG_percent_FH","GNG_percent_CI","GNG_percent_SH"),c("mean_SCR_after_FA","mean_SCR_after_FH","mean_SCR_after_CI","mean_SCR_after_SH")),
                               v.names = c("percent", "SCR"),
                               idvar = "subject",
                               timevar = "condition",
                               times = c("FA","FH","CI","SH")
 )
 
 row.names(df4anova_GNG_scr) <- NULL
 df4anova_GNG_scr <- df4anova_GNG_scr[sort.list(df4anova_GNG_scr$subject),]                        # sort df by subject
 
 df4anova_GNG_scr$subject <- factor(df4anova_GNG_scr$subject)                                          # create factors
 df4anova_GNG_scr$condition <- factor(df4anova_GNG_scr$condition)
 
 
 
 
 # calculate ANOVA with ezANOVA for percent responses GNG
 anova_GNG_percent <- ezANOVA(data = df4anova_GNG_scr, 
                              dv = .(percent), 
                              wid = .(subject), 
                              within = .(condition), 
                              detailed = TRUE)
 
 
 # post hoc comparison: when calculating within-subject ANOVAs, no direct post hoc test is possible; I can only run paitwise t-tests with corresponding adjustment (e.g. Holm, which is better than Bonferroni)
 anova_GNG_percent_posthoc_response_type <- pairwise.t.test(df4anova_GNG_scr$percent,df4anova_GNG_scr$condition,p.adjust.method="holm")
 
 
 
 
 #only execute this because SCR of last subjects not exported yet
 df4anova_GNG_scr <- df4anova_GNG_scr[df4anova_GNG_scr$number < 20,]
 df4anova_GNG_scr$subject <- factor(df4anova_GNG_scr$subject)
 
 
 # calculate ANOVA with ezANOVA for SCR
 anova_GNG_SCR <- ezANOVA(data = df4anova_GNG_scr, 
                          dv = .(SCR), 
                          wid = .(subject), 
                          within = .(condition), 
                          detailed = TRUE)
 
 
 
 # post hoc comparison: when calculating within-subject ANOVAs, no direct post hoc test is possible; I can only run paitwise t-tests with corresponding adjustment (e.g. Holm, which is better than Bonferroni)
 anova_GNG_SCR_posthoc_response_type <- pairwise.t.test(df4anova_GNG_scr$SCR,df4anova_GNG_scr$condition,p.adjust.method="holm")
 

 
 
 ########################################################################################
 #######################   Correlations with    #########################################
 #######################        Traits          #########################################
 ########################################################################################
 
 # Read in questionnaire data, merge with aggregated data, create correlation matrix
 setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study")
 questionnaires <- list.files(pattern = ".sav")                                                         # make sure that only one .sav file (= the current PEQ export) exists there!
 questionnaires <- read.spss(questionnaires, to.data.frame = TRUE)
 questionnaires <- questionnaires[, colSums(is.na(questionnaires)) != nrow(questionnaires)]             # remove columns with only NA
 questionnaires <- questionnaires[,c("CODE", "BD2SUMT0","BASO00T0","BASO01T0","BASO02T0","BASO03T0","BASO04T0","FMPO00T0","FMPO01T0","FMPO02T0","FMPO03T0","FMPO04T0","NNGO00T0","NNGO01T0","OCISUMT0", "OCIO00T0", "OCIO01T0","OCIO02T0","OCIO03T0","OCIO04T0","OCIO05T0","PANO00T0","PANO01T0","PSWSUMT0","STSSUMT0","STTSUMT0","TCISUMT0","TCIO00T0","TCIO01T0","TCIO02T0","TCIO03T0","WSTSUMT0")] # only keep relevant columns (I chose to keep sum scores instead of mean scores, because for BDI and OCi only sum scores are available and missing data are not possible because data aquired with tablet)
 questionnaires <- questionnaires[grep("behav", questionnaires$CODE), ]                                 # only keep participants from behavioral study
 questionnaires$CODE <- as.numeric(substr(questionnaires$CODE,14,15))                                   # convert code to number for joining
 colnames(questionnaires) <-c("CODE", "BDI_Total", "BIS_Total", "BAS_Total", "BAS_Drive", "BAS_Fun_Seeking", "BAS_Reward_Responsiveness", "FMPS_CMD", "FMPS_PEC", "FMPS_PST", "FMPS_ORG", "FMPS_PER/Total", "NEO_Neuroticism", "NEO_Conscientiousness", "OCI_Total", "OCI_Washing", "OCI_Checking", "OCI_Ordering", "OCI_Obsessions", "OCI_Hoarding", "OCI_Neutralising", "PANAS_Pos", "PANAS_Neg", "PSWQ_Total", "STAI_State", "STAI_Trait", "TCI_Total", "TCI_Anticipatory_worry", "TCI_Fear_of_uncertainty", "TCI_Shyness", "TCI_Fatigability", "WST_Total")   # rename columns 
 
 
 df4correlation <- left_join(df4save,questionnaires, by = c("number" = "CODE"))
 df4correlation <- df4correlation[,c("number","mean_priming_overall", "mean_priming_after_FA", "mean_priming_after_FH","mean_SCR_after_FA", "BDI_Total", "BIS_Total", "BAS_Total", "BAS_Drive", "BAS_Fun_Seeking", "BAS_Reward_Responsiveness", "FMPS_CMD", "FMPS_PEC", "FMPS_PST", "FMPS_ORG", "FMPS_PER/Total", "NEO_Neuroticism", "NEO_Conscientiousness", "OCI_Total", "OCI_Washing", "OCI_Checking", "OCI_Ordering", "OCI_Obsessions", "OCI_Hoarding", "OCI_Neutralising", "PANAS_Pos", "PANAS_Neg", "PSWQ_Total", "STAI_State", "STAI_Trait", "TCI_Total", "TCI_Anticipatory_worry", "TCI_Fear_of_uncertainty", "TCI_Shyness", "TCI_Fatigability")]   # only keep relevant columns 
 
 
 
 # When data collection is complete, use shapiro.test to test for normality; if no normality, calculate Kendall's Tau instead of Pearson
 correlation_matrix <- corr.test(df4correlation, adjust = "none")
 p_values <- correlation_matrix[["p"]]
 r_values <- correlation_matrix[["r"]]
 p_values_subset <- data.frame(p_values[2:5,5:(ncol(p_values))])
 r_values_subset <- data.frame(r_values[2:5,5:(ncol(r_values))])
 
 
 
 
 # Template for plotting significant correlations
 ggplot(df4correlation, aes (x= mean_priming_after_FA, y = TCI_Total)) +
   geom_point(color="mediumblue", size=3) +
   geom_smooth(method=lm, se=FALSE, color="limegreen", size = 1) +
   ggtitle("Zusammenhang XXXX und PYYYY") +
   labs(x="Prim", y = "YYYY") +
   theme(axis.title = element_text(size = 15, hjust = 0.5)) +
   theme(axis.text=element_text(size=12)) +
   theme(legend.title=element_text(size=15)) +
   theme(legend.text=element_text(size=12)) +
   theme(plot.title = element_text(size = 17, hjust = 0.5))
 
 
 
 
 
 # #####################    linear mixed models    ####################################
 # #####################      words (mean RT)      #################################### 
 # 
 # # prepare master_words for single trial analysis: exclude incorrect word responses, misses, wrong keys; create column condition; turn relevant variables in factors
 # master_words <- subset(master_words,w_resp <= 52 & t_resp <= 46)
 # master_words$condition <- paste(master_words$valence,"_after_",master_words$response_type)
 # master_words$response_type <- factor(master_words$response_type, levels=c("FA","FH","CI","SH"))   # reorder factor levels! (this is important, because first level is used as reference factor in LMM)
 # master_words$valence <- factor(master_words$valence)                                              # automatic factor order is alphabetical (1 = neg, 2 = pos)
 # master_words$number <- factor(master_words$number)  
 # 
 # # set contrasts
 # contrasts(master_words$response_type) <- contr.sdif(4) 
 # contrasts(master_words$valence) <- contr.sdif(2) 
 # 
 # 
 # 
 # # calculate LMM for rt and plot
 # emm_options(pbkrtest.limit = 9400)
 # emm_options(lmerTest.limit = 9400) # due to warning D.f. calculations have been disabled because the number of observations exceeds 3000.To enable adjustments, set emm_options(pbkrtest.limit = 5866)
 # 
 # LMM_rt <- lmer(w_rt ~ master_words$response_type * valence + (1|number), data=master_words, REML = FALSE)       
 # anova(LMM_rt)                                                                                     # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 # results_LMM_rt <- analyze(LMM_rt)                                                                 # print results
 # print(results_LMM_rt)
 # results_LMM_rt <- get_contrasts(LMM_rt, "response_type * valence")                                # Provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
 # print(results_LMM_rt$contrasts)                                                                   # print contrasts
 # print(results_LMM_rt$means)                                                                       # investigate means
 # 
 # ggplot(results_LMM_rt$means, aes(x=response_type, y=Mean, color=valence, group=valence)) +        # plot
 #   geom_line(position = position_dodge(.3)) +
 #   geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher),
 #                   position = position_dodge(.3)) +
 #   ylab("Response Time in ms") +
 #   xlab("Response Type") +
 #   theme_bw()
 # 
 # 
 # # check normality of residuals
 # qqnorm(resid(LMM_rt))
 # qqline(resid(LMM_rt))
 # shapiro.test(resid(LMM_rt))
 # 
 # 
 # 
 # # post hoc tests for resolving main effects and interaction effects 
 # summary(glht(LMM_rt, linfct=mcp(response_type = "Tukey")), test = adjusted("holm"))  
 # summary(glht(LMM_rt, linfct=mcp(valence = "Tukey")), test = adjusted("holm"))                   # same result can be obtained by lmer_rt_posthoc_response_type <- emmeans(model_lmer_rt, ~ word_valence);;pairs(lmer_rt_posthoc_response_type)
 # LMM_rt_posthoc_priming <- emmeans(LMM_rt, ~ valence|response_type, adjust="tukey")              # compare word valence within each level of response_type
 # pairs(LMM_rt_posthoc_priming) 
 # 
 # 
 # 
 # 
 # 
 # #####################    linear mixed models    ####################################
 # #####################     words (Accuracy)      #################################### 
 # 
 # # accumulate data for LMM accuracy
 # df4accuracy <-   reshape(data = df4save[,c(1:9,29:36)], 
 #                       direction = "long",
 #                       varying = list(c("mean_neg_after_FA","mean_pos_after_FA","mean_neg_after_FH","mean_pos_after_FH","mean_neg_after_CI","mean_pos_after_CI" ,"mean_neg_after_SH","mean_pos_after_SH"),c("percent_correct_neg_after_FA", "percent_correct_pos_after_FA","percent_correct_neg_after_FH","percent_correct_pos_after_FH","percent_correct_neg_after_CI","percent_correct_pos_after_CI","percent_correct_neg_after_SH","percent_correct_pos_after_SH")),
 #                       v.names = c("rt","accuracy"),
 #                       idvar = "number",
 #                       timevar = "condition",
 #                       times = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH","neg_after_CI","pos_after_CI","neg_after_SH","pos_after_SH")
 # )
 # 
 # row.names(df4accuracy) <- NULL
 # df4accuracy <- df4accuracy[sort.list(df4accuracy$number),]                                      # sort df by subject
 # 
 # df4accuracy$response_type <- substr(df4accuracy$condition, 11, 12)                              # create columns needed as factors           
 # df4accuracy$valence <- substr(df4accuracy$condition, 1, 3)
 # 
 # df4accuracy$number <- factor(df4accuracy$number)                                 
 # df4accuracy$response_type <- factor(df4accuracy$response_type, levels=c("FA","FH","CI","SH"))   # create factors response_type and reorder factor levels! (this is important, because first level is used as reference factor in LMM)
 # df4accuracy$valence <- factor(df4accuracy$valence)                                              # automatic factor order is alphabetical (1 = neg, 2 = pos)
 # 
 # 
 # # set contrasts
 # contrasts(df4accuracy$response_type) <- contr.sdif(4) 
 # contrasts(df4accuracy$valence) <- contr.sdif(2) 
 # 
 # 
 # # calculate LMM for accuracy and plot  
 # LMM_accuracy <- lmer(accuracy ~ response_type * valence + (1|number), data=df4accuracy, REML = FALSE)         
 # anova(LMM_accuracy)                                                                             # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 # results_LMM_accuracy <- analyze(LMM_accuracy)                                                   # print results
 # print(results_LMM_accuracy)
 # results_LMM_accuracy <- get_contrasts(LMM_accuracy, "response_type * valence")                  # Provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
 # print(results_LMM_accuracy$contrasts)                                                           # print contrasts
 # print(results_LMM_accuracy$means)                                                               # investigate means
 # 
 # ggplot(results_LMM_accuracy$means, aes(x=response_type, y=Mean, color=valence, group=valence)) +# plot 
 #   geom_line(position = position_dodge(.3)) +
 #   geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher), 
 #                   position = position_dodge(.3)) +
 #   ylab("Accuracy in %") +
 #   xlab("Response Type") +
 #   theme_bw() 
 # 
 # 
 # 
 # # check normality of residuals
 # qqnorm(resid(LMM_accuracy))
 # qqline(resid(LMM_accuracy))
 # shapiro.test(resid(LMM_accuracy))
 # 
 # 
 # 
 # # post hoc tests for resolving main effects and interaction effects 
 # summary(glht(LMM_accuracy, linfct=mcp(response_type = "Tukey")), test = adjusted("holm"))  
 # summary(glht(LMM_accuracy, linfct=mcp(valence = "Tukey")), test = adjusted("holm"))
 # LMM_accuracy_posthoc_priming <- emmeans(LMM_accuracy, ~ valence|response_type, adjust="tukey")
 # pairs(LMM_accuracy_posthoc_priming) 
 # 
 # 
 # 
 # 
 # 
 # 
 # #####################    linear mixed models    ####################################
 # #####################          GNG (RT)         #################################### 
 # 
 # # prepare master_GNG for single trial analysis: ecxlude incorrect responses, misses, wrong keys; turn relevant variables in factors; for LMM_GNG_rt also remove CI
 # master_GNG$response_type <- factor(master_GNG$response_type, levels=c("FA","FH","CI","SH")) # reorder factor levels! (this is important, because first level is used as reference factor in LMM)
 # master_GNG$number <- factor(master_GNG$number)  
 # master_GNG_rt <- subset(master_GNG,t_resp <= 44)
 # 
 # # set contrasts
 # contrasts(master_GNG$response_type) <- contr.sdif(4) 
 # 
 # 
 # 
 # # calculate LMM for rt and plot
 # LMM_GNG_rt <- lmer(t_rt ~ response_type + (1|number), data=master_GNG_rt, REML = FALSE)                     
 # anova(LMM_GNG_rt)                                                                            # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 # results_LMM_GNG_rt <- analyze(LMM_GNG_rt)                                                    # print results
 # print(results_LMM_GNG_rt)
 # results_LMM_GNG_rt <- get_contrasts(LMM_GNG_rt, "response_type")                             # Provide the model and the factors to contrast;; by default, get_contrasts uses the Tukey method for p value adjustment; add ,adjust="none" to turn off automatic p value correction after Tucky
 # print(results_LMM_GNG_rt$contrasts)                                                          # print contrasts
 # print(results_LMM_GNG_rt$means)                                                              # investigate means
 # 
 # ggplot(results_LMM_GNG_rt$means, aes(x=response_type, y=Mean, color=response_type)) +        # plot 
 #   geom_line(position = position_dodge(.3)) +
 #   geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher), 
 #                   position = position_dodge(.3)) +
 #   ylab("Response Time in ms") +
 #   xlab("Response Type") +
 #   theme_bw()
 # 
 # 
 # # post hoc tests for resolving main effect 
 # summary(glht(LMM_GNG_rt, linfct=mcp(response_type = "Tukey")), test = adjusted("holm"))  
 # 
 # 
 # 
 # 
 # 
 # # accumulate data for LMM accuracy
 # df4GNG_percent <-   reshape(data = df4save[,c(1,46:49)], 
 #                                   direction = "long",
 #                                   varying = list(c("GNG_percent_FA","GNG_percent_FH","GNG_percent_CI","GNG_percent_SH")),
 #                                   v.names = c("percent"),
 #                                   idvar = "number",
 #                                   timevar = "condition",
 #                                   times = c("FA","FH","CI","SH")
 # )
 # 
 # row.names(df4GNG_percent) <- NULL
 # df4GNG_percent <- df4GNG_percent[sort.list(df4GNG_percent$number),]                            # sort df by subject
 # 
 # df4GNG_percent$number <- factor(df4GNG_percent$number)                                         # create factors
 # df4GNG_percent$condition <- factor(df4GNG_percent$condition)
 # 
 # 
 # # set contrasts
 # contrasts(df4GNG_percent$condition) <- contr.sdif(4) 
 # 
 # 
 # 
 # # calculate LMM for percent and plot
 # LMM_GNG_percent <- lmer(percent ~ condition + (1|number), data=df4GNG_percent, REML = FALSE)                
 # anova(LMM_GNG_percent)                                                                         # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 # results_LMM_GNG_percent <- analyze(LMM_GNG_percent)                                            # print results
 # print(results_LMM_GNG_percent)
 # results_LMM_GNG_percent <- get_contrasts(LMM_GNG_percent, "condition")                         # Provide the model and the factors to contrast;; by default, get_contrasts uses the Tukey method for p value adjustment; add ,adjust="none" to turn off automatic p value correction after Tucky
 # print(results_LMM_GNG_percent$contrasts)                                                       # print contrasts
 # print(results_LMM_GNG_percent$means)                                                           # investigate means
 # 
 # ggplot(results_LMM_GNG_percent$means, aes(x=condition, y=Mean, color=condition)) +             # plot
 #   geom_line(position = position_dodge(.3)) +
 #   geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher),
 #                   position = position_dodge(.3)) +
 #   ylab("Frequency in %") +
 #   xlab("Response Type") +
 #   theme_bw()
 # 
 # 
 # 
 # # post hoc tests for resolving main effects 
 # summary(glht(LMM_GNG_percent, linfct=mcp(condition = "Tukey")), test = adjusted("holm")) 
 # 
 # 
 # 
 # 
 # 
 # ########################################################################
 # ##############################SCR#######################################
 # ########################################################################
 # ########################################################################
 # 
 # # exclude CI and SH, exclude trials with outlier or error in word or miss/wrong key in GNG 
 # master_scr <- subset(master_scr,w_resp <= 52 & t_resp <= 46 & outlier_words == FALSE & (t_resp == 41 | t_resp == 43 | t_resp == 44))
 # 
 # 
 # # log transform w_rt
 # master_scr$w_rt_log <- log(master_scr$w_rt)
 # 
 # 
 # # set contrasts for fixed effects
 # master_scr$response_type <- as.factor(master_scr$response_type)
 # master_scr$valence <- as.factor(master_scr$valence)
 # contrasts(master_scr$response_type) <- contr.sdif(2) 
 # contrasts(master_scr$valence) <- contr.sdif(2) 
 # 
 # 
 # 
 # # calculate LMM for SCR and plot (does not converge)
 # model_lmer_SCR <- lmer(DDA.AMP.log ~ response_type * valence * w_rt_log + (1+ response_type * valence * w_rt_log|number)  + (w_rt_log|word), data=master_scr, REML = FALSE)       
 # anova(model_lmer_SCR)                                                                          # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 # results_model_lmer_SCR <- analyze(model_lmer_SCR)                                              # print results
 # print(results_model_lmer_SCR)
 # results_model_lmer_SCR <- get_contrasts(model_lmer_SCR, "response_type * valence")             # Provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
 # print(results_model_lmer_SCR$contrasts)                                                        # print contrasts
 # print(results_model_lmer_SCR$means)                                                            # investigate means
 # 
 # 
 # 
 # # try if SCR to errors is predicted by w_rt_log
 # master_scr_FA <- subset(master_scr, t_resp == 43 | t_resp == 44)
 # 
 # # set contrasts
 # contrasts(master_scr_FA$valence) <- contr.sdif(2) 
 # 
 # model_lmer_SCR_FA <- lmer(DDA.AMP.log ~ valence * w_rt_log + (1|number), data=master_scr_FA, REML = FALSE)       
 # anova(model_lmer_SCR_FA)                                                                      # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 # results_model_lmer_SCR_FA <- analyze(model_lmer_SCR_FA)                                       # print results
 # print(results_model_lmer_SCR_FA)
 # results_model_lmer_SCR_FA <- get_contrasts(model_lmer_SCR_FA, "valence*w_rt_log")             # Provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
 # print(results_model_lmer_SCR_FA$contrasts)                                                    # print contrasts
 # print(results_model_lmer_SCR_FA$means)   
 # 
 # 
 # 
 # 
 # # try to calculate correlation (use log transformed and z standardized SCR and log transformed rt)
 # master_scr_FA_neg <- subset(master_scr, (t_resp == 43 | t_resp == 44) & w_resp == 51)
 # master_scr_FA_neg$w_rt_log <- scale(master_scr_FA_neg$w_rt_log, center = TRUE, scale = TRUE)
 # cor.test(master_scr_FA_neg$DDA.AMP.z_scorelog,master_scr_FA_neg$w_rt_log)
 # 
 # ggplot(master_scr_FA_neg, aes (x= DDA.AMP.z_scorelog, y = w_rt_log)) +
 #   geom_point(color="mediumblue", size=3) +
 #   geom_smooth(method=lm, se=FALSE, color="limegreen", size = 1) +
 #   ggtitle("Zusammenhang SCR und w_rt_log") +
 #   labs(x="SCR", y = "w_rt_log") +
 #   theme(axis.title = element_text(size = 15, hjust = 0.5)) + 
 #   theme(axis.text=element_text(size=12)) +
 #   theme(legend.title=element_text(size=15)) +
 #   theme(legend.text=element_text(size=12)) +
 #   theme(plot.title = element_text(size = 17, hjust = 0.5)) 
 # 
 # 
 # 
 # 
 # # Check whether SCR differs between FA and FH calculate LMM for SCR and plot 
 # model_lmer_SCR <- lmer(DDA.AMP.log ~ response_type + (1|number), data=master_scr, REML = FALSE)      
 # anova(model_lmer_SCR)                                                                         # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 # results_model_lmer_SCR <- analyze(model_lmer_SCR)                                             # print results
 # print(results_model_lmer_SCR)
 # results_model_lmer_SCR <- get_contrasts(model_lmer_SCR, "response_type")                      # Provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
 # print(results_model_lmer_SCR$contrasts)                                                       # print contrasts
 # print(results_model_lmer_SCR$means)                                                           # investigate means
 # 
 # 
 # 
 # ggplot(results_model_lmer_SCR$means, aes(x=response_type, y=Mean)) +    # plot
 #   geom_line(position = position_dodge(.3)) +
 #   geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher),
 #                   position = position_dodge(.3)) +
 #   ylab("SCR in log(Microsiemens)") +
 #   xlab("Response Type") +
 #   theme_bw()
 # 
 # 
 # summary(glht(model_lmer_SCR, linfct=mcp(response_type = "Tukey")), test = adjusted("holm"))  
 # 
 # 
 # 
 # 

 
 
 
 
 
