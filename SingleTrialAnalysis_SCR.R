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
library(ez)
#

# clear environment
rm(list=ls())

# force R to not use exponential notation
options(scipen = 999)

# create empty data frames to write aggregated data and single trial data in it
data4mixedmodels <- data.frame()   
df4save <- data.frame()   



####################     load log    #######################################
####################  and scr data   #######################################

# Load logfiles
logfiles <- list.files("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_behav")       

for (subject in logfiles){                                                                              # loop reading txt-files as table, ommit first 58 lines and added lines after trials; use first line as header
  setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_behav")                      # path to folder containing the log files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 
  raw_log <- (read.table(subject, skip = 58, fill = TRUE, header = TRUE, nrows = 516))
  name_of_subject <-  gsub("\\.txt*","",subject)

  
# Load SCR files 
setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/9_SCR_export_preprocessed")               # path to folder containing the scr files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 
filename_scr <- paste0(name_of_subject,"_SCR_Export_era.txt")
raw_scr <- read.table(filename_scr, header = TRUE, sep=",")






####################     merge log and     #######################################
####################       scr data        #######################################

# return all rows of log, merged with matching columns of scr; rows in log with no match in scr will have NAs in the new columns; if there are multiple matches between log and scr, all combinations of the matches are returned
data_log_scr <- left_join(raw_log,raw_scr,by = c("number"= "Participant.ID","count" = "Trial.ID","t_type"="Event.ID.Tar","t_resp"="Event.ID.TarResp","w_type"="Event.ID.Word","w_resp"="Event.ID.WordResp"))
message('Processing ', subject, ' of ', length(logfiles))

# check whether there were rows in log with no match in scr
if (any(is.na(data_log_scr$DDA.AreaSum.TarResp))){
  print("There are rows in log with no match in SCR! Should be only caused by missing SCR trials in subjects 25 (trial 112-137) and 29 (trial 185.5-200)!")
}else{
    print("Matching of logfile and SCR-file is ok. Each row in log has a matching row in SCR")
  }


# check whether there were there are multiple matches between log and scr
if (any(duplicated(data_log_scr[,c('number','count')]))){
  stop("There are rows in log with multiple matches in SCR!")
}else{
  print("Matching of logfile and SCR-file is ok. No row in log has several matching row in SCR")
}




####################   clean merged   ############################
####################    data frame    ############################ 

# select relevant columns
single_trial_data <- subset(data_log_scr, select = c("name","number","count","t_type","t_resp","t_rt","w_type","w_resp","w_rt","word","DDA.AmpSum.TarResp","DDA.AreaSum.TarResp"))


# rename columns
single_trial_data <- rename(single_trial_data,
                            subject = name,
                            subjectID = number,
                            trial = count,
                            gng_condition = t_type,
                            gng_resp = t_resp,
                            gng_rt = t_rt,
                            word_condition = w_type,
                            word_resp = w_resp,
                            word_rt = w_rt,
                            scr_amplitude = DDA.AmpSum.TarResp,
                            scr_area = DDA.AreaSum.TarResp)



####################       define      #####################################
####################      conditions   #####################################

single_trial_data <- mutate(single_trial_data, 
                   valence = ifelse(word_condition > 200, "pos", "neg"),                                                      
                   response_type = ifelse((gng_resp == 43 | gng_resp == 44), "FA", 
                                          ifelse(gng_resp == 41, "FH", 
                                                 ifelse(gng_resp == 45 | gng_resp == 46, "CI", 
                                                        ifelse(gng_resp == 42 , "SH",
                                                               "Miss_or_False_Key"))))
)



####################     define      #####################################
####################  GNG  outliers   #####################################   

# assign TRUE to column GNG_outlier < 150 or > 500 ms (according to Pourtois); here, GNG outliers are excluded only for GNG analyses, not for word categorization in corresponding trial (as in Pourtois)
single_trial_data$gng_invalid_rt <- FALSE
single_trial_data$gng_invalid_rt[(single_trial_data$gng_resp != 45 & single_trial_data$gng_resp != 46 & single_trial_data$gng_resp != 47) & single_trial_data$gng_rt < 150 | single_trial_data$gng_rt > 500]  <- TRUE
count_outlier_GNG <- length(single_trial_data$gng_invalid_rt[single_trial_data$gng_miss_wrong_key_invalid_rt == TRUE])  # count number of outliers




####################   calculate  #####################################
####################    GNG RT    ##################################### 

# FA = False Alarm; FH = Fast Hit; CI = Correctly Inhibited; SH = Slow Hit
# RT not calculated for CI, because there is no rt
GNG_mean_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
GNG_mean_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
GNG_mean_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)




####################   calculate GNG     #####################################
####################  % response types   ##################################### 

GNG_count_FA   <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
GNG_count_FH   <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
GNG_count_CI   <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
GNG_count_SH   <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
GNG_count_miss <- length(single_trial_data[single_trial_data$gng_resp == 47,]$gng_rt)

GNG_all_resp_without_outliers_wrong_keys <- GNG_count_FA + GNG_count_FH + GNG_count_CI + GNG_count_SH + GNG_count_miss 

GNG_percent_FA   <- (GNG_count_FA/GNG_all_resp_without_outliers_wrong_keys)*100
GNG_percent_FH   <- (GNG_count_FH/GNG_all_resp_without_outliers_wrong_keys)*100
GNG_percent_CI   <- (GNG_count_CI/GNG_all_resp_without_outliers_wrong_keys)*100
GNG_percent_SH   <- (GNG_count_SH/GNG_all_resp_without_outliers_wrong_keys)*100
GNG_percent_miss <- (GNG_count_miss/GNG_all_resp_without_outliers_wrong_keys)*100



####################       define      #####################################
####################  word conditions  #####################################

# FA = False Alarm; FH = Fast Hit; CI = Correctly Inhibited; SH = Slow Hit; subsetting only required for outlier detection
neg_after_FA <- subset(single_trial_data, valence == "neg" & response_type == "FA" & word_resp == 51)
pos_after_FA <- subset(single_trial_data, valence == "pos" & response_type == "FA" & word_resp == 52)
neg_after_FH <- subset(single_trial_data, valence == "neg" & response_type == "FH" & word_resp == 51)
pos_after_FH <- subset(single_trial_data, valence == "pos" & response_type == "FH" & word_resp == 52)
neg_after_CI <- subset(single_trial_data, valence == "neg" & response_type == "CI" & word_resp == 51)
pos_after_CI <- subset(single_trial_data, valence == "pos" & response_type == "CI" & word_resp == 52)
neg_after_SH <- subset(single_trial_data, valence == "neg" & response_type == "SH" & word_resp == 51)
pos_after_SH <- subset(single_trial_data, valence == "pos" & response_type == "SH" & word_resp == 52)




####################     define       #####################################
####################  word outliers   ##################################### 

# assign TRUE to rt_values deviating more than 3 median absolute deviations (MAD; better use 2.5?)
single_trial_data <- mutate(single_trial_data,
                   outlier_words = ifelse((valence == "neg" & word_resp == 51),
                                          ifelse(response_type == "FA", 
                                                 (abs(word_rt - median(neg_after_FA$word_rt))/mad(neg_after_FA$word_rt))>3,
                                                 ifelse(response_type == "FH", 
                                                        (abs(word_rt - median(neg_after_FH$word_rt))/mad(neg_after_FH$word_rt))>3,
                                                        ifelse(response_type == "CI", 
                                                               (abs(word_rt - median(neg_after_CI$word_rt))/mad(neg_after_CI$word_rt))>3,
                                                               (abs(word_rt - median(neg_after_SH$word_rt))/mad(neg_after_SH$word_rt))>3)
                                                 )),
                                          # handle pos rows
                                          ifelse((valence == "pos" & word_resp == 52),
                                                 ifelse(response_type == "FA", 
                                                        (abs(word_rt - median(pos_after_FA$word_rt))/mad(pos_after_FA$word_rt))>3,
                                                        ifelse(response_type == "FH", 
                                                               (abs(word_rt - median(pos_after_FH$word_rt))/mad(pos_after_FH$word_rt))>3,
                                                               ifelse(response_type == "CI", 
                                                                      (abs(word_rt - median(pos_after_CI$word_rt))/mad(pos_after_CI$word_rt))>3,
                                                                      (abs(word_rt - median(pos_after_SH$word_rt))/mad(pos_after_SH$word_rt))>3))
                                                 ), FALSE))
)

count_outlier_words_FA_FH <- length(single_trial_data$outlier_words[single_trial_data$outlier_words == TRUE & (single_trial_data$response_type == "FA" | single_trial_data$response_type == "FH" )])  # count number of outliers




####################    calculat priming    ##################################
####################    effect with mean    ##################################

mean_neg_after_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE,]$word_rt)
mean_pos_after_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE,]$word_rt)
mean_neg_after_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE,]$word_rt)
mean_pos_after_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE,]$word_rt)
mean_neg_after_CI <- mean(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE,]$word_rt)
mean_pos_after_CI <- mean(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE,]$word_rt)
mean_neg_after_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE,]$word_rt)
mean_pos_after_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE,]$word_rt)

mean_priming_after_FA <- mean_pos_after_FA - mean_neg_after_FA
mean_priming_after_FH <- mean_neg_after_FH - mean_pos_after_FH
mean_priming_overall  <- (mean_pos_after_FA + mean_neg_after_FH) - (mean_neg_after_FA + mean_pos_after_FH) 




##################   calculate priming   ##################################
##################  effect with median   ##################################

# calculate median without exclusion of word outliers
median_neg_after_FA <- median(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51,]$word_rt)
median_pos_after_FA <- median(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52,]$word_rt)
median_neg_after_FH <- median(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51,]$word_rt)
median_pos_after_FH <- median(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52,]$word_rt)
median_neg_after_CI <- median(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51,]$word_rt)
median_pos_after_CI <- median(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52,]$word_rt)
median_neg_after_SH <- median(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51,]$word_rt)
median_pos_after_SH <- median(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52,]$word_rt)


median_priming_after_FA <- median_pos_after_FA - median_neg_after_FA
median_priming_after_FH <- median_neg_after_FH - median_pos_after_FH
median_priming_overall  <- (median_pos_after_FA + median_neg_after_FH) - (median_neg_after_FA + median_pos_after_FH) 



####################    count number    ##################################
####################     of events      ##################################

# count number of events,  outliers are removed; if I want to work with median, do not exclude outlier
count_neg_after_FA        <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_pos_after_FA        <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_neg_after_FH        <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_pos_after_FH        <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_neg_after_CI        <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_pos_after_CI        <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_neg_after_SH        <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_pos_after_SH        <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE,]$word_rt)

count_incorr_neg_after_FA <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_incorr_pos_after_FA <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_incorr_neg_after_FH <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_incorr_pos_after_FH <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_incorr_neg_after_CI <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_incorr_pos_after_CI <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_incorr_neg_after_SH <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$outlier_words == FALSE,]$word_rt)
count_incorr_pos_after_SH <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$outlier_words == FALSE,]$word_rt)

count_miss_neg_after_FA   <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 55,]$word_rt)
count_miss_pos_after_FA   <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 56,]$word_rt)
count_miss_neg_after_FH   <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 55,]$word_rt)
count_miss_pos_after_FH   <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 56,]$word_rt)
count_miss_neg_after_CI   <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 55,]$word_rt)
count_miss_pos_after_CI   <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 56,]$word_rt)
count_miss_neg_after_SH   <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 55,]$word_rt)
count_miss_pos_after_SH   <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 56,]$word_rt)


# misses are contained in overall number of events so that errors, correct and misses sum up to 100 %
count_all_neg_after_FA    <- count_neg_after_FA + count_incorr_neg_after_FA + count_miss_neg_after_FA
count_all_pos_after_FA    <- count_pos_after_FA + count_incorr_pos_after_FA + count_miss_pos_after_FA
count_all_neg_after_FH    <- count_neg_after_FH + count_incorr_neg_after_FH + count_miss_neg_after_FH
count_all_pos_after_FH    <- count_pos_after_FH + count_incorr_pos_after_FH + count_miss_pos_after_FH
count_all_neg_after_CI    <- count_neg_after_CI + count_incorr_neg_after_CI + count_miss_neg_after_CI
count_all_pos_after_CI    <- count_pos_after_CI + count_incorr_pos_after_CI + count_miss_pos_after_CI
count_all_neg_after_SH    <- count_neg_after_SH + count_incorr_neg_after_SH + count_miss_neg_after_SH
count_all_pos_after_SH    <- count_pos_after_SH + count_incorr_pos_after_SH + count_miss_pos_after_SH




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




####################  log and z transform   ##################################
####################    scr and word rt     ##################################

single_trial_data$scr_amplitude_log <- NA
single_trial_data$word_rt_z_score <- NA
single_trial_data$word_rt_log <- NA
single_trial_data[single_trial_data$gng_invalid_rt == FALSE,]$scr_amplitude_log <- log(single_trial_data[single_trial_data$gng_invalid_rt == FALSE,]$scr_amplitude + 1)                    # log transform SCR (for mean calculation excluding outlier in GNG only, exclude word outlier in single trial part, log-transform and standardize (z-score) SCR data)
single_trial_data$scr_amplitude_z_score_log <- scale(single_trial_data$scr_amplitude_log, center = TRUE, scale = TRUE)                                                                     # z-transform log-SCR
single_trial_data[single_trial_data$gng_invalid_rt == FALSE,]$word_rt_z_score <- scale(single_trial_data[single_trial_data$gng_invalid_rt == FALSE,]$word_rt, center = TRUE, scale = TRUE) # z-transform word rt
single_trial_data[single_trial_data$gng_invalid_rt == FALSE,]$word_rt_log <- log(single_trial_data[single_trial_data$gng_invalid_rt == FALSE,]$word_rt)                                    # log transform word rt


# mark trials followed or preceded by false alarm or wrong key in GNG or by incorrect word categorization or wrong key in word categorization
single_trial_data$followed_or_preceded_by_FA_or_wrong_key <- FALSE

for (i in 2:(nrow(single_trial_data)-1)) {
  current_row <- single_trial_data[i,]
  next_row <- single_trial_data[(i+1),]
  previous_row <- single_trial_data[(i-1),]
  if ((next_row$gng_resp %in% c(43,44,48,49) | next_row$word_resp %in% c(53,54,57,58) | previous_row$gng_resp %in% c(43,44,48,49) | previous_row$word_resp %in% c(53,54,57,58))) {            
    single_trial_data[i,]$followed_or_preceded_by_FA_or_wrong_key <- TRUE
  }     
}



# calculate mean SCR in conditions (some subjects contain NAs, because SCR saving was disrupted -> exclude NA for mean calculation)
# Remember: for CI (and for missing responses), SCR values in table are stimulus-locked to the GNG target; for FA, FH, SH the SCR values are response-locked -> see script "Preprocess_SCR_Export_Files"! there I assigned stimulus-locked SCR values to correctly inhibited trials and trials with missing response (makes no sense to take response-locked SCR value in these, because there is no response and the evaluation trigger is send late)
mean_SCR_after_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & (single_trial_data$word_resp == 51 | single_trial_data$word_resp == 52) & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$scr_amplitude_z_score_log,na.rm=TRUE)
mean_SCR_after_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & (single_trial_data$word_resp == 51 | single_trial_data$word_resp == 52) & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$scr_amplitude_z_score_log,na.rm=TRUE)
mean_SCR_after_CI <- mean(single_trial_data[single_trial_data$response_type == "CI" & (single_trial_data$word_resp == 51 | single_trial_data$word_resp == 52) & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$scr_amplitude_z_score_log,na.rm=TRUE)
mean_SCR_after_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & (single_trial_data$word_resp == 51 | single_trial_data$word_resp == 52) & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$scr_amplitude_z_score_log,na.rm=TRUE)



####################   Add df of one subject to   #####################################
####################  df containing all subjects  #####################################

df4save <- rbind(df4save, data.frame(subject,mean_neg_after_FA,mean_pos_after_FA,mean_neg_after_FH,mean_pos_after_FH,mean_neg_after_CI,mean_pos_after_CI,mean_neg_after_SH,mean_pos_after_SH,mean_priming_overall,mean_priming_after_FA,mean_priming_after_FH,median_neg_after_FA,median_pos_after_FA,median_neg_after_FH,median_pos_after_FH,median_neg_after_CI,median_pos_after_CI,median_neg_after_SH,median_pos_after_SH,median_priming_overall,median_priming_after_FA,median_priming_after_FH,count_neg_after_FA,count_pos_after_FA,count_neg_after_FH,count_pos_after_FH,count_outlier_words_FA_FH,percent_correct_neg_after_FA,percent_correct_pos_after_FA,percent_correct_neg_after_FH,percent_correct_pos_after_FH,percent_correct_neg_after_CI,percent_correct_pos_after_CI,percent_correct_neg_after_SH,percent_correct_pos_after_SH,GNG_mean_FA,GNG_mean_FH,GNG_mean_SH,GNG_percent_FA,GNG_percent_FH,GNG_percent_CI,GNG_percent_SH,GNG_percent_miss, mean_SCR_after_FA,mean_SCR_after_FH,mean_SCR_after_CI,mean_SCR_after_SH))
data4mixedmodels <- rbind(data4mixedmodels,single_trial_data)


}

####################   Prepare data4mixedmodels   #####################################
####################           for LMMs           #####################################
data4mixedmodels$condition     <- paste(data4mixedmodels$valence,"_after_",data4mixedmodels$response_type)
data4mixedmodels$response_type <- factor(data4mixedmodels$response_type, levels=c("FA","FH","CI","SH"))         # reorder factor levels! (this is important, because first level is used as reference factor in LMM)
data4mixedmodels$valence       <- factor(data4mixedmodels$valence)                                              # automatic factor order is alphabetical (1 = neg, 2 = pos)
data4mixedmodels$subjectID     <- factor(data4mixedmodels$subjectID)




 ####################     save df4save    ##################################
 ####################    as excel file    ##################################
 
 
 setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses")    # setting a different folder as working directory to prevent saving stuff into the folder containing the logfiles
 
 date_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
 
 filename <- paste("SummaryStatisticsSingleTrial_For_",length(logfiles),"_subjects_",date_time, ".xlsx", sep = "")
 
 # write.xlsx(df4save, filename)
 

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
 
 # plot rt all conditions
 df4plotRT <- data.frame(
   response = c("false alarm","false alarm","fast hit","fast hit","correctly inhibited","correctly inhibited","slow hit", "slow hit"),
   valence = c("neg","pos","neg","pos","neg","pos","neg","pos"),
   conditions = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH","neg_after_CI","pos_after_CI","neg_after_SH","pos_after_SH"),
   mean = as.numeric(c(descriptive_statistics[2,2],descriptive_statistics[2,3],descriptive_statistics[2,4],descriptive_statistics[2,5],descriptive_statistics[2,6],descriptive_statistics[2,7],descriptive_statistics[2,8],descriptive_statistics[2,9])),
   se = as.numeric(c(descriptive_statistics[3,2],descriptive_statistics[3,3],descriptive_statistics[3,4],descriptive_statistics[3,5],descriptive_statistics[3,6],descriptive_statistics[3,7],descriptive_statistics[3,8],descriptive_statistics[3,9])))
 
 df4plotRT$response <- factor(df4plotRT$response, levels = c("false alarm","fast hit","slow hit","correctly inhibited"))
 
 ggplot(df4plotRT,x = response,y = mean, aes(response, mean, fill = valence))+
   geom_bar(stat="identity", position=position_dodge()) +                                                       # add bars, based on stats values; dodge to avoid stacked bars
   geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.95), width=0.1) +    # add error bars
   ggtitle("Response Time Word Categorization") +                                                               # add title
   xlab("Previous Response") + ylab("Reaction Time (ms)") +                                                     # label axes
   guides(fill=guide_legend(title="Word Valence")) +                                                            # change legend title
   theme(axis.title = element_text(size = 18, hjust = 0.5)) +                                                   # center title
   theme(axis.text=element_text(size=15)) +
   theme(legend.title=element_text(size=15)) +
   theme(legend.text=element_text(size=15)) +
   theme(plot.title = element_text(size = 20, hjust = 0.5)) +
   theme(legend.position="bottom") +
   coord_cartesian(ylim = c(300,900)) +
   scale_fill_manual(values=c("mediumblue", "limegreen"))                                                       # change bar colors
 
 
 # plot accuracy all conditions
 df4plotACC <- data.frame(
   response = c("false alarm","false alarm","fast hit","fast hit","correctly inhibited","correctly inhibited","slow hit", "slow hit"),
   valence = c("neg","pos","neg","pos","neg","pos","neg","pos"),
   conditions = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH","neg_after_CI","pos_after_CI","neg_after_SH","pos_after_SH"),
   mean = as.numeric(c(descriptive_statistics[2,29],descriptive_statistics[2,30],descriptive_statistics[2,31],descriptive_statistics[2,32],descriptive_statistics[2,33],descriptive_statistics[2,34],descriptive_statistics[2,35],descriptive_statistics[2,36])),
   se = as.numeric(c(descriptive_statistics[3,29],descriptive_statistics[3,30],descriptive_statistics[3,31],descriptive_statistics[3,32],descriptive_statistics[3,33],descriptive_statistics[3,34],descriptive_statistics[3,35],descriptive_statistics[3,36])))
 
 df4plotACC$response <- factor(df4plotACC$response, levels = c("false alarm","fast hit","slow hit","correctly inhibited"))
 
 ggplot(df4plotACC,x = response,y = mean, aes(response, mean, fill = valence))+
   geom_bar(stat="identity", position=position_dodge()) +                                                       # add bars, based on stats values; dodge to avoid stacked bars
   geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.95), width=0.1) +    # add error bars
   ggtitle("Accuracy Word Categorization") +                                                                    # add title
   xlab("Previous Response") + ylab("Accuracy (%)") +                                                           # label axes
   guides(fill=guide_legend(title="Word Valence")) +                                                            # change legend title
   theme(axis.title = element_text(size = 18, hjust = 0.5)) +                                                   # center title
   theme(axis.text=element_text(size=15)) +
   theme(legend.title=element_text(size=15)) +
   theme(legend.text=element_text(size=15)) +
   theme(plot.title = element_text(size = 20, hjust = 0.5)) +
   theme(legend.position="bottom") +
   coord_cartesian(ylim = c(30,120)) +
   scale_fill_manual(values=c("mediumblue", "limegreen"))                                                       # change bar colors
 
 
 # plot SCR all conditions
 df4plotSCR <- data.frame(
   response = c("false alarm","fast hit","correctly inhibited","slow hit"),
   valence = c("neg","pos","neg","pos"),
   conditions = c("mean_SCR_after_FA","mean_SCR_after_FH","mean_SCR_after_CI","mean_SCR_after_SH"),
   mean = as.numeric(c(descriptive_statistics[2,45],descriptive_statistics[2,46],descriptive_statistics[2,47],descriptive_statistics[2,48])),
   se = as.numeric(c(descriptive_statistics[3,45],descriptive_statistics[3,46],descriptive_statistics[3,47],descriptive_statistics[3,48])))
 
 df4plotSCR$response <- factor(df4plotSCR$response, levels = c("false alarm","fast hit","slow hit","correctly inhibited"))
 
 ggplot(df4plotSCR,x = response,y = mean, aes(response, mean))+
   geom_bar(stat="identity", position=position_dodge(), fill = "mediumblue") +                                  # add bars, based on stats values; dodge to avoid stacked bars; change bar color
   geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.95), width=0.1) +    # add error bars
   ggtitle("Skin Conductance Response") +                                                                       # add title
   xlab("Response Type") + ylab("Logarithmized SCR (µS)") +                                                     # label axes
   theme(axis.title = element_text(size = 18, hjust = 0.5)) +                                                   # center title
   theme(axis.text=element_text(size=15)) +
   theme(legend.title=element_text(size=18)) +
   theme(legend.text=element_text(size=15)) +
   theme(plot.title = element_text(size = 20, hjust = 0.5)) +
   coord_cartesian(ylim = c(-0.2,0.7)) +
   scale_y_continuous(breaks=seq(0,40,0.2))                                                                       # set specific tic marks                                       
 
 
 
 
 
 
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
                       times = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH","neg_after_CI","pos_after_CI","neg_after_SH","pos_after_SH"))
 
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

 
 # post hoc comparison: when calculating within-subject ANOVAs, no direct post hoc test is possible; I can only run paitwise t-tests with corresponding adjustment (e.g. Holm, which is better than Bonferroni)
 anova_rt_posthoc_response_type       <- pairwise.t.test(df4anova$rt,df4anova$response_type,p.adjust.method="holm",paired=TRUE)
 anova_rt_posthoc_word_valence        <- pairwise.t.test(df4anova$rt,df4anova$word_valence,p.adjust.method="holm",paired=TRUE) 
 anova_rt_posthoc_priming             <- pairwise.t.test(df4anova$rt,df4anova$condition,p.adjust.method="holm",paired=TRUE) 
 
 anova_accuracy_posthoc_response_type <- pairwise.t.test(df4anova$accuracy,df4anova$response_type,p.adjust.method="holm",paired=TRUE)
 anova_accuracy_posthoc_word_valence  <- pairwise.t.test(df4anova$accuracy,df4anova$word_valence,p.adjust.method="holm",paired=TRUE) 
 anova_accuracy_posthoc_priming       <- pairwise.t.test(df4anova$accuracy,df4anova$condition,p.adjust.method="holm",paired=TRUE) 
 
 

 
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
 anova_GNG_rt_posthoc_response_type <- pairwise.t.test(df4anova_GNG_rt$rt,df4anova_GNG_rt$condition,p.adjust.method="holm",paired=TRUE)
 
 
 
 #####################    ANOVA SCR and   ####################################
 #####################     GNG percent    #################################### 
 
 
 # ANOVA requires several rows for each subject, each one row per factor: here I reorder data by reshaping existing data frame
 df4anova_GNG_scr <-   reshape(data = df4save[,c(1,40:43,45:48)], 
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
 anova_GNG_percent_posthoc_response_type <- pairwise.t.test(df4anova_GNG_scr$percent,df4anova_GNG_scr$condition,p.adjust.method="holm",paired=TRUE)
 
 
 
 # calculate ANOVA with ezANOVA for SCR
 anova_GNG_SCR <- ezANOVA(data = df4anova_GNG_scr, 
                          dv = .(SCR), 
                          wid = .(subject), 
                          within = .(condition), 
                          detailed = TRUE)
 
 
 
 # post hoc comparison: when calculating within-subject ANOVAs, no direct post hoc test is possible; I can only run paitwise t-tests with corresponding adjustment (e.g. Holm, which is better than Bonferroni)
 anova_GNG_SCR_posthoc_response_type <- pairwise.t.test(df4anova_GNG_scr$SCR,df4anova_GNG_scr$condition,p.adjust.method="holm",paired=TRUE)
 

 
 
 ########################################################################################
 #######################   Correlations with    #########################################
 #######################        Traits          #########################################
 ########################################################################################
 
 # Read in questionnaire data, merge with aggregated data, create correlation matrix
 setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study")
 questionnaires <- list.files(pattern = ".sav")                                                         # make sure that only one .sav file (= the current PEQ export) exists there!
 questionnaires <- read.spss(questionnaires, to.data.frame = TRUE)
 questionnaires <- questionnaires[,c("CODE", "BD2SUMT0","BASO00T0","BASO01T0","BASO02T0","BASO03T0","BASO04T0","FMPO00T0","FMPO01T0","FMPO02T0","FMPO03T0","FMPO04T0","NNGO00T0","NNGO01T0","OCISUMT0", "OCIO00T0", "OCIO01T0","OCIO02T0","OCIO03T0","OCIO04T0","OCIO05T0","PANO00T0","PANO01T0","PSWSUMT0","STSSUMT0","STTSUMT0","TCISUMT0","TCIO00T0","TCIO01T0","TCIO02T0","TCIO03T0","WSTSUMT0")] # only keep relevant columns (I chose to keep sum scores instead of mean scores, because for BDI and OCi only sum scores are available and missing data are not possible because data aquired with tablet)
 questionnaires <- questionnaires[grep("behav", questionnaires$CODE), ]                                 # only keep participants from behavioral study
 questionnaires <- questionnaires[, colSums(is.na(questionnaires)) != nrow(questionnaires)]             # remove columns with only NA
 questionnaires$CODE <- paste0(questionnaires$CODE,".txt")                                              # convert code for joining (make it identical to code in df4save)                               
 colnames(questionnaires) <-c("CODE", "BDI_Total", "BIS_Total", "BAS_Total", "BAS_Drive", "BAS_Fun_Seeking", "BAS_Reward_Responsiveness", "FMPS_CMD", "FMPS_PEC", "FMPS_PST", "FMPS_ORG", "FMPS_PER/Total", "NEO_Neuroticism", "NEO_Conscientiousness", "OCI_Total", "OCI_Washing", "OCI_Checking", "OCI_Ordering", "OCI_Obsessions", "OCI_Hoarding", "OCI_Neutralising", "PANAS_Pos", "PANAS_Neg", "PSWQ_Total", "STAI_State", "STAI_Trait", "TCI_Total", "TCI_Anticipatory_worry", "TCI_Fear_of_uncertainty", "TCI_Shyness", "TCI_Fatigability", "WST_Total")   # rename columns 
 
 
 df4correlation_matrix <- left_join(df4save,questionnaires, by = c("subject" = "CODE"))
 # possibly exclude subject 14 (great outlier, that causes most correlations/trends)!!!!
 
 
 # test for normality
 normality <- do.call(rbind, lapply(df4correlation_matrix[,-1], function(x) shapiro.test(x)[c("statistic", "p.value")]))
 
 
 # create correlation matrices using my function "function_correlation_matrix_with_significance_levels.R"
 source("P:/LuisaBalzus/1_PhD_Project/5_ModERN_Vorstudie/10_Single_Trial_Analysis_EEG_SCR_GNG/function_correlation_matrix_with_significance_levels.R")
 
 
 df4Pearson_Correlation <- df4correlation_matrix[,c("mean_priming_overall", "mean_priming_after_FA", "mean_priming_after_FH","mean_SCR_after_FA", "BIS_Total","BAS_Fun_Seeking", "BAS_Reward_Responsiveness", "FMPS_CMD", "FMPS_PST", "FMPS_PER/Total", "NEO_Neuroticism", "PANAS_Pos", "PSWQ_Total", "STAI_Trait", "TCI_Total", "TCI_Fear_of_uncertainty", "TCI_Shyness", "TCI_Fatigability")]   # only keep relevant columns 
 method <- "pearson"
 df4correlation <- df4Pearson_Correlation
 correlations_pearson <- correlations_with_significance_level(df4Pearson_Correlation)
 
 
 df4Kendall_Correlation <- df4correlation_matrix[,c("mean_priming_overall", "mean_priming_after_FA", "mean_priming_after_FH","mean_SCR_after_FA","BDI_Total", "BAS_Total", "BAS_Drive", "FMPS_PEC","FMPS_ORG","NEO_Conscientiousness", "OCI_Total", "OCI_Washing", "OCI_Checking", "OCI_Ordering", "OCI_Obsessions", "OCI_Hoarding", "OCI_Neutralising","PANAS_Neg","STAI_State","TCI_Anticipatory_worry")]   # only keep relevant columns 
 method <- "kendall"
 df4correlation <- df4Kendall_Correlation
 correlations_kendall <- correlations_with_significance_level(df4Kendall_Correlation)
 
 

 # Template for plotting significant correlations
 ggplot(df4correlation_matrix, aes (x= TCI_Fear_of_uncertainty, y = mean_SCR_after_FA)) +
   geom_point(color="mediumblue", size=3) +
   geom_smooth(method=lm, se=FALSE, color="limegreen", size = 1) +
   ggtitle("TCI_Fear_of_uncertainty and mean_SCR_after_FA") +
   labs(x= "TCI_Fear_of_uncertainty", y = "mean_SCR_after_FA") +
   theme(axis.title = element_text(size = 18, hjust = 0.5)) +
   theme(axis.text=element_text(size=15)) +
   theme(legend.title=element_text(size=18)) +
   theme(legend.text=element_text(size=15)) +
   theme(plot.title = element_text(size = 20, hjust = 0.5))
 

 
 
 
 #####################    linear mixed models    ####################################
 #####################      words (mean RT)      #################################### 
 
 # prepare master_words for single trial analysis: exclude incorrect word responses, misses, wrong keys; create column condition; turn relevant variables in factors
 data4mixedmodels_words               <- subset(data4mixedmodels,word_resp <= 52 & gng_resp <= 46 & outlier_words == FALSE)

 # set contrasts
 contrasts(data4mixedmodels_words$response_type) <- contr.sdif(4)
 contrasts(data4mixedmodels_words$valence) <- contr.sdif(2)



 # calculate LMM for rt and plot
 emm_options(pbkrtest.limit = 14000)
 emm_options(lmerTest.limit = 14000) # due to warning D.f. calculations have been disabled because the number of observations exceeds 3000.To enable adjustments, set emm_options(pbkrtest.limit = 5866)

 LMM_rt <- lmer(word_rt ~ response_type * valence + (1|subjectID), data=data4mixedmodels_words, REML = FALSE)
 summary(LMM_rt)
 anova(LMM_rt)                                                                                     # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 results_LMM_rt <- analyze(LMM_rt)                                                                 # print results
 print(results_LMM_rt)
 results_LMM_rt <- get_contrasts(LMM_rt, "response_type * valence")                                # Provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
 print(results_LMM_rt$contrasts)                                                                   # print contrasts
 print(results_LMM_rt$means)                                                                       # investigate means

 ggplot(results_LMM_rt$means, aes(x=response_type, y=Mean, color=valence, group=valence)) +        # plot
   geom_line(position = position_dodge(.3)) +
   geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher),
                   position = position_dodge(.3)) +
   ylab("Response Time in ms") +
   xlab("Response Type") +
   theme_bw()


 # check normality of residuals
 qqnorm(resid(LMM_rt))
 qqline(resid(LMM_rt))
 shapiro.test(resid(LMM_rt))



 # post hoc tests for resolving main effects and interaction effects
 summary(glht(LMM_rt, linfct=mcp(response_type = "Tukey")), test = adjusted("holm"))
 summary(glht(LMM_rt, linfct=mcp(valence = "Tukey")), test = adjusted("holm"))                   # same result can be obtained by lmer_rt_posthoc_response_type <- emmeans(model_lmer_rt, ~ word_valence);;pairs(lmer_rt_posthoc_response_type)
 LMM_rt_posthoc_priming <- emmeans(LMM_rt, ~ valence|response_type, adjust="tukey")              # compare word valence within each level of response_type
 pairs(LMM_rt_posthoc_priming)

 


 #####################    linear mixed models    ####################################
 #####################          GNG (RT)         ####################################

 # prepare master_GNG for single trial analysis: ecxlude CI, misses, wrong keys
 data4mixedmodels_GNG <- subset(data4mixedmodels,gng_resp <= 44 & gng_invalid_rt == FALSE)

 # set contrasts
 contrasts(data4mixedmodels_GNG$response_type) <- contr.sdif(3)  # 3 levels: FA, FH, SH



 # calculate LMM for rt and plot
 LMM_GNG_rt <- lmer(gng_rt ~ response_type + (1|subjectID), data=data4mixedmodels_GNG, REML = FALSE)
 summary(LMM_GNG_rt)
 anova(LMM_GNG_rt)                                                                            # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 results_LMM_GNG_rt <- analyze(LMM_GNG_rt)                                                    # print results
 print(results_LMM_GNG_rt)
 results_LMM_GNG_rt <- get_contrasts(LMM_GNG_rt, "response_type")                             # Provide the model and the factors to contrast;; by default, get_contrasts uses the Tukey method for p value adjustment; add ,adjust="none" to turn off automatic p value correction after Tucky
 print(results_LMM_GNG_rt$contrasts)                                                          # print contrasts
 print(results_LMM_GNG_rt$means)                                                              # investigate means

 ggplot(results_LMM_GNG_rt$means, aes(x=response_type, y=Mean, color=response_type)) +        # plot
   geom_line(position = position_dodge(.3)) +
   geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher),
                   position = position_dodge(.3)) +
   ylab("Response Time in ms") +
   xlab("Response Type") +
   theme_bw()


 # post hoc tests for resolving main effect
 summary(glht(LMM_GNG_rt, linfct=mcp(response_type = "Tukey")), test = adjusted("holm"))



 
 
 
 
 ########################################################################
 ##############################SCR#######################################
 ########################################################################
 ########################################################################

 # exclude CI, exclude trials with outlier or error in word or miss/wrong key in GNG, exclude trials with invalid GNG rt, exclude trials followed or preceded by false alarm or wrong key in GNG or by incorrect word categorization or wrong key in word categorization
 data4mixedmodels_scr <- subset(data4mixedmodels,word_resp <= 52 & gng_resp <= 44 & outlier_words == FALSE & gng_invalid_rt == FALSE & followed_or_preceded_by_FA_or_wrong_key == FALSE)


 # set contrasts for fixed effects
 contrasts(data4mixedmodels_scr$response_type) <- contr.sdif(3)
 contrasts(data4mixedmodels_scr$valence) <- contr.sdif(2)

 
 # calculate LMM for SCR and plot (does not converge)
 emm_options(pbkrtest.limit = 7000)
 emm_options(lmerTest.limit = 7000) # due to warning D.f. calculations have been disabled because the number of observations exceeds 3000.To enable adjustments, set emm_options(pbkrtest.limit = 5866)
 
 LMM_SCR <- lmer(scr_amplitude_log ~ response_type * valence * word_rt_log + (1+ response_type * valence * word_rt_log|subjectID)  + (word_rt_log|word), data=data4mixedmodels_scr, REML = FALSE)
 summary(LMM_SCR)
 anova(LMM_SCR)                                                                          # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 results_LMM_SCR <- analyze(LMM_SCR)                                                     # print results
 print(results_LMM_SCR)
 results_LMM_SCR <- get_contrasts(LMM_SCR, "response_type * valence")                    # Provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
 print(results_LMM_SCR$contrasts)                                                        # print contrasts
 print(results_LMM_SCR$means)                                                            # investigate means

 
 # test reduced model: n.s.
 LMM_SCR <- lmer(scr_amplitude_log ~ response_type * valence * word_rt_log + (1|subjectID), data=data4mixedmodels_scr, REML = FALSE)
 summary(LMM_SCR)
 anova(LMM_SCR)                                                                          # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 results_LMM_SCR <- analyze(LMM_SCR)                                                     # print results
 print(results_LMM_SCR)
 results_model_lmer_SCR <- get_contrasts(LMM_SCR, "response_type * valence")             # Provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
 print(results_LMM_SCR$contrasts)                                                        # print contrasts
 print(results_LMM_SCR$means)                                                            # investigate means
 

 
 # try if SCR to errors is predicted by w_rt_log -> interaction valence*word_rt n.s.
 data4mixedmodels_scr_FA <- subset(data4mixedmodels_scr, gng_resp == 43 | gng_resp == 44)
 contrasts(data4mixedmodels_scr_FA$valence) <- contr.sdif(2)

 LMM_SCR_FA <- lmer(scr_amplitude_log ~ valence * word_rt_log + (1|subjectID), data=data4mixedmodels_scr_FA, REML = FALSE)
 summary(LMM_SCR_FA)
 anova(LMM_SCR_FA)                                                                      # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
 results_LMM_SCR_FA <- analyze(LMM_SCR_FA)                                       # print results
 print(results_LMM_SCR_FA)
 results_LMM_SCR_FA <- get_contrasts(LMM_SCR_FA, "valence*word_rt_log")             # Provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
 print(results_LMM_SCR_FA$contrasts)                                                    # print contrasts
 print(results_LMM_SCR_FA$means)



 
 
 
 
 
