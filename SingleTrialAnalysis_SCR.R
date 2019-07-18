  ##### Proprocessing of Evaluative GNG Task in Behavioral Experiment on Affective Priming Effect 2018 ####
  ##### by Luisa Balzus
  ##### 09.11.2018
  
  
  library(dplyr)

  
  # clear environment
  rm(list=ls())
  
  # force R to not use exponential notation
  options(scipen = 999)
  
  # create empty data frames to write aggregated data and single trial data in it
  data4mixedmodels <- data.frame()   
  df4save <- data.frame()   
  
  
  
  ###################   Loop over Subjects to Preprocess Logfile Data and SCR Data for Statistical Analyses   ####################
  
  # Load logfiles and start loop
  logfiles <- list.files("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_Behavioral", pattern = ".txt")       
  
  for (subject in logfiles){                                                                               # loop reading txt-files as table, ommit first 58 lines and added lines after trials; use first line as header
    setwd("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_Behavioral")                 # path to folder containing the log files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 
    raw_log <- read.table(subject, skip = 58, fill = TRUE, header = TRUE, nrows = 516)
    name_of_subject <-  gsub("\\.txt*","",subject)
  
    
  # Load SCR files 
  setwd("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/9_SCR_Export_Preprocessed")               # path to folder containing the scr files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 
  filename_scr <- paste0(name_of_subject,"_SCR_Export_era.txt")
  raw_scr <- read.table(filename_scr, header = TRUE, sep=",")
  
  
  
  
  ####################   merge log and scr data   ####################
 
  # show progress
  message('Processing ', subject, ' of ', length(logfiles))
  
  
  # return all rows of log, merged with matching columns of scr; rows in log with no match in scr will have NAs in the new columns; if there are multiple matches between log and scr, all combinations of the matches are returned
  data_log_scr <- left_join(raw_log,raw_scr,by = c("number"= "Participant.ID","count" = "Trial.ID","t_type"="Event.ID.Tar","t_resp"="Event.ID.TarResp","w_type"="Event.ID.Word","w_resp"="Event.ID.WordResp"))
  
  
  # check whether there were rows in logfile with no match in scr (is true for subject 25 and 29, because SCR recording was interrupted)
  if (any(is.na(data_log_scr$DDA.AreaSum.TarResp))){
    print("There are rows in logfile with no match in SCR file! Should be only caused by missing SCR trials in subjects 25 (trial 112-137) and 29 (trial 185.5-200)!")
  }else{
      print("Matching of logfile and SCR file is ok. Each row in logfile has a matching row in SCR file")
    }
  
  
  # check whether there were there are multiple matches between logfile and scr file
  if (any(duplicated(data_log_scr[,c('number','count')]))){
    stop("There are rows in logfile with multiple matches in SCR file!")
  }else{
    print("Matching of logfile and SCR-file is ok. No row in log has several matching row in SCR")
  }
  
  
  
  
  ####################   clean merged data frame   #################### 
  
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
  
  
  
  
  ####################   define conditions   ####################
  
  single_trial_data <- mutate(single_trial_data, 
                     valence = ifelse(word_condition > 200, "pos", "neg"),                                                      
                     response_type = ifelse((gng_resp == 43 | gng_resp == 44), "FA", 
                                            ifelse(gng_resp == 41, "FH", 
                                                   ifelse(gng_resp == 45 | gng_resp == 46, "CI", 
                                                          ifelse(gng_resp == 42 , "SH",
                                                                 "Miss_or_False_Key"))))
  )
  
  
  
  ####################   define GNG  outliers   ####################   
  
  # assign TRUE to column GNG_outlier < 100 or > 700 ms; here, trials with GNG outliers are excluded ALSO for word categorization in corresponding trial! (Pourtois used 150 and 500 ms, but unclear whether trials with GNG outliers were excluded for word classification there, too; I do not want to loose too many words, so I use 100 and 700 ms) (priming after FA is 119 if no gng outliers are excluded, 116 if 100/700 ms, 104 if 150/500 ms; priming after FH remains always the same)
  single_trial_data$gng_invalid_rt <- FALSE
  single_trial_data$gng_invalid_rt[(single_trial_data$gng_resp != 45 & single_trial_data$gng_resp != 46 & single_trial_data$gng_resp != 47) & single_trial_data$gng_rt < 100 | single_trial_data$gng_rt > 700]  <- TRUE
  count_outlier_GNG <- length(single_trial_data$gng_invalid_rt[single_trial_data$gng_invalid_rt == TRUE])  # count number of outliers
  
  
  
  
  
  ####################   calculate GNG RT   #################### 
  
  # FA = False Alarm; FH = Fast Hit; CI = Correctly Inhibited; SH = Slow Hit
  # RT not calculated for CI, because there is no rt
  GNG_mean_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
  GNG_mean_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
  GNG_mean_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
  
  
  
  
  
  ####################   calculate GNG % response types   #################### 
  
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
  
  
  
  
  
  ####################   define word conditions  ####################
  
  # FA = False Alarm; FH = Fast Hit; CI = Correctly Inhibited; SH = Slow Hit; subsetting only required for outlier detection
  neg_after_FA <- subset(single_trial_data, valence == "neg" & response_type == "FA" & word_resp == 51 & gng_invalid_rt == FALSE)
  pos_after_FA <- subset(single_trial_data, valence == "pos" & response_type == "FA" & word_resp == 52 & gng_invalid_rt == FALSE)
  neg_after_FH <- subset(single_trial_data, valence == "neg" & response_type == "FH" & word_resp == 51 & gng_invalid_rt == FALSE)
  pos_after_FH <- subset(single_trial_data, valence == "pos" & response_type == "FH" & word_resp == 52 & gng_invalid_rt == FALSE)
  neg_after_CI <- subset(single_trial_data, valence == "neg" & response_type == "CI" & word_resp == 51 & gng_invalid_rt == FALSE)
  pos_after_CI <- subset(single_trial_data, valence == "pos" & response_type == "CI" & word_resp == 52 & gng_invalid_rt == FALSE)
  neg_after_SH <- subset(single_trial_data, valence == "neg" & response_type == "SH" & word_resp == 51 & gng_invalid_rt == FALSE)
  pos_after_SH <- subset(single_trial_data, valence == "pos" & response_type == "SH" & word_resp == 52 & gng_invalid_rt == FALSE)
  
  
  
  
  
  ####################   define word outliers   #################### 
  
  # assign TRUE to rt_values deviating more than 3 median absolute deviations (MAD; value of 3 is quite conservative and leads to exclusion of fewer events)
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
  
  
  
  
  
  ####################   calculat priming effect  word classification with mean    ####################
  
  mean_neg_after_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  mean_pos_after_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  mean_neg_after_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  mean_pos_after_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  mean_neg_after_CI <- mean(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  mean_pos_after_CI <- mean(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  mean_neg_after_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  mean_pos_after_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  
  mean_priming_after_FA <-  mean_pos_after_FA - mean_neg_after_FA
  mean_priming_after_FH <-  mean_neg_after_FH - mean_pos_after_FH
  mean_priming_overall  <- (mean_pos_after_FA + mean_neg_after_FH) - (mean_neg_after_FA + mean_pos_after_FH) 
 
 
  
 
  
  ####################   calculate priming effect word classification with median   ####################
  
  # calculate median without exclusion of word outliers and gng outliers
  median_neg_after_FA <- median(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51,]$word_rt)
  median_pos_after_FA <- median(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52,]$word_rt)
  median_neg_after_FH <- median(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51,]$word_rt)
  median_pos_after_FH <- median(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52,]$word_rt)
  median_neg_after_CI <- median(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51,]$word_rt)
  median_pos_after_CI <- median(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52,]$word_rt)
  median_neg_after_SH <- median(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51,]$word_rt)
  median_pos_after_SH <- median(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52,]$word_rt)
  
  
  median_priming_after_FA <-  median_pos_after_FA - median_neg_after_FA
  median_priming_after_FH <-  median_neg_after_FH - median_pos_after_FH
  median_priming_overall  <- (median_pos_after_FA + median_neg_after_FH) - (median_neg_after_FA + median_pos_after_FH) 
  
  
  
  
  
  ####################   count number of events word classification   ####################
  
  # count number of events,  outliers are removed; if I want to work with median, do not exclude outlier
  count_neg_after_FA        <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_pos_after_FA        <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_neg_after_FH        <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_pos_after_FH        <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_neg_after_CI        <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_pos_after_CI        <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_neg_after_SH        <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_pos_after_SH        <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  
  count_incorr_neg_after_FA <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_incorr_pos_after_FA <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_incorr_neg_after_FH <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_incorr_pos_after_FH <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_incorr_neg_after_CI <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_incorr_pos_after_CI <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_incorr_neg_after_SH <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  count_incorr_pos_after_SH <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  
  count_miss_neg_after_FA   <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 55 & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  count_miss_pos_after_FA   <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 56 & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  count_miss_neg_after_FH   <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 55 & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  count_miss_pos_after_FH   <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 56 & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  count_miss_neg_after_CI   <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 55 & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  count_miss_pos_after_CI   <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 56 & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  count_miss_neg_after_SH   <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 55 & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  count_miss_pos_after_SH   <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 56 & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  
  
 
  
  
  ####################   calculate accuracy word classification     ####################
  
  # misses are contained in overall number of events so that errors, correct and misses sum up to 100 %
  percent_correct_neg_after_FA <- count_neg_after_FA / (count_neg_after_FA + count_incorr_neg_after_FA + count_miss_neg_after_FA) * 100 
  percent_correct_pos_after_FA <- count_pos_after_FA / (count_pos_after_FA + count_incorr_pos_after_FA + count_miss_pos_after_FA) * 100  
  percent_correct_neg_after_FH <- count_neg_after_FH / (count_neg_after_FH + count_incorr_neg_after_FH + count_miss_neg_after_FH) * 100  
  percent_correct_pos_after_FH <- count_pos_after_FH / (count_pos_after_FH + count_incorr_pos_after_FH + count_miss_pos_after_FH) * 100  
  percent_correct_neg_after_CI <- count_neg_after_CI / (count_neg_after_CI + count_incorr_neg_after_CI + count_miss_neg_after_CI) * 100 
  percent_correct_pos_after_CI <- count_pos_after_CI / (count_pos_after_CI + count_incorr_pos_after_CI + count_miss_pos_after_CI) * 100  
  percent_correct_neg_after_SH <- count_neg_after_SH / (count_neg_after_SH + count_incorr_neg_after_SH + count_miss_neg_after_SH) * 100  
  percent_correct_pos_after_SH <- count_pos_after_SH / (count_pos_after_SH + count_incorr_pos_after_SH + count_miss_pos_after_SH) * 100  
  
  
  
  
  
  ####################   log and z transform scr and word rt     ####################
  
  single_trial_data$scr_amplitude_z_score_log <- NA    # z transformation depends on mean -> exclude outliers for that
  single_trial_data$word_rt_z_score <- NA
  
  single_trial_data$scr_amplitude_log <- log(single_trial_data$scr_amplitude + 1)                    # log transform SCR (for mean calculation excluding outlier in GNG only, exclude word outlier in single trial part, log-transform and standardize (z-score) SCR data)
  single_trial_data[single_trial_data$gng_invalid_rt == FALSE,]$scr_amplitude_z_score_log <- scale(single_trial_data[single_trial_data$gng_invalid_rt == FALSE,]$scr_amplitude_log, center = TRUE, scale = TRUE)                                                                     # z-transform log-SCR
  single_trial_data[(single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE),]$word_rt_z_score <- scale(single_trial_data[(single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE),]$word_rt, center = TRUE, scale = TRUE) # z-transform word rt
  single_trial_data$word_rt_log <- log(single_trial_data$word_rt)                                    # log transform word rt
  
  
  
  
  
  ####################   mark trials followed or preceded by false alarm or wrong key in GNG or by incorrect word categorization or wrong key in word categorization   ####################

  single_trial_data$followed_or_preceded_by_FA_or_wrong_key <- FALSE
  
  for (i in 2:(nrow(single_trial_data)-1)) {
    current_row <- single_trial_data[i,]
    next_row <- single_trial_data[(i+1),]
    previous_row <- single_trial_data[(i-1),]
    if ((next_row$gng_resp %in% c(43,44,48,49) | next_row$word_resp %in% c(53,54,57,58) | previous_row$gng_resp %in% c(43,44,48,49) | previous_row$word_resp %in% c(53,54,57,58))) {            
      single_trial_data[i,]$followed_or_preceded_by_FA_or_wrong_key <- TRUE
    }     
  }
  
  
  
  ####################   calculate priming facilitation score for LMM   #####################
  
  single_trial_data$facilitation_score <- NA      # score shall be higher the more priming a person shows in a trial (i.e. fast for neg words after FA, slow for pos words; fast for pos words after FH, slow for neg words)
  
  for (j in 1:nrow(single_trial_data)) {
    if (single_trial_data[j,]$response_type == "FH" & single_trial_data[j,]$valence == "neg")  {
    single_trial_data[j,]$facilitation_score <- single_trial_data[j,]$word_rt - mean_neg_after_FH
    }
    if (single_trial_data[j,]$response_type == "FH" & single_trial_data[j,]$valence == "pos")  {
      single_trial_data[j,]$facilitation_score <- mean_pos_after_FH - single_trial_data[j,]$word_rt
    }
    if (single_trial_data[j,]$response_type == "FA" & single_trial_data[j,]$valence == "neg")  {
      single_trial_data[j,]$facilitation_score <- mean_neg_after_FA - single_trial_data[j,]$word_rt
    }
    if (single_trial_data[j,]$response_type == "FA" & single_trial_data[j,]$valence == "pos")  {
      single_trial_data[j,]$facilitation_score <- single_trial_data[j,]$word_rt - mean_pos_after_FA
    }
  }
  # I don't z-standardize facilitation score, because subjects are entered as random intercept in LMM anyway; I don't log-transform facilitation score, because it is normally distributed -> hist(single_trial_data$facilitation_score)
  
  
  
  
  
  
  ####################   calculate mean SCR in conditions (some subjects contain NAs, because SCR recording was disrupted -> exclude NA for mean calculation)   ####################
  
  # for CI (and for missing responses), SCR values in table are stimulus-locked to the GNG target; for FA, FH, SH the SCR values are response-locked -> see script "Preprocess_SCR_Export_Files"! there I assigned stimulus-locked SCR values to correctly inhibited trials and trials with missing response (makes no sense to take response-locked SCR value in these, because there is no response and the evaluation trigger is send late)
  mean_SCR_after_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & (single_trial_data$word_resp == 51 | single_trial_data$word_resp == 52) & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$scr_amplitude_z_score_log,na.rm=TRUE)
  mean_SCR_after_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & (single_trial_data$word_resp == 51 | single_trial_data$word_resp == 52) & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$scr_amplitude_z_score_log,na.rm=TRUE)
  mean_SCR_after_CI <- mean(single_trial_data[single_trial_data$response_type == "CI" & (single_trial_data$word_resp == 51 | single_trial_data$word_resp == 52) & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$scr_amplitude_z_score_log,na.rm=TRUE)
  mean_SCR_after_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & (single_trial_data$word_resp == 51 | single_trial_data$word_resp == 52) & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$scr_amplitude_z_score_log,na.rm=TRUE)
  
  
  
  
  
  
  
  ####################   add df of one subject to df containing all subjects  ####################
  
  df4save <- rbind(df4save, data.frame(subject,mean_neg_after_FA,mean_pos_after_FA,mean_neg_after_FH,mean_pos_after_FH,mean_neg_after_CI,mean_pos_after_CI,mean_neg_after_SH,mean_pos_after_SH,mean_priming_overall,mean_priming_after_FA,mean_priming_after_FH,median_neg_after_FA,median_pos_after_FA,median_neg_after_FH,median_pos_after_FH,median_neg_after_CI,median_pos_after_CI,median_neg_after_SH,median_pos_after_SH,median_priming_overall,median_priming_after_FA,median_priming_after_FH,count_neg_after_FA,count_pos_after_FA,count_neg_after_FH,count_pos_after_FH,count_outlier_words_FA_FH,percent_correct_neg_after_FA,percent_correct_pos_after_FA,percent_correct_neg_after_FH,percent_correct_pos_after_FH,percent_correct_neg_after_CI,percent_correct_pos_after_CI,percent_correct_neg_after_SH,percent_correct_pos_after_SH,GNG_mean_FA,GNG_mean_FH,GNG_mean_SH,GNG_percent_FA,GNG_percent_FH,GNG_percent_CI,GNG_percent_SH,GNG_percent_miss, mean_SCR_after_FA,mean_SCR_after_FH,mean_SCR_after_CI,mean_SCR_after_SH,count_outlier_GNG))
  data4mixedmodels <- rbind(data4mixedmodels,single_trial_data)
  
  
  }
  
  
  
  
  ###################   Prepare data4mixedmodels for LMMs: Define Factors   ####################
 
  data4mixedmodels$condition     <- paste(data4mixedmodels$valence,"_after_",data4mixedmodels$response_type)
  data4mixedmodels$response_type <- factor(data4mixedmodels$response_type, levels=c("FA","FH","CI","SH"))         # reorder factor levels! (this is important, because first level is used as reference factor in LMM)
  data4mixedmodels$valence       <- factor(data4mixedmodels$valence)                                              # automatic factor order is alphabetical (1 = neg, 2 = pos)
  data4mixedmodels$subjectID     <- factor(data4mixedmodels$subjectID)
 
 
  
  
  
  
  
  
  
  ###################   Read in questionnaire data, merge with aggregated data   ####################
  
  setwd("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses")
  questionnaires <- list.files(pattern = ".sav")                                                         # make sure that only one .sav file (= the current PEQ export) exists there!
  questionnaires <- read.spss(questionnaires, to.data.frame = TRUE, add.undeclared.levels = "no")        # add.undeclared.levels required to prevent irrelevant warning message
  questionnaires <- questionnaires[,c("CODE", "BD2SUMT0","BASO00T0","BASO01T0","BASO02T0","BASO03T0","BASO04T0","FMPO00T0","FMPO01T0","FMPO02T0","FMPO03T0","FMPO04T0","NNGO00T0","NNGO01T0","OCISUMT0", "OCIO00T0", "OCIO01T0","OCIO02T0","OCIO03T0","OCIO04T0","OCIO05T0","PANO00T0","PANO01T0","PSWSUMT0","STSSUMT0","STTSUMT0","TCISUMT0","TCIO00T0","TCIO01T0","TCIO02T0","TCIO03T0","WSTSUMT0")] # only keep relevant columns (I chose to keep sum scores instead of mean scores, because for BDI and OCI only sum scores are available and missing data are not possible because data aquired with tablet; except for NEO and BIS/BAS, according to manuals, scoring is based on sum scores; JK also used sum scores for NEO, BIS/BAS)
  # questionnaires <- questionnaires[grep("behav", questionnaires$CODE), ]                               # only keep participants from behavioral study
  questionnaires$CODE <- paste0(questionnaires$CODE,".txt")                                              # convert code for joining (make it identical to code in df4save)                               
  colnames(questionnaires) <-c("CODE", "BDI_Total", "BIS_Total", "BAS_Total", "BAS_Drive", "BAS_Fun_Seeking", "BAS_Reward_Responsiveness", "FMPS_CMD", "FMPS_PEC", "FMPS_PST", "FMPS_ORG", "FMPS_PER/Total", "NEO_Neuroticism", "NEO_Conscientiousness", "OCI_Total", "OCI_Washing", "OCI_Checking", "OCI_Ordering", "OCI_Obsessions", "OCI_Hoarding", "OCI_Neutralising", "PANAS_Pos", "PANAS_Neg", "PSWQ_Total", "STAI_State", "STAI_Trait", "TCI_Harm_Avoidance_Total", "TCI_Anticipatory_worry", "TCI_Fear_of_uncertainty", "TCI_Shyness", "TCI_Fatigability", "WST_Total")   # rename columns 
  

  df4save <- left_join(df4save,questionnaires, by = c("subject" = "CODE"))
  
  
  
  
  
  
  
 
  ###################   Save df4save and data4mixedmodels   ####################

  setwd("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses")    # setting a different folder as working directory to prevent saving stuff into the folder containing the logfiles
  date_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  filename_aggregated <- paste("Data_Aggregated_For_",length(logfiles),"_subjects_",date_time, ".rda", sep = "")
 # save(df4save,file=filename_aggregated)  
  filename_singletrial <- paste("Data_Single_Trial_For_",length(logfiles),"_subjects_",date_time, ".rda", sep = "")
 # save(data4mixedmodels,file=filename_singletrial)    