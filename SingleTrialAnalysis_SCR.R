  ##### Proprocessing of Evaluative GNG Task in Behavioral Experiment on Affective Priming Effect 2018 ####
  ##### by Luisa Balzus
  ##### 09.11.2018
  
  
  library(dplyr)
  library(foreign)
  library(openxlsx)

  
  # clear environment
  rm(list=ls())
  
  # force R to not use exponential notation
  options(scipen = 999)
  
  # create empty data frames to write aggregated data and single trial data in it
  data4mixedmodels <- data.frame()   
  df4save <- data.frame()   
  
  
  
  ###################   Loop over Subjects to Preprocess Logfile Data and SCR Data for Statistical Analyses   ####################
  
  # Load logfiles and start loop
  logfiles <- list.files("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_Behavioral", pattern = ".txt")       
  
  for (subject in logfiles){                                                                               # loop reading txt-files as table, ommit first 58 lines and added lines after trials; use first line as header
    setwd("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_Behavioral")                 # path to folder containing the log files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 
    raw_log <- read.table(subject, skip = 58, fill = TRUE, header = TRUE, nrows = 516)
    name_of_subject <-  gsub("\\.txt*","",subject)
    rating <- read.table(subject, skip = 575, fill = TRUE, header = TRUE, sep = ":", stringsAsFactors = FALSE)
    
    
  # Load SCR files 
  setwd("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/9_SCR_Export_Preprocessed")               # path to folder containing the scr files (use of forward slashes instead of backward slashes is required); should contain ONLY logfiles 
  filename_scr <- paste0(name_of_subject,"_SCR_Export_era.xlsx")
  raw_scr      <- read.xlsx(filename_scr)
  
  

  
  ####################   merge log and scr data   ####################
 
  # show progress
  message('Processing ', subject, ' of ', length(logfiles))
  
  
  # return all rows of log, merged with matching columns of scr; rows in log with no match in scr will have NAs in the new columns; if there are multiple matches between log and scr, all combinations of the matches are returned
  raw_scr$Participant.ID <- as.numeric(gsub("\\D", "", raw_scr$Participant.ID))   # necessary for merging that types of number and Participant.ID are identical (extract number from character Participant.ID)
  data_log_scr           <- left_join(raw_log,raw_scr,by = c("number"= "Participant.ID","count" = "Trial.ID","t_type"="Event.ID.GNG","t_resp"="Event.ID.GNG_Resp","w_type"="Event.ID.Word","w_resp"="Event.ID.Word_Resp"))
  
  
  # check whether there were rows in logfile with no match in scr (is true for subject 25 and 29, because SCR recording was interrupted)
  if (any(is.na(data_log_scr$CDA.ISCR.GNG_Resp))){
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
  single_trial_data <- subset(data_log_scr, select = c("number","count","t_type","t_resp","t_rt","w_type","w_resp","w_rt","word","CDA.ISCR.GNG","CDA.ISCR.GNG_Resp","CDA.ISCR.Word","CDA.ISCR.Word_Resp"))
  
  
  # rename columns
  single_trial_data <- rename(single_trial_data,
                              subjectID          = number,
                              trial              = count,
                              gng_condition      = t_type,
                              gng_resp           = t_resp,
                              gng_rt             = t_rt,
                              word_condition     = w_type,
                              word_resp          = w_resp,
                              word_rt            = w_rt,
                              iscr_gng_stimulus  = CDA.ISCR.GNG,
                              iscr_gng_resp      = CDA.ISCR.GNG_Resp,
                              iscr_word_stimulus = CDA.ISCR.Word,
                              iscr_word_resp     = CDA.ISCR.Word_Resp)
  

  
    ####################   create new word accuracy variable   ####################
  
  single_trial_data$word_accuracy  <- single_trial_data$word_resp
  single_trial_data$word_accuracy[single_trial_data$word_accuracy == 51 | single_trial_data$word_accuracy == 52]  <- 1 # correct responses
  single_trial_data$word_accuracy[single_trial_data$word_accuracy == 53 | single_trial_data$word_accuracy == 54]  <- 0 # incorrect responses
  
  
  
  ####################   define conditions   ####################
  
  single_trial_data <- mutate(single_trial_data, 
                     valence = ifelse(word_condition > 200, "pos", "neg"),                                                      
                     response_type = ifelse((gng_resp == 43 | gng_resp == 44), "FA", 
                                            ifelse(gng_resp == 41, "FH", 
                                                   ifelse(gng_resp == 45 | gng_resp == 46, "CI", 
                                                          ifelse(gng_resp == 42 , "SH",
                                                                 "Miss_or_False_Key"))))
  )
  
 
  single_trial_data$condition     <- paste0(single_trial_data$valence,"_after_",single_trial_data$response_type)
  
  
  
  ####################   inverse transformation word rt and gng_rt (for normality)     ####################  
  
  
  single_trial_data$word_rt_inverse <- -1000/single_trial_data$word_rt                               # inverse transformation word rt
  single_trial_data$gng_rt_inverse  <- -1000/single_trial_data$gng_rt                                # inverse transformation gng rt
  
  
  
  
  
  ####################   define GNG  outliers   ####################   
  
  # assign TRUE to column GNG_outlier < 100 or > 700 ms; here, trials with GNG outliers are excluded ALSO for word categorization in corresponding trial! (Pourtois used 150 and 500 ms, but unclear whether trials with GNG outliers were excluded for word classification there, too; I do not want to loose too many words, so I use 100 and 700 ms) (priming after FA is 119 if no gng outliers are excluded, 116 if 100/700 ms, 104 if 150/500 ms; priming after FH remains always the same)
  single_trial_data$gng_invalid_rt <- FALSE
  single_trial_data$gng_invalid_rt[(single_trial_data$gng_resp < 45) & single_trial_data$gng_rt < 100 | single_trial_data$gng_rt > 700]  <- TRUE
  
  
  
  
  
  
  ####################   calculate GNG RT   #################### 
  
  # FA = False Alarm; FH = Fast Hit; CI = Correctly Inhibited; SH = Slow Hit
  # RT not calculated for CI, because there is no rt
  GNG_rt_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
  GNG_rt_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
  GNG_rt_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
  
  
  
  
  
  ####################   calculate GNG % response types   #################### 
  
  GNG_count_FA   <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
  GNG_count_FH   <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
  GNG_count_CI   <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)
  GNG_count_SH   <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$gng_invalid_rt == FALSE,]$gng_rt)

  GNG_all_resp_without_misses_outliers_wrong_keys <- GNG_count_FA + GNG_count_FH + GNG_count_CI + GNG_count_SH 
  
  GNG_percent_FA   <- (GNG_count_FA/GNG_all_resp_without_misses_outliers_wrong_keys)*100
  GNG_percent_FH   <- (GNG_count_FH/GNG_all_resp_without_misses_outliers_wrong_keys)*100
  GNG_percent_CI   <- (GNG_count_CI/GNG_all_resp_without_misses_outliers_wrong_keys)*100
  GNG_percent_SH   <- (GNG_count_SH/GNG_all_resp_without_misses_outliers_wrong_keys)*100

  GNG_percent_outlier <- ((length(single_trial_data$gng_invalid_rt[single_trial_data$gng_invalid_rt == TRUE]))/(GNG_count_FA + GNG_count_FH + GNG_count_SH))*100  
  
  
  
  
  ####################   define word outliers   ####################   

  # make subset of data to define the outlier threshold (MAD) separately for all categories
  # FA = False Alarm; FH = Fast Hit; CI = Correctly Inhibited; SH = Slow Hit; subsetting only required for outlier detection
  neg_after_FA <- subset(single_trial_data, valence == "neg" & response_type == "FA" & word_resp == 51 & gng_invalid_rt == FALSE)
  pos_after_FA <- subset(single_trial_data, valence == "pos" & response_type == "FA" & word_resp == 52 & gng_invalid_rt == FALSE)
  neg_after_FH <- subset(single_trial_data, valence == "neg" & response_type == "FH" & word_resp == 51 & gng_invalid_rt == FALSE)
  pos_after_FH <- subset(single_trial_data, valence == "pos" & response_type == "FH" & word_resp == 52 & gng_invalid_rt == FALSE)
  neg_after_CI <- subset(single_trial_data, valence == "neg" & response_type == "CI" & word_resp == 51 & gng_invalid_rt == FALSE)
  pos_after_CI <- subset(single_trial_data, valence == "pos" & response_type == "CI" & word_resp == 52 & gng_invalid_rt == FALSE)
  neg_after_SH <- subset(single_trial_data, valence == "neg" & response_type == "SH" & word_resp == 51 & gng_invalid_rt == FALSE)
  pos_after_SH <- subset(single_trial_data, valence == "pos" & response_type == "SH" & word_resp == 52 & gng_invalid_rt == FALSE)
  
  
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
  
  count_outlier_words_FA_FH     <- length(single_trial_data$outlier_words[single_trial_data$outlier_words == TRUE & (single_trial_data$response_type == "FA" | single_trial_data$response_type == "FH" )])  # count number of outliers
  percent_outlier_words_overall <- ((length(single_trial_data$outlier_words[single_trial_data$outlier_words == TRUE & (single_trial_data$response_type == "FA" | single_trial_data$response_type == "FH" | single_trial_data$response_type == "SH" | single_trial_data$response_type == "CI" )]))/516)*100
  
  
  # NOTE: this whole approach for outlier denifition may not be the best strategy, for 3 reasons:
  #   - for defining the outlier threshold (MAD) I here excluded word responses after invalid gng RT; this could be criticized, but is not so problematic that I need to change it now after all analysis were reported (we could argue that after GNG invalid RTs the word response is also not normal) -> do not change this anymore, as this would change all outputs
  #   - the incorrectly categorized words are automatically defined as non-outliers (not sure whether this is good or not) 
  #       - this is not problematic for all (RT) analyses, as I exclude these trials for these analyses anyway 
  #       - but I do NOT exclude incorrect words for Accuracy analysis (GLMM) - so it may DO a bit of a difference there, right? But I still think that defining outliers among incorrect responses isn't really necessary either, as these responses are already incorrect, no matter what their RT is
  #       - for calculation of proportion word outliers it might actually be good, because there I want to see how many among the correctly categorized words were outliers, right?
  #   - I here apply the threshold (MAD) calculated for words after SHs also to the category words after Miss_or_False_Key 
  #       - this is not problematic, as I exclude these trials for all analyses
  #       - for calculation of proportion word outliers it might actually be good, because in the Miss_or_Wrong_Key category there are not enough events to calculate a reliable threshold, but in the SH there are enough events, as this is the biggest category
   
  # a better approach would be: 
  #     single_trial_data <- single_trial_data %>% group_by(subjectID,condition,word_accuracy) %>% mutate(outlier_words = ifelse(abs(word_rt - median(word_rt))/mad(word_rt)>3, TRUE, FALSE))
  
  # I did not apply this new approach here retrospectively, because I want the analysis outputs remain unchanged, but excluding gng_invalid_rt trials with the new approach to define the same outliers as with the old approach causes problems 
  #     1) to get the same outlier thresholds (MADs) with this approach as with the old approach, I would need to exclude the gng_invalid_rt trials before marking the outliers and then merge the overall table with the table including the outlier definition 
  #         - the good news is that then, among the correctly categorized words after valid GNG responses, the outliers are defined correctly both in this approach and in the newer one (see below) -> this is all what matters for all subsequent analyses
  #         - but: I would not get an outlier decision for the gng_invalid rt trials (would enter NA there) because I excluded these before identifying word outliers, which would be a disadvantage when calculating the proportion of word outliers (and the NAs may seem weird when publishing the data)
  #     2) with this approach, I would also get NAs for all specific categories with only one event (combination subject+condition+word_accuracy, e.g. incorrect classified neg word after Miss_or_False_Key in Subject 3)

  # the two methods come to different outlier decisions for words after Miss_or_False_Key gng responses and for incorrectly categorized words (due to reasons stated above) 
  #   - this makes no difference for my analyses, as I exclude these events threre anyway (at least for RT analyses, not for Accuracy analysis)
  #   - but it makes a difference for the calculation of proportion of word outliers
  #   - additionally, I would get NAs in the outlier column when using the new approach but try to use the old outlier threshold (calculated excluding gng_invalid_rt trials) which might be a disadvantage when publishing the data

  # CONCLUSION: so I go with the old approach for publishing the code and for my next project, I will use the new approach (but be aware that than outliers are defined among incorrect responses as well - not sure whether this is good or not; and after Miss_or_Wrong_Key there may be some NAs)
  
  
  
  

  ####################   calculat priming effect  word classification with mean    ####################

  rt_neg_after_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  rt_pos_after_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  rt_neg_after_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  rt_pos_after_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  rt_neg_after_CI <- mean(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  rt_pos_after_CI <- mean(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  rt_neg_after_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)
  rt_pos_after_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$outlier_words == FALSE & single_trial_data$gng_invalid_rt == FALSE,]$word_rt)

  rt_priming_after_FA <-  rt_pos_after_FA - rt_neg_after_FA
  rt_priming_after_FH <-  rt_neg_after_FH - rt_pos_after_FH
  rt_priming_overall  <- (rt_pos_after_FA + rt_neg_after_FH) - (rt_neg_after_FA + rt_pos_after_FH)


  # when calculate median: without exclusion of word outliers and gng outliers




  ####################   count number of events word classification   ####################

  # count number of events,  outliers are removed; if I want to work with median, do not exclude outlier
  correct_responses_neg_after_FA        <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  correct_responses_pos_after_FA        <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  correct_responses_neg_after_FH        <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  correct_responses_pos_after_FH        <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  correct_responses_neg_after_CI        <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  correct_responses_pos_after_CI        <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  correct_responses_neg_after_SH        <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  correct_responses_pos_after_SH        <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)

  incorrect_responses_neg_after_FA <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  incorrect_responses_pos_after_FA <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  incorrect_responses_neg_after_FH <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  incorrect_responses_pos_after_FH <- length(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  incorrect_responses_neg_after_CI <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  incorrect_responses_pos_after_CI <- length(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  incorrect_responses_neg_after_SH <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "neg" & single_trial_data$word_resp == 53 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)
  incorrect_responses_pos_after_SH <- length(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$valence == "pos" & single_trial_data$word_resp == 54 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE,]$word_rt)





  ####################   calculate accuracy word classification     ####################

  # misses are not contained in overall number of events -> errors & correct sum up to 100 %
  percent_correct_responses_neg_after_FA <- correct_responses_neg_after_FA / (correct_responses_neg_after_FA + incorrect_responses_neg_after_FA) * 100
  percent_correct_responses_pos_after_FA <- correct_responses_pos_after_FA / (correct_responses_pos_after_FA + incorrect_responses_pos_after_FA) * 100
  percent_correct_responses_neg_after_FH <- correct_responses_neg_after_FH / (correct_responses_neg_after_FH + incorrect_responses_neg_after_FH) * 100
  percent_correct_responses_pos_after_FH <- correct_responses_pos_after_FH / (correct_responses_pos_after_FH + incorrect_responses_pos_after_FH) * 100
  percent_correct_responses_neg_after_CI <- correct_responses_neg_after_CI / (correct_responses_neg_after_CI + incorrect_responses_neg_after_CI) * 100
  percent_correct_responses_pos_after_CI <- correct_responses_pos_after_CI / (correct_responses_pos_after_CI + incorrect_responses_pos_after_CI) * 100
  percent_correct_responses_neg_after_SH <- correct_responses_neg_after_SH / (correct_responses_neg_after_SH + incorrect_responses_neg_after_SH) * 100
  percent_correct_responses_pos_after_SH <- correct_responses_pos_after_SH / (correct_responses_pos_after_SH + incorrect_responses_pos_after_SH) * 100

  percent_correct_overall      <- (correct_responses_neg_after_FA+correct_responses_pos_after_FA+correct_responses_neg_after_FH+correct_responses_pos_after_FH+
                                   correct_responses_neg_after_CI+correct_responses_pos_after_CI+correct_responses_neg_after_SH+correct_responses_pos_after_SH)/
                                  (correct_responses_neg_after_FA + incorrect_responses_neg_after_FA+
                                   correct_responses_pos_after_FA + incorrect_responses_pos_after_FA+
                                   correct_responses_neg_after_FH + incorrect_responses_neg_after_FH+
                                   correct_responses_pos_after_FH + incorrect_responses_pos_after_FH+
                                   correct_responses_neg_after_CI + incorrect_responses_neg_after_CI+
                                   correct_responses_pos_after_CI + incorrect_responses_pos_after_CI+
                                   correct_responses_neg_after_SH + incorrect_responses_neg_after_SH+
                                   correct_responses_pos_after_SH + incorrect_responses_pos_after_SH) * 100

  # calculate priming effect for accuracy
  accuracy_priming_after_FA <- percent_correct_responses_neg_after_FA - percent_correct_responses_pos_after_FA
  accuracy_priming_after_FH <- percent_correct_responses_pos_after_FH - percent_correct_responses_neg_after_FH
  accuracy_priming_overall  <- accuracy_priming_after_FA + accuracy_priming_after_FH
  
  
  
  
  
  ####################   reorder columns     ####################
  
  single_trial_data <- single_trial_data[,c("subjectID","trial","gng_condition","gng_resp","gng_rt","gng_rt_inverse","gng_invalid_rt","word_condition","word_resp","word_rt","word_rt_inverse","outlier_words","word","word_accuracy","valence","response_type","condition","iscr_gng_stimulus","iscr_gng_resp","iscr_word_stimulus","iscr_word_resp")] 
  
  
  
  

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
  


  ####################   log/square root transform and z-standardize scr     ####################

 # previously, I did the scaling based on all trials. but then I have an almost sign. correlation between SCR after FA and rt_priming_after_FA (p = 0.054 for log transformed data / P = 0.073 for square root transformed data)
 # single_trial_data$iscr_gng_resp_log_z_score      <- scale((log(single_trial_data$iscr_gng_resp + 1)), center = TRUE, scale = TRUE)
 # single_trial_data$iscr_gng_resp_sqrt_z_score     <- scale((sqrt(single_trial_data$iscr_gng_resp)), center = TRUE, scale = TRUE)

  
 # scaling based only on trials that are used for SCR analysis later; if I change exclusion criteria there, I also have to change them here as well
  single_trial_data$iscr_gng_resp_log_z_score  <- NA 
  single_trial_data$iscr_gng_resp_sqrt_z_score <- NA 
 
  single_trial_data[single_trial_data$gng_resp <= 44  & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$word_resp <= 52 & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score    <-  scale((log(single_trial_data[single_trial_data$gng_resp <= 44  & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$word_resp <= 52 & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp + 1)), center = TRUE, scale = TRUE)
  single_trial_data[single_trial_data$gng_resp <= 44  & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$word_resp <= 52 & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_sqrt_z_score   <- scale((sqrt(single_trial_data[single_trial_data$gng_resp <= 44  & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$word_resp <= 52 & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp)), center = TRUE, scale = TRUE)
  
  # NOTE: The z-standardization of the SCR (variable iscr_gng_resp_sqrt_z_score) is based on only trials that also enter the analyses, 
  # not on all trials. This means, that e.g. word and gng outliers, trials with miss or wrong key in word/gng response, incorrect word responses, 
  # and trials followed_or_preceded_by_FA_or_wrong_key were not included when scaling the SCR. I once tried basing it on all trials. The only things 
  # that change are the values of the SCR ANOVAs and Raw-Barplots based on this variable. The results of the LMM on the Relation Priming - SCR remain 
  # the same, only some decimals are slightly different. I thus decided to stick with the scaling based on only the trials that enter the analysis, 
  # because distribution of the variable is then a it nicer (when basing it on all trials, I again get a peak at the value of 0) and basing the scaling 
  # on all trials does not reflect the true value better. For the manuscript, I calculate the variable iscr_gng_resp_sqrt_z_score first when calculating 
  # the LMM Priming - SCR, to avoid having the variable containing many NAs in the Data I publish.*
    
  
  
  #################### for SCR Analyses, check whether number of correctly classified pos / neg words after FA are equal ####################

  # I decided to remove GNG outliers here and word outliers (consistent for all analyses); SCR in a trial may be influenced by invalid GNG response times but not by invalid word response times?
  for_SCR_count_neg_after_FA <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE & !is.na(single_trial_data$iscr_gng_resp),]$word_resp)
  for_SCR_count_pos_after_FA <- length(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE & !is.na(single_trial_data$iscr_gng_resp),]$word_resp)

  for_SCR_count_neg_before_FA <- 0
  for_SCR_count_pos_before_FA <- 0
  for_SCR_count_neg_after_FA_next_trial <- 0
  for_SCR_count_pos_after_FA_next_trial <- 0

  for (i in 2:(nrow(single_trial_data)-1)) {
    current_row <- single_trial_data[i,]
    next_row <- single_trial_data[(i+1),]
    previous_row <- single_trial_data[(i-1),]
    if (current_row$response_type == "FA" & previous_row$word_resp == 51 & current_row$gng_invalid_rt == FALSE & current_row$outlier_words == FALSE & current_row$followed_or_preceded_by_FA_or_wrong_key == FALSE & !is.na(current_row$iscr_gng_resp)) {
      for_SCR_count_neg_before_FA <- for_SCR_count_neg_before_FA + 1
    }
    if (current_row$response_type == "FA" & previous_row$word_resp == 52 & current_row$gng_invalid_rt == FALSE & current_row$outlier_words == FALSE & current_row$followed_or_preceded_by_FA_or_wrong_key == FALSE & !is.na(current_row$iscr_gng_resp)) {
      for_SCR_count_pos_before_FA <- for_SCR_count_pos_before_FA + 1
    }
    if (current_row$response_type == "FA" & next_row$word_resp == 51     & current_row$gng_invalid_rt == FALSE & current_row$outlier_words == FALSE & current_row$followed_or_preceded_by_FA_or_wrong_key == FALSE & !is.na(current_row$iscr_gng_resp)) {
      for_SCR_count_neg_after_FA_next_trial <- for_SCR_count_neg_after_FA_next_trial + 1
    }
    if (current_row$response_type == "FA" & next_row$word_resp == 52     & current_row$gng_invalid_rt == FALSE & current_row$outlier_words == FALSE & current_row$followed_or_preceded_by_FA_or_wrong_key == FALSE & !is.na(current_row$iscr_gng_resp)) {
      for_SCR_count_pos_after_FA_next_trial <- for_SCR_count_pos_after_FA_next_trial + 1
    }
  }





  # ####################   calculate priming facilitation score for LMM   #####################
  # 
  # # makes not much sense, because it does not account for effects of word length and frequency on rt...
  # #
  # # single_trial_data$facilitation_score <- NA      # score shall be higher the more priming a person shows in a trial (i.e. fast for neg words after FA, slow for pos words; fast for pos words after FH, slow for neg words)
  # # 
  # # for (j in 1:nrow(single_trial_data)) {
  # #   if (single_trial_data[j,]$response_type == "FH" & single_trial_data[j,]$valence == "neg")  {
  # #   single_trial_data[j,]$facilitation_score <- single_trial_data[j,]$word_rt - rt_neg_after_FH
  # #   }
  # #   if (single_trial_data[j,]$response_type == "FH" & single_trial_data[j,]$valence == "pos")  {
  # #     single_trial_data[j,]$facilitation_score <- rt_pos_after_FH - single_trial_data[j,]$word_rt
  # #   }
  # #   if (single_trial_data[j,]$response_type == "FA" & single_trial_data[j,]$valence == "neg")  {
  # #     single_trial_data[j,]$facilitation_score <- rt_neg_after_FA - single_trial_data[j,]$word_rt
  # #   }
  # #   if (single_trial_data[j,]$response_type == "FA" & single_trial_data[j,]$valence == "pos")  {
  # #     single_trial_data[j,]$facilitation_score <- single_trial_data[j,]$word_rt - rt_pos_after_FA
  # #   }
  # # }
  # # # I don't z-standardize facilitation score, because subjects are entered as random intercept in LMM anyway; I don't log-transform facilitation score, because it is normally distributed -> hist(single_trial_data$facilitation_score)
  # 
  # 
  # 
  # 
  # 
  # 
  ####################   calculate mean SCR in conditions (some subjects contain NAs, because SCR recording was disrupted -> exclude NA for mean calculation)   ####################

  # for CI (and for missing responses), SCR values in table are stimulus-locked to the GNG target; for FA, FH, SH the SCR values are response-locked -> see script "Preprocess_SCR_Export_Files"! there I assigned stimulus-locked SCR values to correctly inhibited trials and trials with missing response (makes no sense to take response-locked SCR value in these, because there is no response and the evaluation trigger is send late)
  # exclude trials with incorrect word categorization or wrong key or miss in word categorization
  # exclude trials that were preceded by false alarm or wrong key in GNG or by incorrect word categorization or wrong key in word categorization 
  # exclude GNG and word outliers (consistent for all analyses) (after a premature response, the subject may not notice whether the answer was correct or incorrect)

  SCR_after_FA     <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$word_resp <= 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)
  SCR_after_FH     <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$word_resp <= 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)
  SCR_after_CI     <- mean(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$word_resp <= 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)
  SCR_after_SH     <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$word_resp <= 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)

  SCR_neg_after_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)
  SCR_neg_after_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)
  SCR_neg_after_CI <- mean(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)
  SCR_neg_after_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$word_resp == 51 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)

  SCR_pos_after_FA <- mean(single_trial_data[single_trial_data$response_type == "FA" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)
  SCR_pos_after_FH <- mean(single_trial_data[single_trial_data$response_type == "FH" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)
  SCR_pos_after_CI <- mean(single_trial_data[single_trial_data$response_type == "CI" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)
  SCR_pos_after_SH <- mean(single_trial_data[single_trial_data$response_type == "SH" & single_trial_data$word_resp == 52 & single_trial_data$gng_invalid_rt == FALSE & single_trial_data$outlier_words == FALSE & single_trial_data$followed_or_preceded_by_FA_or_wrong_key == FALSE,]$iscr_gng_resp_log_z_score,na.rm=TRUE)



  ####################   add rating variables  #################### 
  
  effort            <- as.numeric(rating[1,1])
  error_avoidance   <- as.numeric(rating[2,1])
  error_frustration <- as.numeric(rating[3,1])
  fatigue           <- as.numeric(rating[4,1])
  
  
  
  ####################   add df of one subject to df containing all subjects  ####################
  
  df4save <- rbind(df4save, data.frame(subject,rt_neg_after_FA,rt_pos_after_FA,rt_neg_after_FH,rt_pos_after_FH,rt_neg_after_SH,rt_pos_after_SH,rt_neg_after_CI,rt_pos_after_CI,rt_priming_overall,rt_priming_after_FA,rt_priming_after_FH,correct_responses_neg_after_FA,correct_responses_pos_after_FA,correct_responses_neg_after_FH,correct_responses_pos_after_FH,count_outlier_words_FA_FH,percent_outlier_words_overall,percent_correct_responses_neg_after_FA,percent_correct_responses_pos_after_FA,percent_correct_responses_neg_after_FH,percent_correct_responses_pos_after_FH,percent_correct_responses_neg_after_SH,percent_correct_responses_pos_after_SH,percent_correct_responses_neg_after_CI,percent_correct_responses_pos_after_CI,percent_correct_overall,accuracy_priming_overall,accuracy_priming_after_FA,accuracy_priming_after_FH,GNG_rt_FA,GNG_rt_FH,GNG_rt_SH,GNG_percent_FA,GNG_percent_FH,GNG_percent_CI,GNG_percent_SH,GNG_percent_outlier,SCR_after_FA,SCR_after_FH,SCR_after_SH,SCR_after_CI,SCR_neg_after_FA,SCR_neg_after_FH,SCR_neg_after_SH,SCR_neg_after_CI,SCR_pos_after_FA,SCR_pos_after_FH,SCR_pos_after_SH,SCR_pos_after_CI,for_SCR_count_neg_after_FA,for_SCR_count_pos_after_FA,for_SCR_count_neg_before_FA,for_SCR_count_pos_before_FA,for_SCR_count_neg_after_FA_next_trial,for_SCR_count_pos_after_FA_next_trial,effort,error_avoidance,error_frustration,fatigue))
  data4mixedmodels <- rbind(data4mixedmodels,single_trial_data)
  
  }
  
  
  
  
  
  
  
  ###################   Read in questionnaire data, merge with aggregated data   ####################
  
  setwd("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses")
  questionnaires           <- list.files(pattern = ".sav")                                                         # make sure that only one .sav file (= the current PEQ export) exists there!
  questionnaires           <- read.spss(questionnaires, to.data.frame = TRUE, add.undeclared.levels = "no")        # add.undeclared.levels required to prevent irrelevant warning message
  questionnaires           <- questionnaires[,c("CODE", "BD2SUMT0","BASO00T0","BASO01T0","BASO02T0","BASO03T0","BASO04T0","FMPO00T0","FMPO01T0","FMPO02T0","FMPO03T0","FMPO04T0","NNGO00T0","NNGO01T0","OCISUMT0", "OCIO00T0", "OCIO01T0","OCIO02T0","OCIO03T0","OCIO04T0","OCIO05T0","PANO00T0","PANO01T0","PSWSUMT0","STSSUMT0","STTSUMT0","TCISUMT0","TCIO00T0","TCIO01T0","TCIO02T0","TCIO03T0","WSTSUMT0")] # only keep relevant columns (I chose to keep sum scores instead of mean scores, because for BDI and OCI only sum scores are available and missing data are not possible because data aquired with tablet; except for NEO and BIS/BAS, according to manuals, scoring is based on sum scores; JK also used sum scores for NEO, BIS/BAS)
  questionnaires$CODE      <- paste0(questionnaires$CODE, ".txt")                                                  # convert code for joining (make it identical to code in df4save)                               
  colnames(questionnaires) <- c("CODE", "BDI_Total", "BIS_Total", "BAS_Total", "BAS_Drive", "BAS_Fun_Seeking", "BAS_Reward_Responsiveness", "FMPS_CMD", "FMPS_PEC", "FMPS_PST", "FMPS_ORG", "FMPS_PER/Total", "NEO_Neuroticism", "NEO_Conscientiousness", "OCI_Total", "OCI_Washing", "OCI_Checking", "OCI_Ordering", "OCI_Obsessions", "OCI_Hoarding", "OCI_Neutralising", "PANAS_Pos", "PANAS_Neg", "PSWQ_Total", "STAI_State", "STAI_Trait", "TCI_Harm_Avoidance_Total", "TCI_Anticipatory_worry", "TCI_Fear_of_uncertainty", "TCI_Shyness", "TCI_Fatigability", "WST_Total")   # rename columns 
  df4save                  <- left_join(df4save,questionnaires, by = c("subject" = "CODE"))
  


 
  ###################   Save df4save and data4mixedmodels   ####################

  setwd("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses")    # setting a different folder as working directory to prevent saving stuff into the folder containing the logfiles
  date_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

  filename_singletrial <- paste("Data_Single_Trial_For_",length(logfiles),"_subjects_",date_time, ".rda", sep = "")
  save(data4mixedmodels,file=filename_singletrial)

  filename_traits <- paste("Data_Trait_Variables_For_",length(logfiles),"_subjects_",date_time, ".rda", sep = "")
  save(questionnaires,file=filename_traits)

  filename_aggregated <- paste("Data_only_for_overview.rda")
  save(df4save,file=filename_aggregated) 
  
  
  