##### Implement necessarry changes in SCR export files for further Analysis (Experiment on Affective Priming Effect 2018) ##### 
##### by Luisa Balzus
##### 11.11.2018


# clear environment
rm(list=ls())



####################     load     #######################################
####################   scr data   #######################################

# Load SCR files (Ledalab output exported as .txt) 
scrfiles <- list.files("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/8_SCR_Export", pattern = ".txt")             # path to folder containing the scr files (use of forward slashes instead of backward slashes is required); should contain ONLY scr-files 
setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/8_SCR_Export")                              


for (subject in scrfiles){                                                                              # loop reading csv-file
  raw_scr <- read.table(subject, header = TRUE)
  subset_scr <- subset(raw_scr, (Event.NID >= 21 & Event.NID <= 58) |(Event.NID >= 141 & Event.NID <= 149) |(Event.NID >= 241 & Event.NID <= 249))   # only keep required markers (responses to GNG target and words)
  
  
  
  ####################   rearrange   #######################################
  ####################   SCR data    #######################################
  
  # create new column that represents correct order of triggers within one trial
  subset_scr$order[subset_scr$Event.NID <=23]                            <-1
  subset_scr$order[subset_scr$Event.NID >=40 & subset_scr$Event.NID <=49]<-2
  subset_scr$order[subset_scr$Event.NID >=140 ]                           <-3
  subset_scr$order[subset_scr$Event.NID >=50 & subset_scr$Event.NID <=59]<-4
  
  number_chunks <- nrow(subset_scr)/4                                                              # one chunk is one trial; one trial contains 4 triggers in scr
  scr_ordered <- data.frame()  
  
  for (i in 0:(number_chunks-1)){
    chunk_start <- 1+4*i
    chunk_end <- 4+4*i
    chunk <- subset_scr[chunk_start:chunk_end,]                                              # extract rows 1-4 of current chunk from scr data frame        
    chunk_ordered <- chunk[order(chunk$order),]                                              # order rows in chunk by the column "order"
    chunk_line <- c(chunk_ordered[1,],chunk_ordered[2,],chunk_ordered[3,],chunk_ordered[4,]) # write correctly ordered rows in single row
    scr_ordered <- rbind(scr_ordered, chunk_line)                                            # add line to new, ordered df                
    names(scr_ordered) <- names(chunk_line)}                                                  # necessary to prevent error that column names of the data frames do not match; if column names don't match, line is not added to the df
  

  
  # rename columns; previous command (names(scr_ordered) <- names(chunk_line)) resulted in repetition of identical column names
  colnames(scr_ordered) <-c("Event.Nr.Tar","DDA.nSCR.Tar","DDA.Latency.Tar","DDA.AmpSum.Tar","DDA.AreaSum.Tar","DDA.Tonic.Tar","TTP.nSCR.Tar","TTP.Latency.Tar","TTP.AmpSum.Tar","Global.Mean.Tar","Global.MaxDeflection.Tar","Event.ID.Tar","Event.Name.Tar","order.Tar",                                     
                            "Event.Nr.TarResp","DDA.nSCR.TarResp","DDA.Latency.TarResp","DDA.AmpSum.TarResp","DDA.AreaSum.TarResp","DDA.Tonic.TarResp","TTP.nSCR.TarResp","TTP.Latency.TarResp","TTP.AmpSum.TarResp","Global.Mean.TarResp","Global.MaxDeflection.TarResp","Event.ID.TarResp","Event.Name.TarResp","order.TarResp",                                    
                            "Event.Nr.Word","DDA.nSCR.Word","DDA.Latency.Word","DDA.AmpSum.Word","DDA.AreaSum.Word","DDA.Tonic.Word","TTP.nSCR.Word","TTP.Latency.Word","TTP.AmpSum.Word","Global.Mean.Word","Global.MaxDeflection.Word","Event.ID.Word","Event.Name.Word","order.Word",                                         
                            "Event.Nr.WordResp","DDA.nSCR.WordResp","DDA.Latency.WordResp","DDA.AmpSum.WordResp","DDA.AreaSum.WordResp","DDA.Tonic.WordResp","TTP.nSCR.WordResp","TTP.Latency.WordResp","TTP.AmpSum.WordResp","Global.Mean.WordResp","Global.MaxDeflection.WordResp","Event.ID.WordResp","Event.Name.WordResp","order.WordResp")                       
  

  
  
  # add participant ID and trial ID
  scr_with_ID <- data.frame(Participant.ID=substr(subject,14,15),Trial.ID =1:nrow(scr_ordered),scr_ordered)   # create new df to add ID for each subject in first column; ID derived from letter 14 and 15 of file name
  
  
  
  # keep only columns of interest
  scr_final <- scr_with_ID[ , c("Participant.ID","Trial.ID","Event.ID.Tar","Event.ID.TarResp","DDA.AmpSum.TarResp","DDA.AreaSum.TarResp", "Event.ID.Word","Event.ID.WordResp" )]
                                                                    

  
  # display progress and abort if number of trials is not 516 -> in subject 5 line added between 1400 and 1401 (trigger 57 was missing)
  message("Preprocessing done for file: ",subject,appendLF=TRUE)
  
  if (nrow(scr_with_ID) != 516){
    stop("Incorrect number of trials in last subject. Check which trigger is missing and add manually in export file!")}
  
  
  
  
  ########################     save preprocessed   ##############################
  ########################        scr data         ##############################  
  
  # setwd("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/SCR_Export_Preprocessed")                              
  write.csv(scr_final, paste("P:/LuisaBalzus/1_PhD_Project/6_ModERN_Behavioral_Study/9_SCR_Export_Preprocessed", subject, sep="/"))
  
}


