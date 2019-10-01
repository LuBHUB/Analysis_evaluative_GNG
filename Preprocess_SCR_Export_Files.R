  ##### Implement necessary changes in SCR export files for Signle Trial Analysis of evaluative GNG task (Experiment on Affective Priming Effect 2018) ##### 
  ##### by Luisa Balzus
  ##### 11.11.2018 
  
  
  library(readxl)
  library(openxlsx)
  
  
  # clear environment
  rm(list=ls())


  ####################   load scr data   ####################
  
  # load SCR files (Ledalab output exported as .xls -> opening it requires readxl) 
  scrfiles <- list.files("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/8_SCR_Export", pattern = ".xls")    # path to folder containing the scr files (use of forward slashes instead of backward slashes is required); should contain ONLY scr files 
  setwd("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/8_SCR_Export")                              
  
  # loop reading SCR file
  for (subject in scrfiles){                                                                                         
    raw_scr <- read_excel(subject,sheet="CDA")
        
 
    
  ####################   clean and rearrange SCR data   ####################
  
  # only keep events (GNG and Word presentation and responses) and columns of interest 
  subset_scr <- subset(raw_scr, (Event.NID >= 21 & Event.NID <= 58) |(Event.NID >= 141 & Event.NID <= 149) |(Event.NID >= 241 & Event.NID <= 249))  
  subset_scr <- subset_scr[,c("Event.Nr","Event.NID","CDA.ISCR [muSxs]")] 
    
   
  # create new column that represents correct order of triggers within one trial
  subset_scr$order[subset_scr$Event.NID <=23]                            <-1
  subset_scr$order[subset_scr$Event.NID >=41 & subset_scr$Event.NID <=49]<-2
  subset_scr$order[subset_scr$Event.NID >=141]                           <-3
  subset_scr$order[subset_scr$Event.NID >=51 & subset_scr$Event.NID <=59]<-4
  
  
  # one chunk is one trial, containing 4 triggers: GNG stimulus, GNG response, word stimulus, word response
  number_chunks <- nrow(subset_scr)/4                                                                       
  scr_ordered   <- data.frame()  
  
  for (i in 0:(number_chunks-1)){
    chunk_start        <- 1+4*i
    chunk_end          <- 4+4*i
    chunk              <- subset_scr[chunk_start:chunk_end,]                                         # extract rows 1-4 of current chunk from scr data frame        
    chunk_ordered      <- chunk[order(chunk$order),]                                                 # order rows in chunk by the column "order" (this step may not be necessary anymore, because the downsampling is now done in Ledalab; before, it was done in BVA, so that the timing information was blurred and the CI triggers occurred after the word triggers in the exported file)
    chunk_line         <- c(chunk_ordered[1,],chunk_ordered[2,],chunk_ordered[3,],chunk_ordered[4,]) # write correctly ordered rows in single row
    scr_ordered        <- rbind(scr_ordered, chunk_line)                                             # add line to new, ordered df                
    names(scr_ordered) <- names(chunk_line)}                                                         # make sure that column names of the data frames do match; if column names don't match, line is not added to the df
  

  
  # rename columns; previous command (names(scr_ordered) <- names(chunk_line)) resulted in repetition of identical column names
  colnames(scr_ordered) <-c("Event.Nr.GNG","Event.ID.GNG","CDA.ISCR.GNG","order.GNG",                                     
                            "Event.Nr.GNG_Resp","Event.ID.GNG_Resp","CDA.ISCR.GNG_Resp","order.GNG_Resp",
                            "Event.Nr.Word","Event.ID.Word","CDA.ISCR.Word","order.Word",
                            "Event.Nr.Word_Resp","Event.ID.Word_Resp","CDA.ISCR.Word_Resp","order.Word_Resp")
                            
               
  # add participant ID and trial ID
  scr_with_ID <- data.frame(Participant.ID=substr(subject,14,15),Trial.ID =1:nrow(scr_ordered),scr_ordered)   # create new df to add ID for each subject in first column; ID derived as numeric from file name
  

  # IMPORTANT: assign stimulus-locked SCR values to correctly inhibited trials and trials with missing response (makes no sense to take response-locked SCR value here, because there is no response and the evaluation trigger is sent late)
  for (ind in 1:nrow(scr_with_ID)) {
    if (scr_with_ID$Event.ID.GNG_Resp[ind] == 45 | scr_with_ID$Event.ID.GNG_Resp[ind]  == 46 | scr_with_ID$Event.ID.GNG_Resp[ind]  == 47)
    {scr_with_ID$CDA.ISCR.GNG_Resp[ind] <- scr_with_ID$CDA.ISCR.GNG[ind]}
  }
  
  
  # keep only columns of interest
  scr_final <- scr_with_ID[ , c("Participant.ID","Trial.ID","Event.ID.GNG","CDA.ISCR.GNG","Event.ID.GNG_Resp","CDA.ISCR.GNG_Resp","Event.ID.Word","CDA.ISCR.Word","Event.ID.Word_Resp","CDA.ISCR.Word_Resp")]
                                                                    

  
  
  ####################   check file   ####################
  
  # display progress and abort if number of trials is not 516 -> in subject 5, one line manually added between 1400 and 1401 (trigger 57 was not send)
  message("Preprocessing done for file: ",subject,appendLF=TRUE)
  
  if (nrow(scr_with_ID) != 516){
    stop("Incorrect number of trials in last subject. Check which trigger is missing and add manually in export file!")}
  
  
  # for subject 25 and 29, SCR recording was interrupted; the missing trials were added manually to the files from the ledalab export (except for trigger number, these rows were left empty, because filling them with NA or NAN does nor work)

  
  
  ####################   save preprocessed scr data as xlsx   ####################

  write.xlsx(scr_final, paste0("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/9_SCR_Export_Preprocessed","/",subject,"x"))
}