Some notes on the files in the analysis folder.


PREPROCESSING SCR

The script Preprocess_SCR_Export_Files.R reads in the files that were exported from the Ledalab 
command (these are stored in folder 8_SCR_Export). The script preprocesses these files and stored 
the proprocessed files in the folder 9_SCR_Export_preprocessed.


PREPROCESSING GENERAL

The script SingleTiralAnalysis_SCR.R reads in raw behavioral data, preprocessed SCR data, and PEQ questionnaire data ("psychoEQExport_18.7.2019_12.25.sav") and does some data cleaning / preparation. Three new data tables are saved as an output: 
- a master table with the single trial data ("Data_Single_Trial_For_[...].rda" (= data4mixedmodels); this is read in for subsequent statistical analyses)
- a table with the questionnaire data ("Data_Trait_Variables_For_[...].rda" (= questionnaires); this is read in for subsequent statistical analyses)
- a table with some aggregated data ("Data_only_for_overview.rda" (= df4save); this is only used to have an overview, it is not used for subsequent analyses)


STATISTICAL ANALYSES

The script Script_(G)LMMs.R was used to build the final LMMs/GLMMs. 

All analyses are done in the Rmd html. The scripts read in the master table and the questionnaire table. To render the whole html, I have to run the script 
RUN_THIS_CODE_TO_RENDER_WEBSITE.R. This needs the following scripts to then render the website:
- _site.yml
- index.Rmd
- Demographic_Data.Rmd
- Word_Categorization_Data.Rmd	- This script creates the file "Data_Aggregated_For_[...].rda" (= df_wide). This is read in and used in the script Trait_Data_Correlations.Rmd. 
- SCR_Data.Rmd			- This script creates the file "Data_SCR_Aggregated_For_[...].rda" (= df_wide_scr). This is read in and used in the script Trait_Data_Correlations.Rmd. 
- GNG_Date.Rmd
- Trait_Data_Correlations.Rmd

The content in the folder website_Modern_behavioral is created when the html is rendered. 

Estimated_Means_[...].rda files (3 of these exist) save the output of the emmeans 
command which was used to backtransform the transformed model estimations to raw 
values (running this script took quite a while, so I saved the output and read it in 
in the Rmd files