setwd("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses")

if (file.exists("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses/website_ModERN_behavioral")){
  unlink("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses/website_ModERN_behavioral", recursive = TRUE)
}
  

rmarkdown::render_site()

