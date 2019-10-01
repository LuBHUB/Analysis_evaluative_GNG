setwd("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses")

if (file.exists("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses/website_ModERN_behavioral")){
  unlink("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses/website_ModERN_behavioral", recursive = TRUE)
}
  

rmarkdown::render_site()

