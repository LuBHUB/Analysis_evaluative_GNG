##### Statistical Anaylsis of Evaluative GNG Task in Behavioral Experiment on Affective Priming Effect 2018 ####
##### by Luisa Balzus
##### 09.11.2018


#library(dplyr)
#library(openxlsx) 
library(pastecs)
library(ggplot2)
#library(lmerTest)
#library(multcomp)
#library(emmeans)
#library(nlme)
#library(tidyr)
#library(psycho)
#library(car)
#library(foreign)
#library(psych)
#library(ez)
#library(jtools)




#setting wd
setwd("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_Behavioral/5_Analyses")

# clear environment
rm(list=ls())

# force R to not use exponential notation
options(scipen = 999)

#reading in datafile
load(file = ".rda")


datafiles <- list.files("P:/Luisa_Balzus/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses", pattern = ".rda")       
for (datafile in datafiles){  
  load(file = datafile)}




###################   Descriptive Statistics     ####################

descriptive_statistics <- stat.desc(df4save,basic=F)




###################   Bar Plots    ####################

# plot word classification rt all conditions
df4plotRT <- data.frame(
  response = c("false alarm","false alarm","fast hit","fast hit","correctly inhibited","correctly inhibited","slow hit", "slow hit"),
  valence = c("neg","pos","neg","pos","neg","pos","neg","pos"),
  conditions = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH","neg_after_CI","pos_after_CI","neg_after_SH","pos_after_SH"),
  mean = as.numeric(c(descriptive_statistics[2,2],descriptive_statistics[2,3],descriptive_statistics[2,4],descriptive_statistics[2,5],descriptive_statistics[2,6],descriptive_statistics[2,7],descriptive_statistics[2,8],descriptive_statistics[2,9])),
  se = as.numeric(c(descriptive_statistics[3,2],descriptive_statistics[3,3],descriptive_statistics[3,4],descriptive_statistics[3,5],descriptive_statistics[3,6],descriptive_statistics[3,7],descriptive_statistics[3,8],descriptive_statistics[3,9])))

df4plotRT$response <- factor(df4plotRT$response, levels = c("false alarm","fast hit","slow hit","correctly inhibited"))  # define order of conditions in plot

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
#ggsave("plot_rt.png", dpi=2000)

# plot word classification accuracy all conditions
df4plotACC <- data.frame(
  response = c("false alarm","false alarm","fast hit","fast hit","correctly inhibited","correctly inhibited","slow hit", "slow hit"),
  valence = c("neg","pos","neg","pos","neg","pos","neg","pos"),
  conditions = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH","neg_after_CI","pos_after_CI","neg_after_SH","pos_after_SH"),
  mean = as.numeric(c(descriptive_statistics[2,29],descriptive_statistics[2,30],descriptive_statistics[2,31],descriptive_statistics[2,32],descriptive_statistics[2,33],descriptive_statistics[2,34],descriptive_statistics[2,35],descriptive_statistics[2,36])),
  se = as.numeric(c(descriptive_statistics[3,29],descriptive_statistics[3,30],descriptive_statistics[3,31],descriptive_statistics[3,32],descriptive_statistics[3,33],descriptive_statistics[3,34],descriptive_statistics[3,35],descriptive_statistics[3,36])))

df4plotACC$response <- factor(df4plotACC$response, levels = c("false alarm","fast hit","slow hit","correctly inhibited")) # define order of conditions in plot

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
#ggsave("plot_acc.png", dpi=2000)

# plot SCR all conditions
df4plotSCR <- data.frame(
  response = c("false alarm","fast hit","correctly inhibited","slow hit"),
  valence = c("neg","pos","neg","pos"),
  conditions = c("mean_SCR_after_FA","mean_SCR_after_FH","mean_SCR_after_CI","mean_SCR_after_SH"),
  mean = as.numeric(c(descriptive_statistics[2,45],descriptive_statistics[2,46],descriptive_statistics[2,47],descriptive_statistics[2,48])),
  se = as.numeric(c(descriptive_statistics[3,45],descriptive_statistics[3,46],descriptive_statistics[3,47],descriptive_statistics[3,48])))

df4plotSCR$response <- factor(df4plotSCR$response, levels = c("false alarm","fast hit","slow hit","correctly inhibited"))  # define order of conditions in plot

ggplot(df4plotSCR,x = response,y = mean, aes(response, mean))+
  geom_bar(stat="identity", position=position_dodge(), fill = "mediumblue") +                                  # add bars, based on stats values; dodge to avoid stacked bars; change bar color
  geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.95), width=0.1) +    # add error bars
  ggtitle("Skin Conductance Response") +                                                                       # add title
  xlab("Response Type") + ylab("Logarithmized SCR (?S)") +                                                     # label axes
  theme(axis.title = element_text(size = 18, hjust = 0.5)) +                                                   # center title
  theme(axis.text=element_text(size=15)) +
  theme(legend.title=element_text(size=18)) +
  theme(legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  coord_cartesian(ylim = c(-0.25,0.7)) +
  scale_y_continuous(breaks=seq(0,40,0.2))                                                                     # set specific tic marks                                       
#ggsave("plot_scr.png", width = 6, height = 5.5, dpi=2000)






###################   Aggregate Data and Calculate ANOVAs   ####################


options(contrasts=c("contr.sum","contr.poly")) # to adhere to the sum-to-zero convention for effect weights, always do this before running ANOVAs in R. This matters sometimes (not always). If I don’t do it, the sum of squares calculations may not match what I get e.g. in SPSS


###################   ANOVA Word Classification   #################### 

# ANOVA requires several rows for each subject, each one row per factor: here I reorder data by reshaping existing data frame
df4anova <-   reshape(data = df4save[,c(1:9,29:36)], 
                      direction = "long",
                      varying = list(c("mean_neg_after_FA","mean_pos_after_FA","mean_neg_after_FH","mean_pos_after_FH","mean_neg_after_CI","mean_pos_after_CI" ,"mean_neg_after_SH","mean_pos_after_SH"),c("percent_correct_neg_after_FA", "percent_correct_pos_after_FA","percent_correct_neg_after_FH","percent_correct_pos_after_FH","percent_correct_neg_after_CI","percent_correct_pos_after_CI","percent_correct_neg_after_SH","percent_correct_pos_after_SH")),
                      v.names = c("rt","accuracy"),
                      idvar = "subject",
                      timevar = "condition",
                      times = c("neg_after_FA","pos_after_FA","neg_after_FH","pos_after_FH","neg_after_CI","pos_after_CI","neg_after_SH","pos_after_SH"))

row.names(df4anova) <- NULL
df4anova <- df4anova[sort.list(df4anova$subject),]                                      # sort df by subject

df4anova$response_type <- substr(df4anova$condition, 11, 12)                            # create columns needed as factors           
df4anova$word_valence  <- substr(df4anova$condition, 1, 3)

df4anova$subject       <- factor(df4anova$subject)                              
df4anova$response_type <- factor(df4anova$response_type, levels=c("FA","FH","CI","SH")) # create factor response_type and reorder factor levels! (this is important, because first level is used as reference factor in LMM)
df4anova$word_valence  <- factor(df4anova$word_valence)                                 # automatic factor order is alphabetical (1 = neg, 2 = pos)


# calculate ANOVA for rt -> gives same result as SPSS :)
anova_rt <- ezANOVA(data = df4anova, 
                    dv = .(rt), 
                    wid = .(subject), 
                    within = .(response_type, word_valence), 
                    detailed = TRUE)

# calculate ANOVA for accuracy -> gives same result as SPSS :)
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





###################   ANOVA GNG (for rt)      #################### 

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

df4anova_GNG_rt$subject   <- factor(df4anova_GNG_rt$subject)                                   # create factors
df4anova_GNG_rt$condition <- factor(df4anova_GNG_rt$condition)


# calculate ANOVA for rt GNG -> gives same result as SPSS :)
anova_GNG_rt <- ezANOVA(data = df4anova_GNG_rt, 
                        dv = .(rt), 
                        wid = .(subject), 
                        within = .(condition), 
                        detailed = TRUE)


# post hoc comparison: when calculating within-subject ANOVAs, no direct post hoc test is possible; I can only run paitwise t-tests with corresponding adjustment (e.g. Holm, which is better than Bonferroni)
anova_GNG_rt_posthoc_response_type <- pairwise.t.test(df4anova_GNG_rt$rt,df4anova_GNG_rt$condition,p.adjust.method="holm",paired=TRUE)



###################   ANOVA SCR and GNG Percent of Responses    #################### 


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

df4anova_GNG_scr$subject   <- factor(df4anova_GNG_scr$subject)                                    # create factors
df4anova_GNG_scr$condition <- factor(df4anova_GNG_scr$condition)




# calculate ANOVA for percent responses GNG
anova_GNG_percent <- ezANOVA(data = df4anova_GNG_scr, 
                             dv = .(percent), 
                             wid = .(subject), 
                             within = .(condition), 
                             detailed = TRUE)


# post hoc comparison: when calculating within-subject ANOVAs, no direct post hoc test is possible; I can only run paitwise t-tests with corresponding adjustment (e.g. Holm, which is better than Bonferroni)
anova_GNG_percent_posthoc_response_type <- pairwise.t.test(df4anova_GNG_scr$percent,df4anova_GNG_scr$condition,p.adjust.method="holm",paired=TRUE)




# calculate ANOVA for SCR
anova_GNG_SCR <- ezANOVA(data = df4anova_GNG_scr, 
                         dv = .(SCR), 
                         wid = .(subject), 
                         within = .(condition), 
                         detailed = TRUE)


# post hoc comparison: when calculating within-subject ANOVAs, no direct post hoc test is possible; I can only run paitwise t-tests with corresponding adjustment (e.g. Holm, which is better than Bonferroni)
anova_GNG_SCR_posthoc_response_type <- pairwise.t.test(df4anova_GNG_scr$SCR,df4anova_GNG_scr$condition,p.adjust.method="holm",paired=TRUE)




###################   Correlations with Traits   ####################

# possibly exclude subject 14 (great outlier, that causes most correlations/trends)!!!!


# test for normality
normality <- do.call(rbind, lapply(df4save[,-1], function(x) shapiro.test(x)[c("statistic", "p.value")]))


# create correlation matrices using my function "function_correlation_matrix_with_significance_levels.R"
source("P:/Luisa_Balzus/1_PhD_Project/5_ModERN_Vorstudie/10_Single_Trial_Analysis_EEG_SCR_GNG/function_correlation_matrix_with_significance_levels.R")


df4Pearson_Correlation <- df4save[,c("mean_priming_overall", "mean_priming_after_FA", "mean_priming_after_FH","mean_SCR_after_FA", "BIS_Total","BAS_Fun_Seeking", "BAS_Reward_Responsiveness", "FMPS_CMD", "FMPS_PST", "FMPS_PER/Total", "NEO_Neuroticism", "PANAS_Pos", "PSWQ_Total", "STAI_Trait", "TCI_Harm_Avoidance_Total", "TCI_Fear_of_uncertainty", "TCI_Shyness", "TCI_Fatigability")]   # only keep normally distributed variables 
method <- "pearson"
df4correlation <- df4Pearson_Correlation
correlations_pearson <- correlations_with_significance_level(df4Pearson_Correlation)


df4Kendall_Correlation <- df4save[,c("mean_priming_overall", "mean_priming_after_FA", "mean_priming_after_FH","mean_SCR_after_FA","BDI_Total", "BAS_Total", "BAS_Drive", "FMPS_PEC","FMPS_ORG","NEO_Conscientiousness", "OCI_Total", "OCI_Washing", "OCI_Checking", "OCI_Ordering", "OCI_Obsessions", "OCI_Hoarding", "OCI_Neutralising","PANAS_Neg","STAI_State","TCI_Anticipatory_worry")]   # only keep non-normally distributed variables (and those of major interest that are normally distributed to calculate correlation with those)
method <- "kendall"
df4correlation <- df4Kendall_Correlation
correlations_kendall <- correlations_with_significance_level(df4Kendall_Correlation)



# Template for plotting significant correlations
ggplot(df4save, aes (x= TCI_Fear_of_uncertainty, y = mean_SCR_after_FA)) +
  geom_point(color="mediumblue", size=3) +
  geom_smooth(method=lm, se=FALSE, color="limegreen", size = 1) +
  ggtitle("TCI_Fear_of_uncertainty and mean_SCR_after_FA") +
  labs(x= "TCI_Fear_of_uncertainty", y = "mean_SCR_after_FA") +
  theme(axis.title = element_text(size = 18, hjust = 0.5)) +
  theme(axis.text=element_text(size=15)) +
  theme(legend.title=element_text(size=18)) +
  theme(legend.text=element_text(size=15)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5))





###################   Linear Mixed Models for Mean RT to Words   #################### 

# prepare master_words for single trial analysis: exclude incorrect word responses, misses, wrong keys
data4mixedmodels_words  <- subset(data4mixedmodels,word_resp <= 52 & gng_resp <= 46 & outlier_words == FALSE)


# set contrasts
contrasts(data4mixedmodels_words$response_type) <- contr.sdif(4)
contrasts(data4mixedmodels_words$valence)       <- contr.sdif(2)


# calculate LMM for rt and plot
emm_options(pbkrtest.limit = 14000)
emm_options(lmerTest.limit = 14000) # due to warning D.f. calculations have been disabled because the number of observations exceeds 3000.To enable adjustments, set emm_options(pbkrtest.limit = 5866)

LMM_rt <- lmer(word_rt ~ response_type * valence + (1|subjectID), data=data4mixedmodels_words, REML = FALSE)
summary(LMM_rt)
anova(LMM_rt)                                                                                     # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
results_LMM_rt <- analyze(LMM_rt)                                                                 # print results
print(results_LMM_rt)
results_LMM_rt <- get_contrasts(LMM_rt, "response_type * valence")                                # provide the model and the factors to contrast; add ,adjust="none" to turn off automatic p value correction after Tucky
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




###################   Linear Mixed Models for Mean RT to GNG Stimulus   ####################

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
results_LMM_GNG_rt <- get_contrasts(LMM_GNG_rt, "response_type")                             # provide the model and the factors to contrast; by default, get_contrasts uses the Tukey method for p value adjustment; add ,adjust="none" to turn off automatic p value correction after Tucky
print(results_LMM_GNG_rt$means)                                                              # investigate means
print(results_LMM_GNG_rt$contrasts)                                                          # print contrasts



ggplot(results_LMM_GNG_rt$means, aes(x=response_type, y=Mean, color=response_type)) +        # plot
  geom_line(position = position_dodge(.3)) +
  geom_pointrange(aes(ymin=CI_lower, ymax=CI_higher),
                  position = position_dodge(.3)) +
  ylab("Response Time in ms") +
  xlab("Response Type") +
  theme_bw()


# post hoc tests for resolving main effect
summary(glht(LMM_GNG_rt, linfct=mcp(response_type = "Tukey")), test = adjusted("holm"))







###################   Linear Mixed Models for SCR   ####################

# exclude CI and SH (SH excluded to only have 2 levels and make interpretation easier), exclude trials with outlier or error in word or miss/wrong key in GNG, exclude trials with invalid GNG rt, exclude trials followed or preceded by false alarm or wrong key in GNG or by incorrect word categorization or wrong key in word categorization
data4mixedmodels_scr <- subset(data4mixedmodels,word_resp <= 52 & (gng_resp == 41 | gng_resp == 43 | gng_resp == 44)  & outlier_words == FALSE & gng_invalid_rt == FALSE & followed_or_preceded_by_FA_or_wrong_key == FALSE)
data4mixedmodels_scr$response_type <- droplevels(data4mixedmodels_scr$response_type) # drop unused levels in response type, because conditions (CI and SH) were excluded


# set contrasts for fixed effects
contrasts(data4mixedmodels_scr$response_type) <- contr.sdif(2)
contrasts(data4mixedmodels_scr$valence)       <- contr.sdif(2)


# calculate LMM for SCR and plot 
emm_options(pbkrtest.limit = 4000)
emm_options(lmerTest.limit = 4000) # due to warning D.f. calculations have been disabled because the number of observations exceeds 3000.To enable adjustments, set emm_options(pbkrtest.limit = 5866)


LMM_SCR <- lmer(scr_amplitude_log ~ response_type*valence*facilitation_score + (1 + response_type | subjectID) + (1|word), data=data4mixedmodels_scr, REML = FALSE)
summary(LMM_SCR)
anova(LMM_SCR)                                                                          # get ANOVA Output from LMM (results differs a bit from that obtained by using aov/ezANOVA; https://stackoverflow.com/questions/20959054/why-is-there-a-dramatic-difference-between-aov-and-lmer)
results_LMM_SCR <- analyze(LMM_SCR)                                                     # print results
print(results_LMM_SCR)                                                                  # reports incorrect beta estimates
results_LMM_SCR <- get_contrasts(LMM_SCR, "response_type * valence")                    # provide the model and the factors to contrast;; add ,adjust="none" to turn off automatic p value correction after Tucky
print(results_LMM_SCR$contrasts)                                                        # print contrasts
print(results_LMM_SCR$means)                                                            # investigate means

stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}

stdCoef.merMod(LMM_SCR)

# significant main effect of response type: FA lead to higher SCR then FH
# significant interaction response_type*facilitation: driven by neg words after FA, high priming magnitude (fast response to neg words) related to low SCR -> see plot of 3-way interaction
# significant interaction response_type*valence*facilitation_score: only for neg words after FA, high priming magnitude (fast response to neg words) related to low SCR
# (hypothesis was that high facilitation for neg words after incorrect responses predicts higher SCR)
interact_plot(LMM_SCR, pred = "facilitation_score", modx = "response_type", mod2 = "valence", interval = TRUE,int.width = 0.8,
              x.label = "Magnitude of Priming (ms)", y.label = "Logarithmized SCR (?S)", modx.labels = c("false alarm","fast hit"), mod2.labels = c("Word Valence = neg","Word Valence = pos"),
              main.title = "SCR Depending on Priming Effect, Response Type, and Word Valence", legend.main = "Response Type")

ggsave("plot_lmm.png", width = 6, height = 3.5, dpi=2000)