##### (G)LMMs for Evaluative GNG Task in Behavioral Experiment on Affective Priming Effect 2018 ####
##### by Luisa Balzus
##### 22.08.2019

library(plyr)                    # for ddply()
library(dplyr)                   # for mutating
library(pastecs)                 # for descriptive stats
library(ggplot2)                 # for plots
library(lme4)                    # for (G)LMMs
library(lmerTest)                # for p values for tests for fixed effects (Satterthwaite's method for approximating degrees of freedom for the t and F tests)
library(MASS)                    # for box cox and contrast definition
library(car)                     # for qqp function
library("GeneralizedHyperbolic") # for fitting inverse gaussian distribution
library(e1071)                   # for functions skewness and kurtosis
library(fitdistrplus)            # for finding distribution that fits SCR data
library(sjPlot)                  # for tab_model

# clear environment
rm(list=ls())

# force R to not use exponential notation
options(scipen = 999)

# reading in datafile
datafiles <- list.files("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses", pattern = ".rda") 
for (datafile in datafiles){  
  # appending full path to filename is necessary to open files in Rmd
  filename <- paste0("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/5_Analyses/",datafile) 
  load(file = filename)}

# prepare labels for tables
pl <- c(
  "(Intercept)" = "Intercept",
  "valence2-1" = "Positive - Negative",
  "response_type2-1" = "FH - SH", 
  "response_type3-2" = "FA - FH",
  "response_type4-3" = "CI - FA",
  "response_type2-1:valence2-1" = "FH - SH x Positive - Negative", 
  "response_type3-2:valence2-1" = "FA - FH x Positive - Negative",
  "response_type4-3:valence2-1" = "CI - FA x Positive - Negative",
  "response_typeSH:valence2-1" = "SH: Positive - Negative",
  "response_typeFH:valence2-1" = "FH: Positive - Negative",
  "response_typeFA:valence2-1" = "FA: Positive - Negative",
  "response_typeCI:valence2-1" = "CI: Positive - Negative")


# define function to test convergence
didLmerConverge = function(lmerModel){
  relativeMaxGradient=signif(max(abs(with(lmerModel@optinfo$derivs,solve(Hessian,gradient)))),3)
  if (relativeMaxGradient < 0.001) {
    cat(sprintf("\tThe relative maximum gradient of %s is less than our 0.001 criterion.\n\tYou can safely ignore any warnings about a claimed convergence failure.\n\n", relativeMaxGradient))
  }
  else {
    cat(sprintf("The relative maximum gradient of %s exceeds our 0.001 criterion.\nThis looks like a real convergence failure; maybe try simplifying your model?\n\n", relativeMaxGradient))
  }
}


# lme4 version 18 and newer seem to cause problems; optimization default was changed from bobyqa to nloptwrap; this led to non convergence of most of my models; I tried to install version 17 but failed. So I am now using version 21 but with the optimizer bobyqa

# model comparison in this script might not be correct and is commented out (I did not really know how to interpret output (LRT can go also in different direction then AIC and BIC) and I chose .05 as cutoff (this is not recommended, see Matuschek 2017)). Schad recommends to only do rePCA in order to exclude zero variance components. Leaving random effects in although they do not improve model fit helps to better control the type I error and allows generalization over these effects.

###################   Preparation for Priming    ####################    


###### STEP 1: PREPARE DATA 

# exclude word and GNG responses with misses or wrong keys, outliers in word RT or GNG RT; incorrect responses for RT LMMs are excluded when the LMM is specified
data4mixedmodels_words <- data4mixedmodels[data4mixedmodels$outlier_words == FALSE & data4mixedmodels$gng_resp <= 46 & data4mixedmodels$gng_invalid_rt == FALSE & data4mixedmodels$word_resp <= 54,]
# 14130 of 15480 trials left

# create variable error rate
df_error_rate <-   ddply(data4mixedmodels_words, 
                             .(subjectID), 
                             summarise, 
                             error_rate = 100-sum(word_accuracy)/length(subjectID)*100)
data4mixedmodels_words <- left_join(data4mixedmodels_words,df_error_rate, by = "subjectID")


# make categorical variables factors
data4mixedmodels_words$response_type <- factor(data4mixedmodels_words$response_type, levels=c("SH","FH","FA","CI"))
data4mixedmodels_words$valence       <- factor(data4mixedmodels_words$valence)
data4mixedmodels_words$condition     <- factor(data4mixedmodels_words$condition, levels=c("neg_after_SH","pos_after_SH","neg_after_FH","pos_after_FH","neg_after_FA","pos_after_FA","neg_after_CI","pos_after_CI"))
data4mixedmodels_words$subjectID     <- factor(data4mixedmodels_words$subjectID)
data4mixedmodels_words$word          <- factor(data4mixedmodels_words$word)
data4mixedmodels_words$word_accuracy <- factor(data4mixedmodels_words$word_accuracy)




###### STEP 2: CHECk DISTRIBUTION 

# transformed rt -> normally distributed, also in the different conditions
plot(density(data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,]$word_rt_inverse), main = "Histogram Word RT Inverse")
ggplot(data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], aes(x = word_rt_inverse)) + geom_density() + facet_wrap(response_type ~ valence) + ggtitle("Histogram Word RT Inverse in Conditions")
qqp(data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,]$word_rt_inverse, "norm", main = "Q-Q Plot Word RT Inverse", ylab = "sample quantiles")

# raw rt -> not normally distributed
plot(density(data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,]$word_rt), main = "Histogram Word RT Raw")
ggplot(data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], aes(x = word_rt)) + geom_density() + facet_wrap(response_type ~ valence) + ggtitle("Histogram Word RT Raw in Conditions")
qqp(data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,]$word_rt, "norm", main = "Q-Q Plot Word RT Raw", ylab = "sample quantiles")

# raw rt -> fit gamma
gamma <- fitdistr(data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,]$word_rt, "gamma")
qqp(data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,]$word_rt, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]], ylab = "sample quantiles",main = "Q-Q Plot Word RT gamma")

# raw rt -> fit inverse gaussian -> inverse gaussian distribution fits quite well (for RTs, Lo and Anderson 2015 recommend gamma or inverse gaussian, Conny used Gamma)
inv_gauss <- nigFit(data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,]$word_rt, plots = FALSE, printOut = FALSE) 

par(mfrow = c(1,2))  # use par function in base graphics approach to arrange plots side by side here; grid.arrange only works for ggplot grahpics
plot(inv_gauss, which = 1, plotTitles = "Histogram of Word RT inv. gaussian")
plot(inv_gauss, which = 3, plotTitles = "Q-Q Plot of Word RT inv. gaussian")
par(mfrow = c(1, 1)) # reset "par" parameter

# I would choose inverse gaussian distribution. Lo & Andrews (2015) recommend gamma or inverse gaussian distribution and identity link. For my data, inverse gaussian distribution fits better.
# Lo & Anderson (2015) recommend gamma or inverse gaussian distribution for RTs. For the word RTs, inverse gaussian distribution fits very well and better than gamma distribution. Hence, I would choose inverse gaussian distribution and the identity link funktion (recommended by Lo & Andrews, 2015) for my GLMM. 


###### STEP 3: DEFINE CONTRASTS 

# I use sliding contrast (reasons for decision see Evernote entry 2019_08_06)


####### uncomment, if I want to customize contrasts, here example with sliding contrast
# # I checked whether the predefined way gives the same result. It does.
#
# # generalized inversed function
# ginv2 <- function(x) # define a function to make the output nicer, but is equal to ginv
#   MASS::fractions(provideDimnames(MASS::ginv(x),base=dimnames(x)[2:1]))
# 
# # specify hypothesis matrix response type  
# t(contr_response_type_inv <- t(cbind(FH_minus_SH=c(XSH = -1,XFH = 1, XFA= 0, XCI = 0),
#                                      FA_minus_FH=c(XSH =  0,XFH =-1, XFA= 1, XCI = 0),
#                                      CI_minus_FA=c(XSH =  0,XFH = 0, XFA=-1, XCI = 1))))
# contr_response_type_inv
# contr_response_type <- ginv2(contr_response_type_inv) # get contrast matrix
# contr_response_type # contrast matix
# contrasts(data4mixedmodels_words$response_type) <- contr_response_type # set contr_response_type as contrast matrix
# 
# # specify hypothesis matrix word valence 
# t(contr_word_valence_inv <- t(cbind(pos_minus_neg=c(Xneg=-1,Xpos=1))))
# contr_word_valence_inv
# contr_word_valence <- ginv2(contr_word_valence_inv) # get contrast matrix
# contr_word_valence # contrast matix
# contrasts(data4mixedmodels_words$valence) <- contr_word_valence # set contr_word_valence as contrast matrix
# 
# # check hypotheses (show hypothesis matrix again)
# ginv2(contr_response_type)
# ginv2(contr_word_valence)
# 
# str(data4mixedmodels_words)
# 
# # add contrast as numerical covariate via model matrix (in order to enter contrasts individually in model and hence to use ||)
# mm_c    <- model.matrix( ~ response_type*valence, data4mixedmodels_words) # stick in model matrix 8 columns
# glimpse(mm_c)                                                             # requires dplyr
# glimpse(data4mixedmodels_words)                                           # at the moment has 25 columns
# 
# # attach to dataframe 
# data4mixedmodels_words[,(ncol(data4mixedmodels_words)+1):(ncol(data4mixedmodels_words)+8)] <- mm_c
# names(data4mixedmodels_words)[(ncol(data4mixedmodels_words)-7):ncol(data4mixedmodels_words)] <- c("Grand Mean", "FH_minus_SH","FA_minus_FH","CI_minus_FA","pos_minus_neg", "FH_minus_SH:pos_minus_neg", "FA_minus_FH:pos_minus_neg", "CI_minus_FA:pos_minus_neg")
# glimpse(data4mixedmodels_words)



####### uncomment, if I want to customize contrasts, here example with custom contrast for testing directly nested interaction 
# 
# # I checked whether the custom contrast investigating directly the nested interaction  gives same result as nested model derived from the full interaction model with repeated contrasts. It does.
# 
# # generalized inversed function
# ginv2 <- function(x) # define a function to make the output nicer, but is equal to ginv
#   MASS::fractions(provideDimnames(MASS::ginv(x),base=dimnames(x)[2:1]))
# #
# # specify hypothesis matrix
# t(contr_inv <-   rbind(FH_minus_SH=c(F1 =-0.5, F2 =-0.5, F3= 0.5, F4 = 0.5, F5 = 0,   F6 = 0,   F7 = 0,   F8 = 0),
#                        FA_minus_FH=c(F1 = 0,   F2 = 0,   F3=-0.5, F4 =-0.5, F5 = 0.5, F6 = 0.5, F7 = 0,   F8 = 0),
#                        CI_minus_FA=c(F1 = 0,   F2 = 0,   F3= 0,   F4 = 0,   F5 =-0.5, F6 =-0.5, F7 = 0.5, F8 = 0.5),
#                 SH_pos_minus_neg = c(F1 =-1,   F2 = 1,   F3= 0,   F4 = 0,   F5 = 0,   F6 = 0,   F7 = 0,   F8 = 0),
#                 FH_pos_minus_neg = c(F1 = 0,   F2 = 0,   F3=-1,   F4 = 1,   F5 = 0,   F6 = 0,   F7 = 0,   F8 = 0),
#                 FA_pos_minus_neg = c(F1 = 0,   F2 = 0,   F3= 0,   F4 = 0,   F5 =-1,   F6 = 1,   F7 = 0,   F8 = 0),
#                 CI_pos_minus_neg = c(F1 = 0,   F2 = 0,   F3= 0,   F4 = 0,   F5 = 0,   F6 = 0,   F7 =-1,   F8 = 1)))
# 
# contr_inv
# contr <- ginv2(contr_inv) # get contrast matrix
# contr # contrast matix
# contrasts(data4mixedmodels_words$condition) <- contr # set contr as contrast matrix
# 
# # check hypotheses (show hypothesis matrix again)
# ginv2(contr)
# 
# str(data4mixedmodels_words)
# 
# # add contrast as numerical covariate via model matrix (in order to enter contrasts individually in model and hence to use ||)
# mm_c    <- model.matrix( ~ condition, data4mixedmodels_words) # stick in model matrix 8 columns
# glimpse(mm_c)                                                 # requires dplyr
# glimpse(data4mixedmodels_words)                               # at the moment has 25 columns
# 
# # attach to dataframe
# data4mixedmodels_words[,(ncol(data4mixedmodels_words)+1):(ncol(data4mixedmodels_words)+8)] <- mm_c
# names(data4mixedmodels_words)[(ncol(data4mixedmodels_words)-7):ncol(data4mixedmodels_words)] <- c("Grand Mean", "FH_minus_SH","FA_minus_FH","CI_minus_FA","pos_minus_neg", "FH_minus_SH:pos_minus_neg", "FA_minus_FH:pos_minus_neg", "CI_minus_FA:pos_minus_neg")
# glimpse(data4mixedmodels_words)
# 
# # example LMM
# example <- lmer(word_rt ~ FH_minus_SH + FA_minus_FH + CI_minus_FA + SH_pos_minus_neg + FH_pos_minus_neg + FA_pos_minus_neg + CI_pos_minus_neg + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA | subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE)
# summary(example)
# example2 <- lmer(word_rt ~ condition + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA | subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE)





####### uncomment, if I want to use predefined contrasts, here example with repeated contrast
contrasts(data4mixedmodels_words$response_type) <- contr.sdif(4)
contrasts(data4mixedmodels_words$valence)       <- contr.sdif(2)

# add contrast as numerical covariate via model matrix (in order to enter contrasts individually in model and hence to use ||)
mm_c    <- model.matrix( ~ response_type*valence, data4mixedmodels_words) # stick in model matrix 8 columns
glimpse(mm_c)                                                             # requires dplyr
glimpse(data4mixedmodels_words)                                           # at the moment has 25 columns

# attach to dataframe 
data4mixedmodels_words[,(ncol(data4mixedmodels_words)+1):(ncol(data4mixedmodels_words)+8)] <- mm_c
names(data4mixedmodels_words)[(ncol(data4mixedmodels_words)-7):ncol(data4mixedmodels_words)] <- c("Grand Mean", "FH_minus_SH","FA_minus_FH","CI_minus_FA","pos_minus_neg", "FH_minus_SH:pos_minus_neg", "FA_minus_FH:pos_minus_neg", "CI_minus_FA:pos_minus_neg")
glimpse(data4mixedmodels_words)


###### STEP 4: SCALE CONTINUOUS PREDICTORS / COVARIATES 

# goal: all PREDICTORS are on similar scale (not independent variable! my constrasts for the categorical predictors are between -1 and 1 -> covariate shall approximately be there as well)
# achieve this by: 1) scaling OR 2) centering and dividing by constant (important is, that transformation is linear, otherwise it may be problematic for interpretation; e.g. 1/RT is non-linear)


# alternative 1: center and divide by standard deviation (= scaling) - my prefered alternative, because it is more commonly used than centering and dividing by a constant and convergence problems for the models I tested are for both approaches the same
data4mixedmodels_words$gng_rt_scaled                                                <- data4mixedmodels_words$gng_rt
data4mixedmodels_words$gng_rt_scaled[data4mixedmodels_words$gng_rt_scaled == 0]     <- NaN # in CI, RT is 0; I reassign NaN to these trials, because they shall not be included in scaling (this would alter mean and sd)
data4mixedmodels_words$gng_rt_scaled[!is.nan(data4mixedmodels_words$gng_rt_scaled)] <- scale(data4mixedmodels_words$gng_rt_scaled[!is.nan(data4mixedmodels_words$gng_rt_scaled)], scale=TRUE, center=TRUE)
data4mixedmodels_words$gng_rt_scaled[is.nan(data4mixedmodels_words$gng_rt_scaled)]  <- 0 # reassign 0, bec. otherwise CI main effects and interactions are not estimated, and other factors are not estimated according to their defined contrast as well, bec. CI is included in contrasts
# with this scaling, the covariate comes quite close to other predictors (constrasts for categorical predictors are between -1 and 1, scaled RT is between -2.5 and 4)


# alternative 2: centering and divide by constant (e.g. 300, wich is approximate mean of GNG RT) to achieve values in range between -1 and 1
# data4mixedmodels_words$gng_rt_scaled2                                                 <- data4mixedmodels_words$gng_rt
# data4mixedmodels_words$gng_rt_scaled2[data4mixedmodels_words$gng_rt_scaled2 == 0]     <- NaN
# data4mixedmodels_words$gng_rt_scaled2[!is.nan(data4mixedmodels_words$gng_rt_scaled2)] <- scale(data4mixedmodels_words$gng_rt_scaled2[!is.nan(data4mixedmodels_words$gng_rt_scaled2)], scale=FALSE, center=TRUE)
# data4mixedmodels_words$gng_rt_scaled2[!is.nan(data4mixedmodels_words$gng_rt_scaled2)] <- data4mixedmodels_words$gng_rt_scaled2[!is.nan(data4mixedmodels_words$gng_rt_scaled2)] / 300
# data4mixedmodels_words$gng_rt_scaled2[is.nan(data4mixedmodels_words$gng_rt_scaled2)]  <- 0 # otherwise CI main effects and interactions are not estimated

 
 # Things I find problematic when including gng_rt as covariate
 #  - when I want to keep my contrasts including CI, there is no gng_rt for these trials
 #      - assigning NaN to these trials leads to dropping of estimation for all CI main effects and interactions
 #      - if I replace NaNs with 0, I get normal estimates (also for CI), but I am not sure in how far this messes up interpretability (now there is a real value for an unreal event)

 # Decision: I favor model without gng_rt as covariate, because reliably testing priming (also in CI) is more important to me than gng_rt covariate
 #  - also, the results don't change when gng_rt is included as covariate
 

# error rate is not continuous, and thus there is no need to scale it (additionally, a magnitude of zero is meaningfull here)


###################   Linear Mixed Models for Priming RT inverse    #################### 

###### STEP 1 to 4: RUN CODE ABOVE (PREPARATION FOR PRIMING)


###### STEP 5: SPECIFY THE MAXIMAL MODEL 

LMM_rt_max <- lmer(word_rt_inverse ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg | subjectID) + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_max)                         # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
summary(rePCA(LMM_rt_max))  
print(VarCorr(LMM_rt_max),comp="Variance")  
# I did not check whether according to didLmerConverge the max model does converge, because as seen below, I have to reduce it for some zero variance components and removing the correlations would have been the first step anyway



###### STEP 6: SIMPLIFICATION OF THE RANDOM EFFECTS STRUCTURE 

# start with zero-correlation parameter model by using ||
LMM_rt_red1 <- lmer(word_rt_inverse ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_red1)                        # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
summary(rePCA(LMM_rt_red1))                 # for word there is one item that does not explain variance
print(VarCorr(LMM_rt_red1),comp="Variance") # it is CI_minus_FA

# the model still converges -> it would also be ok to use the maximal model, as for the full model vs the reduced model all estimates are quite identical effect remain the same 
didLmerConverge(LMM_rt_red1) # 	The relative maximum gradient of 0.00000196 is less than our 0.001 criterion. You can safely ignore any warnings about a claimed convergence failure.

# remove CI_minus_FA from words
LMM_rt_red2 <- lmer(word_rt_inverse ~ FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_red2)        # model does converge, no singular fit -> model is identified
summary(rePCA(LMM_rt_red2)) # all items explain variance
print(VarCorr(LMM_rt_red2),comp="Variance") 


### now that we have an identified model, test non-significant variance components using likelihood ratio tests (starting with component explaning least variance) - Schad recommends drop 1 strategy, but he would rather only do rePCA
# 
# # test whether random slope reaction_type for word improves fit (explains least variance)
#LMM_rt_red2 <- lmer(word_rt_inverse ~ FH_minus_SH + FA_minus_F(1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
#LMM_rt_red3 <- lmer(word_rt_inverse ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_red3)           # model does converge, no singular fit -> model is identified
# anova(LMM_rt_red2,LMM_rt_red3) # no sign. difference, random slope reaction_type for word does not improve fit -> model 3 wins
# 
# # test whether random slope valence for subject improves fit (explains less variance than other random slopes for subjectID)
# LMM_rt_red4 <- lmer(word_rt_inverse ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_red4)           # model does converge, no singular fit -> model is identified
# anova(LMM_rt_red3,LMM_rt_red4) # random slope valence for subject improves fit -> model 3 wins
# 
# # test whether random slope response_type for subject improves fit (slope that explains less variance than random slope of interaction for subjectID)
# LMM_rt_red5 <- lmer(word_rt_inverse ~ response_type*valence + (1 + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_red5)           # model does converge, no singular fit -> model is identified
# anova(LMM_rt_red3,LMM_rt_red5) # random slope response_type for subject improves fit -> model 3 wins
# 
# # test whether random slope of interaction response type*valence for subject improves fit
# LMM_rt_red6 <- lmer(word_rt_inverse ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_red6)           # model does converge, no singular fit -> model is identified
# anova(LMM_rt_red3,LMM_rt_red6) # random slope of interaction response type*valence for subject improves fit -> model 3 wins
# 
# # test whether random intercept for word improves fit
# LMM_rt_red7 <- lmer(word_rt_inverse ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_red7)           # model does converge, no singular fit -> model is identified
# anova(LMM_rt_red3,LMM_rt_red7) # random intercept for word improves fit -> model 3 wins
# 
# # test whether random intercept for subject improves fit
# LMM_rt_red8 <- lmer(word_rt_inverse ~ response_type*valence + (0 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_red8)           # model does converge, no singular fit -> model is identified
# anova(LMM_rt_red3,LMM_rt_red8) # random intercept for subject improves fit -> model 3 wins
# 
# # test whether reintroduction of correlation improves fit
# LMM_rt_red9 <- lmer(word_rt_inverse ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg | subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_red9)           # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
# 
# 
# 
# ###### STEP 7: TEST OMNIBUS SIGNIFICANCE OF FIXED EFFECTS (based on model reduced by rePCA and LRT)
# 
# # test whether fixed effect of interaction response_type*valence improves fit
# LMM_rt_red10 <- lmer(word_rt_inverse ~ response_type + valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_red10)           # model does converge -> identified model
# anova(LMM_rt_red3,LMM_rt_red10) # fixed effect of interaction response_type*valence improves fit  
# 
# # test whether fixed effect of valence improves fit
# LMM_rt_red11 <- lmer(word_rt_inverse ~ response_type + response_type:valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_red11)           # model does converge -> identified model
# anova(LMM_rt_red3,LMM_rt_red11) # no sign. difference; fixed effect of valence does not improve fit 
# 
# # test whether fixed effect of response_type improves fit
# LMM_rt_red12 <- lmer(word_rt_inverse ~ valence + response_type:valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_red12)           # model does converge -> identified model
# anova(LMM_rt_red3,LMM_rt_red12) # sign. difference, but AIC identical; fixed effect of response_type does not improve fit??



###### STEP 8: PRESENT FINAL MODEL BY USING REML ESTIMATIONS 

LMM_rt_final <- lmer(word_rt_inverse ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = TRUE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_final)        # model does converge, no singular fit -> model is identified




###### STEP 9: BREAK DOWN INTERAKTIONS BY USING NESTED MODEL 

# valence within response condition
LMM_rt_nested <- lmer(word_rt_inverse ~ response_type/valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = TRUE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_nested)        # model does converge, no singular fit -> model is identified


# to compute pairwise comparisons to understand the interaction (but see alternative approach with nested models above):
difflsmeans(LMM_rt_final, test.effs = "response_type:valence")
# note that the results provided by difflsmeans and nested models are identical


###### STEP 10: CHECK ASSUMPTIONS 

# LINEARITY & HOMOGENEITY OF VARIANCE
plot(fitted(LMM_rt_final), residuals(LMM_rt_final), xlab = "Fitted Values", ylab = "Residuals")  # plot residuals against the fitted values
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(LMM_rt_final), residuals(LMM_rt_final)))
# scatter looks like a random pattern ->  residuals seem to meet model assumption of normality
# solid line covers the dashed line -> best-fit line fits well
# if the scatter is not random that means there's some other random or fixed effect that explains variation in your data that you're missing
# alternative: plot(LMM_rt_final)


# NORMALITY OF RESIDUALS
qqPlot(resid(LMM_rt_final)) 
qqnorm(resid(LMM_rt_final))
qqline(resid(LMM_rt_final)) # a bit off at the extremes, but that's often the case; again doesn't look too bad


###################   Test Covariates ####
###### Include GNG_RT as covariate??

# For different reasons I do not really want to include GNG RT as covariate (see notes in STEP 3)
# But I want to check whether inclusion of this covariate may change the results
# so  far I did not include it as a random slope; I just include the Covariate as a fixed effect in the final LMM -> not including it as random slope is fine according to Schad, as I do not want to generalize over it / this effect is of no interest 
# including gng rt as covariate improves the fit. However, it causes problems with regard to CI factor level. The results of all main effects and interaction are the same, independent of whether the covariate is included!!! 
# also, the 3-way interaction gng_rt_scaled:valence:response type is not signifcant, suggesting that gng_rt does not affect the priming effect


# warning: fixed-effect model matrix is rank deficient so dropping 2 columns / coefficients -> because for CI there is only one RT
# read here: https://stackoverflow.com/questions/37090722/lme4lmer-reports-fixed-effect-model-matrix-is-rank-deficient-do-i-need-a-fi  AND https://stats.stackexchange.com/questions/242540/r-lm-could-anyone-give-me-an-example-of-the-misleading-case-on-prediction-fr/356083#356083
# Rank-deficiency is not a problem for valid model estimation and comparison
# I also removed interaction CI and gng_rt_inverse from fixed and random effects, as this causes rank deficiency - does not change anything, because it is not estimated anyway - I can use rank-deficient model

LMM_rt_cov <- lmer(word_rt_inverse ~ response_type*valence*gng_rt_scaled + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov)
anova(LMM_rt_final,LMM_rt_cov)         # fit is better with covariate (but models are not nested!!!)

# test whether 3-way interaction improves fit
LMM_rt_cov_red1 <- lmer(word_rt_inverse ~ response_type*valence + gng_rt_scaled + gng_rt_scaled:response_type + gng_rt_scaled:valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_red1)               # model does converge -> identified model
anova(LMM_rt_cov,LMM_rt_cov_red1)      # 3-way interaction gng_rt_scaled:valence:response type does not improve fit (no sign. diff.) -> priming not modulated by gng_rt!!!! -> model LMM_rt_cov_red1 wins

# test whether interaction gng_rt_scaled:valence improves fit
LMM_rt_cov_red2 <- lmer(word_rt_inverse ~ response_type*valence + gng_rt_scaled + gng_rt_scaled:response_type + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_red2)               # model does converge -> identified model
anova(LMM_rt_cov_red1,LMM_rt_cov_red2) # interaction gng_rt_scaled:valence does not improve fit (no sign. diff.) -> model LMM_rt_cov_red2 wins

# test whether interaction gng_rt_scaled:response_type improves fit
LMM_rt_cov_red3 <- lmer(word_rt_inverse ~ response_type*valence + gng_rt_scaled + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_red3)               # model does converge -> identified model
anova(LMM_rt_cov_red2,LMM_rt_cov_red3) # interaction gng_rt_scaled:response_type does not improve fit -> model LMM_rt_cov_red3 wins

# final model
LMM_rt_cov_final <- lmer(word_rt_inverse ~ response_type*valence + gng_rt_scaled + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = TRUE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_final)              # model does converge -> identified model

# nested model
LMM_rt_cov_nested <- lmer(word_rt_inverse ~ response_type/valence + gng_rt_scaled + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = TRUE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_nested)              # model does converge -> identified model
anova(LMM_rt_nested, LMM_rt_cov_nested) # fit is better with covariate



###### Include Error Rate as covariate??
# Center Covariate
data4mixedmodels_words$error_rate <- scale(data4mixedmodels_words$error_rate, center = TRUE, scale = FALSE)

# Full model: response type:valence remains significant; there is no significant 3 way interaction response_type4-3:valence2-1:error_rate
LMM_rt_cov_error_rate <- lmer(word_rt_inverse ~ response_type*valence*error_rate + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_error_rate)

# Nested model: priming after FA, FH, SH remains (for no response type there is a significant interaction valence x error rate)
LMM_rt_cov_error_rate_nested <- lmer(word_rt_inverse ~ (response_type/valence)*error_rate + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_error_rate_nested)

# 3-way interaction error_rate:valence:response type does not improve fit (no sign. diff.) -> priming not modulated by error rate!!!!
LMM_rt_cov_error_rate_red1 <- lmer(word_rt_inverse ~ response_type*valence + error_rate + error_rate:response_type + error_rate:valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_error_rate_red1)                          
anova(LMM_rt_cov_error_rate,LMM_rt_cov_error_rate_red1)      





###### Include Error Frustration / Avoidance as covariate?? (just had a quick glimpse at it)

# Load ratings and merge to data4mixedmodels
ratings  <- data.frame()  
logfiles <- list.files("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_Behavioral", pattern = ".txt")       

for (subject in logfiles){                                                                               
  setwd("C:/Users/Luisa/PhD/1_PhD_Project/6_ModERN_Behavioral_Study/6_Raw_Data_Behavioral")                 
  rating <- read.table(subject, skip = 575, fill = TRUE, header = TRUE, sep = ":", stringsAsFactors = FALSE)
  subjectID <-  factor(as.numeric(substr(subject,14,15)))
  
  effort            <- as.numeric(rating[1,1])
  error_avoidance   <- as.numeric(rating[2,1])
  error_frustration <- as.numeric(rating[3,1])
  fatigue           <- as.numeric(rating[4,1])
  
  ratings <- rbind(ratings,data.frame(subjectID,effort,error_avoidance,error_frustration,fatigue))
} 
data4mixedmodels_words <- left_join(data4mixedmodels_words,ratings,by = "subjectID")


### error avoidance
# Center Covariate
data4mixedmodels_words$error_avoidance <- scale(data4mixedmodels_words$error_avoidance, center = TRUE, scale = FALSE)

# Full model: response type:valence remains significant; there is also a significant 3 way interaction response_type4-3:valence2-1:error_avoidance
LMM_rt_cov_error_avoidance <- lmer(word_rt_inverse ~ response_type*valence*error_avoidance + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_error_avoidance) # interestingly, there is no main effect of error_avoidance on word_rt
 
# Nested model: priming after FA, FH, SH remains (for no response type there is a significant interaction valence x frustration)
LMM_rt_cov_error_avoidance_nested <- lmer(word_rt_inverse ~ (response_type/valence)*error_avoidance + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_error_avoidance_nested)

# 3-way interaction error_avoidance:valence:response type does not improve fit (no sign. diff.) -> priming not modulated by error avoidance!!!!
LMM_rt_cov_error_avoidance_red1 <- lmer(word_rt_inverse ~ response_type*valence + error_avoidance + error_avoidance:response_type + error_avoidance:valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_error_avoidance_red1)                          
anova(LMM_rt_cov_error_avoidance,LMM_rt_cov_error_avoidance_red1)      


### error frustration
# Center Covariate
data4mixedmodels_words$error_frustration <- scale(data4mixedmodels_words$error_frustration, center = TRUE, scale = FALSE)

# Full model: response type:valence remains significant, there is no  3-way interaction error_frustration:valence:response type
LMM_rt_cov_error_frustration <- lmer(word_rt_inverse ~ response_type*valence*error_frustration + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_error_frustration)

# Nested model: priming after FA, FH, and SH remains (for no response type there is a significant interaction valence x frustration)
LMM_rt_cov_error_frustration_nested <- lmer(word_rt_inverse ~ (response_type/valence)*error_frustration + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_error_frustration_nested)

# 3-way interaction error_frustration:valence:response type does not improve fit (no sign. diff.) -> priming not modulated by error frustration!!!!
LMM_rt_cov_error_frustration_red1 <- lmer(word_rt_inverse ~ response_type*valence + error_frustration + error_frustration:response_type + error_frustration:valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_error_frustration_red1)                          
anova(LMM_rt_cov_error_frustration,LMM_rt_cov_error_frustration_red1)      


### effort
# Center Covariate
data4mixedmodels_words$effort <- scale(data4mixedmodels_words$effort, center = TRUE, scale = FALSE)

# Full model: response type:valence remains significant, there is no  3-way interaction effort:valence:response type
LMM_rt_cov_effort <- lmer(word_rt_inverse ~ response_type*valence*effort + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_effort)

# Nested model: priming after FA, FH, SH remains (for no response type there is a significant interaction valence x frustration)
LMM_rt_cov_effort_nested <- lmer(word_rt_inverse ~ (response_type/valence)*effort + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_effort_nested)

# 3-way interaction effort:valence:response type does not improve fit (no sign. diff.) -> priming not modulated by error frustration!!!!
LMM_rt_cov_effort_red1 <- lmer(word_rt_inverse ~ response_type*valence + effort + effort:response_type + effort:valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_cov_effort_red1)                          
anova(LMM_rt_cov_effort,LMM_rt_cov_effort_red1)      




###################   Generalized Linear Mixed Models for Priming RT   #################### 

###### STEP 1 to 4: RUN CODE ABOVE (PREPARATION FOR PRIMING)


###### STEP 5: SPECIFY THE MAXIMAL MODEL 

GLMM_rt_max <- glmer(word_rt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg | subjectID) + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_max)                         # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
summary(rePCA(GLMM_rt_max))  
print(VarCorr(GLMM_rt_max),comp="Variance")  
# I did not check whether according to didLmerConverge the max model does converge, because as seen below, I have to reduce random structure to achieve identified model and removing the correlations would have been the first step anyway



###### STEP 6: SIMPLIFICATION OF THE RANDOM EFFECTS STRUCTURE 

# I increased number of evaluations to achieve convergence

# start with zero-correlation parameter model by using ||
GLMM_rt_red1 <- glmer(word_rt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
summary(GLMM_rt_red1)                        # model failed to converge (but this is a false positive; I can trust in the estimations)
didLmerConverge((GLMM_rt_red1))              # The relative maximum gradient of 0.000569 is less than our 0.001 criterion. You can safely ignore any warnings about a claimed convergence failure.
summary(rePCA(GLMM_rt_red1))                 # all items explained variance
print(VarCorr(GLMM_rt_red1),comp="Variance") 



# ### now that we have an identified model, test non-significant variance components using likelihood ratio tests (starting with component explaning least variance) - Schad recommends drop 1 strategy, but he would rather only do rePCA

# previously I reduce random structure to achieve identified model (I did not know I can safely ignore the warning here); remove random slope that explains the least variance. It is random slope for response type for words 
# GLMM_rt_red2 <- glmer(word_rt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# summary(GLMM_rt_red2)            # model converged -> is identified
# But: If I use this reduced model as final model, I get problems with the corresponding nested model, which would not converge anymore
# 
# 
# # test whether random slope for valence for subjects improves fit
# GLMM_rt_red3 <- glmer(word_rt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_red3)             # model failed to converge 
# anova(GLMM_rt_red2,GLMM_rt_red3)  # I keep random slope for valence for subjects; I also have it for the LMM and model comparison (may not be valid here) also favors model 2
# 
# # test whether random slope for response_type for subjects improves fit
# GLMM_rt_red4 <- glmer(word_rt ~ response_type*valence + (1 + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# summary(GLMM_rt_red4)             # model does converge -> is identified
# anova(GLMM_rt_red2,GLMM_rt_red4)  # sign. difference -> random slope for response_type for subjects improves fit -> model 2 wins
# 
# # test whether random slope for interaction for subjects improves fit
# GLMM_rt_red5 <- glmer(word_rt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# summary(GLMM_rt_red5)             # model does converge -> identified model
# anova(GLMM_rt_red2,GLMM_rt_red5)  # sign. difference -> random slope for interaction for subjects improves fit -> model 2 wins
# 
# # test whether random intercept for word improves fit
# GLMM_rt_red6 <- glmer(word_rt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# summary(GLMM_rt_red6)             # model failed to converge
# anova(GLMM_rt_red2,GLMM_rt_red6)  # I keep random intercept for word; I also have it for the LMM and model comparison (may not be valid here) also favors model 2
# 
# # test whether random intercept for subject improves fit
# GLMM_rt_red7 <- glmer(word_rt ~ response_type*valence + (0 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# summary(GLMM_rt_red7)             # error: step-halvings failed to reduce deviance in pwrssUpdate
# anova(GLMM_rt_red2,GLMM_rt_red7)  # I keep random intercept for subject; I also have it for the LMM and model comparison (may not be valid here) also favors model 2
# 
# # test whether reintroduction of correlation improves fit 
# GLMM_rt_red8 <- glmer(word_rt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg | subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# summary(GLMM_rt_red8)             # model failed to converge
# 
# 
# 
# ###### STEP 7: TEST OMNIBUS SIGNIFICANCE OF FIXED EFFECTS (based on model reduced by rePCA and LRT)
# 
# # test whether fixed effect of interaction response_type*valence improves fit
# GLMM_rt_red9 <- glmer(word_rt ~ response_type + valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# summary(GLMM_rt_red9)            # model failed to converge
# 
# # test whether fixed effect of valence improves fit
# GLMM_rt_red10 <- glmer(word_rt ~ response_type + response_type:valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# summary(GLMM_rt_red10)           # model failed to converge
# 
# # test whether fixed effect of response_type improves fit
# GLMM_rt_red11 <- glmer(word_rt ~ valence + response_type:valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# summary(GLMM_rt_red11)           # model failed to converge



###### STEP 8: PRESENT FINAL MODEL (BY USING REML ESTIMATIONS FOR LMMS) 

GLMM_rt_final <- glmer(word_rt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
summary(GLMM_rt_final)           # model failed to converge (but this is a false positive; I can trust in the estimations)

didLmerConverge(GLMM_rt_final)
# The relative maximum gradient of 0.000569 is less than our 0.001 criterion. You can safely ignore any warnings about a claimed convergence failure.


###### STEP 9: BREAK DOWN INTERAKTIONS BY USING NESTED MODEL 

# valence within response condition
GLMM_rt_nested <- glmer(word_rt ~ response_type/valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA || word), data=data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
summary(GLMM_rt_nested)           # model failed to converge (but this is a false positive; I can trust in the estimations)
# difflsmeans does not work for glmm

didLmerConverge(GLMM_rt_nested)
# The relative maximum gradient of 0.000865 is less than our 0.001 criterion. You can safely ignore any warnings about a claimed convergence failure.




###### STEP 10: CHECK ASSUMPTIONS 

# for GLMMs, testing assumptions is not straightforward, see: https://www.theanalysisfactor.com/regression-diagnostics-glmm/
# Assumption: The chosen link function is appropriate - Is ok, see inverse gaussian distribution fit 
# Assumption: Random effects (intercepts and slopes) are normally distributed - Is only partly ok. As can be seen in the Q-Q plot, there is 1 outlier that leads to deviation from normal distribution of random intercept of subjects. Inspection of the random effect revealed that subject 28 (very high intercept, very slow responses to words in general) is this outlier. Excluding this subject may solve this problem, but there is no other reason to exclude this subject. For this random intercept, the Shapiro-Wilk test is significant. 
# Assumption: Appropriate estimation of variance / No overdispersion -  Overdispersion is irrelevant for models that estimate a scale parameter (i.e. almost anything but Poisson or binomial: Gaussian, Gamma, negative binomial); see: https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html

r_int_subjects <- ranef(GLMM_rt_final_converge)$subjectID$`(Intercept)`
qqnorm(r_int_subjects, main = "Q-Q Plot Random Intercept Subjects")
qqline(r_int_subjects)
shapiro_int_subjects <- shapiro.test(r_int_subjects)["p.value"]
#ranef(GLMM_rt_final_converge)$subjectID$`(Intercept)`: 1 subject (Ss 28) is an extreme outlier; rest is normally distributed, but Shapiro is sign.

r_int_words <- ranef(GLMM_rt_final_converge)$word$`(Intercept)`
qqnorm(r_int_words, main = "Q-Q Plot Random Intercept Words")
qqline(r_int_words)
shapiro_int_words <- shapiro.test(r_int_words)["p.value"]

r_slope1_subjects <- ranef(GLMM_rt_final_converge)$subjectID[,'pos_minus_neg:CI_minus_FA']
qqnorm(r_slope1_subjects, main = "Q-Q Plot Slope CI_minus_FA:pos_minus_neg")
qqline(r_slope1_subjects)
shapiro_slope1_subjects <- shapiro.test(r_slope1_subjects)["p.value"] 

r_slope2_subjects <- ranef(GLMM_rt_final_converge)$subjectID[,'pos_minus_neg:FA_minus_FH']
qqnorm(r_slope2_subjects, main = "Q-Q Plot Slope FA_minus_FH:pos_minus_neg")
qqline(r_slope2_subjects)
shapiro_slope2_subjects <- shapiro.test(r_slope2_subjects)["p.value"] 

r_slope3_subjects <- ranef(GLMM_rt_final_converge)$subjectID[,'FH_minus_SH:pos_minus_neg']
qqnorm(r_slope3_subjects, main = "Q-Q Plot Slope FH_minus_SH:pos_minus_neg")
qqline(r_slope3_subjects)
shapiro_slope3_subjects <- shapiro.test(r_slope3_subjects)["p.value"] 








###################   Generalized Linear Mixed Models for Priming Accuracy    #################### 

###### STEP 1 to 4: RUN CODE ABOVE (PREPARATION FOR PRIMING)


###### STEP 5: SPECIFY THE MAXIMAL MODEL 

GLMM_acc_max <- glmer(word_accuracy ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg | subjectID) + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA | word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
summary(GLMM_acc_max) # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
# I did not check whether according to didLmerConverge the max model does converge, because as seen below, I have to reduce it for some zero variance components and removing the correlations would have been the first step anyway


###### STEP 6: SIMPLIFICATION OF THE RANDOM EFFECTS STRUCTURE 

# start with zero-correlation parameter model by using ||
GLMM_acc_red1 <- glmer(word_accuracy ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA || word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
summary(GLMM_acc_red1)                          # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
summary(rePCA(GLMM_acc_red1))                   # for subject there is one item that does not explain variance
print(VarCorr(GLMM_acc_red1),comp="Variance")   # it is FA_minus_FH

# the model still converges -> it would also be ok to use the maximal model, as for the full model vs the reduced model all estimates are quite identical effect remain the same (see argument Singman and Kellen: lower-order effects should not be removed when higher-order effects are still included; they recommend to rather keep zero variance random effect, if the estimates are quite the same)
didLmerConverge(GLMM_acc_red1) # 	The relative maximum gradient of 0.0000142 is less than our 0.001 criterion. You can safely ignore any warnings about a claimed convergence failure.
# strategy to take out corresponding higher-order effect FA_minus_FH:pos_minus_neg as well does not work, because it has much variance

# remove FA_minus_FH from subject
GLMM_acc_red2 <- glmer(word_accuracy ~ response_type*valence + (1 + FH_minus_SH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA || word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
summary(GLMM_acc_red2)                          # model does converge, no singular fit -> model is identified
summary(rePCA(GLMM_acc_red2))                   # all items explain variance
print(VarCorr(GLMM_acc_red2),comp="Variance")  


# ### now that we have an identified model, test non-significant variance components using likelihood ratio tests (starting with component explaning least variance) - Schad recommends drop 1 strategy, but he would rather only do rePCA
# 
# # test whether random slope reaction_type for word improves fit (explains least variance)
# GLMM_acc_red3 <- glmer(word_accuracy ~ response_type*valence + (1 + FH_minus_SH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_acc_red3)                        # model does converge, no singular fit -> model is identified
# anova(GLMM_acc_red2,GLMM_acc_red3)            # no sign. difference, random slope reaction_type for word does not improve fit -> model 3 wins
# 
# # test whether random slope reaction_type for subject improves fit (now explains least variance)
# GLMM_acc_red4 <- glmer(word_accuracy ~ response_type*valence + (1 + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_acc_red4)                        # model does converge, no singular fit -> model is identified
# anova(GLMM_acc_red3,GLMM_acc_red4)            # no sign. difference, random slope reaction_type for subject does not improve fit -> model 4 wins
# 
# # test whether random slope valence for subject improves fit (now explains least variance)
# GLMM_acc_red5 <- glmer(word_accuracy ~ response_type*valence + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_acc_red5)                        # model does converge, no singular fit -> model is identified
# anova(GLMM_acc_red4,GLMM_acc_red5)            # random slope valence for subject improves fit -> model 4 wins
# 
# # test whether random slope interaction valence*reaction_type for subject improves fit 
# GLMM_acc_red6 <- glmer(word_accuracy ~ response_type*valence + (1 + pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_acc_red6)                        # model does converge, no singular fit -> model is identified
# anova(GLMM_acc_red4,GLMM_acc_red6)            # random slope interaction response_type*valence for subject improves fit -> model 4 wins
# 
# # test whether random intercept of word improves fit
# GLMM_acc_red7 <- glmer(word_accuracy ~ response_type*valence + (1 + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_acc_red7)                        # model does converge, no singular fit -> model is identified
# anova(GLMM_acc_red4,GLMM_acc_red7)            # random intercept for word does improve fit -> model 4 wins
# 
# # test whether random intercept of subject improves fit
# GLMM_acc_red8 <- glmer(word_accuracy ~ response_type*valence + (0 + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_acc_red8)                        # model does converge, no singular fit -> model is identified
# anova(GLMM_acc_red4,GLMM_acc_red8)            # random intercept for subject improves fit -> model 4 wins 
# 
# # test whether reintroduction of correlation improves fit
# GLMM_acc_red9 <- glmer(word_accuracy ~ response_type*valence + (1 + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg | subjectID) + (1 | word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_acc_red9)                        # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits
# 
# 
# 
# ###### STEP 7: TEST OMNIBUS SIGNIFICANCE OF FIXED EFFECTS (based on model reduced by rePCA and LRT)
# 
# # test whether fixed effect of interaction response_type*valence improves fit
# GLMM_acc_red10 <- glmer(word_accuracy ~ response_type + valence + (1 + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_acc_red10)                      # model does converge, no singular fit -> model is identified
# anova(GLMM_acc_red4,GLMM_acc_red10)          # fixed effect of interaction response_type*valence improves fit 
# 
# # test whether fixed effect of valence improves fit
# GLMM_acc_red11 <- glmer(word_accuracy ~ response_type + response_type:valence + (1 + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_acc_red11)                      # model does converge, no singular fit -> model is identified
# anova(GLMM_acc_red4,GLMM_acc_red11)          # fixed effect of valence does not improve fit (AIC identical, but sign. difference)
# 
# # test whether fixed effect of response_type improves fit
# GLMM_acc_red12 <- glmer(word_accuracy ~ valence + response_type:valence + (1 + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_acc_red12)                      # model does converge, no singular fit -> model is identified
# anova(GLMM_acc_red4,GLMM_acc_red12)          # fixed effect of response_type does not improve fit 



###### STEP 8: PRESENT FINAL MODEL (BY USING REML ESTIMATIONS FOR LMMS) 
GLMM_acc_final <- glmer(word_accuracy ~ response_type*valence + (1 + FH_minus_SH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA || word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
summary(GLMM_acc_final)                        



###### STEP 9: BREAK DOWN INTERAKTIONS BY USING NESTED MODEL 

# valence within response condition
GLMM_acc_nested <- glmer(word_accuracy ~ response_type/valence + (1 + FH_minus_SH + CI_minus_FA + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + CI_minus_FA:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH + CI_minus_FA || word), data=data4mixedmodels_words, family = binomial, control=glmerControl(optimizer="bobyqa"))
summary(GLMM_acc_final)      # model fails to converge, unless optimizer="bobyqa" is used              
# difflsmeans does not work for glmm




###### STEP 10: CHECK ASSUMPTIONS 

# for GLMMs, testing assumptions is not straightforward, see: https://www.theanalysisfactor.com/regression-diagnostics-glmm/
# Assumption 1: The chosen link function is appropriate -> Is ok
# Assumption 2: Random effects (intercepts and slopes) are normally distributed -> Is ok only for most effects. 
                # As can be seen in the Q-Q plot, there are 3 outliers that lead to deviation from normal distribution of random intercept of words.
                # Inspection of the random effect revealed that the words "einsam", "herzlos", and "lieblos" are these outliers (these words were related to slow responses). 
                # For the random intercept for words, the Shapiro-Wilk test is significant. It is also significant for the random slope pos_minus_neg for subjects, but this ransom slope only seems to mildly deviate from normal distribution (subjects 3 and 4 are mild outliers here).
# Assumption 3: Appropriate estimation of variance (no overdispersion) -> Overdispersion is not a problem here. 
                # The empirical variance in data does not exceed the nominal variance under the presumed model. The chi-square test of the ratio of the empirical variance in data and the nominal variance under the presumed model is not significant:
  
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(GLMM_acc_final) # If the p-value is < 0.05, the data are overdispersed. Here p > 0.05. Overdispersion is not a problem here. 

# Assumption: Random effects come from a normal distribution - ok for most effects
r_int_subjects <- ranef(GLMM_acc_final)$subjectID$`(Intercept)`
qqnorm(r_int_subjects, main = "Q-Q Plot Random Intercept Subjects")
qqline(r_int_subjects)
shapiro_int_subjects <- shapiro.test(r_int_subjects)["p.value"]

r_slope1_subjects <- ranef(GLMM_acc_final)$subjectID$'FH_minus_SH'
qqnorm(r_slope1_subjects, main = "Q-Q Plot Slope FH_minus_SH")
qqline(r_slope1_subjects)
shapiro_slope1_subjects <- shapiro.test(r_slope1_subjects)["p.value"] 

r_slope2_subjects <- ranef(GLMM_acc_final)$subjectID$'CI_minus_FA'
qqnorm(r_slope2_subjects, main = "Q-Q Plot Slope CI_minus_FA")
qqline(r_slope2_subjects)
shapiro_slope2_subjects <- shapiro.test(r_slope2_subjects)["p.value"] 

r_slope3_subjects <- ranef(GLMM_acc_final)$subjectID$pos_minus_neg
qqnorm(r_slope3_subjects, main = "Q-Q Plot Slope pos_minus_neg")
qqline(r_slope3_subjects)
shapiro_slope3_subjects <- shapiro.test(r_slope3_subjects)["p.value"] 
# ranef(GLMM_acc_final)$subjectID$pos_minus_neg: 2 subjects (Ss 3 and 4) are mild outliers; rest is normally distributed, but Shapiro is sign., but only with 0.04

r_slope4_subjects <- ranef(GLMM_acc_final)$subjectID$'FH_minus_SH:pos_minus_neg'
qqnorm(r_slope4_subjects, main = "Q-Q Plot Slope FH_minus_SH:pos_minus_neg")
qqline(r_slope4_subjects)
shapiro_slope4_subjects <- shapiro.test(r_slope4_subjects)["p.value"] 

r_slope5_subjects <- ranef(GLMM_acc_final)$subjectID$'pos_minus_neg:FA_minus_FH'
qqnorm(r_slope5_subjects, main = "Q-Q Plot Slope FA_minus_FH:pos_minus_neg")
qqline(r_slope5_subjects)
shapiro_slope5_subjects <- shapiro.test(r_slope5_subjects)["p.value"] 

r_slope6_subjects <- ranef(GLMM_acc_final)$subjectID$'CI_minus_FA:pos_minus_neg'
qqnorm(r_slope6_subjects, main = "Q-Q Plot Slope CI_minus_FA:pos_minus_neg")
qqline(r_slope6_subjects)
shapiro_slope6_subjects <- shapiro.test(r_slope6_subjects)["p.value"] 
# ranef(GLMM_acc_final)$subjectID$'CI_minus_FA:pos_minus_neg': 1 subject (Ss 5) is outlier; rest is normally distributed, but Shapiro is sign., but only with 0.04

r_int_words <- ranef(GLMM_acc_final)$word$`(Intercept)`
qqnorm(r_int_words, main = "Q-Q Plot Random Intercept Words")
qqline(r_int_words)
shapiro_int_words <- shapiro.test(r_int_words)["p.value"] 
#ranef(GLMM_acc_final)$word$`(Intercept)`: 3 words (einsam, herzlos, lieblos) are outliers; rest is normally distributed, but Shapiro is sign.

r_slope1_words <- ranef(GLMM_acc_final)$word$FH_minus_SH
qqnorm(r_slope1_words, main = "Q-Q Plot Slope FH_minus_SH")
qqline(r_slope1_words)
shapiro_slope1_words <- shapiro.test(r_slope1_words)["p.value"] 

r_slope2_words <- ranef(GLMM_acc_final)$word$FA_minus_FH
qqnorm(r_slope2_words, main = "Q-Q Plot Slope FA_minus_FH")
qqline(r_slope2_words)
shapiro_slope2_words <- shapiro.test(r_slope2_words)["p.value"] 

r_slope3_words <- ranef(GLMM_acc_final)$word$CI_minus_FA
qqnorm(r_slope3_words, main = "Q-Q Plot Slope CI_minus_FA")
qqline(r_slope3_words)
shapiro_slope3_words <- shapiro.test(r_slope3_words)["p.value"] 




###################   Preparation for Priming and SCR (and for pure SCR)   #################### 


###### STEP 1: PREPARE DATA 

# I decided to report 2 separate models for word rt and word rt with SCR. For the first I will include CI to show that there is no priming effect after CI, for the latter I will exclude CI, because there is no valid SCR and I want to exclude CI from all SCR analyses (starting from calculating z standardized SCR values which are calculated based on only valid trials)
# use df with only valid scr trials: I only included FA, FH, SH trials that had no invalid gng rt, are followed by correct word classification and are not preceded or followed by an incorrect response or wrong key in gng or word classification.
# all SCR analyses excluded CIs, because this is the only condition not confounded by motor response
# exclude word and GNG responses with misses or wrong keys, outliers in word RT or GNG RT; incorrect responses for RT LMMs are excluded when the LMM is specified
data4mixedmodels_words_scr <- data4mixedmodels[data4mixedmodels$outlier_words == FALSE & data4mixedmodels$gng_resp <= 44 & data4mixedmodels$gng_invalid_rt == FALSE & data4mixedmodels$word_resp <= 54 & data4mixedmodels$followed_or_preceded_by_FA_or_wrong_key == FALSE & !is.na(data4mixedmodels$iscr_gng_resp),]
# 7672 of 15480 trials left
data4mixedmodels_words_scr$iscr_gng_resp_sqrt <- sqrt(data4mixedmodels_words_scr$iscr_gng_resp) # this is not needed for the priming-scr analysis, but as dependent variable for the pure scr analysis (see below)


# make categorical variables factors
data4mixedmodels_words_scr$response_type <- factor(data4mixedmodels_words_scr$response_type, levels=c("SH","FH","FA"))
data4mixedmodels_words_scr$valence       <- factor(data4mixedmodels_words_scr$valence)
data4mixedmodels_words_scr$condition     <- factor(data4mixedmodels_words_scr$condition, levels=c("neg_after_SH","pos_after_SH","neg_after_FH","pos_after_FH","neg_after_FA","pos_after_FA"))
data4mixedmodels_words_scr$subjectID     <- factor(data4mixedmodels_words_scr$subjectID)
data4mixedmodels_words_scr$word          <- factor(data4mixedmodels_words_scr$word)
data4mixedmodels_words_scr$word_accuracy <- factor(data4mixedmodels_words_scr$word_accuracy)



###### STEP 2: CHECk DISTRIBUTION 

# transformed rt -> normally distributed, also in the different conditions
plot(density(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$word_rt_inverse), main = "Histogram Word RT Inverse")
ggplot(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], aes(x = word_rt_inverse)) + geom_density() + facet_wrap(response_type ~ valence) + ggtitle("Histogram Word RT Inverse in Conditions")
qqp(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$word_rt_inverse, "norm", main = "Q-Q Plot Word RT Inverse", ylab = "sample quantiles")

# raw rt -> not normally distributed
plot(density(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$word_rt), main = "Histogram Word RT Raw")
ggplot(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], aes(x = word_rt)) + geom_density() + facet_wrap(response_type ~ valence) + ggtitle("Histogram Word RT Raw in Conditions")
qqp(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$word_rt, "norm", main = "Q-Q Plot Word RT Raw", ylab = "sample quantiles")

# raw rt -> fit gamma
gamma <- fitdistr(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$word_rt, "gamma")
qqp(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$word_rt, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]], ylab = "sample quantiles",main = "Q-Q Plot Word RT gamma")

# raw rt -> fit inverse gaussian -> inverse gaussian distribution fits quite well (for RTs, Lo and Anderson 2015 recommend gamma or inverse gaussian, Conny used Gamma)
inv_gauss <- nigFit(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$word_rt, plots = FALSE, printOut = FALSE) 

par(mfrow = c(1,2))  # use par function in base graphics approach to arrange plots side by side here; grid.arrange only works for ggplot grahpics
plot(inv_gauss, which = 1, plotTitles = "Histogram of Word RT inv. gaussian")
plot(inv_gauss, which = 3, plotTitles = "Q-Q Plot of Word RT inv. gaussian")
par(mfrow = c(1, 1)) # reset "par" parameter

# I would choose inverse gaussian distribution. Lo & Andrews (2015) recommend gamma or inverse gaussian distribution and identity link. For my data, inverse gaussian distribution fits better.


###### STEP 3: DEFINE CONTRASTS 

contrasts(data4mixedmodels_words_scr$response_type) <- contr.sdif(3)
contrasts(data4mixedmodels_words_scr$valence)       <- contr.sdif(2)

# add contrast as numerical covariate via model matrix (in order to enter contrasts individually in model and hence to use ||)
mm_c    <- model.matrix( ~ response_type*valence, data4mixedmodels_words_scr) # stick in model matrix 6 columns
glimpse(mm_c)                                                                 # requires dplyr
glimpse(data4mixedmodels_words_scr)                                           # at the moment has 25 columns

# attach to dataframe 
data4mixedmodels_words_scr[,(ncol(data4mixedmodels_words_scr)+1):(ncol(data4mixedmodels_words_scr)+6)] <- mm_c
names(data4mixedmodels_words_scr)[(ncol(data4mixedmodels_words_scr)-5):ncol(data4mixedmodels_words_scr)] <- c("Grand Mean", "FH_minus_SH","FA_minus_FH","pos_minus_neg", "FH_minus_SH:pos_minus_neg", "FA_minus_FH:pos_minus_neg")
glimpse(data4mixedmodels_words_scr)


###### STEP 4: SCALE CONTINUOUS PREDICTORS / COVARIATES 

# not necessary: scr is already scaled :)



###################   Linear Mixed Models for Priming and SCR   #################### 

###### STEP 1 to 4: RUN CODE ABOVE (PREPARATION FOR PRIMING AND SCR)

###### STEP 5: SPECIFY THE MAXIMAL MODEL 
# It would be ok to not include covariate in random structure, because this effect is not of major interest and we only want to control for it, not generalize it over subjects / words (Kliegl often does not include it, Schad prefers including it)
LMM_rt_scr_max <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  | subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_scr_max)                         # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
summary(rePCA(LMM_rt_scr_max))                  # for word there are three items and for subjectID there are five items that do not explain variance
print(VarCorr(LMM_rt_scr_max),comp="Variance")  

# the model still converges (for the fuller model, the trend of main effect SCR is less pronounced) 
didLmerConverge(LMM_rt_scr_max) # The relative maximum gradient of 0.000149 is less than our 0.001 criterion. You can safely ignore any warnings about a claimed convergence failure.



###### STEP 6: SIMPLIFICATION OF THE RANDOM EFFECTS STRUCTURE 

# start with zero-correlation parameter model by using ||
LMM_rt_scr_red1 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_scr_red1)                         # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
summary(rePCA(LMM_rt_scr_red1))                  # for word there are two items and for subjectID there are three items that do not explain variance
print(VarCorr(LMM_rt_scr_red1),comp="Variance")  # it is FA_minus_FH:iscr_gng_resp_sqrt_z_score, FH_minus_SH:iscr_gng_resp_sqrt_z_score for word and FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score, pos_minus_neg:iscr_gng_resp_sqrt_z_score, FA_minus_FH:iscr_gng_resp_sqrt_z_score for subjectID

# remove FH_minus_SH:iscr_gng_resp_sqrt_z_score and FA_minus_FH:iscr_gng_resp_sqrt_z_score for word and FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score, pos_minus_neg:iscr_gng_resp_sqrt_z_score, FA_minus_FH:iscr_gng_resp_sqrt_z_score for subjectID (is same result as done in separate steps)
LMM_rt_scr_red2 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_scr_red2)                         # model does converge -> identified model  
summary(rePCA(LMM_rt_scr_red2))                  # all items explain variance          
print(VarCorr(LMM_rt_scr_red2),comp="Variance")  

# I could also leave out FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score for subject as they also explain < 0.5% of variance; then I do not habe problem that higher-order effects are included and lower-order effects are excluded...
LMM_rt_scr_red3 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score || subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_scr_red3)                         # model does converge -> identified model                    
summary(rePCA(LMM_rt_scr_red3))                  # all items explain variance      
print(VarCorr(LMM_rt_scr_red3),comp="Variance")  



# ### now that we have an identified model, test non-significant variance components using likelihood ratio tests (starting with component explaning least variance) - Schad recommends drop 1 strategy, but he would rather only do rePCA
# 
# # test whether random slope scr for word improves fit (explains least variance)
# LMM_rt_scr_red4 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red4)                         # model does converge -> identified model
# anova(LMM_rt_scr_red3,LMM_rt_scr_red4)           # no sign. difference, random slope scr for word does not improve fit -> model 4 wins
# 
# # test whether random slope FH_minus_SH:iscr_gng_resp_sqrt_z_score for subject improves fit (explains least variance)
# LMM_rt_scr_red5 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red5)                         # model does converge -> identified model
# anova(LMM_rt_scr_red4,LMM_rt_scr_red5)           # no sign. difference, random slope FH_minus_SH:iscr_gng_resp_sqrt_z_score for subject does not improve fit -> model 5 wins
# 
# # test whether random slope scr for subject improves fit 
# LMM_rt_scr_red6 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red6)                         # model does converge -> identified model
# anova(LMM_rt_scr_red5,LMM_rt_scr_red6)           # sign. difference, random slope scr for subject does improve fit  -> model 5 wins
# 
# # test whether random slope FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score for subject improves fit 
# LMM_rt_scr_red7 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red7)                         # model does converge -> identified model
# anova(LMM_rt_scr_red5,LMM_rt_scr_red7)           # no sign. difference, random slope FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score for subject does not improve fit -> model 7 wins
# 
# # test whether random slope response_type for word improves fit 
# LMM_rt_scr_red8 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1  | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red8)                         # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red8)           # sign. difference, random slope response_type for word does improve fit -> model 7 wins
# 
# # test whether random slope valence for subject improves fit 
# LMM_rt_scr_red9 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red9)                         # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red9)           # sign. difference, random slope valence for subject does improve fit -> model 7 wins
# 
# # test whether random slope response_type for subject improves fit 
# LMM_rt_scr_red10 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red10)                        # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red10)          # sign. difference, random slope response_type for subject does improve fit -> model 7 wins
# 
# # test whether random slope interaction valence*response_type for subject improves fit 
# LMM_rt_scr_red11 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red11)                        # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red11)          # sign. difference, random interaction valence*response_type for subject does improve fit -> model 7 wins
# 
# # test whether random intercept for word improves fit
# LMM_rt_scr_red12 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (0 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red12)                        # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red12)          # sign. difference, random intercept for word does improve fit -> model 7 wins
# 
# # test whether random intercept for subject improves fit
# LMM_rt_scr_red13 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (0 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red13)                        # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red13)          # sign. difference, random intercept for subject does improve fit -> model 7 wins
# 
# # test whether reintroduction of correlation improves fit
# LMM_rt_scr_red14 <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  | subjectID) + (1 + FH_minus_SH + FA_minus_FH  | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red14)                        # model failed to converge
# 
# 
# 
# ###### STEP 7: TEST OMNIBUS SIGNIFICANCE OF FIXED EFFECTS (based on model reduced by rePCA and LRT) 
# 
# # test whether fixed effect of interaction response_type*valence*scr improves fit
# LMM_rt_scr_red15 <- lmer(word_rt_inverse ~ FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score +  (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red15)                         # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red15)           # no sign. difference, interaction response_type*valence*scr does not improve fit!!!!
# 
# # test whether fixed effect of interaction valence*scr improves fit
# LMM_rt_scr_red16 <- lmer(word_rt_inverse ~ FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red16)                         # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red16)           # no sign. difference; fixed effect of interaction valence*scr does not improve fit 
# 
# # test whether fixed effect of interaction response_type*scr improves fit
# LMM_rt_scr_red17 <- lmer(word_rt_inverse ~ FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red17)                         # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red17)           # sign. difference; fixed effect of interaction response_type*scr does improve fit 
# 
# # test whether fixed effect of interaction response_type*valence improves fit
# LMM_rt_scr_red18 <- lmer(word_rt_inverse ~ FH_minus_SH + FA_minus_FH + pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red18)                         # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red18)           # sign. difference; fixed effect of interaction response_type*valence does improve fit 
# 
# # test whether fixed effect of scr improves fit
# LMM_rt_scr_red19 <- lmer(word_rt_inverse ~ FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red19)                         # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red19)           # no sign. difference; fixed effect of scr does not improve fit
# 
# # test whether fixed effect of valence improves fit
# LMM_rt_scr_red20 <- lmer(word_rt_inverse ~ FH_minus_SH + FA_minus_FH + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red20)                         # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red20)           # no sign. difference; fixed effect of valence does not improve fit
# 
# # test whether fixed effect of response_type improves fit
# LMM_rt_scr_red21 <- lmer(word_rt_inverse ~ pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_scr_red21)                         # model does converge -> identified model
# anova(LMM_rt_scr_red7,LMM_rt_scr_red21)           # sign. difference; fixed effect of response_type does improve fit




###### STEP 8: PRESENT FINAL MODEL (BY USING REML ESTIMATIONS FOR LMMS) 

LMM_rt_scr_final <- lmer(word_rt_inverse ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score || subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = TRUE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_scr_final)                        # model does converge, no singular fit -> model is identified



###### STEP 9: BREAK DOWN INTERACTIONS BY USING NESTED MODEL 

# valence, scr, valence*scr within response condition 
LMM_rt_scr_nested <- lmer(word_rt_inverse ~ response_type/(valence*iscr_gng_resp_sqrt_z_score) + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = TRUE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_scr_nested)                       # model does converge, no singular fit -> model is identified



###### STEP 10: CHECK ASSUMPTIONS 

# LINEARITY & HOMOGENEITY OF VARIANCE
plot(fitted(LMM_rt_scr_final), residuals(LMM_rt_scr_final), xlab = "Fitted Values", ylab = "Residuals")  # plot residuals against the fitted values
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(LMM_rt_scr_final), residuals(LMM_rt_scr_final)))
# scatter looks like a random pattern ->  residuals seem to meet model assumption of normality
# solid line covers the dashed line -> best-fit line fits well
# if the scatter is not random that means there's some other random or fixed effect that explains variation in your data that you're missing
# alternative: plot(LMM_rt_scr_final)


# NORMALITY OF RESIDUALS
qqPlot(resid(LMM_rt_scr_final)) 
qqnorm(resid(LMM_rt_scr_final))
qqline(resid(LMM_rt_scr_final)) # a bit off at the extremes, but that's often the case; again doesn't look too bad





###################   Generalized Linear Mixed Models for Priming and SCR   #################### 

###### STEP 1 to 4: RUN CODE ABOVE (PREPARATION FOR PRIMING AND SCR)


###### STEP 5: SPECIFY THE MAXIMAL MODEL 

GLMM_rt_scr_max <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  | subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_max)                         # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
summary(rePCA(GLMM_rt_scr_max))  
print(VarCorr(GLMM_rt_scr_max),comp="Variance")  


###### STEP 6: SIMPLIFICATION OF THE RANDOM EFFECTS STRUCTURE 

# start with zero-correlation parameter model by using ||
GLMM_rt_scr_red1 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_red1)                          # model failed to converge 
summary(rePCA(GLMM_rt_scr_red1))                   # all items explain variance
print(VarCorr(GLMM_rt_scr_red1),comp="Variance") 

# the model still converges (I did not try whether the function would say that the correlation model also converges, because it already took ages to run this model and the correlation model is possibly too complex for a GLMM to converge anyway)
didLmerConverge(GLMM_rt_scr_red1) # The relative maximum gradient of 0.000138 is less than our 0.001 criterion. You can safely ignore any warnings about a claimed convergence failure.

# but the corresponding nested model would really not converge. So I simplify it to achieve a converging nested model (I did not quite get there).
GLMM_rt_scr_red_for_nested <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg  || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_red_for_nested)                          # model failed to converge ("This looks like a real convergence failure; maybe try simplifying your model?") AND estimates are strange.
summary(rePCA(GLMM_rt_scr_red_for_nested))                   # all items explain variance
print(VarCorr(GLMM_rt_scr_red_for_nested),comp="Variance") 

# ... I give up

# reduce random structure to achieve identified model; remove random slope that explain least variance. It is pos_minus_neg:iscr_gng_resp_sqrt_z_score for subjectID.
GLMM_rt_scr_red2 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_red2)                          # model failed to converge

# reduce random structure to achieve identified model; remove random slope that explain least variance. It is iscr_gng_resp_sqrt_z_score for word and for subjectID.
GLMM_rt_scr_red3 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_red3)                          # model failed to converge

# reduce random structure to achieve identified model; remove random slope that explain least variance. It is response_type*scr for word.
GLMM_rt_scr_red4 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_red4)                          # model failed to converge

# reduce random structure to achieve identified model; remove random slope that explain least variance. It is valence for subjectID.
GLMM_rt_scr_red5 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_red5)                          # model failed to converge

# reduce random structure to achieve identified model; remove random slope that explain least variance. It is response_type*scr for subjectID.
GLMM_rt_scr_red6 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_red6)                          # model failed to converge

# reduce random structure to achieve identified model; remove random slope that explain least variance. It is response_type for word.
GLMM_rt_scr_red7 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_red7)                          # model failed to converge

# reduce random structure to achieve identified model; remove random slope that explain least variance. It is response_type*scr*valence for subjectID.
GLMM_rt_scr_red8 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH + FA_minus_FH + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_red8)                          # model failed to converge

# reduce random structure to achieve identified model; remove random slope that explain least variance. It is response_type for subjectID.
# increasing number of number of evaluations allows the model to converge up from this point (did not work in models before)
GLMM_rt_scr_red9 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(GLMM_rt_scr_red9)                          # model does converge -> identified model


### now that we have an identified model, test non-significant variance components using likelihood ratio tests (starting with component explaning least variance) - Schad recommends drop 1 strategy, but he would rather only do rePCA

# # test whether response_type for subjectID improves fit
# GLMM_rt_scr_red10 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 | subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_scr_red10)                         # model does converge -> identified model
# anova(GLMM_rt_scr_red9,GLMM_rt_scr_red10)          # sign. difference; response_type for subjectID improves fit -> model 9 wins
# 
# # test whether intercept for subjectID improves fit
# GLMM_rt_scr_red11 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (0 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_scr_red11)                         # model failed to converge
# anova(GLMM_rt_scr_red9,GLMM_rt_scr_red11)          # sign. difference; intercept for subjectID improves fit -> model 9 wins
# 
# # test whether intercept for word improves fit
# GLMM_rt_scr_red12 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_scr_red12)                         # model does converge -> identified model
# anova(GLMM_rt_scr_red9,GLMM_rt_scr_red12)          # sign. difference; intercept for word improves fit -> model 9 wins
# 
# # test whether reintroduction of correlation improves fit 
# GLMM_rt_scr_red13 <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg | subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_scr_red13)                         # model fails to converge 
# 
# 
# 
# 
# ###### STEP 7: TEST OMNIBUS SIGNIFICANCE OF FIXED EFFECTS (based on model reduced by rePCA and LRT)
# 
# # test whether fixed effect of interaction response_type*valence*scr improves fit
# GLMM_rt_scr_red14 <- glmer(word_rt ~ FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_scr_red14)                         # model does converge -> identified model
# anova(GLMM_rt_scr_red9,GLMM_rt_scr_red14)          # no sign. difference, interaction response_type*valence*scr does not improve fit!!!!
# 
# # test whether fixed effect of interaction valence*scr improves fit
# GLMM_rt_scr_red15 <- glmer(word_rt ~ FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_scr_red15)                         # model does converge -> identified model
# anova(GLMM_rt_scr_red9,GLMM_rt_scr_red15)          # no sign. difference; fixed effect of interaction valence*scr does not improve fit 
# 
# # test whether fixed effect of interaction response_type*scr improves fit
# GLMM_rt_scr_red16 <- glmer(word_rt ~ FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_scr_red16)                         # model failed to converge 
# anova(GLMM_rt_scr_red9,GLMM_rt_scr_red16)         
# 
# # test whether fixed effect of interaction response_type*valence improves fit
# GLMM_rt_scr_red17 <- glmer(word_rt ~ FH_minus_SH + FA_minus_FH + pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_scr_red17)                         # model does converge -> identified model
# anova(GLMM_rt_scr_red9,GLMM_rt_scr_red17)          # sign. difference; fixed effect of interaction response_type*valence does improve fit 
# 
# # test whether fixed effect of scr improves fit
# GLMM_rt_scr_red18 <- glmer(word_rt ~ FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_scr_red18)                         # model failed to converge 
# anova(GLMM_rt_scr_red9,GLMM_rt_scr_red18)          
# 
# # test whether fixed effect of valence improves fit
# GLMM_rt_scr_red19 <- glmer(word_rt ~ FH_minus_SH + FA_minus_FH + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_scr_red19)                         # model failed to converge 
# anova(GLMM_rt_scr_red9,GLMM_rt_scr_red19)          
# 
# # test whether fixed effect of response_type improves fit
# GLMM_rt_scr_red20 <- glmer(word_rt ~ pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_rt_scr_red20)                         # model does converge -> identified model
# anova(GLMM_rt_scr_red9,GLMM_rt_scr_red20)          # sign. difference; fixed effect of response_type does improve fit



###### STEP 8: PRESENT FINAL MODEL (BY USING REML ESTIMATIONS FOR LMMS) 

GLMM_rt_scr_final <- glmer(word_rt ~ response_type*valence*iscr_gng_resp_sqrt_z_score + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(GLMM_rt_scr_final)                        # model does converge -> identified model



###### STEP 9: BREAK DOWN INTERAKTIONS BY USING NESTED MODEL 

# valence, scr, valence*scr within response condition 
GLMM_rt_scr_nested <- glmer(word_rt ~ response_type/valence  + response_type/iscr_gng_resp_sqrt_z_score + response_type/(valence*iscr_gng_resp_sqrt_z_score) + (1 + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e15)))
summary(GLMM_rt_scr_nested)           # model does converge, no singular fit -> model is identified
# difflsmeans does not work for glmm

# Here I tried whether a fuller nested model would also work. It won't. And the corresponding full model would not work either. 
GLMM_rt_scr_nested <- glmer(word_rt ~ response_type/(valence*iscr_gng_resp_sqrt_z_score) + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_nested)                          # model failed to converge: This looks like a real convergence failure; maybe try simplifying your model?

GLMM_rt_scr_nested2 <- glmer(word_rt ~ response_type/(valence*iscr_gng_resp_sqrt_z_score) + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score + FH_minus_SH:iscr_gng_resp_sqrt_z_score + FA_minus_FH:iscr_gng_resp_sqrt_z_score + pos_minus_neg:iscr_gng_resp_sqrt_z_score + FH_minus_SH:pos_minus_neg:iscr_gng_resp_sqrt_z_score + FA_minus_FH:pos_minus_neg:iscr_gng_resp_sqrt_z_score  || subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_nested2)                          # model failed to converge: The estimates look good (and effects are like in the LMM_rt_scr_nested, but here FH - SH is also sign.) but I get the message:  This looks like a real convergence failure; maybe try simplifying your model?

GLMM_rt_scr_nested3 <- glmer(word_rt ~ response_type/(valence*iscr_gng_resp_sqrt_z_score) + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg + iscr_gng_resp_sqrt_z_score || subjectID) + (1 + FH_minus_SH + FA_minus_FH + iscr_gng_resp_sqrt_z_score || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_nested3)                          # model failed to converge: The estimates are strange, even though I get the message: The relative maximum gradient of 0.000109 is less than our 0.001 criterion. You can safely ignore any warnings about a claimed convergence failure.

GLMM_rt_scr_nested4 <- glmer(word_rt ~ response_type/(valence*iscr_gng_resp_sqrt_z_score) + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_nested4)                          # model failed to converge: The estimates are strange, even though I get the message: The relative maximum gradient of 0.000413 is less than our 0.001 criterion. You can safely ignore any warnings about a claimed convergence failure.

GLMM_rt_scr_nested5 <- glmer(word_rt ~ response_type/(valence*iscr_gng_resp_sqrt_z_score) + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_scr_nested5)                          # model failed to converge: The estimates look good (and all effects are like in the LMM_rt_scr_nested) but I get the message: The relative maximum gradient of 0.468 exceeds our 0.001 criterion. This looks like a real convergence failure; maybe try simplifying your model?
summary(rePCA(GLMM_rt_scr_nested5))                   # all items explain variance
print(VarCorr(GLMM_rt_scr_nested5),comp="Variance") 
didLmerConverge(GLMM_rt_scr_nested5) # The relative maximum gradient of 0.468 exceeds our 0.001 criterion. This looks like a real convergence failure; maybe try simplifying your model?

# I do not want to further reduce the model. The last one yields good estimates (in line with LMM_rt_scr_nested), but apparently does not converge. Other models seem to converge bur yield weird estimates. This is why I favor LMMs. 



###### STEP 10: CHECK ASSUMPTIONS 

# for GLMMs, testing assumptions is not straightforward, see: https://www.theanalysisfactor.com/regression-diagnostics-glmm/
# Assumption: The chosen link function is appropriate - given
# Assumption: Appropriate estimation of variance - irrelevant here
# Assumption: Random effects come from a normal distribution - seems ok (shapiro only sign for Intercept subjectID die to subject 28)
r_int_subjects <- ranef(GLMM_rt_scr_final)$subjectID$`(Intercept)`
qqnorm(r_int_subjects, main = "Q-Q Plot Random Intercept Subjects")
qqline(r_int_subjects)
shapiro_int_subjects <- shapiro.test(r_int_subjects)["p.value"]
#ranef(GLMM_rt_scr_final)$subjectID$`(Intercept)`: shapiro sign, bec. subject 28 is outlier

r_int_words <- ranef(GLMM_rt_scr_final)$word$`(Intercept)`
qqnorm(r_int_words, main = "Q-Q Plot Random Intercept Words")
qqline(r_int_words)
shapiro_int_words <- shapiro.test(r_int_words)["p.value"]

r_slope1_subjects <- ranef(GLMM_rt_scr_final)$subjectID$'FH_minus_SH:pos_minus_neg'
qqnorm(r_slope1_subjects, main = "Q-Q Plot Slope FH_minus_SH:pos_minus_neg")
qqline(r_slope1_subjects)
shapiro_slope1_subjects <- shapiro.test(r_slope1_subjects)["p.value"] 

r_slope2_subjects <- ranef(GLMM_rt_scr_final)$subjectID$'pos_minus_neg:FA_minus_FH'
qqnorm(r_slope2_subjects, main = "Q-Q Plot Slope FA_minus_FH:pos_minus_neg")
qqline(r_slope2_subjects)
shapiro_slope2_subjects <- shapiro.test(r_slope2_subjects)["p.value"] 





###################   Preparation for GNG   #################### 

###### STEP 1: PREPARE DATA 

# exclude GNG responses with misses or wrong keys, outliers in GNG or Word RT (consistantly excluded everywhere), incorrect responses and CI trials
data4mixedmodels_gng <- data4mixedmodels[data4mixedmodels$gng_resp <= 44 & data4mixedmodels$gng_invalid_rt == FALSE & data4mixedmodels$outlier_words == FALSE,]
# 11028 of 15480 trials left

# make categorical variables factors
data4mixedmodels_gng$response_type <- factor(data4mixedmodels_gng$response_type, levels=c("FH","FA","SH")) # comparing FA with FH and SH; comparison FH and SH makes not much sense, because these differ per definition
data4mixedmodels_gng$subjectID     <- factor(data4mixedmodels_gng$subjectID)




###### STEP 2: CHECk DISTRIBUTION 

# transformed rt -> almost normally distributed, also in the different conditions
plot(density(data4mixedmodels_gng$gng_rt_inverse))
ggplot(data4mixedmodels_gng, aes(x = gng_rt_inverse)) + geom_density() + facet_wrap(response_type ~ valence)
qqp(data4mixedmodels_gng$gng_rt_inverse, "norm")
#qqnorm(data4mixedmodels_gng$gng_rt_inverse)

skewness(data4mixedmodels_gng$gng_rt_inverse) # -0.6 -> is ok
kurtosis(data4mixedmodels_gng$gng_rt_inverse) #  3.4  -> is not ok


# raw rt 
plot(density(data4mixedmodels_gng$gng_rt))
ggplot(data4mixedmodels_gng, aes(x = gng_rt)) + geom_density() + facet_wrap(response_type ~ valence)

# raw rt -> fit gamma
gamma <- fitdistr(data4mixedmodels_gng$gng_rt, "gamma")
qqp(data4mixedmodels_gng$gng_rt, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

# raw rt -> fit inverse gaussian; inverse gaussian distribution fits quite well (for RTs, Lo and Anderson 2015 recommend gamma or inverse gaussian, Conny used Gamma)
inv_gauss <- nigFit(data4mixedmodels_gng$gng_rt, plots = FALSE, printOut = FALSE) 
par(mfrow = c(1,2))  # use par function in base graphics approach to arrange plots side by side here; grid.arrange only works for ggplot grahpics
plot(inv_gauss, which = 1, plotTitles = "Histogram of GNG RT inv. gaussian")
plot(inv_gauss, which = 3, plotTitles = "Q-Q Plot of GNG RT inv. gaussian")
par(mfrow = c(1, 1)) # reset "par" parameter



###### STEP 3: DEFINE CONTRASTS 

contrasts(data4mixedmodels_gng$response_type) <- contr.sdif(3)

# add contrast as numerical covariate via model matrix (in order to enter contrasts individually in model and hence to use ||)
mm_c    <- model.matrix( ~ response_type, data4mixedmodels_gng) # stick in model matrix 3 columns
glimpse(mm_c)                                                   # requires dplyr
glimpse(data4mixedmodels_gng)                                   # at the moment has 24 columns

# attach to dataframe 
data4mixedmodels_gng[,((ncol(data4mixedmodels_gng)+1):(ncol(data4mixedmodels_gng)+3))] <- mm_c
names(data4mixedmodels_gng)[(ncol(data4mixedmodels_gng)-2):ncol(data4mixedmodels_gng)] <- c("Grand Mean", "FA_minus_FH","SH_minus_FA")
glimpse(data4mixedmodels_gng)




###### STEP 4: SCALE CONTINUOUS PREDICTORS / COVARIATES 

# not included here



###################   Linear Mixed Models for GNG   #################### 

###### STEP 1 to 4: RUN CODE ABOVE (PREPARATION FOR GNG)

###### STEP 5: SPECIFY THE MAXIMAL MODEL 

LMM_rt_gng_max <- lmer(gng_rt_inverse ~ SH_minus_FA + FA_minus_FH + (1 + SH_minus_FA + FA_minus_FH  | subjectID), data=data4mixedmodels_gng, REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_gng_max)                         # model does converge -> identified model
summary(rePCA(LMM_rt_gng_max))                  # all items explained variance
print(VarCorr(LMM_rt_gng_max),comp="Variance")  



###### STEP 6: SIMPLIFICATION OF THE RANDOM EFFECTS STRUCTURE 

# ### now that we have an identified model, test non-significant variance components using likelihood ratio tests (starting with component explaning least variance) - Schad recommends drop 1 strategy, but he would rather only do rePCA
# 
# # test whether correlation between random slope and random intercept for subject improves fit
# LMM_rt_gng_red1 <- lmer(gng_rt_inverse ~ SH_minus_FA + FA_minus_FH + (1 + SH_minus_FA + FA_minus_FH  || subjectID), data=data4mixedmodels_gng, REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_gng_red1)                         # model does converge -> identified model
# anova(LMM_rt_gng_max,LMM_rt_gng_red1)            # correlation between random slope and random intercept for subject improves fit (a bit: AIC is lower, BIC is higher - BIC penalizes model complexity more) -> model LMM_rt_gng_max wins
# 
# # test whether random slope for response_type improves fit
# LMM_rt_gng_red2 <- lmer(gng_rt_inverse ~ SH_minus_FA + FA_minus_FH + (1 | subjectID), data=data4mixedmodels_gng, REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_gng_red2)                         # model does converge -> identified model
# anova(LMM_rt_gng_max,LMM_rt_gng_red2)            # random slope for response_type improves fit -> model LMM_rt_gng_max wins
# 
# # test whether random intercept for subject improves fit
# LMM_rt_gng_red3 <- lmer(gng_rt_inverse ~ SH_minus_FA + FA_minus_FH + (0 + SH_minus_FA + FA_minus_FH  | subjectID), data=data4mixedmodels_gng, REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_gng_red3)                         # model failed to converge
# anova(LMM_rt_gng_max,LMM_rt_gng_red3)            # random intercept for subject improves fit -> model LMM_rt_gng_max wins
# 
# # test whether reintroduction of correlation improves fit -> not necessary, because this is included already
# 
# 
# 
# ###### STEP 7: TEST OMNIBUS SIGNIFICANCE OF FIXED EFFECTS (based on model reduced by rePCA and LRT)
# 
# # test whether fixed effect of response_type improves fit
# LMM_rt_gng_red4 <- lmer(gng_rt_inverse ~ + (1 + SH_minus_FA + FA_minus_FH  | subjectID), data=data4mixedmodels_gng, REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_rt_gng_red4)                         # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
# anova(LMM_rt_gng_max,LMM_rt_gng_red4)            # fixed effect for response_type improves fit 
# 


###### STEP 8: PRESENT FINAL MODEL BY USING REML ESTIMATIONS for LMMs

# final model is
LMM_rt_gng_final <- lmer(gng_rt_inverse ~ SH_minus_FA + FA_minus_FH + (1 + SH_minus_FA + FA_minus_FH  | subjectID), data=data4mixedmodels_gng, REML = TRUE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_rt_gng_final)                         # model does converge -> identified model
# results: SH > FA, FH = FA


###### STEP 9: BREAK DOWN INTERAKTIONS BY USING NESTED MODEL 

# not necessary, because no interactions included



###### STEP 10: CHECK ASSUMPTIONS 

# LINEARITY & HOMOGENEITY OF VARIANCE
plot(fitted(LMM_rt_gng_final), residuals(LMM_rt_gng_final), xlab = "Fitted Values", ylab = "Residuals")  # plot residuals against the fitted values
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(LMM_rt_gng_final), residuals(LMM_rt_gng_final)))
# scatter looks not really like a random pattern ...?
# solid line covers the dashed line -> best-fit line fits well
# if the scatter is not random that means there's some other random or fixed effect that explains variation in your data that you're missing
# alternative: plot(LMM_rt_gng_final)


# NORMALITY OF RESIDUALS
qqPlot(resid(LMM_rt_gng_final)) 
qqnorm(resid(LMM_rt_gng_final))
qqline(resid(LMM_rt_gng_final)) # quite off at the extremes ...

skewness(resid(LMM_rt_gng_final)) # -0.5 -> is ok
kurtosis(resid(LMM_rt_gng_final)) #  6.3 -> is not ok




###################   Genralized Linear Mixed Models for GNG   #################### 

###### STEP 1 to 4: RUN CODE ABOVE (PREPARATION FOR GNG)

###### STEP 5: SPECIFY THE MAXIMAL MODEL 

GLMM_rt_gng_max <- glmer(gng_rt ~ SH_minus_FA + FA_minus_FH + (1 + SH_minus_FA + FA_minus_FH  | subjectID), data=data4mixedmodels_gng, family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(GLMM_rt_gng_max)                         # model failed to converge
summary(rePCA(GLMM_rt_gng_max))                  # all items explained variance
print(VarCorr(GLMM_rt_gng_max),comp="Variance")  

# the model still converges -> it would also be ok to use the maximal model
didLmerConverge(GLMM_rt_gng_max)

###### STEP 6: SIMPLIFICATION OF THE RANDOM EFFECTS STRUCTURE 

# start with zero-correlation parameter model by using ||
GLMM_rt_gng_red1 <- glmer(gng_rt ~ SH_minus_FA + FA_minus_FH + (1 + SH_minus_FA + FA_minus_FH  || subjectID), data=data4mixedmodels_gng, family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_gng_red1)                         # model does converge -> identified model

### now that we have an identified model, test non-significant variance components using likelihood ratio tests (starting with component explaning least variance) - Schad recommends drop 1 strategy, but he would rather only do rePCA

# # test whether random slope for response_type improves fit
# GLMM_rt_gng_red2 <- glmer(gng_rt ~ SH_minus_FA + FA_minus_FH + (1 | subjectID), data=data4mixedmodels_gng, family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_rt_gng_red2)                         # model does converge -> identified model
# anova(GLMM_rt_gng_red1,GLMM_rt_gng_red2)          # random slope for response_type improves fit -> model GLMM_rt_gng_red1 wins
# 
# # test whether random intercept for subject improves fit
# GLMM_rt_gng_red3 <- glmer(gng_rt ~ SH_minus_FA + FA_minus_FH + (0 + SH_minus_FA + FA_minus_FH  || subjectID), data=data4mixedmodels_gng, family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_rt_gng_red3)                         # model failed to converge
# anova(GLMM_rt_gng_red1,GLMM_rt_gng_red3)          # random intercept for subject improves fit -> model GLMM_rt_gng_red1 wins
# 
# 
# 
# 
# ###### STEP 7: TEST OMNIBUS SIGNIFICANCE OF FIXED EFFECTS (based on model reduced by rePCA and LRT)
# 
# # test whether fixed effect of response_type improves fit
# GLMM_rt_gng_red4 <- glmer(gng_rt ~ + (1 + SH_minus_FA + FA_minus_FH  || subjectID), data=data4mixedmodels_gng, family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
# summary(GLMM_rt_gng_red4)                         # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
# anova(GLMM_rt_gng_red1,GLMM_rt_gng_red4)          # fixed effect for response_type improves fit 
# 


###### STEP 8: PRESENT FINAL MODEL (BY USING REML ESTIMATIONS for LMMs)

# final model is
GLMM_rt_gng_final <- glmer(gng_rt ~ SH_minus_FA + FA_minus_FH + (1 + SH_minus_FA + FA_minus_FH  || subjectID), data=data4mixedmodels_gng, family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
summary(GLMM_rt_gng_final)                         # model does converge -> identified model
# results: SH > FA, FH = FA



###### STEP 9: BREAK DOWN INTERAKTIONS BY USING NESTED MODEL 

# not necessary, because no interactions included



###### STEP 10: CHECK ASSUMPTIONS 

# for GLMMs, testing assumptions is not straightforward, see: https://www.theanalysisfactor.com/regression-diagnostics-glmm/

# Assumption 1: The chosen link function is appropriate -> Is ok, see inverse gaussian distribution fit 
  
# Assumption 2: Random effects (intercepts and slopes) are normally distributed ->  Is only partly ok. 
                   # As can be seen in the Q-Q plots, there are 2 outliers that lead to deviation from normal distribution of random intercept and random slope SH_minus_FA. 
                   # Inspection of the random effects revealed that subject 29 (very high intercept and SH_minus_FA) and subject 26 (very high SH_minus_FA) are these outliers. Excluding these subjects may solve this problem. But as gng RT is of no real interest to me, I will not do that. 
                   # For the random intercept and random slope SH_minus_FA, the Shapiro-Wilk test is significant. 

# Assumption 2: Random effects (intercepts and slopes) are normally distributed - partly given
r_int<- ranef(GLMM_rt_gng_final)$subjectID$`(Intercept)`
qqnorm(r_int, main = "Q-Q Plot Random Intercept")
qqline(r_int)
shapiro_int <- shapiro.test(r_int)["p.value"]
#ranef(GLMM_rt_gng_final)$subjectID$`(Intercept)`: 1 subject (Ss 29) is an extreme outlier; rest is normally distributed, but Shapiro is sign.

r_slope1<- ranef(GLMM_rt_gng_final)$subjectID$SH_minus_FA
qqnorm(r_slope1, main = "Q-Q Plot Ran. Slope SH_minus_FA")
qqline(r_slope1)
shapiro_slope1 <- shapiro.test(r_slope1)["p.value"] 

r_slope2<- ranef(GLMM_rt_gng_final)$subjectID$FA_minus_FH
qqnorm(r_slope2, main = "Q-Q Plot Ran. Slope FA_minus_FH")
qqline(r_slope2)
shapiro_slope2 <- shapiro.test(r_slope2)["p.value"] 

# Assumption 3: Appropriate estimation of variance / No overdispersion - not relevant 
# Overdispersion: empirical variance in data exceeds nominal variance under some presumed model
# irrelevant for models that estimate a scale parameter (i.e. almost anything but Poisson or binomial: Gaussian, Gamma, negative binomial .)
# see how to proceed in case of overdispersion here: http://rstudio-pubs-static.s3.amazonaws.com/263877_d811720e434d47fb8430b8f0bb7f7da4.html












###################   Check Distribution for SCR    ####################    

###### STEP 2: CHECk DISTRIBUTION 

# raw scr -> not normally distributed
plot(density(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp), main = "Histogram SCR Raw")
ggplot(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], aes(x = iscr_gng_resp)) + geom_density() + facet_wrap(response_type ~ valence) + ggtitle("Histogram SCR Raw in Conditions")

# square root transformed scr -> not normally distributed (distribution is cut at 0)
plot(density(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp_sqrt), main = "Histogram Square Root Transformed SCR")
ggplot(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], aes(x = iscr_gng_resp_sqrt)) + geom_density() + facet_wrap(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$response_type) + ggtitle("Histogram Square Root Transformed SCR in Conditions")
qqp(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp_sqrt, "norm", main = "Q-Q Plot Square Root Transformed SCR", ylab = "sample quantiles")

# square root transformed z-standardized scr -> normally distributed, also in the different conditions
plot(density(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp_sqrt_z_score), main = "Histogram Square Root Transformed Z-Standardized SCR")
ggplot(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], aes(x = iscr_gng_resp_sqrt_z_score)) + geom_density() + facet_wrap(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$response_type) + ggtitle("Histogram Square Root Transformed Z-Standardized SCR in Conditions")
qqp(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp_sqrt_z_score, "norm", main = "Q-Q Plot Square Root Transformed Z-Standardized SCR", ylab = "sample quantiles")


### which distribution fits to raw scr data?

# overview - which distribution may fit?
descdist(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp, discrete = FALSE) # suggests gamma or beta

# try different distributions
normal <- fitdist(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp+1, "norm")
plot(normal)

gamma <- fitdist(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp+1, "gamma")
plot(gamma)

lognormal <- fitdist(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp+1, "lnorm")
plot(lognormal)

weibull <- fitdist(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp+1, "weibull")
plot(weibull)

inv_gauss <- nigFit(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp, plots = FALSE, printOut = FALSE) 
par(mfrow = c(1,2))  # use par function in base graphics approach to arrange plots side by side here; grid.arrange only works for ggplot grahpics
plot(inv_gauss, which = 1, plotTitles = "Histogram of Word RT inv. gaussian")
plot(inv_gauss, which = 3, plotTitles = "Q-Q Plot of Word RT inv. gaussian")
par(mfrow = c(1, 1)) # reset "par" parameter


normal$aic
gamma$aic   
lognormal$aic # AIC lowest for lognormal
weibull$aic

# distribution would be Box-Cox Power Exponential Distribution - not implemented in glmer (families ar gaussian, binomial, poisson, gamma, inverse.gaussian, abd quasi)
# library(gamlss)
# library(gamlss.dist)
# library(gamlss.add)
# fit <- fitDist(data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,]$iscr_gng_resp+1, k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)
# summary(fit)





###################   Linear Mixed Models for SCR ###################

###### STEP 1 to 4: RUN CODE ABOVE (PREPARATION FOR PRIMING AND SCR (AND FOR PURE SCR))

# I use df that was prepared for the priming-scr analysis; for pure SCR analysis, I actually would not need to exclude word outliers; but I do not want to exclude different trials for each analysis, so I keep it consistant here


###### STEP 5: SPECIFY THE MAXIMAL MODEL (HERE INCLUDING WORD VALENCE - see below for model without valence) 

LMM_scr_max <- lmer(iscr_gng_resp_sqrt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg | subjectID) + (1 + FH_minus_SH + FA_minus_FH | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_scr_max)                         # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits -> didLmerConverge showed that this is a real convergence failure
summary(rePCA(LMM_scr_max))  
print(VarCorr(LMM_scr_max),comp="Variance")  


###### STEP 6: SIMPLIFICATION OF THE RANDOM EFFECTS STRUCTURE 

# start with zero-correlation parameter model by using ||
LMM_scr_red1 <- lmer(iscr_gng_resp_sqrt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_scr_red1)                        # warning that model may be singular; post-fitting convergence checks not performed for (nearly-)singular fits 
summary(rePCA(LMM_scr_red1))                 # for word there are 2 items that do not explain variance, for subject ID one
print(VarCorr(LMM_scr_red1),comp="Variance") # for word it is intercept and FA_minus_FH; for subject ID it is FA_minus_FH:pos_minus_neg


# the model still converges -> it would also be ok to use the maximal model, as for the full model vs the reduced model all estimates are quite identical effect remain the same (see argument Singman and Kellen: lower-order effects should not be removed when higher-order effects are still included; they recommend to rather keep zero variance random effect, if the estimates are quite the same)
didLmerConverge(LMM_scr_red1)

# remove intercept and FA_minus_FH from words, FA_minus_FH:pos_minus_neg from subjectID (same result if this is done in subsequent steps)
LMM_scr_red2 <- lmer(iscr_gng_resp_sqrt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg || subjectID) + (0 + FH_minus_SH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_scr_red2)                        # model does converge, no singular fit -> model is identified
summary(rePCA(LMM_scr_red2))                 # all items explain variance
print(VarCorr(LMM_scr_red2),comp="Variance") 


### now that we have an identified model, test non-significant variance components using likelihood ratio tests (starting with component explaning least variance) - Schad recommends drop 1 strategy, but he would rather only do rePCA
# 
# # test whether random slope valence for subject improves fit (explains least variance)
# LMM_scr_red3 <- lmer(iscr_gng_resp_sqrt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + FH_minus_SH:pos_minus_neg || subjectID) + (0 + FH_minus_SH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_red3)            # model does converge, no singular fit -> model is identified
# anova(LMM_scr_red2,LMM_scr_red3) # no sign. difference, random slope valence for subject does not improve fit -> model 3 wins
# 
# # test whether random slope interaction reaction_type*valence for subjectID improves fit (explains least variance)
# LMM_scr_red4 <- lmer(iscr_gng_resp_sqrt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH || subjectID) + (0 + FH_minus_SH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_red4)            # model does converge, no singular fit -> model is identified
# anova(LMM_scr_red3,LMM_scr_red4) # no sign. difference, random slope interaction reaction_type*valence for subjectID does not improve fit -> model 4 wins
# 
# # test whether random slope reaction_type for word improves fit (explains least variance)
# LMM_scr_red5 <- lmer(iscr_gng_resp_sqrt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH || subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_red5)            # model does converge, no singular fit -> model is identified
# anova(LMM_scr_red4,LMM_scr_red5) # random slope reaction_type for word does not improve fit -> model 5 wins
# 
# # test whether random slope reaction_type for subject improves fit (explains least variance)
# LMM_scr_red6 <- lmer(iscr_gng_resp_sqrt ~ response_type*valence + (1 | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_red6)            # model does converge, no singular fit -> model is identified
# anova(LMM_scr_red5,LMM_scr_red6) # random slope reaction_type for subject does improve fit -> model 5 wins
# 
# # test whether random intercept for subject improves fit
# LMM_scr_red7 <- lmer(iscr_gng_resp_sqrt ~ response_type*valence + (0 + FH_minus_SH + FA_minus_FH || subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_red7)            # model does converge, no singular fit -> model is identified
# anova(LMM_scr_red5,LMM_scr_red7) # random intercept for subjectID improves fit -> model 5 wins
# 
# # test whether reintroduction of correlation improves fit
# LMM_scr_red8 <- lmer(iscr_gng_resp_sqrt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_red8)            # model does converge, no singular fit -> model is identified
# anova(LMM_scr_red5,LMM_scr_red8) # sign. difference (but small difference; AIC favors model 8)
# 
# 
# 
# ###### STEP 7: TEST OMNIBUS SIGNIFICANCE OF FIXED EFFECTS (based on model reduced by rePCA and LRT)
# 
# # test whether fixed effect of interaction response_type*valence improves fit
# LMM_scr_red9 <- lmer(iscr_gng_resp_sqrt ~ response_type+valence + (1 + FH_minus_SH + FA_minus_FH | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_red9)            # model does converge, no singular fit -> model is identified
# anova(LMM_scr_red8,LMM_scr_red9) # no sign. difference, fixed effect interaction reaction_type*valence does not improve fit 
# 
# # test whether fixed effect of valence improves fit
# LMM_scr_red10 <- lmer(iscr_gng_resp_sqrt ~ response_type + response_type:valence + (1 + FH_minus_SH + FA_minus_FH | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_red10)            # model does converge, no singular fit -> model is identified
# anova(LMM_scr_red8,LMM_scr_red10) # sign. difference but AIC identical; fixed effect of valence does not improve fit?
# 
# # test whether fixed effect of response_type improves fit
# LMM_scr_red11 <- lmer(iscr_gng_resp_sqrt ~ valence + response_type:valence + (1 + FH_minus_SH + FA_minus_FH | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_red11)            # model does converge, no singular fit -> model is identified
# anova(LMM_scr_red8,LMM_scr_red11) # no sign. difference; fixed effect of response type does not improve fit 




###### STEP 8: PRESENT FINAL MODEL BY USING REML ESTIMATIONS for LMMs

LMM_scr_final <- lmer(iscr_gng_resp_sqrt ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg || subjectID) + (0 + FH_minus_SH  || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = TRUE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_scr_final)            # model does converge, no singular fit -> model is identified
# I also computed it with iscr_gng_resp_sqrt_z_score and results are the same


###### STEP 9: BREAK DOWN INTERAKTIONS BY USING NESTED MODEL 

# not necessary, bec. no sign. interaction



###### STEP 10: CHECK ASSUMPTIONS 

# LINEARITY & HOMOGENEITY OF VARIANCE
plot(fitted(LMM_scr_final), residuals(LMM_scr_final), xlab = "Fitted Values", ylab = "Residuals")  # plot residuals against the fitted values
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(LMM_scr_final), residuals(LMM_scr_final)))
# scatter looks not like a random pattern (looks better when using the z score, but I actually do not want to do that)
# homogeneity of variance seems to be given but not linearity
# solid line covers the dashed line -> best-fit line fits well


# NORMALITY OF RESIDUALS
qqPlot(resid(LMM_scr_final)) 
qqnorm(resid(LMM_scr_final))
qqline(resid(LMM_scr_final)) # a bit off at the extremes, but that's often the case; again doesn't look too bad

skewness(resid(LMM_scr_final)) # 0.96 -> is ok
kurtosis(resid(LMM_scr_final)) # 1.85 -> is ok



####### Model without valence

# maximal model
LMM_scr_no_val_max <- lmer(iscr_gng_resp_sqrt ~ response_type + (1 + FH_minus_SH + FA_minus_FH  | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_scr_no_val_max)                         # model does converge, no singular fit -> model is identified
summary(rePCA(LMM_scr_no_val_max))                  # all items explain variance
print(VarCorr(LMM_scr_no_val_max),comp="Variance") 

# # test whether correlation random effects improves fit
# LMM_scr_no_val_red1 <- lmer(iscr_gng_resp_sqrt ~ response_type + (1 + FH_minus_SH + FA_minus_FH  || subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_no_val_red1)                        # model does converge, no singular fit -> model is identified
# anova(LMM_scr_no_val_max,LMM_scr_no_val_red1)       # sign. difference, AIC slightly smaller for max model
# 
# # test whether random slope response type improves fit
# LMM_scr_no_val_red2 <- lmer(iscr_gng_resp_sqrt ~ response_type + (1 | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_no_val_red2)                        # model does converge, no singular fit -> model is identified
# anova(LMM_scr_no_val_max,LMM_scr_no_val_red2)       # sign. difference, max model wins
# 
# # test whether random intercept improves fit
# LMM_scr_no_val_red3 <- lmer(iscr_gng_resp_sqrt ~ response_type + (0 + FH_minus_SH + FA_minus_FH  | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_no_val_red3)                        # model does converge, no singular fit -> model is identified
# anova(LMM_scr_no_val_max,LMM_scr_no_val_red3)       # sign. difference, max model wins
# 
# # test whether fixed effect response type improves fit
# LMM_scr_no_val_red4 <- lmer(iscr_gng_resp_sqrt ~ 1 + (1 + FH_minus_SH + FA_minus_FH  | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = FALSE, control=lmerControl(optimizer="bobyqa"))
# summary(LMM_scr_no_val_red4)                        # model does converge, no singular fit -> model is identified
# anova(LMM_scr_no_val_max,LMM_scr_no_val_red4)       # sign. difference

# the final model is 
LMM_scr_no_val_final <- lmer(iscr_gng_resp_sqrt ~ response_type + (1 + FH_minus_SH + FA_minus_FH  | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], REML = TRUE, control=lmerControl(optimizer="bobyqa"))
summary(LMM_scr_no_val_final)                         # model does converge, no singular fit -> model is identified
# this is also the final model I get when I initially include (1 + FH_minus_SH + FA_minus_FH  | word) as random effect in the maximal model (rePCA shows that for word all random effects including intercept explain zero variance) ...



# check assumptions
# LINEARITY & HOMOGENEITY OF VARIANCE
plot(fitted(LMM_scr_no_val_final), residuals(LMM_scr_no_val_final), xlab = "Fitted Values", ylab = "Residuals")  # plot residuals against the fitted values
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(LMM_scr_no_val_final), residuals(LMM_scr_no_val_final)))
# scatter looks not like a random pattern (looks better when using the z score, but I actually do not want to do that)
# homogeneity of variance seems to be given but not linearity
# solid line covers the dashed line -> best-fit line fits well


# NORMALITY OF RESIDUALS
qqPlot(resid(LMM_scr_no_val_final)) 
qqnorm(resid(LMM_scr_no_val_final))
qqline(resid(LMM_scr_no_val_final)) # a bit off at the extremes, but that's often the case; again doesn't look too bad

skewness(resid(LMM_scr_no_val_final)) # 0.96 -> is ok
kurtosis(resid(LMM_scr_no_val_final)) # 1.89 -> is ok




###################   Generalized Linear Mixed Model for SCR ###################

###### STEP 1 to 4: RUN CODE ABOVE (PREPARATION FOR SCR)

###### STEP 5: SPECIFY THE MAXIMAL MODEL (HERE INCLUDING WORD VALENCE) 

GLMM_scr_max <- glmer(iscr_gng_resp+1 ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg | subjectID) + (1 + FH_minus_SH + FA_minus_FH | word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(GLMM_scr_max)                         # warning: maxfun < 10 * length(par)^2 is not recommended. -> too many parameters to be estimated
summary(rePCA(GLMM_scr_max))  
print(VarCorr(GLMM_scr_max),comp="Variance")  



###### STEP 6: SIMPLIFICATION OF THE RANDOM EFFECTS STRUCTURE 

# start with zero-correlation parameter model by using ||
GLMM_scr_red1 <- glmer(iscr_gng_resp+1 ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (1 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(GLMM_scr_red1)                       # model fails to converge (according to didLmerConverge, I could ignore the warning, but the estimates and statistics are really strange)   
summary(rePCA(GLMM_scr_red1))                # all items explain variance
print(VarCorr(GLMM_scr_red1),comp="Variance") 


# remove intercept from words
GLMM_scr_red2 <- glmer(iscr_gng_resp+1 ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + pos_minus_neg + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (0 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(GLMM_scr_red2)                       # model fails to converge    

# remove random slope valence from subjectID
GLMM_scr_red3 <- glmer(iscr_gng_resp+1 ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID) + (0 + FH_minus_SH + FA_minus_FH || word), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(GLMM_scr_red3)                       # model fails to converge    

# remove random slope response type from word
GLMM_scr_red4 <- glmer(iscr_gng_resp+1 ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH + FH_minus_SH:pos_minus_neg + FA_minus_FH:pos_minus_neg || subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(GLMM_scr_red4)                       # model fails to converge   

# remove random slope interaction reaction_type*valence from subjectID 
GLMM_scr_red5 <- glmer(iscr_gng_resp+1 ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH  || subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(GLMM_scr_red5)                       # model does converge, no singular fit -> model is identified


# ### now that we have an identified model, test non-significant variance components using likelihood ratio tests (starting with component explaning least variance) - Schad recommends drop 1 strategy, but he would rather only do rePCA
# 
# # test whether random slope reaction_type for subjectID improves fit (explains least variance)
# GLMM_scr_red6 <- glmer(iscr_gng_resp+1 ~ response_type*valence + (1 | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_scr_red6)                       # model does converge, no singular fit -> model is identified
# anova(GLMM_scr_red5,GLMM_scr_red6)           # random slope reaction_type for subjectID does improve fit -> model 5 wins
# 
# # test whether random intercept for subject improves fit
# GLMM_scr_red7 <- glmer(iscr_gng_resp+1 ~ response_type*valence + (0 + FH_minus_SH + FA_minus_FH || subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_scr_red7)                       # model fails to converge
# anova(GLMM_scr_red5,GLMM_scr_red7)           # random intercept for subjectID improves fit -> model 5 wins

# test whether reintroduction of correlation improves fit
GLMM_scr_red8 <- glmer(iscr_gng_resp+1 ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(GLMM_scr_red8)                       # model does converge, no singular fit -> model is identified
anova(GLMM_scr_red5,GLMM_scr_red8)           # sign. difference (but small difference; AIC favors model 8, BIC model 5 -> go with model 8, because also identical to LMM_final_model)



# ###### STEP 7: TEST OMNIBUS SIGNIFICANCE OF FIXED EFFECTS (based on model reduced by rePCA and LRT)
# 
# # test whether fixed effect of interaction response_type*valence improves fit
# GLMM_scr_red9 <- glmer(iscr_gng_resp+1 ~ response_type + valence + (1 + FH_minus_SH + FA_minus_FH | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_scr_red9)                      # model does converge, no singular fit -> model is identified
# anova(GLMM_scr_red8,GLMM_scr_red9)          # no sign. difference, fixed effect interaction reaction_type*valence does not improve fit 
# 
# # test whether fixed effect of valence improves fit
# GLMM_scr_red10 <- glmer(iscr_gng_resp+1 ~ response_type + response_type:valence + (1 + FH_minus_SH + FA_minus_FH | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_scr_red10)                     # model does converge, no singular fit -> model is identified
# anova(GLMM_scr_red8,GLMM_scr_red10)         # sign. difference but AIC and BIC are identical
# 
# # test whether fixed effect of response_type improves fit
# GLMM_scr_red11 <- glmer(iscr_gng_resp+1 ~ valence + response_type:valence + (1 + FH_minus_SH + FA_minus_FH | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
# summary(GLMM_scr_red11)                     # model does converge, no singular fit -> model is identified
# anova(GLMM_scr_red8,GLMM_scr_red11)         # no sign. difference and AIC and BIC are identical



###### STEP 8: PRESENT FINAL MODEL (BY USING REML ESTIMATIONS for LMMs)

GLMM_scr_final <- glmer(iscr_gng_resp+1 ~ response_type*valence + (1 + FH_minus_SH + FA_minus_FH | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(GLMM_scr_final)                    # model does converge, no singular fit -> model is identified
# results the same as in LMM_scr_final



###### STEP 9: BREAK DOWN INTERAKTIONS BY USING NESTED MODEL 

# not necessary, bec. no sign, interaction



###### STEP 10: CHECK ASSUMPTIONS 

# for GLMMs, testing assumptions is not straightforward, see: https://www.theanalysisfactor.com/regression-diagnostics-glmm/
# Assumption: The chosen link function is appropriate - I hope so
# Assumption: Appropriate estimation of variance - irrelevant here

# Assumption: Random effects come from a normal distribution - not really given
r_int_subjects <- ranef(GLMM_scr_final)$subjectID$`(Intercept)`
qqnorm(r_int_subjects, main = "Q-Q Plot Random Intercept Subjects")
qqline(r_int_subjects)
shapiro_int_subjects <- shapiro.test(r_int_subjects)["p.value"]
#ranef(GLMM_scr_final)$subjectID$`(Intercept)`: 1 subject (Ss 16) is an extreme outlier; rest is normally distributed, but Shapiro is sign.

r_slope1_subjects <- ranef(GLMM_scr_final)$subjectID$FH_minus_SH
qqnorm(r_slope1_subjects, main = "Q-Q Plot Slope FH_minus_SH")
qqline(r_slope1_subjects)
shapiro_slope1_subjects <- shapiro.test(r_slope1_subjects)["p.value"] 

r_slope2_subjects <- ranef(GLMM_scr_final)$subjectID$FA_minus_FH
qqnorm(r_slope2_subjects, main = "Q-Q Plot Slope FA_minus_FH")
qqline(r_slope2_subjects)
shapiro_slope2_subjects <- shapiro.test(r_slope2_subjects)["p.value"] 
#ranef(GLMM_scr_final)$subjectID$FA_minus_FH: Shapiro is sign.





####### Model without valence

# maximal model
GLMM_scr_vo_val_final <- glmer(iscr_gng_resp+1 ~ response_type + (1 + FH_minus_SH + FA_minus_FH | subjectID), data=data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], family = inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(GLMM_scr_vo_val_final)                    # model does converge, no singular fit -> model is identified






###################   Transformation Decisions     ####################

data4mixedmodels_words     <- data4mixedmodels[data4mixedmodels$outlier_words == FALSE & data4mixedmodels$gng_resp <= 46 & data4mixedmodels$gng_invalid_rt == FALSE & data4mixedmodels$word_resp <= 54,]
data4mixedmodels_gng       <- data4mixedmodels[data4mixedmodels$gng_resp <= 44 & data4mixedmodels$gng_invalid_rt == FALSE,]
data4mixedmodels_words_scr <- data4mixedmodels[data4mixedmodels$outlier_words == FALSE & data4mixedmodels$gng_resp <= 44 & data4mixedmodels$gng_invalid_rt == FALSE & data4mixedmodels$word_resp <= 54 & data4mixedmodels$followed_or_preceded_by_FA_or_wrong_key == FALSE & !is.na(data4mixedmodels$iscr_gng_resp),]


### For RTs (word and gng)

# estimate optimal lambda word_rts: - 0.83
bc <- boxcox(word_rt ~ 1, data = data4mixedmodels_words[data4mixedmodels_words$word_accuracy==1,])
optlambda <- bc$x[which.max(bc$y)]

# confidence interval of optimal lambda is betw -0.79 and -0.86, so again arguing for inverse transformation (equals lambda -1; log transformation equals lambda = 0; https://robjhyndman.com/talks/RevolutionR/7-Transformations.pdf)
CI_optlambda <- range(bc$x[bc$y > max(bc$y)-qchisq(0.95,1)/2])


# estimate optimal lambda gng_rts: - 0.59, not far from -1, so go for inverse transformation here as well (dont do methods too complicated and these results to not really matter)
bc_gng <- boxcox(gng_rt ~ 1, data = data4mixedmodels_gng)
optlambda <- bc_gng$x[which.max(bc_gng$y)]

# confidence interval of optimal lambda is betw -0.54 and -0.63, so again arguing for inverse transformation (equals lambda -1)
CI_optlambda <- range(bc_gng$x[bc_gng$y > max(bc_gng$y)-qchisq(0.95,1)/2])



### For Accuracy

# estimate optimal lambda word accuracy: 2
bc <- boxcox(percent_correct_overall ~ 1, data = df4save,lambda = seq(-14, 14, 1/10))
optlambda <- bc$x[which.max(bc$y)]

# confidence interval of optimal lambda is betw 1.31 and 2 so again arguing for square transformation (equals lambda 2)
CI_optlambda <- range(bc$x[bc$y > max(bc$y)-qchisq(0.95,1)/2])



### For SCR

# estimate optimal lambda scr: -3.7 (default margins are -2/2 -> changed that) 
bc <- boxcox((iscr_gng_resp+1) ~ 1, data = data4mixedmodels_words_scr[data4mixedmodels_words_scr$word_accuracy==1,], lambda = seq(-6, 6, 1/10))
optlambda <- bc$x[which.max(bc$y)]

# confidence interval of optimal lambda is betw -3.8 and -3.6, so again arguing for transformation x = y^-4 (equals lambda -4)
CI_optlambda <- range(bc$x[bc$y > max(bc$y)-qchisq(0.95,1)/2])