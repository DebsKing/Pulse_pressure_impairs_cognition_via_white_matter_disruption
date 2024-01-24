### Multiple linear regression models of vascular factors predicting PSMD

# Latent vascular factors (LVF) 1-3 were produced in a prior analysis. 
# https://www.nature.com/articles/s41598-022-27252-1. 
# https://github.com/DebsKing/Distinct_Vascular_Components_Relate_To_Cognition. 

# Here, we relate them to the marker of white matter microstructure, 
# the peak width of skeletonised mean diffusivity (PSMD).

#-----------------------------------------------------------------------------#

rm(list=ls()) # clear all

#-----------------------------------------------------------------------------#
### Load packages

## General packages
library(psych) # for describe(), to explore data.
library(base) # for scale(), to z-score
library(plyr)

## Packages for robust linear regression
library(MASS) # for robust linear regression with rlm()
library(sjPlot) # for HTML tables of regression models

# Visualisations
library(GGally) # for correlation plot
library(ggplot2)

#-----------------------------------------------------------------------------#
### Load data

mydata <- readRDS("/.../2_mydata.rds")

# #-----------------------------------------------------------------------------#
# ### Check variable distributions before running regressions
# 
# # Skew should be in the bounds of positive and negative 2. 
# # Kurtosis should be in the bounds of positive and negative 7. 
# # (Kline 2011; Aminu et al 2014 EJBM) 
# 
# # Note that deviations from normality in terms of skewness and kurtosis often do not 
# # make a substantive difference to analysis when the samples are n>200. 
# # (Tabachnick and Fidell, 2013). 
# 
# # References about skew and kurtosis are on:
# # https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/
# # This includes the reference:  Aminu et al., 2014, European Journal of Business and Management. 
# # This is a paper on normality, see section 4.4.3 for with further discussion & references on cut off values.
# 
# 
# ### Vascular factors 
# par(mfrow=c(3,2))
# qqnorm(mydata$LVF1); hist(mydata$LVF1, breaks=40)
# qqnorm(mydata$LVF2); hist(mydata$LVF2, breaks=40)
# qqnorm(mydata$LVF3); hist(mydata$LVF3, breaks=40)
# describe(mydata[,cbind("LVF1", "LVF2", "LVF3")])
# # all very comfortably within bounds of skew +-2, kurtosis +- 7. 
# # LVFs are peaky around 0, which we now from previous analysis and have
# # investigated; it is not problematic or due to an error.
# 
# ### Cognitive factors 
# qqnorm(mydata$Speed)
# hist(mydata$Speed, breaks=40)
# describe(mydata$Speed) # minimal skew and kurtosis.
# 
# ### PSMD
# hist(mydata$PSMD_k)
# describe(mydata$PSMD_k)
# # skew = 1.59, kurtosis = 3.59. 
# # Not ideal, but within acceptable bounds, 
# # especially when both in a large sample and using robust regression.
# # For comparison, I log transformed speed raw observations bcse of skew >2, and kurtosis of >7 and up to 42.
# # Comparatively, raw PSMD is ok, less visible outliers:
# par(mfrow=c(1,3))
# hist(mydata$PSMD_k, breaks=40)
# hist(as.numeric(mydata$RTchoice_RTsd_all), breaks=40)
# hist(as.numeric(mydata$RTsimple_RTmean), breaks=40)
# 
# #-----------------------------------------------------------------#
# ### Notes on regression
# 
# # In the previous analysis (code at https://github.com/DebsKing/Distinct_Vascular_Components_Relate_To_Cognition), 
# # I tested for the need to use robust regression, and found that I should use it. 
# # Here, I use robust regression for consistency. 
# # I also check the assumptions of regression on line ~300 below and find that
# # use of robust regression is justified. 

#-----------------------------------------------------------------------------#
### Prepare covariates specifically for regression

## Sex
mydata$Sex <- as.factor(mydata$Sex) 
# head(mydata$Sex) # does have factor levels

## Medications
# coding binary variables as factors ensures they are treated categorically.
# medication classes explained in prior project:
# https://github.com/DebsKing/Distinct_Vascular_Components_Relate_To_Cognition 
# https://www.nature.com/articles/s41598-022-27252-1 

mydata$AntiH_Binary <- as.factor(mydata$AntiH_Binary) 
# head(mydata$AntiH_Binary) # check it, and yes it does have factor levels

mydata$drugclass2BB <- as.factor(mydata$drugclass2BB) 
# head(mydata$drugclass2BB) 

mydata$drugclass2Diu <- as.factor(mydata$drugclass2Diu) 
# head(mydata$drugclass2Diu) 

mydata$drugclass2Stat <- as.factor(mydata$drugclass2Stat) 
# head(mydata$drugclass2Stat) 

#-----------------------------------------------------------------#
### Z-score before regression

# Only z-score on participants in analysis, with complete data in vars of interest.
# cannot do this directly within each rlm() command, because we don't z-score all, eg Sex

mydata_reg <- data.frame( na.omit(subset(mydata, 
                                     select=c(PSMD_k, 
                                              LVF1,LVF2, LVF3,
                                              Stripe_index, 
                                              Age,Sex ,
                                              AntiH_Binary, drugclass2BB, drugclass2Diu, drugclass2Stat 
                                     ) )))

# n=611 for regressions

## Z-score
mydata_reg$LVF1_z <- scale(mydata_reg$LVF1, center=TRUE, scale=TRUE)
mydata_reg$LVF2_z <- scale(mydata_reg$LVF2, center=TRUE, scale=TRUE)
mydata_reg$LVF3_z <- scale(mydata_reg$LVF3, center=TRUE, scale=TRUE)

mydata_reg$PSMD_k_z <- scale(mydata_reg$PSMD_k, center=TRUE, scale=TRUE)

mydata_reg$Stripe_index_z <- scale(mydata_reg$Stripe_index, center=TRUE, scale=TRUE)

mydata_reg$Age_z <-scale(mydata_reg$Age, scale=TRUE, center=TRUE)

# #-----------------------------------------------------------------#
# ### Check LVF2 quadratic is as we expect it to look
# 
# plot((poly(mydata_reg$LVF2,2))[,2])
# 
# mydata_reg$LVF2_z_q <- (poly(mydata_reg$LVF2_z,2))[,2]
# 
# ggplot(mydata_reg, aes(y=LVF2_z_q, x=Age)) +
#   geom_point(alpha=0.5)

#-----------------------------------------------------------------#
### Models

# Note: in paper and in PhD thesis, models 1-5 here are written up as Models 1a-e.

# 1. without Age
lm_1 <- rlm(PSMD_k_z ~ scale(poly(LVF1_z,2)) + 
                      scale(poly(LVF2_z,2)) + 
                        scale(poly(LVF3_z,2)) + 
                        Stripe_index_z +  Sex ,
                        data=mydata_reg)
tab_model(lm_1, show.se=TRUE) #LVF2 linear & quadratic are significant

# 2. With linear age
lm_2 <- rlm(PSMD_k_z ~ scale(poly(LVF1_z,2)) + 
              scale(poly(LVF2_z,2)) + 
              scale(poly(LVF3_z,2)) +  
              Stripe_index_z +  Sex + Age_z,
            data=mydata_reg)
tab_model(lm_2, show.se=TRUE) # LVF2 linear & quadratic are significant above Linear Age, Additive.
AIC(lm_1,lm_2)
BIC(lm_1,lm_2) # Model 2 with Linear Age, Additive wins on AIC & BIC.
# Take linear age into Model 3, with quadratic Age also.

# 3. With quadratic age
lm_3 <- rlm(PSMD_k_z ~ scale(poly(LVF1_z,2)) + 
              scale(poly(LVF2_z,2)) + 
              scale(poly(LVF3_z,2)) +   
              Stripe_index_z +  Sex + scale(poly(Age_z,2)),
            data=mydata_reg) 
tab_model(lm_3, show.se=TRUE) # LVF2 quadratic is significant above Quadratic Age. :)
AIC(lm_1,lm_2, lm_3) 
BIC(lm_1,lm_2,lm_3) # with quadratic age wins

### We are no longer testing the regression with interactions, as of 2023-07-20. 
### If we were to test and include it, it would make no difference bcse
### it loses model comparisons to the model without age interactions.  
#
# # 4. With quadratic age interactions
# lm_4 <- rlm(PSMD_k_z ~ scale(poly(LVF1_z,2)) * scale(poly(Age_z,2)) +
#               scale(poly(LVF2_z,2)) * scale(poly(Age_z,2))+ 
#               scale(poly(LVF3_z,2)) * scale(poly(Age_z,2)) +  
#               Stripe_index_z +  Sex ,
#             data=mydata_reg) 
# 
# AIC(lm_3, lm_4) # model 3 wins 
# BIC(lm_3, lm_4) # model 3 wins 

### 4. With medications. 
lm_4 <- rlm(PSMD_k_z ~
              # LVF1  * each medication
              scale( poly(LVF1_z,2)) *AntiH_Binary +
              scale(poly(LVF1_z,2)) *drugclass2BB +
              scale(poly(LVF1_z,2)) *drugclass2Diu +
              scale(poly(LVF1_z,2)) *drugclass2Stat +
              
              # LVF2 * each medication
              scale(poly(LVF2_z,2)) *AntiH_Binary +
              scale(poly(LVF2_z,2)) *drugclass2BB +
              scale(poly(LVF2_z,2)) *drugclass2Diu +
              scale(poly(LVF2_z,2)) *drugclass2Stat +
              
              # LVF3 * each medication
              scale( poly(LVF3_z,2)) *AntiH_Binary +
              scale(poly(LVF3_z,2)) *drugclass2BB +
              scale(poly(LVF3_z,2)) *drugclass2Diu +
              scale( poly(LVF3_z,2)) *drugclass2Stat +
              
              scale(poly(Age_z,2)) +
              
              Sex +Stripe_index_z,  data=mydata_reg)
AIC(lm_3, lm_4)
BIC(lm_3,lm_5)
# Similar on AIC, but Model 3 clearly wins on BIC. 
# Model 3 - without medications and without age interactions - wins.
# don't take medications into final Model 5. 

# 5. With sex interactions 
# This was motivated by Raes et al 2020, Hypertension.
lm_5 <- rlm(PSMD_k_z ~ scale(poly(LVF1_z,2))* Sex + 
              scale(poly(LVF2_z,2))* Sex + 
              scale(poly(LVF3_z,2))* Sex +   
              Stripe_index_z  + scale(poly(Age_z,2)),
            data=mydata_reg) 
tab_model(lm_5, show.se=TRUE) 
# LVF2 quadratic is significant above Quadratic Age
# and sex interactions:)

### Conclusion of model comparisons: 
AIC(lm_1,lm_2,lm_3,lm_4,lm_5) # Models 3 and 4 have lowest AIC
BIC(lm_1,lm_2,lm_3,lm_4,lm_5) # Model 3 has lowest BIC, Model 4 is not close
# so Model 3 wins.

#-----------------------------------------------------------------#
### Save all models, for write up

setwd("/.../Results_Regressions")
tab_model(lm_1,show.se=TRUE, digits=2, file="results_linearmodel1.xls") 
tab_model(lm_2,show.se=TRUE, digits=2, file="results_linearmodel2.xls") 
tab_model(lm_3,show.se=TRUE, digits=2, file="results_linearmodel3.xls") 
tab_model(lm_4,show.se=TRUE, digits=2, file="results_linearmodel4.xls") 
tab_model(lm_5,show.se=TRUE, digits=2, file="results_linearmodel5.xls") 

#-----------------------------------------------------------------#
### Save a table of model comparisons, for write up

# first, extract fit metrics
aic <- AIC(lm_1,lm_2,lm_3, lm_4, lm_5)
bic <- BIC(lm_1,lm_2,lm_3, lm_4, lm_5)


## Create table to store results in
comparisons <- as.data.frame(matrix (nrow=4, ncol=3))
rownames(comparisons)<- c("Model.1v2",
                          "Model.2v3",
                          "Model.3v4",
                          "Model.3v5")
colnames(comparisons) <- c("Difference in AIC",
                           "Difference in BIC",
                           "Difference in Sum of Squares")

### Calculations to report:

# Model 1 vs 2
comparisons$`Difference in AIC`[1] <- round((aic$AIC[1]-aic$AIC[2]),2)
comparisons$`Difference in BIC`[1] <- round((bic$BIC[1]-bic$BIC[2]),2)
temp <- round((anova(lm_1, lm_2)),2)
comparisons$`Difference in Sum of Squares`[1] <- temp$`Sum of Sq`[2]

# Model 2 vs 3
comparisons$`Difference in AIC`[2] <-round((aic$AIC[2]-aic$AIC[3]),2)
comparisons$`Difference in BIC`[2] <- round((bic$BIC[2]-bic$BIC[3]),2)
temp <- round((anova(lm_2, lm_3)),2)
comparisons$`Difference in Sum of Squares`[2] <- temp$`Sum of Sq`[2]

# Model 3 vs 4
comparisons$`Difference in AIC`[3] <-round((aic$AIC[3]-aic$AIC[4]),2)
comparisons$`Difference in BIC`[3] <- round((bic$BIC[3]-bic$BIC[4]),2)
temp <-round((anova(lm_3, lm_4)),2)
comparisons$`Difference in Sum of Squares`[3] <- temp$`Sum of Sq`[2]


# Model 3 vs 5
comparisons$`Difference in AIC`[4] <-round((aic$AIC[3]-aic$AIC[5]),2)
comparisons$`Difference in BIC`[4] <-round((bic$BIC[3]-bic$BIC[5]),2)
temp <-round((anova(lm_3, lm_5)),2)
comparisons$`Difference in Sum of Squares`[4] <- temp$`Sum of Sq`[2]

rm(temp)

# Save model comparisons
saveRDS(comparisons, file="linear_model_comparisons_PSMD.rds")

View(comparisons)

#------------------------------------------------------------------#
### Check assumptions of winning Model 3

# Plots to check assumptions usually run on non-robust regression, i.e. lm(). 
# We ran robust regression, i.e. rlm(). 
# The plots do not account for the 'robust' application of this function. 
# Therefore if these plots show that the assumptions are violated (eg by outliers),
# it motivates the use of robust regression.

## Residuals v fitted
# Checks the linearity of the data. Aim: no pattern & approx horizontal red line.
par(mfrow=c(1,1))
plot(lm_3,1) # linear. acceptable...but some points high up.

## QQ for normal distribution
# Aim: points follow dashed line
plot(lm_3,2) # outliers with points not following line

## Spread-location plot
# Checks homogeneity of variance. Aim: horizontal line & equally spread points,
# showing that the variance of residuals does not change in line with the predictors.
# Outliers are okay if they are the tail end of a normal distribution, but not if
# they have undue influence and differ from the pattern of the data.
plot(lm_3,3) 
# not a horizontal line, and spread of points (variance) increases from L to R
# this is called minor homoscedacity: 
# I tried log(PSMD) to fix this, and it is not a good solution.

## Residuals vs leverage plot.
# Aim: standardised residuals <3 absolute, or considered outliers.
# Aim: leverage <0.02 & no points in upper/lower R corners, or considered influential.
plot(lm_3,5) # Indicates possible influential outlier.

# This confirms that use of robust regression was appropriate, 
# and it is also consistent with my thesis Chapter 3 // Scientific Reports paper.

rm(aic, bic, lm_1, lm_2, lm_3, lm_4, lm_5, lm_6)

#-----------------------------------------------------------------------------#
## Investigating Effects in Model 3, and of LVF2 onto PSMD

# Reminder that LVF2 is the latent vascular factor expressing predominantly pulse pressure.
# And that PSMD is a marker of white matter microstructure.

# In regressions, we made LVF2^2 on the fly
# So here make it again, for use in the plots
mydata_reg$LVF2_z_q_z <- scale(poly(mydata_reg$LVF2_z,2)[,2], center=TRUE, scale=TRUE)

# Plot LVF2 linear against PSMD, with a quadratic line of fit. 
ggplot(mydata_reg, aes(y=PSMD_k_z, x=LVF2_z)) +
  geom_point(alpha=0.5)+
  stat_smooth(method="lm",se=FALSE, formula = y ~ x + I(x^2)) +
  xlab("Pulse Pressure Latent Factor (standardized)") +  
  ylab("PSMD (standardized)") +
  scale_x_continuous(breaks=c(-3,-2,-1,0,1,2,3)) +
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3,4,5)) +
  # make text bold:
  theme(axis.text=element_text(face="bold")) + 
  theme(axis.title=element_text(face="bold"))

# This plot was not reported. 

#-----------------------------------------------------------------------------#

### Ends





