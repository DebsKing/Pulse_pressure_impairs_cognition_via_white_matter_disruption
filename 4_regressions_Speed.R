### Multiple linear regression models of PSMD predicting Speed

# Here, a marker of white matter microstructure, the peak width of skeletonised
# mean diffusivity (PSMD), is related to processing speed. 

# Speed is a latent variable, produced in confirmatory factor analysis, 
# in script 2 of this analysis.

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
library(ggpubr) # for plotting muitple GGplots on one pane with ggarrange()

#-----------------------------------------------------------------------------#
### Load data

mydata <- readRDS("/.../2_mydata.rds")

# #-----------------------------------------------------------------#
# ### Notes on regression
# 
# # In my thesis chapter 1 (code at https:://github/DebsKing), 
# # I tested for the need to use robust regression, and found that I should use it. 
# # Here, I use robust regression for consistency. 
# # I also check the assumptions of regression on line ~330 below and find that
# # use of robust regression is justified. 

#-----------------------------------------------------------------------------#
### Prepare for regression

## Encode binary Sex as a factor
mydata$Sex <- as.factor(mydata$Sex) 
# head(mydata$Sex) # does have factor levels

### Z-score 
# Only z-score on participants in analysis, with complete data in vars of interest.
mydata_reg <- data.frame( na.omit(subset(mydata, 
                                         select=c(Speed, 
                                                  PSMD_k,
                                                  Age,Sex  
                                         ) )))
# n=611 for regressions

## Z-score
mydata_reg$Speed_z <- scale(mydata_reg$Speed, center=TRUE, scale=TRUE)
mydata_reg$PSMD_k_z <- scale(mydata_reg$PSMD_k, center=TRUE, scale=TRUE)
mydata_reg$Age_z <- scale(mydata_reg$Age, center=TRUE, scale=TRUE)

#-----------------------------------------------------------------#
### Multiple Linear Regression Models

# Models 1-5 are in script 3.

# Here, Models 6-8 are tested. 
# In the paper and thesis chapter, these models are labelled 2a-c. 

# 6. Without Age
lm_6 <- rlm(Speed_z ~ PSMD_k_z  + Sex,
            data=mydata_reg)
tab_model(lm_6,show.se=TRUE, digits=2)  



# 7. With Age linear
lm_7 <- rlm(Speed_z ~ PSMD_k_z  + Sex + Age_z ,
            data=mydata_reg)
tab_model(lm_7,show.se=TRUE, digits=2)  

AIC(lm_6, lm_7)
BIC(lm_6, lm_7)
# Model 7 wins on AIC and BIC
# so we take linear age into the next model:



# 8. With Age quadratic
lm_8 <- rlm(Speed_z ~ PSMD_k_z + scale(poly(Age_z,2))   + Sex,
            data=mydata_reg)
tab_model(lm_8,show.se=TRUE, digits=2)  

AIC(lm_7, lm_8)
BIC(lm_7, lm_8)
# Model 7 wins on AIC and BIC, 
# don't take Age quadratic into next model 

#-----------------------------------------------------------------#
### Save all models, for thesis write up

setwd("/.../Results_Regressions")
tab_model(lm_6,show.se=TRUE, digits=2, file="results_linearmodel6.xls") 
tab_model(lm_7,show.se=TRUE, digits=2, file="results_linearmodel7.xls") 
tab_model(lm_8,show.se=TRUE, digits=2, file="results_linearmodel8.xls") 

#-----------------------------------------------------------------#
### Save a table of model comparisons, for thesis write up

# first, extract fit metrics
aic <- AIC(lm_6,lm_7,lm_8)
bic <- BIC(lm_6,lm_7,lm_8)


## Create table to store results in
comparisons <- as.data.frame(matrix (nrow=2, ncol=3))
rownames(comparisons)<- c("Model.6v7",
                          "Model.7v8",
                          "Model.7v9")
colnames(comparisons) <- c("Difference in AIC",
                           "Difference in BIC",
                           "Difference in Sum of Squares")

### Calculations to report:

# Model 6 vs 7
comparisons$`Difference in AIC`[1] <- round((aic$AIC[1]-aic$AIC[2]),2)
comparisons$`Difference in BIC`[1] <- round((bic$BIC[1]-bic$BIC[2]),2)
temp <- round((anova(lm_6, lm_7)),2)
comparisons$`Difference in Sum of Squares`[1] <- temp$`Sum of Sq`[2]

# Model 7 vs 8
comparisons$`Difference in AIC`[2] <-round((aic$AIC[2]-aic$AIC[3]),2)
comparisons$`Difference in BIC`[2] <- round((bic$BIC[2]-bic$BIC[3]),2)
temp <- round((anova(lm_7, lm_8)),2)
comparisons$`Difference in Sum of Squares`[2] <- temp$`Sum of Sq`[2]

rm(temp)

# Save model comparisons
saveRDS(comparisons, file="linear_model_comparisons_SPEED.rds")

View(comparisons)

#-----------------------------------------------------------------#

### Ends


