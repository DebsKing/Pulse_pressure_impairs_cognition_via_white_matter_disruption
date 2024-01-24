### Mediation Models, developed in the Framework of Structural Equation Modelling 

# Models asking whether the relationship between the pulse pressure latent
# vascular factor (LVF2) and processing speed, is fully or partially
# accounted for through a marker of white matter microstructure, called
# the peak width of skeletonised mean diffusivity (PSMD).

#-----------------------------------------------------------------------------#
### Updates to analysis, following further discussions: 

# Not rejecting models based on chi-sq p<0.05, because it is biased in n>200. 
# Instead, assess fit on indices such as RMSEA for the first model. 
# This model is next expanded in subsequent models through the addition of
# extra measures (e.g. Age). Assess the benefit of each addition to the model
# by comparing the newly expanded model, to a version of itself where the 
# new path is constrained to be equal to zero. Only take the newly added path forwards
# into subsequent models if shown to benefit model fit.

# Taking quadratic LVFs into SEM. 
# The effect size and significant of mediation is reported as a combined total
# of linear and quadratic predictors. 

# Covarying Age on outcome variable only.
# Not moderating the 'a' path with Age (as presented at BNA and Rotman conferences).

#-----------------------------------------------------------------------------#

rm(list=ls()) # clear all

#-----------------------------------------------------------------------------#
### Load packages

library(base) # for scale(), to z-score
library(lavaan) # for SEM
library(semPlot) # for plotting SEM diagram.

#-----------------------------------------------------------------------------#
### Create function to run SEM

# This allows the user to quickly and iteratively expand model structures, below.

fit_sem <- function(mymodel, mydata){ 
  
  set.seed(2022)  
  modelfit <- sem(mymodel, 
                  data=mydata, 
                  estimator='ml', # default, not robust but bootstrap is robust 
                  se = "bootstrap", 
                  bootstrap = 5000,    
                  fixed.x=FALSE, # Rogier: this matters only for factor models -  it sets the first indicator to 1 by default/
                  meanstructure=FALSE # Rogier: this line says you don't care about means/intercepts - usually true in mediations
  )
  
  semPaths(modelfit, 
           what='est', 
           rotation = 2, # default rotation = 1 with four options
           curve = 2, # pull covariances' curves out a little
           nCharNodes = 4,
           nCharEdges = 5, # don't limit variable name lengths
           sizeMan = 8, # font size of manifest variable names
           style = "lisrel", # single-headed arrows vs. # "ram"'s 2-headed for variances
           edge.label.cex=0.8, 
           curvePivot = TRUE, 
           fade=FALSE, 
           residuals = FALSE, intercepts = FALSE,
           sizeInt = 7,sizeLat=5)
  
  fits <- round(fitMeasures(modelfit, c("cfi", "rmsea", "srmr")),digits=2)  # Aim: CFI>0.9; RMSEA<0.08; SRMR <0.08.
  View(fits)
  
  results_model <- as.data.frame (parameterEstimates (modelfit, ci=TRUE, standardized = TRUE) )
  is.num <- sapply (results_model, is.numeric)
  results_model[is.num] <- lapply (results_model[is.num], round,3)
  # keep only results data needed for results table in paper
  #results_model <- results_model[c("label", "est", "se", "z", "pvalue", "ci.lower", "ci.upper")]
  #colnames(results_model)<- c("Parameter", "Estimates", "SE", "Z", "P-value", "CI.lower", "CI.upper")
  View(results_model)
  
  write.csv(results_model, "/.../results_model.xls" , row.names=FALSE)
  modelfit
}

#-----------------------------------------------------------------------------#
### Load data

# Load the data output from Script 2, which created the Speed latent variable
mydata <- readRDS("/.../2_mydata.rds")

#-----------------------------------------------------------------#
### Z-score 

# Z-score only on subset of participants input to models

mydata <- data.frame( na.omit(subset(mydata, select=c(Age, 
                                                      LVF1, LVF2, LVF3, 
                                                      PSMD_k, Stripe_index, 
                                                      Speed, Sex  ) )))

mydata$Age_z <-scale(mydata$Age, scale=TRUE, center=TRUE)

mydata$LVF1_z <- scale(mydata$LVF1, center=TRUE, scale=TRUE)
mydata$LVF2_z <- scale(mydata$LVF2, center=TRUE, scale=TRUE)
mydata$LVF3_z <- scale(mydata$LVF3, center=TRUE, scale=TRUE)

mydata$PSMD_k_z <- scale(mydata$PSMD_k, center=TRUE, scale=TRUE)  # new PSMD data
mydata$Stripe_index_z <- scale(mydata$Stripe_index, center=TRUE, scale=TRUE)   
mydata$Speed_z <- scale(mydata$Speed, center=TRUE, scale=TRUE)

### Define quadratic terms

mydata$LVF1_q <- poly(mydata$LVF1,2)[,2]
mydata$LVF1_q_z <-scale(mydata$LVF1_q, scale=TRUE, center=TRUE)

mydata$LVF2_q <- poly(mydata$LVF2,2)[,2]
mydata$LVF2_q_z <-scale(mydata$LVF2_q, scale=TRUE, center=TRUE)

mydata$LVF3_q <- poly(mydata$LVF3,2)[,2]
mydata$LVF3_q_z <-scale(mydata$LVF3_q, scale=TRUE, center=TRUE)

mydata$Age_q <- poly(mydata$Age,2)[,2]
mydata$Age_q_z <-scale(mydata$Age_q, scale=TRUE, center=TRUE)

#-----------------------------------------------------------------------------#
### SEM models in set 1

### 1a. LVF2 to PSMD to Speed
mymodel_1a <- ' 
## Regressions
    Speed_z ~ c1 * LVF2_z
    Speed_z ~ c2 * LVF2_q_z


    PSMD_k_z  ~ a1 * LVF2_z
    PSMD_k_z  ~ a2 * LVF2_q_z

    Speed_z  ~ b1 * PSMD_k_z

## Covariates
    PSMD_k_z  ~ Sex + Stripe_index_z
    Speed_z ~ Sex  

## Terms
    Total := (abs(a1*b1)) + (abs(a2*b1)) + (abs(c1)) + (abs(c2))
    Mediation_PP := (abs(a1*b1)) + (abs(a2*b1))
    Proportion_PP := Mediation_PP / Total
    Direct_effects := (abs(c1)) + (abs(c2))  # ** newly added.
'
mymodel_1a <- fit_sem(mymodel_1a, mydata)

summary(mymodel_1a, fit.measures=TRUE)  

setwd("/.../Results_SEM_2023-10-03/")
saveRDS(mymodel_1a, file="mymodel_1a.rds")


# ### To comapre to a version where PP and PSMD are reversed:
mymodel_1a_reverse <- '
## Regressions
    Speed_z ~ c1 * PSMD_k_z

    LVF2_z    ~ a1 * PSMD_k_z
    LVF2_q_z  ~ a2 * PSMD_k_z

    Speed_z  ~ b1 * LVF2_z
    Speed_z  ~ b1 * LVF2_q_z

## Covariates
    PSMD_k_z  ~ Sex + Stripe_index_z
    Speed_z ~ Sex
'
mymodel_1a_reverse <- fit_sem(mymodel_1a_reverse, mydata)

setwd("/.../Results_SEM_2023-10-03/")
saveRDS(mymodel_1a_reverse, file="mymodel_1a_reverse.rds")

anova(mymodel_1a, mymodel_1a_reverse)

AIC(mymodel_1a,  mymodel_1a_reverse)
BIC(mymodel_1a, mymodel_1a_reverse)


#-----------------------------------------------------------------------------------#
### SEM 1b - add linear age on dependent variable only

### 1b. LVF2 to PSMD to Speed
mymodel_1b <- ' 
## Regressions
Speed_z ~ c1 * LVF2_z
Speed_z ~ c2 * LVF2_q_z

PSMD_k_z  ~ a1 * LVF2_z
PSMD_k_z  ~ a2 * LVF2_q_z

Speed_z  ~ b1 * PSMD_k_z

## Covariates
PSMD_k_z  ~ Sex + Stripe_index_z
Speed_z ~ Sex  + Age_z

## Terms
    Total := (abs(a1*b1)) + (abs(a2*b1)) + (abs(c1)) + (abs(c2))
    Mediation_PP := (abs(a1*b1)) + (abs(a2*b1))
    Proportion_PP := Mediation_PP / Total
    Direct_effects := (abs(c1)) + (abs(c2))  # ** newly added.


'
mymodel_1b <- fit_sem(mymodel_1b, mydata)

summary(mymodel_1b, fit.measures=TRUE) 

setwd("/.../Results_SEM_2023-10-03/")
saveRDS(mymodel_1b, file="mymodel_1b.rds")



### Compare SEM 1b with Age to a version where the path connecting Age is zero-ed.
mymodel_1b_zero <- ' 
## Regressions
Speed_z ~ c1 * LVF2_z
Speed_z ~ c2 * LVF2_q_z

PSMD_k_z  ~ a1 * LVF2_z
PSMD_k_z  ~ a2 * LVF2_q_z

Speed_z  ~ b1 * PSMD_k_z

## Covariates
PSMD_k_z  ~ Sex + Stripe_index_z
Speed_z ~ Sex  + 0*Age_z

'
mymodel_1b_zero <- fit_sem(mymodel_1b_zero, mydata)
saveRDS(mymodel_1b_zero, file="mymodel_1b_zero.rds")

anova(mymodel_1b, mymodel_1b_zero) 

#-----------------------------------------------------------------------------------#
### SEM 1c - add quadratic Age

### 1c. LVF2 to PSMD to Speed
mymodel_1c <- ' 
## Regressions
Speed_z ~ c1 * LVF2_z
Speed_z ~ c2 * LVF2_q_z

PSMD_k_z  ~ a1 * LVF2_z
PSMD_k_z  ~ a2 * LVF2_q_z

Speed_z  ~ b1 * PSMD_k_z

## Covariates
PSMD_k_z  ~ Sex + Stripe_index_z
Speed_z ~ Sex  + Age_z + Age_q_z

## Terms
Total := (abs(a1*b1)) + (abs(a2*b1)) + (abs(c1)) + (abs(c2))
Mediation_PP := (abs(a1*b1)) + (abs(a2*b1))
Proportion_PP := Mediation_PP / Total
    Direct_effects := (abs(c1)) + (abs(c2))  # ** newly added.

'
mymodel_1c <- fit_sem(mymodel_1c, mydata)

setwd("/.../Results_SEM_2023-10-03/")
saveRDS(mymodel_1c, file="mymodel_1c.rds")

summary(mymodel_1c, fit.measures=TRUE) 



### Compare SEM 1c with Age to a version where Age is zero-ed.
mymodel_1c_zero <- ' 
## Regressions
Speed_z ~ c1 * LVF2_z
Speed_z ~ c2 * LVF2_q_z

PSMD_k_z  ~ a1 * LVF2_z
PSMD_k_z  ~ a2 * LVF2_q_z

Speed_z  ~ b1 * PSMD_k_z

## Covariates
PSMD_k_z  ~ Sex + Stripe_index_z
Speed_z ~ Sex  + Age_z + 0*Age_q_z

'
mymodel_1c_zero <- fit_sem(mymodel_1c_zero, mydata)
saveRDS(mymodel_1c_zero, file="mymodel_1c_zero.rds")

anova(mymodel_1c, mymodel_1c_zero) 


#-----------------------------------------------------------------------------------#
### SEM 1d - add all LVFs, covary Age linear only

# This is a sensitivity analysis, to see if the effects of pulse pressure (LVF2)
# on white matter microstructure (PSMD) and processing speed, hold after
# accounting for other latent vascular factors, of BP and HRV. 
# Here, latent vascular factors (LVFs) are named LVF1-3.
# Formation of all LVFs is outlined in the prior paper https://www.nature.com/articles/s41598-022-27252-1 

mymodel_1d <- ' 
## Regressions
Speed_z ~ c1 * LVF1_z
Speed_z ~ c2 * LVF1_q_z

Speed_z ~ c3 * LVF2_z
Speed_z ~ c4 * LVF2_q_z

Speed_z ~ c5 * LVF3_z
Speed_z ~ c6 * LVF3_q_z

PSMD_k_z  ~ a1 * LVF1_z
PSMD_k_z  ~ a2 * LVF1_q_z

PSMD_k_z  ~ a3 * LVF2_z
PSMD_k_z  ~ a4 * LVF2_q_z

PSMD_k_z  ~ a5 * LVF3_z
PSMD_k_z  ~ a6 * LVF3_q_z

Speed_z  ~ b1 * PSMD_k_z

## Covariates
PSMD_k_z  ~ Sex + Stripe_index_z  
Speed_z ~ Sex  + Age_z

## Terms
Total := (abs(a1*b1)) + (abs(a2*b1)) +  # LVF1
(abs(a3*b1)) + (abs(a4*b1)) +  # LVF2
(abs(a5*b1)) + (abs(a6*b1)) +  # LVF3
(abs(c1)) + (abs(c2)) + # LVF1
(abs(c3)) + (abs(c4)) + # LVF2
(abs(c5)) + (abs(c6))  # LVF3

Mediation_LVF1 := abs(a1*b1) +  abs(a2*b1) # combining linear and quadratic
Mediation_LVF2 := abs(a3*b1) + abs(a4*b1)
Mediation_LVF3 := abs(a5*b1) + abs(a6*b1)

Proportion_LVF1 := Mediation_LVF1 / Total
Proportion_m_LVF2 := Mediation_LVF2 / Total
Proportion_m_LVF3 := Mediation_LVF3 / Total

Direct_LVF1 := (abs(c1)) + (abs(c2)) 
Direct_LVF2 := (abs(c3)) + (abs(c4))
Direct_LVF3 := (abs(c5)) + (abs(c6)) # ** newly added.


'
mymodel_1d <- fit_sem(mymodel_1d, mydata)

setwd("/.../Results_SEM_2023-10-03/")
saveRDS(mymodel_1d, file="mymodel_1d.rds")

#-----------------------------------------------------------------------------------#

# The model structure was next tested in sub-groups split by age. 
# In script 6. 

### ENDS 
