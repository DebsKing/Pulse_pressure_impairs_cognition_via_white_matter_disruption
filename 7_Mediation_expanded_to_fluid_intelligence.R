### Mediation Models, developed in the Framework of Structural Equation Modelling 


# This builds on the initial and main models in script 5:
# Pulse pressure > white matter > processing speed.

# These are here expanded to the cognitive ability discrepancy, 
# and to fluid intellignece, spearately. 

# Giving two separate dual mediation models: 

# Pulse pressure > white matter > speed > ability discrepancy 
# Pulse pressure > white matter > speed > fluid intelligence


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

mydata <- readRDS("/.../2_mydata.rds")

#-----------------------------------------------------------------#
### Z-score 

# Z-score only on subset of participants input to models

mydata <- data.frame( na.omit(subset(mydata, select=c(Age, 
                                                      LVF1, LVF2, LVF3, 
                                                      PSMD_k, Stripe_index, 
                                                      LCF2, DS, Speed, Sex  ) )))

mydata$Age_z <-scale(mydata$Age, scale=TRUE, center=TRUE)
mydata$LVF2_z <- scale(mydata$LVF2, center=TRUE, scale=TRUE)
mydata$PSMD_k_z <- scale(mydata$PSMD_k, center=TRUE, scale=TRUE)  # new PSMD data
mydata$Stripe_index_z <- scale(mydata$Stripe_index, center=TRUE, scale=TRUE)   
mydata$Speed_z <- scale(mydata$Speed, center=TRUE, scale=TRUE)
mydata$Fluid_z <- scale(mydata$LCF2, center=TRUE, scale=TRUE)
mydata$DS_z <- scale(mydata$DS, center=TRUE, scale=TRUE)

### Define quadratic terms
mydata$LVF2_q <- poly(mydata$LVF2,2)[,2]
mydata$LVF2_q_z <-scale(mydata$LVF2_q, scale=TRUE, center=TRUE)

#-----------------------------------------------------------------------------------#
### 3a. LVF2 to PSMD to Speed to ability discrepancy

mymodel_3a <- ' 
## Regressions
DS_z ~ c1 * LVF2_z
DS_z ~ c2 * LVF2_q_z

PSMD_k_z  ~ a1 * LVF2_z
PSMD_k_z  ~ a2 * LVF2_q_z

Speed_z  ~ b1 * PSMD_k_z
DS_z  ~ b2 * Speed_z

## Covariates
PSMD_k_z  ~ Sex + Stripe_index_z
Speed_z ~ Sex + Age_z
DS_z ~ Sex + Age_z

## Terms
Total := (abs(a1*b1*b2)) + (abs(a2*b1*b2)) + (abs(c1)) + (abs(c2))
Mediation_DS := (abs(a1*b1*b2) ) + (abs(a2*b1*b2) )
Proportion_DS := Mediation_DS / Total
'
mymodel_3a <- fit_sem(mymodel_3a, mydata)

setwd("/.../Results_SEM_2023-10-03/")
saveRDS(mymodel_3a, file="mymodel_3b.rds")

summary(mymodel_3a, fit.measures=TRUE) 


## compare 3a to zero-ing extra path
mymodel_3a_zero <- ' 
## Regressions
DS_z ~ c1 * LVF2_z
DS_z ~ c2 * LVF2_q_z

PSMD_k_z  ~ a1 * LVF2_z
PSMD_k_z  ~ a2 * LVF2_q_z

Speed_z  ~ b1 * PSMD_k_z
DS_z  ~ 0 * Speed_z

## Covariates
PSMD_k_z  ~ Sex + Stripe_index_z
Speed_z ~ Sex + Age_z
DS_z ~ Sex + Age_z 
'

mymodel_3a_zero <- fit_sem(mymodel_3a_zero, mydata)
anova(mymodel_3a, mymodel_3a_zero)
setwd("/.../Results_SEM_2023-10-03/")
saveRDS(mymodel_3a_zero, file="mymodel_3a_zero.rds")

#-----------------------------------------------------------------------------------#
### 3b. LVF2 to PSMD to Speed to fluid

# covary linear age onto Speed and Fluid

mymodel_3b <- ' 
## Regressions
Fluid_z ~ c1 * LVF2_z
Fluid_z ~ c2 * LVF2_q_z

PSMD_k_z  ~ a1 * LVF2_z
PSMD_k_z  ~ a2 * LVF2_q_z

Speed_z  ~ b1 * PSMD_k_z
Fluid_z  ~ b2 * Speed_z

## Covariates
PSMD_k_z  ~ Sex + Stripe_index_z
Speed_z ~ Sex + Age_z
Fluid_z ~ Sex + Age_z

## Terms
Total := (abs(a1*b1*b2)) + (abs(a2*b1*b2)) + (abs(c1)) + (abs(c2))
Mediation_PP := (abs(a1*b1*b2) ) + ( abs(a2*b1*b2))
ProportionPP := Mediation_PP / Total
Direct_LVF1 := (abs(c1)) + (abs(c2)) 

'
mymodel_3b <- fit_sem(mymodel_3b, mydata)

setwd("/.../Results_SEM_2023-10-03/")
saveRDS(mymodel_3b, file="mymodel_3b.rds")

summary(mymodel_3b, fit.measures=TRUE) 



## compare 3b to zero-ing extra path
mymodel_3b_zero <- ' 
## Regressions
Fluid_z ~ c1 * LVF2_z
Fluid_z ~ c2 * LVF2_q_z

PSMD_k_z  ~ a1 * LVF2_z
PSMD_k_z  ~ a2 * LVF2_q_z

Speed_z  ~ b1 * PSMD_k_z
Fluid_z  ~ 0 * Speed_z

## Covariates
PSMD_k_z  ~ Sex + Stripe_index_z
Speed_z ~ Sex + Age_z
Fluid_z ~ Sex + Age_z 
'
mymodel_3b_zero <- fit_sem(mymodel_3b_zero, mydata)
anova(mymodel_3b, mymodel_3b_zero)
setwd("/.../Results_SEM_2023-10-03/")
saveRDS(mymodel_3a_zero, file="mymodel_3b_zero.rds")



### reverse SEM 3b on fluid
mymodel_3b_rev <- ' 
## Regressions
Speed_z ~ c1 * LVF2_z
Speed_z ~ c2 * LVF2_q_z

PSMD_k_z  ~ a1 * LVF2_z
PSMD_k_z  ~ a2 * LVF2_q_z

Fluid_z  ~ b1 * PSMD_k_z
Speed_z  ~ b2 * Speed_z

## Covariates
PSMD_k_z  ~ Sex + Stripe_index_z
Speed_z ~ Sex + Age_z
Fluid_z ~ Sex + Age_z
'
mymodel_3b_rev <- fit_sem(mymodel_3b_rev, mydata)

setwd("/.../Results_SEM_2023-10-03/")
saveRDS(mymodel_3b_rev, file="mymodel_3b_rev.rds")


## This comaprison is necessary following watershed logic.
anova(mymodel_3b, mymodel_3b_rev)

#-----------------------------------------------------------------------------#
### ENDS