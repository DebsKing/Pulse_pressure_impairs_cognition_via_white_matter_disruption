### Multi-group SEM

# The model SEM 1a did not include Age, and was:
# LVF2 > PSMD > Speed. 

# Here, SEM 1a is repeated such that the path estimates were allowed to vary, 
# across sub-groups of young, middle and old aged participants. 

# We run one model where paths can vary, across the 3 age groups. 
# And a second model where paths are constrained to be equal across the 3 age groups. 

# Models are compared. 
# If the model where paths are constrained wins, this shows that the strength
# of the mediation pathway does not vary with age. 

# Here, models are called 1 and 2. 
# These are reported in the write up as SEM 2a and 2b.

#-----------------------------------------------------------------------------#
### Notes

# In other SEM scripts of this analysis, I begin script by defining 
# the function: fit_sem(), given to me by Rik.
# This nicely stores and views in tables the fit metrics and the
# results. I don't need this for the multigroup comparison.
# All I need is anova on the two models, to see if the model where
# paths vary with age fits better.
# therefore I'm not using the fit_sem() defined function, 
# instead I'm using the standard lavaan sem(), and applying my
# bootstrap etc arguments to it. Results are the same.

#-----------------------------------------------------------------------------#

rm(list=ls()) # clear all

#-----------------------------------------------------------------------------#
### Load packages

library(base) # for scale(), to z-score
library(lavaan) # for SEM
library(semPlot) # for plotting SEM diagram.

#-----------------------------------------------------------------------------#
### Load data

mydata <- readRDS("/.../2_mydata.rds")

#-----------------------------------------------------------------------------#
### Z-score 

mydata <- data.frame( na.omit(subset(mydata, select=c(LVF2, PSMD_k, Speed,Sex, Stripe_index, Age) )))

mydata$Age_z <-scale(mydata$Age, scale=TRUE, center=TRUE)
mydata$LVF2_z <- scale(mydata$LVF2, center=TRUE, scale=TRUE)
mydata$PSMD_k_z <- scale(mydata$PSMD_k, center=TRUE, scale=TRUE)  # new PSMD data
mydata$Stripe_index_z <- scale(mydata$Stripe_index, center=TRUE, scale=TRUE) 
mydata$Speed_z <- scale(mydata$Speed, center=TRUE, scale=TRUE)

## Define quadratic terms
mydata$LVF2_q <- poly(mydata$LVF2,2)[,2]
mydata$LVF2_q_z <-scale(mydata$LVF2_q, scale=TRUE, center=TRUE)

#------------------------------------------------------------------------------#
### Split into age groups

# Rank order by Age 
attach(mydata)
mydata <- mydata[order(Age),]

# with speed -- only 570 people, n=190 per group equally.  

## Speed:
mydata$Age_Group <- "Young"
mydata[191:380,]$Age_Group  <- "Middle"
mydata[381:570,]$Age_Group  <- "Old"

mydata$Age_Group <- as.factor(mydata$Age_Group)
head(mydata$Age_Group)
library(plyr); count(mydata$Age_Group)

# this is bad coding:
range(mydata[1:190,]$Age) # 18-44
range(mydata[191:380,]$Age) # 44-65
range(mydata[381:570,]$Age) # 65-87

#------------------------------------------------------------------------------#
### SEM models

# Remove all constants / letters ("c1" "a1" etc), apart from the paths you do want to compare.  

# In Model 1, constrain all paths to be equal (i.e. nothing varies with Age)
# In Model 2, constrain all paths APART from LVF2 > PSMD, which can vary with Age. 

# Constrain paths by using: 
# c(a1, a1, a1) 

set.seed(2022) 
# for SEM imputations
# These SEM models are slow to fit, because fitting each with n=5,000 bootstrap imputations, over 3 age groups.


### Model 1 - no paths vary with age
mymodel_1 <- ' 
## Regressions 
    PSMD_k_z  ~ c(a1, a1, a1) * LVF2_z
    PSMD_k_z  ~ c(b1, b1, b1) * LVF2_q_z
    Speed_z ~ c(c1, c1, c1) * LVF2_z
    Speed_z ~ c(d1, d1, d1) * LVF2_q_z
    Speed_z  ~ c(e1, e1, e1) * PSMD_k_z

## Covariates
    PSMD_k_z  ~ Sex + Stripe_index_z
    Speed_z ~ Sex 
'
mymodel_1 <- sem(mymodel_1, data=mydata, group="Age_Group",
                 estimator='ml', # default, not robust but bootstrap is robust 
                 se = "bootstrap",  bootstrap = 5000, 
                 fixed.x=FALSE, # Rogier: this matters only for factor models -  it sets the first indicator to 1 by default/
                 meanstructure=FALSE ) # Rogier: this line says you don't care about means/intercepts - usually true in mediations)

setwd("/.../Results_SEM_multigroup__2023-10-03/")
saveRDS(mymodel_1, file="mymodel_1.rds")



### Model 2 - path of LVF2 to PSMD can vary across 3 age groups
mymodel_2 <- ' 
## Regressions 
    PSMD_k_z  ~  LVF2_z       ## remove constants: can vary with age.
    PSMD_k_z  ~ LVF2_q_z      ## remove constants: can vary with age.
    Speed_z ~ c(c1, c1, c1) * LVF2_z
    Speed_z ~ c(d1, d1, d1) * LVF2_q_z
    Speed_z  ~ c(e1, e1, e1) * PSMD_k_z

## Covariates
    PSMD_k_z  ~ Sex + Stripe_index_z
    Speed_z ~ Sex 
'
mymodel_2 <- sem(mymodel_2, data=mydata, group="Age_Group",
                 estimator='ml',se = "bootstrap",  bootstrap = 5000, 
                 fixed.x=FALSE, meanstructure=FALSE ) 

setwd("/.../Results_SEM_multigroup__2023-10-03//")
saveRDS(mymodel_2, file="mymodel_2.rds")


#-----------------------------------------------------------------------------#
### Compare models

anova(mymodel_1, mymodel_2)

# this comparison asks: 
# Does PSMD ~ LVF2 path need to vary with Age (Model 2), 
# when other paths are NOT allowed to vary with Age?
# If PSMD ~ LVF2 path does need to vary with Age, Model 2 will win, 
# over Model 1 where all paths are constrained.

## rename models to help anova interpretation:
mymodel_1_constrain_all <- mymodel_1
mymodel_2_PP_varies <- mymodel_2
anova(mymodel_1_constrain_all, mymodel_2_PP_varies)

#-----------------------------------------------------------------------------#

### ENDS