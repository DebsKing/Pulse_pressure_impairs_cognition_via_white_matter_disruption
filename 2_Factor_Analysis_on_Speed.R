### CFA on cognitive variables 

#-----------------------------------------------------------------------------#
###  LOAD

rm(list=ls()) # clear all

library(corrplot)
library(lavaan)
library(plyr) # count() function
library(psych) 
library(ggplot2)
library(GGally)
library(semPlot)

mydata <- readRDS("/.../1_mydata.rds")

#------------------------------------------------------------------------------#
### Prepare reaction time variables: recode data entry errors 

# One participant has calculated standard deviation = zero for simple RT, 
# and a second participant has the same for choice RT 
# SD = 0. Indicates that all values are the same.
# It seems improbable that reaction times would exactly be the same, across the 50 trials in each task.
# Therefore it seems improbable that standard deviation is zero. And this suggests a data entry error.
# I therefore re-code the two instances of SD=0 to be NaN.

## Recode SD=0 as NaN. Their data will be imputed in CFA below. 
mydata$RTsimple_RTsd [mydata$RTsimple_RTsd == "0"] <- "NaN"
mydata$RTchoice_RTsd_all [mydata$RTchoice_RTsd_all == "0"] <- "NaN"

#------------------------------------------------------------------------------------------------#
### Only run CFA on participants with =>2 cognitive observations

# In order to quantify the number of cognitive observations per participant, 
# I first record whether each cognitive varaible is missing or present, for each participant.
# 1 = present;  0 = missing data. Then re-format it to be in numeric format. 
# This data of present/missing is encoded in a "mask" variable.

mydata$RTsimple_RTtrim3mean_mask [mydata$RTsimple_RTtrim3mean != "NaN"] <- "1"
mydata$RTsimple_RTtrim3mean_mask [mydata$RTsimple_RTtrim3mean == "NaN"] <- "0"
mydata$RTsimple_RTtrim3mean_mask <- as.numeric(unlist(mydata$RTsimple_RTtrim3mean_mask))

mydata$RTsimple_RTsd_mask [mydata$RTsimple_RTsd != "NaN"] <- "1"
mydata$RTsimple_RTsd_mask [mydata$RTsimple_RTsd == "NaN"] <- "0"
mydata$RTsimple_RTsd_mask <- as.numeric(unlist(mydata$RTsimple_RTsd_mask))

mydata$RTchoice_RTtrim3mean_all_mask [mydata$RTchoice_RTtrim3mean_all != "NaN"] <- "1"
mydata$RTchoice_RTtrim3mean_all_mask [mydata$RTchoice_RTtrim3mean_all == "NaN"] <- "0"
mydata$RTchoice_RTtrim3mean_all_mask <- as.numeric(unlist(mydata$RTchoice_RTtrim3mean_all_mask))

mydata$RTchoice_RTsd_all_mask [mydata$RTchoice_RTsd_all != "NaN"] <- "1"
mydata$RTchoice_RTsd_all_mask [mydata$RTchoice_RTsd_all == "NaN"] <- "0"
mydata$RTchoice_RTsd_all_mask <- as.numeric(unlist(mydata$RTchoice_RTsd_all_mask))

### Quantify complete data per participant
# This identifie participants with <2 observations, across speed observations
# Loop calculates the total number of variables recorded per participant
for (i in 1:708) {
  mydata$Observations_speed [i] <- (sum(
    # as.numeric(
    mydata$RTsimple_RTtrim3mean_mask[i], mydata$RTsimple_RTsd_mask[i],
    mydata$RTchoice_RTtrim3mean_all_mask[i], mydata$RTchoice_RTsd_all_mask[i]
  ))
}

# Across the 4 speed measures, what is the range of total observations per participant?
range(as.numeric(mydata$Observations_speed)) # 0-4

sum(mydata$Observations_speed <2 ) 
# 44 participants have 0 observations, across the 4 measures of speed; 
# these participants will be automatically omitted from CFA. 

plot(mydata$Observations_speed)
# zero participants have 1 observations. 
sum(mydata$Observations_speed >=2 )
# n=664 to be included in CFA. 


#------------------------------------------------------------------------------#
### Check distributions

# Acceptable bounds of +-2 skew, and +-7 kurtosis (Kline 2011; Aminu et al 2014 EJBM) .
# https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/. 


### RT simple SD -- distribution
mydata$RTsimple_RTsd <- as.numeric(mydata$RTsimple_RTsd) # unsure why this necessary
plot(mydata$RTsimple_RTsd); 
qqnorm((mydata$RTsimple_RTsd))
hist(((mydata$RTsimple_RTsd)), breaks=40)
describe(mydata[,cbind("RTsimple_RTsd")]) 
# skew =2.1, kurtosis= 6.0. 
# outside bounds of acceptable skew and kurtosis. 
# range doesn't cross 0 and is continuous; log transform. 


# RT simple mean -- distribution
plot(mydata$RTsimple_RTtrim3mean) 
qqnorm(mydata$RTsimple_RTtrim3mean)
hist(mydata$RTsimple_RTtrim3mean, breaks=40)
describe(mydata[,cbind("RTsimple_RTtrim3mean")]) 
# skew = 1.92, kurtosis =6.97. 
# outside bounds of acceptable skew and kurtosis. 
# range doesn't cross 0 and is continuous; log transform. 


# RT choice SD 
mydata$RTchoice_RTsd_all <- as.numeric(mydata$RTchoice_RTsd_all)
plot(mydata$RTchoice_RTsd_all) 
qqnorm(mydata$RTchoice_RTsd_all)
hist(mydata$RTchoice_RTsd_all, breaks=40)
describe(mydata[,cbind("RTchoice_RTsd_all")]) 
# skew = 4.13, kurtosis =42.29 
# outside bounds of acceptable skew and kurtosis. 
# range doesn't cross 0 and is continuous; log transform. 


### RT choice mean 
plot(mydata$RTchoice_RTtrim3mean_all); 
qqnorm((mydata$RTchoice_RTtrim3mean_all))
hist(((mydata$RTchoice_RTtrim3mean_all)), breaks=40)
describe(mydata[,cbind("RTchoice_RTtrim3mean_all")]) 
# skew =1.45, kurtosis=3.3 
# range doesn't cross 0.
# reasonable distribution, but log transform as consistent with other RT variables.


# Now apply: log(x)
mydata$RTsimple_RTtrim3mean_log <- log(mydata$RTsimple_RTtrim3mean )
mydata$RTsimple_RTsd_log <- log(as.numeric(mydata$RTsimple_RTsd)) # initially "as character", so put as numeric.
mydata$RTchoice_RTtrim3mean_all_log <- log(as.numeric(mydata$RTchoice_RTtrim3mean_all)) 
mydata$RTchoice_RTsd_all_log <- log(as.numeric(mydata$RTchoice_RTsd_all)) 


# check if distributions improved after log(x)
par(mfrow=c(2,2))
hist(mydata$RTsimple_RTtrim3mean_log, breaks=40) 
hist(mydata$RTsimple_RTsd_log, breaks=40)
hist(mydata$RTchoice_RTtrim3mean_all_log, breaks=40) 
hist(mydata$RTchoice_RTsd_all_log, breaks=40)

describe(mydata$RTsimple_RTtrim3mean_log)
describe(mydata$RTsimple_RTsd_log)
describe(mydata$RTchoice_RTtrim3mean_all_log)
describe(mydata$RTchoice_RTsd_all_log)

# After log transform, distributions are improved and within bounds of skew (+-2), kurtosis (+-7).

## Invert reaction times
# Invert to follow expected trajectory of decreases in performance with age
mydata$RTsimple_RTtrim3mean_log_inv <- mydata$RTsimple_RTtrim3mean_log * -1
mydata$RTsimple_RTsd_log_inv  <- mydata$RTsimple_RTsd_log * -1
mydata$RTchoice_RTtrim3mean_all_log_inv <- mydata$RTchoice_RTtrim3mean_all_log * -1
mydata$RTchoice_RTsd_all_log_inv <- mydata$RTchoice_RTsd_all_log * -1


# # Check inversion shows a performance decrease with age
par(mfrow=c(2,2))

plot(mydata$Age, mydata$RTsimple_RTtrim3mean_log_inv, 
     ylab="RT Simple mean (inverted log)")

plot(mydata$Age, mydata$RTsimple_RTsd_log_inv, 
     ylab="RT Simple SD (inverted log)")

plot(mydata$Age, mydata$RTchoice_RTtrim3mean_all_log_inv, 
     ylab="RT Choice mean (inverted log)")

plot(mydata$Age, mydata$RTchoice_RTsd_all_log_inv, 
     ylab="RT Choice SD (inverted log)")

#------------------------------------------------------------------------------------------------#
### For methods of paper, count observations in each measure

sum(mydata$RTsimple_RTtrim3mean_log_inv != "NaN") # n = 658
sum(mydata$RTsimple_RTsd_log_inv != "NaN")# n = 658

sum(mydata$RTchoice_RTtrim3mean_all_log_inv != "NaN") # n=656
sum(mydata$RTchoice_RTsd_all_log_inv != "NaN") # n = 656

#------------------------------------------------------------------------------#
###  correlation plot of cognitive factors input to CFA
 
# Figure B.1 in Appendix of Chapter 4 in PhD thesis. 
# Figure 1 in Appendix of paper.
 
 mydata_cor <-  ( data.frame(mydata$Age,
                             mydata$RTsimple_RTtrim3mean_log_inv,mydata$RTsimple_RTsd_log_inv,
                             mydata$RTchoice_RTtrim3mean_all_log_inv, mydata$RTchoice_RTsd_all_log_inv
 ) )
 
 colnames(mydata_cor) <- c("Age", 
                           "Simple_Mean", "Simple_SD", 
                           "Choice_Mean", "Choice_SD")
 par(mfrow=c(1,1))
 
 my_fn <- function (data, mapping ) {
   p <- ggplot(data=mydata_cor,mapping=mapping)+
     geom_jitter(alpha=0.1) +
     geom_smooth (method=lm, color = "blue" ,se=FALSE) 
   
   p
 }
 
 g = ggpairs(mydata_cor,
             lower=list(continuous=my_fn) )

 g
 
 myplot_cor <- g
 
 setwd("/.../")
 ggsave("myplot_cor_speed_cfa_inputs.png")
 
#------------------------------------------------------------------------------------------------#
### CFA

# Note that in coding cfa(), the additional arguments to standardize and estimate missing
# data produce in the results print out "Intercepts" which can have negative values. This 
# is not to report. We do however care if the "Variances" are negative.

model_CF_1 <- 'Speed =~ RTsimple_RTtrim3mean_log_inv + RTsimple_RTsd_log_inv + 
RTchoice_RTtrim3mean_all_log_inv + RTchoice_RTsd_all_log_inv '

fit_CF_1 <- cfa(model_CF_1, data=mydata, std.lv=T,std.ov=T,missing='fiml', estimator='mlr') 

### Try with 1 vs 2 factors in CFA, and in EFA
# Actually, this doesn't work. the 2-factor model is saturated with DoF=-1, and negative variances. 
# This is because it has 4 observations, which each have a path leading to a latent variable, i.e. 4 paths,
# and a additional path connects (covaries) the latent variables. This gives 5 paths in total, but
# only 4 input observations. Hence DoF=-1. 


# Visualise the 1-factor CFA o speed to check it looks as it should:
par(mfrow=c(1,1))
semPaths(fit_CF_1, # change to each fit
         "par",
         sizeMan = 7, sizeInt = 7, sizeLat=7, fade=FALSE,
         edge.label.cex=.7,# Make edge labels readable
         residuals = FALSE,  # Remove the residual variance, within each variable
         intercepts = FALSE,# Remove intercepts
         rotation = 3,
         nCharNodes=9)

# We don't report fit metrics for this theory-driven CFA. 

### Get the path  loadings to create a figure in powerpoint: 
# # The summary shows the estimates which are the path loadings.
summary(fit_CF_1, standardized=TRUE)

# Extract values for participant loadings onto each latent var, created in CFA:
CF <- predict(fit_CF_1) 

# Put participant loadings onto each latent variable into the main dataframe:
mydata <- data.frame (mydata[,],CF)
tail(colnames(mydata))


# Tidy workspace
rm(model_CF_1, fit_CF_1, fit_CF_2,  CF, i , model_CF_1, model_CF_2, my_fn, 
   fit_efa_1, fit_efa_2,  mydata_cor)

#------------------------------------------------------------------------------#
### Sanity checks 

### Check imputation worked 
count(mydata$Speed != "NaN" )
temp<-subset(mydata, Speed != "NaN", c="Speed" )
rm(temp)
# n=664 have complete data after CFA

# #------------------------------------------------------------------------------#
### Save data

setwd("/.../")

# Save as RDS for use in future scripts. 
saveRDS(mydata, file="2.2_mydata.rds")

# Save as CSV for viewing externally.
write.csv(mydata, "2.2_mydata.csv")

#---------------------------------------------------------------------------------------------------------------------------#
### correlation plot of linear LVFs, PSMD, latent cognitive factors, Age

# Figure 4.3 in PhD thesis, 
# Figure 3 in paper.

mydata_cor <-  ( data.frame(mydata$Age,
                            mydata$LVF1, mydata$LVF2, mydata$LVF3,
                            mydata$PSMD_k,
                            mydata$Stripe_index,
                            mydata$Speed , mydata$LCF2
) )


colnames(mydata_cor) <- c("Age", 
                          "SSBP", "PP" , "HRV",
                          "PSMD",
                          "Motion",
                          "Speed",  "Fluid" #, "Discrepancy", 
)

# Missing values are not allowed with poly(), so remove NaNs
mydata_cor <- na.omit(mydata_cor)

par(mfrow=c(1,1))

my_fn <- function (data, mapping ) {
  p <- ggplot(data=mydata_cor,mapping=mapping)+
    # geom_point()  +  
    geom_jitter(alpha=0.1) +
    geom_smooth (method=lm, color = "red" ,se=FALSE)+
    geom_smooth (method=lm, formula = y ~ x + I(x^2), size=1,
                 color = "blue",se=FALSE) 
  # plot the linear line second, so we only see red quadratic if they're different
  p
}

g = ggpairs(mydata_cor,
            lower=list(continuous=my_fn) )
g

myplot_cor <- g
setwd("/.../")
ggsave("myplot_cor.png")


#---------------------------------------------------------------------------------------------------------------------------#
### Repeat with quadratic LVFs 

# Figure B.2 in the thesis.
# Figure 2 of Appendix in paper.

mydata_cor$Age_q <- poly(mydata_cor$Age,2)[,2]
mydata_cor$SSBP_q <- poly(mydata_cor$SSBP,2)[,2]
mydata_cor$PP_q <- poly(mydata_cor$PP,2)[,2]
mydata_cor$HRV_q <- poly(mydata_cor$HRV,2)[,2]

# Check newly named quadratic cols
colnames(mydata_cor)

# only select some variables
mydata_cor <- subset(mydata_cor, select=-c(SSBP, PP, HRV))
colnames(mydata_cor)

# Re-order by col names
col_order <- c("Age", "Age_q",
               "SSBP_q", "PP_q","HRV_q",
               "PSMD", "Motion", # meaning head motion
               "Speed", "Fluid")
mydata_cor<-mydata_cor[,col_order]
colnames(mydata_cor)

# Missing values are not allowed with poly(), so remove NaNs
mydata_cor <- na.omit(mydata_cor)

par(mfrow=c(1,1))

my_fn <- function (data, mapping ) {
  p <- ggplot(data=mydata_cor,mapping=mapping)+
    # geom_point()  +  
    geom_jitter(alpha=0.1) +
    geom_smooth (method=lm, color = "red" ,se=FALSE) +
    geom_smooth (method=lm, formula = y ~ x + I(x^2), size=1,
                 color = "blue",se=FALSE)
  p
}

g = ggpairs(mydata_cor,
            lower=list(continuous=my_fn) )
g

myplot_cor <- g

setwd("/.../")
ggsave("myplot_cor_LVFquadratic.png")

rm(mydata)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

### ENDS





