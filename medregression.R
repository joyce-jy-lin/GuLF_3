# Median regression

setwd("/Users/joycelin/Desktop/Gulf/Aim3/GuLF_3")
library(tidyverse)
library(dplyr)
library(readxl)
library(Hmisc) 
library(gtsummary)
library(summarytools)
library(ggplot2)

## descriptive plots metals and neuro by quartile of metal -----------------------------------------------
data<- read.csv("NEW_CE_neuro_quantiles.csv")

data$AlQ <- factor(data$AlQ)
data$AsQ <- factor(data$AsQ)
data$CaQ <- factor(data$CaQ)
data$CrQ <- factor(data$CrQ)
data$CuQ <- factor(data$CuQ)
data$FeQ <- factor(data$FeQ)
data$PbQ <- factor(data$PbQ)
data$MgQ <- factor(data$MgQ)
data$MnQ <- factor(data$MnQ)
data$HgQ <- factor(data$HgQ)
data$NiQ <- factor(data$NiQ)
data$SeQ <- factor(data$SeQ)
data$ZnQ <- factor(data$ZnQ)
data$CdBinary <- factor(data$CdBinary)
data$CoBinary <- factor(data$CoBinary)
data$MoBinary <- factor(data$MoBinary)
data$SbTertile <- factor(data$SbTertile)
data$VTertile <- factor(data$VTertile)
data$race <- relevel(factor(data$race), ref = "White")
data$QTHC <- factor(data$QTHC)
data$EN_FORMERSMOKER <- factor(data$EN_FORMERSMOKER)
data$THC_CUMULATIVE1<- as.numeric(data$THC_CUMULATIVE1)
data$headinjuryever <- factor(data$headinjuryever)
data$drinklot <- factor(data$drinklot)
data$Caffiene2h <- factor(data$Caffiene2h)

data<- data %>% mutate(marital = case_when(EN_MARITAL == 1 ~ 'Married or with partner',
                                           EN_MARITAL == 6 ~ 'Married or with partner',
                                           EN_MARITAL == 2 ~ 'Separated or single',
                                           EN_MARITAL == 4 ~ 'Separated or single',
                                           EN_MARITAL == 3 ~ 'Separated or single',
                                           EN_MARITAL == 5 ~ 'Separated or single')) 

data$marital <- factor(data$marital)

data<- data %>% mutate(BMI = case_when(CE_BMI <= 24.9 ~ 'healthy',
                                       CE_BMI >=25 & CE_BMI <=30 ~ 'overweight',
                                       CE_BMI >=30 ~ 'obese')) # end function

data$BMI <- factor(data$BMI)
data <- data %>% 
  rename("Race" = "race")

blacksubgroup<-filter(data, Race =="Black")
whitesubgroup<-filter(data, Race =="White")

data_log <- data %>% mutate(across(c(Mg, Al, Ca, Cr, Mn, Fe, Ni, Cu, Zn, As, Se, Hg, Pb), log10))

## dichotomize outcome for fractional responses (CPT CR Fraction, CPT hit fraction)
data_log <- data_log %>% mutate(CPTcr1 = case_when(CE_BARS_CPT_CR_FRACTION < 0.80 ~ '1',
                                                   CE_BARS_CPT_CR_FRACTION >= 0.8 ~ '0')) 
data_log <- data_log %>% mutate(CPThit1 = case_when(CE_BARS_CPT_HIT_FRACTION < 0.80 ~ '1',
                                                    CE_BARS_CPT_HIT_FRACTION >= 0.8 ~ '0')) 
data_log$CPTcr1 <- as.numeric(data_log$CPTcr1)
data_log$CPThit1 <- as.numeric(data_log$CPThit1)

# Regular linear regression ----------------------------
mod1 <- lm(CE_BARS_DST_REVERSE_CNT ~ MnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + headinjuryever + drinklot + marital, data = data_log)
summary(mod1)
confint(mod1)

# Quantile regression with bootstrap SE ----Mn DST reverse-------------
library(quantreg)
mod1 <- rq(CE_BARS_DST_REVERSE_CNT ~ MnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + headinjuryever + drinklot + marital, data = data_log)
mod1<-summary.rq(mod1, se="boot")
coef=mod1$coefficients[,1]
err=mod1$coefficients[,2]
ci<- list()
for (i in 1:length(coef)){
  ci[[i]] <- coef[i] + c(-1,1)*err[i]*qnorm(0.975)}
ci_dat <- data.frame(t(sapply(ci,c)))
ci_dat %>% 
  rename(
    lower = X1,
    upper = X2)

wmod1 <- rq(CE_BARS_DST_REVERSE_CNT ~ MnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + headinjuryever + drinklot + marital, data = whitesubgroup)
wmod1<-summary.rq(wmod1, se="boot")
wcoef=wmod1$coefficients[,1]
werr=wmod1$coefficients[,2]
wci<- list()
for (i in 1:length(wcoef)){
  wci[[i]] <- wcoef[i] + c(-1,1)*err[i]*qnorm(0.975)}
wci
wci_dat <- data.frame(t(sapply(wci,c)))
wci_dat %>% 
  rename(
    lower = X1,
    upper = X2)

bmod1 <- rq(CE_BARS_DST_REVERSE_CNT ~ MnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + headinjuryever + drinklot + marital, data = blacksubgroup)
bmod1<-summary.rq(bmod1, se="boot")
bcoef=bmod1$coefficients[,1]
berr=bmod1$coefficients[,2]
bci<- list()
for (i in 1:length(bcoef)){
  bci[[i]] <- bcoef[i] + c(-1,1)*err[i]*qnorm(0.975)}
bci
bci_dat <- data.frame(t(sapply(bci,c)))
bci_dat %>% 
  rename(
    lower = X1,
    upper = X2)


Mn <- data.frame(x =c(coef["MnQ2"], coef["MnQ3"], coef["MnQ4"], wcoef["MnQ2"], wcoef["MnQ3"], wcoef["MnQ4"], bcoef["MnQ2"], bcoef["MnQ3"], bcoef["MnQ4"]),
                 y = c("Mn Q2", "Mn Q3", "Mn Q4","Mn Q2", "Mn Q3", "Mn Q4", "Mn Q2", "Mn Q3", "Mn Q4"),
                 lower = c(ci_dat[2,1], ci_dat[3,1], ci_dat[4,1], wci_dat[2,1], wci_dat[3,1], wci_dat[4,1], bci_dat[2,1], bci_dat[3,1], bci_dat[4,1]),
                 upper = c(ci_dat[2,2], ci_dat[3,2], ci_dat[4,2], wci_dat[2,2], wci_dat[3,2], wci_dat[4,2], bci_dat[2,2], bci_dat[3,2], bci_dat[4,2]),
                 Race = c("All", "All", "All", "White", "White", "White","Black", "Black","Black"))

Mn$y <- factor(Mn$y, c("Mn Q4", "Mn Q3", "Mn Q2"))
Mn$Race <- factor(Mn$Race, c("Black", "White", "All"))

library("ggplot2")
p<-ggplot(Mn, aes(y, x, color=Race, group=Race)) +     
  geom_point(position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0.2, position=position_dodge(width=0.7)) + 
  coord_flip() + geom_vline(xintercept = 0)
p

# Quantile regression with bootstrap SE --Cr D Prime---------------
library(quantreg)
mod1 <- rq(CE_BARS_CPT_D_PRIME ~ CrQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + headinjuryever + drinklot + marital, data = data_log)
mod1<-summary.rq(mod1, se="boot")
coef=mod1$coefficients[,1]
err=mod1$coefficients[,2]
ci<- list()
for (i in 1:length(coef)){
  ci[[i]] <- coef[i] + c(-1,1)*err[i]*qnorm(0.975)}
ci_dat <- data.frame(t(sapply(ci,c)))
ci_dat %>% 
  rename(
    lower = X1,
    upper = X2)

wmod1 <- rq(CE_BARS_CPT_D_PRIME ~ CrQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + headinjuryever + drinklot + marital, data = whitesubgroup)
wmod1<-summary.rq(wmod1, se="boot")
wcoef=wmod1$coefficients[,1]
werr=wmod1$coefficients[,2]
wci<- list()
for (i in 1:length(wcoef)){
  wci[[i]] <- wcoef[i] + c(-1,1)*err[i]*qnorm(0.975)}
wci
wci_dat <- data.frame(t(sapply(wci,c)))
wci_dat %>% 
  rename(
    lower = X1,
    upper = X2)

bmod1 <- rq(CE_BARS_CPT_D_PRIME ~ CrQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + headinjuryever + drinklot + marital, data = blacksubgroup)
bmod1<-summary.rq(bmod1, se="boot")
bcoef=bmod1$coefficients[,1]
berr=bmod1$coefficients[,2]
bci<- list()
for (i in 1:length(bcoef)){
  bci[[i]] <- bcoef[i] + c(-1,1)*err[i]*qnorm(0.975)}
bci
bci_dat <- data.frame(t(sapply(bci,c)))
bci_dat %>% 
  rename(
    lower = X1,
    upper = X2)


Cr <- data.frame(x =c(coef["CrQ2"], coef["CrQ3"], coef["CrQ4"], wcoef["CrQ2"], wcoef["CrQ3"], wcoef["CrQ4"], bcoef["CrQ2"], bcoef["CrQ3"], bcoef["CrQ4"]),
                 y = c("Cr Q2", "Cr Q3", "Cr Q4","Cr Q2", "Cr Q3", "Cr Q4", "Cr Q2", "Cr Q3", "Cr Q4"),
                 lower = c(ci_dat[2,1], ci_dat[3,1], ci_dat[4,1], wci_dat[2,1], wci_dat[3,1], wci_dat[4,1], bci_dat[2,1], bci_dat[3,1], bci_dat[4,1]),
                 upper = c(ci_dat[2,2], ci_dat[3,2], ci_dat[4,2], wci_dat[2,2], wci_dat[3,2], wci_dat[4,2], bci_dat[2,2], bci_dat[3,2], bci_dat[4,2]),
                 Race = c("All", "All", "All", "White", "White", "White","Black", "Black","Black"))

Cr$y <- factor(Cr$y, c("Cr Q4", "Cr Q3", "Cr Q2"))
Cr$Race <- factor(Cr$Race, c("Black", "White", "All"))

library("ggplot2")
p<-ggplot(Cr, aes(y, x, color=Race, group=Race)) +     
  geom_point(position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0.2, position=position_dodge(width=0.7)) + 
  coord_flip() + geom_vline(xintercept = 0)
p

