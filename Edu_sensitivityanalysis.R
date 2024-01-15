setwd("/Users/joycelin/Desktop/Gulf/Aim3/GuLF_3")
library(tidyverse)
library(dplyr)
library(readxl)
library(Hmisc) 
library(gtsummary)
library(summarytools)
library(ggplot2)
library(dotwhisker)

## Loading data and data management  -----------------------------------------------
data<- read.csv("NEW_CE_neuro_quantiles.csv")
data$race <- relevel(factor(data$race), ref = "White")
data$QTHC <- factor(data$QTHC)
data$EN_FORMERSMOKER <- factor(data$EN_FORMERSMOKER)
data$THC_CUMULATIVE1<- as.numeric(data$THC_CUMULATIVE1)
data$headinjuryever <- factor(data$headinjuryever)
data<- data %>% mutate(drinklot = case_when(drinksweek <= 14 ~ '0',
                                            drinksweek >14 ~ '1')) # more than 14 drinks/week
data$drinklot <- factor(data$drinklot)
data$Caffiene2h <- factor(data$Caffiene2h)

data<- data %>% mutate(marital = case_when(EN_MARITAL == 1 ~ 'Married or with partner',
                                           EN_MARITAL == 6 ~ 'Married or with partner',
                                           EN_MARITAL == 2 ~ 'Separated or single',
                                           EN_MARITAL == 4 ~ 'Separated or single',
                                           EN_MARITAL == 3 ~ 'Separated or single',
                                           EN_MARITAL == 5 ~ 'Separated or single')) 
data$marital <- factor(data$marital)

data<- data %>% mutate(BMI = case_when(CE_BMI <18.5 ~ 'underweight',
                                       CE_BMI >=18.5 & CE_BMI < 25 ~ 'healthy',
                                       CE_BMI >=25 & CE_BMI <=30 ~ 'overweight',
                                       CE_BMI >=30 ~ 'obese')) # end function
data$BMI <- factor(data$BMI)

#log10 transform metal concentrations
data_log <- data %>% mutate(across(c(Mg, Al, Ca, Cr, Mn, Fe, Ni, Cu, Zn, As, Se, Hg, Pb), log10))

#dichotomize outcome for fractional responses (CPT CR Fraction, CPT hit fraction)
data_log <- data_log %>% mutate(CPTcr1 = case_when(CE_BARS_CPT_CR_FRACTION < 0.80 ~ '0',
                                                   CE_BARS_CPT_CR_FRACTION >= 0.8 ~ '1')) 
data_log <- data_log %>% mutate(CPThit1 = case_when(CE_BARS_CPT_HIT_FRACTION < 0.80 ~ '0',
                                                    CE_BARS_CPT_HIT_FRACTION >= 0.8 ~ '1')) 
data_log$CPTcr1 <- as.numeric(data_log$CPTcr1)
data_log$CPThit1 <- as.numeric(data_log$CPThit1)

data_log <- data_log %>% 
  rename("Race" = "race")
blacksubgroup<-filter(data_log, Race =="Black")
whitesubgroup<-filter(data_log, Race =="White")
othersubgroup<-filter(data_log, Race =="Other")

## education grouping
lowed<-filter(data_log, CE_C1 <=15)
highed<-filter(data_log, CE_C1 >=16)

summaryhighed <- highed %>%
  group_by(Race) %>%
  summarise(Count = n())
summaryhighed

## Stratified by race -------------------------------------------------------
## adjusted for other metals continuously
mod1<- lm(CE_BARS_CPT_D_PRIME ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME  ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME  ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

summary(mod1)
confint(mod1)
summary(mod2)
confint(mod2)
summary(mod3)
confint(mod3)

mod1<- glm(CPTcr1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPTcr1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPTcr1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

summary(mod1)
confint(mod1)
summary(mod2)
summary(mod3)
confint(mod3)

## Adjusted scaled by IQR
#Cr
Criqr<-IQR(data_log$Cr) 
WCriqr<-IQR(data_log$Cr) 
BCriqr<-IQR(data_log$Cr) 
data_log<-data_log %>% mutate(CrIQR = Cr/Criqr)
whitesubgroup<-whitesubgroup %>% mutate(WCrIQR = Cr/WCriqr)
blacksubgroup<-blacksubgroup %>% mutate(BCrIQR = Cr/BCriqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ CrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ CrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WCrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WCrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BCrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BCrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

#Cu
Cuiqr<-IQR(data_log$Cu) 
WCuiqr<-IQR(data_log$Cu) 
BCuiqr<-IQR(data_log$Cu) 
data_log<-data_log %>% mutate(CuIQR = Cu/Cuiqr)
whitesubgroup<-whitesubgroup %>% mutate(WCuIQR = Cu/WCuiqr)
blacksubgroup<-blacksubgroup %>% mutate(BCuIQR = Cu/BCuiqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ CuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ CuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WCuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WCuIQR + As + Mn + Cr + Zn + Pb + Hg + Se+ CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BCuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BCuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# As
Asiqr<-IQR(data_log$As) 
WAsiqr<-IQR(data_log$As) 
BAsiqr<-IQR(data_log$As) 
data_log<-data_log %>% mutate(AsIQR = As/Asiqr)
whitesubgroup<-whitesubgroup %>% mutate(WAsIQR = As/WAsiqr)
blacksubgroup<-blacksubgroup %>% mutate(BAsIQR = As/BAsiqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ AsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ AsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# Hg
Hgiqr<-IQR(data_log$Hg) 
WHgiqr<-IQR(data_log$Hg) 
BHgiqr<-IQR(data_log$Hg) 
data_log<-data_log %>% mutate(HgIQR = Hg/Hgiqr)
whitesubgroup<-whitesubgroup %>% mutate(WHgIQR = Hg/WHgiqr)
blacksubgroup<-blacksubgroup %>% mutate(BHgIQR = Hg/BHgiqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ HgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ HgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WHgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WHgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BHgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BHgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# Mn
Mniqr<-IQR(data_log$Mn) 
WMniqr<-IQR(data_log$Mn) 
BMniqr<-IQR(data_log$Mn) 
data_log<-data_log %>% mutate(MnIQR = Mn/Mniqr)
whitesubgroup<-whitesubgroup %>% mutate(WMnIQR = Mn/WMniqr)
blacksubgroup<-blacksubgroup %>% mutate(BMnIQR = Mn/BMniqr)

mod1<- lm(CE_BARS_DST_REVERSE_CNT ~ MnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log)
mod1<- glm(CPTcr1 ~ MnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm(CE_BARS_DST_REVERSE_CNT ~ WMnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WMnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm(CE_BARS_DST_REVERSE_CNT ~ BMnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BMnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# Pb
Pbiqr<-IQR(data_log$Pb) 
WPbiqr<-IQR(data_log$Pb) 
BPbiqr<-IQR(data_log$Pb) 
data_log<-data_log %>% mutate(PbIQR = Pb/Pbiqr)
whitesubgroup<-whitesubgroup %>% mutate(WPbIQR = Pb/WPbiqr)
blacksubgroup<-blacksubgroup %>% mutate(BPbIQR = Pb/BPbiqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ PbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log)
mod1<- glm(CPTcr1 ~ PbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WPbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WPbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BPbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BPbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# Se
Seiqr<-IQR(data_log$Se) 
WSeiqr<-IQR(data_log$Se) 
BSeiqr<-IQR(data_log$Se) 
data_log<-data_log %>% mutate(SeIQR = Se/Seiqr)
whitesubgroup<-whitesubgroup %>% mutate(WSeIQR = Se/WSeiqr)
blacksubgroup<-blacksubgroup %>% mutate(BSeIQR = Se/BSeiqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ SeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ SeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WSeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WSeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BSeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BSeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# Zn
Zniqr<-IQR(data_log$Zn) 
WZniqr<-IQR(data_log$Zn) 
BZniqr<-IQR(data_log$Zn) 
data_log<-data_log %>% mutate(ZnIQR = Zn/Zniqr)
whitesubgroup<-whitesubgroup %>% mutate(WZnIQR = Zn/WZniqr)
blacksubgroup<-blacksubgroup %>% mutate(BZnIQR = Zn/BZniqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ ZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ ZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb+ CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

## dot whisker plot IQR scaled ---------------------------------------------
mod1<- lm(CE_BARS_CPT_D_PRIME ~ AsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ WAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ BAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod4<- lm(CE_BARS_CPT_D_PRIME ~ CrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod5<- lm(CE_BARS_CPT_D_PRIME ~ WCrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod6<- lm(CE_BARS_CPT_D_PRIME ~ BCrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod7<- lm(CE_BARS_CPT_D_PRIME ~ CuIQR + Cr + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod8<- lm(CE_BARS_CPT_D_PRIME ~ WCuIQR + Cr + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod9<- lm(CE_BARS_CPT_D_PRIME ~ BCuIQR + Cr + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod10<- lm(CE_BARS_CPT_D_PRIME ~ HgIQR + Cu + As + Cr + Zn + Pb + Mn + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod11<- lm(CE_BARS_CPT_D_PRIME ~ WHgIQR + Cu + As + Cr + Zn + Pb + Mn + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod12<- lm(CE_BARS_CPT_D_PRIME ~ BHgIQR + Cu + As + Cr + Zn + Pb + Mn + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod13<- lm(CE_BARS_CPT_D_PRIME ~ MnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod14<- lm(CE_BARS_CPT_D_PRIME ~ WMnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod15<- lm(CE_BARS_CPT_D_PRIME ~ BMnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod16<- lm(CE_BARS_CPT_D_PRIME ~ PbIQR + Cu + As + Cr + Zn + Mn + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod17<- lm(CE_BARS_CPT_D_PRIME ~ WPbIQR + Cu + As + Cr + Zn + Mn + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod18<- lm(CE_BARS_CPT_D_PRIME ~ BPbIQR + Cu + As + Cr + Zn + Mn + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod19<- lm(CE_BARS_CPT_D_PRIME ~ SeIQR + Cu + As + Cr + Zn + Pb + Hg + Mn + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod20<- lm(CE_BARS_CPT_D_PRIME ~ WSeIQR + Cu + As + Cr + Zn + Pb + Hg + Mn + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod21<- lm(CE_BARS_CPT_D_PRIME ~ BSeIQR + Cu + As + Cr + Zn + Pb + Hg + Mn + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod22<- lm(CE_BARS_CPT_D_PRIME ~ ZnIQR + Cu + As + Cr + Mn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod23<- lm(CE_BARS_CPT_D_PRIME ~ WZnIQR + Cu + As + Cr + Mn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod24<- lm(CE_BARS_CPT_D_PRIME ~ BZnIQR + Cu + As + Cr + Mn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)

mod1 <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod2 <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod3 <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod4 <-broom::tidy(mod4) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod5 <-broom::tidy(mod5) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod6 <-broom::tidy(mod6) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod7 <-broom::tidy(mod7) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod8 <-broom::tidy(mod8) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod9 <-broom::tidy(mod9) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod10 <-broom::tidy(mod10) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod11 <-broom::tidy(mod11) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod12 <-broom::tidy(mod12) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod13 <-broom::tidy(mod13) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod14 <-broom::tidy(mod14) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod15 <-broom::tidy(mod15) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod16 <-broom::tidy(mod16) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod17 <-broom::tidy(mod17) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod18 <-broom::tidy(mod18) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod19 <-broom::tidy(mod19) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod20 <-broom::tidy(mod20) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod21 <-broom::tidy(mod21) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod22 <-broom::tidy(mod22) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod23 <-broom::tidy(mod23) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod24 <-broom::tidy(mod24) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")

all_models <- combine(mod1, mod2, mod3, mod4, mod5, mod6,mod7, mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75)) +
  scale_x_continuous(limits=c(-0.45,0.32),breaks=c( -0.4, -0.2, 0, 0.2))
DSTr_Mn


## Stratified by education level --------------------------------------------------------
## adjusted for other metals continuously
mod1<- lm(CE_BARS_CPT_D_PRIME ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME  ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = highed)
mod3<- lm(CE_BARS_CPT_D_PRIME  ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = lowed)

summary(mod1)
confint(mod1)
summary(mod2)
confint(mod2)
summary(mod3)
confint(mod3)

mod1<- glm(CPTcr1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPTcr1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = highed, family = binomial)
mod3<- glm(CPTcr1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = lowed, family = binomial)

summary(mod1)
confint(mod1)
summary(mod2)
summary(mod3)
confint(mod3)

mod1<- glm(CPThit1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPThit1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = highed, family = binomial)
mod3<- glm(CPThit1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = lowed, family = binomial)

summary(mod1)
confint(mod1)
summary(mod2)
summary(mod3)
confint(mod3)

## Adjusted scaled by IQR, should remove CE_C1 as a variable from models? yes?
#Cr
Criqr<-IQR(data_log$Cr) 
HCriqr<-IQR(data_log$Cr) 
LCriqr<-IQR(data_log$Cr) 
data_log<-data_log %>% mutate(CrIQR = Cr/Criqr)
highed<-highed %>% mutate(HCrIQR = Cr/HCriqr)
lowed<-lowed %>% mutate(LCrIQR = Cr/LCriqr)

mod1<- lm( CE_BARS_CPT_D_PRIME ~ CrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_CPT_D_PRIME ~ HCrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_CPT_D_PRIME ~ LCrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
summary(mod3)
confint(mod3)

#Cu
Cuiqr<-IQR(data_log$Cu) 
HCuiqr<-IQR(data_log$Cu) 
LCuiqr<-IQR(data_log$Cu) 
data_log<-data_log %>% mutate(CuIQR = Cu/Cuiqr)
highed<-highed %>% mutate(HCuIQR = Cu/HCuiqr)
lowed<-lowed %>% mutate(LCuIQR = Cu/LCuiqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ CuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WCuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BCuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE+ EN_FORMERSMOKER +   marital, data = lowed)
summary(mod3)
confint(mod3)

# As
Asiqr<-IQR(data_log$As) 
HAsiqr<-IQR(data_log$As) 
LAsiqr<-IQR(data_log$As) 
data_log<-data_log %>% mutate(AsIQR = As/Asiqr)
highed<-highed %>% mutate(HAsIQR = As/HAsiqr)
lowed<-lowed %>% mutate(LAsIQR = As/LAsiqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ AsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
summary(mod3)
confint(mod3)

# Hg
Hgiqr<-IQR(data_log$Hg) 
HHgiqr<-IQR(data_log$Hg) 
LHgiqr<-IQR(data_log$Hg) 
data_log<-data_log %>% mutate(HgIQR = Hg/Hgiqr)
highed<-highed %>% mutate(HHgIQR = Hg/HHgiqr)
lowed<-lowed %>% mutate(LHgIQR = Hg/LHgiqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ HgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WHgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BHgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
summary(mod3)
confint(mod3)

# Mn
Mniqr<-IQR(data_log$Mn) 
HMniqr<-IQR(data_log$Mn) 
LMniqr<-IQR(data_log$Mn) 
data_log<-data_log %>% mutate(MnIQR = Mn/Mniqr)
highed<-highed %>% mutate(HMnIQR = Mn/HMniqr)
lowed<-lowed %>% mutate(LMnIQR = Mn/LMniqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ MnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ HMnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ LMnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
summary(mod3)
confint(mod3)

# Pb
Pbiqr<-IQR(data_log$Pb) 
HPbiqr<-IQR(data_log$Pb) 
LPbiqr<-IQR(data_log$Pb) 
data_log<-data_log %>% mutate(PbIQR = Pb/Pbiqr)
highed<-highed %>% mutate(HPbIQR = Pb/HPbiqr)
lowed<-lowed %>% mutate(LPbIQR = Pb/LPbiqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ PbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WPbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE  + EN_FORMERSMOKER +   marital, data = highed)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BPbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
summary(mod3)
confint(mod3)

# Se
Seiqr<-IQR(data_log$Se) 
HSeiqr<-IQR(data_log$Se) 
LSeiqr<-IQR(data_log$Se) 
data_log<-data_log %>% mutate(SeIQR = Se/Seiqr)
highed<-highed %>% mutate(HSeIQR = Se/HSeiqr)
lowed<-lowed %>% mutate(LSeIQR = Se/LSeiqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ SeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ SeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WSeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BSeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
summary(mod3)
confint(mod3)

# Zn
Zniqr<-IQR(data_log$Zn) 
HZniqr<-IQR(data_log$Zn) 
LZniqr<-IQR(data_log$Zn) 
data_log<-data_log %>% mutate(ZnIQR = Zn/Zniqr)
highed<-highed %>% mutate(HZnIQR = Zn/HZniqr)
lowed<-lowed %>% mutate(LZnIQR = Zn/LZniqr)

mod1<- lm( CE_BARS_DST_REVERSE_CNT ~ ZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ ZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_DST_REVERSE_CNT ~ WZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_DST_REVERSE_CNT ~ BZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
summary(mod3)
confint(mod3)

## dot whisker plot IQR scaled ---------------------------------------------
mod1<- lm(CE_BARS_CPT_D_PRIME ~ AsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ HAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ LAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
mod4<- lm(CE_BARS_CPT_D_PRIME ~ CrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = data_log)
mod5<- lm(CE_BARS_CPT_D_PRIME ~ HCrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
mod6<- lm(CE_BARS_CPT_D_PRIME ~ LCrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
mod7<- lm(CE_BARS_CPT_D_PRIME ~ CuIQR + Cr + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = data_log)
mod8<- lm(CE_BARS_CPT_D_PRIME ~ HCuIQR + Cr + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
mod9<- lm(CE_BARS_CPT_D_PRIME ~ LCuIQR + Cr + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
mod10<- lm(CE_BARS_CPT_D_PRIME ~ HgIQR + Cu + As + Cr + Zn + Pb + Mn + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = data_log)
mod11<- lm(CE_BARS_CPT_D_PRIME ~ HHgIQR + Cu + As + Cr + Zn + Pb + Mn + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
mod12<- lm(CE_BARS_CPT_D_PRIME ~ LHgIQR + Cu + As + Cr + Zn + Pb + Mn + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
mod13<- lm(CE_BARS_CPT_D_PRIME ~ MnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = data_log)
mod14<- lm(CE_BARS_CPT_D_PRIME ~ HMnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
mod15<- lm(CE_BARS_CPT_D_PRIME ~ LMnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
mod16<- lm(CE_BARS_CPT_D_PRIME ~ PbIQR + Cu + As + Cr + Zn + Mn + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = data_log)
mod17<- lm(CE_BARS_CPT_D_PRIME ~ HPbIQR + Cu + As + Cr + Zn + Mn + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
mod18<- lm(CE_BARS_CPT_D_PRIME ~ LPbIQR + Cu + As + Cr + Zn + Mn + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
mod19<- lm(CE_BARS_CPT_D_PRIME ~ SeIQR + Cu + As + Cr + Zn + Pb + Hg + Mn + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = data_log)
mod20<- lm(CE_BARS_CPT_D_PRIME ~ HSeIQR + Cu + As + Cr + Zn + Pb + Hg + Mn + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
mod21<- lm(CE_BARS_CPT_D_PRIME ~ LSeIQR + Cu + As + Cr + Zn + Pb + Hg + Mn + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)
mod22<- lm(CE_BARS_CPT_D_PRIME ~ ZnIQR + Cu + As + Cr + Mn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = data_log)
mod23<- lm(CE_BARS_CPT_D_PRIME ~ HZnIQR + Cu + As + Cr + Mn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = highed)
mod24<- lm(CE_BARS_CPT_D_PRIME ~ LZnIQR + Cu + As + Cr + Mn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +   marital, data = lowed)

mod1 <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod2 <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod3 <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod4 <-broom::tidy(mod4) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod5 <-broom::tidy(mod5) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod6 <-broom::tidy(mod6) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod7 <-broom::tidy(mod7) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod8 <-broom::tidy(mod8) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod9 <-broom::tidy(mod9) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod10 <-broom::tidy(mod10) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod11 <-broom::tidy(mod11) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod12 <-broom::tidy(mod12) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod13 <-broom::tidy(mod13) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod14 <-broom::tidy(mod14) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod15 <-broom::tidy(mod15) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod16 <-broom::tidy(mod16) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod17 <-broom::tidy(mod17) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod18 <-broom::tidy(mod18) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod19 <-broom::tidy(mod19) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod20 <-broom::tidy(mod20) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod21 <-broom::tidy(mod21) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod22 <-broom::tidy(mod22) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod23 <-broom::tidy(mod23) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod24 <-broom::tidy(mod24) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")

all_models <- combine(mod1, mod2, mod3, mod4, mod5, mod6,mod7, mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75)) +
  scale_x_continuous(limits=c(-0.75,0.6),breaks=c( -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6))
DSTr_Mn

# EDU and RACE combined plots for Cr, Mn, Cu, Zn only ------------------------
mod1<- lm(CE_BARS_CPT_D_PRIME ~ CrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER + CE_C1 + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ WCrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ BCrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = blacksubgroup)
mod4<- lm(CE_BARS_CPT_D_PRIME ~ HCrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +  marital, data = highed)
mod5<- lm(CE_BARS_CPT_D_PRIME ~ LCrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +  marital, data = lowed)

mod6<- lm(CE_BARS_CPT_D_PRIME ~ MnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +  CE_C1 + marital, data = data_log)
mod7<- lm(CE_BARS_CPT_D_PRIME ~ WMnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +  marital, data = whitesubgroup)
mod8<- lm(CE_BARS_CPT_D_PRIME ~ BMnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +  marital, data = blacksubgroup)
mod9<- lm(CE_BARS_CPT_D_PRIME ~ HMnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER + marital, data = highed)
mod10<- lm(CE_BARS_CPT_D_PRIME ~ LMnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER + marital, data = lowed)

mod11<- lm(CE_BARS_CPT_D_PRIME ~ CuIQR + Mn + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +  CE_C1 + marital, data = data_log)
mod12<- lm(CE_BARS_CPT_D_PRIME ~ WCuIQR + Mn + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +  marital, data = whitesubgroup)
mod13<- lm(CE_BARS_CPT_D_PRIME ~ BCuIQR + Mn + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +  marital, data = blacksubgroup)
mod14<- lm(CE_BARS_CPT_D_PRIME ~ HCuIQR + Mn + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER + marital, data = highed)
mod15<- lm(CE_BARS_CPT_D_PRIME ~ LCuIQR + Mn + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER + marital, data = lowed)

mod16<- lm(CE_BARS_CPT_D_PRIME ~ ZnIQR + Mn + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER +  CE_C1 + marital, data = data_log)
mod17<- lm(CE_BARS_CPT_D_PRIME ~ WZnIQR + Mn + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +  marital, data = whitesubgroup)
mod18<- lm(CE_BARS_CPT_D_PRIME ~ BZnIQR + Mn + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +  marital, data = blacksubgroup)
mod19<- lm(CE_BARS_CPT_D_PRIME ~ HZnIQR + Mn + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER + marital, data = highed)
mod20<- lm(CE_BARS_CPT_D_PRIME ~ LZnIQR + Mn + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + EN_FORMERSMOKER + marital, data = lowed)

mod1 <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod2 <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod3 <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod4 <-broom::tidy(mod4) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod5 <-broom::tidy(mod5) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod6 <-broom::tidy(mod6) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod7 <-broom::tidy(mod7) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod8 <-broom::tidy(mod8) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod9 <-broom::tidy(mod9) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod10 <-broom::tidy(mod10) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod11 <-broom::tidy(mod11) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod12 <-broom::tidy(mod12) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod13 <-broom::tidy(mod13) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod14 <-broom::tidy(mod14) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod15 <-broom::tidy(mod15) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod16 <-broom::tidy(mod16) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod17 <-broom::tidy(mod17) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod18 <-broom::tidy(mod18) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod19 <-broom::tidy(mod19) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")
mod20 <-broom::tidy(mod20) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Se")%>% filter(term != "Hg")%>% filter(term != "Zn")

all_models <- combine(mod1, mod2, mod3, mod4, mod5, mod6,mod7, mod8,mod9,mod10, mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75)) +
  scale_x_continuous(limits=c(-0.45,0.3),breaks=c(-0.4, -0.2, 0,0.2, 0.4))
DSTr_Mn
