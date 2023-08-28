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

data<- data %>% mutate(BMI = case_when(CE_BMI <= 24.9 ~ 'healthy',
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


## adjusted for other metals continuously -------------------------------
mod1<- lm(CE_BARS_DST_REVERSE_CNT ~ Mn +As + Cr + Cu + Zn + Pb + Hg + Se + headinjury + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME  ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME  ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = blacksubgroup)

summary(mod1)
confint(mod1)
summary(mod2)
confint(mod2)
summary(mod3)
confint(mod3)

mod1<- glm(CPTcr1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log, family = binomial)
mod2<- glm(CPTcr1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPTcr1 ~ Cr + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = blacksubgroup, family = binomial)

summary(mod1)
confint(mod1)
summary(mod2)
summary(mod3)
confint(mod3)

## Adjusted scaled by IQR
#Cr
Criqr<-IQR(data_log$Cr) 
WCriqr<-IQR(whitesubgroup$Cr) 
BCriqr<-IQR(blacksubgroup$Cr) 
data_log<-data_log %>% mutate(CrIQR = Cr/Criqr)
whitesubgroup<-whitesubgroup %>% mutate(WCrIQR = Cr/WCriqr)
blacksubgroup<-blacksubgroup %>% mutate(BCrIQR = Cr/BCriqr)

mod1<- lm( CE_BARS_CPT_D_PRIME ~ CrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ CrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_CPT_D_PRIME ~ WCrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WCrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_CPT_D_PRIME ~ BCrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BCrIQR + As + Mn + Cu + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

#Cu
Cuiqr<-IQR(data_log$Cu) 
WCuiqr<-IQR(whitesubgroup$Cu) 
BCuiqr<-IQR(blacksubgroup$Cu) 
data_log<-data_log %>% mutate(CuIQR = Cu/Cuiqr)
whitesubgroup<-whitesubgroup %>% mutate(WCuIQR = Cu/WCuiqr)
blacksubgroup<-blacksubgroup %>% mutate(BCuIQR = Cu/BCuiqr)

mod1<- lm( CE_BARS_CPT_D_PRIME ~ CuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ CuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_CPT_D_PRIME ~ WCuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WCuIQR + As + Mn + Cr + Zn + Pb + Hg + Se+ CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_CPT_D_PRIME ~ BCuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BCuIQR + As + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# As
Asiqr<-IQR(data_log$As) 
WAsiqr<-IQR(whitesubgroup$As) 
BAsiqr<-IQR(blacksubgroup$As) 
data_log<-data_log %>% mutate(AsIQR = As/Asiqr)
whitesubgroup<-whitesubgroup %>% mutate(WAsIQR = As/WAsiqr)
blacksubgroup<-blacksubgroup %>% mutate(BAsIQR = As/BAsiqr)

mod1<- lm( CE_BARS_CPT_D_PRIME ~ AsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ AsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_CPT_D_PRIME ~ WAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_CPT_D_PRIME ~ BAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# Hg
Hgiqr<-IQR(data_log$Hg) 
WHgiqr<-IQR(whitesubgroup$Hg) 
BHgiqr<-IQR(blacksubgroup$Hg) 
data_log<-data_log %>% mutate(HgIQR = Hg/Hgiqr)
whitesubgroup<-whitesubgroup %>% mutate(WHgIQR = Hg/WHgiqr)
blacksubgroup<-blacksubgroup %>% mutate(BHgIQR = Hg/BHgiqr)

mod1<- lm( CE_BARS_CPT_D_PRIME ~ HgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ HgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_CPT_D_PRIME ~ WHgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WHgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_CPT_D_PRIME ~ BHgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BHgIQR + Cu + Mn + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# Mn
Mniqr<-IQR(data_log$Mn) 
WMniqr<-IQR(whitesubgroup$Mn) 
BMniqr<-IQR(blacksubgroup$Mn) 
data_log<-data_log %>% mutate(MnIQR = Mn/Mniqr)
whitesubgroup<-whitesubgroup %>% mutate(WMnIQR = Mn/WMniqr)
blacksubgroup<-blacksubgroup %>% mutate(BMnIQR = Mn/BMniqr)

mod1<- lm( CE_BARS_CPT_D_PRIME ~ MnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log)
mod1<- glm(CPTcr1 ~ MnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_CPT_D_PRIME ~ WMnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WMnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_CPT_D_PRIME ~ BMnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BMnIQR + Cu + Hg + Cr + Zn + Pb + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# Pb
Pbiqr<-IQR(data_log$Pb) 
WPbiqr<-IQR(whitesubgroup$Pb) 
BPbiqr<-IQR(blacksubgroup$Pb) 
data_log<-data_log %>% mutate(PbIQR = Pb/Pbiqr)
whitesubgroup<-whitesubgroup %>% mutate(WPbIQR = Pb/WPbiqr)
blacksubgroup<-blacksubgroup %>% mutate(BPbIQR = Pb/BPbiqr)

mod1<- lm( CE_BARS_CPT_D_PRIME ~ PbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log)
mod1<- glm(CPTcr1 ~ PbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_CPT_D_PRIME ~ WPbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WPbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_CPT_D_PRIME ~ BPbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BPbIQR + Cu + Hg + Cr + Zn + Mn + As + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# Se
Seiqr<-IQR(data_log$Se) 
WSeiqr<-IQR(whitesubgroup$Se) 
BSeiqr<-IQR(blacksubgroup$Se) 
data_log<-data_log %>% mutate(SeIQR = Se/Seiqr)
whitesubgroup<-whitesubgroup %>% mutate(WSeIQR = Se/WSeiqr)
blacksubgroup<-blacksubgroup %>% mutate(BSeIQR = Se/BSeiqr)

mod1<- lm( CE_BARS_CPT_D_PRIME ~ SeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ SeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_CPT_D_PRIME ~ WSeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WSeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_CPT_D_PRIME ~ BSeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BSeIQR + Cu + Hg + Cr + Zn + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

# Zn
Zniqr<-IQR(data_log$Zn) 
WZniqr<-IQR(whitesubgroup$Zn) 
BZniqr<-IQR(blacksubgroup$Zn) 
data_log<-data_log %>% mutate(ZnIQR = Zn/Zniqr)
whitesubgroup<-whitesubgroup %>% mutate(WZnIQR = Zn/WZniqr)
blacksubgroup<-blacksubgroup %>% mutate(BZnIQR = Zn/BZniqr)

mod1<- lm( CE_BARS_CPT_D_PRIME ~ ZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod1<- glm(CPTcr1 ~ ZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log, family = binomial)
summary(mod1)
confint(mod1)

mod2<- lm( CE_BARS_CPT_D_PRIME ~ WZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod2<- glm(CPTcr1 ~ WZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb+ CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup, family = binomial)
summary(mod2)
confint(mod2)

mod3<- lm( CE_BARS_CPT_D_PRIME ~ BZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod3<- glm(CPTcr1 ~ BZnIQR + Cu + Hg + Cr + Se + Mn + As + Pb + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup, family = binomial)
summary(mod3)
confint(mod3)

## dot whisker plot IQR scaled ---------------------------------------------
mod1<- lm(CE_BARS_DST_REVERSE_CNT ~ AsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod2<- lm(CE_BARS_DST_REVERSE_CNT ~ WAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_REVERSE_CNT ~ BAsIQR + Cu + Mn + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod4<- lm(CE_BARS_DST_REVERSE_CNT ~ CrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod5<- lm(CE_BARS_DST_REVERSE_CNT ~ WCrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod6<- lm(CE_BARS_DST_REVERSE_CNT ~ BCrIQR + Cu + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod7<- lm(CE_BARS_DST_REVERSE_CNT ~ CuIQR + Cr + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod8<- lm(CE_BARS_DST_REVERSE_CNT ~ WCuIQR + Cr + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod9<- lm(CE_BARS_DST_REVERSE_CNT ~ BCuIQR + Cr + As + Mn + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod10<- lm(CE_BARS_DST_REVERSE_CNT ~ HgIQR + Cu + As + Cr + Zn + Pb + Mn + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod11<- lm(CE_BARS_DST_REVERSE_CNT ~ WHgIQR + Cu + As + Cr + Zn + Pb + Mn + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod12<- lm(CE_BARS_DST_REVERSE_CNT ~ BHgIQR + Cu + As + Cr + Zn + Pb + Mn + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod13<- lm(CE_BARS_DST_REVERSE_CNT ~ MnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod14<- lm(CE_BARS_DST_REVERSE_CNT ~ WMnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod15<- lm(CE_BARS_DST_REVERSE_CNT ~ BMnIQR + Cu + As + Cr + Zn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod16<- lm(CE_BARS_DST_REVERSE_CNT ~ PbIQR + Cu + As + Cr + Zn + Mn + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod17<- lm(CE_BARS_DST_REVERSE_CNT ~ WPbIQR + Cu + As + Cr + Zn + Mn + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod18<- lm(CE_BARS_DST_REVERSE_CNT ~ BPbIQR + Cu + As + Cr + Zn + Mn + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod19<- lm(CE_BARS_DST_REVERSE_CNT ~ SeIQR + Cu + As + Cr + Zn + Pb + Hg + Mn + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod20<- lm(CE_BARS_DST_REVERSE_CNT ~ WSeIQR + Cu + As + Cr + Zn + Pb + Hg + Mn + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod21<- lm(CE_BARS_DST_REVERSE_CNT ~ BSeIQR + Cu + As + Cr + Zn + Pb + Hg + Mn + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)
mod22<- lm(CE_BARS_DST_REVERSE_CNT ~ ZnIQR + Cu + As + Cr + Mn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = data_log)
mod23<- lm(CE_BARS_DST_REVERSE_CNT ~ WZnIQR + Cu + As + Cr + Mn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = whitesubgroup)
mod24<- lm(CE_BARS_DST_REVERSE_CNT ~ BZnIQR + Cu + As + Cr + Mn + Pb + Hg + Se + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER +   marital, data = blacksubgroup)

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


# supplemental fig race differences in toenail metal conc
datalog<-data_log
mod1 <- lm(Mn ~ CE_AGE + THC_CUMULATIVE1 + EN_FORMERSMOKER + passivesmoke + race + CE_BMI + EN_STATE, data = datalog)
mod2 <- lm(Se ~ CE_AGE + THC_CUMULATIVE1 + EN_FORMERSMOKER + passivesmoke + race + CE_BMI + EN_STATE, data = datalog)
mod3 <- lm(As ~ CE_AGE + THC_CUMULATIVE1 + EN_FORMERSMOKER + passivesmoke + race + CE_BMI + EN_STATE, data = datalog)
mod4 <- lm(Cr ~ CE_AGE + THC_CUMULATIVE1 + EN_FORMERSMOKER + passivesmoke + race + CE_BMI + EN_STATE, data = datalog)
mod5 <- lm(Hg ~ CE_AGE + THC_CUMULATIVE1 + EN_FORMERSMOKER + passivesmoke + race + CE_BMI + EN_STATE, data = datalog)
mod6 <- lm(Cu ~ CE_AGE + THC_CUMULATIVE1 + EN_FORMERSMOKER + passivesmoke + race + CE_BMI + EN_STATE, data = datalog)
mod7 <- lm(Pb ~ CE_AGE + THC_CUMULATIVE1 + EN_FORMERSMOKER + passivesmoke + race + CE_BMI + EN_STATE, data = datalog)
mod8 <- lm(Zn ~ CE_AGE + THC_CUMULATIVE1 + EN_FORMERSMOKER + passivesmoke + race + CE_BMI + EN_STATE, data = datalog)

Mn <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term != "THC_CUMULATIVE1") %>% filter(term != "EN_FORMERSMOKER") %>% filter(term != "CE_C1") %>% filter(term != "EN_FORMERSMOKER1") %>% filter(term != "passivesmoke1")%>% filter(term != "passivesmoke")%>% filter(term != "raceOther")%>% filter(term != "CE_BMI")%>% filter(term != "EN_STATETX")%>% filter(term != "EN_STATEFL")%>% filter(term != "EN_STATELA")%>% filter(term != "EN_STATEMS")
Se <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term != "THC_CUMULATIVE1") %>% filter(term != "EN_FORMERSMOKER") %>% filter(term != "CE_C1") %>% filter(term != "EN_FORMERSMOKER1") %>% filter(term != "passivesmoke1")%>% filter(term != "passivesmoke")%>% filter(term != "raceOther")%>% filter(term != "CE_BMI")%>% filter(term != "EN_STATETX")%>% filter(term != "EN_STATEFL")%>% filter(term != "EN_STATELA")%>% filter(term != "EN_STATEMS")
As <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term != "THC_CUMULATIVE1") %>% filter(term != "EN_FORMERSMOKER") %>% filter(term != "CE_C1") %>% filter(term != "EN_FORMERSMOKER1") %>% filter(term != "passivesmoke1")%>% filter(term != "passivesmoke")%>% filter(term != "raceOther")%>% filter(term != "CE_BMI")%>% filter(term != "EN_STATETX")%>% filter(term != "EN_STATEFL")%>% filter(term != "EN_STATELA")%>% filter(term != "EN_STATEMS")
Cr <-broom::tidy(mod4) %>% filter(term != "CE_AGE") %>% filter(term != "THC_CUMULATIVE1") %>% filter(term != "EN_FORMERSMOKER") %>% filter(term != "CE_C1") %>% filter(term != "EN_FORMERSMOKER1") %>% filter(term != "passivesmoke1")%>% filter(term != "passivesmoke")%>% filter(term != "raceOther")%>% filter(term != "CE_BMI")%>% filter(term != "EN_STATETX")%>% filter(term != "EN_STATEFL")%>% filter(term != "EN_STATELA")%>% filter(term != "EN_STATEMS")
Hg <-broom::tidy(mod5) %>% filter(term != "CE_AGE") %>% filter(term != "THC_CUMULATIVE1") %>% filter(term != "EN_FORMERSMOKER") %>% filter(term != "CE_C1") %>% filter(term != "EN_FORMERSMOKER1") %>% filter(term != "passivesmoke1")%>% filter(term != "passivesmoke")%>% filter(term != "raceOther")%>% filter(term != "CE_BMI")%>% filter(term != "EN_STATETX")%>% filter(term != "EN_STATEFL")%>% filter(term != "EN_STATELA")%>% filter(term != "EN_STATEMS")
Cu <-broom::tidy(mod6) %>% filter(term != "CE_AGE") %>% filter(term != "THC_CUMULATIVE1") %>% filter(term != "EN_FORMERSMOKER") %>% filter(term != "CE_C1") %>% filter(term != "EN_FORMERSMOKER1") %>% filter(term != "passivesmoke1")%>% filter(term != "passivesmoke")%>% filter(term != "raceOther")%>% filter(term != "CE_BMI")%>% filter(term != "EN_STATETX")%>% filter(term != "EN_STATEFL")%>% filter(term != "EN_STATELA")%>% filter(term != "EN_STATEMS")
Pb<-broom::tidy(mod7) %>% filter(term != "CE_AGE") %>% filter(term != "THC_CUMULATIVE1") %>% filter(term != "EN_FORMERSMOKER") %>% filter(term != "CE_C1") %>% filter(term != "EN_FORMERSMOKER1") %>% filter(term != "passivesmoke1")%>% filter(term != "passivesmoke")%>% filter(term != "raceOther")%>% filter(term != "CE_BMI")%>% filter(term != "EN_STATETX")%>% filter(term != "EN_STATEFL")%>% filter(term != "EN_STATELA")%>% filter(term != "EN_STATEMS")
Zn<-broom::tidy(mod8) %>% filter(term != "CE_AGE") %>% filter(term != "THC_CUMULATIVE1") %>% filter(term != "EN_FORMERSMOKER") %>% filter(term != "CE_C1") %>% filter(term != "EN_FORMERSMOKER1") %>% filter(term != "passivesmoke1")%>% filter(term != "passivesmoke")%>% filter(term != "raceOther")%>% filter(term != "CE_BMI")%>% filter(term != "EN_STATETX")%>% filter(term != "EN_STATEFL")%>% filter(term != "EN_STATELA")%>% filter(term != "EN_STATEMS")


summary(mod1)
(10^(coef(mod8))-1)*100
(10^(confint(mod8))-1)*100

library(dotwhisker)
library(gdata)
library(dplyr)
essential_models <- combine(Zn, Se, Pb, Mn, Hg, Cu, Cr, As)
colnames(essential_models)[6] ="model"

plotRace_essential<- dwplot(essential_models, dodge_size = 1, dot_args = list(aes(color = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% relabel_predictors(c("raceBlack" = "1"))
plotRace_essential
plotRace_essential1 <- plotRace_essential + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + ggtitle("NEI Metals") + theme(text = element_text(size = 12)) +  labs(x = bquote('Percent difference compared to White')) + 
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) +  
  scale_x_continuous(breaks = log10(pretty(10^(essential_models$estimate),n=500)), labels = pretty(10^(essential_models$estimate),n=500)) + theme_bw() + theme(panel.border = element_rect(fill=NA, colour = "black", size=1))
plotRace_essential1 
