setwd("/Users/joycelin/Desktop/Gulf/Aim3/")
library(tidyverse)
library(readxl)
library(Hmisc) 
library(gtsummary)
library(summarytools)
library(ggplot2)

# tidy long data for facet wrap plot
data<- read.csv("CE_neuroscaled_quantiles.csv")
data_log <- data %>% mutate(across(c(Mg, Al, Ca, Cr, Mn, Fe, Ni, Cu, Zn, As, Se, Hg, Pb), log10))

data_log$CdBinary <- factor(data_log$CdBinary)
data_log$CoBinary <- factor(data_log$CoBinary)
data_log$MoBinary <- factor(data_log$MoBinary)
data_log$SbTertile <- factor(data_log$SbTertile)
data_log$VTertile <- factor(data_log$VTertile)
data_log$Batch <- factor(data_log$Batch)
data_log$EN_FORMERSMOKER<- factor(data_log$EN_FORMERSMOKER)
data_log$passivesmoke<- factor(data_log$passivesmoke)
data_log$headinjury <- factor(data_log$headinjury)
data_log$headinjuryever <- factor(data_log$headinjuryever)
data_log$drinklot <- factor(data_log$drinklot)
data_log$EDU <- relevel(factor(data_log$EDU), ref = "College or more")
data_log$race <- relevel(factor(data_log$race), ref = "White")

data_log<- data_log %>% mutate(headinjuryever = case_when(headinjury == 2 ~ '1',
                                                          headinjury ==0 ~ '0',
                                                          headinjury ==1 ~ '0')) # head injury where loss consciousness
data_log<- data_log %>% mutate(drinklot = case_when(drinksweek <= 14 ~ '0',
                                                    drinksweek >14 ~ '1')) # more than 14 drinks/week
## plot descriptive
ggplot(data_log, aes(x=(Pb), y=CE_V4B_TEST_A_TIME)) + geom_point() + geom_smooth(method=lm) + labs(y = bquote('CE_V4B_TEST_A_TIME'), x = bquote('log10Pb'))
## probably need to be modeled as counts
hist(data_log$CE_BARS_CPT_CR_FRACTION)
hist(data_log$CE_BARS_CPT_FA_FRACTION)
hist(data_log$CE_BARS_CPT_HIT_FRACTION)
hist(data_log$CE_BARS_MTS_COR_CNT)

#can be modeled as linear regression
hist(data_log$CE_BARS_CPT_D_PRIME)
hist(data_log$CE_BARS_CPT_FA_LATENCY)
hist(data_log$CE_BARS_DST_REVERSE_CNT)

## try neg binomial model for count data where mean does not equal variance
require(MASS)
require(foreign)

#neg binomial
m1 <- glm.nb(CE_BARS_CPT_CR_FRACTION ~ Pb + CE_AGE + EDU + Batch + headinjuryever + drinklot + EN_FORMERSMOKER + passivesmoke, data = data_log)
summary(m1)
(est <- cbind(Estimate = coef(m1), confint(m1)))
exp(est)

#poisson
m3 <- glm(CE_BARS_CPT_CR_FRACTION ~ Pb + CE_AGE + EDU + Batch + headinjuryever + drinklot + EN_FORMERSMOKER + passivesmoke, family = "poisson", data = data_log)
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)

# regular linear regression
mod1 <- lm(CE_BARS_CPT_CR_FRACTION ~ Pb + CE_AGE + EDU + Batch + headinjuryever + drinklot + EN_FORMERSMOKER + passivesmoke, data = data_log)
summary(mod1)
confint(mod1)
plot(mod1)
AIC(mod1)
