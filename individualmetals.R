setwd("/Users/joycelin/Desktop/Gulf/Aim3")
library(tidyverse)
library(readxl)
library(Hmisc) 
library(gtsummary)
library(summarytools)
library(ggplot2)

# tidy long data for facet wrap plot
data<- read_excel("CE_perfneuroscaled_tertiles.xlsx")
data_log <- data %>% mutate(across(c(Mg, Al, Ca, Cr, Mn, Fe, Ni, Cu, Zn, As, Se, Hg, Pb), log10))


## group for metals below 60% detect ----------------------
data_log$CE_BARS_MTS_COR_CNT <- factor(data_log$CE_BARS_MTS_COR_CNT)
data_log$CE_BARS_DST_REVERSE_CNT <- factor(data_log$CE_BARS_DST_REVERSE_CNT)
data_log$CE_BARS_DST_FORWARD_CNT <- factor(data_log$CE_BARS_DST_FORWARD_CNT)
data_log$CdBinary <- factor(data_log$CdBinary)
data_log$CoBinary <- factor(data_log$CoBinary)
data_log$MoBinary <- factor(data_log$MoBinary)
data_log$SbTertile <- factor(data_log$SbTertile)
data_log$VTertile <- factor(data_log$VTertile)
data_log$Batch <- factor(data_log$Batch)
data_log$EN_FORMERSMOKER<- factor(data_log$EN_FORMERSMOKER)
data_log$headinjury <- factor(data_log$headinjury)
data_log$headinjuryever <- factor(data_log$headinjuryever)
data_log$drinklot <- factor(data_log$drinklot)

data_log<- data_log %>% mutate(headinjuryever = case_when(headinjury == 2 ~ '1',
                                          headinjury ==0 ~ '0')) # end function
data_log<- data_log %>% mutate(drinklot = case_when(drinksweek <= 14 ~ '0',
                                                          drinksweek >14 ~ '1')) # end function



## for fully adjusted regressions ------------------------------------------------
data_log$EDU <- relevel(factor(data_log$EDU), ref = "College or more")
data_log$race <- relevel(factor(data_log$race), ref = "White")
data_log$EN_FORMERSMOKER <- as.factor(data_log$EN_FORMERSMOKER)
data_log$EN_MARITAL <- factor(data_log$EN_MARITAL)
data_log$quartileTHC <- factor(data_log$quartileTHC)

data_log <- data_log[-(414:429),]
data <- data[-(414:429),]

data_log$CE_V4B_TEST_A_TIME <- as.numeric(data_log$CE_V4B_TEST_A_TIME)

###plot relationship between covariates and outcomes
ggplot(data_log, aes(x=(CE_AGE), y=CE_BARS_MTS_COR_CNT)) + geom_point(size=1) + geom_smooth(method=loess)
ggplot(data_log, aes(x=(EDU), y=CE_BARS_MTS_COR_CNT)) + geom_boxplot(aes(fill=EDU))
ggplot(data_log, aes(x=(Batch), y=CE_BARS_MTS_COR_CNT)) + geom_boxplot(aes(fill=Batch))
ggplot(data_log, aes(x=(headinjuryever), y=CE_BARS_MTS_COR_CNT)) + geom_boxplot(aes(fill=headinjuryever))
ggplot(data_log, aes(x=(EN_FORMERSMOKER), y=CE_BARS_MTS_COR_CNT)) + geom_boxplot(aes(fill=EN_FORMERSMOKER))
ggplot(data_log, aes(x=(drinklot), y=CE_BARS_MTS_COR_CNT)) + geom_boxplot(aes(fill=drinklot))
## split by race
ggplot(data_log, aes(x=(drinklot), y=CE_V4B_TEST_A_TIME)) + geom_boxplot(aes(fill=race))
ggplot(data_log, aes(x=(EN_FORMERSMOKER), y=CE_V4B_TEST_A_TIME)) + geom_boxplot(aes(fill=race))
ggplot(data_log, aes(x=(headinjuryever), y=CE_V4B_TEST_A_TIME)) + geom_boxplot(aes(fill=race))
ggplot(data_log, aes(x=(EDU), y=CE_V4B_TEST_A_TIME)) + geom_boxplot(aes(fill=race))
ggplot(data_log, aes(x=(CE_AGE), y=CE_V4B_TEST_A_TIME)) + geom_boxplot(aes(fill=race))
ggplot(data_log, aes(x=(CE_AGE), y=CE_V4B_TEST_A_TIME, color = race)) + geom_point() + geom_smooth(method=lm) + labs(y = bquote('CE_V4B_TEST_A_TIME'), x = bquote('Age'))



#### GAM to visualize relationship
if (!require("pacman")) install.packages("pacman", INSTALL_opts = '--no-lock')
pacman::p_load(qgcomp, knitr, ggplot2, mgcv, ggcorr,plot)

gam1 <- gam(CE_V4B_TEST_A_TIME ~ s(As, bs="cr") + CE_AGE + EDU + Batch + headinjuryever +drinklot + EN_FORMERSMOKER, data = data_log)
summary(gam1)
plot(gam1)

gam1 <- gam(CE_BARS_MTS_COR_CNT ~ s(Hg, bs="cr") + CE_AGE + EDU + Batch + headinjuryever +drinklot + EN_FORMERSMOKER, data = data_log)
summary(gam1)
plot(gam1)

gam1 <- gam(CE_BARS_CPT_HIT_FRACTION ~ s(Mn, bs="cr") + CE_AGE + EDU + Batch + headinjuryever +drinklot + EN_FORMERSMOKER, data = data_log)
summary(gam1)
plot(gam1)

gam1 <- gam(CE_BARS_CPT_FA_FRACTION ~ s(Mn, bs="cr") + CE_AGE + EDU + Batch + headinjuryever +drinklot + EN_FORMERSMOKER, data = data_log)
summary(gam1)
plot(gam1)

gam1 <- gam(CE_BARS_CPT_COR_HIT_FRACTION ~ s(Mn, bs="cr") + CE_AGE + EDU + Batch + headinjuryever +drinklot + EN_FORMERSMOKER, data = data_log)
summary(gam1)
plot(gam1)

gam1 <- gam(CE_BARS_CPT_D_PRIME ~ s(Mn, bs="cr") + CE_AGE + EDU + Batch + headinjuryever +drinklot + EN_FORMERSMOKER, data = data_log)
summary(gam1)
plot(gam1)

gam1 <- gam(CE_BARS_DST_FORWARD_CNT ~ s(Pb, bs="cr") + CE_AGE + EDU + Batch + headinjuryever +drinklot + EN_FORMERSMOKER, data = data_log)
summary(gam1)
plot(gam1)

gam1 <- gam(CE_BARS_DST_REVERSE_CNT ~ s(Pb, bs="cr") + CE_AGE + EDU + Batch + headinjuryever +drinklot + EN_FORMERSMOKER, data = data_log)
summary(gam1)
plot(gam1)

## regular linear models and for toxic metals of interest only (Pb, Mn, Cd)

mod1 <- lm(CE_BARS_CPT_D_PRIME ~ CdBinary + CE_AGE + EDU + Batch + headinjuryever + drinklot + EN_FORMERSMOKER, data = data_log)
summary(mod1)
confint(mod1)
plot(mod1)

mod1 <- lm(CE_BARS_CPT_COR_HIT_FRACTION ~ CdBinary + CE_AGE + EDU + Batch + headinjuryever + drinklot + EN_FORMERSMOKER, data = data_log)
summary(mod1)
confint(mod1)
plot(mod1)

mod1 <- lm(CE_BARS_CPT_HIT_FRACTION ~ CdBinary + CE_AGE + EDU + Batch + headinjuryever + drinklot + EN_FORMERSMOKER, data = data_log)
summary(mod1)
confint(mod1)
plot(mod1)

mod1 <- lm(CE_BARS_CPT_FA_FRACTION ~ CdBinary + CE_AGE + EDU + Batch + headinjuryever + drinklot + EN_FORMERSMOKER, data = data_log)
summary(mod1)
plot(mod1)
AIC(mod1)
confint(mod1)
(exp(coef(mod1))-1)*100


mod1 <- lm(CE_V4B_TEST_A_TIME ~ CdBinary + CE_AGE + EDU + Batch+ headinjuryever + drinklot + EN_FORMERSMOKER, data = data_log)
summary(mod1)
plot(mod1)
confint(mod1)
(exp(coef(mod1))-1)*100

mod1 <- glm(CE_BARS_DST_REVERSE_CNT ~Mn + CE_AGE + EDU + Batch + headinjuryever + drinklot + EN_FORMERSMOKER, data = data_log, family = binomial)
summary(mod1)
plot(mod1)

mod1 <- glm(CE_BARS_DST_FORWARD_CNT ~Mn + CE_AGE + EDU + Batch + headinjuryever + drinklot + EN_FORMERSMOKER, data = data_log, family = binomial)
summary(mod1)
plot(mod1)

mod1 <- glm(CE_BARS_MTS_COR_CNT ~Pb + CE_AGE + EDU + Batch + headinjuryever + drinklot + EN_FORMERSMOKER, data = data_log, family = binomial)
summary(mod1)
plot(mod1)

## controlling for other metals
mod1 <- lm(CE_BARS_CPT_FA_FRACTION ~ Mn + Pb + CdBinary + CE_AGE + EDU + Batch, data = data_log)
summary(mod1)
plot(mod1)