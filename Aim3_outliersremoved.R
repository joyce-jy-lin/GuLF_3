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

#### GAM to visualize relationships ---------------------------------------------------
# visually assess for potential outliers
orighist <- data %>% ggplot(aes(x = As)) +geom_histogram(position="identity", bins = 50) 
orighist

library(gridExtra)
boxplot(data_log$Mg, main = "Boxplot")
qqnorm(data_log$Mg, main = "Normal Q-Q plot")

# first subset data to exclude extreme values 
#Rosners method, ESD many outliers test
library(EnvStats)
Asout<-rosnerTest(data_log$As, k = 70)$all.stats #55 outliers
Asout <- Asout%>% filter (Outlier==TRUE)
outids<-Asout$Obs.Num
As_data <- filter(data_log, !(row.names(data_log) %in% outids))
# characteristics of outliers
As_out<-filter(data_log,(row.names(data_log) %in% outids))
As_out %>%
  group_by(race) %>%
  dplyr::summarise(Count = n())

orighist <- data_log %>% ggplot(aes(x = As)) +geom_histogram(position="identity", bins = 50) 
orighist
noouthist <- As_data %>% ggplot(aes(x = As)) +geom_histogram(position="identity", bins = 50) 
noouthist

Crout<-rosnerTest(data_log$Cr, k = 70)$all.stats #5 outliers
Crout <- Crout%>% filter (Outlier==TRUE)
outids<-Crout$Obs.Num
Cr_data <- filter(data_log, !(row.names(data_log) %in% outids))
# characteristics of outliers
Cr_out<-filter(data_log,(row.names(data_log) %in% outids))
Cr_out %>%
  group_by(race) %>%
  dplyr::summarise(Count = n())

orighist <- data_log %>% ggplot(aes(x = Cr)) +geom_histogram(position="identity", bins = 50) 
orighist
noouthist <- Cr_data %>% ggplot(aes(x = Cr)) +geom_histogram(position="identity", bins = 50) 
noouthist

Cuout<-rosnerTest(data_log$Cu, k = 70)$all.stats #7 outliers
Cuout <- Cuout%>% filter (Outlier==TRUE)
outids<-Cuout$Obs.Num
Cu_data <- filter(data_log, !(row.names(data_log) %in% outids))
# characteristics of outliers
Cu_out<-filter(data_log,(row.names(data_log) %in% outids))
Cu_out %>%
  group_by(race) %>%
  dplyr::summarise(Count = n())
orighist <- data_log %>% ggplot(aes(x = Cu)) +geom_histogram(position="identity", bins = 50) 
orighist
noouthist <- Cu_data %>% ggplot(aes(x = Cu)) +geom_histogram(position="identity", bins = 50) 
noouthist


Feout<-rosnerTest(data_log$Fe, k = 70)$all.stats #6 outliers
Feout <- Feout%>% filter (Outlier==TRUE)
outids<-Feout$Obs.Num
Fe_data <- filter(data_log, !(row.names(data_log) %in% outids))
# characteristics of outliers
Fe_out<-filter(data_log,(row.names(data_log) %in% outids))
Fe_out %>%
  group_by(race) %>%
  dplyr::summarise(Count = n())
orighist <- data_log %>% ggplot(aes(x = Fe)) +geom_histogram(position="identity", bins = 50) 
orighist
noouthist <- Fe_data %>% ggplot(aes(x = Fe)) +geom_histogram(position="identity", bins = 50) 
noouthist


Hgout<-rosnerTest(data_log$Hg, k = 70)$all.stats #9 outliers
Hgout <- Hgout%>% filter (Outlier==TRUE)
outids<-Hgout$Obs.Num
Hg_data <- filter(data_log, !(row.names(data_log) %in% outids))
# characteristics of outliers
Hg_out<-filter(data_log,(row.names(data_log) %in% outids))
Hg_out %>%
  group_by(race) %>%
  dplyr::summarise(Count = n())
orighist <- data_log %>% ggplot(aes(x = Hg)) +geom_histogram(position="identity", bins = 50) 
orighist
noouthist <- Hg_data %>% ggplot(aes(x = Hg)) +geom_histogram(position="identity", bins = 50) 
noouthist


Mgout<-rosnerTest(data_log$Mg, k = 70)$all.stats #9 outliers
Mgout <- Mgout%>% filter (Outlier==TRUE)
outids<-Mgout$Obs.Num
Mg_data <- filter(data_log, !(row.names(data_log) %in% outids))
# characteristics of outliers
Mg_out<-filter(data_log,(row.names(data_log) %in% outids))
Mg_out %>%
  group_by(race) %>%
  dplyr::summarise(Count = n())
orighist <- data_log %>% ggplot(aes(x = Mg)) +geom_histogram(position="identity", bins = 50) 
orighist
noouthist <- Mg_data %>% ggplot(aes(x = Mg)) +geom_histogram(position="identity", bins = 50) 
noouthist


Mnout<-rosnerTest(data_log$Mn, k = 70)$all.stats #10 outliers
Mnout <- Mnout%>% filter (Outlier==TRUE)
outids<-Mnout$Obs.Num
Mn_data <- filter(data_log, !(row.names(data_log) %in% outids))
# characteristics of outliers
Mn_out<-filter(data_log,(row.names(data_log) %in% outids))
Mn_out %>%
  group_by(race) %>%
  dplyr::summarise(Count = n())
orighist <- data_log %>% ggplot(aes(x = Cu)) +geom_histogram(position="identity", bins = 50) 
orighist
noouthist <- Cu_data %>% ggplot(aes(x = Cu)) +geom_histogram(position="identity", bins = 50) 
noouthist



Pbout<-rosnerTest(data_log$Pb, k = 70)$all.stats #26 outliers
Pbout <- Pbout%>% filter (Outlier==TRUE)
outids<-Pbout$Obs.Num
Pb_data <- filter(data_log, !(row.names(data_log) %in% outids))
# characteristics of outliers
Pb_out<-filter(data_log,(row.names(data_log) %in% outids))
Pb_out %>%
  group_by(race) %>%
  dplyr::summarise(Count = n())
orighist <- data_log %>% ggplot(aes(x = Pb)) +geom_histogram(position="identity", bins = 50) 
orighist
noouthist <- Pb_data %>% ggplot(aes(x = Pb)) +geom_histogram(position="identity", bins = 50) 
noouthist


Seout<-rosnerTest(data_log$Se, k = 70)$all.stats #15 outliers
Seout <- Seout%>% filter (Outlier==TRUE)
outids<-Seout$Obs.Num
Se_data <- filter(data_log, !(row.names(data_log) %in% outids))
# characteristics of outliers
Se_out<-filter(data_log,(row.names(data_log) %in% outids))
Se_out %>%
  group_by(race) %>%
  dplyr::summarise(Count = n())
orighist <- data_log %>% ggplot(aes(x = Se)) +geom_histogram(position="identity", bins = 50) 
orighist
noouthist <- Se_data %>% ggplot(aes(x = Se)) +geom_histogram(position="identity", bins = 50) 
noouthist


Znout<-rosnerTest(data_log$Zn, k = 70)$all.stats #14 outliers
Znout <- Znout%>% filter (Outlier==TRUE)
outids<-Znout$Obs.Num
Zn_data <- filter(data_log, !(row.names(data_log) %in% outids))
# characteristics of outliers
Zn_out<-filter(data_log,(row.names(data_log) %in% outids))
Zn_out %>%
  group_by(race) %>%
  dplyr::summarise(Count = n())

orighist <- data_log %>% ggplot(aes(x = Zn)) +geom_histogram(position="identity", bins = 50) 
orighist
noouthist <- Zn_data %>% ggplot(aes(x = Zn)) +geom_histogram(position="identity", bins = 50) 
noouthist

# (5 and 95th %iles) -------
quantile(data_log$As, 0.95)
quantile(data_log$As, 0.05)

As_data<-data_log %>% filter(As<quantile(data_log$As, 0.975))%>% filter(As> quantile(data_log$As, 0.025))
Ashist <- As_data %>% ggplot(aes(x = As)) +geom_histogram(position="identity", bins = 50) 
Ashist

Al_data<-data_log %>% filter(Al<quantile(data_log$Al, 0.95))%>% filter(Al> quantile(data_log$Al, 0.05))
Pb_data<-data_log %>% filter(Pb<quantile(data_log$Pb, 0.95))%>% filter(Pb> quantile(data_log$Pb, 0.05))
Mn_data<-data_log %>% filter(Mn<quantile(data_log$Mn, 0.95))%>% filter(Mn> quantile(data_log$Mn, 0.05))
Hg_data<-data_log %>% filter(Hg<quantile(data_log$Hg, 0.95))%>% filter(Hg> quantile(data_log$Hg, 0.05))
Cr_data<-data_log %>% filter(Cr<quantile(data_log$Cr, 0.95))%>% filter(Cr> quantile(data_log$Cr, 0.05))
Zn_data<-data_log %>% filter(Zn<quantile(data_log$Zn, 0.95))%>% filter(Zn> quantile(data_log$Zn, 0.05))
Se_data<-data_log %>% filter(Se<quantile(data_log$Se, 0.95))%>% filter(Se> quantile(data_log$Se, 0.05))
Fe_data<-data_log %>% filter(Fe<quantile(data_log$Fe, 0.95))%>% filter(Fe> quantile(data_log$Fe, 0.05))
Cu_data<-data_log %>% filter(Cu<quantile(data_log$Cu, 0.95))%>% filter(Cu> quantile(data_log$Cu, 0.05))
Ca_data<-data_log %>% filter(Ca<quantile(data_log$Ca, 0.95))%>% filter(Ca> quantile(data_log$Ca, 0.05))

# GAMs
if (!require("pacman")) install.packages("pacman", INSTALL_opts = '--no-lock')
pacman::p_load(qgcomp, knitr, ggplot2, mgcv, ggcorr,plot)

gam1 <- gam(CE_BARS_CPT_D_PRIME ~ s(Mg, bs="cr") + CE_AGE + EDU +drinklot + EN_FORMERSMOKER + marital, data = Mg_data)
summary(gam1)
plot(gam1)

gam2 <- gam(CE_BARS_CPT_HIT_FRACTION ~s(Mg, bs="cr") + CE_AGE + EDU +drinklot + EN_FORMERSMOKER + marital,
            data = Mg_data,
            family = binomial,
            method = "REML")
summary(gam2)
plot(gam2, pages = 1)

gam2 <- gam(CE_BARS_CPT_CR_FRACTION ~s(Mg, bs="cr") + CE_AGE + EDU +drinklot + EN_FORMERSMOKER + marital,
                data = Mg_data,
                family = binomial,
                method = "REML")
summary(gam2)
plot(gam2, pages = 1)
#plot(gam2, pages = 1, trans = plogis, shade = TRUE, shade.col = "lightgreen", col = "purple")

gam1 <- gam(CE_BARS_DST_FORWARD_CNT ~ s(Mg, bs="cr") + CE_AGE + EDU +drinklot + EN_FORMERSMOKER + marital, data = Mg_data)
summary(gam1)
plot(gam1)

gam1 <- gam(CE_BARS_DST_REVERSE_CNT ~ s(Mg, bs="cr") + CE_AGE + EDU +drinklot + EN_FORMERSMOKER + marital, data = Mg_data)
summary(gam1)
plot(gam1)

## Reorganizing quartiles without outliers
# #Split data into quartiles by metal log transformed
library(data.table)
setDT(As_data)
setDT(Pb_data)
setDT(Mn_data)
setDT(Hg_data)
setDT(Cr_data)
setDT(Ca_data)
setDT(Cu_data)
setDT(Fe_data)
setDT(Se_data)
setDT(Zn_data)
setDT(Al_data)

As_data[ , AsQ := cut(As,
                   breaks=quantile(As,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]
Al_data[ , AlQ := cut(Al,
                   breaks=quantile(Al,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]

Ca_data[ , CaQ := cut(Ca,
                   breaks=quantile(Ca,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]

Cr_data[ , CrQ := cut(Cr,
                   breaks=quantile(Cr,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]

Cu_data[ , CuQ := cut(Cu,
                   breaks=quantile(Cu,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]

Fe_data[ , FeQ := cut(Fe,
                   breaks=quantile(Fe,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]

Pb_data[ , PbQ := cut(Pb,
                   breaks=quantile(Pb,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]

Mg_data[ , MgQ := cut(Mg,
                   breaks=quantile(Mg,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]


Mn_data[ , MnQ := cut(Mn,
                   breaks=quantile(Mn,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]

Hg_data[ , HgQ := cut(Hg,
                   breaks=quantile(Hg,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]

Ni_data[ , NiQ := cut(Ni,
                   breaks=quantile(Ni,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]

Se_data[ , SeQ := cut(Se,
                   breaks=quantile(Se,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]

Zn_data[ , ZnQ := cut(Zn,
                   breaks=quantile(Zn,
                   probs=seq(0, 1, by=0.25), na.rm=T),
                   include.lowest= TRUE, labels=1:4) ]

#### plot dot whisker for all, black, and white together ----------------------------------
library(dotwhisker)
library(dplyr)
library(gdata)

blacksubgroup<-filter(As_data, race =="Black")
whitesubgroup<-filter(As_data, race =="White")

mod1<- lm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = As_data)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

plotRace_As<- dwplot(all_models, dodge_size = 1, dot_args = list(aes(colour = model),size = 5), whisker_args = list(aes(colour=model), size = 3)) %>% relabel_predictors(c("raceBlack" = "Black", "AsQuartile2" = "As Q2", "AsQuartile3" = "As Q3", "AsQuartile4" = "As Q4"))
plotRace_As
plotRace_As1 <- plotRace_As + theme_bw(base_size = 4) + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + ggtitle("CPT D Prime") + theme(text = element_text(size = 16)) + labs(x = bquote('Coefficent Estimate with 95% CIs')) + theme(axis.text=element_text(size=17)) + theme(axis.title.x = element_text(margin = margin(t = 15))) + theme(panel.border = element_rect(fill=NA, colour = "black", size=1))  
plotRace_As1 

describe(whitesubgroup$CrQ)
describe(blacksubgroup$CrQ)

##plot black white all whisker plot for dichotomized outcomes log reg
#### plot dot whisker for all, black, and white together ----------------------------------

blacksubgroup<-filter(Pb_data, race =="Black")
whitesubgroup<-filter(Pb_data, race =="White")

mod1<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = Pb_data, family = binomial)
mod2<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, All)
colnames(all_models)[6] ="model"

plotRace_As<- dwplot(all_models, dodge_size = 1, dot_args = list(aes(colour = model),size = 5), whisker_args = list(aes(colour=model), size = 3)) %>% relabel_predictors(c("raceBlack" = "Black", "AsQuartile2" = "As Q2", "AsQuartile3" = "As Q3", "AsQuartile4" = "As Q4"))
plotRace_As
plotRace_As1 <- plotRace_As + theme_bw(base_size = 4) + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + ggtitle("CPT CR Fraction (dichotomous)") + theme(text = element_text(size = 16)) + labs(x = bquote('Coefficent Estimate with 95% CIs')) + theme(axis.text=element_text(size=17)) + theme(axis.title.x = element_text(margin = margin(t = 15))) + theme(panel.border = element_rect(fill=NA, colour = "black", size=1))  
plotRace_As1 

## outliers kept in for comparison---------------------------------
data_log$AsQ <- factor(data_log$AsQ)
data_log$CrQ <- factor(data_log$CrQ)

blacksubgroup<-filter(data_log, race =="Black")
whitesubgroup<-filter(data_log, race =="White")

mod1<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "Cr")%>% filter(term != "Pb")%>% filter(term != "Cu")%>% filter(term != "Hg")%>% filter(term != "As")%>% filter(term != "Mg")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "Cr")%>% filter(term != "Pb")%>% filter(term != "Cu")%>% filter(term != "Hg")%>% filter(term != "As")%>% filter(term != "Mg")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "Cr")%>% filter(term != "Pb")%>% filter(term != "Cu")%>% filter(term != "Hg")%>% filter(term != "As")%>% filter(term != "Mg")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

plotRace_As<- dwplot(all_models, dodge_size = 1, dot_args = list(aes(colour = model),size = 5), whisker_args = list(aes(colour=model), size = 3)) %>% relabel_predictors(c("raceBlack" = "Black", "AsQuartile2" = "As Q2", "AsQuartile3" = "As Q3", "AsQuartile4" = "As Q4"))
plotRace_As
plotRace_As1 <- plotRace_As + theme_bw(base_size = 4) + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + ggtitle("CPT D Prime") + theme(text = element_text(size = 16)) + labs(x = bquote('Coefficent Estimate with 95% CIs')) + theme(axis.text=element_text(size=17)) + theme(axis.title.x = element_text(margin = margin(t = 15))) + theme(panel.border = element_rect(fill=NA, colour = "black", size=1))  
plotRace_As1 



