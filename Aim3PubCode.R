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

data_log$AlQ <- factor(data_log$AlQ)
data_log$AsQ <- factor(data_log$AsQ)
data_log$CaQ <- factor(data_log$CaQ)
data_log$CrQ <- factor(data_log$CrQ)
data_log$CuQ <- factor(data_log$CuQ)
data_log$FeQ <- factor(data_log$FeQ)
data_log$PbQ <- factor(data_log$PbQ)
data_log$MgQ <- factor(data_log$MgQ)
data_log$MnQ <- factor(data_log$MnQ)
data_log$HgQ <- factor(data_log$HgQ)
data_log$NiQ <- factor(data_log$NiQ)
data_log$SeQ <- factor(data_log$SeQ)
data_log$ZnQ <- factor(data_log$ZnQ)
data_log$CdBinary <- factor(data_log$CdBinary)
data_log$CoBinary <- factor(data_log$CoBinary)
data_log$MoBinary <- factor(data_log$MoBinary)
data_log$SbTertile <- factor(data_log$SbTertile)
data_log$VTertile <- factor(data_log$VTertile)

data_log$AlQ <- ordered(data_log$AlQ)
data_log$AsQ <- ordered(data_log$AsQ)
data_log$CaQ <- ordered(data_log$CaQ)
data_log$CrQ <- ordered(data_log$CrQ)
data_log$CuQ <- ordered(data_log$CuQ)
data_log$FeQ <- ordered(data_log$FeQ)
data_log$PbQ <- ordered(data_log$PbQ)
data_log$MgQ <- ordered(data_log$MgQ)
data_log$MnQ <- ordered(data_log$MnQ)
data_log$HgQ <- ordered(data_log$HgQ)
data_log$NiQ <- ordered(data_log$NiQ)
data_log$SeQ <- ordered(data_log$SeQ)
data_log$ZnQ <- ordered(data_log$ZnQ)

data_log <- data_log %>% 
  rename("Race" = "race")

blacksubgroup<-filter(data_log, Race =="Black")
whitesubgroup<-filter(data_log, Race =="White")

## Overall regressions and outcomes for all metals ---------------------------------------
mod_CPTdprime <- lm(CE_BARS_CPT_D_PRIME~ MgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
summary(mod_CPTdprime)
confint(mod_CPTdprime)

mod_DSTfor <- lm(CE_BARS_DST_FORWARD_CNT~ MgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
summary(mod_DSTfor)
confint(mod_DSTfor)

mod_DSTrev <- lm(CE_BARS_DST_REVERSE_CNT~ MgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
summary(mod_DSTrev)
confint(mod_DSTrev)

mod_CPThit <- glm(CPThit1~ MgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family="binomial")
summary(mod_CPThit)
confint(mod_CPThit)

mod_CPTcr <- glm(CPTcr1~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family="binomial")
summary(mod_CPTcr)
confint(mod_CPTcr)

## Overall regressions and outcomes controlling for all other metals ---------------------------------------
mod_CPTcr <- glm(CPTcr1~ AsQ + MnQ + PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family="binomial")
summary(mod_CPTcr)
confint(mod_CPTcr)

## Comparing minimally to maximally adjusted models ---------------------------------------
mod_CPTdprime <- lm(CE_BARS_CPT_D_PRIME~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
summary(mod_CPTdprime)
confint(mod_CPTdprime)
plot(mod_CPTdprime)
BIC(mod_CPTdprime)

mod_CPTdprime_thc <- lm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital + THC_CUMULATIVE1, data = data_log)
summary(mod_CPTdprime_thc)
BIC(mod_CPTdprime_thc)

mod_CPTdprime_caf <- glm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital + Caffiene2h, data = data_log)
summary(mod_CPTdprime_caf)
BIC(mod_CPTdprime_caf)

mod_CPTdprime_head <- glm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital + headinjuryever, data = data_log)
summary(mod_CPTdprime_head)
BIC(mod_CPTdprime_head)

# show regression outcomes in table
library(sjPlot)
library(sjmisc)
library(sjlabelled)
tab_model(mod_CPTdprime, mod_CPTdprime_thc, mod_CPTdprime_caf, mod_CPTdprime_head)

## plot dot whisker for all, black, and white together ----------------------------------
library(dotwhisker)
library(dplyr)
library(gdata)

col <- c("All" ="#4DAF4A", "White"= "#377EB8", "Black" = "#E41A1C")

# CPT D Prime
mod1<- lm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

CPTdprime_As<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
                relabel_predictors(c( "AsQ2" = "As Q2", "AsQ3" = "As Q3", "AsQ4" = "As Q4"))
CPTdprime_As
CPTdprime_As <- CPTdprime_As + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
                ggtitle("D Prime") + theme(plot.title = element_text(face="bold", size=12)) +
                scale_color_brewer(palette="Set1") + theme(text = element_text(size = 12)) + 
                theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
                theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") 
CPTdprime_As 

# CPT CR
mod1<- glm(CPTcr1 ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPTcr1 ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPTcr1 ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, All)
colnames(all_models)[6] ="model"

CPTcr_As<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "AsQ2" = "As Q2", "AsQ3" = "As Q3", "AsQ4" = "As Q4"))
CPTcr_As
CPTcr_As <- CPTcr_As + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("CR Fraction") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_manual(values=c("#E41A1C","#4DAF4A")) + theme(text = element_text(size = 12)) + 
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") + theme(axis.text.y=element_blank()) 
CPTcr_As 

# CPT Hit
mod1<- glm(CPThit1 ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPThit1 ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPThit1 ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, All)
colnames(all_models)[6] ="model"

CPThit_As<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "AsQ2" = "As Q2", "AsQ3" = "As Q3", "AsQ4" = "As Q4"))
CPThit_As
CPThit_As <- CPThit_As + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("Hit Fraction") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_manual(values=c("#E41A1C","#4DAF4A")) + theme(text = element_text(size = 12)) + 
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") + theme(axis.text.y=element_blank()) 
CPThit_As 

# DST forward
mod1<- lm(CE_BARS_DST_FORWARD_CNT ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_DST_FORWARD_CNT ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_FORWARD_CNT ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

DSTfor_As<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "AsQ2" = "As Q2", "AsQ3" = "As Q3", "AsQ4" = "As Q4"))
DSTfor_As
DSTfor_As <- DSTfor_As + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("Forward") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_brewer(palette="Set1") + theme(text = element_text(size = 12)) +
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") + theme(axis.text.y=element_blank()) 
DSTfor_As

# DST reverse
mod1<- lm(CE_BARS_DST_REVERSE_CNT ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_DST_REVERSE_CNT ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_REVERSE_CNT ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

DSTrev_As<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "AsQ2" = "As Q2", "AsQ3" = "As Q3", "AsQ4" = "As Q4"))
DSTrev_As
DSTrev_As <- DSTrev_As + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("Reverse") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_brewer(palette="Set1") + theme(text = element_text(size = 12)) +
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(axis.text.y=element_blank()) + labs(color=NULL) 
DSTrev_As

# combine plots
library(patchwork) 
Asplot<- CPTdprime_As + CPTfa_As + CPTcr_As + CPThit_As + DSTfor_As + DSTrev_As + plot_layout(ncol=6)
Asplot
ggsave("asplot.png", Asplot, bg='transparent')
png("asplot.png", width = 15, height = 5, units = 'in', res = 400) 
Asplot
dev.off()

## compare groupings of metals by race by quartile ------------------------------------
othersubgroup<-filter(data, race =="Other")
ggsave("Qhist.png", Qhist, bg='transparent')
png("Qhist.png", width = 8, height = 5, units = 'in', res = 400) 
Qhist
dev.off()

# Plots subset for main text --------------------------------

library(dotwhisker)
library(dplyr)
library(gdata)

col <- c("All" = "#4DAF4A", "White"="#377EB8", "Black" = "#E41A1C")
# CPT D Prime -----------------------------------------------------------------------
#As
mod1<- lm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

CPTdprime_As<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "AsQ2" = "Q2", "AsQ3" = "Q3", "AsQ4" = "Q4"))
CPTdprime_As
CPTdprime_As <- CPTdprime_As + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("As") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_brewer(palette="Set1") + theme(text = element_text(size = 12)) + 
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") 
CPTdprime_As 

# Cr
mod1<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

CPTdprime_Cr<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "CrQ2" = "Q2", "CrQ3" = "Q3", "CrQ4" = "Q4"))
CPTdprime_Cr
CPTdprime_Cr <- CPTdprime_Cr + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("Cr") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_brewer(palette="Set1") + theme(text = element_text(size = 12)) +
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") + theme(axis.text.y=element_blank()) 
CPTdprime_Cr

# Hg
mod1<- lm(CE_BARS_CPT_D_PRIME ~ HgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ HgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ HgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

CPTdprime_Hg<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "HgQ2" = "Q2", "HgQ3" = "Q3", "HgQ4" = "Q4"))
CPTdprime_Hg
CPTdprime_Hg <- CPTdprime_Hg + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("Hg") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_brewer(palette="Set1") + theme(text = element_text(size = 12)) +
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") + theme(axis.text.y=element_blank()) 
CPTdprime_Hg

# combine plots
library(patchwork) 
Dprime<- CPTdprime_As + CPTdprime_Cr + CPTdprime_Hg + plot_layout(ncol=3)
Dprime
ggsave("Dprime.png", Dprime, bg='transparent')
png("Dprime.png", width = 7, height = 5, units = 'in', res = 400) 
Dprime
dev.off()

## CR Fraction plots --------------------------------------------------
#As
mod1<- glm(CPTcr1 ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPTcr1 ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPTcr1 ~ AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, All)
colnames(all_models)[6] ="model"

CPTcr_As<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "AsQ2" = "Q2", "AsQ3" = "Q3", "AsQ4" = "Q4"))
CPTcr_As
CPTcr_As <- CPTcr_As + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("As") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_manual(values=c("#E41A1C","#4DAF4A")) + theme(text = element_text(size = 12)) + 
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") 
CPTcr_As 

# Pb
mod1<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, All)
colnames(all_models)[6] ="model"

CPTcr_Pb<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "PbQ2" = "Q2", "PbQ3" = "Q3", "PbQ4" = "Q4"))
CPTcr_Pb
CPTcr_Pb <- CPTcr_Pb + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("Pb") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_manual(values=c("#E41A1C","#4DAF4A")) + theme(text = element_text(size = 12)) + 
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") + theme(axis.text.y=element_blank()) 
CPTcr_Pb

# combine plots
library(patchwork) 
CRfrac<- CPTcr_As + CPTcr_Pb + plot_layout(ncol=2)
CRfrac
ggsave("CRfrac.png", CRfrac, bg='transparent')
png("CRfrac.png", width = 4.7, height = 5, units = 'in', res = 400) 
CRfrac
dev.off()

## Hit Fraction plots --------------------------------------------------
#Cr
mod1<- glm(CPThit1 ~ CrQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPThit1 ~ CrQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPThit1 ~ CrQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, All)
colnames(all_models)[6] ="model"

CPThit_Cr<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "CrQ2" = "Q2", "CrQ3" = "Q3", "CrQ4" = "Q4"))
CPThit_Cr
CPThit_Cr <- CPThit_Cr + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("Cr") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_manual(values=c("#E41A1C","#4DAF4A")) + theme(text = element_text(size = 12)) + 
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") 
CPThit_Cr 

ggsave("Hitfrac.png", CPThit_Cr, bg='transparent')
png("Hitfrac.png", width = 2.3333, height = 5, units = 'in', res = 400) 
CPThit_Cr
dev.off()


## DST forward -----------------------------------------------------------
#Mg
mod1<- lm(CE_BARS_DST_FORWARD_CNT ~ MgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_DST_FORWARD_CNT ~ MgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_FORWARD_CNT ~ MgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

DSTf_Mg<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "MgQ2" = "Q2", "MgQ3" = "Q3", "MgQ4" = "Q4"))
DSTf_Mg
DSTf_Mg <- DSTf_Mg + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("Mg") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_brewer(palette="Set1") + theme(text = element_text(size = 12)) + 
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") 
DSTf_Mg 

ggsave("DSTf.png", DSTf_Mg, bg='transparent')
png("DSTf.png", width = 2.3333, height = 5, units = 'in', res = 400) 
DSTf_Mg
dev.off()

## DST reverse --------------------------------------------------------
#Mg
mod1<- lm(CE_BARS_DST_REVERSE_CNT ~ MgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_DST_REVERSE_CNT ~ MgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_REVERSE_CNT ~ MgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

DSTr_Mg<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "MgQ2" = "Q2", "MgQ3" = "Q3", "MgQ4" = "Q4"))
DSTr_Mg
DSTr_Mg <- DSTr_Mg + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("Mg") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_brewer(palette="Set1") + theme(text = element_text(size = 12)) + 
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5)) + theme(legend.position = "none") 
DSTr_Mg 

#Mn
mod1<- lm(CE_BARS_DST_REVERSE_CNT ~ MnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_DST_REVERSE_CNT ~ MnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_REVERSE_CNT ~ MnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.35, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% 
  relabel_predictors(c( "MnQ2" = "Q2", "MnQ3" = "Q3", "MnQ4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw(base_size = 4) +
  ggtitle("Mn") + theme(plot.title = element_text(face="bold", size=12)) +
  scale_color_brewer(palette="Set1") + theme(text = element_text(size = 12)) + 
  theme(axis.text=element_text(size=12)) + theme(axis.title.x = element_text(margin = margin(t = 12))) + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5))+ theme(axis.text.y=element_blank()) 
DSTr_Mn 

# combine plots
library(patchwork) 
DSTr<- DSTr_Mg + DSTr_Mn + plot_layout(ncol=2)
DSTr
ggsave("DSTr.png", DSTr, bg='transparent')
png("DSTr.png", width = 4.7, height = 5, units = 'in', res = 400) 
DSTr
dev.off()

### Median regression plots ---------------------------------
# CPT D Prime
library(quantreg)
mod1 <- rq(CE_BARS_CPT_D_PRIME ~ HgQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
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

wmod1 <- rq(CE_BARS_CPT_D_PRIME ~ HgQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
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

bmod1 <- rq(CE_BARS_CPT_D_PRIME ~ HgQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)
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
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0, position=position_dodge(width=0.5), linewidth=1)+
  geom_hline(yintercept=0, colour = "grey60", linetype = 2) + scale_fill_discrete(breaks=c('All', 'White', 'Black'))+
  scale_y_continuous(limits=c(-0.5, 1),breaks=c(-0.5, 0, 0.5, 1)) + coord_flip() + 
  theme(axis.title.y = element_text(size = 16)) + labs(y = bquote('Difference compared to Q1')) + theme_bw()+ 
  theme(panel.border = element_rect(fill=NA, colour = "black", linewidth=1)) + 
  theme(axis.text=element_text(size=12)) + theme(legend.text = element_text(size = 12), legend.title = element_text(size = 12))
p

## DST Forward
library(quantreg)
mod1 <- rq(CE_BARS_DST_FORWARD_CNT ~ SeQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
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

wmod1 <- rq(CE_BARS_DST_FORWARD_CNT ~ SeQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
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

bmod1 <- rq(CE_BARS_DST_FORWARD_CNT ~ SeQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)
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


Se <- data.frame(x =c(coef["SeQ2"], coef["SeQ3"], coef["SeQ4"], wcoef["SeQ2"], wcoef["SeQ3"], wcoef["SeQ4"], bcoef["SeQ2"], bcoef["SeQ3"], bcoef["SeQ4"]),
                 y = c("Se Q2", "Se Q3", "Se Q4","Se Q2", "Se Q3", "Se Q4", "Se Q2", "Se Q3", "Se Q4"),
                 lower = c(ci_dat[2,1], ci_dat[3,1], ci_dat[4,1], wci_dat[2,1], wci_dat[3,1], wci_dat[4,1], bci_dat[2,1], bci_dat[3,1], bci_dat[4,1]),
                 upper = c(ci_dat[2,2], ci_dat[3,2], ci_dat[4,2], wci_dat[2,2], wci_dat[3,2], wci_dat[4,2], bci_dat[2,2], bci_dat[3,2], bci_dat[4,2]),
                 Race = c("All", "All", "All", "White", "White", "White","Black", "Black","Black"))

Se$y <- factor(Se$y, c("Se Q4", "Se Q3", "Se Q2"))
Se$Race <- factor(Se$Race, c("Black", "White", "All"))

library("ggplot2")
p<-ggplot(Se, aes(y, x, color=Race, group=Race)) +     
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0, position=position_dodge(width=0.5), linewidth=1)+
  geom_hline(yintercept=0, colour = "grey60", linetype = 2) + scale_fill_discrete(breaks=c('All', 'White', 'Black'))+
  scale_y_continuous(limits=c(-1.65, 1.65),breaks=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) + coord_flip() + 
  theme(axis.title.y = element_text(size = 16)) + labs(y = bquote('Difference compared to Q1')) + theme_bw()+ 
  theme(panel.border = element_rect(fill=NA, colour = "black", linewidth=1)) + 
  theme(axis.text=element_text(size=12)) + theme(legend.text = element_text(size = 12), legend.title = element_text(size = 12))
p


## DST Reverse
library(quantreg)
mod1 <- rq(CE_BARS_DST_REVERSE_CNT ~ ZnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
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

wmod1 <- rq(CE_BARS_DST_REVERSE_CNT ~ ZnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
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

bmod1 <- rq(CE_BARS_DST_REVERSE_CNT ~ ZnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)
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


Zn <- data.frame(x =c(coef["ZnQ2"], coef["ZnQ3"], coef["ZnQ4"], wcoef["ZnQ2"], wcoef["ZnQ3"], wcoef["ZnQ4"], bcoef["ZnQ2"], bcoef["ZnQ3"], bcoef["ZnQ4"]),
                 y = c("Zn Q2", "Zn Q3", "Zn Q4","Zn Q2", "Zn Q3", "Zn Q4", "Zn Q2", "Zn Q3", "Zn Q4"),
                 lower = c(ci_dat[2,1], ci_dat[3,1], ci_dat[4,1], wci_dat[2,1], wci_dat[3,1], wci_dat[4,1], bci_dat[2,1], bci_dat[3,1], bci_dat[4,1]),
                 upper = c(ci_dat[2,2], ci_dat[3,2], ci_dat[4,2], wci_dat[2,2], wci_dat[3,2], wci_dat[4,2], bci_dat[2,2], bci_dat[3,2], bci_dat[4,2]),
                 Race = c("All", "All", "All", "White", "White", "White","Black", "Black","Black"))

Zn$y <- factor(Zn$y, c("Zn Q4", "Zn Q3", "Zn Q2"))
Zn$Race <- factor(Zn$Race, c("Black", "White", "All"))

library("ggplot2")
p<-ggplot(Zn, aes(y, x, color=Race, group=Race)) +     
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0, position=position_dodge(width=0.5), linewidth=1)+
  geom_hline(yintercept=0, colour = "grey60", linetype = 2) + scale_fill_discrete(breaks=c('All', 'White', 'Black'))+
  scale_y_continuous(limits=c(-1.5, 1.25),breaks=c(-1.5, -1, -0.5, 0, 0.5, 1)) + coord_flip() + 
  theme(axis.title.y = element_text(size = 16)) + labs(y = bquote('Difference compared to Q1')) + theme_bw()+ 
  theme(panel.border = element_rect(fill=NA, colour = "black", linewidth=1)) + 
  theme(axis.text=element_text(size=12)) + theme(legend.text = element_text(size = 12), legend.title = element_text(size = 12))
p

## CPT Hit
library(dotwhisker)
library(dplyr)
library(gdata)

mod1<- glm(CPThit1 ~ CdBinary + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPThit1 ~ CdBinary + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPThit1 ~ CdBinary + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "Cr")%>% filter(term != "Pb")%>% filter(term != "Cu")%>% filter(term != "Hg")%>% filter(term != "As")%>% filter(term != "Mg")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "Cr")%>% filter(term != "Pb")%>% filter(term != "Cu")%>% filter(term != "Hg")%>% filter(term != "As")%>% filter(term != "Mg")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "Cr")%>% filter(term != "Pb")%>% filter(term != "Cu")%>% filter(term != "Hg")%>% filter(term != "As")%>% filter(term != "Mg")

all_models <- combine(Black, White,All)
colnames(all_models)[6] ="model"

plotRace_As<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% relabel_predictors(c("raceBlack" = "Black", "AsQuartile2" = "As Q2", "AsQuartile3" = "As Q3", "AsQuartile4" = "As Q4"))
plotRace_As
plotRace_As1 <- plotRace_As + theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  ggtitle("CPT Hit Fraction (dichotomous)") + theme(text = element_text(size = 12)) 
  labs(x = bquote('Log odds compared to Q1 ')) + theme(axis.text=element_text(size=12)) + scale_x_continuous(limits=c(-4.1, 3.2),breaks=c(-4, -2, 0, 2))+
  theme(axis.title.x = element_text(margin = margin(t = 12))) + theme(panel.border = element_rect(fill=NA, colour = "black", size=1))  
plotRace_As1 

## CPT Correct Response
col <- c("All" = "#619CFF", "White"="#00BA38", "Black" = "#F8766D")

mod1<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "Cr")%>% filter(term != "Pb")%>% filter(term != "Cu")%>% filter(term != "Hg")%>% filter(term != "As")%>% filter(term != "Mg")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "Cr")%>% filter(term != "Pb")%>% filter(term != "Cu")%>% filter(term != "Hg")%>% filter(term != "As")%>% filter(term != "Mg")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "Cr")%>% filter(term != "Pb")%>% filter(term != "Cu")%>% filter(term != "Hg")%>% filter(term != "As")%>% filter(term != "Mg")

all_models <- combine(Black, All)
colnames(all_models)[6] ="model"

plotRace_As<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 2), whisker_args = list(aes(colour=model), size = 1)) %>% relabel_predictors(c("raceBlack" = "Black", "AsQuartile2" = "As Q2", "AsQuartile3" = "As Q3", "AsQuartile4" = "As Q4"))
plotRace_As
plotRace_As1 <- plotRace_As + theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
 theme(text = element_text(size = 12)) + scale_color_manual(values=c("#F8766D","#00BA38"))+labs(x = bquote('Log odds compared to Q1 ')) + theme(axis.text=element_text(size=12)) + 
  scale_x_continuous(limits=c(-4.6, 3.2),breaks=c(-4, -2, 0, 2))+theme(axis.title.x = element_text(margin = margin(t = 12))) + theme(panel.border = element_rect(fill=NA, colour = "black", size=1))  
plotRace_As1 
