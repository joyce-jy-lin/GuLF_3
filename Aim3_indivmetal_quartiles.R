setwd("/Users/joycelin/Desktop/Gulf/Aim3")
library(tidyverse)
library(dplyr)
library(readxl)
library(Hmisc) 
library(gtsummary)
library(summarytools)
library(ggplot2)


## descriptive plots metals and neuro by quartile of metal -----------------------------------------------

data<- read_excel("CE_neuroscaled_quantiles.xlsx")

data$AlQuartile <- factor(data$AlQuartile)
data$AsQuartile <- factor(data$AsQuartile)
data$CaQuartile <- factor(data$CaQuartile)
data$CrQuartile <- factor(data$CrQuartile)
data$CuQuartile <- factor(data$CuQuartile)
data$FeQuartile <- factor(data$FeQuartile)
data$PbQuartile <- factor(data$PbQuartile)
data$MgQuartile <- factor(data$MgQuartile)
data$MnQuartile <- factor(data$MnQuartile)
data$HgQuartile <- factor(data$HgQuartile)
data$NiQuartile <- factor(data$NiQuartile)
data$SeQuartile <- factor(data$SeQuartile)
data$ZnQuartile <- factor(data$ZnQuartile)
data$CdBinary <- factor(data$CdBinary)
data$CoBinary <- factor(data$CoBinary)
data$MoBinary <- factor(data$MoBinary)
data$SbTertile <- factor(data$SbTertile)
data$VTertile <- factor(data$VTertile)
data$race <- relevel(factor(data$race), ref = "White")
data$quartileTHC <- factor(data$quartileTHC)
data$EN_FORMERSMOKER <- factor(data$EN_FORMERSMOKER)
data$CE_BARS_TAP_RIGHT_PREF_AVG <- as.numeric(data$CE_BARS_TAP_RIGHT_PREF_AVG)
data$CE_V4B_TEST_A_TIME <- as.numeric(data$CE_V4B_TEST_A_TIME)
data$CE_V4B_TEST_B_TIME <- as.numeric(data$CE_V4B_TEST_B_TIME)


data<- data %>% mutate(headinjuryever = case_when(headinjury == 2 ~ '1',
                                                          headinjury ==0 ~ '0',
                                                          headinjury ==1 ~ '0')) # head injury where loss consciousness
data<- data %>% mutate(drinklot = case_when(drinksweek <= 14 ~ '0',
                                                    drinksweek >14 ~ '1')) # more than 14 drinks/week

data$headinjuryever <- factor(data$headinjuryever)
data$drinklot <- factor(data$drinklot)

data<- data %>% mutate(marital = case_when(EN_MARITAL == 1 ~ 'Married or with partner',
                                              EN_MARITAL == 6 ~ 'Married or with partner',
                                              EN_MARITAL == 2 ~ 'Separated or single',
                                              EN_MARITAL == 3 ~ 'Separated or single',
                                              EN_MARITAL == 5 ~ 'Separated or single')) 

data$marital <- factor(data$marital)

data_log <- data %>% mutate(across(c(Mg, Al, Ca, Cr, Mn, Fe, Ni, Cu, Zn, As, Se, Hg, Pb), log10))

## plot descriptive
library(dotwhisker)

## Tests to include in aim 3 main text
## Attention
hist(data_log$CE_V4B_TEST_A_TIME)
hist(data_log$CE_BARS_CPT_HIT_FRACTION)
hist(data_log$CE_BARS_CPT_CR_FRACTION)
hist(data_log$CE_BARS_CPT_FA_FRACTION)
hist(data_log$CE_BARS_CPT_FA_LATENCY)
hist(data_log$CE_BARS_CPT_D_PRIME)
##Memory
hist(data_log$CE_BARS_DST_FORWARD_CNT)
hist(data_log$CE_BARS_MTS_COR_CNT)

# regular linear regression all data
mod1 <- lm(CE_BARS_CPT_FA_FRACTION ~ AsQuartile + CE_AGE + EDU + EN_FORMERSMOKER + headinjuryever + drinklot + marital + CE_BMI, data = data)
summary(mod1)
confint(mod1)
plot(mod1)
AIC(mod1)

#regular linear regression race subset
blacksubgroup<-filter(data, race =="Black")
whitesubgroup<-filter(data, race =="White")

mod1 <- lm(CE_BARS_CPT_FA_FRACTION ~ AsQuartile + CE_AGE + EDU + EN_FORMERSMOKER + headinjuryever + drinklot + marital + CE_BMI, data = blacksubgroup)
summary(mod1)
confint(mod1)
plot(mod1)
AIC(mod1)

# start plots switch out neuro test using cntrlf -------------------------------------------------------

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ CdBinary + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital, data = data)
Cd <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")%>%
  filter(term != "maritalSeparated or single")
plotmod<- dwplot(Cd)
plotCd<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotCd

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ AsQuartile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital + CE_BMI, data = blacksubgroup)
As <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "passivesmoke1")%>%
  filter(term != "passivesmokeNA")%>%
  filter(term != "CE_BMI")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")%>%
  filter(term != "maritalSeparated or single")
plotmod<- dwplot(As)
plotAs<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotAs

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ CuQuartile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot, data = data)
Cu <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")

plotmod<- dwplot(Cu)
plotCu<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotCu

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ FeQuartile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot, data = data)
Fe <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")
plotmod<- dwplot(Fe)
plotFe<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotFe

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ PbQuartile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital + CE_BMI, data = data)
Pb <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "passivesmoke1")%>%
  filter(term != "passivesmokeNA")%>%
  filter(term != "CE_BMI")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")%>%
  filter(term != "maritalSeparated or single")
plotmod<- dwplot(Pb)
plotPb<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotPb

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ MnQuartile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital, data = data)
Mn <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "passivesmoke1")%>%
  filter(term != "passivesmokeNA")%>%
  filter(term != "CE_BMI")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")%>%
  filter(term != "maritalSeparated or single")
plotmod<- dwplot(Mn)
plotMn<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotMn

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ HgQuartile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital, data = data)
Hg <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "passivesmoke1")%>%
  filter(term != "passivesmokeNA")%>%
  filter(term != "CE_BMI")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")%>%
  filter(term != "maritalSeparated or single")
plotmod<- dwplot(Hg)
plotHg<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotHg


mod1<- lm( CE_BARS_CPT_FA_FRACTION~ SeQuartile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital, data = data)
Se <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")%>%
  filter(term != "maritalSeparated or single")
plotmod<- dwplot(Se)
plotSe<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotSe

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ ZnQuartile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital, data = data)
Zn <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")%>%
  filter(term != "maritalSeparated or single")

plotmod<- dwplot(Zn)
plotZn<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotZn


library(patchwork)
 CE_BARS_CPT_FA_FRACTION<-(plotCd| plotAs|plotZn)/
(plotCu| plotFe| plotPb)/
(plotSe| plotMn| plotHg)


 CE_BARS_CPT_FA_FRACTION<-  CE_BARS_CPT_FA_FRACTION+ plot_annotation(title = 'Vibrotactile Hold')
 CE_BARS_CPT_FA_FRACTION

png(" CE_BARS_CPT_FA_FRACTION.png", width = 6, height = 10, units = 'in', res = 300) 
 CE_BARS_CPT_FA_FRACTION
dev.off() 

#### other metals not in main text
mod1<- lm( CE_BARS_CPT_FA_FRACTION~ MoBinary + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital, data = data)
Mo<-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")%>%
  filter(term != "maritalSeparated or single")

plotmod<- dwplot(Mo)
plotMo<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotMo

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ SbTertile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital, data = data)
Sb <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")
plotmod<- dwplot(Sb)
plotSb<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2)  + theme(legend.position = "none") 
plotSb

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ VTertile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot, data = data)
V <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")
plotmod<- dwplot(V)
plotV<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotV

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ MgQuartile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital, data = data)
Mg <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")%>%
  filter(term != "maritalSeparated or single")

plotmod<- dwplot(Mg)
plotMg<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotMg


mod1<- lm( CE_BARS_CPT_FA_FRACTION~ CaQuartile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot, data = data)
Ca <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")
plotmod<- dwplot(Ca)
plotCa<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2)  + theme(legend.position = "none") 
plotCa

mod1<- lm( CE_BARS_CPT_FA_FRACTION~ CrQuartile + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot, data = data)
Cr <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")

plotmod<- dwplot(Cr)
plotCr<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2)  + theme(legend.position = "none") 
plotCr


mod1<- lm( CE_BARS_CPT_FA_FRACTION~ CoBinary + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital, data = data)
Co <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% 
  filter(term != "Batch") %>% 
  filter(term != "Mass") %>% 
  filter(term !="EDUHighschool") %>% 
  filter(term !="EDULess than Highschool") %>% 
  filter(term !="EDUSome College") %>% 
  filter(term != "quartileTHC2")%>%
  filter(term != "quartileTHC3")%>%
  filter(term != "quartileTHC4")%>%
  filter(term != "quartileTHCNA")%>%
  filter(term != "EN_FORMERSMOKER")%>%
  filter(term != "EN_FORMERSMOKER1")%>%
  filter(term != "passivesmoke")%>%
  filter(term != "headinjuryever1")%>%
  filter(term != "drinklot1")%>%
  filter(term != "maritalSeparated or single")

plotmod<- dwplot(Co)
plotCo<- plotmod+ theme_bw() + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme(legend.position = "none") 
plotCo
                    
           
library(patchwork)
CE_BARS_CPT_FA_FRACTION<-(plotCo| plotMo|plotSb|plotV)/
  (plotZn | plotAl| plotCa)/ 
  (plotMg| plotCr|plotNi)

CE_BARS_CPT_FA_FRACTION<-  CE_BARS_CPT_FA_FRACTION+ plot_annotation(title = 'Vibrotactile Hold')
CE_BARS_CPT_FA_FRACTION

png(" CE_BARS_CPT_FA_FRACTION.png", width = 6, height = 10, units = 'in', res = 300) 
CE_BARS_CPT_FA_FRACTION
dev.off()          





#### plot dot whisker for all, black, and white together ----------------------------------
mod1<- lm( CE_BARS_CPT_FA_FRACTION~ CdBinary + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital + CE_BMI, data = data)
mod2<- lm( CE_BARS_CPT_FA_FRACTION~ CdBinary + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital + CE_BMI, data = whitesubgroup)
mod3<- lm( CE_BARS_CPT_FA_FRACTION~ CdBinary + CE_AGE + Batch + EDU + EN_FORMERSMOKER + passivesmoke+ headinjuryever + drinklot + marital + CE_BMI, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term != "Batch") %>% filter(term != "Mass") %>% filter(term !="EDUHighschool") %>% filter(term !="EDULess than Highschool") %>% filter(term !="EDUSome College") %>% filter(term != "quartileTHC2")%>%filter(term != "quartileTHC3")%>%filter(term != "quartileTHC4")%>% filter(term != "quartileTHCNA")%>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term != "Batch") %>% filter(term != "Mass") %>% filter(term !="EDUHighschool") %>% filter(term !="EDULess than Highschool") %>% filter(term !="EDUSome College") %>% filter(term != "quartileTHC2")%>%filter(term != "quartileTHC3")%>%filter(term != "quartileTHC4")%>% filter(term != "quartileTHCNA")%>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term != "Batch") %>% filter(term != "Mass") %>% filter(term !="EDUHighschool") %>% filter(term !="EDULess than Highschool") %>% filter(term !="EDUSome College") %>% filter(term != "quartileTHC2")%>%filter(term != "quartileTHC3")%>%filter(term != "quartileTHC4")%>% filter(term != "quartileTHCNA")%>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")

library(dotwhisker)
library(dplyr)
library(gdata)
all_models <- combine(All, White, Black)
colnames(all_models)[6] ="model"

plotRace_As<- dwplot(all_models, dodge_size = 1, dot_args = list(aes(colour = model),size = 5), whisker_args = list(aes(colour=model), size = 3)) %>% relabel_predictors(c("raceBlack" = "Black"))
plotRace_As
plotRace_As1 <- plotRace_As + theme_bw(base_size = 4) + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + ggtitle("CPT FA Fraction") + scale_color_brewer(palette="Set2") + theme(text = element_text(size = 16)) + labs(x = bquote('Coefficent Estimate with 95% CIs')) + theme(axis.text=element_text(size=17)) + theme(axis.title.x = element_text(margin = margin(t = 15))) + theme(panel.border = element_rect(fill=NA, colour = "black", size=1))  
plotRace_As1 
ggsave("metals.png", plotRace_essential1, bg='transparent')
png("AsRace.png", width = 9.9, height = 5.4, units = 'in', res = 300) 
plotRace_essential1 
dev.off()
                    