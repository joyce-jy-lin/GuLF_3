setwd("/Users/joycelin/Desktop/Gulf/Aim3")
library(sas7bdat)
library(tidyverse)
library(ggplot2)
library(broom)       # for tidy model output using tidy() function
library(Hmisc) 
library(gtsummary)
library(summarytools)
library(readxl)

# Table1 by metal conc and neuro scores -----------------------------------
library(tidyverse)
library(dplyr)
library(readxl)
data<- read_excel("CE_neuroscaled_quantiles.xlsx")

# make <60% detect metals binary detect non-detect, values are imputed LOD/SQRT2 for each metal
data$V[data$V < 0.00284] <-0
data$V[data$V > 0.00284] <-1
data$V <- factor(data$V)

data$Co[data$Co < 0.0000273] <-0
data$Co[data$Co > 0.0000273] <-1
data$Co <- factor(data$Co)

data$Cd[data$Cd < 0.00031] <-0
data$Cd[data$Cd > 0.00031] <-1
data$Cd <- factor(data$Cd)

data$Mo[data$Mo < 0.00046] <-0
data$Mo[data$Mo > 0.00046] <-1
data$Mo <- factor(data$Mo)

data$Sb[data$Sb < 0.00044] <-0
data$Sb[data$Sb > 0.00044] <-1
data$Sb <- factor(data$Sb)


# Age ---------------------------------------------------------------------
data <- data %>% mutate(CEagegroup = case_when(CE_AGE >= 70 ~ '5',
                                             CE_AGE >= 60 & CE_AGE <=69 ~ '4',
                                             CE_AGE >= 40 & CE_AGE<= 59 ~ '3',
                                             CE_AGE >= 20 & CE_AGE <= 39 ~ '2',
                                             CE_AGE < 20 ~ '1')) # end function
summaryCEage <- data %>%
  group_by(CEagegroup) %>%
  summarise(Count = n())
summaryCEage


summaryCEage<- data %>% group_by(CEagegroup) %>% summarise(median = median(As), n=n())

summaryCEage <- data %>% group_by(CEagegroup) %>% summarise(Se = quantile(Se, c(0.5, 0.2, 0.75)), q = c(0.5, 0.25, 0.75))
summaryCEage

## race subset
datablack <-filter(data, EN_RACE2 =="2")
datawhite <-filter(data, EN_RACE2 =="1")

summaryCEage <- datawhite %>%
  group_by(CEagegroup) %>%
  summarise(Count = n())
summaryCEage


# Race --------------------------------------------------------------------
describe(data$EN_RACE2)

summaryrace <- data %>% group_by(EN_RACE2) %>% summarise(Se = quantile(Se, c(0.5, 0.2, 0.75)), q = c(0.5, 0.25, 0.75))
(summaryrace)

## race subset smoking
summarysmok <- data %>%
  group_by(EN_FORMERSMOKER) %>%
  summarise(Count = n())
summarysmok

# Education ---------------------------------------------------------------
as.numeric(data$EN_EDU)
data%>% count(data$EN_EDU)

data <- data %>% mutate(EDUgroup = case_when(EN_EDU >= 18 ~ '4',
                                               EN_EDU >= 15 & EN_EDU <=17 ~ '3',
                                               EN_EDU >= 13 & EN_EDU<= 14 ~ '2',
                                               EN_EDU <=12 ~ '1')) # end function

summaryEDU <- data %>% group_by(EDUgroup) %>% summarise(Se = quantile(Se, c(0.5, 0.2, 0.75)), q = c(0.5, 0.25, 0.75))
summaryEDU

## race subset
datablack <-filter(data, EN_RACE2 =="2")
datawhite <-filter(data, EN_RACE2 =="1")

summaryedu <- datawhite %>%
  group_by(EDUgroup) %>%
  summarise(Count = n())
summaryedu

# Income ---------------------------------------------------------------
data%>% count(data$EN_TOTINCOME)
dat<-data%>%drop_na(EN_TOTINCOME)
as.numeric(dat$EN_TOTINCOME)
sum(is.na(dat$EN_TOTINCOME))


dat<- dat %>% mutate(incomegroup = case_when(EN_TOTINCOME >= 5 ~ '3',
                                                EN_TOTINCOME > 1 & EN_TOTINCOME<= 4 ~ '2',
                                                EN_TOTINCOME ==1 ~ '1')) # end function

summaryincome <- dat %>% group_by(incomegroup) %>% summarise(Se = quantile(Se, c(0.5, 0.2, 0.75)), q = c(0.5, 0.25, 0.75))
summaryincome

## race subset
datablack <-filter(dat, EN_RACE2 =="2")
datawhite <-filter(dat, EN_RACE2 =="1")

summaryinc <- datablack %>%
  group_by(incomegroup) %>%
  summarise(Count = n())
summaryinc


# State -------------------------------------------------------------------
data%>% count(data$EN_STATE)
summarystate <- data %>% group_by(EN_STATE) %>% summarise(Se = quantile(Se, c(0.5, 0.2, 0.75)), q = c(0.5, 0.25, 0.75))
summarystate

# Cum_THC -------------------------------------------------------------------
as.numeric(data$THC_CUMULATIVE1)
class(data$THC_CUMULATIVE1)
quantile(data$THC_CUMULATIVE1, na.rm= TRUE)
data <- data%>% mutate(quartilTHCCum = ntile(THC_CUMULATIVE1, 4)) 
describe(data$quartilTHCCum)

summaryTHC <- data %>% group_by(quartilTHCCum) %>% summarise(Se = quantile(Se, c(0.5, 0.2, 0.75)), q = c(0.5, 0.25, 0.75))
summaryTHC

datablack <-filter(data, EN_RACE2 =="2")
datawhite <-filter(data, EN_RACE2 =="1")

summaryTHC <- datawhite %>%
  group_by(quartilTHCCum) %>%
  summarise(Count = n())
summaryTHC


# Descriptive analysis of metals and neuro --------------------------------
library(ggplot2)
library(scales)
require(scales)

# tidy long data for facet wrap plot
data <- read_excel("GulfData_perfnerv.xlsx")
low_metals <- c("Cd", "Co", "Mo" ,"V", "Sb") # metals below 60% detect
`%notin%` <- Negate(`%in%`)

datlong<-pivot_longer(data, 6:23, names_to = "Metal", values_to = "Concentration")
datlongmain<-datlong %>% filter(Metal %notin% low_metals) 
datlonglow<- datlong %>% filter(Metal %in% low_metals)

# metal concentrations by race, try to produce facet wrap for all metals
datlongmain$EN_RACE2 <- factor(datlong$EN_RACE2, 
                        levels = c( 1, 2, 3, 4, 9),
                        labels = c("White", "Black", "Asian", "Other", "Multiracial"))

p2<- ggplot(datlongmain, aes(x=EN_RACE2, y=log10(Concentration), fill=EN_RACE2)) + geom_boxplot() + facet_wrap(~Metal, scale="free")
p2

# metal concentrations by state
datlongmain$EN_STATE <- factor(data$EN_STATE)
p2<- ggplot(datlongmain, aes(x=EN_STATE, y=log10(Concentration), fill=EN_STATE)) + geom_boxplot() + facet_wrap(~Metal, scale="free")
p2

#metal concentrations by THC Cumulative
datlongmain <- datlongmain %>% mutate(quartileTHC = ntile(THC_CUMULATIVE1, 4))
datlongmain$quartileTHC <- factor(datlongmain$quartileTHC)

p2<- ggplot(datlongmain, aes(x=quartileTHC, y=log10(Concentration), fill=quartileTHC)) + geom_boxplot() + facet_wrap(~Metal, scale="free")
p2

## GM and IQR
library(EnvStats)
data %>% dplyr:: summarize(n=n(), gm = geoMean(Se), mean = mean(Ni), SD = sd(Ni), GSD = geoSD(Se))

## race subset
datablack <-filter(data, EN_RACE2 =="2")
datawhite <-filter(data, EN_RACE2 =="1")

datawhite %>% dplyr:: summarize(n=n(), gm = geoMean(Al), mean = mean(Ni), SD = sd(Ni), GSD = geoSD(Al))


# summary stats for neuro tests
library(EnvStats)
data %>% dplyr:: summarize(n=n(), gm = geoMean(CE_BARS_CPT_COR_HIT_FRACTION), GSD = geoSD(CE_BARS_CPT_COR_HIT_FRACTION))
data  %>% dplyr:: group_by(EN_RACE2)%>% dplyr:: summarize(n=n(), gm = geoMean(CE_BARS_CPT_COR_HIT_FRACTION), GSD = geoSD(CE_BARS_CPT_COR_HIT_FRACTION))

## characteristics by race
data$CE_D1D_AVG_SYS_NUM 

describe(data$CE_C18_YN)
describe(datablack$CE_C19_YN)

## descriptive scatterplots of metals and neuro
setwd("/Users/joycelin/Desktop/Gulf/Aim3")
data <- read.csv("CE_perfneuroscaled_tertiles.csv")
data_mainlong<-pivot_longer(data, 6:23, names_to = "Metal", values_to = "Concentration")

data_mainlong <- data_mainlong %>% mutate(Concentration = scale(Concentration))
p1<- ggplot(data_mainlong, aes(x=(Concentration), y=CE_BARS_CPT_COR_HIT_FRACTION)) + geom_point(size=1) + facet_wrap(~Metal, scale="free") +  geom_smooth(method=loess)
p1

