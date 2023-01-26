setwd("/Users/joycelin/Desktop/Gulf/Aim3")
library(tidyverse)
library(readxl)
library(Hmisc) 
library(gtsummary)
library(summarytools)
library(ggplot2)

data <- read_excel("GulfData_perfnerv.xlsx")
low_metals <- c("Cd", "Co", "Mo" ,"V", "Sb") # metals below 60% detect
`%notin%` <- Negate(`%in%`)

datlong<-pivot_longer(data, 6:23, names_to = "Metal", values_to = "Concentration")
cols.num <- c(181:269)
datlong[cols.num] <- sapply(datlong[cols.num],as.numeric)
data[cols.num] <- sapply(data[cols.num],as.numeric)

# summary stats by race for neuro tests
library(EnvStats)
data %>% dplyr:: summarize(n=n(), mean = mean(CE_BARS_PRT_TOTAL_TAPS, na.rm=T), SD = sd(CE_BARS_PRT_TOTAL_TAPS, na.rm=T))
data  %>% dplyr:: group_by(EN_RACE2)%>% dplyr:: summarize(n=n(), mean = mean(CE_BARS_PRT_TOTAL_TAPS, na.rm=T), SD = sd(CE_BARS_PRT_TOTAL_TAPS, na.rm=T))
datlong$CE_V4B_TEST_A_TIME <- as.numeric(as.character(datlong$CE_V4B_TEST_A_TIME))
datlong$CE_V4B_TEST_B_TIME <- as.numeric(as.character(datlong$CE_V4B_TEST_B_TIME))
data$CE_V4B_TEST_A_TIME <- as.numeric(as.character(data$CE_V4B_TEST_A_TIME))
data$CE_V4B_TEST_B_TIME <- as.numeric(as.character(data$CE_V4B_TEST_B_TIME))

#negate neuro tests so for all tests lower score indicates worse cog function. results in inverse relationship between neuro scores and metals = adverse effect
datlong1<- datlong %>%  mutate(CE_BARS_CPT_FA_FRACTION = CE_BARS_CPT_FA_FRACTION *-1,
                               CE_BARS_CPT_COR_FA_FRACTION = CE_BARS_CPT_COR_FA_FRACTION *-1,
                               CE_BARS_CPT_FA_LATENCY = CE_BARS_CPT_FA_LATENCY *-1,
                               CE_BARS_CPT_HIT_LATENCY= CE_BARS_CPT_HIT_LATENCY *-1,
                               CE_BARS_CPT_MISSES= CE_BARS_CPT_MISSES *-1,
                               CE_BARS_MTS_AVE_COR_LAT = CE_BARS_MTS_AVE_COR_LAT*-1,
                               CE_BARS_SDT_AVE_COR_LAT = CE_BARS_SDT_AVE_COR_LAT*-1,
                               CE_BARS_SDT_TOTAL_ERRS = CE_BARS_SDT_TOTAL_ERRS*-1,
                               CE_BARS_SRT_AVE_COR_LAT = CE_BARS_SRT_AVE_COR_LAT*-1,
                               CE_BARS_SRT_TOTAL_ERRS = CE_BARS_SRT_TOTAL_ERRS*-1,
                               CE_V4B_TEST_A_TIME = CE_V4B_TEST_A_TIME*-1,
                               CE_V4B_TEST_B_TIME = CE_V4B_TEST_B_TIME*-1,)

dat1<- data %>%  mutate(CE_BARS_CPT_FA_FRACTION = CE_BARS_CPT_FA_FRACTION *-1,
                        CE_BARS_CPT_COR_FA_FRACTION = CE_BARS_CPT_COR_FA_FRACTION *-1,
                        CE_BARS_CPT_FA_LATENCY = CE_BARS_CPT_FA_LATENCY *-1,
                        CE_BARS_CPT_HIT_LATENCY= CE_BARS_CPT_HIT_LATENCY *-1,
                        CE_BARS_CPT_MISSES= CE_BARS_CPT_MISSES *-1,
                        CE_BARS_MTS_AVE_COR_LAT = CE_BARS_MTS_AVE_COR_LAT*-1,
                        CE_BARS_SDT_AVE_COR_LAT = CE_BARS_SDT_AVE_COR_LAT*-1,
                        CE_BARS_SDT_TOTAL_ERRS = CE_BARS_SDT_TOTAL_ERRS*-1,
                        CE_BARS_SRT_AVE_COR_LAT = CE_BARS_SRT_AVE_COR_LAT*-1,
                        CE_BARS_SRT_TOTAL_ERRS = CE_BARS_SRT_TOTAL_ERRS*-1,
                        CE_V4B_TEST_A_TIME = CE_V4B_TEST_A_TIME*-1,
                        CE_V4B_TEST_B_TIME = CE_V4B_TEST_B_TIME*-1)

#scale neuro outcomes so comparable
datlong2 <- datlong1 %>% mutate(s_CE_BARS_CPT_HITS = scale(CE_BARS_CPT_HITS),
                                s_CE_BARS_CPT_CR_FRACTION = scale(CE_BARS_CPT_CR_FRACTION),
                                s_CE_BARS_CPT_MISSES = scale(CE_BARS_CPT_MISSES),
                                s_CE_BARS_CPT_HIT_FRACTION = scale(CE_BARS_CPT_HIT_FRACTION),
                                s_CE_BARS_CPT_COR_HIT_FRACTION = scale(CE_BARS_CPT_COR_HIT_FRACTION),
                                s_CE_BARS_CPT_HIT_LATENCY  = scale(CE_BARS_CPT_HIT_LATENCY),
                                s_CE_BARS_CPT_FA_FRACTION = scale(CE_BARS_CPT_FA_FRACTION),
                                s_CE_BARS_CPT_COR_FA_FRACTION = scale(CE_BARS_CPT_COR_FA_FRACTION),
                                s_CE_BARS_CPT_FA_LATENCY = scale(CE_BARS_CPT_FA_LATENCY),
                                s_CE_BARS_CPT_D_PRIME = scale(CE_BARS_CPT_D_PRIME),
                                s_CE_BARS_DST_FORWARD_CNT = scale(CE_BARS_DST_FORWARD_CNT),
                                s_CE_BARS_DST_REVERSE_CNT = scale(CE_BARS_DST_REVERSE_CNT),
                                s_CE_BARS_MTS_AVE_COR_LAT = scale(CE_BARS_MTS_AVE_COR_LAT),
                                s_CE_BARS_MTS_COR_CNT = scale(CE_BARS_MTS_COR_CNT),
                                s_CE_BARS_SDT_AVE_COR_LAT = scale(CE_BARS_SDT_AVE_COR_LAT),
                                s_CE_BARS_SDT_TOTAL_ERRS = scale(CE_BARS_SDT_TOTAL_ERRS),
                                s_CE_BARS_SRT_AVE_COR_LAT = scale(CE_BARS_SRT_AVE_COR_LAT),
                                s_CE_BARS_SRT_TOTAL_ERRS = scale(CE_BARS_SRT_TOTAL_ERRS),
                                #s_CE_BARS_TAP_ALT_BOTH_AVG = scale(CE_BARS_TAP_ALT_BOTH_AVG),
                                #s_CE_BARS_TAP_RIGHT_PREF_AVG = scale(CE_BARS_TAP_RIGHT_PREF_AVG),
                                s_CE_BARS_PRT_TOTAL_TAPS = scale(CE_BARS_PRT_TOTAL_TAPS),
                                s_CE_V4B_TEST_A_TIME = scale(CE_V4B_TEST_A_TIME),
                                s_CE_V4B_TEST_B_TIME = scale(CE_V4B_TEST_B_TIME))

dat2 <- dat1 %>% mutate(s_CE_BARS_CPT_HITS = scale(CE_BARS_CPT_HITS),
                        s_CE_BARS_CPT_CR_FRACTION = scale(CE_BARS_CPT_CR_FRACTION),
                        s_CE_BARS_CPT_MISSES = scale(CE_BARS_CPT_MISSES),
                        s_CE_BARS_CPT_HIT_FRACTION = scale(CE_BARS_CPT_HIT_FRACTION),
                        s_CE_BARS_CPT_COR_HIT_FRACTION = scale(CE_BARS_CPT_COR_HIT_FRACTION),
                        s_CE_BARS_CPT_HIT_LATENCY  = scale(CE_BARS_CPT_HIT_LATENCY),
                        s_CE_BARS_CPT_FA_FRACTION = scale(CE_BARS_CPT_FA_FRACTION),
                        s_CE_BARS_CPT_COR_FA_FRACTION = scale(CE_BARS_CPT_COR_FA_FRACTION),
                        s_CE_BARS_CPT_FA_LATENCY = scale(CE_BARS_CPT_FA_LATENCY),
                        s_CE_BARS_CPT_D_PRIME = scale(CE_BARS_CPT_D_PRIME),
                        s_CE_BARS_DST_FORWARD_CNT = scale(CE_BARS_DST_FORWARD_CNT),
                        s_CE_BARS_DST_REVERSE_CNT = scale(CE_BARS_DST_REVERSE_CNT),
                        s_CE_BARS_MTS_AVE_COR_LAT = scale(CE_BARS_MTS_AVE_COR_LAT),
                        s_CE_BARS_MTS_COR_CNT = scale(CE_BARS_MTS_COR_CNT),
                        s_CE_BARS_SDT_AVE_COR_LAT = scale(CE_BARS_SDT_AVE_COR_LAT),
                        s_CE_BARS_SDT_TOTAL_ERRS = scale(CE_BARS_SDT_TOTAL_ERRS),
                        s_CE_BARS_SRT_AVE_COR_LAT = scale(CE_BARS_SRT_AVE_COR_LAT),
                        s_CE_BARS_SRT_TOTAL_ERRS = scale(CE_BARS_SRT_TOTAL_ERRS),
                        #s_CE_BARS_TAP_ALT_BOTH_AVG = scale(CE_BARS_TAP_ALT_BOTH_AVG),
                        #s_CE_BARS_TAP_RIGHT_PREF_AVG = scale(CE_BARS_TAP_RIGHT_PREF_AVG),
                        s_CE_BARS_PRT_TOTAL_TAPS = scale(CE_BARS_PRT_TOTAL_TAPS),
                        s_CE_V4B_TEST_A_TIME = scale(CE_V4B_TEST_A_TIME),
                        s_CE_V4B_TEST_B_TIME = scale(CE_V4B_TEST_B_TIME))

# quartiles ---------------------------------------------------------------
# #Split data into quartiles by metal log transformed
dat2 <- dat2 %>% mutate(AlQuartile = ntile(Al, 4))
# make into factor variable
dat2$AlQuartile <- factor(dat2$AlQuartile)

dat2 <- dat2 %>% mutate(AsQuartile = ntile(As, 4))
dat2$AsQuartile <- factor(dat2$AsQuartile)

dat2 <- dat2 %>% mutate(CaQuartile = ntile(Ca, 4))
dat2$CaQuartile <- factor(dat2$CaQuartile)

dat2 <- dat2 %>% mutate(CrQuartile = ntile(Cr, 4))
dat2$CrQuartile <- factor(dat2$CrQuartile)

dat2 <- dat2 %>% mutate(CuQuartile = ntile(Cu, 4))
dat2$CuQuartile <- factor(dat2$CuQuartile)

dat2 <- dat2 %>% mutate(FeQuartile = ntile(Fe, 4))
dat2$FeQuartile <- factor(dat2$FeQuartile)

dat2 <- dat2 %>% mutate(PbQuartile = ntile(Pb, 4))
dat2$PbQuartile <- factor(dat2$PbQuartile)

dat2 <- dat2 %>% mutate(MgQuartile = ntile(Mg, 4))
dat2$MgQuartile <- factor(dat2$MgQuartile)

dat2 <- dat2 %>% mutate(MnQuartile = ntile(Mn, 4))
dat2$MnQuartile <- factor(dat2$MnQuartile)

dat2 <- dat2 %>% mutate(HgQuartile = ntile(Hg, 4))
dat2$HgQuartile <- factor(dat2$HgQuartile)

dat2 <- dat2 %>% mutate(NiQuartile = ntile(Ni, 4))
dat2$NiQuartile <- factor(dat2$NiQuartile)

dat2 <- dat2 %>% mutate(SeQuartile = ntile(Se, 4))
dat2$SeQuartile <- factor(dat2$SeQuartile)

dat2 <- dat2 %>% mutate(ZnQuartile = ntile(Zn, 4))
dat2$ZnQuartile <- factor(dat2$ZnQuartile)

## group for metals below 60% detect ----------------------
dat2 <- dat2 %>% mutate(CdBinary = Cd)
dat2$CdBinary[dat2$CdBinary < 0.00031] <-0
dat2$CdBinary[dat2$CdBinary > 0.00031] <-1
dat2$CdBinary <- factor(dat2$CdBinary)

dat2 <- dat2 %>% mutate(CoBinary = Co)
dat2$CoBinary[dat2$CoBinary < 0.0000273] <-0
dat2$CoBinary[dat2$CoBinary > 0.0000273] <-1
dat2$CoBinary <- factor(dat2$CoBinary)

dat2 <- dat2 %>% mutate(MoBinary = Mo)
dat2$MoBinary[dat2$MoBinary <  0.00046] <-0
dat2$MoBinary[dat2$MoBinary >  0.00046] <-1
dat2$MoBinary <- factor(dat2$MoBinary)

Sb <- sort(dat2$Sb)
Sb
nth(Sb, 286)
dat2 <- dat2 %>% mutate(SbTertile = case_when(Sb <= 0.00044 ~ '0',
                                                      Sb > 0.00044 & Sb <= 0.0343302 ~ '1',
                                                      Sb >0.0343302 ~ '2'))
dat2$SbTertile <- factor(dat2$SbTertile)

V <- sort(dat2$V)
V
nth(V, 289)
dat2 <- dat2 %>% mutate(VTertile = case_when(V <= 0.00284 ~ '0',
                                                     V > 0.00284 & V <= 0.01947112 ~ '1',
                                                     V >0.01947112 ~ '2'))
dat2$VTertile <- factor(dat2$VTertile)

dat2 <- dat2 %>% mutate(quartileTHC = ntile(THC_CUMULATIVE1, 4))
dat2$quartileTHC <- factor(dat2$quartileTHC)

## for fully adjusted regressions ------------------------------------------------
dat2 <- dat2 %>% mutate(EDU = case_when(EN_EDU >= 18 ~ 'College or more',
                                                EN_EDU >= 15 & EN_EDU <=17 ~ 'Some College',
                                                EN_EDU >= 13 & EN_EDU<= 14 ~ 'Highschool',
                                                EN_EDU <=12 ~ 'Less than Highschool')) # end function

dat2$EDU <- relevel(factor(dat2$EDU), ref = "College or more")

dat2 <- dat2 %>% mutate(race = case_when(EN_RACE2 >=3 ~ 'Other', 
                                                 EN_RACE2 == 2 ~ 'Black',
                                                 EN_RACE2 == 1 ~ 'White')) # end function

dat2$race <- relevel(factor(dat2$race), ref = "White")

dat2<- dat2 %>% mutate(headinjuryever = case_when(headinjury == 2 ~ '1',
                                                  headinjury ==0 ~ '0',
                                                  headinjury ==1 ~ '0')) # head injury where loss consciousness
dat2<- dat2 %>% mutate(drinklot = case_when(drinksweek <= 14 ~ '0',
                                            drinksweek >14 ~ '1')) # more than 14 drinks/week

dat2$EN_FORMERSMOKER <- as.factor(dat2$EN_FORMERSMOKER)
dat2$headinjury <- factor(dat2$headinjury)
dat2$headinjuryever <- factor(dat2$headinjuryever)
dat2$drinklot <- factor(dat2$drinklot)


write.csv(dat2, file = "CE_neuroscaled_quantiles.csv")



