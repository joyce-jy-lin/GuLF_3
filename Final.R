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

## QUARTILE regression table-----------------------------------
#CPT D Prime
metal <- "AsQ"
predictors <- c( "CE_AGE", "CE_BMI", "CE_C1", "EN_FORMERSMOKER","marital")
outcome <- "CE_BARS_CPT_D_PRIME"
fit_models <- function(metal, predictors, outcome){
  tmp <- data_log %>% dplyr::select(outcome, all_of(predictors), metal)
  mod1 <- lm(CE_BARS_CPT_D_PRIME ~ ., data = tmp)
  tmp <- whitesubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod2 <- lm(CE_BARS_CPT_D_PRIME ~ ., data = tmp)
  tmp <- blacksubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod3 <- lm(CE_BARS_CPT_D_PRIME ~ ., data = tmp)
  
  summary1 <- summary(mod1)$coefficients[7,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "full"
  )
  summary2 <- summary(mod2)$coefficients[7,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "white"
  )
  summary3 <- summary(mod3)$coefficients[7,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "black"
  )
  
  bind_rows(summary1, summary2, summary3) %>% mutate(metal = metal)
  
}
fit_models(metal = "AsQ", predictors = predictors, outcome=outcome)

list_metals <- c("AsQ", "CdBinary", "CrQ", "CuQ", "FeQ", "HgQ", "MgQ", "MnQ", "PbQ","SeQ","ZnQ")
Allmetal <- purrr::map_dfr(.x = list_metals, .f = fit_models, predictors=predictors,
                           outcome=outcome)
print(Allmetal,digits=2)
Allmetal<- Allmetal %>% 
  mutate_if(is.numeric, round, digits = 4)

# DST forward
metal <- "AsQ"
predictors <- c( "CE_AGE", "CE_BMI", "CE_C1", "EN_FORMERSMOKER","drinklot","marital")
outcome <- "CE_BARS_DST_FORWARD_CNT"
fit_models <- function(metal, predictors, outcome){
  tmp <- data_log %>% dplyr::select(outcome, all_of(predictors), metal)
  mod1 <- lm(CE_BARS_DST_FORWARD_CNT ~ ., data = tmp)
  tmp <- whitesubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod2 <- lm(CE_BARS_DST_FORWARD_CNT ~ ., data = tmp)
  tmp <- blacksubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod3 <- lm(CE_BARS_DST_FORWARD_CNT ~ ., data = tmp)
  
  summary1 <- summary(mod1)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "full"
  )
  summary2 <- summary(mod2)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "white"
  )
  summary3 <- summary(mod3)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "black"
  )
  
  bind_rows(summary1, summary2, summary3) %>% mutate(metal = metal)
  
}
fit_models(metal = "AsQ", predictors = predictors, outcome=outcome)

list_metals <- c("AsQ", "CdBinary", "CrQ", "CuQ", "FeQ", "HgQ", "MgQ", "MnQ", "PbQ","SeQ","ZnQ")
Allmetal <- purrr::map_dfr(.x = list_metals, .f = fit_models, predictors=predictors,
                           outcome=outcome)
print(Allmetal,digits=2)
Allmetal<- Allmetal %>% 
  mutate_if(is.numeric, round, digits = 4)

# DST reverse
metal <- "AsQ"
predictors <- c( "CE_AGE", "CE_BMI", "CE_C1", "EN_FORMERSMOKER","drinklot","marital")
outcome <- "CE_BARS_DST_REVERSE_CNT"
fit_models <- function(metal, predictors, outcome){
  tmp <- data_log %>% dplyr::select(outcome, all_of(predictors), metal)
  mod1 <- lm(CE_BARS_DST_REVERSE_CNT ~ ., data = tmp)
  tmp <- whitesubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod2 <- lm(CE_BARS_DST_REVERSE_CNT ~ ., data = tmp)
  tmp <- blacksubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod3 <- lm(CE_BARS_DST_REVERSE_CNT ~ ., data = tmp)
  
  summary1 <- summary(mod1)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "full"
  )
  summary2 <- summary(mod2)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "white"
  )
  summary3 <- summary(mod3)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "black"
  )
  
  bind_rows(summary1, summary2, summary3) %>% mutate(metal = metal)
  
}
fit_models(metal = "AsQ", predictors = predictors, outcome=outcome)

list_metals <- c("AsQ", "CdBinary", "CrQ", "CuQ", "FeQ", "HgQ", "MgQ", "MnQ", "PbQ","SeQ","ZnQ")
Allmetal <- purrr::map_dfr(.x = list_metals, .f = fit_models, predictors=predictors,
                           outcome=outcome)
print(Allmetal,digits=2)
Allmetal<- Allmetal %>% 
  mutate_if(is.numeric, round, digits = 4)


# CPT cr
metal <- "AsQ"
predictors <- c( "CE_AGE", "CE_BMI", "CE_C1", "EN_FORMERSMOKER","drinklot","marital")
outcome <- "CPTcr1"
fit_models <- function(metal, predictors, outcome){
  tmp <- data_log %>% dplyr::select(outcome, all_of(predictors), metal)
  mod1 <- glm(CPTcr1 ~ ., data = tmp, family = binomial)
  tmp <- blacksubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod3 <- glm(CPTcr1 ~ ., data = tmp, family = binomial)
  
  summary1 <- summary(mod1)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "full"
  )
  summary3 <- summary(mod3)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "black"
  )
  
  bind_rows(summary1, summary3) %>% mutate(metal = metal)
  
}
fit_models(metal = "AsQ", predictors = predictors, outcome=outcome)


list_metals <- c("AsQ", "CdBinary", "CrQ", "CuQ", "FeQ", "HgQ", "MgQ", "MnQ", "PbQ","SeQ","ZnQ")
Allmetal <- purrr::map_dfr(.x = list_metals, .f = fit_models, predictors=predictors,
                           outcome=outcome)
print(Allmetal,digits=2)
Allmetal<- Allmetal %>% 
  mutate_if(is.numeric, round, digits = 4)

# CPT hit
metal <- "AsQ"
predictors <- c( "CE_AGE", "CE_BMI", "CE_C1", "EN_FORMERSMOKER","drinklot","marital")
outcome <- "CPThit1"
fit_models <- function(metal, predictors, outcome){
  tmp <- data_log %>% dplyr::select(outcome, all_of(predictors), metal)
  mod1 <- glm(CPThit1 ~ ., data = tmp, family = binomial)
  tmp <- whitesubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod2 <- glm(CPThit1 ~ ., data = tmp, family = binomial)
  tmp <- blacksubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod3 <- glm(CPThit1 ~ ., data = tmp, family = binomial)
  
  summary1 <- summary(mod1)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "full"
  )
  summary2 <- summary(mod2)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "white"
  )
  summary3 <- summary(mod3)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "black"
  )
  
  bind_rows(summary1, summary2, summary3) %>% mutate(metal = metal)
  
}
fit_models(metal = "AsQ", predictors = predictors, outcome=outcome)


list_metals <- c("AsQ", "CdBinary", "CrQ", "CuQ", "FeQ", "HgQ", "MgQ", "MnQ", "PbQ","SeQ","ZnQ")
Allmetal <- purrr::map_dfr(.x = list_metals, .f = fit_models, predictors=predictors,
                           outcome=outcome)
print(Allmetal,digits=2)
Allmetal<- Allmetal %>% 
  mutate_if(is.numeric, round, digits = 4)

## Continuous -----------------------------------------------------------------------------

## QUARTILE regression table-----------------------------------
#CPT D Prime
metal <- "As"
predictors <- c( "CE_AGE", "CE_BMI", "CE_C1", "EN_FORMERSMOKER","drinklot","marital")
outcome <- "CE_BARS_CPT_D_PRIME"
fit_models <- function(metal, predictors, outcome){
  tmp <- data_log %>% dplyr::select(outcome, all_of(predictors), metal)
  mod1 <- lm(CE_BARS_CPT_D_PRIME ~ ., data = tmp)
  tmp <- whitesubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod2 <- lm(CE_BARS_CPT_D_PRIME ~ ., data = tmp)
  tmp <- blacksubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod3 <- lm(CE_BARS_CPT_D_PRIME ~ ., data = tmp)
  
  summary1 <- summary(mod1)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "full"
  )
  summary2 <- summary(mod2)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "white"
  )
  summary3 <- summary(mod3)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "black"
  )
  
  bind_rows(summary1, summary2, summary3) %>% mutate(metal = metal)
  
}
fit_models(metal = "As", predictors = predictors, outcome=outcome)

list_metals <- c("Al", "As", "CdBinary", "Cr", "Cu", "Fe", "Hg", "Mg", "Mn", "Pb","Se","Zn")
Allmetal <- purrr::map_dfr(.x = list_metals, .f = fit_models, predictors=predictors,
                           outcome=outcome)
print(Allmetal,digits=2)
Allmetal<- Allmetal %>% 
  mutate_if(is.numeric, round, digits = 4)

# DST forward
metal <- "As"
predictors <- c( "CE_AGE", "CE_BMI", "CE_C1", "EN_FORMERSMOKER","drinklot","marital")
outcome <- "CE_BARS_DST_FORWARD_CNT"
fit_models <- function(metal, predictors, outcome){
  tmp <- data_log %>% dplyr::select(outcome, all_of(predictors), metal)
  mod1 <- lm(CE_BARS_DST_FORWARD_CNT ~ ., data = tmp)
  tmp <- whitesubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod2 <- lm(CE_BARS_DST_FORWARD_CNT ~ ., data = tmp)
  tmp <- blacksubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod3 <- lm(CE_BARS_DST_FORWARD_CNT ~ ., data = tmp)
  
  summary1 <- summary(mod1)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "full"
  )
  summary2 <- summary(mod2)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "white"
  )
  summary3 <- summary(mod3)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "black"
  )
  
  bind_rows(summary1, summary2, summary3) %>% mutate(metal = metal)
  
}
fit_models(metal = "As", predictors = predictors, outcome=outcome)

list_metals <- c("Al","As", "CdBinary", "Cr", "Cu", "Fe", "Hg", "Mg", "Mn", "Pb","Se","Zn")
Allmetal <- purrr::map_dfr(.x = list_metals, .f = fit_models, predictors=predictors,
                           outcome=outcome)
print(Allmetal,digits=2)
Allmetal<- Allmetal %>% 
  mutate_if(is.numeric, round, digits = 4)

# DST reverse
metal <- "As"
predictors <- c( "CE_AGE", "CE_BMI", "CE_C1", "EN_FORMERSMOKER","drinklot","marital")
outcome <- "CE_BARS_DST_REVERSE_CNT"
fit_models <- function(metal, predictors, outcome){
  tmp <- data_log %>% dplyr::select(outcome, all_of(predictors), metal)
  mod1 <- lm(CE_BARS_DST_REVERSE_CNT ~ ., data = tmp)
  tmp <- whitesubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod2 <- lm(CE_BARS_DST_REVERSE_CNT ~ ., data = tmp)
  tmp <- blacksubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod3 <- lm(CE_BARS_DST_REVERSE_CNT ~ ., data = tmp)
  
  summary1 <- summary(mod1)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "full"
  )
  summary2 <- summary(mod2)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "white"
  )
  summary3 <- summary(mod3)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "black"
  )
  
  bind_rows(summary1, summary2, summary3) %>% mutate(metal = metal)
  
}
fit_models(metal = "As", predictors = predictors, outcome=outcome)

list_metals <- c("Al","As", "CdBinary", "Cr", "Cu", "Fe", "Hg", "Mg", "Mn", "Pb","Se","Zn")
Allmetal <- purrr::map_dfr(.x = list_metals, .f = fit_models, predictors=predictors,
                           outcome=outcome)
print(Allmetal,digits=2)
Allmetal<- Allmetal %>% 
  mutate_if(is.numeric, round, digits = 4)


# CPT cr
metal <- "As"
predictors <- c( "CE_AGE", "CE_BMI", "CE_C1", "EN_FORMERSMOKER","drinklot","marital")
outcome <- "CPTcr1"
fit_models <- function(metal, predictors, outcome){
  tmp <- data_log %>% dplyr::select(outcome, all_of(predictors), metal)
  mod1 <- glm(CPTcr1 ~ ., data = tmp, family = binomial)
  tmp <- blacksubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod3 <- glm(CPTcr1 ~ ., data = tmp, family = binomial)
  
  summary1 <- summary(mod1)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "full"
  )
  summary3 <- summary(mod3)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "black"
  )
  
  bind_rows(summary1, summary3) %>% mutate(metal = metal)
  
}
fit_models(metal = "As", predictors = predictors, outcome=outcome)


list_metals <- c("Al", "As", "CdBinary", "Cr", "Cu", "Fe", "Hg", "Mg", "Mn", "Pb","Se","Zn")
Allmetal <- purrr::map_dfr(.x = list_metals, .f = fit_models, predictors=predictors,
                           outcome=outcome)
print(Allmetal,digits=2)
Allmetal<- Allmetal %>% 
  mutate_if(is.numeric, round, digits = 4)

# CPT hit
metal <- "As"
predictors <- c( "CE_AGE", "CE_BMI", "CE_C1", "EN_FORMERSMOKER","drinklot","marital")
outcome <- "CPThit1"
fit_models <- function(metal, predictors, outcome){
  tmp <- data_log %>% dplyr::select(outcome, all_of(predictors), metal)
  mod1 <- glm(CPThit1 ~ ., data = tmp, family = binomial)
  tmp <- whitesubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod2 <- glm(CPThit1 ~ ., data = tmp, family = binomial)
  tmp <- blacksubgroup %>% dplyr::select(outcome, all_of(predictors), metal)
  mod3 <- glm(CPThit1 ~ ., data = tmp, family = binomial)
  
  summary1 <- summary(mod1)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "full"
  )
  summary2 <- summary(mod2)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "white"
  )
  summary3 <- summary(mod3)$coefficients[8,] %>% t() %>% data.frame() %>% mutate(
    lb = Estimate - (1.96*Std..Error),
    ub = Estimate + (1.96 * Std..Error),
    model = "black"
  )
  
  bind_rows(summary1, summary2, summary3) %>% mutate(metal = metal)
  
}
fit_models(metal = "As", predictors = predictors, outcome=outcome)


list_metals <- c("Al", "As", "CdBinary", "Cr", "Cu", "Fe", "Hg", "Mg", "Mn", "Pb","Se","Zn")
Allmetal <- purrr::map_dfr(.x = list_metals, .f = fit_models, predictors=predictors,
                           outcome=outcome)
print(Allmetal,digits=2)
Allmetal<- Allmetal %>% 
  mutate_if(is.numeric, round, digits = 4)


## Linear trend assessment -----------------------------------------
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

blacksubgroup<-filter(data_log, Race =="Black")
whitesubgroup<-filter(data_log, Race =="White")

mod1 <- lm(CE_BARS_CPT_D_PRIME ~ ZnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2 <- lm(CE_BARS_CPT_D_PRIME ~ ZnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3 <- lm(CE_BARS_CPT_D_PRIME ~ ZnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)
summary(mod1)
summary(mod2)
summary(mod3)

mod1 <- lm(CE_BARS_DST_FORWARD_CNT ~ ZnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2 <- lm(CE_BARS_DST_FORWARD_CNT ~ ZnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3 <- lm(CE_BARS_DST_FORWARD_CNT ~ ZnQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)
summary(mod1)
summary(mod2)
summary(mod3)

mod1 <- lm(CE_BARS_DST_REVERSE_CNT ~ HgQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2 <- lm(CE_BARS_DST_REVERSE_CNT ~ HgQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3 <- lm(CE_BARS_DST_REVERSE_CNT ~ HgQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)
summary(mod1)
summary(mod2)
summary(mod3)

mod1 <- glm(CPTcr1 ~ PbQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family=binomial)
mod2 <- glm(CPTcr1 ~ PbQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family=binomial)
mod3 <- glm(CPTcr1 ~ PbQ + CE_AGE + CE_BMI + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family=binomial)
vif(mod1)
summary(mod1)
summary(mod2)
summary(mod3)

## Dot whisker plots --------------------------------------------------------------------
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

blacksubgroup<-filter(data_log, Race =="Black")
whitesubgroup<-filter(data_log, Race =="White")

# CPT D Prime
mod1<- lm(CE_BARS_CPT_D_PRIME ~ HgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ HgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ HgQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- vctrs::vec_c(Black, White, All)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  scale_x_continuous(limits=c(-1, 0.7),breaks=c(-1, -0.5, 0, 0.5)) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75))
DSTr_Mn

# DST Forward
mod1<- lm(CE_BARS_DST_FORWARD_CNT ~ MnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_DST_FORWARD_CNT ~ MnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_FORWARD_CNT ~ MnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  scale_x_continuous(limits=c(-1.75, 1.6),breaks=c( -1, 0, 1)) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75))
DSTr_Mn

# DST Reverse
mod1<- lm(CE_BARS_DST_REVERSE_CNT ~ ZnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_DST_REVERSE_CNT ~ ZnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_REVERSE_CNT ~ ZnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  scale_x_continuous(limits=c(-1.1, 1.6),breaks=c(-1, 0, 1)) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75))
DSTr_Mn

# CPT Hit
mod1<- glm(CPThit1 ~ MnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPThit1 ~ MnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPThit1 ~ MnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75))+
  scale_x_continuous(limits = c(-4.6, 3.2), breaks = log(pretty(exp(Black$estimate),n=5)), labels = (pretty(exp(Black$estimate),n=5)))

DSTr_Mn

# CPT CR
mod1<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPTcr1 ~ PbQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")

all_models <- combine(Black, All)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75)) +
  scale_x_continuous(breaks = log(pretty(exp(Black$estimate),n=2)), labels = (pretty(exp(Black$estimate),n=2)))
DSTr_Mn


## Plots adjusted for other metals continuously -------------------------------
# CPT D Prime
mod1<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + HgQ + AsQ + CuQ + SeQ + ZnQ + PbQ + MnQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + HgQ + AsQ + CuQ + SeQ + ZnQ + PbQ + MnQ  + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + HgQ + AsQ + CuQ + SeQ + ZnQ + PbQ + MnQ  + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

library(car)
vif(mod2)

library(dotwhisker)
library(dplyr)
library(gdata)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "MnQ")%>% filter(term != "AsQ")%>% filter(term != "CuQ")%>% filter(term != "MgQ")%>% filter(term != "PbQ")%>% filter(term != "FeQ")%>% filter(term != "CrQ")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "MnQ")%>% filter(term != "AsQ")%>% filter(term != "CuQ")%>% filter(term != "MgQ")%>% filter(term != "PbQ")%>% filter(term != "FeQ")%>% filter(term != "CrQ")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "MnQ")%>% filter(term != "AsQ")%>% filter(term != "CuQ")%>% filter(term != "MgQ")%>% filter(term != "PbQ")%>% filter(term != "FeQ")%>% filter(term != "CrQ")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75))
DSTr_Mn

# DST Forward
mod1<- lm(CE_BARS_DST_FORWARD_CNT ~ SeQ + ZnQ + HgQ + CrQ + PbQ + CuQ + MnQ + AsQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_DST_FORWARD_CNT ~ SeQ + ZnQ + HgQ + CrQ + PbQ + CuQ + MnQ + AsQ  + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_FORWARD_CNT ~ SeQ + ZnQ + HgQ + CrQ + PbQ + CuQ + MnQ + AsQ  + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

vif(mod3)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "As")%>% filter(term != "Cu")%>% filter(term != "Mg")%>% filter(term != "Pb")%>% filter(term != "Fe")%>% filter(term != "Cr")%>% filter(term != "Zn")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "As")%>% filter(term != "Cu")%>% filter(term != "Mg")%>% filter(term != "Pb")%>% filter(term != "Fe")%>% filter(term != "Cr")%>% filter(term != "Zn")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Mn")%>% filter(term != "As")%>% filter(term != "Cu")%>% filter(term != "Mg")%>% filter(term != "Pb")%>% filter(term != "Fe")%>% filter(term != "Cr")%>% filter(term != "Zn")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  scale_x_continuous(limits=c(-2, 1.6),breaks=c( -1, 0, 1)) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75))
DSTr_Mn

# DST Reverse
mod1<- lm(CE_BARS_DST_REVERSE_CNT ~ MnQ + AsQ+ CuQ + ZnQ + HgQ + CrQ + PbQ + SeQ  + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_DST_REVERSE_CNT ~ MnQ + AsQ+ CuQ + ZnQ + HgQ + CrQ + PbQ + SeQ  + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_REVERSE_CNT ~ MnQ + AsQ+ CuQ + ZnQ + HgQ + CrQ + PbQ + SeQ  + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

vif(mod3)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Fe")%>% filter(term != "Mg")%>% filter(term != "Mn")%>% filter(term != "Zn")%>% filter(term != "As")%>% filter(term != "Cu")%>% filter(term != "Cr")%>% filter(term != "Pb")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Fe")%>% filter(term != "Mg")%>% filter(term != "Mn")%>% filter(term != "Zn")%>% filter(term != "As")%>% filter(term != "Cu")%>% filter(term != "Cr")%>% filter(term != "Pb")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "Fe")%>% filter(term != "Mg")%>% filter(term != "Mn")%>% filter(term != "Zn")%>% filter(term != "As")%>% filter(term != "Cu")%>% filter(term != "Cr")%>% filter(term != "Pb")

all_models <- combine(Black, White, All)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75))
DSTr_Mn

# CPT Hit # colinearity issue
mod1<- glm(CPThit1 ~ SeQ + ZnQ + HgQ + CrQ + PbQ + CuQ + MnQ + AsQ  + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPThit1 ~ SeQ + ZnQ + HgQ + CrQ + PbQ + CuQ + MnQ + AsQ  + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPThit1 ~ SeQ + ZnQ + HgQ + CrQ + PbQ + CuQ + MnQ + AsQ  + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

vif(mod3)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Mg")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Fe")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Mg")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Fe")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Mg")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Fe")

all_models <- combine(Black, All)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75))+
  scale_x_continuous(limits = c(-4.6, 4),breaks=c(-4,-2, 0, 2, 4))

DSTr_Mn

# CPT CR
mod1<- glm(CPTcr1 ~ SeQ + ZnQ + HgQ + CrQ + PbQ + CuQ + MnQ + AsQ  +CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPTcr1 ~ SeQ + ZnQ + HgQ + CrQ + PbQ + CuQ + MnQ + AsQ  +CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPTcr1 ~ SeQ + ZnQ + HgQ + CrQ + PbQ + CuQ + MnQ + AsQ  +CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

All <-broom::tidy(mod1) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Mg")%>% filter(term != "Fe")%>% filter(term != "Zn")
White <-broom::tidy(mod2) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Mg")%>% filter(term != "Fe")%>% filter(term != "Zn")
Black <-broom::tidy(mod3) %>% filter(term != "CE_AGE") %>% filter(term !="CE_C1") %>% filter(term != "EN_FORMERSMOKER")%>%filter(term != "EN_FORMERSMOKER1")%>%filter(term != "passivesmoke")%>% filter(term != "passivesmoke1")%>% filter(term != "passivesmokeNA")%>%filter(term != "CE_BMI")%>% filter(term != "headinjuryever1")%>%filter(term != "drinklot1")%>% filter(term != "maritalSeparated or single")%>% filter(term != "Caffiene2h1")%>% filter(term != "Caffiene2h")%>% filter(term != "As")%>% filter(term != "Mn")%>% filter(term != "Pb")%>% filter(term != "Cr")%>% filter(term != "Cu")%>% filter(term != "Mg")%>% filter(term != "Fe")%>% filter(term != "Zn")

vif(mod1)
all_models <- combine(All)
colnames(all_models)[6] ="model"

DSTr_Mn<- dwplot(all_models, dodge_size = 0.5, dot_args = list(aes(colour = model),size = 1.1), whisker_args = list(aes(colour=model), size = 0.7)) %>% 
  relabel_predictors(c( "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4"))
DSTr_Mn
DSTr_Mn <- DSTr_Mn + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw() +
  theme(text = element_text(size = 11)) + theme(axis.text=element_text(size=11)) + 
  theme(axis.title.x = element_text(margin = margin(t = 11))) + 
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75)) +
  scale_x_continuous(limits=c(-2,11),breaks=c( 0, 2,4, 6, 8, 10))
DSTr_Mn

## Adjusted regression table ---------------------------------
## Linear trend assessment
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

blacksubgroup<-filter(data_log, Race =="Black")
whitesubgroup<-filter(data_log, Race =="White")
# CPT D Prime
mod1<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_CPT_D_PRIME ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

summary(mod1)
confint(mod1)
summary(mod2)
confint(mod2)
summary(mod3)
confint(mod3)

# CPT Hit
mod1<- glm(CPThit1 ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ+ CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPThit1 ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ+ CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPThit1 ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ+ CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

summary(mod1)
confint(mod1)
summary(mod2)
confint(mod2)
summary(mod3)
confint(mod3)

# CPT CR
mod1<- glm(CPTcr1 ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log, family = binomial)
mod2<- glm(CPTcr1 ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ+ CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup, family = binomial)
mod3<- glm(CPTcr1 ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ+ CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup, family = binomial)

summary(mod1)
confint(mod1)
summary(mod2)
summary(mod3)
confint(mod3)

# DST Forward
mod1<- lm(CE_BARS_DST_FORWARD_CNT ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_DST_FORWARD_CNT ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_FORWARD_CNT ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

summary(mod1)
confint(mod1)
summary(mod2)
confint(mod2)
summary(mod3)
confint(mod3)

# DST Reverse
mod1<- lm(CE_BARS_DST_REVERSE_CNT ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = data_log)
mod2<- lm(CE_BARS_DST_REVERSE_CNT ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = whitesubgroup)
mod3<- lm(CE_BARS_DST_REVERSE_CNT ~ CrQ + AsQ + MnQ + CuQ + ZnQ + PbQ + HgQ + SeQ + CE_BMI + CE_AGE + CE_C1 + EN_FORMERSMOKER + drinklot + marital, data = blacksubgroup)

summary(mod1)
confint(mod1)
summary(mod2)
summary(mod3)
confint(mod3)

## adjusted for other metals continuously -------------------------------
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


## Plotting zoom in from BKMR output

CuCr<- pred.resp.bivar.levels.data %>%
  filter(variable1 == "Cu", variable2 == "Cr")
MnCr<- pred.resp.bivar.levels.data %>%
  filter(variable1 == "Mn", variable2 == "Cr")

p<- ggplot(MnCr, aes(z1, est))+ geom_smooth(aes(col = quantile), stat = "identity")+
  geom_hline(yintercept=0, linetype= "dashed") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Outcome") +
  xlab("Exposure")
p


