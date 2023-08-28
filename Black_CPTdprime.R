
if (!require("pacman")) install.packages("pacman", INSTALL_opts = '--no-lock')
pacman::p_load(bkmr, ggplot2, bkmrhat, rstan, coda)
Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet") # for future package
future::plan(strategy = future::multisession, workers=10, .skip=TRUE)

library(tidyverse)
library(readxl)
library(ggplot2)

setwd("/users/jyijoyce")
data <- read.csv("NEW_CE_neuro_quantiles.csv")
data <- data %>% mutate(across(c(Mg, Al, Ca, Cr, Mn, Fe, Ni, Cu, Zn, As, Se, Hg, Pb), log10))

data <- data %>% mutate(Al = scale(Al),
                        As = scale(As),
                        Ca = scale(Ca),
                        Cr = scale(Cr),
                        Cu = scale(Cu),
                        Fe = scale(Fe),
                        Pb = scale(Pb),
                        Mg = scale(Mg),
                        Mn = scale(Mn),
                        Hg = scale(Hg),
                        Ni = scale(Ni),
                        Se = scale(Se),
                        Zn = scale(Zn))

data <- data %>% mutate(CE_BARS_DST_REVERSE_CNT = scale(CE_BARS_DST_REVERSE_CNT))
data <- data %>% mutate(CE_BARS_CPT_D_PRIME = scale(CE_BARS_CPT_D_PRIME))
data <- data %>% mutate(CE_AGE = scale(CE_AGE))
data <- data %>% mutate(CE_BMI = scale(CE_BMI))

data<- data %>% mutate(headinjuryever = case_when(headinjury == 2 ~ '1',
                                                  headinjury == 0 ~ '0',
                                                  headinjury == 1 ~ '0')) # head injury where loss consciousness
data<- data %>% mutate(drinklot = case_when(drinksweek <= 14 ~ '0',
                                            drinksweek >14 ~ '1')) # more than 14 drinks/week
data<- data %>% mutate(marital = case_when(EN_MARITAL == 1 ~ 'Married or with partner',
                                           EN_MARITAL == 6 ~ 'Married or with partner',
                                           EN_MARITAL == 2 ~ 'Separated or single',
                                           EN_MARITAL == 4 ~ 'Separated or single',
                                           EN_MARITAL == 3 ~ 'Separated or single',
                                           EN_MARITAL == 5 ~ 'Separated or single')) 
data<- data %>% mutate(BMI = case_when(CE_BMI <= 24.9 ~ 'healthy',
                                       CE_BMI >=25 & CE_BMI <=30 ~ 'overweight',
                                       CE_BMI >=30 ~ 'obese')) # end function
data$BMI <- as.numeric(as.factor((data$BMI)))
data$marital <- as.numeric(as.factor((data$marital)))
data$drinklot <- as.numeric(as.factor((data$drinklot)))
data$headinjuryever <- as.numeric(as.factor((data$headinjuryever)))
data$EDU <- as.numeric(as.factor((data$EDU)))

data<-filter(data, race =="Black")

## create R object with all exposures for regression models
xnm <- c('As', 'Cr', 'Cu', 'Hg', 'Mn','Pb','Se','Zn')
## create R object with all covariates
covars = c('CE_AGE', 'CE_C1', 'EN_FORMERSMOKER', 'marital','CE_BMI')

bkmr.data <- na.omit(data[,c(xnm, covars, 'CE_BARS_CPT_D_PRIME')])
z.data = as.matrix(bkmr.data[,c(xnm)])
x.data = as.matrix(bkmr.data[,c(covars)])
y.data= bkmr.data$CE_BARS_CPT_D_PRIME

set.seed(1234)
system.time(fitkm.model1 <- kmbayes_parallel(nchains=5, y = y.data, Z = z.data, X = x.data, iter = 40000, verbose = FALSE, varsel = TRUE))


## model diagnostics
multidiag = kmbayes_diagnose(fitkm.model1, warmup=20000, digits_summary=2)
fitkm.model1.coda <- as.mcmc.list(fitkm.model1, iterstart=20001)
trace_plot <- coda::traceplot(fitkm.model1.coda)
crosscorr_plot <- coda::crosscorr(fitkm.model1.coda)
autocorr_plot <- autocorr.plot(fitkm.model1.coda)
eff_size <- effectiveSize(fitkm.model1.coda) 
dens_plot <- densplot(fitkm.model1.coda)

#Posterior inclusion probabilities#
lapply(fitkm.model1, function(x) t(ExtractPIPs(x)))
fitkm.model1.comb = kmbayes_combine(fitkm.model1)
summary(fitkm.model1.comb)
fitkm.model1.pips<- ExtractPIPs(fitkm.model1.comb)
z_label <- as_labeller(c(z1='Al', z2='As', z3='Ca',z4='Cr',z5='Cu', z6='Fe', z7='Pb', z8='Mg', z9='Mn',z10='Hg',z11='Ni',z12='Se',z13='Zn'))

#Independent effects, holding all other exposures at 50th percentile (i.e., median)#
pred.resp.univar.data <- PredictorResponseUnivar(fit=fitkm.model1.comb, sel=seq(100001, 200000, by=1), ngrid=100, center=TRUE)

p1 <- ggplot(pred.resp.univar.data, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity", color = "black") + 
  facet_wrap(~ variable, ncol = 4) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black")) + 
  ylab("Outcome") + 
  xlab("Exposure")

#Bivariate associations, holding second exposure at different quantiles and all other exposures at 50th percentile (i.e., median)#
pred.resp.bivar.data <- PredictorResponseBivar(fit=fitkm.model1.comb, sel=seq(100001, 200000, by=1), ngrid=100, center=TRUE, min.plot.distance=1)

pred.resp.bivar.levels.data <- PredictorResponseBivarLevels(pred.resp.df = pred.resp.bivar.data, Z=z.data, qs=c(0.25, 0.50, 0.75))

p2 <- ggplot(pred.resp.bivar.levels.data, aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) + 
  geom_hline(yintercept=0, linetype= "dashed") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Outcome") +
  xlab("Exposure") 


#Overall mixture effect#
overall.estimates.5chain <- OverallRiskSummaries(fit=fitkm.model1.comb, y=y.data, Z=z.data, X=x.data, qs=seq(0.10, 0.90, by=0.05), q.fixed=0.50, method="approx", sel=seq(100001, 200000, by=1))
overall.estimates.5chain

p3 <- ggplot(overall.estimates.5chain, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange(size = 0.4) + 
  geom_hline(yintercept = 0.00, linetype = "dashed") + 
  theme_bw() + scale_y_continuous(limits=c(-1,1.3),breaks=c( -0.8, -0.4,  0, 0.4, 0.8,1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Quantile of Exposure") + 
  ylab("Difference in outcome") 


#Overall mixture effect, but decreasing/increasing one exposure while holding all other exposures at different quantiles#
risks.singvar <- SingVarRiskSummaries(fit=fitkm.model1.comb, y=y.data, Z=z.data, X=x.data, qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75), method="approx", sel=seq(100001, 200000, by=1))

p4 <- ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd, col =  q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept=0, linetype = "dashed") +
  ylab("Expected difference in outcome") +
  xlab("Exposure")

save.image(file = "BlackCPTdprimeall_results.RData")
