# Load libraries ----------------------------------------------------------

library(data.table)
library(lubridate)
library(lme4)
library(DHARMa)
library(nasapower)
library(performance)
library(glmmTMB)
library(car)
library(bbmle)
library(dplyr)
library(stringr)
library(tidyr)
library(ggeffects)
library(ggplot2)
library(readr)



# Read in data ------------------------------------------------------------

Gibbon_chainsaw_data_combined <- read.csv('data/Gibbon_chainsaw_data_combined.csv')


# Modeling script ---------------------------------------------------------

str(Gibbon_chainsaw_data_combined)

 # Gibbon_chainsaw_data_combined <-
 #   subset(Gibbon_chainsaw_data_combined,Time < 7)

Gibbon_chainsaw_data_combined <- Gibbon_chainsaw_data_combined %>%
  mutate(Chainsaw_binary = if_else(Chainsaw_frequency > 0, 1, 0))

##Coding the model

gibbon_chainsaw_model_1 <- glmmTMB(Gibbon_pres ~
                                     Chainsaw_binary+
                                     (1| Site),
                                   family = binomial,
                                   data = Gibbon_chainsaw_data_combined)
summary(gibbon_chainsaw_model_1)


gibbon_bin_model <- glmmTMB(Gibbon_pres ~ Chainsaw_binary + PRECTOTCORR + (1 | Site),
                            family = binomial,
                            data = Gibbon_chainsaw_data_combined)

summary(gibbon_bin_model)

test <- simulateResiduals(gibbon_bin_model) # Check the distribution of the residuals
#plot(test)

check_overdispersion(gibbon_chainsaw_model_1) # If overdispersed, can use negative binomial

## P-value testing

null_model <- glmmTMB(Gibbon_pres ~  (1| Site), family = binomial, data = Gibbon_chainsaw_data_combined)

full_model <- glmmTMB(Gibbon_pres ~ Chainsaw_binary + PRECTOTCORR +
                        Wind.Speed + RH2M+ T2M+
                      + (1 | Site),
                            family = binomial,
                            data = Gibbon_chainsaw_data_combined)

AICtab(null_model,
       gibbon_chainsaw_model_1,
       gibbon_bin_model,
       full_model,
       weights=T)

sjPlot::plot_model(full_model,type='pred')
