
########### -- Script setup
#setwd("C:/Users/cagger/Documents/Cornell WHI Acoustics/Gibbon illegal paper/Analysis/data")

library(dplyr)
library(stringr)
library(tidyr)
library(ggeffects)
library(ggplot2)
library(readr)
library(data.table)
library(lubridate)
library(lme4)
library(DHARMa)
library(nasapower)
library(performance)
library(glmmTMB)
library(car)


##################################### --- Does chainsaw activity have an effect on gibbons calls --- ##############################################



####################################### -- Pulling the data from selection tables

### Chainsaws


# Read all selection table locations in
Chainsaw_selection_tables <- data.frame(File = list.files(
  path = "data/S1139_Dep01_FLAC",
  recursive = TRUE, pattern = "txt$", full.names = TRUE
))

# Create a metadata sheet from selection tables

Chainsaw_selection_table_dataset <- Chainsaw_selection_tables %>%
  mutate(
    FileName = basename(File),
    Site = str_sub(FileName, 19, 21),
    Recorder = str_sub(FileName, 23, 27),
    SDCard = str_sub(FileName, 29, 34),
    Date = as.numeric(str_sub(FileName, 36, 43)),
    Hour = as.numeric(str_sub(FileName, 45, 46)),
    Month = as.numeric(str_sub(FileName, 40, 41)),
  )


# Create the dataframe from positive detection to be joined to the metadata

chainsaw_positive_files <- list.files(path = "data/Positive",
                                      pattern = ".wav", full.names = TRUE)

chainsaw_positive_dataset <- data.frame(File = basename(chainsaw_positive_files)) %>%
  mutate(
    Recorder = str_extract(File, "R\\d+"), # extract recorder
    SDCard = str_extract(File, "SD\\d+"), # extract SDcard
    Date = as.numeric(str_extract(File, "\\d{8}")), #extract date
    Hour = as.numeric(str_sub(str_extract(File, "_\\d{6}\\+"), 2, 3)), #extract hour
    Presence_chainsaw = 1 # create Presence_chainsaw when file exists
  )


# Join both metadata sheet and positive detection

chainsaw_dataset_together <- Chainsaw_selection_table_dataset %>%
  left_join(chainsaw_positive_dataset %>%
              select(Recorder, SDCard, Date, Hour, Presence_chainsaw),
            by = c("Recorder", "SDCard", "Date", "Hour")) %>%
  mutate(Presence_chainsaw = replace_na(Presence_chainsaw, 0))

chainsaw_dataset_together$Date <- as.Date(format(parse_date_time(chainsaw_dataset_together$Date, "%Y/%m/%d"), "%Y-%m-%d"))


### Gibbons

# Read in detections over 96% confidence interval that ave been manually verified

Gibbon_selection_tables <- read.csv("data/gibbonverifieddetections0.96.csv")

nrow(Gibbon_selection_tables) # 2438 manually verified detection

# Changing some filenames and column types so compatable with chainsaw dataset

Gibbon_selection_tables <- Gibbon_selection_tables %>%
  rename(Site = Plot)

Gibbon_selection_tables <- Gibbon_selection_tables %>%
  mutate(Gibbon_pres = as.numeric(c("1")))

Gibbon_selection_tables$Date <- as.Date(format(parse_date_time(Gibbon_selection_tables$Date, "%Y/%m/%d"), "%Y-%m-%d"))


################################################### -- Combining gibbon and chainsaw detections

# Want to combine so that Chainsaw is numbered by the amount of detections per day
# Gibbon can be treated as a yes of no call for the day, not sure in the accuracy when we are picking up multiple calls over a day
# chainsaws should be treated as frequency. Need to produce this frequency for every day, per transect
# In one dataset I need to retain the number of days and transects in a total list if chainsaw or gibbon present or not, day is the repeating factor/ repetition


# Collapsing chainsaw pres/abs columns so number of chainsaw events detected per day is the main determanent

Chainsaw_dataset_frequency <- chainsaw_dataset_together %>%
  group_by(Date, Site) %>%
  summarise(Chainsaw_frequency = sum(Presence_chainsaw))

# Joining gibbon dataset to chainsaw dataset, chainsaw as the main dataset as it contains all days listed out from the data collection period.
# We need all days listed here so that gibbon call is either 0 or 1 for that day, then have associated chainsaw detection frequecny for that day

Gibbon_chainsaw_data_combined <-
  left_join(Chainsaw_dataset_frequency, Gibbon_selection_tables,
            by = c("Date", "Site"))

# Replacing NA values in gibbon_pres with 0 as they did not call that day at that transect

Gibbon_chainsaw_data_combined <- mutate(Gibbon_chainsaw_data_combined, Gibbon_pres = replace_na(Gibbon_pres, 0))

# Collapse the dataset so their is only one repeat of each transect for each day
# At the moment there are multiple gibbons call detection per transect per day, all I need is presence and absence

Gibbon_chainsaw_data_combined <- Gibbon_chainsaw_data_combined %>%
  distinct(Date, Site, .keep_all = TRUE)

#Changing format of Date column to be a date
Gibbon_chainsaw_data_combined$Date <- as.Date(format(parse_date_time(Gibbon_chainsaw_data_combined$Date, "%Y/%m/%d"), "%Y-%m-%d"))


############################################################# -- Adding weather into combined dataset

### Weather data from weather station Jahoo

#Read in the data
weather_data <- read.csv("data/CUMULATIVE DATA_20230605-ongoing_Jahoo weather station.csv")

# subset for only the useful comments
weather_data <- subset(weather_data, select = c("Date", "Solar.Radiation", "Wind.Speed"))

# Getting date column into the format needed to subset
weather_data$Date <- as.Date(format(parse_date_time(weather_data$Date, "%d/%m/%Y"), "%Y-%m-%d"))

# subsetting to take out the days we don't need
weather_data <- weather_data[weather_data$Date > "2024-02-19" &
                               weather_data$Date < "2024-08-27",]

weather_data$Solar.Radiation <- as.numeric(weather_data$Solar.Radiation)
weather_data$Wind.Speed <- as.numeric(weather_data$Wind.Speed)

# Average the data over a day so it can be combined
weather_data <- weather_data %>%
  group_by(Date) %>%
  summarise(Solar.Radiation = mean(Solar.Radiation), Wind.Speed = mean(Wind.Speed))

# Joining weather data to main dataset

Gibbon_chainsaw_data_combined <-
  left_join(Gibbon_chainsaw_data_combined, weather_data,
            by = c("Date"))

### Weather data that we do not have from station has been downloaded from NASA

# https://docs.ropensci.org/nasapower/articles/nasapower.html -- Place where weather data was downloaded
#Below is the code used to fetch this data, but now it has been written, there is no need to run again
#Kept for archive purposes

# NASA_weather_variables <- get_power(
#  community = "ag",
#  lonlat = c(106.92, 12.15),
#  pars = c("RH2M", "T2M", "PRECTOTCORR"),
#  dates = c("2024-02-19", "2024-08-27"),
#  temporal_api = "daily"
# )
#
# write.csv(NASA_weather_variables, file = "data/NASA_weather_variables.csv")

NASA_weather_variables <- read.csv("data/NASA_weather_variables.csv")

## Combining weather data with full dataset

# Stripping down NASA so only have the things needed

NASA_weather_variables <- subset(NASA_weather_variables, select = c("YYYYMMDD", "RH2M", "T2M", "PRECTOTCORR"))

# Renaming date collumn in NASA dataset so can left join to main dataset

names(NASA_weather_variables)[names(NASA_weather_variables) == 'YYYYMMDD'] <- 'Date'

# Making the Date column in NASA variable as date format so it can left join

NASA_weather_variables$Date <- as.Date(format(parse_date_time(NASA_weather_variables$Date, "%Y-%m-%d"), "%Y-%m-%d"))

# Joinging to the full dataset

Gibbon_chainsaw_data_combined <-
  left_join(Gibbon_chainsaw_data_combined, NASA_weather_variables,
            by = c("Date"))


##################################################-- Look through detection history on each transect, remove transect where gibbons were never detected


Gibbon_chainsaw_data_combined <- Gibbon_chainsaw_data_combined %>%
  group_by(Site) %>%
  filter(sum(Gibbon_pres) >= 1)

# Apply correction only to T19
Gibbon_chainsaw_data_combined <- Gibbon_chainsaw_data_combined %>%
  mutate(
   Time = if_else(
      Site == "T19",
      as.numeric(Time) - 13,
      as.numeric(Time)
    )
  )

Gibbon_chainsaw_data_combined$Time <-
  if_else(Gibbon_chainsaw_data_combined$Time < 0, Gibbon_chainsaw_data_combined$Time + 24, Gibbon_chainsaw_data_combined$Time)


########## -- Data exploration

# Basic stats and figures -- Of the recorders where gibbons were detected

nrow(Gibbon_chainsaw_data_combined) # 2286 recorder days
sum(Gibbon_chainsaw_data_combined$Gibbon_pres) # 367 days that include gibbon call
sum(Gibbon_chainsaw_data_combined$Chainsaw_frequency == 0) # 1724 recorder days without chainsaw noise
min(Gibbon_chainsaw_data_combined$Date) # 20240220 is the first date of the dataset
max(Gibbon_chainsaw_data_combined$Date) # 0240827 the last date of the dataset

# Basic graphs

### The timing of chainsaw noise by hour in a day

chainsaw_noise_per_hour <- chainsaw_dataset_together %>%
  group_by(Hour) %>%
  summarise(chainsaw_count = sum(Presence_chainsaw, na.rm = T))

ggplot(chainsaw_noise_per_hour, aes(x = Hour, y = chainsaw_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Hour of day", y = "Frequency of chainsaw detections",
       title = "Distribution of chainsaw detections per hour") +
  theme_minimal()

### The timing of gibbon calls by hour in a day

gibbon_call_per_hour <- Gibbon_chainsaw_data_combined %>%
  group_by(Time) %>%
  summarise(gibbon_call_count = sum(Gibbon_pres, na.rm = T))

ggplot(gibbon_call_per_hour, aes(x = Time, y = gibbon_call_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Hour of day", y = "Frequency of gibbon call detections",
       title = "Distribution of gibbon call detections per hour") +
  theme_minimal()

chainsaw_noise_per_transect <- chainsaw_dataset_together %>%
  group_by(Site) %>%
  summarise(chainsaw_count = sum(Presence_chainsaw, na.rm = T))

ggplot(chainsaw_noise_per_transect, aes(x = Site, y = chainsaw_count)) +
  geom_bar(stat = "identity", fill = "red") +
  labs(x = "Transect number", y = "Frequency of chainsaw detections",
       title = "Distribution of chainsaw detections per transect") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


## The amount of gibbon calls per transect

gibbon_call_per_transect <- Gibbon_selection_tables %>%
  group_by(Site) %>%
  summarise(gibbon_call_count = sum(Gibbon_pres, na.rm = T))

ggplot(gibbon_call_per_transect, aes(x = Site, y = gibbon_call_count)) +
  geom_bar(stat = "identity", fill = "red") +
  labs(x = "Site", y = "Frequency of gibbon call detections",
       title = "Distribution of gibbon call detections per hour") +
  theme_minimal()

#write.csv(Gibbon_chainsaw_data_combined,'data/Gibbon_chainsaw_data_combined.csv')
########################## ---------- GLMM modelling

##Coding the model

gibbon_chainsaw_model_1 <- glmmTMB(Gibbon_pres ~ Chainsaw_frequency + (Solar.Radiation + Wind.Speed + PRECTOTCORR + RH2M | Site), family = binomial, data = Gibbon_chainsaw_data_combined)
summary(gibbon_chainsaw_model_1)

test <- simulateResiduals(gibbon_chainsaw_model_1) # Check the distribution of the residuals
plot(test)

check_overdispersion(gibbon_chainsaw_model_1) # If overdispersed, can use negative binomial

## P-value testing

null_model <- glmmTMB(Gibbon_pres ~ 1 + (Solar.Radiation + Wind.Speed + PRECTOTCORR + RH2M | Site), family = binomial, data = Gibbon_chainsaw_data_combined)

anova(null_model, gibbon_chainsaw_model_1, test = 'LRT')

## Quickly predicting and plotting

#create dataframe to predict onto
basic_test_predcting <- data.frame(sort(unique(Gibbon_chainsaw_data_combined$Chainsaw_frequency)))

basic_test_predcting <- basic_test_predcting %>%
  rename(Chainsaw_frequency = sort.unique.Gibbon_chainsaw_data_combined.Chainsaw_frequency..)

#running the prediction
basic_test_predcting$predicted <- predict(gibbon_chainsaw_model_1, basic_test_predcting, type = "response", re.form = NA)

# Very quickly plotting --

ggplot(basic_test_predcting, aes(x = Chainsaw_frequency, y = predicted)) +
  geom_line(stat = "identity") +
  labs(x = "Number of chainsaw detections within 1 hour bins in a day", y = "Probability of gibbon calling on a given morning")+
  theme_minimal()


#################################### --------------------- Investigating the random effects structure and if it is colinear

## Is any of the weather related predictors co-linear with chainsaw activity
# Use pearson correlation to see if correlated. +1 or -1 indicates strong correlation
# Values close to 0 indicate no linear relationship
cor.test(Gibbon_chainsaw_data_combined$Chainsaw_frequency, Gibbon_chainsaw_data_combined$Solar.Radiation, use = "complete.obs", method = "pearson")
cor.test(Gibbon_chainsaw_data_combined$Chainsaw_frequency, Gibbon_chainsaw_data_combined$Wind.Speed, use = "complete.obs", method = "pearson")
cor.test(Gibbon_chainsaw_data_combined$Chainsaw_frequency, Gibbon_chainsaw_data_combined$PRECTOTCORR, use = "complete.obs", method = "pearson")
cor.test(Gibbon_chainsaw_data_combined$Chainsaw_frequency, Gibbon_chainsaw_data_combined$RH2M, use = "complete.obs", method = "pearson")

plot(Gibbon_chainsaw_data_combined$Chainsaw_frequency, Gibbon_chainsaw_data_combined$Solar.Radiation, main = "Scatterplot", xlab = "Var1", ylab = "Var2")
abline(lm(Solar.Radiation ~ Chainsaw_frequency, data = Gibbon_chainsaw_data_combined), col = "red", lwd = 2)

plot(Gibbon_chainsaw_data_combined$Chainsaw_frequency, Gibbon_chainsaw_data_combined$Wind.Speed, main = "Scatterplot", xlab = "Var1", ylab = "Var2")
abline(lm(Wind.Speed ~ Chainsaw_frequency, data = Gibbon_chainsaw_data_combined), col = "red", lwd = 2)

plot(Gibbon_chainsaw_data_combined$Chainsaw_frequency, Gibbon_chainsaw_data_combined$PRECTOTCORR, main = "Scatterplot", xlab = "Var1", ylab = "Var2")
abline(lm(PRECTOTCORR ~ Chainsaw_frequency, data = Gibbon_chainsaw_data_combined), col = "red", lwd = 2)

plot(Gibbon_chainsaw_data_combined$Chainsaw_frequency, Gibbon_chainsaw_data_combined$RH2M, main = "Scatterplot", xlab = "Var1", ylab = "Var2")
abline(lm(RH2M ~ Chainsaw_frequency, data = Gibbon_chainsaw_data_combined), col = "red", lwd = 2)


###### Looking at the output of this first model and making notes

##### Need to double check all weather variables to make sure they are there all the time and there and no NAs

# Check the random effects summary
# Very small Standard deviations suggests these random effects may not be contrbuting: PRECTOTCORR is very small, potentially remove
# Correlations in the output
# The interpretation of the random intercept shows that there's substantial variability in baseline gibbon presence across sites.

# Solar raditation -- variance is tiny, high correlation with the intercept suggests a computational artifact, remove
# Wind speed -- Moderate variation in how wind speed affects gibbon pres, correlation with other effects is weak, keep
# PRECTOTCORR -- Esentially 0 variance, very high correlation likely due to numerical instability, consider removing.
# RH2M relative humidity - small but none 0 variance, very high negative correlation with other random effects suggests collinearity or overparamatization, remove.

# Look at the r2 to see how much variation the radom effects explains. If it's alot keep as is
r2(gibbon_chainsaw_model_1)## This is suggesting that chainsaw freq is explaining virtually none of the variation in gibbon pres, need to try a leaner model

# Quickly test for singularity
check_singularity(gibbon_chainsaw_model_1) ## Model is singular, meaning we are trying to estimate more parameters than the data can support

############## --------------- Improved GLM modelling based on what has been learned above

# Althought the conclusion of the random effects above is saying to take out rain, I think ecologically it is important so maybe try and keep in for now
# Or investigate the variable and see how it changes over the days, plot is !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Or find another data source to get this data from that might be more accurate!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#NEED TO DO THE ABOVE

gibbon_chainsaw_model_2 <- glmmTMB(Gibbon_pres ~ Chainsaw_frequency + (Wind.Speed | Site), family = binomial, data = Gibbon_chainsaw_data_combined)
summary(gibbon_chainsaw_model_2)




########################################### -------- Timing of calls analysis ----------###############################################################

# https://onlinelibrary.wiley.com/doi/10.1111/btp.13205 -- one paper that might help with this analysis. Good approach, but might have to adapt

# - Other research into how we can emasure timing

# https://www.sciencedirect.com/science/article/pii/S0895435622001779#:~:text=Contemporaneous%20effects%20are%20the%20effects,%2Dtwo%20=%20t%20+%202.

################################################# -------------- CANDIDATE MODEL APPROACH WITH TIME BINS

#### Wrangling data into the correct format.

# Renaming date in gibbon data so can join datasets correctly

Gibbon_selection_tables_timing_analysis <- Gibbon_selection_tables %>%
  rename(Hour = Time)

Gibbon_selection_tables_timing_analysis$Month_gib <- format(Gibbon_selection_tables_timing_analysis$Date, "%m")

# Need to collapse gibbon calling dataset so only includes one detection per hour
# At the moment whenever there is a detection it is in the data, meaning there are multiple detections from the same hour

Gibbon_selection_tables_timing_analysis <- Gibbon_selection_tables_timing_analysis %>%
  distinct(Date, Hour, Month_gib, .keep_all = TRUE)

# Merging chainsaw and gibbon calls data together so it is organised hourly

GLMM_organised_hourly_data <- chainsaw_dataset_together %>%
  left_join(Gibbon_selection_tables_timing_analysis %>%
              select(Site, Date, Hour, Gibbon_pres),
            by = c("Site", "Date", "Hour")) %>%
  mutate(Gibbon_pres = replace_na(Gibbon_pres, 0))

# Removing transect deployment locations where gibbons were never detected

GLMM_organised_hourly_data <- GLMM_organised_hourly_data %>%
  group_by(Site) %>%
  filter(sum(Gibbon_pres) >= 1)

# Check for zero inflation in the dataset -  Likely to change when get new gibbon verified data from Dena

GLMM_organised_hourly_data <- GLMM_organised_hourly_data %>%
  ungroup()

GLMM_organised_hourly_data %>%
  count(Presence_chainsaw) #### 0 = 49737, 1 = 1871

GLMM_organised_hourly_data %>%
  count(Gibbon_pres) #### 0 = 51295, 1 = 313


########################----------------------- Simple check to see if the model is overparamatised without binning hour

Chainsaw_and_hour_interaction_model <- glmmTMB(Gibbon_pres ~ Presence_chainsaw * Hour + (1 | Site), family = binomial, data = GLMM_organised_hourly_data)
summary(Chainsaw_and_hour_interaction_model)

### --- Set the reference level of the interaction and set factors so model knows what type of data it is running

# Reference of chainsaw presence is non/0
GLMM_organised_hourly_data$Presence_chainsaw <- factor(GLMM_organised_hourly_data$Presence_chainsaw, levels = c("0", "1"))

# Treating hour as a factor as I expect nonlinear effect between hours and gibbon calls
GLMM_organised_hourly_data$Hour <- factor(GLMM_organised_hourly_data$Hour)

# Run the model
Chainsaw_and_hour_interaction_model_factors <- glmmTMB(Gibbon_pres ~ Presence_chainsaw * Hour + (1 | Site), family = binomial, data = GLMM_organised_hourly_data)
summary(Chainsaw_and_hour_interaction_model_factors)

### Have alot of NA in the output,
#might be because there is insufficient data in the combination of factors,
# HAve alook at this

table(GLMM_organised_hourly_data$Presence_chainsaw, GLMM_organised_hourly_data$Hour)#  Doesn't seem to be this, each hour contains chainsaws

### Might also be because of overparatisation, too many hours in a day in the model. Group together to get around this?
#Too many parameters vs observations, if this dataset contains more columns than rows, model could be overparamatised

X <- model.matrix(Chainsaw_and_hour_interaction_model_factors)
dim(X) # first number is number of rows, second number is the number of columns. Doesn't look like this is the problem here

### MODEL IS OVERPARAMATIZED.

########################################## ----- Time bin the time of days so the model actually works normally

# manipulating the data to bin hours into the correct format
### Have binned to 4 at this time, but maybe if we do more manual verification we up the number of detections and provide more data to produce smaller bins/more complicated models

GLMM_organised_hourly_data$Hour <- as.numeric(GLMM_organised_hourly_data$Hour)

GLMM_4hour_time_bins_data <- GLMM_organised_hourly_data %>%
  mutate(Hour_bin = cut(Hour,
                        breaks = seq(0, 24, by = 4),
                        labels = c("00–03", "04–07", "08–11", "12–16",
                                   "17–20", "21–24"),
                        right = FALSE))

## Try to model again but with hours into 3 hourly time bins.
# The middle of the night (0-2) is set as the level of comparison here

######## Run the candidate models with the 3 hour timebins

######################################################## ------------ Candidate model situation approach with time bins


#### Data wrangling into the correct format for the analysis

####### dependent variable is binomial if gibbons called in the morning or not. OR THE HOUR OF THEIR CALL?

###### Candidate model 1: Null model of only random effect of recorder.

Candidate_model_1 <- glmmTMB(Gibbon_pres ~ 1 + (1 | Site), family = binomial, data = GLMM_4hour_time_bins_data)
summary(Candidate_model_1)

###### Candidate model 2: random effect of recorder, recorder as predictors.
###### Can then look and see if there is a difference of timing of calls across recorders.
###### Not including hour variable in here means model just looks at each repeat, but doesn't know anything about the time it is happening,
###### Just seeing if recorder location is pushng the trend of when gibbons are calling

Candidate_model_2 <- glmmTMB(Gibbon_pres ~ Site, family = binomial, data = GLMM_4hour_time_bins_data)
summary(Candidate_model_2)

###### Candidate model 3: random effect of recorder and fixed effect of Hour.
###### This model is the check if the hour of calling impacts the most

Candidate_model_3 <- glmmTMB(Gibbon_pres ~ Hour_bin + (1 | Site.x), family = binomial, data = GLMM_4hour_time_bins_data)
summary(Candidate_model_3)

###### Candidate model 4: fixed effects of hour and chainsaw pres/abs, and random effect of recorder.
###### Reflects the prediction there will be differences in gibbon call detection by time of day but also by chainsaw pres/abs

Candidate_model_4 <- glmmTMB(Gibbon_pres ~ Presence_chainsaw + Hour_bin + (1 | Site), family = binomial, data = GLMM_4hour_time_bins_data)
summary(Candidate_model_4)

###### Candidate model 5: Interaction between hour and chainsaw detection pres/abs, as well as random effect of recorder location nested within Date .
###### test the prediction that there was an effect of the timing of chainsaw events
###### When plot see what the predicted gibbon call probability is depending on hour of day.
###### example graph in notes

Candidate_model_5 <- glmmTMB(Gibbon_pres ~ Presence_chainsaw * Hour_bin + (1 | Site), family = binomial, data = GLMM_4hour_time_bins_data)
summary(Candidate_model_5)

# Checking model assumptions

check_overdispersion(Candidate_model_5)

test_hour_bins <- simulateResiduals(Candidate_model_5) # Check the distribution of the residuals
plot(test_hour_bins)

# Checking model significance using LRT

anova(Candidate_model_1, Candidate_model_5, test = 'LRT')

# Predicting and plotting, dataset to predict on

dataframe_for_preds <- expand.grid(
  Presence_chainsaw = c(0, 1),
  Hour_bin = levels(GLMM_4hour_time_bins_data$Hour_bin)
)

# Defining bootstrap function for the prediction
# Number of bootsrap iterations
n_boot <- 100

# Matrix to store the predicted values
pred_matrix_hour_bins <- matrix(NA, nrow = n_boot, ncol = nrow(dataframe_for_preds))

set.seed(123) # For reproduceability

# Bootstrapping function
for (i in 1:n_boot) {
  #Resample the origional data with replacement
  boot_data <- GLMM_3hour_time_bins_data[sample(nrow(GLMM_3hour_time_bins_data), replace = TRUE), ]
  #Refit the model on bootstrap sample
  Candidate_model_5_boot <- glmmTMB(Gibbon_pres ~ Presence_chainsaw * Hour_bin +
                                                           (1 | Site), family = binomial, data = boot_data)
  # Predict and exclude random effects
  pred_matrix_hour_bins[i, ] <- predict(Candidate_model_5_boot,
                                        newdata = dataframe_for_preds, type = "response", re.form = NA)
}

# Summarise confidence intervals, calculate them and then add into the dataset

pred_mean <- apply(pred_matrix_hour_bins, 2, mean)
pred_lwr  <- apply(pred_matrix_hour_bins, 2, quantile, probs = 0.025)
pred_upr  <- apply(pred_matrix_hour_bins, 2, quantile, probs = 0.975)
dataframe_for_preds$pred_mean <- pred_mean
dataframe_for_preds$lwr_CI  <- pred_lwr
dataframe_for_preds$upr_CI  <- pred_upr

# Plot the predictons and CIs


ggplot(dataframe_for_preds, aes(x = Hour_bin, y = pred_mean, fill = as.factor(Presence_chainsaw))) +
  geom_col(position = position_dodge(0.8)) +
  geom_errorbar(position = position_dodge(0.8), aes(ymin = lwr_CI, ymax = upr_CI), width = 0.2) +
  labs(y = "Predicted probability of gibbon calling within a given day", x = "Hour Bin", fill = "Chainsaw Present") +
  theme_minimal()



############################################################ ------------ Restricting the timeframe of the model

#### Data wrangling
#Checking which time the gibbons call in the dataset so I can restrict the analysis just to those times they have called.

calls_x_hour <- GLMM_organised_hourly_data %>%
  group_by(Gibbon_pres, Hour) %>%
  summarise(count = n())

# Mainly call at 5,6,7,8,9, so include this in the analysis timeframe

Time_restricted_dataset <- GLMM_organised_hourly_data %>%
  filter(Hour %in% c("5", "6"))

Time_restricted_dataset$Hour <- factor(Time_restricted_dataset$Hour)
Time_restricted_dataset$Gibbon_pres <- factor(Time_restricted_dataset$Gibbon_pres)

#### Modelling framework

Time_restricted_model <- glmmTMB(Gibbon_pres ~ Presence_chainsaw * Hour + (1 | Site), family = binomial, data = Time_restricted_dataset)
summary(Time_restricted_model)

### No statistical significance here


############################################################ ------------ Impact of chainsaw noise from the day before

### Data wrangling

GLMM_organised_hourly_day_before <- GLMM_organised_hourly_data %>%
  mutate(Date = as.Date(Date),
         Day_before = Date - 1)

# Use the new day before column as an index to pull across if chainsaw was present at that time.
# Can say look at date and time column and pull across if gibbon_pres was 1 or 0 for that time the day before

GLMM_organised_hourly_day_before <- GLMM_organised_hourly_day_before %>%
  left_join(
    GLMM_organised_hourly_day_before %>%
      select(Date, Hour, Presence_chainsaw, Site) %>%
      rename(Chainsaw_prev_day = Presence_chainsaw),
    by = c("Day_before" = "Date", "Hour" = "Hour", "Site" = "Site")
  )

#Organise chainsaw presence day before into time bins or if it was present or not. I.e the evening before. Can categories as the time it is dark
# Do the analysis the same as the time bins above? Just for the day before?
# Put into bins again, I think it depends on if I have enough data or not. Need to do more verfications

