#Data Preparation

data <- readRDS("data\\AIdataset.Rds")
#MV   Mechanical Ventilation
#ICU  intensive care unit
#PVA  Patient-ventilator asynchronies
#SOFA measures the degree of organ dysfunction and assesses 
#the evolution of the patient’s severity over the ICU stay
#AI  measures PVAs and assesses overall severity
#### Data Changes ####
#ai is somehow called di
####Initial Changes####
colnames(data)[6] <- "ai"
library(dplyr)
#get the biggest value for every id
icu_days_obs <- data %>%  group_by(id) %>%  summarise(icu_days_obs = max(day)) 
#mean max min all the same, since it is the same value for the same id
icu <- data %>%group_by(id) %>%summarise(icu_days_said = mean(icustaymv))
#in theory those values should be the same but they are not. So we assume that
#the people just joined the study after they where already in the MV unit for a 
#while.

icu$icu_days_obs <- icu_days_obs$icu_days_obs
icu$difference <-icu$icu_days_said - icu$icu_days_obs
#sometimes the observed days are higher as the proposed ones, I just said that then
#the icu value was one to short and so I basically said icu is one day bigger
#else - shift in the data really not needed
icu$difference <- ifelse(icu$difference == -1 , 0, icu$difference )
u <- as.numeric(data$id)
data$day <- data$day + icu$difference[u]
age_real <- data$age
day_real <- data$day
#now again normalize the data as before 
data$age <- data$age - mean(data$age)
data$day <- data$day - mean(data$day)
saveRDS(data, file = "data\\AIdataset_normalized.Rds")
