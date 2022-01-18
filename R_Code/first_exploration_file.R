getwd()
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
age_real <- data$age
day_real <- data$day
data$age <- data$age - mean(data$age)
data$day <- data$day - mean(data$day)
saveRDS(data, file = "data\\AIdataset_normalized.Rds")
#### some Plots ####
library(ggplot2)
ggplot(data,aes(x=day_real,y=sofa)) + stat_binhex(na.rm = T) + scale_fill_gradient(low = "blue", high = "green")
ggplot(data,aes(x=day_real,y=sofa)) + geom_point(aes(colour = id)) +guides(colour="none")
ggplot(data,aes(x=day_real,y=sofa)) + geom_point(aes(colour = id)) + 
  guides(colour="none")+geom_text(mapping = aes(label = id),hjust=0, vjust=0)


ggplot(data,aes(x=age_real,y=sofa)) + stat_binhex(na.rm = T) + scale_fill_gradient(low = "blue", high = "green")
ggplot(data,aes(x=age_real,y=sofa)) + geom_point(aes(colour = id)) +guides(colour="none")
ggplot(data,aes(x=age_real,y=sofa)) + geom_point(aes(colour = id)) + 
  guides(colour="none")+geom_text(mapping = aes(label = id),hjust=0, vjust=0)

####Frequentist analysis ####
library(lme4)

sofa_freq <- lmer(sofa ~ age + day +(1|id), data, REML = T)
summary(sofa_freq)
