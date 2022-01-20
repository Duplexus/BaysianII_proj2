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

#If observations are meaned by person
library(dplyr)
data <- readRDS("data\\AIdataset.Rds")
data_by_id <- data %>%  group_by(id) %>% summarise(sofa = mean(sofa, na.rm = T), day = mean(day, na.rm = T), age = mean(age))

ggplot(data_by_id,aes(x=day,y=sofa)) + stat_binhex(na.rm = T) + scale_fill_gradient(low = "blue", high = "green")
ggplot(data_by_id,aes(x=day,y=sofa)) + geom_point() +guides(colour="none")
ggplot(data_by_id,aes(x=day,y=sofa)) + geom_point() + 
  guides(colour="none")+geom_text(mapping = aes(label = age),hjust=0, vjust=0)

ggplot(data_by_id,aes(x=age,y=sofa)) + stat_binhex(na.rm = T) + scale_fill_gradient(low = "blue", high = "green")
ggplot(data_by_id,aes(x=age,y=sofa)) + geom_point(aes(colour = id)) +guides(colour="none")
ggplot(data_by_id,aes(x=age,y=sofa)) + geom_point(aes(colour = id)) + 
  guides(colour="none")+geom_text(mapping = aes(label = id),hjust=0, vjust=0)

#grouped by day, therefore evolution age and sofa with spline line
library(splines)
data_by_day <- data %>%  group_by(day) %>% summarise(sofa = mean(sofa, na.rm = T), age = mean(age))

spline_mod1 <- lm(sofa ~ bs(day, knots = c(10,20)),data_by_day)
y_sp_mod1 <- predict(spline_mod1, data.frame(day = c(1:30)))

ggplot(data_by_day,aes(x=day,y=sofa)) + geom_point() + 
  geom_line(aes(x = c(1:30)), y = y_sp_mod1 )

#grouped by age, therefore evolution age and sofa with spline line
data_by_age <- data %>%  group_by(age) %>% summarise(sofa = mean(sofa, na.rm = T), day = mean(day))

summary(data_by_age$age)

spline_mod2 <- lm(sofa ~ bs(age, knots = c(20,40,60)),data_by_age)
age_ests <- seq(17,88,length.out = length(data_by_age$age))
y_sp_mod2 <- predict(spline_mod2, data.frame(age = age_ests))

ggplot(data_by_age,aes(x=age,y=sofa)) + geom_point() + 
  geom_line(aes(x = age_ests), y = y_sp_mod2 )





#same again this time for ai as response
#### some Plots ####
data <- readRDS("data\\AIdataset.Rds")
####Initial Changes####
age_real <- data$age
day_real <- data$day
data$age <- data$age - mean(data$age)
data$day <- data$day - mean(data$day)
library(ggplot2)
library(lattice)
xyplot(data$ai ~ data$day|as.factor(data$id), cex=0.5,
       par.settings=list(superpose.symbol = list(pch=c(1,3,20))),
       as.table=T, auto.key=list(points=T,columns=3), xlab="day", ylab="SOFA")

xyplot(data$ai ~ data$age|as.factor(data$day), cex=0.5,
       par.settings=list(superpose.symbol = list(pch=c(1,3,20))),
       as.table=T, auto.key=list(points=T,columns=3), xlab="day", ylab="AI")




ggplot(data,aes(x=day_real,y=ai)) + stat_binhex(na.rm = T) + scale_fill_gradient(low = "blue", high = "green")
ggplot(data,aes(x=day_real,y=ai)) + geom_point(aes(colour = id)) +guides(colour="none")
ggplot(data,aes(x=day_real,y=ai)) + geom_point(aes(colour = id)) + 
  guides(colour="none")+geom_text(mapping = aes(label = id),hjust=0, vjust=0)


ggplot(data,aes(x=age_real,y=ai)) + stat_binhex(na.rm = T) + scale_fill_gradient(low = "blue", high = "green")
ggplot(data,aes(x=age_real,y=ai)) + geom_point(aes(colour = id)) +guides(colour="none")
ggplot(data,aes(x=age_real,y=ai)) + geom_point(aes(colour = id)) + 
  guides(colour="none")+geom_text(mapping = aes(label = id),hjust=0, vjust=0)

#If observations are meaned by person
library(dplyr)
data <- readRDS("data\\AIdataset.Rds")
colnames(data)[6] <- "ai"

data_by_id <- data %>%  group_by(id) %>% summarise(ai = mean(ai, na.rm = T), day = mean(day, na.rm = T), age = mean(age))

ggplot(data_by_id,aes(x=day,y=ai)) + stat_binhex(na.rm = T) + scale_fill_gradient(low = "blue", high = "green")
ggplot(data_by_id,aes(x=day,y=ai)) + geom_point() +guides(colour="none")
ggplot(data_by_id,aes(x=day,y=ai)) + geom_point() + 
  guides(colour="none")+geom_text(mapping = aes(label = age),hjust=0, vjust=0)

ggplot(data_by_id,aes(x=age,y=ai)) + stat_binhex(na.rm = T) + scale_fill_gradient(low = "blue", high = "green")
ggplot(data_by_id,aes(x=age,y=ai)) + geom_point(aes(colour = id)) +guides(colour="none")
ggplot(data_by_id,aes(x=age,y=ai)) + geom_point(aes(colour = id)) + 
  guides(colour="none")+geom_text(mapping = aes(label = id),hjust=0, vjust=0)

#grouped by day, therefore evolution age and ai with spline line
library(splines)
data_by_day <- data %>%  group_by(day) %>% summarise(ai = mean(ai, na.rm = T), age = mean(age))

spline_mod1 <- lm(ai ~ bs(day, knots = c(10,20)),data_by_day)
y_sp_mod1 <- predict(spline_mod1, data.frame(day = c(1:30)))

ggplot(data_by_day,aes(x=day,y=ai)) + geom_point() + 
  geom_line(aes(x = c(1:30)), y = y_sp_mod1 )

#grouped by age, therefore evolution age and ai with spline line
data_by_age <- data %>%  group_by(age) %>% summarise(ai = mean(ai, na.rm = T), day = mean(day))

summary(data_by_age$age)

spline_mod2 <- lm(ai ~ bs(age, knots = c(20,40,60)),data_by_age)
age_ests <- seq(17,88,length.out = length(data_by_age$age))
y_sp_mod2 <- predict(spline_mod2, data.frame(age = age_ests))

ggplot(data_by_age,aes(x=age,y=ai)) + geom_point() + 
  geom_line(aes(x = age_ests), y = y_sp_mod2 )




####Frequentist analysis ####
library(lme4)

sofa_freq <- lmer(sofa ~ age + day +(1|id), data, REML = T)
summary(sofa_freq)
