### Stats analysis for clutch and hatchling level phenotypes ###

# Created by Chris, edited by Charley, 
# 16th May 2024
# For Evol Apps

########################################################################################################################################

###############################
###### Prep environment #######
###############################

### Load packages ###
library(ggplot2)
library(lme4)
library(lmerTest)
library(MASS)
library(car)
library(emmeans)

### Set directories ###
setwd("/Users/ecyen/Dropbox/Documents/PhD/Project/00_Scripts/04_Peoples_Code/Hatchery/Chris_EvolApps_Reloc3")

### Load data sets ###

# Clutch level dataset
data2<-read.csv("Reloc3_NestLevelData_WithTemperatureInfo_NERCExp2021.csv")

# Hatchling level dataset
data<-read.csv("Reloc3_Hatchling_Fitness_NERC_Exp_2021_WithPatrolDates_WithIDs_T.csv")

########################################################################################################################################

#################################
###### Clutch-level stats #######
#################################

# Retain data with temperature information to merge in models
remove_na_values <- function(dataset, column_name) {
  # Remove NA values from the specified column
  cleaned_dataset <- dataset[!is.na(dataset[[column_name]]), ]
  
  # Return the cleaned dataset
  return(cleaned_dataset)
}

data2 <- remove_na_values(data2, "Full_Incubation_MeanTemp")

################################################
###### (A) Linear model to test if mean sub-clutch incubation temperature is associated with depth treatment ######

model = lm(data2$Full_Incubation_MeanTemp ~ data2$Treatment)
summary(model)
anova(model)

### Plots ###
# Plot of temperature vs treatment 
ggplot(data2, aes(x = Treatment, y = Full_Incubation_MeanTemp, fill = Treatment)) + 
  geom_boxplot()  +
  labs(x="Treatment", y="Mean incubation temperature") +
  scale_fill_manual(values=c("#FA9E3B","#0D0887")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_classic() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))

# Plot of temperature vs clutch size
# To show known effect of metabolic heating
ggplot(data2, aes(x = Split_Clutch_Size, y = Full_Incubation_MeanTemp, color = Treatment)) + 
  geom_point(size=2)  +
  geom_smooth(method = "lm") +
  labs(x="Sub-clutch size (eggs)", y="Mean incubation temperature") +
  scale_color_manual(values=c("#FA9E3B","#0D0887")) +
  theme_classic() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))

model = lm(data2$Full_Incubation_MeanTemp ~ data2$Treatment*data2$Split_Clutch_Size)
summary(model)
anova(model) # Effects of treatment and temperature, but both independent from one another

################################################
###### (B) Linear model to test if mean sub-clutch incubation duration is associated with depth treatment ######

# Testing incubation duration against temperature to verify classic assumption of a correlation
model = lm(data2$Incubation_Duration ~ data2$Treatment*data2$Full_Incubation_MeanTemp)
summary(model)
anova(model) # Effects of treatment and temperature, but both independent from one another

### Plots ###
# Plot of incubation duration vs treatment 
ggplot(data2, aes(x = Treatment, y = Incubation_Duration, fill = Treatment)) + 
  geom_boxplot()  +
  labs(x="Treatment", y="Incubation duration (days)") +
  scale_fill_manual(values=c("#FA9E3B","#0D0887")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_classic() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))

# Plot incubation duration vs temperature 
ggplot(data2, aes(x = Full_Incubation_MeanTemp, y = Incubation_Duration, color = Treatment)) + 
  geom_point(size=2)  +
  geom_smooth(method = "lm") +
  labs(x="Mean incubation temperature", y="Incubation duration (days)") +
  scale_color_manual(values=c("#FA9E3B","#0D0887")) +
  theme_classic() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))

################################################
###### (C) Linear model to test for correlates of hatching success ######

# Calculate hatching success:total split-clutch size minus total dead hatchlings, divided by total split clutch size
data2["Success_Rate"] <- (data2$Split_Clutch_Size - (data2$Dead_Full_Term+data2$Dead_Pip+data2$Unhatched_Egg_No_Embryo+data2$Unhatched_Egg_Embryo+data2$Unhatched_Egg_Fully_Formed_Embryo))/data2$Split_Clutch_Size

# Calculate clutch size squared:
# from visual inspection there is a quadratic relationship, which helps us to formulate model (see plot)
data2$Clutch_Size2=data2$Clutch_Size^2 

## Use residuals of a linear model between incubation temperature and treatment (model A), rather than incubation temperature alone
# to disentangle effects of temperature and treatment on response
mod=lm(data2$Full_Incubation_MeanTemp ~ data2$Treatment)

# Full, starting model
model=lm(data2$Success_Rate ~ data2$Treatment + mod$residuals + data2$Clutch_Size + data2$Clutch_Size2 +
           mod$residuals:data2$Treatment + data2$Clutch_Size:data2$Treatment + data2$Clutch_Size2:data2$Treatment)

# Stepwise select model with AIC
stepAIC(model)

# Stepwise selected model: 
final_model <- lm(data2$Success_Rate ~ data2$Treatment + data2$Clutch_Size + data2$Clutch_Size2)
summary(final_model)
anova(final_model)

### Plots ###
# Plot success rate vs clutch size
ggplot(data2, aes(x = Clutch_Size, y = Success_Rate, color = Treatment, group = Treatment)) +
  geom_point(size = 3) +
  geom_line() + 
  labs(x = "Clutch size", y = "Hatching success rate") +
  scale_color_manual(values = c("Deep" = "#0D0887", "Shallow" = "#FA9E3B")) +
  theme_classic() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))

# Plot success rate vs treatment
ggplot(data2, aes(x = Treatment, y = Success_Rate, fill = Treatment)) + 
  geom_boxplot()  +
  labs(x="Treatment", y="Hatching success rate") +
  scale_fill_manual(values=c("#FA9E3B","#0D0887")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_classic() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))

### Response to reviewer ###
##1. remove temperature
model=lm(data2$Success_Rate ~ data2$Treatment + data2$Clutch_Size + data2$Clutch_Size2 + 
           data2$Clutch_Size:data2$Treatment)
anova(model)
stepAIC(model)

model2=lm(formula = data2$Success_Rate ~ data2$Clutch_Size + data2$Treatment + 
            data2$Clutch_Size2)
anova(model2)

##2. remove clutch size
model=lm(data2$Success_Rate ~ data2$Treatment + mod$residuals + 
           mod$residuals:data2$Treatment)
anova(model)
stepAIC(model)


# Plot success rate vs sub-clutch size to show there is no difference with clutch size ####
ggplot(data2, aes(x = Split_Clutch_Size, y = Success_Rate, color = Treatment, group = Treatment)) +
  geom_point(size = 3) +  # Uses color aesthetic, increased point size for better visibility
  geom_line() +   # Lines also use the color aesthetic
  labs(x = "Sub-Clutch Size", y = "Success Rate") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15)) +
  scale_color_manual(values = c("Deep" = "#0D0887", "Shallow" = "#FA9E3B"))

########################################################################################################################################

####################################
###### Hatchling-level stats #######
####################################

# NB temperature is available for 18 of the 20 nests

attach(data)

#Remove hatchlings for which there is no temperatures
remove_na_values <- function(dataset, column_name) {
  # Remove NA values from the specified column
  cleaned_dataset <- dataset[!is.na(dataset[[column_name]]), ]
  
  # Return the cleaned dataset
  return(cleaned_dataset)
}

data1 <- remove_na_values(data, "Full_Incubation_MeanTemp")
print(data1)

attach(data1)

################################################
###### Hatchling size (mass and SCL) ###### 

###### (D) Linear mixed effects model to test for correlates of straight carapace length (SCL) ###### 

# take residuals of correlation between Temperature and Treatment. More important here than at clutch level because pseudo-points 
model=lm(data1$Full_Incubation_MeanTemp~data1$Treatment)
anova(model)

# Full starting model
model1 = lmer(data1$SCL_mm ~ data1$Treatment + data1$Clutch_size + model$residuals + 
                data1$Treatment*data1$Clutch_size + data1$Treatment*model$residuals +
              + (1| factor(data1$Maternal_ID)))
summary(model1)
anova(model1)

step(model1) ## simplify complex model

# Final, stepwise-selected model
model2=lmer(data1$SCL_mm ~ data1$Treatment + data1$Clutch_size + data1$Treatment:data1$Clutch_size+ (1 | factor(data1$Maternal_ID)))
summary(model2)
anova(model2) # Final model retains only a significant interaction of Treatment and clutch size on SCL. 

### Post hoc ###
# Visualise fitted lines
emmip(model2, Clutch_size ~ Treatment, cov.reduce = range)
# Summary of contrasts at max and min of range, by treatment
# cov.reduce required for continuous variable
emm <- emmeans(model2, pairwise ~ Clutch_size * Treatment, type = "response", cov.reduce = range)
pairs(emm, simple = "Treatment")

### Plots ###
hist (SCL_mm)

# SCL vs clutch size
ggplot(data1, aes(x = Clutch_size, y = SCL_mm, color= Treatment)) + 
  geom_point(size=2)  +
  geom_smooth(method = lm)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Clutch size", y="Hatchling straight carapace length (mm)") +
  theme(text = element_text(size = 15))+
  scale_color_manual(values = c("Deep" = "#0D0887", "Shallow" = "#FA9E3B"))

###### (E) Linear mixed effects model to test size-mass relationship ###### 

# Known allometric relationship between mass and SCL

# take residuals of correlation between Temperature and Treatment. More important here than at clutch level because pseudo-points 
model=lm(data1$Full_Incubation_MeanTemp~data1$Treatment)
anova(model)

# Full starting model
model1 = lmer(data1$SCL_mm ~ data1$Treatment + data1$Clutch_size + model$residuals + data1$Weight_g +
                data1$Treatment:data1$Clutch_size + data1$Treatment:model$residuals + data1$Treatment:data1$Weight_g +
                (1| factor(data1$Maternal_ID)))
step(model1)

# Final, stepwise-selected model
model2 = lmer(data1$SCL_mm ~ data1$Treatment + data1$Clutch_size + data1$Weight_g + (1 | factor(data1$Maternal_ID)) + 
                data1$Weight_g:data1$Treatment + 
                data1$Treatment:data1$Clutch_size)
              
summary(model2)
anova(model2)

### Post hoc ###
# Visualise fitted lines
emmip(model2, Clutch_size ~ Treatment, cov.reduce = range)
emmip(model2, Weight_g ~ Treatment, cov.reduce = range)
# Summary of contrasts at max and min of range, by treatment
# cov.reduce required for continuous variable
# Clutch size
emm <- emmeans(model2, pairwise ~ Clutch_size * Treatment, type = "response", cov.reduce = range)
pairs(emm, simple = "Treatment")
# Weight
emm <- emmeans(model2, pairwise ~ Weight_g * Treatment, type = "response", cov.reduce = range )
pairs(emm, simple = "Treatment")

### Plots ###
# Interaction between weight and treatment
ggplot(data, aes(x = Weight_g, y = SCL_mm, color= Treatment)) + 
  geom_point(size=1)  +
  geom_smooth(method = lm)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Mass (g)", y="Straight carapce length (mm)") +
  theme(text = element_text(size = 15))+
  scale_color_manual(values = c("Deep" = "#0D0887", "Shallow" = "#FA9E3B"))

### Response to reviewer ###

model=lm(data1$Full_Incubation_MeanTemp~data1$Treatment)
anova(model)

model2=lm(data1$Weight_g~data1$Treatment)
anova(model2)

model3=lmer(data1$SCL_mm ~ model2$residuals + model$residuals* data1$Treatment*data1$ Clutch_size+ (1| factor(data1$Maternal_ID)))
step(model3)

model4=lmer(data1$SCL_mm ~ model2$residuals + data1$Treatment + data1$Clutch_size + (1 | factor(data1$Maternal_ID)) + data1$Treatment:data1$Clutch_size)
anova(model4)

model5=lmer(data1$SCL_mm ~  data1$Treatment*data1$ Clutch_size+ (1| factor(data1$Maternal_ID)))
anova(model5)

model6=lmer(data1$SCL_mm ~ model2$residuals + model$residuals+ data1$Treatment+data1$ Clutch_size+ (1| factor(data1$Maternal_ID)))
anova(model6)

model7=lmer(data1$SCL_mm ~ model$residuals+ data1$Treatment+data1$ Clutch_size+ (1| factor(data1$Maternal_ID)))
anova(model7)

# Interaction between clutch size and temperature
library("plot3D")
## plotted for response to reviewer
scatter3D(model$residuals, data1$Clutch_size, data1$SCL_mm, clab = c("Size", "(mm)"), phi=0, bty="b2", xlab= "Mean Temperature (corrected by Treatment)", ylab="Nest size", zlab="Size in mm")
scatter3D(model$residuals, data1$Clutch_size, data1$SCL_mm, clab = c("SCL", "(mm)"), phi=0, bty="b2", xlab= "Mean Temperature", ylab="Clutch Size", zlab="SCL (mm)") ## plotted for response to reviewer

################################################
###### (F) Linear mixed effects model to test for correlates of mean run time ###### 

# Focusing on run time, taking the mean of the 2 runs
# First remove hatchlings that failed
data2 <- remove_na_values(data1, "Run_average")
print(data2)

# Take residuals of correlation between Temperature and Treatment. More important here than at clutch level because pseudo-points 
model=lm(data2$Full_Incubation_MeanTemp~data2$Treatment)
anova(model)

hist(data2$Run_average) # note Run_average needs to be log + 1 transformed
hist(log(data2$Run_average+1))

# Full starting model
# SCL used as a covariable
model1=lmer(log(data2$Run_average+1) ~ data2$Treatment + data2$Clutch_size + model$residuals + data2$SCL_mm + 
              data2$Treatment:data2$Clutch_size + data2$Treatment:model$residuals + data2$Treatment:data2$SCL_mm +
              (1| factor(data2$Maternal_ID)))

step(model1)

# Final, stepwise-selected model
model2=lmer(log(data2$Run_average + 1) ~ data2$Treatment + data2$Clutch_size + model$residuals + (1 | factor(data2$Maternal_ID)) + 
              data2$Clutch_size:data2$Treatment + data2$Treatment:model$residuals)

summary(model2)
anova(model2)

### Post hoc ###
# Visualise fitted lines
emmip(model2, Clutch_size ~ Treatment, cov.reduce = range)
emmip(model2, residuals ~ Treatment, cov.reduce = range)
# Summary of contrasts at max and min of range, by treatment
# cov.reduce required for continuous variable
# Clutch size
emm <- emmeans(model2, pairwise ~ Clutch_size * Treatment, cov.reduce = range)
pairs(emm, simple = "Treatment")
# Weight
emm <- emmeans(model2, pairwise ~ residuals * Treatment, cov.reduce = range )
pairs(emm, simple = "Treatment")


### Plots ###

#plotting run time by temperature split for treatment, corrected so they overlap
ggplot(data2, aes(x = (Full_Incubation_MeanTemp), y = log(Run_average+1), color = Treatment)) + 
  geom_point(size=2)  +
  geom_smooth(method = lm)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Mean incubation temperature (°C)", y="Average run time (s)") +
  scale_color_manual(values = c("Deep" = "#0D0887", "Shallow" = "#FA9E3B"))

#plotting run time by temperature split for treatment 
ggplot(data2, aes(x = Full_Incubation_MeanTemp, y = log(Run_average+1), color = Treatment)) + 
  geom_point(size=2)  +
  geom_smooth(method = lm)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Mean incubation temperatures (°C)", y="Average run time (s)") +
  scale_color_manual(values = c("Deep" = "#0D0887", "Shallow" = "#FA9E3B"))

### Response to reviewer ###
# exclude temperature
model1=lmer(log(data2$Run_average+1) ~ data2$SCL_mm+ data2$Clutch_size * data2$Treatment+ (1| factor(data2$Maternal_ID))) #### note here SCL is simply added as a covariate
step(model1)

model2=lmer(log(data2$Run_average + 1) ~ data2$Clutch_size + data2$Treatment + (1 | factor(data2$Maternal_ID)) + data2$Clutch_size:data2$Treatment)
anova(model2)

#run success vs failure
#two runs, so values can be 0 if both fails, 1 if one failed, 2 if all successful.
#Note, in this data it is only 0 or 2 so transforming into 0 for failed, and 1 for success. 
data$Run_success <- ifelse(data$Run_success == 2, 1, 0)

ggplot(data, aes(x = factor(Run_success), y = SCL_mm, color = Treatment)) + 
  geom_boxplot(size=2)  +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Run Success", y="SCL") 

model_Run_success <- glmer(data$Run_success ~ data$SCL_mm * data$Treatment+ (1| factor(data$Maternal_ID)), data = data, family = binomial) #### note here SCL is simply added as a covariate
Anova(model_Run_success, type=3)

pred_data <- with(data, expand.grid(
  Treatment = levels(factor(Treatment)),
  SCL_mm = seq(min(SCL_mm), max(SCL_mm), length.out = 408),
  Maternal_ID = levels(factor(Maternal_ID))[1]  # Use the first level of Family as a placeholder
))

# Predict probabilities
pred_data$probability <- predict(model_Run_success, newdata = pred_data, type = "response")

# Plotting the data
ggplot(pred_data, aes(x = SCL_mm, y = probability, colour = Treatment)) +
  geom_line() +
  labs(title = "Interaction Effect of Treatment and Hatchling Size on Success",
       x = "Straight Carapace Length (mm)",
       y = "Probability of Success",
       colour = "Treatment") +
  theme_minimal()


ggplot(data, aes(x = SCL_mm, color=factor(Run_success))) +
  geom_histogram(position = "identity", bins = 30, alpha = 0.1) +
  scale_fill_manual(values = c("Failure" = "#FF9999", "Success" = "#9999FF")) +  # Optional color coding
  labs(title = "Histogram of Hatchling Size by Treatment and Success",
       x = "Straight Carapace Length (mm)",
       y = "Frequency",
       fill = "Run Success") +
  theme_minimal()

data$Run_success <- factor(data$Run_success, levels = c(0, 2), labels = c("Failure", "Success"))

library(dplyr)
summary_data <- data %>%
  group_by(Treatment, Run_success) %>%
  summarise(Count = n(), .groups = 'drop')

# Plotting the data
ggplot(summary_data, aes(x = Treatment, y = Count, fill = Run_success)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of Successes and Failures by Treatment",
       x = "Treatment",
       y = "Count of Outcomes",
       fill = "Outcome") +
  theme_minimal()

#plotting run time by temperature split for treatment 
ggplot(data2, aes(x = Full_Incubation_MeanTemp, y = log(Run_average+1), color = Treatment)) + 
  geom_point(size=2)  +
  geom_smooth(method = lm)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Mean incubation temperatures (°C))", y="Average run time (s)") 

#plotting run time by clutch size for treatment 
ggplot(data2, aes(x = Clutch_size, y = log(Run_average+1), color = Treatment)) + 
  geom_point(size=2)  +
  geom_smooth(method = lm)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Clutch Size (eggs))", y="Average run time (log+1, in sec)") 


################################################
###### (G) Linear mixed effects model to test for correlates of mean self-righting time ###### 

data3 <- remove_na_values(data1, "Flip_average")
print(data3)

# Take residuals of correlation between Temperature and Treatment. More important here than at clutch level because pseudo-points 
model=lm(data3$Full_Incubation_MeanTemp~data3$Treatment)
anova(model)

hist(data2$Flip_average) # note Flip average needs to be log + 1 transformed
hist(log(data3$Flip_average+1))

# NB. Model with 3 way interactions: none get dropped by stepwise selection
# The model assumes many independent variables but they are correlated. Need simpler approach
# When almost everything is significant, this suggests there is something wrong and probably overparametrisation or non independence of variables
# For run time, it works out likely because it is a simpler outcome
# Therefore need to simplify - just include 2 way interactions with treatment, as done on clutch-level

# Full, starting model
# SCL used as a covariable
model1=lmer(log(data3$Flip_average+1) ~ data3$Treatment + model$residuals + data3$Clutch_size + data3$SCL_mm +
              data3$Treatment:model$residuals + data3$Treatment:data3$Clutch_size + data3$Treatment:data3$SCL_mm +
              (1| factor(data3$Maternal_ID)))

step(model1)

# Final, step-selected model
model2=lmer(log(data3$Flip_average + 1) ~ data3$Treatment + (1 | factor(data3$Maternal_ID)))
summary(model2)
anova(model2)

### Plots ###
ggplot(data3, aes(x = Treatment, y = log(Flip_average+1), fill = Treatment)) + 
  geom_boxplot()  +
  labs(x="Treatment", y="Time to self-righting (log+1) in seconds") +
  scale_fill_manual(values=c("#FA9E3B","#0D0887")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_classic() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))



