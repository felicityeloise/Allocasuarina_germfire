# Written by Felicity Charles
# Caveat emptor
# Date: 12th March 2025

# Analysis of full experiment germination data and recruitment data


#1. Load packages ----
library(dplyr)
library(gamm4)
library(terra)
library(sf)
library(AICcmodavg)
library(mgcv)




# 2. Read in data ----
# For analyses, we want to analyse the species separately as we expect them to have different responses. 
# Germination data with seed characteristics and environmental data
dat_cum.prop <- read.csv('./00_Data/Full_experiment/Full_experiment_cumulative_germ.csv', header = T, stringsAsFactors = T)
head(dat_cum.prop); dim(dat_cum.prop)
dat_cum.prop$Treatment <- factor(dat_cum.prop$Treatment, levels = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
str(dat_cum.prop) # We also need to specify set as a factor variable
dat_cum.prop$Set <- factor(dat_cum.prop$Set, levels = c('1', '2', '3'))

tor_cum.prop <- dat_cum.prop[dat_cum.prop$Species == "torulosa", ]
lit_cum.prop <- dat_cum.prop[dat_cum.prop$Species == "littoralis", ]
head(tor_cum.prop); dim(tor_cum.prop)
head(lit_cum.prop); dim(lit_cum.prop)


dat_transects <- read.csv('./00_Data/Transect_data/Transect_level_enviro.csv', header = T, stringsAsFactors = T)
head(dat_transects); dim(dat_transects)
tor_transects <- dat_transects[dat_transects$Species == "torulosa", ]
lit_transects <- dat_transects[dat_transects$Species == "littoralis", ]
# Gatton transect was reduced to 50m in line with sampling for subsequent transects.

dat_tree <- read.csv('./00_Data/Transect_data/Tree_level_enviro.csv', header = T, stringsAsFactors = T)
head(dat_tree); dim(dat_tree)
tor_tree <- dat_tree[dat_tree$Species == "torulosa", ]
lit_tree <- dat_tree[dat_tree$Species == "littoralis", ]
# For GB1, this transect was not sampled as comprehensively as other transects as this was the first transect surveyed so much tree level data is missing. It has been limited to 50m, where it was originally measured out to 100m.

# 3. Data exploration ----
# Q1 : How are germination rates influenced by seed treatment, seed attributes and/or fire frequency?

# Proportion germination and time to 50% as response
# Treatment
# Seed weight
# Fire frequency
# Latitude
# Treatment * seed weight
# Treatment * fire frequency
# Seed weight * fire frequency 
# Treatment * latitude
# Fire frequency * latitude
# seed weight * latitude



# Q2 : How does recent fire activity influence population age structure and female fecundity?
# Proportions of seedling, saplings, recruits and number of cones as response
# TSF
# TSF * latitude
# TSF * FPC


### For fecundity only- but we have test whether there height and fire frequency are correlated
# Height as additive - https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2435.14619 
# TSF + height 
# Would this also include environmental variables?

# May include nested effects of location and transect


# Q3 :  How does contemporary fire history (i.e., fire frequency) and environmental attributes influence reproductive traits?
# Proportions of seedling, saplings, recruits and number of cones, seed size as response
# Fire frequency 
# Fire frequency * latitude
# Fire frequency * FPC



### For fecundity only - but we have test whether there height and fire frequency are correlated
# Height as additive - https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2435.14619 
# Fire frequency + height
# Would this also include environmental variables?

# May include nested effects of location and transect



# Question 1: Do seed treatment, seed attributes, or fire frequency influence germination rates? 
#### How do we incorporate environment here? Is this to do with seed weight * latitude interaction with decreasing latitude decreasing seed weight? Increasing latitude = increased fire frequency?

# The random effect of individual and set (i.e., replicate) are crossed as all individuals are in each set. 

# Proportion germination, time to 50%
par(mfrow = c(2,2))
#~ Treatment
boxplot(tor_cum.prop$Proportion_germ ~ tor_cum.prop$Treatment, xlab = "Treatment", ylab = "Proportion germination (%)", las = 1, main = expression(italic("Allocasuarina torulosa")))
boxplot(lit_cum.prop$Proportion_germ ~ lit_cum.prop$Treatment, xlab = "Treatment", ylab = "Proportion germination (%)", las = 1, main = expression(italic("Allocasuarina littoralis")))
mtext(expression(bold("Proportion germination")), line = 2, at = -0.5)


boxplot(tor_cum.prop$t50 ~ tor_cum.prop$Treatment, xlab = "Treatment", ylab = "Time to 50%", las = 1, ylim = c(0,25))
boxplot(lit_cum.prop$t50 ~ lit_cum.prop$Treatment, xlab = "Treatment", ylab = "Time to 50%", las = 1, ylim = c(0,25))
mtext(expression(bold("Time to 50%")), line = 1.5, at = -0.5)




#~ seed size - plot as xy
plot(tor_cum.prop$seed_weight, tor_cum.prop$Proportion_germ, xlab = "Seed weight (g)", ylab = "Proportion germination (%)", las = 1, main = expression(italic("Allocasuarina torulosa")))
plot(lit_cum.prop$seed_weight, lit_cum.prop$Proportion_germ, xlab = "Seed weight (g)", ylab = "Proportion germination (%)", las = 1, main = expression(italic("Allocasuarina littoralis")))
mtext(expression(bold("Proportion germination")), line = 2, at = -0.5)


plot(tor_cum.prop$seed_weight, tor_cum.prop$t50, xlab = "Seed weight (g)", ylab = "Time to 50%", las = 1)
plot(lit_cum.prop$seed_weight, lit_cum.prop$t50, xlab = "Seed weight (g)", ylab = "Time to 50%", las = 1)
mtext(expression(bold("Time to 50%")), line = 1.5, at = -0.5)





#~ fire frequency
boxplot(tor_cum.prop$Proportion_germ ~ tor_cum.prop$Fire_freq, xlab = "Fire frequency", ylab = "Proportion germination (%)", las = 1, main = expression(italic("Allocasuarina torulosa")))
boxplot(lit_cum.prop$Proportion_germ ~ lit_cum.prop$Fire_freq, xlab = "Fire frequency", ylab = "Proportion germination (%)", las = 1, main = expression(italic("Allocasuarina littoralis")))
mtext(expression(bold("Proportion germination")), line = 2, at = -0.5)


boxplot(tor_cum.prop$t50 ~ tor_cum.prop$Fire_freq, xlab = "Fire frequency", ylab = "Time to 50%", las = 1)
boxplot(lit_cum.prop$t50 ~ lit_cum.prop$Fire_freq, xlab = "Fire frequency", ylab = "Time to 50%", las = 1)
mtext(expression(bold("Time to 50%")), line = 1.5, at = -0.5)



# Fire frequency * seed weight
boxplot(tor_cum.prop$seed_weight ~ tor_cum.prop$Fire_freq)
boxplot(lit_cum.prop$seed_weight ~ lit_cum.prop$Fire_freq)

# Seed weight * latitude
boxplot(round(tor_cum.prop$seed_weight, 4) ~ round(tor_cum.prop$Latitude, 2))
boxplot(round(lit_cum.prop$seed_weight, 4) ~ round(lit_cum.prop$Latitude, 2))

# Fire frequency * latitude
boxplot(round(tor_cum.prop$Latitude) ~ tor_cum.prop$Fire_freq)
boxplot(round(lit_cum.prop$Latitude, 2) ~ lit_cum.prop$Fire_freq)


# ~ latitude
boxplot(tor_cum.prop$Proportion_germ ~ round(tor_cum.prop$Latitude, 2))
boxplot(lit_cum.prop$Proportion_germ ~ round(lit_cum.prop$Latitude, 2))

boxplot(tor_cum.prop$t50 ~ round(tor_cum.prop$Latitude, 2))
boxplot(lit_cum.prop$t50 ~ round(lit_cum.prop$Latitude, 2))









# Question 2: How does recent fire activity influence population age structure and female fecundity?
# Seedlings:adults ~ TSF
# Sapling:adults ~ TSF
# Seed and sap:adults ~ TSF
# Fecundity ~ TSF
par(mfrow = c(4,2))

boxplot(tor_transects$Proportion_seedlings ~ tor_transects$TSF, xlab = "Time since fire (years)", ylab = "Prop. seedlings", las = 1, main = expression(italic("Allocasuarina torulosa")))
boxplot(lit_transects$Proportion_seedlings ~ lit_transects$TSF, xlab = "Time since fire (years)", ylab = "Prop. seedlings", las = 1, main = expression(italic("Allocasuarina littoralis")))
mtext(expression(bold("Proportion seedlings")), line = 2,at = -0.2)

boxplot(tor_transects$Proportion_saplings ~ tor_transects$TSF, xlab = "Time since fire (years)", ylab = "Prop. saplings", las = 1)
boxplot(lit_transects$Proportion_saplings ~ lit_transects$TSF, xlab = "Time since fire (years)", ylab = "Prop. saplings", las = 1)
mtext(expression(bold("Proportion saplings")), line = 2,at = -0.2)


boxplot(tor_transects$Proportion_recruits ~ tor_transects$TSF, xlab = "Time since fire (years)", ylab = "Prop. recruits", las = 1)
boxplot(lit_transects$Proportion_recruits ~ lit_transects$TSF, xlab = "Time since fire (years)", ylab = "Prop. recruits", las = 1)
mtext(expression(bold("Proportion recruits")), line = 2,at = -0.2)


boxplot(tor_transects$Average_cones ~ tor_transects$TSF, xlab = "Time since fire (years)", ylab = "Female fecundity", las = 1)
boxplot(lit_transects$Average_cones ~ lit_transects$TSF, xlab = "Time since fire (years)", ylab = "Female fecundity", las = 1)
mtext(expression(bold("Female fecundity")), line = 2,at = -0.2)





# Question 3: How does contemporary fire history (i.e., fire frequency) and environmental attributes influence reproductive traits?
# Looking at proportions of seedlings, saplings and recruit, fecundity and seed size
# Seedlings:adults ~ FF
# Sapling:adults ~ FF
# Seed and sap:adults ~ FF
# Fecundity ~ FF
# Seedlings:adults ~ FF * Latitude * FPC
# Sapling:adults ~ FF * Latitude * FPC
# Seed and sap:adults ~ FF * Latitude * FPC
# Fecundity ~ FF * Latitude * FPC


par(mfrow = c(5,2))

boxplot(tor_transects$Proportion_seedlings ~ tor_transects$Fire_freq, xlab = "Fire frequency", ylab = "Prop. seedlings", las = 1, main = expression(italic("Allocasuarina torulosa")))
boxplot(lit_transects$Proportion_seedlings ~ lit_transects$Fire_freq, xlab = "Fire frequency", ylab = "Prop. seedlings", las = 1, main = expression(italic("Allocasuarina littoralis")))
mtext(expression(bold("Proportion seedlings")), line = 2,at = -0.2)

boxplot(tor_transects$Proportion_saplings ~ tor_transects$Fire_freq, xlab = "Fire frequency", ylab = "Prop. saplings", las = 1)
boxplot(lit_transects$Proportion_saplings ~ lit_transects$Fire_freq, xlab = "Fire frequency", ylab = "Prop. saplings", las = 1)
mtext(expression(bold("Proportion saplings")), line = 2,at = -0.2)


boxplot(tor_transects$Proportion_recruits ~ tor_transects$Fire_freq, xlab = "Fire frequency", ylab = "Prop. recruits", las = 1)
boxplot(lit_transects$Proportion_recruits ~ lit_transects$Fire_freq, xlab = "Fire frequency", ylab = "Prop. recruits", las = 1)
mtext(expression(bold("Proportion recruits")), line = 2,at = -0.2)


boxplot(tor_transects$Average_cones ~ tor_transects$Fire_freq, xlab = "Fire frequency", ylab = "Female fecundity", las = 1)
boxplot(lit_transects$Average_cones ~ lit_transects$Fire_freq, xlab = "Fire frequency", ylab = "Female fecundity", las = 1)
mtext(expression(bold("Female fecundity")), line = 2,at = -0.2)



# Seed size ~ FF
boxplot(tor_cum.prop$Fire_freq ~ round(tor_cum.prop$seed_weight, 4), xlab = "", ylab = "Fire frequency", las = 2, ylim = c(0,7))
mtext("Seed weight (g)", line = 4, side = 1, cex = 0.7)

boxplot(lit_cum.prop$Fire_freq ~ round(lit_cum.prop$seed_weight, 4), xlab = "", ylab = "Fire frequency", las = 2)
mtext(expression(bold("Seed size")), line = 2,at = -0.2)
mtext("Seed weight (g)", line = 4, side = 1, cex = 0.7)


# FF * Latitude
boxplot(round(tor_transects$Latitude, 2) ~ tor_transects$Fire_freq)
boxplot(round(lit_transects$Latitude, 2) ~ lit_transects$Fire_freq)

# FF * FPC
boxplot(tor_transects$FPC ~ tor_transects$Fire_freq)
boxplot(lit_transects$FPC ~ lit_transects$Fire_freq)











