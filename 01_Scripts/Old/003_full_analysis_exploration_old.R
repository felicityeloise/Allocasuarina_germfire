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

# 1.1 Load custom functions ----
invisible(lapply(paste("./04_Functions/", dir ("04_Functions"), sep = ""), function(x) source (x)))


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


### For fecundity only
# TSF + height 
# Would this also include environmental variables?

# May include nested effects of location and transect


# Q3 :  How does contemporary fire history (i.e., fire frequency) and environmental attributes influence reproductive traits?
# Proportions of seedling, saplings, recruits and number of cones, seed size as response
# Fire frequency 
# Fire frequency * latitude
# Fire frequency * FPC



### For fecundity only
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




#~ seed size 
plot(tor_cum.prop$seed_wt_mg, tor_cum.prop$Proportion_germ, xlab = "Seed weight (mg)", ylab = "Proportion germination (%)", las = 1, main = expression(italic("Allocasuarina torulosa")))
plot(lit_cum.prop$seed_wt_mg, lit_cum.prop$Proportion_germ, xlab = "Seed weight (mg)", ylab = "Proportion germination (%)", las = 1, main = expression(italic("Allocasuarina littoralis")))
mtext(expression(bold("Proportion germination")), line = 2, at = -0.5)


plot(tor_cum.prop$seed_wt_mg, tor_cum.prop$t50, xlab = "Seed weight (mg)", ylab = "Time to 50%", las = 1)
plot(lit_cum.prop$seed_wt_mg, lit_cum.prop$t50, xlab = "Seed weight (mg)", ylab = "Time to 50%", las = 1)
mtext(expression(bold("Time to 50%")), line = 1.5, at = -0.5)





#~ fire frequency
boxplot(tor_cum.prop$Proportion_germ ~ tor_cum.prop$Fire_freq, xlab = "Fire frequency", ylab = "Proportion germination (%)", las = 1, main = expression(italic("Allocasuarina torulosa")))
boxplot(lit_cum.prop$Proportion_germ ~ lit_cum.prop$Fire_freq, xlab = "Fire frequency", ylab = "Proportion germination (%)", las = 1, main = expression(italic("Allocasuarina littoralis")))
mtext(expression(bold("Proportion germination")), line = 2, at = -0.5)


boxplot(tor_cum.prop$t50 ~ tor_cum.prop$Fire_freq, xlab = "Fire frequency", ylab = "Time to 50%", las = 1)
boxplot(lit_cum.prop$t50 ~ lit_cum.prop$Fire_freq, xlab = "Fire frequency", ylab = "Time to 50%", las = 1)
mtext(expression(bold("Time to 50%")), line = 1.5, at = -0.5)



# Fire frequency * seed weight
boxplot(tor_cum.prop$seed_wt_mg ~ tor_cum.prop$Fire_freq)
boxplot(lit_cum.prop$seed_wt_mg ~ lit_cum.prop$Fire_freq)

# Seed weight * latitude
boxplot(round(tor_cum.prop$seed_wt_mg, 4) ~ round(tor_cum.prop$Latitude, 2))
boxplot(round(lit_cum.prop$seed_wt_mg, 4) ~ round(lit_cum.prop$Latitude, 2))

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
boxplot(tor_cum.prop$Fire_freq ~ round(tor_cum.prop$seed_wt_mg, 4), xlab = "", ylab = "Fire frequency", las = 2, ylim = c(0,7))
mtext("Seed weight (mg)", line = 4, side = 1, cex = 0.7)

boxplot(lit_cum.prop$Fire_freq ~ round(lit_cum.prop$seed_wt_mg, 4), xlab = "", ylab = "Fire frequency", las = 2)
mtext(expression(bold("Seed size")), line = 2,at = -0.2)
mtext("Seed weight (mg)", line = 4, side = 1, cex = 0.7)


# FF * Latitude
boxplot(round(tor_transects$Latitude, 2) ~ tor_transects$Fire_freq)
boxplot(round(lit_transects$Latitude, 2) ~ lit_transects$Fire_freq)

# FF * FPC
boxplot(tor_transects$FPC ~ tor_transects$Fire_freq)
boxplot(lit_transects$FPC ~ lit_transects$Fire_freq)






# 4. Preliminary analyses to determine modelling method ----
# Scale variables and rerun glmer - scale(). Do we get non-significant effects similar to gam. Get predictions from glm
glmn <- glmer(Proportion_germ ~ 1 + (1|Individual) + (1|Set), family = binomial, data = tor_cum.prop)

mod1 <- glmer(Proportion_germ ~ seed_wt_mg * Fire_freq + (1|Individual) + (1|Set), family = binomial, data = tor_cum.prop)
summary(mod1)


# Rescale
r_seed_wt <- tor_cum.prop$seed_wt_mg
r_seed_wt <- scale(r_seed_wt)
tor_cum.prop$r_seed_wt <- scale(tor_cum.prop$seed_wt_mg)
tor_cum.prop$r_Fire_freq <- scale(tor_cum.prop$Fire_freq)
head(tor_cum.prop)
table(tor_cum.prop$r_Fire_freq, tor_cum.prop$Fire_freq)


mod2 <- glmer(Proportion_germ ~ r_seed_wt * r_Fire_freq + (1|Individual) + (1|Set), family = binomial, data = tor_cum.prop)
summary(mod2) # Rescaling removes the significant effect for fire frequency but retains all others. Only left with the binomial glm warning


# glmmTMB allegedly can fit linear mixed effect model with beta family but failed to run so we can only investigate beta regression for GAM
gbn <- gamm4(Proportion_germ ~ 1, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))

mod3 <- gamm4(Proportion_germ ~ ti(seed_wt_mg, Fire_freq) + ti(seed_wt_mg) + ti(Fire_freq), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(mod3$gam)
par(mfrow = c(2,2))
gam.check(mod3$gam) # Could change k' for seed weight and fire freq. Significant effect of interaction and seed weight.
plot.gam(mod3$gam) 

mod3.1 <- gamm4(Proportion_germ ~ ti(seed_wt_mg, Fire_freq) + ti(seed_wt_mg) + ti(Fire_freq, k = 7), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
par(mfrow = c(2,2))
gam.check(mod3.1$gam)
plot.gam(mod3.1$gam)

# mgcv::betar for GAM beta regression
gbb <- gam(Proportion_germ ~ 1, family = betar, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))

mod4 <- gam(Proportion_germ ~ ti(seed_wt_mg, Fire_freq) + ti(seed_wt_mg) + ti(Fire_freq), family = betar, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
par(mfrow = c(2,2))
gam.check(mod4) # edf is not 1 for seed weight and fire frequency as it is for the gamm4 with binomial family.
plot.gam(mod4) # As suggested by the check.gam results, we could improve model fit by adjusting k


mod4.1 <- gam(Proportion_germ ~ ti(seed_wt_mg, Fire_freq) + ti(seed_wt_mg, k = 10) + ti(Fire_freq, k = 4), family = betar, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
par(mfrow = c(2,2))
gam.check(mod4.1)
plot.gam(mod4.1)
# We can only improve the wiggliness for fire frequency by reducing k but this then means we cannot achieve the same smoothness for seed weight

mod4.2 <- gam(Proportion_germ ~ ti(seed_wt_mg, Fire_freq) + ti(seed_wt_mg, k = 6) + ti(Fire_freq), family = betar, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
par(mfrow = c(2,2))
gam.check(mod4.2)
plot.gam(mod4.2)


mod4.3 <- gam(Proportion_germ ~ ti(seed_wt_mg, Fire_freq) + ti(seed_wt_mg, k = 8) + ti(Fire_freq, k = 4), family = betar, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
par(mfrow = c(2,2))
gam.check(mod4.3)
plot.gam(mod4.3) # What if we keep seed weight as smooth as possible, without the extra wigglinness that is allowed if we specify k as 10. 



# Now we want to compare AICc for these models
glm_set <- list(mod1, mod2, glmn) # Better than null
aictab(glm_set) # No change in AIC but we should use the rescaled variables. 


gam_set <- list(mod3$mer, mod3.1$mer, gbn$mer)
aictab(gam_set) # Even though this results in confidence intervals crossing, the second model is considered a better model. Null model is worse but only marginally than mod3.1


# Need to calculate AICc manually, then create a table with all GAM models
gam.aic <- as.data.frame(AICc(mod4))
gam.aic$Model <- "mod4"
colnames(gam.aic) <- c("AICc", "Model")

gam.aic2 <- as.data.frame(AICc(mod4.1))
gam.aic2$Model <- "mod4.1"
colnames(gam.aic2) <- c("AICc", "Model")

gam.aic3 <- as.data.frame(AICc(mod4.2))
gam.aic3$Model <- "mod4.2"
colnames(gam.aic3) <- c("AICc", "Model")

gam.aic4 <- as.data.frame(AICc(mod4.3))
gam.aic4$Model <- "mod4.3"
colnames(gam.aic4) <- c("AICc", "Model")

gam.aci5 <- as.data.frame(AICc(gbb))
gam.aci5$Model <- "null"
colnames(gam.aci5) <- c("AICc", "Model")

gam.aic <- rbind(gam.aic, gam.aic2, gam.aic3, gam.aic4, gam.aci5)
gam.aic
gam.aic <- gam.aic[order(gam.aic$AICc),] 
gam.aic
# This suggests that the model which includes wiggliness for seed weight would be the best beta regression model. This does seem reasonable as seed weight was quite variable and we could see this wiggliness in the data. And that this is better than the null model. Something hasn't gone right with the later models.

# Compare the modelling methods by taking the best model from each set and plotting predictions from these models
mod2 <- glmer(Proportion_germ ~ r_seed_wt * r_Fire_freq + (1|Individual) + (1|Set), family = binomial, data = tor_cum.prop)
summary(mod2)

mod3.1 <- gamm4(Proportion_germ ~ ti(seed_wt_mg, Fire_freq) + ti(seed_wt_mg) + ti(Fire_freq, k = 7), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))

mod3.2<- gamm4(Proportion_germ ~ seed_wt_mg*Fire_freq, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))


mod4.1 <- gam(Proportion_germ ~ ti(seed_wt_mg, Fire_freq) + ti(seed_wt_mg, k = 10) + ti(Fire_freq, k = 4), family = betar, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))

mod4.2 <- gam(Proportion_germ ~ seed_wt_mg*Fire_freq, family = betar, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))

# Generate predictions
# As we have interactions in all models, we need to create new data in a different way than usual
# GLMER
# NOTE: For GLMER we need to use the predictSE() as predict() otherwise expects us to provide new data for the random effects which is incorrect

new_m2 <- expand.grid(
  r_Fire_freq_min = min(tor_cum.prop$r_Fire_freq),
  r_Fire_freq_max = max(tor_cum.prop$r_Fire_freq),
  r_Fire_freq_mean = mean(tor_cum.prop$r_Fire_freq),
  
  r_seed_wt = seq(min(tor_cum.prop$r_seed_wt), max(tor_cum.prop$r_seed_wt), length = 50)
)
m2_pred <- predictSE(mod2, newdata = new_m2, se.fit = T, type = 'response')
new_m2$fit <- m2_pred$fit
new_m2$se.fit <- m2_pred$se.fit
new_m2$lci <- new_m2$fit - (new_m2$se.fit * 1.96)
new_m2$uci <- new_m2$fit + (new_m2$se.fit * 1.96)
head(new_m2)

# Need to rescale fire frequency??
plot(new_m2$r_seed_wt, new_m2$fit,type="l")
plot(new_m2$r_Fire_freq, new_m2$fit)


# GAM binomial
new_m3 <- expand.grid(
  Fire_freq = seq(min(tor_cum.prop$Fire_freq), max(tor_cum.prop$Fire_freq), length = 50),
  seed_wt_mg = seq(min(tor_cum.prop$seed_wt_mg), max(tor_cum.prop$seed_wt_mg), length = 50)
)
m3_pred <- predict(mod3.1$gam, newdata = new_m3, se.fit = T, type = 'response')
new_m3$fit <- m3_pred$fit
new_m3$se.fit <- m3_pred$se.fit
new_m3$lci <- new_m3$fit - (new_m3$se.fit * 1.96)
new_m3$uci <- new_m3$fit + (new_m3$se.fit * 1.96)
head(new_m3)

plot(new_m3$seed_wt_mg, new_m3$fit)

plot(new_m3$Fire_freq, new_m3$fit)

# If we consider the binomial GLM and GAM, the log likelihood of the best model is the same but the GAM does not require us to rescale the variablesas it automatically takes into account that fire frequency and seed weight are on different scales. The prediction plots are also identical so we would recommend selecting the binomial GAM over the binomial GLM for ease of modelling.



# GAM beta regression
new_m4 <- expand.grid(
  Fire_freq = seq(min(tor_cum.prop$Fire_freq), max(tor_cum.prop$Fire_freq), length = 50),
  seed_wt_mg = seq(min(tor_cum.prop$seed_wt_mg), max(tor_cum.prop$seed_wt_mg), length = 50)
)
m4_pred <- predict(mod4.1, newdata = new_m4, se.fit = T, type = 'response')
new_m4$fit <- m4_pred$fit
new_m4$se.fit <- m4_pred$se.fit
new_m4$lci <- new_m4$fit - (new_m4$se.fit * 1.96)
new_m4$uci <- new_m4$fit + (new_m4$se.fit * 1.96)
head(new_m4)

plot(new_m4$seed_wt_mg, new_m4$fit)
plot(new_m4$Fire_freq, new_m4$fit)


# If we consider the binomial GAM and beta regression GAM the results are quite different, likely due to when we tried to improve model fit for the binomial GAM which meant that confidence intervals crossed at a certain point. Have a quick look at what happpens if we did not adjust k, which was a slightly worse model just to see how this would change the plotting


new_m3.1 <- expand.grid(
  Fire_freq = seq(min(tor_cum.prop$Fire_freq), max(tor_cum.prop$Fire_freq), length = 50),
  seed_wt_mg = seq(min(tor_cum.prop$seed_wt_mg), max(tor_cum.prop$seed_wt_mg), length = 50)
)
m3.1_pred <- predict(mod3$gam, newdata = new_m3.1, se.fit = T, type = 'response')
new_m3.1$fit <- m3.1_pred$fit
new_m3.1$se.fit <- m3.1_pred$se.fit
new_m3.1$lci <- new_m3.1$fit - (new_m3.1$se.fit * 1.96)
new_m3.1$uci <- new_m3.1$fit + (new_m3.1$se.fit * 1.96)

plot(new_m3.1$seed_wt_mg, new_m3.1$fit)

plot(new_m3.1$Fire_freq, new_m3.1$fit)


# This introduces some more variability that we likely expect to see. Take a look at the original plot of raw data.

plot(tor_cum.prop$seed_wt_mg, tor_cum.prop$Proportion_germ, xlab = "Seed weight (mg)", ylab = "Proportion germination (%)", las = 1, main = expression(italic("Allocasuarina torulosa")))

# Considering the raw data points for seed weight, the beta regression GAM is capturing the relationship of seed weight with proportion germination better as we see that for really low seed weights, proportion germination may be 0 which we see in the raw data. Beta regression GAMs are going to be the most appropriate for modelling proportion germination. Also better captures the variability of low or high germination rates at different fire frequencies. 


save.image('./02_Workspaces/003_full_analysis_exploration.RData')










# Notes on analyses:
# 1. We use a beta regression family as proportion germination is bounded by 0 and 1. 

# 2. We use Akaike's information criterion corrected for small sample sizes to rank models with the best model chosen as the model with the lowest AICc, which improves model fit over the null model by a change in AICc > 2. 

# 3. AIC does not work well with gam() from mgcv but beta regression is not available for gamm4 which is compatible with AIC, so we manually calculate AIC tables similar to that produced by AICtab.


# 4. Smooth terms in gamm4 only are 's' or 't2/te' or ' ti'. When fitting interactions we want to use tensor prduct smooths as this accounts for variables being on different scales (i.e., seed weight = grams, germination rates = % or days), thus, not requiring any variable rescaling. We use s smooth terms for main effects of treatment. We use 'ti' tensor product smooths for interaction terms, as this is suitable for specifying models with main effects and interactions, and t2 for main effects.


# 5. The marginal basis function (the basis for the smooth term) can equal a few things; e.g. fs (factor smooth), tp (default numeric), cr (cubic regression spline) or cc (cyclic regression spline). We use 'fs' for factor variables such as Treatment.

# 6. The order of predictor variables matter, we want to specify interaction first then continuous variables then factors.


# Question 1: How are germination rates influenced by seed treatment, seed attributes and/or fire frequency? ----
# This is in terms of proportion germination and time to 50% germination
# Treatment
# Seed weight
# Fire frequency
# Treatment * seed weight
# Treatment * fire frequency
# Seed weight * fire frequency 

pt_n <- gam(Proportion_germ ~ 1, family = betar, data = tor_cum.prop, random = ~ (1 | Individual) + (1 | Set))

pl_n <- gam(Proportion_germ ~ 1, family = betar, data = lit_cum.prop, random = ~ (1 | Individual) + (1 | Set))

pt_m1 <- gam(Proportion_germ ~ s(Treatment, bs = 'fs'), family = betar, data = tor_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(2,2))
gam.check(pt_m1)
plot.gam(pt_m1)

pl_m1 <- gam(Proportion_germ ~ s(Treatment, bs = 'fs'), family = betar, data = lit_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(2,2))
gam.check(pl_m1)
plot.gam(pl_m1)


pt_m2 <- gam(Proportion_germ ~ t2(seed_wt_mg), family = betar, data = tor_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(3,2))
gam.check(pt_m2)
plot.gam(pt_m2) # Changing k here does not improve model fit

pl_m2 <- gam(Proportion_germ ~ t2(seed_wt_mg), family = betar, data = lit_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(3,2))
gam.check(pl_m2)
plot.gam(pl_m2) # Changing k doesn't change fit of the GAM



pt_m3 <- gam(Proportion_germ ~ t2(Fire_freq), family = betar, data = tor_cum.prop, random = ~ (1 | Individual) + (1 | Set))
par(mfrow = c(3,2))
gam.check(pt_m3)
plot.gam(pt_m3)

pl_m3 <- gam(Proportion_germ ~ t2(Fire_freq, k = 3), family = betar, data = lit_cum.prop, random = ~ (1 | Individual) + (1 | Set))
par(mfrow = c(3,2))
gam.check(pl_m3)
plot.gam(pl_m3)



pt_m4 <- gam(Proportion_germ ~ ti(seed_wt_mg, Treatment, bs = c('tp', 'fs')) + t2(seed_wt_mg, k = 6) + s(Treatment, bs = 'fs') , random = ~ (1 | Individual) + (1 | Set), family = betar, data = tor_cum.prop)
par(mfrow = c(3,2))
gam.check(pt_m4)
plot.gam(pt_m4)

pl_m4 <- gam(Proportion_germ ~ ti(seed_wt_mg, Treatment, bs = c('tp', 'fs')) + t2(seed_wt_mg) + s(Treatment, bs = 'fs') , random = ~ (1 | Individual) + (1 | Set), family = betar, data = lit_cum.prop)
par(mfrow = c(3,2))
gam.check(pl_m4)
plot.gam(pl_m4)


pt_m5 <- gam(Proportion_germ ~ ti(Fire_freq, Treatment, bs = c('tp', 'fs')) + t2(Fire_freq, k = 6) + s(Treatment, bs = 'fs'), random = ~ (1 | Individual) + (1 | Set), family = betar, data = tor_cum.prop)
par(mfrow = c(3,2))
gam.check(pt_m5)
plot(pt_m5)

pl_m5 <- gam(Proportion_germ ~ ti(Fire_freq, Treatment, bs = c('tp', 'fs'), k = 3) + t2(Fire_freq, k = 3) + s(Treatment, bs = 'fs'), random = ~ (1 | Individual) + (1 | Set), family = betar, data = lit_cum.prop)
par(mfrow = c(3,2))
gam.check(pl_m5) # For this model we need to specify k for any predictor variables including fire frequency
plot(pl_m5) 



pt_m6 <- gam(formula = Proportion_germ ~ ti(seed_wt_mg, Fire_freq) +  t2(seed_wt_mg, k = 6) + t2(Fire_freq), family = betar, data = tor_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(2,2))
gam.check(pt_m6)
plot(pt_m6)

pl_m6 <- gam(formula = Proportion_germ ~ ti(seed_wt_mg, Fire_freq, k = 3) + t2(seed_wt_mg, k = 7) + t2(Fire_freq, k = 3), family = betar, data = lit_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(2,2))
gam.check(pl_m6)
plot(pl_m6)



pt_m7 <- gam(formula = Proportion_germ ~ t2(Latitude, k = 7), family = betar, data = tor_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(3,2))
gam.check(pt_m7)
plot(pt_m7)

pl_m7 <- gam(formula = Proportion_germ ~ t2(Latitude), family = betar, data = lit_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(3,2))
gam.check(pl_m7)
plot(pl_m7)



pt_m8 <- gam(formula = Proportion_germ ~ ti(Latitude, Treatment, bs = c('tp', 'fs')) + t2(Latitude, k = 6) + t2(Treatment, bs = 'fs'), family = betar, data = tor_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(3,2))
gam.check(pt_m8)
plot(pt_m8)

pl_m8 <- gam(formula = Proportion_germ ~ ti(Latitude, Treatment, bs = c('tp', 'fs')) + t2(Latitude) + t2(Treatment, bs = 'fs'), family = betar, data = lit_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(3,2))
gam.check(pl_m8)
plot(pl_m8)



pt_m9 <- gam(formula = Proportion_germ ~ ti(Fire_freq, Latitude, k = 3) + t2(Fire_freq, k = 3) + t2(Latitude, k = 8), family = betar, data = tor_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(3,2))
gam.check(pt_m9)
plot(pt_m9)


pl_m9 <- gam(formula = Proportion_germ ~ ti(Fire_freq, Latitude, k = 3) + t2(Fire_freq, k = 3) + t2(Latitude), family = betar, data = lit_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(3,3))
gam.check(pl_m9)
plot(pl_m9)



pt_m10 <- gam(formula = Proportion_germ ~ ti(seed_wt_mg, Latitude) + t2(seed_wt_mg, k = 6) + t2(Latitude, k = 6), family = betar, data = tor_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(3,3))
gam.check(pt_m10)
plot(pt_m10)


pl_m10 <- gam(formula = Proportion_germ ~ ti(seed_wt_mg, Latitude) + t2(seed_wt_mg, k = 6) + t2(Latitude), family = betar, data = lit_cum.prop, random = ~(1 | Individual) + (1 | Set))
par(mfrow = c(2,2))
gam.check(pl_m10)
plot(pl_m10)

##### For each model we get a warning that we have saturated likelihood which may be inaccurate, this is due to the data being one inflated (lots of 100% germination). Recommended to use bayesian regression models using 'Stan' as it allows for zero/one inflated beta


# Compare AICc
tp_aic <- as.data.frame(1:11)
tp_aic$AICc <- "NA"
tp_aic$Model <- "NA"
tp_aic$LL <- "NA"
tp_aic$AICc[1] <- AICc(pt_n)
tp_aic$Model[1] <- "Null"
tp_aic$LL[1] <- logLik(pt_n)
tp_aic$AICc[2] <- AICc(pt_m1)
tp_aic$Model[2] <- "m1"
tp_aic$LL[2] <- logLik(pt_m1)
tp_aic$AICc[3] <- AICc(pt_m2)
tp_aic$Model[3] <- "m2"
tp_aic$LL[3] <- logLik(pt_m3)
tp_aic$AICc[4] <- AICc(pt_m3)
tp_aic$Model[4] <- "m3"
tp_aic$LL[4] <- logLik(pt_m3)
tp_aic$AICc[5] <- AICc(pt_m4)
tp_aic$Model[5] <- "m4"
tp_aic$LL[5] <- logLik(pt_m4)
tp_aic$AICc[6] <- AICc(pt_m5)
tp_aic$Model[6] <- "m5"
tp_aic$LL[6] <- logLik(pt_m5)
tp_aic$AICc[7] <- AICc(pt_m6)
tp_aic$Model[7] <- "m6"
tp_aic$LL[7] <- logLik(pt_m6)
tp_aic$AICc[8] <- AICc(pt_m7)
tp_aic$Model[8] <- "m7"
tp_aic$LL[8] <- logLik(pt_m7)
tp_aic$AICc[9] <- AICc(pt_m8)
tp_aic$Model[9] <- "m8"
tp_aic$LL[9] <- logLik(pt_m8)
tp_aic$AICc[10] <- AICc(pt_m9)
tp_aic$Model[10] <- "m9"
tp_aic$LL[10] <- logLik(pt_m9)
tp_aic$AICc[11] <- AICc(pt_m10)
tp_aic$Model[11] <- "m10"
tp_aic$LL[11] <- logLik(pt_m10)
tp_aic <- tp_aic[, 2:ncol(tp_aic)]



lp_aic <- as.data.frame(1:11)
lp_aic$AICc <- "NA"
lp_aic$Model <- "NA"
lp_aic$LL <- "NA"
lp_aic$AICc[1] <- AICc(pl_n)
lp_aic$Model[1] <- "Null"
lp_aic$LL[1] <- logLik(pl_n)
lp_aic$AICc[2] <- AICc(pl_m1)
lp_aic$Model[2] <- "m1"
lp_aic$LL[2] <- logLik(pl_m1)
lp_aic$AICc[3] <- AICc(pl_m2)
lp_aic$Model[3] <- "m2"
lp_aic$LL[3] <- logLik(pl_m3)
lp_aic$AICc[4] <- AICc(pl_m3)
lp_aic$Model[4] <- "m3"
lp_aic$LL[4] <- logLik(pl_m3)
lp_aic$AICc[5] <- AICc(pl_m4)
lp_aic$Model[5] <- "m4"
lp_aic$LL[5] <- logLik(pl_m4)
lp_aic$AICc[6] <- AICc(pl_m5)
lp_aic$Model[6] <- "m5"
lp_aic$LL[6] <- logLik(pl_m5)
lp_aic$AICc[7] <- AICc(pl_m6)
lp_aic$Model[7] <- "m6"
lp_aic$LL[7] <- logLik(pl_m6)
lp_aic$AICc[8] <- AICc(pl_m7)
lp_aic$Model[8] <- "m7"
lp_aic$LL[8] <- logLik(pl_m7)
lp_aic$AICc[9] <- AICc(pl_m8)
lp_aic$Model[9] <- "m8"
lp_aic$LL[9] <- logLik(pl_m8)
lp_aic$AICc[10] <- AICc(pl_m9)
lp_aic$Model[10] <- "m9"
lp_aic$LL[10] <- logLik(pl_m9)
lp_aic$AICc[11] <- AICc(pl_m10)
lp_aic$Model[11] <- "m10"
lp_aic$LL[11] <- logLik(pl_m10)
lp_aic <- lp_aic[, 2:ncol(lp_aic)]


# Re-order and calculate Delta AICc
str(tp_aic)
tp_aic$AICc <- as.numeric(tp_aic$AICc)
tp_aic$LL <- as.numeric(tp_aic$LL)
str(tp_aic)
tp_aic <- tp_aic[order(tp_aic$AICc), ]
tp_aic # The best model is the model with interactive effect of fire frequency and seed weight.
tp_aic$Delta_AICc <- "0.00"
tp_aic$Delta_AICc[2] <- round(tp_aic$AICc[1]-tp_aic$AICc[2], 2)
tp_aic$Delta_AICc[3] <- round(tp_aic$AICc[1]-tp_aic$AICc[3], 2)
tp_aic$Delta_AICc[4] <- round(tp_aic$AICc[1]-tp_aic$AICc[4], 2)
tp_aic$Delta_AICc[5] <- round(tp_aic$AICc[1]-tp_aic$AICc[5], 2)
tp_aic$Delta_AICc[6] <- round(tp_aic$AICc[1]-tp_aic$AICc[6], 2)
tp_aic$Delta_AICc[7] <- round(tp_aic$AICc[1]-tp_aic$AICc[7], 2)
tp_aic$Delta_AICc[8] <- round(tp_aic$AICc[1]-tp_aic$AICc[8], 2)
tp_aic$Delta_AICc[9] <- round(tp_aic$AICc[1]-tp_aic$AICc[9], 2)
tp_aic$Delta_AICc[10] <- round(tp_aic$AICc[1]-tp_aic$AICc[10], 2)
tp_aic$Delta_AICc[11] <- round(tp_aic$AICc[1]-tp_aic$AICc[11],2)
tp_aic # All models are better than the null model, m6 with the interactive effect of seed weight and fire frequency is the best model.


str(lp_aic)
lp_aic$AICc <- as.numeric(lp_aic$AICc)
lp_aic$LL <- as.numeric(lp_aic$LL)
str(lp_aic)
lp_aic <- lp_aic[order(lp_aic$AICc), ]
lp_aic # The best model is the model with interactive effect of fire frequency and seed weight.
lp_aic$Delta_AICc <- "0.00"
lp_aic$Delta_AICc[2] <- round(lp_aic$AICc[1]-lp_aic$AICc[2], 2)
lp_aic$Delta_AICc[3] <- round(lp_aic$AICc[1]-lp_aic$AICc[3], 2)
lp_aic$Delta_AICc[4] <- round(lp_aic$AICc[1]-lp_aic$AICc[4], 2)
lp_aic$Delta_AICc[5] <- round(lp_aic$AICc[1]-lp_aic$AICc[5], 2)
lp_aic$Delta_AICc[6] <- round(lp_aic$AICc[1]-lp_aic$AICc[6], 2)
lp_aic$Delta_AICc[7] <- round(lp_aic$AICc[1]-lp_aic$AICc[7], 2)
lp_aic$Delta_AICc[8] <- round(lp_aic$AICc[1]-lp_aic$AICc[8], 2)
lp_aic$Delta_AICc[9] <- round(lp_aic$AICc[1]-lp_aic$AICc[9], 2)
lp_aic$Delta_AICc[10] <- round(lp_aic$AICc[1]-lp_aic$AICc[10], 2)
lp_aic$Delta_AICc[11] <- round(lp_aic$AICc[1]-lp_aic$AICc[11], 2)
lp_aic # Most models are better than the null model. Model 10 with latitude and seed weight is the best model with no model ranked within AICc <2.


# Take a better look at the best models
dev.off()
vis.gam(x = pl_m10,
        view = c("seed_wt_mg", "Latitude"),
        plot.type = 'contour')

par(mfrow = c(1,2))
vis.gam(x = pl_m10,
        view = c("seed_wt_mg", "Latitude"),
        plot.type = 'persp')
vis.gam(x = pl_m10,
        view = c("seed_wt_mg", "Latitude"),
        plot.type = 'persp',
        se = 2)

dev.off()
vis.gam(x = pt_m6,
        view = c("seed_wt_mg", "Fire_freq"),
        plot.type = 'contour')

par(mfrow = c(1,2))
vis.gam(x = pt_m6,
        view = c("seed_wt_mg", "Fire_freq"),
        plot.type = 'persp')
vis.gam(x = pt_m6,
        view = c("seed_wt_mg", "Fire_freq"),
        plot.type = 'persp',
        se = 2)
dev.off()

# Predict from the best model
tp_new <- expand.grid(
  seed_wt_mg = seq(min(tor_cum.prop$seed_wt_mg), max(tor_cum.prop$seed_wt_mg), length = 50),
  Fire_freq = seq(min(tor_cum.prop$Fire_freq), max(tor_cum.prop$Fire_freq), length = 50)
)
tp_pred <- predict(pt_m6, newdata = tp_new, se.fit = T, type = 'response')
tp_new$fit <- tp_pred$fit
tp_new$se.fit <- tp_pred$se.fit
tp_new$lci <- tp_new$fit - (tp_new$se.fit * 1.96)
tp_new$uci <- tp_new$fit + (tp_new$se.fit * 1.96)
head(tp_new)
tail(tp_new)
plot(tp_new$seed_wt_mg, tp_new$fit)
# We can treat fire as a categorical variable, colouring by fire frequency but not sure how to do this, probably for loop

plot(tp_new$seed_wt_mg, tp_new$fit, ylim = c(0,1), ylab = expression(bold("Proportion germination")), xlim = c(2.3, 7.1), xaxt = "n", las = 1, xlab = "", pch = 19, col = tp_new$Fire_freq)
axis(side = 1, at = seq(2.3, 7.1, 0.4), las = 1)
axis(side = 1, at = seq(2.3, 7.1, 0.1), las = 1, labels = F)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)

# Or should we be ignoring the fact that this is an interaction and producing new predictions for seed weight and fire frequency separately?

tp_sw <- data.frame(
  seed_wt_mg = seq(min(tor_cum.prop$seed_wt_mg), max(tor_cum.prop$seed_wt_mg), length = 50),
  Fire_freq = mean(tor_cum.prop$Fire_freq)
)
tp_swp <- predict(pt_m6, newdata = tp_sw, se.fit = T, type = 'response')
tp_sw$fit <- tp_swp$fit
tp_sw$se.fit <- tp_swp$se.fit
tp_sw$lci <- tp_sw$fit - (tp_sw$se.fit * 1.96)
tp_sw$uci <- tp_sw$fit + (tp_sw$se.fit * 1.96)
head(tp_sw)
plot(tp_sw$seed_wt_mg, tp_sw$fit)

par(mfrow = c(1,2))
plot(tp_sw$seed_wt_mg, tp_sw$fit, ylim = c(0,1), ylab = expression(bold("Poportion germination")), xlim = c(2.3, 7.1), xaxt = "n", las = 1, xlab = "", type = 'l')
axis(side = 1, at = seq(2.3, 7.1, 0.4), las = 1)
axis(side = 1, at = seq(2.3, 7.1, 0.1), labels = F)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "seed_wt_mg", data = "tp_sw", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")


tp_ff <- data.frame(
  seed_wt_mg = mean(tor_cum.prop$seed_weight),
  Fire_freq = seq(min(tor_cum.prop$Fire_freq), max(tor_cum.prop$Fire_freq), length = 50)
)
tp_ffp <- predict(pt_m6, newdata = tp_ff, se.fit = T, type = 'response')
tp_ff$fit <- tp_ffp$fit
tp_ff$se.fit <- tp_ffp$se.fit
tp_ff$lci <- tp_ff$fit - (tp_ff$se.fit * 1.96)
tp_ff$uci <- tp_ff$fit + (tp_ff$se.fit * 1.96)

plot(tp_ff$Fire_freq, tp_ff$fit, ylim = c(0,1), ylab = expression(bold("Proportion germination")), las = 1, xlab = expression(bold("Fire frequency")), type = 'l')
pg.ci(x = "Fire_freq", data = 'tp_ff', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")



lp_new <- expand.grid(
  seed_wt_mg = seq(min(lit_cum.prop$seed_wt_mg), max(lit_cum.prop$seed_wt_mg), length = 50),
  Latitude = seq(min(lit_cum.prop$Latitude), max(lit_cum.prop$Latitude), length = 50)
)
lp_pred <- predict(pl_m10, newdata = lp_new, se.fit = T, type = 'response')
lp_new$fit <- lp_pred$fit
lp_new$se.fit <- lp_pred$se.fit
lp_new$lci <- lp_new$fit - (lp_new$se.fit * 1.96)
lp_new$uci <- lp_new$fit + (lp_new$se.fit * 1.96)
head(lp_new)


plot(lp_new$Latitude, lp_new$fit)

lp_sw <- data.frame(
  seed_wt_mg = seq(min(lit_cum.prop$seed_wt_mg), max(lit_cum.prop$seed_wt_mg), length = 50),
  Latitude = mean(lit_cum.prop$Latitude)
)
lp_swp <- predict(pl_m10, newdata = lp_sw, se.fit = T, type = 'response')
lp_sw$fit <- lp_swp$fit
lp_sw$se.fit <- lp_swp$se.fit
lp_sw$lci <- lp_sw$fit - (lp_sw$se.fit * 1.96)
lp_sw$uci <- lp_sw$fit + (lp_sw$se.fit * 1.96)
head(lp_sw)
plot(lp_sw$seed_wt_mg, lp_sw$fit)

par(mfrow = c(1,2))
plot(lp_sw$seed_wt_mg, lp_sw$fit, ylim = c(0,1), ylab = expression(bold("Poportion germination")), xlim = c(1.3, 3.1), xaxt = "n", las = 1, xlab = "", type = 'l')
axis(side = 1, at = seq(1.3, 3.1, 0.2), las = 1)
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "seed_wt_mg", data = "lp_sw", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")


lp_l <- data.frame(
  seed_wt_mg = mean(lit_cum.prop$seed_wt_mg),
  Latitude = seq(min(lit_cum.prop$Latitude), max(lit_cum.prop$Latitude), length = 50)
)
lp_lp <- predict(pl_m10, newdata = lp_l, se.fit = T, type = 'response')
lp_l$fit <- lp_lp$fit
lp_l$se.fit <- lp_lp$se.fit
lp_l$lci <- lp_l$fit - (lp_l$se.fit * 1.96)
lp_l$uci <- lp_l$fit + (lp_l$se.fit * 1.96)
head(lp_l)

plot(lp_l$Latitude, lp_l$fit, ylim = c(0,1), ylab = expression(bold("Proportion germination")), las = 1, xlab = expression(bold("Latitude")), type = 'l')
pg.ci(x = "Latitude", data = 'lp_l', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")

save.image('./02_Workspaces/003_full_analysis_exploration.RData')



