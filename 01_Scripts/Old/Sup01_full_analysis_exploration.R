# Written by Felicity Charles
# Caveat emptor
# Date: 12th March 2025

# Analysis of full factorial germination experiment


#1. Load packages ----
library(dplyr)
library(gamm4)
library(terra)
library(sf)
library(AICcmodavg)
library(mgcv)
library(arm)

# 1.1 Load custom functions ----
invisible(lapply(paste("./04_Functions/", dir ("04_Functions"), sep = ""), function(x) source (x)))
load('./02_Workspaces/002_full_analysis_exploration.RData')


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
# TSF
# TSF + height 
# TSF * latitude
# TSF * FPC
# TSF * precipitation seasonality
# TSF * temperature seasonality

# May include nested effects of location and transect


# Q3 :  How does contemporary fire history (i.e., fire frequency) and environmental attributes influence reproductive traits?
# Proportions of seedling, saplings, recruits and number of cones, seed size as response
# Fire frequency 
# Fire frequency * latitude
# Fire frequency * FPC



### For fecundity only
# TSF
# TSF + height 
# TSF * latitude
# TSF * FPC
# TSF * precipitation seasonality
# TSF * temperature seasonality

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







# 4. QUESTION 1: How are germination rates influenced by seed treatment, seed attributes and/or fire frequency?
# 4.1 Proportion germination preliminary analyses -----
# 4.1.1 GLMER ----
# We need to rescale fire frequency and seed weight prior to modelling, but we need to wrap this in c() for its dimension-stripping properties to be able to predict from gamms with these rescaled variables.
tor_cum.prop$r_seed_wt <- c(scale(tor_cum.prop$seed_wt_mg))
lit_cum.prop$r_seed_wt <- c(scale(lit_cum.prop$seed_wt_mg))
tor_cum.prop$r_fire_freq <- c(scale(tor_cum.prop$Fire_freq))
lit_cum.prop$r_fire_freq <- c(scale(lit_cum.prop$Fire_freq))

mnull_t <- glmer(Proportion_germ ~ 1 + (1|Individual) + (1|Set), family = binomial, data = tor_cum.prop)
mnull_l <- glmer(Proportion_germ ~ 1 + (1|Individual) + (1|Set), family = binomial, data = lit_cum.prop)

m1_t <- glmer(Proportion_germ ~ Treatment + (1|Individual) + (1|Set), family = binomial, data = tor_cum.prop)
summary(m1_t)
m1_l <- glmer(Proportion_germ ~ Treatment + (1|Individual) + (1|Set), family = binomial, data = lit_cum.prop)
summary(m1_l)


m2_t <- glmer(Proportion_germ ~ r_seed_wt + (1|Individual) + (1|Set), family = binomial, data = tor_cum.prop)
summary(m2_t)
m2_l <- glmer(Proportion_germ ~ r_seed_wt + (1|Individual) + (1|Set), family = binomial, data = lit_cum.prop)
summary(m2_l)


m3_t <- glmer(Proportion_germ ~ r_fire_freq + (1|Individual) + (1|Set), family = binomial, data = tor_cum.prop)
summary(m3_t)
m3_l <- glmer(Proportion_germ ~ r_fire_freq + (1|Individual) + (1|Set), family = binomial, data = lit_cum.prop)
summary(m3_l)


m4_t <- glmer(Proportion_germ ~ Treatment * r_seed_wt + (1|Individual) + (1|Set), family = binomial, data = tor_cum.prop)
summary(m4_t)
m4_l <- glmer(Proportion_germ ~ Treatment * r_seed_wt + (1|Individual) + (1|Set), family = binomial, data = lit_cum.prop)
summary(m4_l)


m5_t <- glmer(Proportion_germ ~ Treatment * r_fire_freq + (1|Individual) + (1|Set), family = binomial, data = tor_cum.prop)
summary(m5_t)
m5_l <- glmer(Proportion_germ ~ Treatment * r_fire_freq + (1|Individual) + (1|Set), family = binomial, data = lit_cum.prop)
summary(m5_l)


m6_t <- glmer(Proportion_germ ~ r_seed_wt * r_fire_freq + (1|Individual) + (1|Set), family = binomial, data = tor_cum.prop)
summary(m6_t)
m6_l <- glmer(Proportion_germ ~ r_seed_wt * r_fire_freq + (1|Individual) + (1|Set), family = binomial, data = lit_cum.prop)
summary(m6_l)


# Compare model fit
tp_aic <- list(mnull_t, m1_t, m2_t, m3_t, m4_t, m5_t, m6_t)
lp_aic <- list(mnull_l, m1_l, m2_l, m3_l, m4_l, m5_l, m6_l)
aictab(tp_aic)
aictab(lp_aic)

# Predict from best GLMER models ----
# Littoralis
new_lp <- data.frame(r_seed_wt = seq(min(lit_cum.prop$r_seed_wt), max(lit_cum.prop$r_seed_wt), length = 50))
p_lp <- predictSE(m2_l, newdata = new_lp, se.fit = T, type = 'link')
new_lp$fit.link <- p_lp$fit
new_lp$se.link <- p_lp$se.fit
new_lp$lci.link <- new_lp$fit.link - (new_lp$se.link * 1.96)
new_lp$uci.link <- new_lp$fit.link + (new_lp$se.link * 1.96)

new_lp$fit <- invlogit(new_lp$fit.link)
new_lp$se <- invlogit(new_lp$se.link)
new_lp$lci <- invlogit(new_lp$lci.link)
new_lp$uci <- invlogit(new_lp$uci.link)



# Torulosa
# For predictions, we will choose 3 fire frequencies within the range of the data - min, mean, max
summary(tor_cum.prop$r_fire_freq) # Note the mean here is different to the mean provided by mean(tor_cum.prop$r_fire_freq)


new_tpl <- data.frame(r_seed_wt = seq(min(tor_cum.prop$r_seed_wt), max(tor_cum.prop$r_seed_wt), length = 50),
                      r_fire_freq = min(tor_cum.prop$r_fire_freq))
p_tpl <- predictSE(m6_t, newdata = new_tpl, se.fit = T, type = 'link')
new_tpl$fit.link <- p_tpl$fit
new_tpl$se.link <- p_tpl$se.fit
new_tpl$lci.link <- new_tpl$fit.link - (new_tpl$se.link * 1.96)
new_tpl$uci.link <- new_tpl$fit.link + (new_tpl$se.link * 1.96)

# Predict on the link scale (logit scale)
new_tpl$fit <- invlogit(new_tpl$fit.link)
new_tpl$se <- invlogit(new_tpl$se.link)
new_tpl$lci <- invlogit(new_tpl$lci.link)
new_tpl$uci <- invlogit(new_tpl$uci.link)

new_tpa <- data.frame(r_seed_wt = seq(min(tor_cum.prop$r_seed_wt), max(tor_cum.prop$r_seed_wt), length = 50),
                      r_fire_freq = 0)
p_tpa <- predictSE(m6_t, newdata = new_tpa, se.fit = T, type = 'link')
new_tpa$fit.link <- p_tpa$fit
new_tpa$se.link <- p_tpa$se.fit
new_tpa$lci.link <- new_tpa$fit.link - (new_tpa$se.link * 1.96)
new_tpa$uci.link <- new_tpa$fit.link + (new_tpa$se.link * 1.96)

# Predict on the link scale (logit scale)
new_tpa$fit <- invlogit(new_tpa$fit.link)
new_tpa$se <- invlogit(new_tpa$se.link)
new_tpa$lci <- invlogit(new_tpa$lci.link)
new_tpa$uci <- invlogit(new_tpa$uci.link)

new_tph <- data.frame(r_seed_wt = seq(min(tor_cum.prop$r_seed_wt), max(tor_cum.prop$r_seed_wt), length = 50),
                      r_fire_freq = max(tor_cum.prop$r_fire_freq))
p_tph <- predictSE(m6_t, newdata = new_tph, se.fit = T, type = 'link')
new_tph$fit.link <- p_tph$fit
new_tph$se.link <- p_tph$se.fit
new_tph$lci.link <- new_tph$fit.link - (new_tph$se.link * 1.96)
new_tph$uci.link <- new_tph$fit.link + (new_tph$se.link * 1.96)

# Predict on the link scale (logit scale)
new_tph$fit <- invlogit(new_tph$fit.link)
new_tph$se <- invlogit(new_tph$se.link)
new_tph$lci <- invlogit(new_tph$lci.link)
new_tph$uci <- invlogit(new_tph$uci.link)


# Plot predictions
# Black blue red
par(mfrow = c(1,1), mar = c(4,4,2,8), mgp = c(2.7,1,0))
plot(new_tpl$r_seed_wt, new_tpl$fit, ylim = c(0, 1), ylab = expression(bold("Proportion germination")), las = 1, type = 'l', xlab = "", xlim = c(-2.2, 2.6), col = "blue", xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = F)
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = seq(2.2, 7.2, 0.2))
axis(side = 2, at = seq(0, 1, 0.1,), las = 1)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "r_seed_wt", data = "new_tpl", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(new_tpa$r_seed_wt, new_tpa$fit)
pg.ci(x = "r_seed_wt", data = "new_tpa", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')

lines(new_tph$r_seed_wt, new_tph$fit, col = 'red')
pg.ci(x = 'r_seed_wt', data = 'new_tph', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')

par(xpd = NA)
legend(x = 3, y = 1, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", 'black', 'red'), title = expression(bold("Fire frequency")), lty = 1)
par(xpd = F)


dev.off()
plot(new_lp$r_seed_wt, new_lp$fit, ylim = c(0,1), ylab = expression(bold("Proportion germination")), las = 1, type = 'l', xlab = "", xaxt = "n", xlim = c(-2.8, 2.6), yaxt = "n")
axis(side = 1, at = seq(-2.7, 2.7, 0.3), labels = F)
axis(side = 1, at = seq(-2.7, 2.7, 0.3), labels = seq(1.3, 3.1, 0.1))
axis(side = 2, at = seq(0,1, 0.1), las = 1)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "r_seed_wt", data = "new_lp", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")



# Binomial GAMs no smoothing ----
# Test out binomial GAMs with no smooth functions - without smooth functions we need to use the rescaled variables.
# NOTE this is really just a GLM so we should just use GLMER instead

gnt <- gamm4(Proportion_germ ~ 1, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
gnl <- gamm4(Proportion_germ ~ 1, family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))

gt1 <- gamm4(Proportion_germ ~ Treatment, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt1$gam)
par(mfrow = c(2,2))
gam.check(gt1$gam)
gl1 <- gamm4(Proportion_germ ~ Treatment, family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl1$gam)
par(mfrow = c(2,2))
gam.check(gl1$gam)


gt2 <- gamm4(Proportion_germ ~ r_seed_wt, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt2$gam)
par(mfrow = c(2,2))
gam.check(gt2$gam)
gl2 <- gamm4(Proportion_germ ~ r_seed_wt, family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl2$gam)
par(mfrow = c(2,2))
gam.check(gl2$gam)

gt3 <- gamm4(Proportion_germ ~ r_fire_freq, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt3$gam)
par(mfrow = c(2,2))
gam.check(gt3$gam)
gl3 <- gamm4(Proportion_germ ~ r_fire_freq, family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl3$gam)
par(mfrow = c(2,2))
gam.check(gl3$gam)

gt4 <- gamm4(Proportion_germ ~ Treatment * r_seed_wt, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt4$gam)
par(mfrow = c(2,2))
gam.check(gt4$gam)
gl4 <- gamm4(Proportion_germ ~ Treatment * r_seed_wt, family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl4$gam)
par(mfrow = c(2,2))
gam.check(gl4$gam)


gt5 <- gamm4(Proportion_germ ~ Treatment * r_fire_freq, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt5$gam)
par(mfrow = c(2,2))
gam.check(gt5$gam)
gl5 <- gamm4(Proportion_germ ~ Treatment * r_fire_freq, family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl5$gam)
par(mfrow = c(2,2))
gam.check(gl5$gam)

gt6 <- gamm4(Proportion_germ ~ r_seed_wt * r_fire_freq, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt6$gam)
par(mfrow = c(2,2))
gam.check(gt6$gam)
gl6 <- gamm4(Proportion_germ ~ r_seed_wt * r_fire_freq, family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl6$gam)
par(mfrow = c(2,2))
gam.check(gl6$gam)



# Compare model fit
gt_aic <- list(gnt, gt1, gt2, gt3, gt4, gt5, gt6)
gl_aic <- list(gnl, gl1, gl2, gl3, gl4, gl5, gl6)
aictab(tp_aic)
aictab(lp_aic)

#  As with the GLMMERs, the same models are best. 


# Predict from best GAM no smoothing models ----
# Littoralis
newl <- data.frame(r_seed_wt = seq(min(lit_cum.prop$r_seed_wt), max(lit_cum.prop$r_seed_wt), length = 50))
pl <- predict(gl2$gam, newdata = newl, se.fit = T, type = 'link')
newl$fit.link <- pl$fit
newl$se.link <- pl$se.fit
newl$lci.link <- newl$fit.link - (newl$se.link * 1.96)
newl$uci.link <- newl$fit.link + (newl$se.link * 1.96)

newl$fit <- invlogit(newl$fit.link)
newl$se <- invlogit(newl$se.link)
newl$lci <- invlogit(newl$lci.link)
newl$uci <- invlogit(newl$uci.link)



# Torulosa
# For predictions, we will choose 3 fire frequencies within the range of the data - min, mean, max
summary(tor_cum.prop$r_fire_freq) # Note the mean here is different to the mean provided by mean(tor_cum.prop$r_fire_freq)


newtl <- data.frame(r_seed_wt = seq(min(tor_cum.prop$r_seed_wt), max(tor_cum.prop$r_seed_wt), length = 50),
                      r_fire_freq = min(tor_cum.prop$r_fire_freq))

ptl <- predict(gt6$gam, newdata = newtl, se.fit = T, type = 'link')
newtl$fit.link <- ptl$fit
newtl$se.link <- ptl$se.fit
newtl$lci.link <- newtl$fit.link - (newtl$se.link * 1.96)
newtl$uci.link <- newtl$fit.link + (newtl$se.link * 1.96)

# Predict on the link scale (logit scale)
newtl$fit <- invlogit(newtl$fit.link)
newtl$se <- invlogit(newtl$se.link)
newtl$lci <- invlogit(newtl$lci.link)
newtl$uci <- invlogit(newtl$uci.link)

newta <- data.frame(r_seed_wt = seq(min(tor_cum.prop$r_seed_wt), max(tor_cum.prop$r_seed_wt), length = 50),
                      r_fire_freq = 0)
pta <- predictSE(m6_t, newdata = newta, se.fit = T, type = 'link')
newta$fit.link <- pta$fit
newta$se.link <- pta$se.fit
newta$lci.link <- newta$fit.link - (newta$se.link * 1.96)
newta$uci.link <- newta$fit.link + (newta$se.link * 1.96)

# Predict on the link scale (logit scale)
newta$fit <- invlogit(newta$fit.link)
newta$se <- invlogit(newta$se.link)
newta$lci <- invlogit(newta$lci.link)
newta$uci <- invlogit(newta$uci.link)

newth <- data.frame(r_seed_wt = seq(min(tor_cum.prop$r_seed_wt), max(tor_cum.prop$r_seed_wt), length = 50),
                      r_fire_freq = max(tor_cum.prop$r_fire_freq))
pth <- predictSE(m6_t, newdata = newth, se.fit = T, type = 'link')
newth$fit.link <- pth$fit
newth$se.link <- pth$se.fit
newth$lci.link <- newth$fit.link - (newth$se.link * 1.96)
newth$uci.link <- newth$fit.link + (newth$se.link * 1.96)

# Predict on the link scale (logit scale)
newth$fit <- invlogit(newth$fit.link)
newth$se <- invlogit(newth$se.link)
newth$lci <- invlogit(newth$lci.link)
newth$uci <- invlogit(newth$uci.link)


# Plot predictions
# Black blue red
par(mfrow = c(1,1), mar = c(4,4,2,8), mgp = c(2.7,1,0))
plot(new_tpl$r_seed_wt, new_tpl$fit, ylim = c(0, 1), ylab = expression(bold("Proportion germination")), las = 1, type = 'l', xlab = "", xlim = c(-2.2, 2.6), col = "blue", xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = F)
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = seq(2.2, 7.1, 0.2))
axis(side = 2, at = seq(0, 1, 0.1,), las = 1)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "r_seed_wt", data = "new_tpl", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(new_tpa$r_seed_wt, new_tpa$fit)
pg.ci(x = "r_seed_wt", data = "new_tpa", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')

lines(new_tph$r_seed_wt, new_tph$fit, col = 'red')
pg.ci(x = 'r_seed_wt', data = 'new_tph', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')

par(xpd = NA)
legend(x = 3, y = 1, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", 'black', 'red'), title = expression(bold("Fire frequency")), lty = 1)
par(xpd = F)
# Identical to the GLMMER plot



dev.off()
plot(new_lp$r_seed_wt, new_lp$fit, ylim = c(0,1), ylab = expression(bold("Proportion germination")), las = 1, type = 'l', xlab = "", xaxt = "n", xlim = c(-2.8, 2.6), yaxt = "n")
axis(side = 1, at = seq(-2.7, 2.7, 0.3), labels = F)
axis(side = 1, at = seq(-2.7, 2.7, 0.3), labels = seq(1.3, 3.1, 0.1))
axis(side = 2, at = seq(0,1, 0.1), las = 1)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "r_seed_wt", data = "new_lp", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
# Also identical to the GLMER plot.



# 4.1.2 GAMs with default smoothing ----
# Test binomial gam with smooth terms, try only s() with rescale variables, then try te() or ti()
# Something to keep in mind with using s() and ti() terms is that the basis functions are different with s() using bs = 'tp' but ti() using bs = 'cr', and k is also different wit the default being k = 10 for s() and k = 5 for ti(). Best practice would be using a combination of univaritate main effect s() and interactive effect ti(). Simon wood has previously suggested that the functoinality allowing univariate tensor product smooths would be removed. Another suggestion was to remove the functionality of using te() or t2() when fitting models where main effects and interactive terms occur simultaneously, as ti() should be used in these instances.

# Best practice would dictate fitting models with main effects as s() and interactive effects as ti(). In these cases we can probably use the unscaled variables as ti() will model the interaction taking into account that the variables are on different scales. 
# So, we need to run three separate tests, one where we use defaults for these smooths, and two others where we adjust bs and k arguments to match for each smooth (i.e., bs = 'tp' and k = 10 or bs = 'cr' and k = 5). 
  # For the using the defaults for s(), we need to change k for fire frequency to the maximum for each dataset as 10 is too many knots for this variable. 

# See how this changes model prediction plots, then decide if there are any other types of smooths we would like to test. 


# We cannot fit gamm4 with ti() terms, I'm not sure whether t2() would work in an equivalent manner. Another issue I have encountered is that I cannot fit the main effect only model for littoralis fire frequency, even when adjusting k to 3 or less. 

gnt <- gam(Proportion_germ ~ 1, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
gnl <- gam(Proportion_germ ~ 1, family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))


# NOTE: This is just a linear model as there is no smoothing occurring. It may be more robust to model this as such. However, this would change the assumptions of what is the best model as this model has a lower AIC than model 6 which is the best GAM.
gt1_s <- gam(Proportion_germ ~ Treatment, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt1_s)
par(mfrow = c(2,2))
gam.check(gt1_s)
plot(gt1_s)
gl1_s <- gam(Proportion_germ ~ Treatment, family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl1_s)
par(mfrow = c(2,2))
gam.check(gl1_s)
plot(gl1_s)


gt2_s <- gam(Proportion_germ ~ s(seed_wt_mg), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt2_s)
par(mfrow = c(2,2))
gam.check(gt2_s)
plot(gt2_s)
gl2_s <- gam(Proportion_germ ~ s(seed_wt_mg), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl2_s)
par(mfrow = c(2,2))
gam.check(gl2_s)
plot(gl2_s)

gt3_s <- gam(Proportion_germ ~ s(Fire_freq, k = 7), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt3_s)
par(mfrow = c(2,2))
gam.check(gt3_s)
plot(gt3_s)
gl3_s <- gam(Proportion_germ ~ s(Fire_freq, k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set)) 
summary(gl3_s)
par(mfrow = c(2,2))
gam.check(gl3_s)
plot(gl3_s)

gt4_s <- gam(Proportion_germ ~ s(Treatment, bs = 're') + s(seed_wt_mg) + ti(seed_wt_mg, by = Treatment), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt4_s)
par(mfrow = c(2,2))
gam.check(gt4_s)
plot(gt4_s)
gl4_s <- gam(Proportion_germ ~ s(Treatment, bs ='re') + s(seed_wt_mg) + ti(seed_wt_mg, by = Treatment), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl4_s)
par(mfrow = c(2,2))
gam.check(gl4_s)
plot(gl4_s)


gt5_s <- gam(Proportion_germ ~ s(Treatment, bs ='re') + s(Fire_freq, k = 7) + ti(Fire_freq, by = Treatment), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt5_s)
par(mfrow = c(2,2))
gam.check(gt5_s)
plot(gt5_s)
gl5_s <- gam(Proportion_germ ~ s(Treatment, bs ='re') + s(Fire_freq, k = 3) + ti(Fire_freq, by = Treatment, k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl5_s)
par(mfrow = c(2,2))
gam.check(gl5_s)
plot(gl5_s)

gt6_s <- gam(Proportion_germ ~ s(seed_wt_mg) + s(Fire_freq,k = 7) + ti(seed_wt_mg, Fire_freq), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt6_s)
par(mfrow = c(2,2))
gam.check(gt6_s)
plot(gt6_s)
gl6_s <- gam(Proportion_germ ~ s(seed_wt_mg) + s(Fire_freq, k = 3) + ti(seed_wt_mg, Fire_freq, k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl6_s)
par(mfrow = c(2,2))
gam.check(gl6_s)
plot(gl6_s)



# Compare model fit
tg_aic <-  as.data.frame(1:7)
tg_aic$AICc <- "NA"
tg_aic$Model <- "NA"
tg_aic$LL <- "NA"
tg_aic$AICc[1] <- AICc(gnt)
tg_aic$Model[1] <- "Null"
tg_aic$LL[1] <- logLik(gnt)
tg_aic$AICc[2] <- AICc(gt1_s)
tg_aic$Model[2] <- "m1"
tg_aic$LL[2] <- logLik(gt1_s)
tg_aic$AICc[3] <- AICc(gt2_s)
tg_aic$Model[3] <- "m2"
tg_aic$LL[3] <- logLik(gt2_s)
tg_aic$AICc[4] <- AICc(gt3_s)
tg_aic$Model[4] <- "m3"
tg_aic$LL[4] <- logLik(gt3_s)
tg_aic$AICc[5] <- AICc(gt4_s)
tg_aic$Model[5] <- "m4"
tg_aic$LL[5] <- logLik(gt4_s)
tg_aic$AICc[6] <- AICc(gt5_s)
tg_aic$Model[6] <- "m5"
tg_aic$LL[6] <- logLik(gt5_s)
tg_aic$AICc[7] <- AICc(gt6_s)
tg_aic$Model[7] <- "m6"
tg_aic$LL[7] <- logLik(gt6_s)
tg_aic <- tg_aic[, 2:ncol(tg_aic)]


lg_aic <-  as.data.frame(1:7)
lg_aic$AICc <- "NA"
lg_aic$Model <- "NA"
lg_aic$LL <- "NA"
lg_aic$AICc[1] <- AICc(gnl)
lg_aic$Model[1] <- "Null"
lg_aic$LL[1] <- logLik(gnl)
lg_aic$AICc[2] <- AICc(gl1_s)
lg_aic$Model[2] <- "m1"
lg_aic$LL[2] <- logLik(gl1_s)
lg_aic$AICc[3] <- AICc(gl2_s)
lg_aic$Model[3] <- "m2"
lg_aic$LL[3] <- logLik(gl2_s)
lg_aic$AICc[4] <- AICc(gl3_s)
lg_aic$Model[4] <- "m3"
lg_aic$LL[4] <- logLik(gl3_s)
lg_aic$AICc[5] <- AICc(gl4_s)
lg_aic$Model[5] <- "m4"
lg_aic$LL[5] <- logLik(gl4_s)
lg_aic$AICc[6] <- AICc(gl5_s)
lg_aic$Model[6] <- "m5"
lg_aic$LL[6] <- logLik(gl5_s)
lg_aic$AICc[7] <- AICc(gl6_s)
lg_aic$Model[7] <- "m6"
lg_aic$LL[7] <- logLik(gl6_s)
lg_aic <- lg_aic[, 2:ncol(lg_aic)]


# Re-order and calculate Delta AICc
str(tg_aic)
tg_aic$AICc <- as.numeric(tg_aic$AICc)
tg_aic$LL <- as.numeric(tg_aic$LL)
str(tg_aic)
tg_aic <- tg_aic[order(tg_aic$AICc), ]
tg_aic # The best model is the model with interactive effect of fire frequency and seed weight.
tg_aic$Delta_AICc <- "0.00"
tg_aic$Delta_AICc[2] <- round(tg_aic$AICc[1]-tg_aic$AICc[2], 2)
tg_aic$Delta_AICc[3] <- round(tg_aic$AICc[1]-tg_aic$AICc[3], 2)
tg_aic$Delta_AICc[4] <- round(tg_aic$AICc[1]-tg_aic$AICc[4], 2)
tg_aic$Delta_AICc[5] <- round(tg_aic$AICc[1]-tg_aic$AICc[5], 2)
tg_aic$Delta_AICc[6] <- round(tg_aic$AICc[1]-tg_aic$AICc[6], 2)
tg_aic$Delta_AICc[7] <- round(tg_aic$AICc[1]-tg_aic$AICc[7], 2)
tg_aic
# As with other model 6 is best

str(lg_aic)
lg_aic$AICc <- as.numeric(lg_aic$AICc)
lg_aic$LL <- as.numeric(lg_aic$LL)
str(lg_aic)
lg_aic <- lg_aic[order(lg_aic$AICc), ]
lg_aic # The best model is the model with interactive effect of fire frequency and seed weight.
lg_aic$Delta_AICc <- "0.00"
lg_aic$Delta_AICc[2] <- round(lg_aic$AICc[1]-lg_aic$AICc[2], 2)
lg_aic$Delta_AICc[3] <- round(lg_aic$AICc[1]-lg_aic$AICc[3], 2)
lg_aic$Delta_AICc[4] <- round(lg_aic$AICc[1]-lg_aic$AICc[4], 2)
lg_aic$Delta_AICc[5] <- round(lg_aic$AICc[1]-lg_aic$AICc[5], 2)
lg_aic$Delta_AICc[6] <- round(lg_aic$AICc[1]-lg_aic$AICc[6], 2)
lg_aic$Delta_AICc[7] <- round(lg_aic$AICc[1]-lg_aic$AICc[7], 2)
lg_aic
# Model 2 is best with no model ranked within delta AIC < 2.



# 4.1.3 GAMs with s() default parameter smoothing ----
# Fit interacions with tp and k = 10, where possible. For fire frequency, set k to maximum closest to 10.

gt4_s_tp <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 10) + s(seed_wt_mg) + ti(seed_wt_mg, by = Treatment, bs = 'tp', k = 10), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt4_s_tp)
par(mfrow = c(2,2))
gam.check(gt4_s_tp)
plot(gt4_s_tp)
gl4_s_tp <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 10) + s(seed_wt_mg) + ti(seed_wt_mg, by = Treatment, bs = 'tp', k = 9), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl4_s_tp)
par(mfrow = c(2,2))
gam.check(gl4_s_tp)
plot(gl4_s_tp)


gt5_s_tp <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 10) + s(Fire_freq, k = 7) + ti(Fire_freq, by = Treatment, bs = 'tp', k = 7), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt5_s_tp)
par(mfrow = c(2,2))
gam.check(gt5_s_tp)
plot(gt5_s_tp)
gl5_s_tp <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 10) + s(Fire_freq, k = 3) + ti(Fire_freq, by = Treatment, bs = 'tp', k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl5_s_tp)
par(mfrow = c(2,2))
gam.check(gl5_s_tp)
plot(gl5_s_tp)

gt6_s_tp <- gam(Proportion_germ ~ s(seed_wt_mg) + s(Fire_freq,k = 7) + ti(seed_wt_mg, Fire_freq, bs = 'tp', k = 7), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt6_s_tp)
par(mfrow = c(2,2))
gam.check(gt6_s_tp)
plot(gt6_s_tp)
gl6_s_tp <- gam(Proportion_germ ~ s(seed_wt_mg) + s(Fire_freq, k = 3) + ti(seed_wt_mg, Fire_freq, bs = 'tp', k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl6_s_tp)
par(mfrow = c(2,2))
gam.check(gl6_s_tp)
plot(gl6_s_tp)



# Compare model fit
tg_s_aic <-  as.data.frame(1:7)
tg_s_aic$AICc <- "NA"
tg_s_aic$Model <- "NA"
tg_s_aic$LL <- "NA"
tg_s_aic$AICc[1] <- AICc(gnt)
tg_s_aic$Model[1] <- "Null"
tg_s_aic$LL[1] <- logLik(gnt)
tg_s_aic$AICc[2] <- AICc(gt1_s)
tg_s_aic$Model[2] <- "m1"
tg_s_aic$LL[2] <- logLik(gt1_s)
tg_s_aic$AICc[3] <- AICc(gt2_s)
tg_s_aic$Model[3] <- "m2"
tg_s_aic$LL[3] <- logLik(gt2_s)
tg_s_aic$AICc[4] <- AICc(gt3_s)
tg_s_aic$Model[4] <- "m3"
tg_s_aic$LL[4] <- logLik(gt3_s)
tg_s_aic$AICc[5] <- AICc(gt4_s_tp)
tg_s_aic$Model[5] <- "m4"
tg_s_aic$LL[5] <- logLik(gt4_s_tp)
tg_s_aic$AICc[6] <- AICc(gt5_s_tp)
tg_s_aic$Model[6] <- "m5"
tg_s_aic$LL[6] <- logLik(gt5_s_tp)
tg_s_aic$AICc[7] <- AICc(gt6_s_tp)
tg_s_aic$Model[7] <- "m6"
tg_s_aic$LL[7] <- logLik(gt6_s_tp)
tg_s_aic <- tg_s_aic[, 2:ncol(tg_s_aic)]


lg_s_aic <-  as.data.frame(1:7)
lg_s_aic$AICc <- "NA"
lg_s_aic$Model <- "NA"
lg_s_aic$LL <- "NA"
lg_s_aic$AICc[1] <- AICc(gnl)
lg_s_aic$Model[1] <- "Null"
lg_s_aic$LL[1] <- logLik(gnl)
lg_s_aic$AICc[2] <- AICc(gl1_s)
lg_s_aic$Model[2] <- "m1"
lg_s_aic$LL[2] <- logLik(gl1_s)
lg_s_aic$AICc[3] <- AICc(gl2_s)
lg_s_aic$Model[3] <- "m2"
lg_s_aic$LL[3] <- logLik(gl2_s)
lg_s_aic$AICc[4] <- AICc(gl3_s)
lg_s_aic$Model[4] <- "m3"
lg_s_aic$LL[4] <- logLik(gl3_s)
lg_s_aic$AICc[5] <- AICc(gl4_s_tp)
lg_s_aic$Model[5] <- "m4"
lg_s_aic$LL[5] <- logLik(gl4_s_tp)
lg_s_aic$AICc[6] <- AICc(gl5_s_tp)
lg_s_aic$Model[6] <- "m5"
lg_s_aic$LL[6] <- logLik(gl5_s_tp)
lg_s_aic$AICc[7] <- AICc(gl6_s_tp)
lg_s_aic$Model[7] <- "m6"
lg_s_aic$LL[7] <- logLik(gl6_s_tp)
lg_s_aic <- lg_s_aic[, 2:ncol(lg_s_aic)]


# Re-order and calculate Delta AICc
str(tg_s_aic)
tg_s_aic$AICc <- as.numeric(tg_s_aic$AICc)
tg_s_aic$LL <- as.numeric(tg_s_aic$LL)
str(tg_s_aic)
tg_s_aic <- tg_s_aic[order(tg_s_aic$AICc), ]
tg_s_aic # The best model is the model with interactive effect of fire frequency and seed weight.
tg_s_aic$Delta_AICc <- "0.00"
tg_s_aic$Delta_AICc[2] <- round(tg_s_aic$AICc[1]-tg_s_aic$AICc[2], 2)
tg_s_aic$Delta_AICc[3] <- round(tg_s_aic$AICc[1]-tg_s_aic$AICc[3], 2)
tg_s_aic$Delta_AICc[4] <- round(tg_s_aic$AICc[1]-tg_s_aic$AICc[4], 2)
tg_s_aic$Delta_AICc[5] <- round(tg_s_aic$AICc[1]-tg_s_aic$AICc[5], 2)
tg_s_aic$Delta_AICc[6] <- round(tg_s_aic$AICc[1]-tg_s_aic$AICc[6], 2)
tg_s_aic$Delta_AICc[7] <- round(tg_s_aic$AICc[1]-tg_s_aic$AICc[7], 2)
tg_s_aic
# As with other model 6 is best

str(lg_s_aic)
lg_s_aic$AICc <- as.numeric(lg_s_aic$AICc)
lg_s_aic$LL <- as.numeric(lg_s_aic$LL)
str(lg_s_aic)
lg_s_aic <- lg_s_aic[order(lg_s_aic$AICc), ]
lg_s_aic # The best model is the model with interactive effect of fire frequency and seed weight.
lg_s_aic$Delta_AICc <- "0.00"
lg_s_aic$Delta_AICc[2] <- round(lg_s_aic$AICc[1]-lg_s_aic$AICc[2], 2)
lg_s_aic$Delta_AICc[3] <- round(lg_s_aic$AICc[1]-lg_s_aic$AICc[3], 2)
lg_s_aic$Delta_AICc[4] <- round(lg_s_aic$AICc[1]-lg_s_aic$AICc[4], 2)
lg_s_aic$Delta_AICc[5] <- round(lg_s_aic$AICc[1]-lg_s_aic$AICc[5], 2)
lg_s_aic$Delta_AICc[6] <- round(lg_s_aic$AICc[1]-lg_s_aic$AICc[6], 2)
lg_s_aic$Delta_AICc[7] <- round(lg_s_aic$AICc[1]-lg_s_aic$AICc[7], 2)
lg_s_aic
# Model 2



# 4.1.4 GAMs with ti() default parameter smoothing ----
# Fit interacions with cr and k = 5, where possible.

gt2_s_cr <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'cr', k = 5), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt2_s_cr)
par(mfrow = c(2,2))
gam.check(gt2_s_cr)
plot(gt2_s_cr)
gl2_s_cr <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'cr', k = 5), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl2_s_cr)
par(mfrow = c(2,2))
gam.check(gl2_s_cr)
plot(gl2_s_cr)

gt3_s_cr <- gam(Proportion_germ ~ s(Fire_freq, bs = 'cr', k = 5), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt3_s_cr)
par(mfrow = c(2,2))
gam.check(gt3_s_cr)
plot(gt3_s_cr)
gl3_s_cr <- gam(Proportion_germ ~ s(Fire_freq, bs = 'cr', k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set)) 
summary(gl3_s_cr)
par(mfrow = c(2,2))
gam.check(gl3_s_cr)
plot(gl3_s_cr)

gt4_s_cr <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 5) + s(seed_wt_mg, bs = 'cr', k = 5) + ti(seed_wt_mg, by = Treatment, bs = 'cr', k = 5), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt4_s_cr)
par(mfrow = c(2,2))
gam.check(gt4_s_cr)
plot(gt4_s_cr)
gl4_s_cr <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 5) + s(seed_wt_mg, bs = 'cr', k = 5) + ti(seed_wt_mg, by = Treatment, bs = 'cr', k = 5), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl4_s_cr)
par(mfrow = c(2,2))
gam.check(gl4_s_cr)
plot(gl4_s_cr)


gt5_s_cr <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 5) + s(Fire_freq, bs = 'cr', k = 5) + ti(Fire_freq, by = Treatment, bs = 'cr', k = 5), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt5_s_cr)
par(mfrow = c(2,2))
gam.check(gt5_s_cr)
plot(gt5_s_cr)
gl5_s_cr <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 5) + s(Fire_freq, bs = 'cr', k = 3) + ti(Fire_freq, by = Treatment, bs = 'cr', k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl5_s_cr)
par(mfrow = c(2,2))
gam.check(gl5_s_cr)
plot(gl5_s_cr)

gt6_s_cr <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'cr', k = 5) + s(Fire_freq, bs = 'cr', k = 5) + ti(seed_wt_mg, Fire_freq, bs = 'cr', k = 5), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gt6_s_cr)
par(mfrow = c(2,2))
gam.check(gt6_s_cr)
plot(gt6_s_cr)
gl6_s_cr <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'cr', k = 5) + s(Fire_freq, bs = 'cr', k = 3) + ti(seed_wt_mg, Fire_freq, bs = 'cr', k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set))
summary(gl6_s_cr)
par(mfrow = c(2,2))
gam.check(gl6_s_cr)
plot(gl6_s_cr)


# Compare model fit
tg_cr_aic <-  as.data.frame(1:7)
tg_cr_aic$AICc <- "NA"
tg_cr_aic$Model <- "NA"
tg_cr_aic$LL <- "NA"
tg_cr_aic$AICc[1] <- AICc(gnt)
tg_cr_aic$Model[1] <- "Null"
tg_cr_aic$LL[1] <- logLik(gnt)
tg_cr_aic$AICc[2] <- AICc(gt1_s)
tg_cr_aic$Model[2] <- "m1"
tg_cr_aic$LL[2] <- logLik(gt1_s)
tg_cr_aic$AICc[3] <- AICc(gt2_s_cr)
tg_cr_aic$Model[3] <- "m2"
tg_cr_aic$LL[3] <- logLik(gt2_s_cr)
tg_cr_aic$AICc[4] <- AICc(gt3_s_cr)
tg_cr_aic$Model[4] <- "m3"
tg_cr_aic$LL[4] <- logLik(gt3_s_cr)
tg_cr_aic$AICc[5] <- AICc(gt4_s_cr)
tg_cr_aic$Model[5] <- "m4"
tg_cr_aic$LL[5] <- logLik(gt4_s_cr)
tg_cr_aic$AICc[6] <- AICc(gt5_s_cr)
tg_cr_aic$Model[6] <- "m5"
tg_cr_aic$LL[6] <- logLik(gt5_s_cr)
tg_cr_aic$AICc[7] <- AICc(gt6_s_cr)
tg_cr_aic$Model[7] <- "m6"
tg_cr_aic$LL[7] <- logLik(gt6_s_cr)
tg_cr_aic <- tg_cr_aic[, 2:ncol(tg_cr_aic)]


lg_cr_aic <-  as.data.frame(1:7)
lg_cr_aic$AICc <- "NA"
lg_cr_aic$Model <- "NA"
lg_cr_aic$LL <- "NA"
lg_cr_aic$AICc[1] <- AICc(gnl)
lg_cr_aic$Model[1] <- "Null"
lg_cr_aic$LL[1] <- logLik(gnl)
lg_cr_aic$AICc[2] <- AICc(gl1_s)
lg_cr_aic$Model[2] <- "m1"
lg_cr_aic$LL[2] <- logLik(gl1_s)
lg_cr_aic$AICc[3] <- AICc(gl2_s_cr)
lg_cr_aic$Model[3] <- "m2"
lg_cr_aic$LL[3] <- logLik(gl2_s_cr)
lg_cr_aic$AICc[4] <- AICc(gl3_s_cr)
lg_cr_aic$Model[4] <- "m3"
lg_cr_aic$LL[4] <- logLik(gl3_s_cr)
lg_cr_aic$AICc[5] <- AICc(gl4_s_cr)
lg_cr_aic$Model[5] <- "m4"
lg_cr_aic$LL[5] <- logLik(gl4_s_cr)
lg_cr_aic$AICc[6] <- AICc(gl5_s_cr)
lg_cr_aic$Model[6] <- "m5"
lg_cr_aic$LL[6] <- logLik(gl5_s_cr)
lg_cr_aic$AICc[7] <- AICc(gl6_s_cr)
lg_cr_aic$Model[7] <- "m6"
lg_cr_aic$LL[7] <- logLik(gl6_s_cr)
lg_cr_aic <- lg_cr_aic[, 2:ncol(lg_cr_aic)]


# Re-order and calculate Delta AICc
str(tg_cr_aic)
tg_cr_aic$AICc <- as.numeric(tg_cr_aic$AICc)
tg_cr_aic$LL <- as.numeric(tg_cr_aic$LL)
str(tg_cr_aic)
tg_cr_aic <- tg_cr_aic[order(tg_cr_aic$AICc), ]
tg_cr_aic # The best model is the model with interactive effect of fire frequency and seed weight.
tg_cr_aic$Delta_AICc <- "0.00"
tg_cr_aic$Delta_AICc[2] <- round(tg_cr_aic$AICc[1]-tg_cr_aic$AICc[2], 2)
tg_cr_aic$Delta_AICc[3] <- round(tg_cr_aic$AICc[1]-tg_cr_aic$AICc[3], 2)
tg_cr_aic$Delta_AICc[4] <- round(tg_cr_aic$AICc[1]-tg_cr_aic$AICc[4], 2)
tg_cr_aic$Delta_AICc[5] <- round(tg_cr_aic$AICc[1]-tg_cr_aic$AICc[5], 2)
tg_cr_aic$Delta_AICc[6] <- round(tg_cr_aic$AICc[1]-tg_cr_aic$AICc[6], 2)
tg_cr_aic$Delta_AICc[7] <- round(tg_cr_aic$AICc[1]-tg_cr_aic$AICc[7], 2)
tg_cr_aic
# As with other model 6 is best

str(lg_cr_aic)
lg_cr_aic$AICc <- as.numeric(lg_cr_aic$AICc)
lg_cr_aic$LL <- as.numeric(lg_cr_aic$LL)
str(lg_cr_aic)
lg_cr_aic <- lg_cr_aic[order(lg_cr_aic$AICc), ]
lg_cr_aic # The best model is the model with interactive effect of fire frequency and seed weight.
lg_cr_aic$Delta_AICc <- "0.00"
lg_cr_aic$Delta_AICc[2] <- round(lg_cr_aic$AICc[1]-lg_cr_aic$AICc[2], 2)
lg_cr_aic$Delta_AICc[3] <- round(lg_cr_aic$AICc[1]-lg_cr_aic$AICc[3], 2)
lg_cr_aic$Delta_AICc[4] <- round(lg_cr_aic$AICc[1]-lg_cr_aic$AICc[4], 2)
lg_cr_aic$Delta_AICc[5] <- round(lg_cr_aic$AICc[1]-lg_cr_aic$AICc[5], 2)
lg_cr_aic$Delta_AICc[6] <- round(lg_cr_aic$AICc[1]-lg_cr_aic$AICc[6], 2)
lg_cr_aic$Delta_AICc[7] <- round(lg_cr_aic$AICc[1]-lg_cr_aic$AICc[7], 2)
lg_cr_aic
# Model 2 is best



# 4.1.5  GAMs with s() default parameter smoothing and method = 'ML' -----
# Test fitting GAMs with tp defaults and method = ML, if ML is as wiggly as the default UBRE (GCV), then the GAM is not overfitting and we don't need to specify method


gnt_m <- gam(Proportion_germ ~ 1, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
gnl_m <- gam(Proportion_germ ~ 1, family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')

gt1_m <- gam(Proportion_germ ~ Treatment, family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt1_m)
par(mfrow = c(2,2))
gam.check(gt1_m)
gl1_m <- gam(Proportion_germ ~ Treatment, family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl1_m)
par(mfrow = c(2,2))
gam.check(gl1_m)


gt2_m <- gam(Proportion_germ ~ s(seed_wt_mg), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt2_m)
par(mfrow = c(2,2))
gam.check(gt2_m)
gl2_m <- gam(Proportion_germ ~ s(seed_wt_mg), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl2_m)
par(mfrow = c(2,2))
gam.check(gl2_m)

gt3_m <- gam(Proportion_germ ~ s(Fire_freq, k = 7), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt3_m)
par(mfrow = c(2,2))
gam.check(gt3_m)
gl3_m <- gam(Proportion_germ ~ s(Fire_freq, k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl3_m)
par(mfrow = c(2,2))
gam.check(gl3_m)

gt4_m <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 10) + s(seed_wt_mg) + ti(seed_wt_mg, by = Treatment, bs = 'tp', k = 10), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt4_m)
par(mfrow = c(2,2))
gam.check(gt4_m)
plot(gt4_m)
gl4_m <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 10) + s(seed_wt_mg) + ti(seed_wt_mg, by = Treatment, bs = 'tp', k = 6), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl4_m)
par(mfrow = c(2,2))
gam.check(gl4_m)
plot(gl4_m)


gt5_m <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 10) + s(Fire_freq, k = 7) + ti(Fire_freq, by = Treatment, bs = 'tp', k = 10), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt5_m)
par(mfrow = c(2,2))
gam.check(gt5_m)
plot(gt5_m)
gl5_m <- gam(Proportion_germ ~ s(Treatment, bs ='re', k = 10) + s(Fire_freq, k = 3) + ti(Fire_freq, by = Treatment, bs = 'tp', k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl5_m)
par(mfrow = c(2,2))
gam.check(gl5_m)
plot(gl5_m)

gt6_m <- gam(Proportion_germ ~ s(seed_wt_mg) + s(Fire_freq,k = 7) + ti(seed_wt_mg, Fire_freq, bs = 'tp', k = 7), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt6_m)
par(mfrow = c(2,2))
gam.check(gt6_m)
plot(gt6_m)
gl6_m <- gam(Proportion_germ ~ s(seed_wt_mg) + s(Fire_freq, k = 3) + ti(seed_wt_mg, Fire_freq, bs = 'tp', k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl6_m)
par(mfrow = c(2,2))
gam.check(gl6_m)
plot(gl6_m)



# Compare model fit
tg_m_aic <-  as.data.frame(1:7)
tg_m_aic$AICc <- "NA"
tg_m_aic$Model <- "NA"
tg_m_aic$LL <- "NA"
tg_m_aic$AICc[1] <- AICc(gnt_m)
tg_m_aic$Model[1] <- "Null"
tg_m_aic$LL[1] <- logLik(gnt_m)
tg_m_aic$AICc[2] <- AICc(gt1_m)
tg_m_aic$Model[2] <- "m1"
tg_m_aic$LL[2] <- logLik(gt1_m)
tg_m_aic$AICc[3] <- AICc(gt2_m)
tg_m_aic$Model[3] <- "m2"
tg_m_aic$LL[3] <- logLik(gt2_m)
tg_m_aic$AICc[4] <- AICc(gt3_m)
tg_m_aic$Model[4] <- "m3"
tg_m_aic$LL[4] <- logLik(gt3_m)
tg_m_aic$AICc[5] <- AICc(gt4_m)
tg_m_aic$Model[5] <- "m4"
tg_m_aic$LL[5] <- logLik(gt4_m)
tg_m_aic$AICc[6] <- AICc(gt5_m)
tg_m_aic$Model[6] <- "m5"
tg_m_aic$LL[6] <- logLik(gt5_m)
tg_m_aic$AICc[7] <- AICc(gt6_m)
tg_m_aic$Model[7] <- "m6"
tg_m_aic$LL[7] <- logLik(gt6_m)
tg_m_aic <- tg_m_aic[, 2:ncol(tg_m_aic)]


lg_m_aic <-  as.data.frame(1:7)
lg_m_aic$AICc <- "NA"
lg_m_aic$Model <- "NA"
lg_m_aic$LL <- "NA"
lg_m_aic$AICc[1] <- AICc(gnl_m)
lg_m_aic$Model[1] <- "Null"
lg_m_aic$LL[1] <- logLik(gnl_m)
lg_m_aic$AICc[2] <- AICc(gl1_m)
lg_m_aic$Model[2] <- "m1"
lg_m_aic$LL[2] <- logLik(gl1_m)
lg_m_aic$AICc[3] <- AICc(gl2_m)
lg_m_aic$Model[3] <- "m2"
lg_m_aic$LL[3] <- logLik(gl2_m)
lg_m_aic$AICc[4] <- AICc(gl3_m)
lg_m_aic$Model[4] <- "m3"
lg_m_aic$LL[4] <- logLik(gl3_m)
lg_m_aic$AICc[5] <- AICc(gl4_m)
lg_m_aic$Model[5] <- "m4"
lg_m_aic$LL[5] <- logLik(gl4_m)
lg_m_aic$AICc[6] <- AICc(gl5_m)
lg_m_aic$Model[6] <- "m5"
lg_m_aic$LL[6] <- logLik(gl5_m)
lg_m_aic$AICc[7] <- AICc(gl6_m)
lg_m_aic$Model[7] <- "m6"
lg_m_aic$LL[7] <- logLik(gl6_m)
lg_m_aic <- lg_m_aic[, 2:ncol(lg_m_aic)]


# Re-order and calculate Delta AICc
str(tg_m_aic)
tg_m_aic$AICc <- as.numeric(tg_m_aic$AICc)
tg_m_aic$LL <- as.numeric(tg_m_aic$LL)
str(tg_m_aic)
tg_m_aic <- tg_m_aic[order(tg_m_aic$AICc), ]
tg_m_aic # The best model is the model with interactive effect of fire frequency and seed weight.
tg_m_aic$Delta_AICc <- "0.00"
tg_m_aic$Delta_AICc[2] <- round(tg_m_aic$AICc[1]-tg_m_aic$AICc[2], 2)
tg_m_aic$Delta_AICc[3] <- round(tg_m_aic$AICc[1]-tg_m_aic$AICc[3], 2)
tg_m_aic$Delta_AICc[4] <- round(tg_m_aic$AICc[1]-tg_m_aic$AICc[4], 2)
tg_m_aic$Delta_AICc[5] <- round(tg_m_aic$AICc[1]-tg_m_aic$AICc[5], 2)
tg_m_aic$Delta_AICc[6] <- round(tg_m_aic$AICc[1]-tg_m_aic$AICc[6], 2)
tg_m_aic$Delta_AICc[7] <- round(tg_m_aic$AICc[1]-tg_m_aic$AICc[7], 2)
tg_m_aic
# As with other model 6 is best

str(lg_m_aic)
lg_m_aic$AICc <- as.numeric(lg_m_aic$AICc)
lg_m_aic$LL <- as.numeric(lg_m_aic$LL)
str(lg_m_aic)
lg_m_aic <- lg_m_aic[order(lg_m_aic$AICc), ]
lg_m_aic # The best model is the model with interactive effect of fire frequency and seed weight.
lg_m_aic$Delta_AICc <- "0.00"
lg_m_aic$Delta_AICc[2] <- round(lg_m_aic$AICc[1]-lg_m_aic$AICc[2], 2)
lg_m_aic$Delta_AICc[3] <- round(lg_m_aic$AICc[1]-lg_m_aic$AICc[3], 2)
lg_m_aic$Delta_AICc[4] <- round(lg_m_aic$AICc[1]-lg_m_aic$AICc[4], 2)
lg_m_aic$Delta_AICc[5] <- round(lg_m_aic$AICc[1]-lg_m_aic$AICc[5], 2)
lg_m_aic$Delta_AICc[6] <- round(lg_m_aic$AICc[1]-lg_m_aic$AICc[6], 2)
lg_m_aic$Delta_AICc[7] <- round(lg_m_aic$AICc[1]-lg_m_aic$AICc[7], 2)
lg_m_aic
# Model 2 but the null is ranked within delta AICc <2.



# 4.1.6 GAM with ti() default parameter smoothing and method = 'ML' ----

gl2_mc <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'cr', k = 5), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl2_mc)
par(mfrow = c(2,2))
gam.check(gl2_mc)
plot(gl2_mc)

gt6_mc <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'cr', k = 5) + s(Fire_freq, bs = 'cr', k = 5) + ti(seed_wt_mg, Fire_freq, bs = 'cr', k = 5), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt6_mc)
par(mfrow = c(2,2))
gam.check(gt6_mc)
plot(gt6_mc)




# GAM littoralis predictions -----
# Model 2 - seed weight main effect only 
# Predict to new data
# Defaults and s('tp', k = 10) are the same model so only ran default gam with smooth terms
newl <- data.frame(seed_wt_mg = seq(min(lit_cum.prop$seed_wt_mg), max(lit_cum.prop$seed_wt_mg), length = 50))
pl_s <- predict(gl2_s, newdata = newl, se.fit = T, type = 'response')
newl_s <- newl
newl_s$fit <- pl_s$fit
newl_s$se <- pl_s$se.fit
newl_s$lci <- newl_s$fit - (newl_s$se * 1.96)
newl_s$uci <- newl_s$fit + (newl_s$se * 1.96)



pl_s_cr <- predict(gl2_s_cr, newdata = newl, se.fit = T, type = 'response')
newl_s_cr <- newl
newl_s_cr$fit <- pl_s_cr$fit
newl_s_cr$se <- pl_s_cr$se.fit
newl_s_cr$lci <- newl_s_cr$fit - (newl_s_cr$se * 1.96)
newl_s_cr$uci <- newl_s_cr$fit + (newl_s_cr$se * 1.96)



pl_m <- predict(gl2_m, newdata = newl, se.fit = T, type = 'response')
newl_m <- newl
newl_m$fit <- pl_m$fit
newl_m$se <- pl_m$se.fit
newl_m$lci <- newl_m$fit - (newl_m$se * 1.96)
newl_m$uci <- newl_m$fit + (newl_m$se * 1.96)


pl_mc <- predict(gl2_mc, newdata = newl, se.fit = T, type = 'response')
newl_mc <- newl
newl_mc$fit <- pl_mc$fit
newl_mc$se <- pl_mc$se.fit
newl_mc$lci <- newl_mc$fit - (newl_mc$se * 1.96)
newl_mc$uci <- newl_mc$fit + (newl_mc$se * 1.96)



# Plot GLMER with GAMs with different bs and k
dev.new(width = 20, height = 5, noRStudioGD = T, dpi = 300)
par(mfrow = c(1,3), mar = c(6,6,3,2))

plot(new_lp$r_seed_wt, new_lp$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xaxt = "n", xlim = c(-2.8, 2.6), yaxt = "n")
axis(side = 1, at = seq(-2.7, 2.7, 0.3), labels = F)
text(seq(-2.7, 2.7, 0.3), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "r_seed_wt", data = "new_lp", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(a) GLMER")), cex = 1.5)

plot(newl_s$seed_wt_mg, newl_s$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_s', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(b) GAM + default s(bs = 'tp', k = 10)")), cex = 1.5)


plot(newl_s_cr$seed_wt_mg, newl_s_cr$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_s_cr', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(c) GAM + default ti(bs = 'cr', k = 5)")), cex = 1.5)
# Only slight difference between tp or cr basis function



# Plot GLMER with GAMS with tp settings but different methods
dev.new(width = 14, height = 18, noRStudioGD = T, dpi = 300)
par(mfrow = c(3,2), mar = c(6,6,3,2))

plot(new_lp$r_seed_wt, new_lp$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xaxt = "n", xlim = c(-2.8, 2.6), yaxt = "n")
axis(side = 1, at = seq(-2.7, 2.7, 0.3), labels = F)
text(seq(-2.7, 2.7, 0.3), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "r_seed_wt", data = "new_lp", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(a) GLMER")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(m2_l), 3), sep = ""), line = -8, cex = 1.2)


plot(lit_cum.prop$seed_wt_mg, lit_cum.prop$Proportion_germ, ylab = "", las = 1, pch = 19, xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
mtext(side = 3, expression(bold("(b) Raw")), cex = 1.5)


plot(newl_s$seed_wt_mg, newl_s$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_s', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(c) GAM + 'tp' k=10, method = 'UBRE'")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_s), 3), sep = ""), line = -8, cex = 1.2)


plot(newl_m$seed_wt_mg, newl_m$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_m', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(d) GAM + 'tp' k=10,method = 'ML'")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_m), 3), sep = ""), line = -8, cex = 1.2)



plot(newl_s_cr$seed_wt_mg, newl_s_cr$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_s_cr', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(e) GAM + 'cr' k=5,method = 'UBRE'")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_s_cr), 3), sep = ""), line = -8, cex = 1.2)


plot(newl_mc$seed_wt_mg, newl_mc$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_mc', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(d) GAM + 'cr' k=5,method = 'ML'")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_mc), 3), sep = ""), line = -8, cex = 1.2)


# UBRE/GCV method results may overestimate the likelihood that small seeds would have 100% germination, ML seems more appropriate,  with the fit more similar to the GLMER for bs = 'tp'. Basis function 'cr' does not seem as approrpraite for the litteroalis data for this question.


# Torulosa GAM predictions ----
# For predictions, we will choose 3 fire frequencies within the range of the data - min, mean, max
summary(tor_cum.prop$r_fire_freq) # Note the mean here is different to the mean provided by mean(tor_cum.prop$r_fire_freq)

gt6_s
gt6_s_tp
gt6_s_cr


# Default smooth parameters
newtl_s <- data.frame(seed_wt_mg = seq(min(tor_cum.prop$seed_wt_mg), max(tor_cum.prop$seed_wt_mg), length = 50),
                    Fire_freq = min(tor_cum.prop$Fire_freq))

ptl_s <- predict(gt6_s, newdata = newtl_s, se.fit = T, type = 'response')
newtl_s_s <- newtl_s
newtl_s_s$fit <- ptl_s$fit
newtl_s_s$se <- ptl_s$se.fit
newtl_s_s$lci <- newtl_s_s$fit - (newtl_s_s$se * 1.96)
newtl_s_s$uci <- newtl_s_s$fit + (newtl_s_s$se * 1.96)


newta_s <- data.frame(seed_wt_mg = seq(min(tor_cum.prop$seed_wt_mg), max(tor_cum.prop$seed_wt_mg), length = 50),
                    Fire_freq = 3)
pta_s <- predict(gt6_s, newdata = newta_s, se.fit = T, type = 'response')
newta_s_s <- newta_s
newta_s_s$fit <- pta_s$fit
newta_s_s$se <- pta_s$se.fit
newta_s_s$lci <- newta_s_s$fit - (newta_s_s$se * 1.96)
newta_s_s$uci <- newta_s_s$fit + (newta_s_s$se * 1.96)



newth_s <- data.frame(seed_wt_mg = seq(min(tor_cum.prop$seed_wt_mg), max(tor_cum.prop$seed_wt_mg), length = 50),
                    Fire_freq = max(tor_cum.prop$Fire_freq))
pth_s <- predict(gt6_s, newdata = newth_s, se.fit = T, type = 'response')
newth_s_s <- newth_s
newth_s_s$fit <- pth_s$fit
newth_s_s$se <- pth_s$se.fit
newth_s_s$lci <- newth_s_s$fit - (newth_s_s$se * 1.96)
newth_s_s$uci <- newth_s_s$fit + (newth_s_s$se * 1.96)



# Default s() smoothing parameters with tp
ptl_tp <- predict(gt6_s_tp, newdata = newtl_s, se.fit = T, type = 'response')
newtl_s_tp <- newtl_s
newtl_s_tp$fit <- ptl_tp$fit
newtl_s_tp$se <- ptl_tp$se.fit
newtl_s_tp$lci <- newtl_s_tp$fit - (newtl_s_tp$se * 1.96)
newtl_s_tp$uci <- newtl_s_tp$fit + (newtl_s_tp$se * 1.96)


pta_tp <- predict(gt6_s_tp, newdata = newta_s, se.fit = T, type = 'response')
newta_s_tp <- newta_s
newta_s_tp$fit <- pta_tp$fit
newta_s_tp$se <- pta_tp$se.fit
newta_s_tp$lci <- newta_s_tp$fit - (newta_s_tp$se * 1.96)
newta_s_tp$uci <- newta_s_tp$fit + (newta_s_tp$se * 1.96)


pth_tp <- predict(gt6_s_tp, newdata = newth_s, se.fit = T, type = 'response')
newth_s_tp <- newth_s
newth_s_tp$fit <- pth_tp$fit
newth_s_tp$se <- pth_tp$se.fit
newth_s_tp$lci <- newth_s_tp$fit - (newth_s_tp$se * 1.96)
newth_s_tp$uci <- newth_s_tp$fit + (newth_s_tp$se * 1.96)





# Default ti() smoothing parameters with cr
ptl_cr <- predict(gt6_s_cr, newdata = newtl_s, se.fit = T, type = 'response')
newtl_s_cr <- newtl_s
newtl_s_cr$fit <- ptl_cr$fit
newtl_s_cr$se <- ptl_cr$se.fit
newtl_s_cr$lci <- newtl_s_cr$fit - (newtl_s_cr$se * 1.96)
newtl_s_cr$uci <- newtl_s_cr$fit + (newtl_s_cr$se * 1.96)


pta_cr <- predict(gt6_s_cr, newdata = newta_s, se.fit = T, type = 'response')
newta_s_cr <- newta_s
newta_s_cr$fit <- pta_cr$fit
newta_s_cr$se <- pta_cr$se.fit
newta_s_cr$lci <- newta_s_cr$fit - (newta_s_cr$se * 1.96)
newta_s_cr$uci <- newta_s_cr$fit + (newta_s_cr$se * 1.96)


pth_cr <- predict(gt6_s_cr, newdata = newth_s, se.fit = T, type = 'response')
newth_s_cr <- newth_s
newth_s_cr$fit <- pth_cr$fit
newth_s_cr$se <- pth_cr$se.fit
newth_s_cr$lci <- newth_s_cr$fit - (newth_s_cr$se * 1.96)
newth_s_cr$uci <- newth_s_cr$fit + (newth_s_cr$se * 1.96)



# Default s() smoothing parameters with tp and method = ML
ptl_m <- predict(gt6_m, newdata = newtl_s, se.fit = T, type = 'response')
newtl_m <- newtl_s
newtl_m$fit <- ptl_m$fit
newtl_m$se <- ptl_m$se.fit
newtl_m$lci <- newtl_m$fit - (newtl_m$se * 1.96)
newtl_m$uci <- newtl_m$fit + (newtl_m$se * 1.96)


pta_m <- predict(gt6_m, newdata = newta_s, se.fit = T, type = 'response')
newta_m <- newta_s
newta_m$fit <- pta_m$fit
newta_m$se <- pta_m$se.fit
newta_m$lci <- newta_m$fit - (newta_m$se * 1.96)
newta_m$uci <- newta_m$fit + (newta_m$se * 1.96)


pth_m <- predict(gt6_m, newdata = newth_s, se.fit = T, type = 'response')
newth_m <- newth_s
newth_m$fit <- pth_m$fit
newth_m$se <- pth_m$se.fit
newth_m$lci <- newth_m$fit - (newth_m$se * 1.96)
newth_m$uci <- newth_m$fit + (newth_m$se * 1.96)




# Default ti() smoothing parameters with cr and method = ML
ptl_mc <- predict(gt6_mc, newdata = newtl_s, se.fit = T, type = 'response')
newtl_mc <- newtl_s
newtl_mc$fit <- ptl_mc$fit
newtl_mc$se <- ptl_mc$se.fit
newtl_mc$lci <- newtl_mc$fit - (newtl_mc$se * 1.96)
newtl_mc$uci <- newtl_mc$fit + (newtl_mc$se * 1.96)


pta_mc <- predict(gt6_mc, newdata = newta_s, se.fit = T, type = 'response')
newta_mc <- newta_s
newta_mc$fit <- pta_mc$fit
newta_mc$se <- pta_mc$se.fit
newta_mc$lci <- newta_mc$fit - (newta_mc$se * 1.96)
newta_mc$uci <- newta_mc$fit + (newta_mc$se * 1.96)


pth_mc <- predict(gt6_mc, newdata = newth_s, se.fit = T, type = 'response')
newth_mc <- newth_s
newth_mc$fit <- pth_mc$fit
newth_mc$se <- pth_mc$se.fit
newth_mc$lci <- newth_mc$fit - (newth_mc$se * 1.96)
newth_mc$uci <- newth_mc$fit + (newth_mc$se * 1.96)



# Torulosa GAM plots -----
dev.new(width = 20, height = 5, noRStudioGD = T, dpi = 300)
par(mfrow = c(1,3), mar = c(6,6,2,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))
labs <- seq(2.2, 7.2, 0.2)

plot(newtl_s_s$seed_wt_mg, newtl_s_s$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_s", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_s_s$seed_wt_mg, newth_s_s$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_s', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) Default s() and ti() parameters")), side = 3, line = 0.3, cex = 1.2)


plot(newtl_s_tp$seed_wt_mg, newtl_s_tp$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_tp", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_s_tp$seed_wt_mg, newth_s_tp$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_tp', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b) bs = 'tp' and k = 10 parameters")), side = 3, line = 0.3, cex = 1.2)



plot(newtl_s_cr$seed_wt_mg, newtl_s_cr$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_cr", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_s_cr$seed_wt_mg, newth_s_cr$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_cr', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) bs = 'cr' and k = 5 parameters")), side = 3, line = 0.3, cex = 1.2)


par(xpd = NA)
legend(x = 7.4, y = 1, legend = c("0 fires", "6 fires"), col = c("blue", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)
# Feel the tp default smoothing parameters perform better than the cr. 


dev.new(width = 20, height = 5, noRStudioGD = T, dpi = 300)
par(mfrow = c(1,3), mar = c(6,6,2,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))
labs <- seq(2.2, 7.2, 0.2)

plot(newtl_s_s$seed_wt_mg, newtl_s_s$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_s", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newta_s_s$seed_wt_mg, newta_s_s$fit, col = 'black')
pg.ci(x = 'seed_wt_mg', data = 'newta_s_s', colour = rgb(0, 0, 0, 0.1), lower = 'lci', upper = 'uci')

lines(newth_s_s$seed_wt_mg, newth_s_s$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_s', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) Default s() and ti() parameters")), side = 3, line = 0.3, cex = 1.2)


plot(newtl_s_tp$seed_wt_mg, newtl_s_tp$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_tp", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newta_s_tp$seed_wt_mg, newta_s_tp$fit, col = 'black')
pg.ci(x = 'seed_wt_mg', data = 'newta_s_tp', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')

lines(newth_s_tp$seed_wt_mg, newth_s_tp$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_tp', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b) bs = 'tp' and k = 10 parameters")), side = 3, line = 0.3, cex = 1.2)



plot(newtl_s_cr$seed_wt_mg, newtl_s_cr$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_cr", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newta_s_cr$seed_wt_mg, newta_s_cr$fit, col = 'black')
pg.ci(x = 'seed_wt_mg', data = 'newta_s_cr', colour = rgb(0, 0, 0, 0.1), lower = 'lci', upper = 'uci')

lines(newth_s_cr$seed_wt_mg, newth_s_cr$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_cr', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) bs = 'cr' and k = 5 parameters")), side = 3, line = 0.3, cex = 1.2)


par(xpd = NA)
legend(x = 7.4, y = 1, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", "black", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)





# Compare GAMs to the GLMMER ----
dev.new(width = 20, height = 5, noRStudioGD = T, dpi = 300)
par(mfrow = c(1,3), mar = c(6,6,2,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))
labs <- seq(2.2, 7.2, 0.2)



plot(new_tpl$r_seed_wt, new_tpl$fit, ylim = c(0, 1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(-2.2, 2.6), col = "blue", xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = F)
text(seq(-2.2, 2.6, 0.19), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0, 1, 0.1,), las = 1, cex.axis = 1.4)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "r_seed_wt", data = "new_tpl", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(new_tph$r_seed_wt, new_tph$fit, col = 'red')
pg.ci(x = 'r_seed_wt', data = 'new_tph', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) GLMMER")), side = 3, line = 0.3, cex = 1.2)


plot(newtl_s_tp$seed_wt_mg, newtl_s_tp$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_tp", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_s_tp$seed_wt_mg, newth_s_tp$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_tp', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b) GAM + default s(bs = 'tp', k = 10)")), side = 3, line = 0.3, cex = 1.2)



plot(newtl_s_cr$seed_wt_mg, newtl_s_cr$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_cr", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_s_cr$seed_wt_mg, newth_s_cr$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_cr', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) GAM + default s(bs = 'cr', k = 5)")), side = 3, line = 0.3, cex = 1.2)


par(xpd = NA)
legend(x = 7.4, y = 1, legend = c("0 fires", "6 fires"), col = c("blue", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)




dev.new(width = 20, height = 5, noRStudioGD = T, dpi = 300)
par(mfrow = c(1,3), mar = c(6,6,2,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))
labs <- seq(2.2, 7.2, 0.2)



plot(new_tpl$r_seed_wt, new_tpl$fit, ylim = c(0, 1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(-2.2, 2.6), col = "blue", xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = F)
text(seq(-2.2, 2.6, 0.19), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0, 1, 0.1,), las = 1, cex.axis = 1.4)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "r_seed_wt", data = "new_tpl", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(new_tpa$r_seed_wt, new_tpa$fit, col = 'black')
pg.ci(x = 'r_seed_wt', data = 'new_tpa', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')

lines(new_tph$r_seed_wt, new_tph$fit, col = 'red')
pg.ci(x = 'r_seed_wt', data = 'new_tph', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) GLMMER")), side = 3, line = 0.3, cex = 1.2)


plot(newtl_s_tp$seed_wt_mg, newtl_s_tp$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_tp", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newta_s_tp$seed_wt_mg, newta_s_tp$fit, col = 'black')
pg.ci(x = 'seed_wt_mg', data = 'newta_s_tp', colour  = rgb(0, 0, 0, 0.1), lower = 'lci', upper = 'uci')

lines(newth_s_tp$seed_wt_mg, newth_s_tp$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_tp', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b) GAM + default s(bs = 'tp', k = 10)")), side = 3, line = 0.3, cex = 1.2)



plot(newtl_s_cr$seed_wt_mg, newtl_s_cr$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_cr", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newta_s_cr$seed_wt_mg, newta_s_cr$fit, col = 'black')
pg.ci(x = 'seed_wt_mg', data = 'newta_s_cr', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')

lines(newth_s_cr$seed_wt_mg, newth_s_cr$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_cr', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) GAM + default s(bs = 'cr', k = 5)")), side = 3, line = 0.3, cex = 1.2)


par(xpd = NA)
legend(x = 7.4, y = 1, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", 'black', 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)
# Feel the tp default smoothing parameters perform better than the cr. 


save.image('./02_Workspaces/002_full_analysis_exploration.RData')



# Compare GLMER and s() default smoothing for methods UBRE and ML
dev.new(width = 18, height = 18, noRStudioGD = T, dpi = 300)
par(mfrow = c(3,2), mar = c(6,6,2,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))
labs <- seq(2.2, 7.2, 0.2)



plot(new_tpl$r_seed_wt, new_tpl$fit, ylim = c(0, 1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(-2.2, 2.6), col = "blue", xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = F)
text(seq(-2.2, 2.6, 0.19), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0, 1, 0.1,), las = 1, cex.axis = 1.4)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "r_seed_wt", data = "new_tpl", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(new_tph$r_seed_wt, new_tph$fit, col = 'red')
pg.ci(x = 'r_seed_wt', data = 'new_tph', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) GLMMER")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(m6_t), 3), sep = ""), line = -4, cex = 1.2)


plot(tor_cum.prop$seed_wt_mg[tor_cum.prop$Fire_freq == 0], tor_cum.prop$Proportion_germ[tor_cum.prop$Fire_freq == 0], pch = 19, col = 'blue', ylim = c(0,1), ylab = "", las = 1, xlab = "",xaxt = "n", yaxt = "n", xlim = c(2.2, 7.2))
points(tor_cum.prop$seed_wt_mg[tor_cum.prop$Fire_freq == 6], tor_cum.prop$Proportion_germ[tor_cum.prop$Fire_freq == 6], pch = 19, col = 'red')
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
mtext(expression(bold("(b) Raw")), side = 3, line = 0.3, cex = 1.2)


par(xpd = NA)
legend(x = 7.4, y = 1, legend = c("0 fires", "6 fires"), col = c("blue", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)


plot(newtl_s_tp$seed_wt_mg, newtl_s_tp$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_tp", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_s_tp$seed_wt_mg, newth_s_tp$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_tp', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) GAM + tp k=10, method = 'UBRE'")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_s_tp), 3), sep = ""), line = -4, cex = 1.2)




plot(newtl_m$seed_wt_mg, newtl_m$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_m", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_m$seed_wt_mg, newth_m$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_m', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(d) GAM + tp k=10 method = 'ML'")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_m), 3), sep = ""), line = -4, cex = 1.2)





plot(newtl_s_cr$seed_wt_mg, newtl_s_cr$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_s_cr", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_s_cr$seed_wt_mg, newth_s_cr$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_s_cr', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(e) GAM + cr k = 5, method = 'UBRE")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_s_cr), 3), sep = ""), line = -4, cex = 1.2)



plot(newtl_mc$seed_wt_mg, newtl_mc$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_mc", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_mc$seed_wt_mg, newth_mc$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_mc', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(f) GAM +  cr k=5, method = 'ML'")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_mc), 3), sep = ""), line = -4, cex = 1.2)

# If we are considering which method to use, 'ML' is more appropriate as UBRE/GCV is definitely overestimating and producing overly wiggly lines. Without optimisation this would also suggest that the bs = 'cr' is more appropriate as this reflects better the response that we would see with the GLMER. 












# Test k and bs changes for best models only ----
# Each time we run these analyses we are told the model with only seed weight for littoralis and the model with seed weight and fire frequency is best for torulosa. So what if we take these models and we test changing k for a tp basis function until we get good fit, then lets investigate using other basis functions with a default k

gl2_k <- gam(Proportion_germ ~ s(seed_wt_mg, k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl2_k)
par(mfrow = c(2,2)); gam.check(gl2_k) # QQ plot isn't very good 
plot(gl2_k) # Crossing of intervals at 0 with default k

gt6_k <- gam(Proportion_germ ~ s(seed_wt_mg, k = 3) + s(Fire_freq, k = 3) + ti(seed_wt_mg, Fire_freq, bs = 'tp', k = 5), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt6_k) 
par(mfrow = c(2,2)); gam.check(gt6_k) # QQ plot not very good
plot(gt6_k) # Crossing of intervals at 0 with default k



gl2_ts <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ts', k = 3), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl2_ts)
par(mfrow = c(2,2)); gam.check(gl2_ts)
plot(gl2_ts) # We cannot optimise smoothing by adjusting k

gt6_ts <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ts') + s(Fire_freq, bs= 'ts', k = 7) + ti(seed_wt_mg, Fire_freq, bs = 'ts'), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt6_ts)
par(mfrow = c(2,2)); gam.check(gt6_ts)
plot(gt6_ts) # We cannot optimise smoothing by adjusting k 



# Default is a cubic spline
gl2_ps <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ps', k = 7), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl2_ps)
par(mfrow = c(2,2)); gam.check(gl2_ps)
plot(gl2_ps) # Default k here is just as good, this is the smallest k which preserves the shape of CI and kept it bounded between -1 and 2

gt6_ps <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ps', k = 4) + s(Fire_freq, bs= 'ps', k = 7) + ti(seed_wt_mg, Fire_freq, bs = 'ps'), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt6_ps)
par(mfrow = c(2,2)); gam.check(gt6_ps)
plot(gt6_ps) # Default is fine here too, but adjusting k for seed weight allowed the smoother to bound between -1 and 1. Adjusting k for fire frequency does not change the smoother but removes a warning that basis dimension is larger than the number of unique covariates. 



gl2_bs <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'bs', k = 7), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl2_bs)
par(mfrow = c(2,2)); gam.check(gl2_bs)
plot(gl2_bs) # Default k here is just as good, this is the smallest k which preserves the shape of CI and kept it bounded between -1 and 2

gt6_bs <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'bs', k = 6) + s(Fire_freq, bs= 'bs', k = 7) + ti(seed_wt_mg, Fire_freq, bs = 'bs'), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt6_bs)
par(mfrow = c(2,2)); gam.check(gt6_bs)
plot(gt6_bs) 



gl2_ck <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'cr', k = 4), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl2_ck)
par(mfrow = c(2,2)); gam.check(gl2_ck)
plot(gl2_ck)


gt6_ck <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'cr', k = 6) + s(Fire_freq, bs = 'cr', k = 4) + ti(seed_wt_mg, Fire_freq, bs = 'cr'), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt6_ck)
par(mfrow = c(2,2)); gam.check(gt6_ck)
plot(gt6_ck) 



gl2_ad <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ad', k = 25), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl2_ad)
par(mfrow = c(2,2)); gam.check(gl2_ad)
plot(gl2_ad) # Default k here is preferable but we get a warning basis dimension is larger than number of unique covariates, reducing k to 25 preserves the CI shape but removes this warning.

gt6_ad <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ad', k = 8) + s(Fire_freq, bs = 'ad', k = 8) + ti(seed_wt_mg, Fire_freq, bs = 'tp', k = 5), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt6_ad)
par(mfrow = c(2,2)); gam.check(gt6_ad)
plot(gt6_ad) # Unable to completely remove the warning that the basis dimension is larger than number of unique covariates. 


gt6_adc <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ad', k = 8) + s(Fire_freq, bs = 'ad', k = 8) + ti(seed_wt_mg, Fire_freq, bs = 'cr', k = 5), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gt6_adc)
par(mfrow = c(2,2)); gam.check(gt6_adc)
plot(gt6_adc) # Unable to completely remove the warning that the basis dimension is larger than number of unique covariates. 
 

# Check if these models are better than the null
AICc(gnl_m)
AICc(gl2_k)
AICc(gl2_ts)
AICc(gl2_ps)
AICc(gl2_bs)
AICc(gl2_ck)
AICc(gl2_ad)
# All better than the null


AICc(gnt_m)
AICc(gt6_k)
AICc(gt6_ts)
AICc(gt6_ps)
AICc(gt6_bs)
AICc(gt6_ck)
AICc(gt6_ad)
AICc(gt6_adc)
# All better than the null


# Predict to new data for different bs -----
# Littoralis
pl_k <- predict(gl2_k, newdata = newl, se.fit = T, type = 'response')
newl_k <- newl
newl_k$fit <- pl_k$fit
newl_k$se <- pl_k$se.fit
newl_k$lci <- newl_k$fit - (newl_k$se * 1.96)
newl_k$uci <- newl_k$fit + (newl_k$se * 1.96)


pl_ck <- predict(gl2_ck, newdata = newl, se.fit = T, type = 'response')
newl_ck <- newl
newl_ck$fit <- pl_ck$fit
newl_ck$se <- pl_ck$se.fit
newl_ck$lci <- newl_ck$fit - (newl_ck$se * 1.96)
newl_ck$uci <- newl_ck$fit + (newl_ck$se * 1.96)

pl_ts <- predict(gl2_ts, newdata = newl, se.fit = T, type = 'response')
newl_ts <- newl
newl_ts$fit <- pl_ts$fit
newl_ts$se <- pl_ts$se.fit
newl_ts$lci <- newl_ts$fit - (newl_ts$se * 1.96)
newl_ts$uci <- newl_ts$fit + (newl_ts$se * 1.96)


pl_ps <- predict(gl2_ps, newdata = newl, se.fit = T, type = 'response')
newl_ps <- newl
newl_ps$fit <- pl_ps$fit
newl_ps$se <- pl_ps$se.fit
newl_ps$lci <- newl_ps$fit - (newl_ps$se * 1.96)
newl_ps$uci <- newl_ps$fit + (newl_ps$se * 1.96)


pl_bs <- predict(gl2_bs, newdata = newl, se.fit = T, type = 'response')
newl_bs <- newl
newl_bs$fit <- pl_bs$fit
newl_bs$se <- pl_bs$se.fit
newl_bs$lci <- newl_bs$fit - (newl_bs$se * 1.96)
newl_bs$uci <- newl_bs$fit + (newl_bs$se * 1.96)


pl_ad <- predict(gl2_ad, newdata = newl, se.fit = T, type = 'response')
newl_ad <- newl
newl_ad$fit <- pl_ad$fit
newl_ad$se <- pl_ad$se.fit
newl_ad$lci <- newl_ad$fit - (newl_ad$se * 1.96)
newl_ad$uci <- newl_ad$fit + (newl_ad$se * 1.96)



# Prediction plots
# Optimised k for tp only
dev.new(width = 20, height = 5, noRStudioGD = T, dpi = 300)
par(mfrow = c(1,3), mar = c(6,6,3,2))

plot(new_lp$r_seed_wt, new_lp$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xaxt = "n", xlim = c(-2.8, 2.6), yaxt = "n")
axis(side = 1, at = seq(-2.7, 2.7, 0.3), labels = F)
text(seq(-2.7, 2.7, 0.3), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "r_seed_wt", data = "new_lp", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(a) GLMER")), cex = 1.5)


plot(newl_m$seed_wt_mg, newl_m$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_m', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(b) GAM + 'tp' k=10, method = 'ML'")), cex = 1.5)
AICc(gt6_s_tp)

plot(newl_k$seed_wt_mg, newl_k$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_k', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(c) GAM 'tp' k optimised, method = 'ML'")), cex = 1.5)
AICc(gt6_k)

# Optimising k is not necessarily producing a better fit, also indicated by AICc. We may be overestimating proportion germination at low seed weights.



# Compare all different basis functions to GLMER and tp default and tp optimised.
dev.new(width = 20, height = 15, noRStudioGD = T, dpi = 300)
par(mfrow = c(3,4), mar = c(6,6,3,2))


plot(new_lp$r_seed_wt, new_lp$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xaxt = "n", xlim = c(-2.8, 2.6), yaxt = "n")
axis(side = 1, at = seq(-2.7, 2.7, 0.3), labels = F)
text(seq(-2.7, 2.7, 0.3), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "r_seed_wt", data = "new_lp", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(a) GLMER")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(m2_l), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_m$seed_wt_mg, newl_m$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_m', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(b) GAM tp default")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_s), 3), sep = ""), line = -10, cex = 1.2)

plot(newl_k$seed_wt_mg, newl_k$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_k', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(c) GAM tp optimised")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_k), 3), sep = ""), line = -10, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gl2_k), 3), sep = ""), line = -10, cex = 1.2)



plot(newl_mc$seed_wt_mg, newl_mc$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_mc', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(d) GAM cr default")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_mc), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_ck$seed_wt_mg, newl_ck$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_ck', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(e) GAM ck optimised")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_ck), 3), sep = ""), line = -10, cex = 1.2)

plot(newl_ts$seed_wt_mg, newl_ts$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_ts', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(f) GAM ts optimised")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_ts), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_ps$seed_wt_mg, newl_ps$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_ps', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(g) GAM ps optimised")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_ps), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_bs$seed_wt_mg, newl_bs$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_bs', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(h) GAM bs optimised")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_bs), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_ad$seed_wt_mg, newl_ad$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_ad', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(i) GAM ad optimised")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_ad), 3), sep = ""), line = -10, cex = 1.2)



plot(lit_cum.prop$seed_wt_mg, lit_cum.prop$Proportion_germ, ylab = "", las = 1, pch = 19, xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
mtext(side = 3, expression(bold("(j) Raw")), cex = 1.5)


# Proportion germination for littoralis is likely best modelled with GLMER or GAM with default tp settings



# Torulosa
# tp optimised
ptl_k <- predict(gt6_k, newdata = newtl_s, se.fit = T, type = 'response')
newtl_k <- newtl_s
newtl_k$fit <- ptl_k$fit
newtl_k$se <- ptl_k$se.fit
newtl_k$lci <- newtl_k$fit - (newtl_k$se * 1.96)
newtl_k$uci <- newtl_k$fit + (newtl_k$se * 1.96)


pta_k <- predict(gt6_k, newdata = newta_s, se.fit = T, type = 'response')
newta_k <- newta_s
newta_k$fit <- pta_k$fit
newta_k$se <- pta_k$se.fit
newta_k$lci <- newta_k$fit - (newta_k$se * 1.96)
newta_k$uci <- newta_k$fit + (newta_k$se * 1.96)


pth_k <- predict(gt6_k, newdata = newth_s, se.fit = T, type = 'response')
newth_k <- newth_s
newth_k$fit <- pth_k$fit
newth_k$se <- pth_k$se.fit
newth_k$lci <- newth_k$fit - (newth_k$se * 1.96)
newth_k$uci <- newth_k$fit + (newth_k$se * 1.96)



# ck optimised
ptl_ck <- predict(gt6_ck, newdata = newtl_s, se.fit = T, type = 'response')
newtl_ck <- newtl_s
newtl_ck$fit <- ptl_ck$fit
newtl_ck$se <- ptl_ck$se.fit
newtl_ck$lci <- newtl_ck$fit - (newtl_ck$se * 1.96)
newtl_ck$uci <- newtl_ck$fit + (newtl_ck$se * 1.96)


pta_ck <- predict(gt6_ck, newdata = newta_s, se.fit = T, type = 'response')
newta_ck <- newta_s
newta_ck$fit <- pta_ck$fit
newta_ck$se <- pta_ck$se.fit
newta_ck$lci <- newta_ck$fit - (newta_ck$se * 1.96)
newta_ck$uci <- newta_ck$fit + (newta_ck$se * 1.96)


pth_ck <- predict(gt6_ck, newdata = newth_s, se.fit = T, type = 'response')
newth_ck <- newth_s
newth_ck$fit <- pth_ck$fit
newth_ck$se <- pth_ck$se.fit
newth_ck$lci <- newth_ck$fit - (newth_ck$se * 1.96)
newth_ck$uci <- newth_ck$fit + (newth_ck$se * 1.96)

# ts
ptl_ts <- predict(gt6_ts, newdata = newtl_s, se.fit = T, type = 'response')
newtl_ts <- newtl_s
newtl_ts$fit <- ptl_ts$fit
newtl_ts$se <- ptl_ts$se.fit
newtl_ts$lci <- newtl_ts$fit - (newtl_ts$se * 1.96)
newtl_ts$uci <- newtl_ts$fit + (newtl_ts$se * 1.96)


pta_ts <- predict(gt6_ts, newdata = newta_s, se.fit = T, type = 'response')
newta_ts <- newta_s
newta_ts$fit <- pta_ts$fit
newta_ts$se <- pta_ts$se.fit
newta_ts$lci <- newta_ts$fit - (newta_ts$se * 1.96)
newta_ts$uci <- newta_ts$fit + (newta_ts$se * 1.96)


pth_ts <- predict(gt6_ts, newdata = newth_s, se.fit = T, type = 'response')
newth_ts <- newth_s
newth_ts$fit <- pth_ts$fit
newth_ts$se <- pth_ts$se.fit
newth_ts$lci <- newth_ts$fit - (newth_ts$se * 1.96)
newth_ts$uci <- newth_ts$fit + (newth_ts$se * 1.96)

# ps
ptl_ps <- predict(gt6_ps, newdata = newtl_s, se.fit = T, type = 'response')
newtl_ps <- newtl_s
newtl_ps$fit <- ptl_ps$fit
newtl_ps$se <- ptl_ps$se.fit
newtl_ps$lci <- newtl_ps$fit - (newtl_ps$se * 1.96)
newtl_ps$uci <- newtl_ps$fit + (newtl_ps$se * 1.96)


pta_ps <- predict(gt6_ps, newdata = newta_s, se.fit = T, type = 'response')
newta_ps <- newta_s
newta_ps$fit <- pta_ps$fit
newta_ps$se <- pta_ps$se.fit
newta_ps$lci <- newta_ps$fit - (newta_ps$se * 1.96)
newta_ps$uci <- newta_ps$fit + (newta_ps$se * 1.96)


pth_ps <- predict(gt6_ps, newdata = newth_s, se.fit = T, type = 'response')
newth_ps <- newth_s
newth_ps$fit <- pth_ps$fit
newth_ps$se <- pth_ps$se.fit
newth_ps$lci <- newth_ps$fit - (newth_ps$se * 1.96)
newth_ps$uci <- newth_ps$fit + (newth_ps$se * 1.96)

# bs
ptl_bs <- predict(gt6_bs, newdata = newtl_s, se.fit = T, type = 'response')
newtl_bs <- newtl_s
newtl_bs$fit <- ptl_bs$fit
newtl_bs$se <- ptl_bs$se.fit
newtl_bs$lci <- newtl_bs$fit - (newtl_bs$se * 1.96)
newtl_bs$uci <- newtl_bs$fit + (newtl_bs$se * 1.96)


pta_bs <- predict(gt6_bs, newdata = newta_s, se.fit = T, type = 'response')
newta_bs <- newta_s
newta_bs$fit <- pta_bs$fit
newta_bs$se <- pta_bs$se.fit
newta_bs$lci <- newta_bs$fit - (newta_bs$se * 1.96)
newta_bs$uci <- newta_bs$fit + (newta_bs$se * 1.96)


pth_bs <- predict(gt6_bs, newdata = newth_s, se.fit = T, type = 'response')
newth_bs <- newth_s
newth_bs$fit <- pth_bs$fit
newth_bs$se <- pth_bs$se.fit
newth_bs$lci <- newth_bs$fit - (newth_bs$se * 1.96)
newth_bs$uci <- newth_bs$fit + (newth_bs$se * 1.96)

# ad with ti('tp' optimised)
ptl_ad <- predict(gt6_ad, newdata = newtl_s, se.fit = T, type = 'response')
newtl_ad <- newtl_s
newtl_ad$fit <- ptl_ad$fit
newtl_ad$se <- ptl_ad$se.fit
newtl_ad$lci <- newtl_ad$fit - (newtl_ad$se * 1.96)
newtl_ad$uci <- newtl_ad$fit + (newtl_ad$se * 1.96)


pta_ad <- predict(gt6_ad, newdata = newta_s, se.fit = T, type = 'response')
newta_ad <- newta_s
newta_ad$fit <- pta_ad$fit
newta_ad$se <- pta_ad$se.fit
newta_ad$lci <- newta_ad$fit - (newta_ad$se * 1.96)
newta_ad$uci <- newta_ad$fit + (newta_ad$se * 1.96)


pth_ad <- predict(gt6_ad, newdata = newth_s, se.fit = T, type = 'response')
newth_ad <- newth_s
newth_ad$fit <- pth_ad$fit
newth_ad$se <- pth_ad$se.fit
newth_ad$lci <- newth_ad$fit - (newth_ad$se * 1.96)
newth_ad$uci <- newth_ad$fit + (newth_ad$se * 1.96)

# ad with ti ('cr' optimised)
ptl_adc <- predict(gt6_adc, newdata = newtl_s, se.fit = T, type = 'response')
newtl_adc <- newtl_s
newtl_adc$fit <- ptl_adc$fit
newtl_adc$se <- ptl_adc$se.fit
newtl_adc$lci <- newtl_adc$fit - (newtl_adc$se * 1.96)
newtl_adc$uci <- newtl_adc$fit + (newtl_adc$se * 1.96)


pta_adc <- predict(gt6_adc, newdata = newta_s, se.fit = T, type = 'response')
newta_adc <- newta_s
newta_adc$fit <- pta_adc$fit
newta_adc$se <- pta_adc$se.fit
newta_adc$lci <- newta_adc$fit - (newta_adc$se * 1.96)
newta_adc$uci <- newta_adc$fit + (newta_adc$se * 1.96)


pth_adc <- predict(gt6_adc, newdata = newth_s, se.fit = T, type = 'response')
newth_adc <- newth_s
newth_adc$fit <- pth_adc$fit
newth_adc$se <- pth_adc$se.fit
newth_adc$lci <- newth_adc$fit - (newth_adc$se * 1.96)
newth_adc$uci <- newth_adc$fit + (newth_adc$se * 1.96)


# Prediction plots



# Compare all different basis functions to GLMER including tp and cr default and tp and cr optimised.
dev.new(width = 25, height = 15, noRStudioGD = T, dpi = 300)
par(mfrow = c(3,4), mar = c(6,6,3,2))
labs <- seq(2.2, 7.2, 0.2)


plot(new_tpl$r_seed_wt, new_tpl$fit, ylim = c(0, 1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(-2.2, 2.6), col = "blue", xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = F)
text(seq(-2.2, 2.6, 0.19), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0, 1, 0.1,), las = 1, cex.axis = 1.4)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "r_seed_wt", data = "new_tpl", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(new_tph$r_seed_wt, new_tph$fit, col = 'red')
pg.ci(x = 'r_seed_wt', data = 'new_tph', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) GLMER")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(m2_t), 3), sep = ""), line = -1.5, cex = 1.2)


plot(newtl_m$seed_wt_mg, newtl_m$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_m", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_m$seed_wt_mg, newth_m$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_m', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b) GAM 'tp' k=10")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_m), 3), sep = ""), line = -1.5, cex = 1.2)



plot(newtl_k$seed_wt_mg, newtl_k$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_k", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_k$seed_wt_mg, newth_k$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_k', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) GAM 'tp' k optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_k), 3), sep = ""), line = -1.5, cex = 1.2)




plot(newtl_mc$seed_wt_mg, newtl_mc$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_mc", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_mc$seed_wt_mg, newth_mc$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_mc', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(d) GAM 'cr' default")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_mc), 3), sep = ""), line = -1.5, cex = 1.2)




plot(newtl_ck$seed_wt_mg, newtl_ck$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_ck", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_ck$seed_wt_mg, newth_ck$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_ck', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(e) GAM 'ck' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_ck), 3), sep = ""), line = -1.5, cex = 1.2)


plot(newtl_ts$seed_wt_mg, newtl_ts$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_ts", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_ts$seed_wt_mg, newth_ts$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_ts', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(f) GAM 'ts' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_ts), 3), sep = ""), line = -1.5, cex = 1.2)


plot(newtl_ps$seed_wt_mg, newtl_ps$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_ps", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_ps$seed_wt_mg, newth_ps$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_ps', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(g) GAM 'ps' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_ps), 3), sep = ""), line = -1.5, cex = 1.2)


plot(newtl_bs$seed_wt_mg, newtl_bs$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_bs", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_bs$seed_wt_mg, newth_bs$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_bs', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(h) GAM 'bs' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_bs), 3), sep = ""), line = -1.5, cex = 1.2)


plot(newtl_ad$seed_wt_mg, newtl_ad$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_ad", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_ad$seed_wt_mg, newth_ad$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_ad', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(i) GAM 'ad' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_ad), 3), sep = ""), line = -1.5, cex = 1.2)



plot(newtl_adc$seed_wt_mg, newtl_adc$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_adc", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_adc$seed_wt_mg, newth_adc$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_adc', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(j) GAM 'adc' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_adc), 3), sep = ""), line = -1.5, cex = 1.2)



plot(tor_cum.prop$seed_wt_mg[tor_cum.prop$Fire_freq == 0], tor_cum.prop$Proportion_germ[tor_cum.prop$Fire_freq == 0], pch = 19, col = 'blue', ylim = c(0,1), ylab = "", las = 1, xlab = "",xaxt = "n", yaxt = "n", xlim = c(2.2, 7.2))
points(tor_cum.prop$seed_wt_mg[tor_cum.prop$Fire_freq == 6], tor_cum.prop$Proportion_germ[tor_cum.prop$Fire_freq == 6], pch = 19, col = 'red')
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
mtext(expression(bold("(k) Raw")), side = 3, line = 0.3, cex = 1.2)



par(xpd = NA)
legend(x = 9, y = 1, legend = c("0 fires", "6 fires"), col = c("blue", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)



# For littoralis, while the optimised bs = 'ps' has the lowest AIC, the model fit is not very good. Our best model fit for a GAM would be bs = 'bs', which has a slightly better AIC than the cr default (i.e., the default basis function and k used by ti for all covariates).  

# Compare different basis functions without optimisation -----


gl2_tsd <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ts'), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
gt6_tsd <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ts') + s(Fire_freq, bs= 'ts', k = 7) + ti(seed_wt_mg, Fire_freq, bs = 'ts'), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
AICc(gl2_ts); AICc(gl2_tsd)
AICc(gt6_ts); AICc(gt6_tsd) # Identical


# Default is a cubic spline
gl2_psd <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ps'), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
gt6_psd <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ps') + s(Fire_freq, bs= 'ps') + ti(seed_wt_mg, Fire_freq, bs = 'ps'), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
AICc(gl2_ps); AICc(gl2_psd)
AICc(gt6_ps); AICc(gt6_psd) # Identical


gl2_bsd <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'bs'), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
gt6_bsd <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'bs') + s(Fire_freq, bs= 'bs') + ti(seed_wt_mg, Fire_freq, bs = 'bs'), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
AICc(gl2_bs); AICc(gl2_bsd)
AICc(gt6_bs); AICc(gt6_bsd)


gl2_add <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ad'), family = binomial, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(gl2_ad)
gt6_add <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ad') + s(Fire_freq, bs = 'ad') + ti(seed_wt_mg, Fire_freq, bs = 'tp'), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
AICc(gl2_ad); AICc(gl2_add)
AICc(gt6_ad); AICc(gt6_add)

gt6_adcd <- gam(Proportion_germ ~ s(seed_wt_mg, bs = 'ad') + s(Fire_freq, bs = 'ad') + ti(seed_wt_mg, Fire_freq, bs = 'cr'), family = binomial, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
AICc(gt6_adc); AICc(gt6_adcd)




# Predict to new data for different bs -----
# Littoralis
pl_tsd <- predict(gl2_tsd, newdata = newl, se.fit = T, type = 'response')
newl_tsd <- newl
newl_tsd$fit <- pl_tsd$fit
newl_tsd$se <- pl_tsd$se.fit
newl_tsd$lci <- newl_tsd$fit - (newl_tsd$se * 1.96)
newl_tsd$uci <- newl_tsd$fit + (newl_tsd$se * 1.96)


pl_psd <- predict(gl2_psd, newdata = newl, se.fit = T, type = 'response')
newl_psd <- newl
newl_psd$fit <- pl_psd$fit
newl_psd$se <- pl_psd$se.fit
newl_psd$lci <- newl_psd$fit - (newl_psd$se * 1.96)
newl_psd$uci <- newl_psd$fit + (newl_psd$se * 1.96)


pl_bsd <- predict(gl2_bsd, newdata = newl, se.fit = T, type = 'response')
newl_bsd <- newl
newl_bsd$fit <- pl_bsd$fit
newl_bsd$se <- pl_bsd$se.fit
newl_bsd$lci <- newl_bsd$fit - (newl_bsd$se * 1.96)
newl_bsd$uci <- newl_bsd$fit + (newl_bsd$se * 1.96)


pl_add <- predict(gl2_add, newdata = newl, se.fit = T, type = 'response')
newl_add <- newl
newl_add$fit <- pl_add$fit
newl_add$se <- pl_add$se.fit
newl_add$lci <- newl_add$fit - (newl_add$se * 1.96)
newl_add$uci <- newl_add$fit + (newl_add$se * 1.96)



# Compare all different basis functions at defaults to GLMER
dev.new(width = 20, height = 15, noRStudioGD = T, dpi = 300)
par(mfrow = c(3,4), mar = c(6,6,3,2))


plot(new_lp$r_seed_wt, new_lp$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xaxt = "n", xlim = c(-2.8, 2.6), yaxt = "n")
axis(side = 1, at = seq(-2.7, 2.7, 0.3), labels = F)
text(seq(-2.7, 2.7, 0.3), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "r_seed_wt", data = "new_lp", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(a) GLMER")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(m2_l), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_m$seed_wt_mg, newl_m$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_m', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(b) GAM tp default")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_s), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_mc$seed_wt_mg, newl_mc$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_mc', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(c) GAM cr default")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_mc), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_ts$seed_wt_mg, newl_ts$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_ts', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(d) GAM ts optimised")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_ts), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_tsd$seed_wt_mg, newl_tsd$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_tsd', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(d) GAM ts default")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_tsd), 3), sep = ""), line = -10, cex = 1.2)



plot(newl_ps$seed_wt_mg, newl_ps$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_ps', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(e) GAM ps optimised")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_ps), 3), sep = ""), line = -10, cex = 1.2)



plot(newl_psd$seed_wt_mg, newl_psd$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_psd', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(e) GAM ps default")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_psd), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_bs$seed_wt_mg, newl_bs$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_bs', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(f) GAM bs optimised")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_bs), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_bsd$seed_wt_mg, newl_bsd$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_bsd', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(g) GAM bs default")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_bsd), 3), sep = ""), line = -10, cex = 1.2)


plot(newl_ad$seed_wt_mg, newl_ad$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_ad', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(h) GAM ad optimised")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_ad), 3), sep = ""), line = -10, cex = 1.2)



plot(newl_add$seed_wt_mg, newl_add$fit, ylab = "", las = 1, type = 'l', xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n", ylim = c(0,1))
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = 'newl_add', colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(side = 3, expression(bold("(i) GAM ad default")), cex = 1.5)
mtext(paste("AICc = ", round(AICc(gl2_add), 3), sep = ""), line = -10, cex = 1.2)



plot(lit_cum.prop$seed_wt_mg, lit_cum.prop$Proportion_germ, ylab = "", las = 1, pch = 19, xlab = "", xlim = c(1.3, 3.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(1.3, 3.1, 0.1), labels = F)
text(seq(1.3, 3.1, 0.1), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3.5,cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
mtext(side = 3, expression(bold("(j) Raw")), cex = 1.5)


# Proportion germination for littoralis is likely best modelled with GLMER or GAM with default tp settings
# Optimised and defaults are pretty similar. Does not change the best models 


# Torulosa
# bs default
ptl_bsd <- predict(gt6_bsd, newdata = newtl_s, se.fit = T, type = 'response')
newtl_bsd <- newtl_s
newtl_bsd$fit <- ptl_bsd$fit
newtl_bsd$se <- ptl_bsd$se.fit
newtl_bsd$lci <- newtl_bsd$fit - (newtl_bsd$se * 1.96)
newtl_bsd$uci <- newtl_bsd$fit + (newtl_bsd$se * 1.96)


pta_bsd <- predict(gt6_bsd, newdata = newta_s, se.fit = T, type = 'response')
newta_bsd <- newta_s
newta_bsd$fit <- pta_bsd$fit
newta_bsd$se <- pta_bsd$se.fit
newta_bsd$lci <- newta_bsd$fit - (newta_bsd$se * 1.96)
newta_bsd$uci <- newta_bsd$fit + (newta_bsd$se * 1.96)


pth_bsd <- predict(gt6_bsd, newdata = newth_s, se.fit = T, type = 'response')
newth_bsd <- newth_s
newth_bsd$fit <- pth_bsd$fit
newth_bsd$se <- pth_bsd$se.fit
newth_bsd$lci <- newth_bsd$fit - (newth_bsd$se * 1.96)
newth_bsd$uci <- newth_bsd$fit + (newth_bsd$se * 1.96)

# add with tp default
ptl_add <- predict(gt6_add, newdata = newtl_s, se.fit = T, type = 'response')
newtl_add <- newtl_s
newtl_add$fit <- ptl_add$fit
newtl_add$se <- ptl_add$se.fit
newtl_add$lci <- newtl_add$fit - (newtl_add$se * 1.96)
newtl_add$uci <- newtl_add$fit + (newtl_add$se * 1.96)


pta_add <- predict(gt6_add, newdata = newta_s, se.fit = T, type = 'response')
newta_add <- newta_s
newta_add$fit <- pta_add$fit
newta_add$se <- pta_add$se.fit
newta_add$lci <- newta_add$fit - (newta_add$se * 1.96)
newta_add$uci <- newta_add$fit + (newta_add$se * 1.96)


pth_add <- predict(gt6_add, newdata = newth_s, se.fit = T, type = 'response')
newth_add <- newth_s
newth_add$fit <- pth_add$fit
newth_add$se <- pth_add$se.fit
newth_add$lci <- newth_add$fit - (newth_add$se * 1.96)
newth_add$uci <- newth_add$fit + (newth_add$se * 1.96)

# add with ti defauld
ptl_adcd <- predict(gt6_adcd, newdata = newtl_s, se.fit = T, type = 'response')
newtl_adcd <- newtl_s
newtl_adcd$fit <- ptl_adcd$fit
newtl_adcd$se <- ptl_adcd$se.fit
newtl_adcd$lci <- newtl_adcd$fit - (newtl_adcd$se * 1.96)
newtl_adcd$uci <- newtl_adcd$fit + (newtl_adcd$se * 1.96)


pta_adcd <- predict(gt6_adcd, newdata = newta_s, se.fit = T, type = 'response')
newta_adcd <- newta_s
newta_adcd$fit <- pta_adcd$fit
newta_adcd$se <- pta_adcd$se.fit
newta_adcd$lci <- newta_adcd$fit - (newta_adcd$se * 1.96)
newta_adcd$uci <- newta_adcd$fit + (newta_adcd$se * 1.96)


pth_adcd <- predict(gt6_adcd, newdata = newth_s, se.fit = T, type = 'response')
newth_adcd <- newth_s
newth_adcd$fit <- pth_adcd$fit
newth_adcd$se <- pth_adcd$se.fit
newth_adcd$lci <- newth_adcd$fit - (newth_adcd$se * 1.96)
newth_adcd$uci <- newth_adcd$fit + (newth_adcd$se * 1.96)


# Prediction plots



# Compare all different basis functions to GLMER including tp and cr default and tp and cr optimised.
dev.new(width = 25, height = 15, noRStudioGD = T, dpi = 300)
par(mfrow = c(3,4), mar = c(6,6,3,2))
labs <- seq(2.2, 7.2, 0.2)


plot(new_tpl$r_seed_wt, new_tpl$fit, ylim = c(0, 1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(-2.2, 2.6), col = "blue", xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = F)
text(seq(-2.2, 2.6, 0.19), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0, 1, 0.1,), las = 1, cex.axis = 1.4)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "r_seed_wt", data = "new_tpl", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(new_tph$r_seed_wt, new_tph$fit, col = 'red')
pg.ci(x = 'r_seed_wt', data = 'new_tph', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) GLMER")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(m2_t), 3), sep = ""), line = -1.5, cex = 1.2)




plot(newtl_k$seed_wt_mg, newtl_k$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_k", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_k$seed_wt_mg, newth_k$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_k', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b) GAM 'tp' k optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_k), 3), sep = ""), line = -1.5, cex = 1.2)




plot(newtl_mc$seed_wt_mg, newtl_mc$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_mc", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_mc$seed_wt_mg, newth_mc$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_mc', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) GAM 'cr' default")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_mc), 3), sep = ""), line = -1.5, cex = 1.2)




plot(newtl_ck$seed_wt_mg, newtl_ck$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_ck", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_ck$seed_wt_mg, newth_ck$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_ck', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(d) GAM 'ck' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_ck), 3), sep = ""), line = -1.5, cex = 1.2)


plot(newtl_ts$seed_wt_mg, newtl_ts$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_ts", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_ts$seed_wt_mg, newth_ts$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_ts', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(e) GAM 'ts' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_ts), 3), sep = ""), line = -1.5, cex = 1.2)


plot(newtl_ps$seed_wt_mg, newtl_ps$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_ps", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_ps$seed_wt_mg, newth_ps$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_ps', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(f) GAM 'ps' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_ps), 3), sep = ""), line = -1.5, cex = 1.2)


plot(newtl_bs$seed_wt_mg, newtl_bs$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_bs", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_bs$seed_wt_mg, newth_bs$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_bs', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(g) GAM 'bs' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_bs), 3), sep = ""), line = -1.5, cex = 1.2)


plot(newtl_bsd$seed_wt_mg, newtl_bsd$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_bsd", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_bsd$seed_wt_mg, newth_bsd$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_bs', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(h) GAM 'bs' default")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_bsd), 3), sep = ""), line = -1.5, cex = 1.2)


plot(newtl_ad$seed_wt_mg, newtl_ad$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_ad", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_ad$seed_wt_mg, newth_ad$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_ad', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(i) GAM 'ad' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_ad), 3), sep = ""), line = -1.5, cex = 1.2)

plot(newtl_add$seed_wt_mg, newtl_add$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_add", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_add$seed_wt_mg, newth_add$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_add', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(j) GAM 'ad' default")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_add), 3), sep = ""), line = -1.5, cex = 1.2)



plot(newtl_adc$seed_wt_mg, newtl_adc$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_adc", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_adc$seed_wt_mg, newth_adc$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_adc', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(k) GAM 'adc' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_adc), 3), sep = ""), line = -1.5, cex = 1.2)


plot(newtl_adcd$seed_wt_mg, newtl_adcd$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_adcd", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newth_adcd$seed_wt_mg, newth_adcd$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_adc', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(l) GAM 'adc' defauld")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_adcd), 3), sep = ""), line = -1.5, cex = 1.2)


# Optimisation here definitely makes a difference, while the AICc of the bs = 'bs' default is better than the optimised, the fit for high fire frequency does not fall well within CIs



# Create plots with the GLMER, GAM bs optimised including the average, cr default and optimised  
dev.new(width = 16, height = 12, noRStudioGD = T, dpi = 300)
par(mfrow = c(2,2), mar = c(6,6,3,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))
labs <- seq(2.2, 7.2, 0.2)


plot(new_tpl$r_seed_wt, new_tpl$fit, ylim = c(0, 1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(-2.2, 2.6), col = "blue", xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = F)
text(seq(-2.2, 2.6, 0.19), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0, 1, 0.1,), las = 1, cex.axis = 1.4)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 2.5)
pg.ci(x = "r_seed_wt", data = "new_tpl", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(new_tpa$r_seed_wt, new_tpa$fit, col = 'black')
pg.ci(x = 'r_seed_wt', data = 'new_tpa', colour = rgb(0,0,0, 0.1), lower = 'lci', upper = 'uci')
lines(new_tph$r_seed_wt, new_tph$fit, col = 'red')
pg.ci(x = 'r_seed_wt', data = 'new_tph', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) GLMER")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(m2_t), 3), sep = ""), line = -1.5, cex = 1.2)



plot(newtl_bs$seed_wt_mg, newtl_bs$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_bs", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newta_bs$seed_wt_mg, newta_bs$fit, col = 'black')
pg.ci(x = 'seed_wt_mg', data = 'newta_bs', colour = rgb(0,0,0, 0.1), lower = 'lci', upper = 'uci')
lines(newth_bs$seed_wt_mg, newth_bs$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_bs', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b) GAM 'bs' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_bs), 3), sep = ""), line = -1.5, cex = 1.2)

par(xpd = NA)
legend(x = 7.4, y = 1, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", "black", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)

plot(newtl_mc$seed_wt_mg, newtl_mc$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_mc", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newta_mc$seed_wt_mg, newta_mc$fit, col = 'black')
pg.ci(x = 'seed_wt_mg', data = 'newta_mc', colour = rgb(0,0,0, 0.1), lower = 'lci', upper = 'uci')
lines(newth_mc$seed_wt_mg, newth_mc$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_mc', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) GAM 'cr' default")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_mc), 3), sep = ""), line = -1.5, cex = 1.2)




plot(newtl_ck$seed_wt_mg, newtl_ck$fit, ylim = c(0,1), ylab = "", las = 1, type = 'l', xlab = "", xlim = c(2.2, 7.2), col = 'blue', xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(2.2, 7.2, 0.2), labels = F)
text(seq(2.2, 7.2, 0.2), par("usr")[3]-0.02, labels = labs, srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(0,1, 0.1), las = 1, cex.axis = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 1.5)
pg.ci(x = "seed_wt_mg", data = "newtl_ck", colour = rgb(0/255, 0/255, 255/255, 0.1), lower = "lci", upper = "uci")

lines(newta_ck$seed_wt_mg, newta_ck$fit, col = 'black')
pg.ci(x = 'seed_wt_mg', data = 'newta_ck', colour = rgb(0,0,0, 0.1), lower = 'lci', upper = 'uci')
lines(newth_ck$seed_wt_mg, newth_ck$fit, col = 'red')
pg.ci(x = 'seed_wt_mg', data = 'newth_ck', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(d) GAM 'cr' optimised")), side = 3, line = 0.3, cex = 1.2)
mtext(paste("AICc = ", round(AICc(gt6_ck), 3), sep = ""), line = -1.5, cex = 1.2)







# 4.2 Time to 50% germination preliminary analyses ----
# 4.2.1 GLMER ----
mn_tt50 <- glmer(t50 ~ 1 + (1|Individual) + (1|Set), family = poisson, data = tor_cum.prop)
mn_lt50 <- glmer(t50 ~ 1 + (1|Individual) + (1|Set), family = poisson, data = lit_cum.prop)

m1_tt50 <- glmer(t50 ~ Treatment + (1|Individual) + (1|Set), family = poisson, data = tor_cum.prop)
summary(m1_tt50)
m1_lt50 <- glmer(t50 ~ Treatment + (1|Individual) + (1|Set), family = poisson, data = lit_cum.prop)
summary(m1_lt50)


m2_tt50 <- glmer(t50 ~ r_seed_wt + (1|Individual) + (1|Set), family = poisson, data = tor_cum.prop)
summary(m2_tt50)
m2_lt50 <- glmer(t50 ~ r_seed_wt + (1|Individual) + (1|Set), family = poisson, data = lit_cum.prop)
summary(m2_lt50)


m3_tt50 <- glmer(t50 ~ r_fire_freq + (1|Individual) + (1|Set), family = poisson, data = tor_cum.prop)
summary(m3_tt50)
m3_lt50 <- glmer(t50 ~ r_fire_freq + (1|Individual) + (1|Set), family = poisson, data = lit_cum.prop)
summary(m3_lt50)


m4_tt50 <- glmer(t50 ~ Treatment * r_seed_wt + (1|Individual) + (1|Set), family = poisson, data = tor_cum.prop)
summary(m4_t)
m4_lt50 <- glmer(t50 ~ Treatment * r_seed_wt + (1|Individual) + (1|Set), family = poisson, data = lit_cum.prop)
summary(m4_lt50)


m5_tt50 <- glmer(t50 ~ Treatment * r_fire_freq + (1|Individual) + (1|Set), family = poisson, data = tor_cum.prop)
summary(m5_tt50)
m5_lt50 <- glmer(t50 ~ Treatment * r_fire_freq + (1|Individual) + (1|Set), family = poisson, data = lit_cum.prop)
summary(m5_lt50)


m6_tt50 <- glmer(t50 ~ r_seed_wt * r_fire_freq + (1|Individual) + (1|Set), family = poisson, data = tor_cum.prop)
summary(m6_tt50)
m6_lt50 <- glmer(t50 ~ r_seed_wt * r_fire_freq + (1|Individual) + (1|Set), family = poisson, data = lit_cum.prop)
summary(m6_lt50)


# Compare model fit
tt_aic <- list(mn_tt50, m1_tt50, m2_tt50, m3_tt50, m4_tt50, m5_tt50, m6_tt50)
lt_aic <- list(mn_lt50, m1_lt50, m2_lt50, m3_lt50, m4_lt50, m5_lt50, m6_lt50)
aictab(tt_aic) # Models 2, 5, and 6 were ranked within delta AICc <2 and are better than the null
aictab(lt_aic) # Model 2 is better than the null but the null is within delta AICc <2. 



# 4.2.2 GAM with 'tp', k = 10 defaults ----
ttn <- gam(t50 ~ 1, family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
ltn <- gam(t50 ~ 1, family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')

tt1_s <- gam(t50 ~ Treatment, family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt1_s)
par(mfrow = c(2,2))
gam.check(tt1_s)
plot(tt1_s)
lt1_s <- gam(t50 ~ Treatment, family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt1_s)
par(mfrow = c(2,2))
gam.check(lt1_s)
plot(lt1_s)


tt2_s <- gam(t50 ~ s(seed_wt_mg), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt2_s)
par(mfrow = c(2,2))
gam.check(tt2_s)
plot(tt2_s)
lt2_s <- gam(t50 ~ s(seed_wt_mg), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt2_s)
par(mfrow = c(2,2))
gam.check(lt2_s)
plot(lt2_s)

tt3_s <- gam(t50 ~ s(Fire_freq, k = 7), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt3_s)
par(mfrow = c(2,2))
gam.check(tt3_s)
plot(tt3_s)
lt3_s <- gam(t50 ~ s(Fire_freq, k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML') 
summary(lt3_s)
par(mfrow = c(2,2))
gam.check(lt3_s)
plot(lt3_s)

tt4_s <- gam(t50 ~ s(Treatment, bs = 're') + s(seed_wt_mg) + ti(seed_wt_mg, by = Treatment, bs = 'tp', k = 10), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt4_s)
par(mfrow = c(2,2))
gam.check(tt4_s)
plot(tt4_s)
lt4_s <- gam(t50 ~ s(Treatment, bs ='re') + s(seed_wt_mg) + ti(seed_wt_mg, by = Treatment, bs = 'tp', k = 9), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt4_s)
par(mfrow = c(2,2))
gam.check(lt4_s)
plot(lt4_s)


tt5_s <- gam(t50 ~ s(Treatment, bs ='re') + s(Fire_freq, k = 7) + ti(Fire_freq, by = Treatment, bs = 'tp', k = 7), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt5_s)
par(mfrow = c(2,2))
gam.check(tt5_s)
plot(tt5_s)
lt5_s <- gam(t50 ~ s(Treatment, bs ='re') + s(Fire_freq, k = 3) + ti(Fire_freq, by = Treatment, bs = 'tp', k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt5_s)
par(mfrow = c(2,2))
gam.check(lt5_s)
plot(lt5_s)

tt6_s <- gam(t50 ~ s(seed_wt_mg) + s(Fire_freq, k = 7) + ti(seed_wt_mg, Fire_freq, bs = 'tp', k = 7), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt6_s)
par(mfrow = c(2,2))
gam.check(tt6_s)
plot(tt6_s)
lt6_s <- gam(t50 ~ s(seed_wt_mg) + s(Fire_freq, k = 3) + ti(seed_wt_mg, Fire_freq, bs = 'tp', k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt6_s)
par(mfrow = c(2,2))
gam.check(lt6_s)
plot(lt6_s)



# Compare model fit
tt_aic <-  as.data.frame(1:7)
tt_aic$AICc <- "NA"
tt_aic$Model <- "NA"
tt_aic$LL <- "NA"
tt_aic$AICc[1] <- AICc(ttn)
tt_aic$Model[1] <- "Null"
tt_aic$LL[1] <- logLik(ttn)
tt_aic$AICc[2] <- AICc(tt1_s)
tt_aic$Model[2] <- "m1"
tt_aic$LL[2] <- logLik(tt1_s)
tt_aic$AICc[3] <- AICc(tt2_s)
tt_aic$Model[3] <- "m2"
tt_aic$LL[3] <- logLik(tt2_s)
tt_aic$AICc[4] <- AICc(tt3_s)
tt_aic$Model[4] <- "m3"
tt_aic$LL[4] <- logLik(tt3_s)
tt_aic$AICc[5] <- AICc(tt4_s)
tt_aic$Model[5] <- "m4"
tt_aic$LL[5] <- logLik(tt4_s)
tt_aic$AICc[6] <- AICc(tt5_s)
tt_aic$Model[6] <- "m5"
tt_aic$LL[6] <- logLik(tt5_s)
tt_aic$AICc[7] <- AICc(tt6_s)
tt_aic$Model[7] <- "m6"
tt_aic$LL[7] <- logLik(tt6_s)
tt_aic <- tt_aic[, 2:ncol(tt_aic)]


lt_aic <-  as.data.frame(1:7)
lt_aic$AICc <- "NA"
lt_aic$Model <- "NA"
lt_aic$LL <- "NA"
lt_aic$AICc[1] <- AICc(ltn)
lt_aic$Model[1] <- "Null"
lt_aic$LL[1] <- logLik(ltn)
lt_aic$AICc[2] <- AICc(lt1_s)
lt_aic$Model[2] <- "m1"
lt_aic$LL[2] <- logLik(lt1_s)
lt_aic$AICc[3] <- AICc(lt2_s)
lt_aic$Model[3] <- "m2"
lt_aic$LL[3] <- logLik(lt2_s)
lt_aic$AICc[4] <- AICc(lt3_s)
lt_aic$Model[4] <- "m3"
lt_aic$LL[4] <- logLik(lt3_s)
lt_aic$AICc[5] <- AICc(lt4_s)
lt_aic$Model[5] <- "m4"
lt_aic$LL[5] <- logLik(lt4_s)
lt_aic$AICc[6] <- AICc(lt5_s)
lt_aic$Model[6] <- "m5"
lt_aic$LL[6] <- logLik(lt5_s)
lt_aic$AICc[7] <- AICc(lt6_s)
lt_aic$Model[7] <- "m6"
lt_aic$LL[7] <- logLik(lt6_s)
lt_aic <- lt_aic[, 2:ncol(lt_aic)]


# Re-order and calculate Delta AICc
str(tt_aic)
tt_aic$AICc <- as.numeric(tt_aic$AICc)
tt_aic$LL <- as.numeric(tt_aic$LL)
str(tt_aic)
tt_aic <- tt_aic[order(tt_aic$AICc), ]
tt_aic # The best model is the model with interactive effect of fire frequency and seed weight.
tt_aic$Delta_AICc <- "0.00"
tt_aic$Delta_AICc[2] <- round(tt_aic$AICc[1]-tt_aic$AICc[2], 2)
tt_aic$Delta_AICc[3] <- round(tt_aic$AICc[1]-tt_aic$AICc[3], 2)
tt_aic$Delta_AICc[4] <- round(tt_aic$AICc[1]-tt_aic$AICc[4], 2)
tt_aic$Delta_AICc[5] <- round(tt_aic$AICc[1]-tt_aic$AICc[5], 2)
tt_aic$Delta_AICc[6] <- round(tt_aic$AICc[1]-tt_aic$AICc[6], 2)
tt_aic$Delta_AICc[7] <- round(tt_aic$AICc[1]-tt_aic$AICc[7], 2)
tt_aic
# Unlike GLMER model 4 is best. No model ranked within delta AICc <2. The best model for the GLMER were model 2 ranked 5th, model 5 ranked 2nd and model 6 ranked 4th here.

str(lt_aic)
lt_aic$AICc <- as.numeric(lt_aic$AICc)
lt_aic$LL <- as.numeric(lt_aic$LL)
str(lt_aic)
lt_aic <- lt_aic[order(lt_aic$AICc), ]
lt_aic # The best model is the model with interactive effect of fire frequency and seed weight.
lt_aic$Delta_AICc <- "0.00"
lt_aic$Delta_AICc[2] <- round(lt_aic$AICc[1]-lt_aic$AICc[2], 2)
lt_aic$Delta_AICc[3] <- round(lt_aic$AICc[1]-lt_aic$AICc[3], 2)
lt_aic$Delta_AICc[4] <- round(lt_aic$AICc[1]-lt_aic$AICc[4], 2)
lt_aic$Delta_AICc[5] <- round(lt_aic$AICc[1]-lt_aic$AICc[5], 2)
lt_aic$Delta_AICc[6] <- round(lt_aic$AICc[1]-lt_aic$AICc[6], 2)
lt_aic$Delta_AICc[7] <- round(lt_aic$AICc[1]-lt_aic$AICc[7], 2)
lt_aic
# The best model here is model 1 with no other model ranked within AICc <2. This is also different to the GLMER which ranked model 2 best, but this is ranked 2nd by the GAMs.



# 4.2.3 GAM with 'cr',  k = 5 defaults -----
tt2_sc <- gam(t50 ~ s(seed_wt_mg, bs = 'cr', k = 5), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt2_sc)
par(mfrow = c(2,2))
gam.check(tt2_sc)
plot(tt2_sc)
lt2_sc <- gam(t50 ~ s(seed_wt_mg, bs = 'cr', k = 5), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt2_sc)
par(mfrow = c(2,2))
gam.check(lt2_sc)
plot(lt2_sc)

tt3_sc <- gam(t50 ~ s(Fire_freq, bs = 'cr', k = 5), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt3_sc)
par(mfrow = c(2,2))
gam.check(tt3_sc)
plot(tt3_sc)
lt3_sc <- gam(t50 ~ s(Fire_freq, bs = 'cr', k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML') 
summary(lt3_sc)
par(mfrow = c(2,2))
gam.check(lt3_sc)
plot(lt3_sc)

tt4_sc <- gam(t50 ~ s(Treatment, bs = 're') + s(seed_wt_mg, bs = 'cr', k = 5) + ti(seed_wt_mg, by = Treatment, bs = 'cr', k = 5), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt4_sc)
par(mfrow = c(2,2))
gam.check(tt4_sc)
plot(tt4_sc)
lt4_sc <- gam(t50 ~ s(Treatment, bs ='re') + s(seed_wt_mg, bs = 'cr', k = 5) + ti(seed_wt_mg, by = Treatment, bs = 'cr', k = 5), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt4_sc)
par(mfrow = c(2,2))
gam.check(lt4_sc)
plot(lt4_sc)


tt5_sc <- gam(t50 ~ s(Treatment, bs ='re') + s(Fire_freq, bs = 'cr', k = 5) + ti(Fire_freq, by = Treatment, bs = 'cr', k = 5), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt5_sc)
par(mfrow = c(2,2))
gam.check(tt5_sc)
plot(tt5_sc)
lt5_sc <- gam(t50 ~ s(Treatment, bs ='re') + s(Fire_freq, bs = 'cr', k = 3) + ti(Fire_freq, by = Treatment, bs = 'cr', k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt5_sc)
par(mfrow = c(2,2))
gam.check(lt5_sc)
plot(lt5_sc)

tt6_sc <- gam(t50 ~ s(seed_wt_mg) + s(Fire_freq, bs = 'cr', k = 5) + ti(seed_wt_mg, Fire_freq, bs = 'cr', k = 5), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt6_sc)
par(mfrow = c(2,2))
gam.check(tt6_sc)
plot(tt6_sc)
lt6_sc <- gam(t50 ~ s(seed_wt_mg) + s(Fire_freq, bs = 'cr', k = 3) + ti(seed_wt_mg, Fire_freq, bs = 'cr', k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt6_sc)
par(mfrow = c(2,2))
gam.check(lt6_sc)
plot(lt6_sc)



# Compare model fit
tt_c_aic <-  as.data.frame(1:7)
tt_c_aic$AICc <- "NA"
tt_c_aic$Model <- "NA"
tt_c_aic$LL <- "NA"
tt_c_aic$AICc[1] <- AICc(ttn)
tt_c_aic$Model[1] <- "Null"
tt_c_aic$LL[1] <- logLik(ttn)
tt_c_aic$AICc[2] <- AICc(tt1_s)
tt_c_aic$Model[2] <- "m1"
tt_c_aic$LL[2] <- logLik(tt1_s)
tt_c_aic$AICc[3] <- AICc(tt2_s)
tt_c_aic$Model[3] <- "m2"
tt_c_aic$LL[3] <- logLik(tt2_s)
tt_c_aic$AICc[4] <- AICc(tt3_sc)
tt_c_aic$Model[4] <- "m3"
tt_c_aic$LL[4] <- logLik(tt3_sc)
tt_c_aic$AICc[5] <- AICc(tt4_sc)
tt_c_aic$Model[5] <- "m4"
tt_c_aic$LL[5] <- logLik(tt4_sc)
tt_c_aic$AICc[6] <- AICc(tt5_sc)
tt_c_aic$Model[6] <- "m5"
tt_c_aic$LL[6] <- logLik(tt5_sc)
tt_c_aic$AICc[7] <- AICc(tt6_sc)
tt_c_aic$Model[7] <- "m6"
tt_c_aic$LL[7] <- logLik(tt6_sc)
tt_c_aic <- tt_c_aic[, 2:ncol(tt_c_aic)]


lt_c_aic <-  as.data.frame(1:7)
lt_c_aic$AICc <- "NA"
lt_c_aic$Model <- "NA"
lt_c_aic$LL <- "NA"
lt_c_aic$AICc[1] <- AICc(ltn)
lt_c_aic$Model[1] <- "Null"
lt_c_aic$LL[1] <- logLik(ltn)
lt_c_aic$AICc[2] <- AICc(lt1_s)
lt_c_aic$Model[2] <- "m1"
lt_c_aic$LL[2] <- logLik(lt1_s)
lt_c_aic$AICc[3] <- AICc(lt2_s)
lt_c_aic$Model[3] <- "m2"
lt_c_aic$LL[3] <- logLik(lt2_s)
lt_c_aic$AICc[4] <- AICc(lt3_sc)
lt_c_aic$Model[4] <- "m3"
lt_c_aic$LL[4] <- logLik(lt3_sc)
lt_c_aic$AICc[5] <- AICc(lt4_sc)
lt_c_aic$Model[5] <- "m4"
lt_c_aic$LL[5] <- logLik(lt4_sc)
lt_c_aic$AICc[6] <- AICc(lt5_sc)
lt_c_aic$Model[6] <- "m5"
lt_c_aic$LL[6] <- logLik(lt5_sc)
lt_c_aic$AICc[7] <- AICc(lt6_sc)
lt_c_aic$Model[7] <- "m6"
lt_c_aic$LL[7] <- logLik(lt6_sc)
lt_c_aic <- lt_c_aic[, 2:ncol(lt_c_aic)]


# Re-order and calculate Delta AICc
str(tt_c_aic)
tt_c_aic$AICc <- as.numeric(tt_c_aic$AICc)
tt_c_aic$LL <- as.numeric(tt_c_aic$LL)
str(tt_c_aic)
tt_c_aic <- tt_c_aic[order(tt_c_aic$AICc), ]
tt_c_aic # The best model is the model with interactive effect of fire frequency and seed weight.
tt_c_aic$Delta_AICc <- "0.00"
tt_c_aic$Delta_AICc[2] <- round(tt_c_aic$AICc[1]-tt_c_aic$AICc[2], 2)
tt_c_aic$Delta_AICc[3] <- round(tt_c_aic$AICc[1]-tt_c_aic$AICc[3], 2)
tt_c_aic$Delta_AICc[4] <- round(tt_c_aic$AICc[1]-tt_c_aic$AICc[4], 2)
tt_c_aic$Delta_AICc[5] <- round(tt_c_aic$AICc[1]-tt_c_aic$AICc[5], 2)
tt_c_aic$Delta_AICc[6] <- round(tt_c_aic$AICc[1]-tt_c_aic$AICc[6], 2)
tt_c_aic$Delta_AICc[7] <- round(tt_c_aic$AICc[1]-tt_c_aic$AICc[7], 2)
tt_c_aic # model 4 is best, with no model ranked within delta AICc <2. 

str(lt_c_aic)
lt_c_aic$AICc <- as.numeric(lt_c_aic$AICc)
lt_c_aic$LL <- as.numeric(lt_c_aic$LL)
str(lt_c_aic)
lt_c_aic <- lt_c_aic[order(lt_c_aic$AICc), ]
lt_c_aic # The best model is the model with interactive effect of fire frequency and seed weight.
lt_c_aic$Delta_AICc <- "0.00"
lt_c_aic$Delta_AICc[2] <- round(lt_c_aic$AICc[1]-lt_c_aic$AICc[2], 2)
lt_c_aic$Delta_AICc[3] <- round(lt_c_aic$AICc[1]-lt_c_aic$AICc[3], 2)
lt_c_aic$Delta_AICc[4] <- round(lt_c_aic$AICc[1]-lt_c_aic$AICc[4], 2)
lt_c_aic$Delta_AICc[5] <- round(lt_c_aic$AICc[1]-lt_c_aic$AICc[5], 2)
lt_c_aic$Delta_AICc[6] <- round(lt_c_aic$AICc[1]-lt_c_aic$AICc[6], 2)
lt_c_aic$Delta_AICc[7] <- round(lt_c_aic$AICc[1]-lt_c_aic$AICc[7], 2)
lt_c_aic # Model 1 is best with no model ranked within delta AICc <2.


# 4.2.4 GAM with 'tp' optimised ----
tt2_so <- gam(t50 ~ s(seed_wt_mg, k = 5), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt2_so)
par(mfrow = c(2,2))
gam.check(tt2_so)
plot(tt2_so)
lt2_so <- gam(t50 ~ s(seed_wt_mg), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt2_so)
par(mfrow = c(2,2))
gam.check(lt2_so)
plot(lt2_so) # Cannot be optimised

tt3_so <- gam(t50 ~ s(Fire_freq, k = 7), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt3_so)
par(mfrow = c(2,2))
gam.check(tt3_so)
plot(tt3_so) # Cannot be optimised
lt3_so <- gam(t50 ~ s(Fire_freq, k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML') 
summary(lt3_so)
par(mfrow = c(2,2))
gam.check(lt3_so)
plot(lt3_so) # Cannot be optimised

tt4_so <- gam(t50 ~ s(Treatment, bs = 're') + s(seed_wt_mg, k = 4) + ti(seed_wt_mg, by = Treatment, bs = 'tp', k = 10), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt4_so)
par(mfrow = c(2,2))
gam.check(tt4_so) # Doesn't really optimise fit
par(mfrow = c(3,3)); plot(tt4_so)
lt4_so <- gam(t50 ~ s(Treatment, bs ='re') + s(seed_wt_mg, k = 4) + ti(seed_wt_mg, by = Treatment, bs = 'tp'), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt4_so)
par(mfrow = c(2,2))
gam.check(lt4_so)
par(mfrow = c(3,3));plot(lt4_so) # Doesn't really optimise fit


tt5_so <- gam(t50 ~ s(Treatment, bs ='re') + s(Fire_freq, k = 5) + ti(Fire_freq, by = Treatment, bs = 'tp', k = c(5, 3, 3, 7, 3, 5)), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt5_so)
par(mfrow = c(2,2))
gam.check(tt5_so)
par(mfrow = c(3,3)); plot(tt5_so)
lt5_so <- gam(t50 ~ s(Treatment, bs ='re') + s(Fire_freq, k = 3) + ti(Fire_freq, by = Treatment, bs = 'tp', k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt5_so)
par(mfrow = c(2,2))
gam.check(lt5_so)
par(mfrow = c(3,3)); plot(lt5_so) # Cannot be optimised as k is already at the lowest it can be set.


tt6_so <- gam(t50 ~ s(seed_wt_mg) + s(Fire_freq, k = 3) + ti(seed_wt_mg, Fire_freq, bs = 'tp', k = 7), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt6_so)
par(mfrow = c(2,2))
gam.check(tt6_so)
plot(tt6_so)
lt6_so <- gam(t50 ~ s(seed_wt_mg) + s(Fire_freq, k = 3) + ti(seed_wt_mg, Fire_freq, bs = 'tp', k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt6_so)
par(mfrow = c(2,2))
gam.check(lt6_so)
plot(lt6_so) # Cannot optimise fit



# Compare model fit
tt_so_aic <-  as.data.frame(1:7)
tt_so_aic$AICc <- "NA"
tt_so_aic$Model <- "NA"
tt_so_aic$LL <- "NA"
tt_so_aic$AICc[1] <- AICc(ttn)
tt_so_aic$Model[1] <- "Null"
tt_so_aic$LL[1] <- logLik(ttn)
tt_so_aic$AICc[2] <- AICc(tt1_s)
tt_so_aic$Model[2] <- "m1"
tt_so_aic$LL[2] <- logLik(tt1_s)
tt_so_aic$AICc[3] <- AICc(tt2_s)
tt_so_aic$Model[3] <- "m2"
tt_so_aic$LL[3] <- logLik(tt2_s)
tt_so_aic$AICc[4] <- AICc(tt3_s)
tt_so_aic$Model[4] <- "m3"
tt_so_aic$LL[4] <- logLik(tt3_s)
tt_so_aic$AICc[5] <- AICc(tt4_so)
tt_so_aic$Model[5] <- "m4"
tt_so_aic$LL[5] <- logLik(tt4_so)
tt_so_aic$AICc[6] <- AICc(tt5_so)
tt_so_aic$Model[6] <- "m5"
tt_so_aic$LL[6] <- logLik(tt5_so)
tt_so_aic$AICc[7] <- AICc(tt6_so)
tt_so_aic$Model[7] <- "m6"
tt_so_aic$LL[7] <- logLik(tt6_so)
tt_so_aic <- tt_so_aic[, 2:ncol(tt_so_aic)]


lt_so_aic <-  as.data.frame(1:7)
lt_so_aic$AICc <- "NA"
lt_so_aic$Model <- "NA"
lt_so_aic$LL <- "NA"
lt_so_aic$AICc[1] <- AICc(ltn)
lt_so_aic$Model[1] <- "Null"
lt_so_aic$LL[1] <- logLik(ltn)
lt_so_aic$AICc[2] <- AICc(lt1_s)
lt_so_aic$Model[2] <- "m1"
lt_so_aic$LL[2] <- logLik(lt1_s)
lt_so_aic$AICc[3] <- AICc(lt2_s)
lt_so_aic$Model[3] <- "m2"
lt_so_aic$LL[3] <- logLik(lt2_s)
lt_so_aic$AICc[4] <- AICc(lt3_s)
lt_so_aic$Model[4] <- "m3"
lt_so_aic$LL[4] <- logLik(lt3_s)
lt_so_aic$AICc[5] <- AICc(lt4_so)
lt_so_aic$Model[5] <- "m4"
lt_so_aic$LL[5] <- logLik(lt4_so)
lt_so_aic$AICc[6] <- AICc(lt5_s)
lt_so_aic$Model[6] <- "m5"
lt_so_aic$LL[6] <- logLik(lt5_s)
lt_so_aic$AICc[7] <- AICc(lt6_s)
lt_so_aic$Model[7] <- "m6"
lt_so_aic$LL[7] <- logLik(lt6_s)
lt_so_aic <- lt_so_aic[, 2:ncol(lt_so_aic)]


# Re-order and calculate Delta AICc
str(tt_so_aic)
tt_so_aic$AICc <- as.numeric(tt_so_aic$AICc)
tt_so_aic$LL <- as.numeric(tt_so_aic$LL)
str(tt_so_aic)
tt_so_aic <- tt_so_aic[order(tt_so_aic$AICc), ]
tt_so_aic # The best model is the model with interactive effect of fire frequency and seed weight.
tt_so_aic$Delta_AICc <- "0.00"
tt_so_aic$Delta_AICc[2] <- round(tt_so_aic$AICc[1]-tt_so_aic$AICc[2], 2)
tt_so_aic$Delta_AICc[3] <- round(tt_so_aic$AICc[1]-tt_so_aic$AICc[3], 2)
tt_so_aic$Delta_AICc[4] <- round(tt_so_aic$AICc[1]-tt_so_aic$AICc[4], 2)
tt_so_aic$Delta_AICc[5] <- round(tt_so_aic$AICc[1]-tt_so_aic$AICc[5], 2)
tt_so_aic$Delta_AICc[6] <- round(tt_so_aic$AICc[1]-tt_so_aic$AICc[6], 2)
tt_so_aic$Delta_AICc[7] <- round(tt_so_aic$AICc[1]-tt_so_aic$AICc[7], 2)
tt_so_aic
# Model 4 is best

str(lt_so_aic)
lt_so_aic$AICc <- as.numeric(lt_so_aic$AICc)
lt_so_aic$LL <- as.numeric(lt_so_aic$LL)
str(lt_so_aic)
lt_so_aic <- lt_so_aic[order(lt_so_aic$AICc), ]
lt_so_aic # The best model is the model with interactive effect of fire frequency and seed weight.
lt_so_aic$Delta_AICc <- "0.00"
lt_so_aic$Delta_AICc[2] <- round(lt_so_aic$AICc[1]-lt_so_aic$AICc[2], 2)
lt_so_aic$Delta_AICc[3] <- round(lt_so_aic$AICc[1]-lt_so_aic$AICc[3], 2)
lt_so_aic$Delta_AICc[4] <- round(lt_so_aic$AICc[1]-lt_so_aic$AICc[4], 2)
lt_so_aic$Delta_AICc[5] <- round(lt_so_aic$AICc[1]-lt_so_aic$AICc[5], 2)
lt_so_aic$Delta_AICc[6] <- round(lt_so_aic$AICc[1]-lt_so_aic$AICc[6], 2)
lt_so_aic$Delta_AICc[7] <- round(lt_so_aic$AICc[1]-lt_so_aic$AICc[7], 2)
lt_so_aic
# Model 1 is best



# 4.2.5 GAM with 'cr' optimised ----
tt2_sco <- gam(t50 ~ s(seed_wt_mg, bs = 'cr', k = 5), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt2_sco)
par(mfrow = c(2,2))
gam.check(tt2_sco)
plot(tt2_sco) # Doesn't need optimisation
lt2_sco <- gam(t50 ~ s(seed_wt_mg, bs = 'cr', k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt2_sco)
par(mfrow = c(2,2))
gam.check(lt2_sco)
plot(lt2_sco)

tt3_sco <- gam(t50 ~ s(Fire_freq, bs = 'cr', k = 5), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt3_sco)
par(mfrow = c(2,2))
gam.check(tt3_sco)
plot(tt3_sco) # Cannot be optimised
lt3_sco <- gam(t50 ~ s(Fire_freq, bs = 'cr', k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML') 
summary(lt3_sco)
par(mfrow = c(2,2))
gam.check(lt3_sco)
plot(lt3_sco) # Cannot be optimised

tt4_sco <- gam(t50 ~ s(Treatment, bs = 're') + s(seed_wt_mg, bs = 'cr', k = 4) + ti(seed_wt_mg, by = Treatment, bs = 'cr', k = 5), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt4_sco)
par(mfrow = c(2,2))
gam.check(tt4_sco)
par(mfrow = c(3,3)); plot(tt4_sco)
lt4_sco <- gam(t50 ~ s(Treatment, bs ='re') + s(seed_wt_mg, bs = 'cr', k = 5) + ti(seed_wt_mg, by = Treatment, bs = 'cr', k = 5), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt4_sco)
par(mfrow = c(2,2))
gam.check(lt4_sco)
par(mfrow = c(3,3)); plot(lt4_sco) # Cannot be optimised


tt5_sco <- gam(t50 ~ s(Treatment, bs ='re') + s(Fire_freq, bs = 'cr', k = 6) + ti(Fire_freq, by = Treatment, bs = 'cr', k = 4), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt5_sco)
par(mfrow = c(2,2))
gam.check(tt5_sco)
par(mfrow = c(3,3)); plot(tt5_sco)
lt5_sco <- gam(t50 ~ s(Treatment, bs ='re') + s(Fire_freq, bs = 'cr', k = 3) + ti(Fire_freq, by = Treatment, bs = 'cr', k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt5_sco)
par(mfrow = c(2,2))
gam.check(lt5_sco)
par(mfrow = c(3,3)); plot(lt5_sco) # Cannot be optimised

tt6_sco <- gam(t50 ~ s(seed_wt_mg) + s(Fire_freq, bs = 'cr', k = 5) + ti(seed_wt_mg, Fire_freq, bs = 'cr', k = 5), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(tt6_sco)
par(mfrow = c(2,2))
gam.check(tt6_sco)
plot(tt6_sco) # Cannot be optimised
lt6_sco <- gam(t50 ~ s(seed_wt_mg) + s(Fire_freq, bs = 'cr', k = 3) + ti(seed_wt_mg, Fire_freq, bs = 'cr', k = 3), family = poisson, data = lit_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
summary(lt6_sco)
par(mfrow = c(2,2))
gam.check(lt6_sco)
plot(lt6_sco) # Cannot be optimised



# Compare model fit
tt_co_aic <-  as.data.frame(1:7)
tt_co_aic$AICc <- "NA"
tt_co_aic$Model <- "NA"
tt_co_aic$LL <- "NA"
tt_co_aic$AICc[1] <- AICc(ttn)
tt_co_aic$Model[1] <- "Null"
tt_co_aic$LL[1] <- logLik(ttn)
tt_co_aic$AICc[2] <- AICc(tt1_s)
tt_co_aic$Model[2] <- "m1"
tt_co_aic$LL[2] <- logLik(tt1_s)
tt_co_aic$AICc[3] <- AICc(tt2_sc)
tt_co_aic$Model[3] <- "m2"
tt_co_aic$LL[3] <- logLik(tt2_sc)
tt_co_aic$AICc[4] <- AICc(tt3_sc)
tt_co_aic$Model[4] <- "m3"
tt_co_aic$LL[4] <- logLik(tt3_sc)
tt_co_aic$AICc[5] <- AICc(tt4_sco)
tt_co_aic$Model[5] <- "m4"
tt_co_aic$LL[5] <- logLik(tt4_sco)
tt_co_aic$AICc[6] <- AICc(tt5_sco)
tt_co_aic$Model[6] <- "m5"
tt_co_aic$LL[6] <- logLik(tt5_sco)
tt_co_aic$AICc[7] <- AICc(tt6_sc)
tt_co_aic$Model[7] <- "m6"
tt_co_aic$LL[7] <- logLik(tt6_sc)
tt_co_aic <- tt_co_aic[, 2:ncol(tt_co_aic)]


lt_co_aic <-  as.data.frame(1:7)
lt_co_aic$AICc <- "NA"
lt_co_aic$Model <- "NA"
lt_co_aic$LL <- "NA"
lt_co_aic$AICc[1] <- AICc(ltn)
lt_co_aic$Model[1] <- "Null"
lt_co_aic$LL[1] <- logLik(ltn)
lt_co_aic$AICc[2] <- AICc(lt1_s)
lt_co_aic$Model[2] <- "m1"
lt_co_aic$LL[2] <- logLik(lt1_s)
lt_co_aic$AICc[3] <- AICc(lt2_sco)
lt_co_aic$Model[3] <- "m2"
lt_co_aic$LL[3] <- logLik(lt2_sco)
lt_co_aic$AICc[4] <- AICc(lt3_sc)
lt_co_aic$Model[4] <- "m3"
lt_co_aic$LL[4] <- logLik(lt3_sc)
lt_co_aic$AICc[5] <- AICc(lt4_sc)
lt_co_aic$Model[5] <- "m4"
lt_co_aic$LL[5] <- logLik(lt4_sc)
lt_co_aic$AICc[6] <- AICc(lt5_sc)
lt_co_aic$Model[6] <- "m5"
lt_co_aic$LL[6] <- logLik(lt5_sc)
lt_co_aic$AICc[7] <- AICc(lt6_sc)
lt_co_aic$Model[7] <- "m6"
lt_co_aic$LL[7] <- logLik(lt6_sc)
lt_co_aic <- lt_co_aic[, 2:ncol(lt_co_aic)]


# Re-order and calculate Delta AICc
str(tt_co_aic)
tt_co_aic$AICc <- as.numeric(tt_co_aic$AICc)
tt_co_aic$LL <- as.numeric(tt_co_aic$LL)
str(tt_co_aic)
tt_co_aic <- tt_co_aic[order(tt_co_aic$AICc), ]
tt_co_aic # The best model is the model with interactive effect of fire frequency and seed weight.
tt_co_aic$Delta_AICc <- "0.00"
tt_co_aic$Delta_AICc[2] <- round(tt_co_aic$AICc[1]-tt_co_aic$AICc[2], 2)
tt_co_aic$Delta_AICc[3] <- round(tt_co_aic$AICc[1]-tt_co_aic$AICc[3], 2)
tt_co_aic$Delta_AICc[4] <- round(tt_co_aic$AICc[1]-tt_co_aic$AICc[4], 2)
tt_co_aic$Delta_AICc[5] <- round(tt_co_aic$AICc[1]-tt_co_aic$AICc[5], 2)
tt_co_aic$Delta_AICc[6] <- round(tt_co_aic$AICc[1]-tt_co_aic$AICc[6], 2)
tt_co_aic$Delta_AICc[7] <- round(tt_co_aic$AICc[1]-tt_co_aic$AICc[7], 2)
tt_co_aic # model 4 is best, with no model ranked within delta AICc <2. 

str(lt_co_aic)
lt_co_aic$AICc <- as.numeric(lt_co_aic$AICc)
lt_co_aic$LL <- as.numeric(lt_co_aic$LL)
str(lt_co_aic)
lt_co_aic <- lt_co_aic[order(lt_co_aic$AICc), ]
lt_co_aic # The best model is the model with interactive effect of fire frequency and seed weight.
lt_co_aic$Delta_AICc <- "0.00"
lt_co_aic$Delta_AICc[2] <- round(lt_co_aic$AICc[1]-lt_co_aic$AICc[2], 2)
lt_co_aic$Delta_AICc[3] <- round(lt_co_aic$AICc[1]-lt_co_aic$AICc[3], 2)
lt_co_aic$Delta_AICc[4] <- round(lt_co_aic$AICc[1]-lt_co_aic$AICc[4], 2)
lt_co_aic$Delta_AICc[5] <- round(lt_co_aic$AICc[1]-lt_co_aic$AICc[5], 2)
lt_co_aic$Delta_AICc[6] <- round(lt_co_aic$AICc[1]-lt_co_aic$AICc[6], 2)
lt_co_aic$Delta_AICc[7] <- round(lt_co_aic$AICc[1]-lt_co_aic$AICc[7], 2)
lt_co_aic # Model 1 is best with no model ranked within delta AICc <2.






# Predict from models for t50 ----
# GLMER  ----
# Littoralis
m2_lt50
nt50_m2_l <- data.frame(r_seed_wt = seq(min(lit_cum.prop$r_seed_wt), max(lit_cum.prop$r_seed_wt), length = 50)) 
pl_t50 <- predictSE(m2_lt50, newdata = nt50_m2_l, se.fit = T, type = 'link')
nt50_m2_l$fit.link <- pl_t50$fit
nt50_m2_l$se.link <- pl_t50$se.fit
nt50_m2_l$lci.link <- nt50_m2_l$fit.link - (nt50_m2_l$se.link * 1.96)
nt50_m2_l$uci.link <- nt50_m2_l$fit.link + (nt50_m2_l$se.link * 1.96)

# Back-transform - Poisson uses a log link function so the inverse is exp()
nt50_m2_l$fit <- exp(nt50_m2_l$fit.link)
nt50_m2_l$se <- exp(nt50_m2_l$se.link)
nt50_m2_l$lci <- exp(nt50_m2_l$lci.link)
nt50_m2_l$uci <- exp(nt50_m2_l$uci.link)

# Plot predictions
plot(nt50_m2_l$r_seed_wt, nt50_m2_l$fit, type = 'l', xaxt = "n", ylab = "", xlab = "", xlim = c(-2.7, 2.7), ylim = c(10, 13))
axis(side = 1, at = seq(-2.7, 2.7, 0.3), labels = F)
text(seq(-2.7, 2.7, 0.3), par("usr")[3]-0.02, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(10,13, 0.5), label = F)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3.5, cex = 1.5)
pg.ci(x = 'r_seed_wt', data = 'nt50_m2_l', colour = rgb(0,0,0, 0.1), lower = 'lci', upper = 'uci')
mtext(paste("AICc = ", round(AICc(m2_lt50), 3), sep = ""), line = -1.5, cex = 1.2)
# Doesn't look very good because the confidence intervals are really blown out

# Torulosa
nt50_m2_t <- data.frame(r_seed_wt = seq(min(tor_cum.prop$r_seed_wt), max(tor_cum.prop$r_seed_wt), length = 50)) 
pt50_t <- predictSE(m2_tt50, newdata = nt50_m2_t, se.fit = T, type = 'link')
nt50_m2_t <- nt50_m2_t
nt50_m2_t$fit.link <- pt50_t$fit
nt50_m2_t$se.link <- pt50_t$se.fit
nt50_m2_t$lci.link <- nt50_m2_t$fit.link - (nt50_m2_t$se.link * 1.96)
nt50_m2_t$uci.link <- nt50_m2_t$fit.link + (nt50_m2_t$se.link * 1.96)

nt50_m2_t$fit <- exp(nt50_m2_t$fit.link)
nt50_m2_t$se <- exp(nt50_m2_t$se.link)
nt50_m2_t$lci <- exp(nt50_m2_t$lci)
nt50_m2_t$uci <- exp(nt50_m2_t$uci)


nt50_m5_tl <- expand.grid(r_fire_freq = min(tor_cum.prop$r_fire_freq),
                         Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
pt50_m5_tl <- predictSE(m5_tt50, newdata = nt50_m5_tl, se.fit = T, type = 'link')
nt50_m5_tl$fit.link <- pt50_m5_tl$fit
nt50_m5_tl$se.link <- pt50_m5_tl$se.fit
nt50_m5_tl$lci.link <- nt50_m5_tl$fit.link - (nt50_m5_tl$se.link * 1.96)
nt50_m5_tl$uci.link <- nt50_m5_tl$fit.link + (nt50_m5_tl$se.link * 1.96)

nt50_m5_tl$fit <- exp(nt50_m5_tl$fit.link)
nt50_m5_tl$se <- exp(nt50_m5_tl$se.link)
nt50_m5_tl$lci <- exp(nt50_m5_tl$lci.link)
nt50_m5_tl$uci <- exp(nt50_m5_tl$uci.link)

nt50_m5_ta <- expand.grid(r_fire_freq = 3,
                          Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
pt50_m5_ta <- predictSE(m5_tt50, newdata = nt50_m5_ta, se.fit = T, type = 'link')
nt50_m5_ta$fit.link <- pt50_m5_ta$fit
nt50_m5_ta$se.link <- pt50_m5_ta$se.fit
nt50_m5_ta$lci.link <- nt50_m5_ta$fit.link - (nt50_m5_ta$se.link * 1.96)
nt50_m5_ta$uci.link <- nt50_m5_ta$fit.link + (nt50_m5_ta$se.link * 1.96)

nt50_m5_ta$fit <- exp(nt50_m5_ta$fit.link)
nt50_m5_ta$se <- exp(nt50_m5_ta$se.link)
nt50_m5_ta$lci <- exp(nt50_m5_ta$lci.link)
nt50_m5_ta$uci <- exp(nt50_m5_ta$uci.link)


nt50_m5_th <- expand.grid(r_fire_freq = max(tor_cum.prop$r_fire_freq),
                          Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
pt50_m5_th <- predictSE(m5_tt50, newdata = nt50_m5_th, se.fit = T, type = 'link')
nt50_m5_th$fit.link <- pt50_m5_th$fit
nt50_m5_th$se.link <- pt50_m5_th$se.fit
nt50_m5_th$lci.link <- nt50_m5_th$fit.link - (nt50_m5_th$se.link * 1.96)
nt50_m5_th$uci.link <- nt50_m5_th$fit.link + (nt50_m5_th$se.link * 1.96)

nt50_m5_th$fit <- exp(nt50_m5_th$fit.link)
nt50_m5_th$se <- exp(nt50_m5_th$se.link)
nt50_m5_th$lci <- exp(nt50_m5_th$lci.link)
nt50_m5_th$uci <- exp(nt50_m5_th$uci.link)





pt50_m6_t <- predictSE(m6_tt50, newdata = new_tpl, se.fit = T, type = 'link')
nt50_m6_tl <- new_tpl
nt50_m6_tl$fit.link <- pt50_m6_t$fit
nt50_m6_tl$se.link <- pt50_m6_t$se.fit
nt50_m6_tl$lci.link <- nt50_m6_tl$fit.link - (nt50_m6_tl$se.link * 1.96)
nt50_m6_tl$uci.link <- nt50_m6_tl$fit.link + (nt50_m6_tl$se.link * 1.96)

nt50_m6_tl$fit <- exp(nt50_m6_tl$fit.link)
nt50_m6_tl$se <- exp(nt50_m6_tl$se.link)
nt50_m6_tl$lci <- exp(nt50_m6_tl$lci.link)
nt50_m6_tl$uci <- exp(nt50_m6_tl$uci.link)


pt50_m6_ta <- predictSE(m6_tt50, newdata = new_tpa, se.fit = T, type = 'link')
nt50_m6_ta <- new_tpa
nt50_m6_ta$fit.link <- pt50_m6_ta$fit
nt50_m6_ta$se.link <- pt50_m6_ta$se.fit
nt50_m6_ta$lci.link <- nt50_m6_ta$fit.link - (nt50_m6_ta$se.link * 1.96)
nt50_m6_ta$uci.link <- nt50_m6_ta$fit.link + (nt50_m6_ta$se.link * 1.96)

nt50_m6_ta$fit <- exp(nt50_m6_ta$fit.link)
nt50_m6_ta$se <- exp(nt50_m6_ta$se.link)
nt50_m6_ta$lci <- exp(nt50_m6_ta$lci.link)
nt50_m6_ta$uci <- exp(nt50_m6_ta$uci.link)


pt50_m6_th <- predictSE(m6_tt50, newdata = new_tph, se.fit = T, type = 'link')
nt50_m6_th <- new_tph
nt50_m6_th$fit.link <- pt50_m6_th$fit
nt50_m6_th$se.link <- pt50_m6_th$se.fit
nt50_m6_th$lci.link <- nt50_m6_th$fit.link - (nt50_m6_th$se.link * 1.96)
nt50_m6_th$uci.link <- nt50_m6_th$fit.link + (nt50_m6_th$se.link * 1.96)

nt50_m6_th$fit <- exp(nt50_m6_th$fit.link)
nt50_m6_th$se <- exp(nt50_m6_th$se.link)
nt50_m6_th$lci <- exp(nt50_m6_th$lci.link)
nt50_m6_th$uci <- exp(nt50_m6_th$uci.link)


dev.new(width = 20, height = 5, noRStudioGD = T, dpi = 300)
par(mfrow = c(1,3), mar = c(8,6,3,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))

plot(nt50_m2_t$r_seed_wt, nt50_m2_t$fit, type = 'l', xaxt = "n", ylab = "", xlab = "", ylim = c(5,18), xlim = c(-2.2, 2.6), las = 1, cex.axis = 1.4)
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = F)
text(seq(-2.2, 2.6, 0.19), par("usr")[3]-0.6, labels =  seq(2.2, 7.2, 0.2), srt = 60, pos = 1, xpd = T, cex = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
pg.ci(x = "r_seed_wt", data = "nt50_m2_t", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(paste("AICc = ", round(AICc(m2_tt50), 3), sep = ""), line = -1.5, cex = 1.2)




plot.default(nt50_m5_ta$Treatment, nt50_m5_ta$fit, pch = 19, type = 'p', xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(5, 18))
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.6, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_m5_ta$lci, c(1:6), nt50_m5_ta$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(m5_tt50), 3), sep = ""), line = -1.5, cex = 1.2)

points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_m5_tl$fit, pch = 19, col = "steelblue2")
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_m5_tl$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_m5_tl$uci, length = 0.05, code = 3, angle = 90, col = 'steelblue2')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_m5_th$fit, pch = 19, col = 'red')
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_m5_th$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_m5_th$uci, length = 0.05, code = 3, angle = 90, col = 'red')




plot(nt50_m6_ta$r_seed_wt, nt50_m6_ta$fit, type = 'l', xaxt = "n", xlab = "", ylab = "", las = 1, col = "black", cex.axis = 1.4, ylim = c(5,18))
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = F)
text(seq(-2.2, 2.6, 0.19), par("usr")[3]-0.6, labels =  seq(2.2, 7.2, 0.2), srt = 60, pos = 1, xpd = T, cex = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
pg.ci(x = "r_seed_wt", data = "nt50_m6_ta", colour = rgb(0/255, 0/255, 0/255, 0.05), lower = "lci", upper = "uci")

lines(nt50_m6_th$r_seed_wt, nt50_m6_th$fit, type = 'l', col = 'red')
pg.ci(x = "r_seed_wt", data = "nt50_m6_th", colour = rgb(1,0,0, 0.1), lower = "lci", upper = "uci")
lines(nt50_m6_tl$r_seed_wt, nt50_m6_tl$fit, type = 'l', col = 'steelblue2', lwd = 1.2)
pg.ci(x = 'r_seed_wt', data = 'nt50_m6_tl', colour = rgb(92/255, 172/255, 238/255, 0.1), lower = 'lci', upper = 'uci')
mtext(paste("AICc = ", round(AICc(m6_tt50), 3), sep = ""), line = -1.5, cex = 1.2)


par(xpd = NA)
legend(x = 2.8, y = 18, legend = c("0 fires", "3 fires", "6 fires"), col = c("steelblue2", "black", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)

# For GLMER focus on the best model only for torulosa. We can still make some comparison for torulosa between GLMER and GAM but the GAM just includes treatment with seed weight. We cannot make any comparisons between the GLMER and GAM for littoralis as GAMs selected the model with treatment only as the best model.



# GAMs ----
# Littoralis - model 1
lt1_s
nt50_m1_l <- data.frame(Treatment = as.factor(c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke")))
p50_m1_ls <- predict(lt1_s, newdata = nt50_m1_l, se.fit = T, type = 'response')
nt50_m1_l$fit <- p50_m1_ls$fit
nt50_m1_l$se <- p50_m1_ls$se.fit
nt50_m1_l$lci <- nt50_m1_l$fit - (nt50_m1_l$se * 1.96) 
nt50_m1_l$uci <- nt50_m1_l$fit + (nt50_m1_l$se * 1.96) 


# Plot
dev.new(width = 15, height = 8, noRStudioGD = T, dpi = 300)
par(mfrow = c(1,2), mar = c(6,6,3,2))

plot(nt50_m2_l$r_seed_wt, nt50_m2_l$fit, type = 'l', xaxt = "n", ylab = "", xlab = "", xlim = c(-2.7, 2.7), ylim = c(10, 13), las = 1, cex.axis = 1.4)
axis(side = 1, at = seq(-2.7, 2.7, 0.3), labels = F)
text(seq(-2.7, 2.7, 0.3), par("usr")[3]-0.1, labels = seq(1.3, 3.1, 0.1), srt = 60, pos = 1, xpd = T, cex = 1.4)
axis(side = 2, at = seq(10,13, 0.5), label = F)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3.5, cex = 1.5)
pg.ci(x = 'r_seed_wt', data = 'nt50_m2_l', colour = rgb(0,0,0, 0.1), lower = 'lci', upper = 'uci')
mtext(paste("AICc = ", round(AICc(m2_lt50), 3), sep = ""), line = -1.5, cex = 1.2)
mtext(expression(bold("(a) GLMER")), cex = 2)

plot.default(nt50_m1_l$Treatment, nt50_m1_l$fit, pch = 19, type = 'p', xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(9, 15))
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
text(x = c(1:6), y = par("usr")[3] - 0.6, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 5, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_m1_l$lci, c(1:6), nt50_m1_l$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(lt1_s), 3), sep = ""), line = -1.5, cex = 1.2)
mtext(expression(bold("(b) GAM")), cex = 2)


# Torulosa
nt50_tl <- expand.grid(seed_wt_mg = c(min(tor_cum.prop$seed_wt_mg)), 
                       Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
nt50_ta <- expand.grid(seed_wt_mg = c(mean(tor_cum.prop$seed_wt_mg)), 
                       Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
nt50_th <- expand.grid(seed_wt_mg = c(max(tor_cum.prop$seed_wt_mg)), 
                       Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))

pt50_tls <- predict(tt4_s, newdata = nt50_tl, se.fit = T, type = 'response')
nt50_tls <- nt50_tl
nt50_tls$fit <- pt50_tls$fit
nt50_tls$se <- pt50_tls$se.fit
nt50_tls$lci <- nt50_tls$fit - (nt50_tls$se * 1.96)
nt50_tls$uci <- nt50_tls$fit + (nt50_tls$se * 1.96)

pt50_tas <- predict(tt4_s, newdata = nt50_ta, se.fit = T, type = 'response')
nt50_tas <- nt50_ta
nt50_tas$fit <- pt50_tas$fit
nt50_tas$se <- pt50_tas$se.fit
nt50_tas$lci <- nt50_tas$fit - (nt50_tas$se * 1.96)
nt50_tas$uci <- nt50_tas$fit + (nt50_tas$se * 1.96)

pt50_ths <- predict(tt4_s, newdata = nt50_th, se.fit = T, type = 'response')
nt50_ths <- nt50_th
nt50_ths$fit <- pt50_ths$fit
nt50_ths$se <- pt50_ths$se.fit
nt50_ths$lci <- nt50_ths$fit - (nt50_ths$se * 1.96)
nt50_ths$uci <- nt50_ths$fit + (nt50_ths$se * 1.96)



pt50_tlso <- predict(tt4_so, newdata = nt50_tl, se.fit = T, type = 'response')
nt50_tlso <- nt50_tl
nt50_tlso$fit <- pt50_tlso$fit
nt50_tlso$se <- pt50_tlso$se.fit
nt50_tlso$lci <- nt50_tlso$fit - (nt50_tlso$se * 1.96)
nt50_tlso$uci <- nt50_tlso$fit + (nt50_tlso$se * 1.96)

pt50_taso <- predict(tt4_so, newdata = nt50_ta, se.fit = T, type = 'response')
nt50_taso <- nt50_ta
nt50_taso$fit <- pt50_taso$fit
nt50_taso$se <- pt50_taso$se.fit
nt50_taso$lci <- nt50_taso$fit - (nt50_taso$se * 1.96)
nt50_taso$uci <- nt50_taso$fit + (nt50_taso$se * 1.96)

pt50_thso <- predict(tt4_so, newdata = nt50_th, se.fit = T, type = 'response')
nt50_thso <- nt50_th
nt50_thso$fit <- pt50_thso$fit
nt50_thso$se <- pt50_thso$se.fit
nt50_thso$lci <- nt50_thso$fit - (nt50_thso$se * 1.96)
nt50_thso$uci <- nt50_thso$fit + (nt50_thso$se * 1.96)






pt50_tlsc <- predict(tt4_sc, newdata = nt50_tl, se.fit = T, type = 'response')
nt50_tlsc <- nt50_tl
nt50_tlsc$fit <- pt50_tlsc$fit
nt50_tlsc$se <- pt50_tlsc$se.fit
nt50_tlsc$lci <- nt50_tlsc$fit - (nt50_tlsc$se * 1.96)
nt50_tlsc$uci <- nt50_tlsc$fit + (nt50_tlsc$se * 1.96)

pt50_tasc <- predict(tt4_sc, newdata = nt50_ta, se.fit = T, type = 'response')
nt50_tasc <- nt50_ta
nt50_tasc$fit <- pt50_tasc$fit
nt50_tasc$se <- pt50_tasc$se.fit
nt50_tasc$lci <- nt50_tasc$fit - (nt50_tasc$se * 1.96)
nt50_tasc$uci <- nt50_tasc$fit + (nt50_tasc$se * 1.96)

pt50_thsc <- predict(tt4_sc, newdata = nt50_th, se.fit = T, type = 'response')
nt50_thsc <- nt50_th
nt50_thsc$fit <- pt50_thsc$fit
nt50_thsc$se <- pt50_thsc$se.fit
nt50_thsc$lci <- nt50_thsc$fit - (nt50_thsc$se * 1.96)
nt50_thsc$uci <- nt50_thsc$fit + (nt50_thsc$se * 1.96)





pt50_tlsco <- predict(tt4_sco, newdata = nt50_tl, se.fit = T, type = 'response')
nt50_tlsco <- nt50_tl
nt50_tlsco$fit <- pt50_tlsco$fit
nt50_tlsco$se <- pt50_tlsco$se.fit
nt50_tlsco$lci <- nt50_tlsco$fit - (nt50_tlsco$se * 1.96)
nt50_tlsco$uci <- nt50_tlsco$fit + (nt50_tlsco$se * 1.96)


pt50_tasco <- predict(tt4_sco, newdata = nt50_ta, se.fit = T, type = 'response')
nt50_tasco <- nt50_ta
nt50_tasco$fit <- pt50_tasco$fit
nt50_tasco$se <- pt50_tasco$se.fit
nt50_tasco$lci <- nt50_tasco$fit - (nt50_tasco$se * 1.96)
nt50_tasco$uci <- nt50_tasco$fit + (nt50_tasco$se * 1.96)


pt50_thsco <- predict(tt4_sco, newdata = nt50_th, se.fit = T, type = 'response')
nt50_thsco <- nt50_th
nt50_thsco$fit <- pt50_thsco$fit
nt50_thsco$se <- pt50_thsco$se.fit
nt50_thsco$lci <- nt50_thsco$fit - (nt50_thsco$se * 1.96)
nt50_thsco$uci <- nt50_thsco$fit + (nt50_thsco$se * 1.96)



# Plots
dev.new(width = 20, height = 12, noRStudioGD = T, dpi = 300)
par(mfrow = c(2,3), mar = c(8,6,3,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))

plot(nt50_m2_t$r_seed_wt, nt50_m2_t$fit, type = 'l', xaxt = "n", ylab = "", xlab = "", ylim = c(5,23), xlim = c(-2.2, 2.6), las = 1, cex.axis = 1.4, yaxt = "n")
axis(side = 1, at = seq(-2.2, 2.6, 0.19), labels = F)
axis(side = 2, at = seq(5,23,2), labels = seq(5,23,2), cex.axis = 1.4, las = 1)
axis(side = 2, at = seq(5,23,1), labels = F)
text(seq(-2.2, 2.6, 0.19), par("usr")[3]-0.6, labels =  seq(2.2, 7.2, 0.2), srt = 60, pos = 1, xpd = T, cex = 1.4)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
pg.ci(x = "r_seed_wt", data = "nt50_m2_t", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext(paste("AICc = ", round(AICc(m2_tt50), 3), sep = ""), line = -1.5, cex = 1.2)
mtext(expression(bold("(a) GLMER - seed weight only")), cex = 2)



plot.default(nt50_tas$Treatment, nt50_tas$fit, type = 'p', pch = 19, xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(5, 23), yaxt = "n")
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
axis(side = 2, at = seq(5,23,2), labels = seq(5,23,2), cex.axis = 1.4, las = 1)
axis(side = 2, at = seq(5,23,1), labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.4, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_tas$lci, c(1:6), nt50_tas$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(tt4_s), 3), sep = ""), line = -1.5, cex = 1.2)
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tls$fit, pch = 19, col = "steelblue2")
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tls$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tls$uci, length = 0.05, code = 3, angle = 90, col = 'steelblue2')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_ths$fit, pch = 19, col = 'red')
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_ths$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_ths$uci, length = 0.05, code = 3, angle = 90, col = 'red')
mtext(expression(bold("(b) GAM - default s(bs = 'tp')")), cex = 2)




plot.default(nt50_taso$Treatment, nt50_taso$fit, type = 'p', pch = 19, xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(5, 23), yaxt = "n")
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
axis(side = 2, at = seq(5,23,2), labels = seq(5,23,2), cex.axis = 1.4, las = 1)
axis(side = 2, at = seq(5,23,1), labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.4, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_taso$lci, c(1:6), nt50_taso$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(tt4_so), 3), sep = ""), line = -1.5, cex = 1.2)
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlso$fit, pch = 19, col = "steelblue2")
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlso$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlso$uci, length = 0.05, code = 3, angle = 90, col = 'steelblue2')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thso$fit, pch = 19, col = 'red')
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thso$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thso$uci, length = 0.05, code = 3, angle = 90, col = 'red')
mtext(expression(bold("(c) GAM - s() optimised")), cex = 2)

par(xpd = NA)
legend(x = 6.2, y = 23, legend = c("0 fires", "3 fires", "6 fires"), col = c("steelblue2", "black", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)


plot.default(nt50_tasc$Treatment, nt50_tasc$fit, type = 'p', pch = 19, xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(5, 23), yaxt = "n")
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
axis(side = 2, at = seq(5,23,2), labels = seq(5,23,2), cex.axis = 1.4, las = 1)
axis(side = 2, at = seq(5,23,1), labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.4, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_tasc$lci, c(1:6), nt50_tasc$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(tt4_sc), 3), sep = ""), line = -1.5, cex = 1.2)
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlsc$fit, pch = 19, col = "steelblue2")
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlsc$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlsc$uci, length = 0.05, code = 3, angle = 90, col = 'steelblue2')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thsc$fit, pch = 19, col = 'red')
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thsc$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thsc$uci, length = 0.05, code = 3, angle = 90, col = 'red')
mtext(expression(bold("(d) GAM - ti(bs = 'cr') default")), cex = 2)



plot.default(nt50_tasco$Treatment, nt50_tasco$fit, type = 'p', pch = 19, xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(5, 23), yaxt = "n")
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
axis(side = 2, at = seq(5,23,2), labels = seq(5,23,2), cex.axis = 1.4, las = 1)
axis(side = 2, at = seq(5,23,1), labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.4, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_tasco$lci, c(1:6), nt50_tasco$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(tt4_sco), 3), sep = ""), line = -1.5, cex = 1.2)
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlsco$fit, pch = 19, col = "steelblue2")
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlsco$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlsco$uci, length = 0.05, code = 3, angle = 90, col = 'steelblue2')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thsco$fit, pch = 19, col = 'red')
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thsco$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thsco$uci, length = 0.05, code = 3, angle = 90, col = 'red')
mtext(expression(bold("(e) GAM - ti(bs = 'cr') optimised")), cex = 2)


# Test some other basis functions for torulosa model 4. As littoralis best model is model 1, treatment is not wrapped in smooth terms so no further optimisatoin is available
tt4_ts <- gam(t50 ~ s(Treatment, bs = 're') + s(seed_wt_mg, bs = 'ts', k = 5) + ti(seed_wt_mg, by = Treatment, bs = 'ts', k = 7), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
par(mfrow = c(2,2)); gam.check(tt4_ts)
par(mfrow = c(3,3));plot(tt4_ts)

tt4_ps <- gam(t50 ~ s(Treatment, bs = 're') + s(seed_wt_mg, bs = 'ps', k = 6) + ti(seed_wt_mg, by = Treatment, bs = 'ps', k = 4), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
par(mfrow = c(2,2)); gam.check(tt4_ps)
par(mfrow = c(3,3));plot(tt4_ps)


tt4_bs <- gam(t50 ~ s(Treatment, bs = 're') + s(seed_wt_mg, bs = 'bs', k = 7) + ti(seed_wt_mg, by = Treatment, bs = 'bs'), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
par(mfrow = c(2,2)); gam.check(tt4_bs)
par(mfrow = c(3,3));plot(tt4_bs)

tt4_ad <- gam(t50 ~ s(Treatment, bs = 're') + s(seed_wt_mg, bs = 'ad', k = 15) + ti(seed_wt_mg, by = Treatment, bs = 'tp', k = 4), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
par(mfrow = c(2,2)); gam.check(tt4_ad)
par(mfrow = c(3,3));plot(tt4_ad)

tt4_adc <- gam(t50 ~ s(Treatment, bs = 're') + s(seed_wt_mg, bs = 'ad', k = 15) + ti(seed_wt_mg, by = Treatment, bs = 'cr', k = 5), family = poisson, data = tor_cum.prop, random = ~ (1|Individual) + (1|Set), method = 'ML')
par(mfrow = c(2,2)); gam.check(tt4_adc)
par(mfrow = c(3,3));plot(tt4_adc)


# Produce predictions from these new models
# ts
pt50_tl_ts <- predict(tt4_s, newdata = nt50_tl, se.fit = T, type = 'response')
nt50_tl_ts <- nt50_tl
nt50_tl_ts$fit <- pt50_tl_ts$fit
nt50_tl_ts$se <- pt50_tl_ts$se.fit
nt50_tl_ts$lci <- nt50_tl_ts$fit - (nt50_tl_ts$se * 1.96)
nt50_tl_ts$uci <- nt50_tl_ts$fit + (nt50_tl_ts$se * 1.96)

pt50_ta_ts <- predict(tt4_s, newdata = nt50_ta, se.fit = T, type = 'response')
nt50_ta_ts <- nt50_ta
nt50_ta_ts$fit <- pt50_ta_ts$fit
nt50_ta_ts$se <- pt50_ta_ts$se.fit
nt50_ta_ts$lci <- nt50_ta_ts$fit - (nt50_ta_ts$se * 1.96)
nt50_ta_ts$uci <- nt50_ta_ts$fit + (nt50_ta_ts$se * 1.96)

pt50_th_ts <- predict(tt4_s, newdata = nt50_th, se.fit = T, type = 'response')
nt50_th_ts <- nt50_th
nt50_th_ts$fit <- pt50_th_ts$fit
nt50_th_ts$se <- pt50_th_ts$se.fit
nt50_th_ts$lci <- nt50_th_ts$fit - (nt50_th_ts$se * 1.96)
nt50_th_ts$uci <- nt50_th_ts$fit + (nt50_th_ts$se * 1.96)



#ps
pt50_tl_ps <- predict(tt4_s, newdata = nt50_tl, se.fit = T, type = 'response')
nt50_tl_ps <- nt50_tl
nt50_tl_ps$fit <- pt50_tl_ps$fit
nt50_tl_ps$se <- pt50_tl_ps$se.fit
nt50_tl_ps$lci <- nt50_tl_ps$fit - (nt50_tl_ps$se * 1.96)
nt50_tl_ps$uci <- nt50_tl_ps$fit + (nt50_tl_ps$se * 1.96)

pt50_ta_ps <- predict(tt4_s, newdata = nt50_ta, se.fit = T, type = 'response')
nt50_ta_ps <- nt50_ta
nt50_ta_ps$fit <- pt50_ta_ps$fit
nt50_ta_ps$se <- pt50_ta_ps$se.fit
nt50_ta_ps$lci <- nt50_ta_ps$fit - (nt50_ta_ps$se * 1.96)
nt50_ta_ps$uci <- nt50_ta_ps$fit + (nt50_ta_ps$se * 1.96)

pt50_th_ps <- predict(tt4_s, newdata = nt50_th, se.fit = T, type = 'response')
nt50_th_ps <- nt50_th
nt50_th_ps$fit <- pt50_th_ps$fit
nt50_th_ps$se <- pt50_th_ps$se.fit
nt50_th_ps$lci <- nt50_th_ps$fit - (nt50_th_ps$se * 1.96)
nt50_th_ps$uci <- nt50_th_ps$fit + (nt50_th_ps$se * 1.96)



# bs
pt50_tl_bs <- predict(tt4_s, newdata = nt50_tl, se.fit = T, type = 'response')
nt50_tl_bs <- nt50_tl
nt50_tl_bs$fit <- pt50_tl_bs$fit
nt50_tl_bs$se <- pt50_tl_bs$se.fit
nt50_tl_bs$lci <- nt50_tl_bs$fit - (nt50_tl_bs$se * 1.96)
nt50_tl_bs$uci <- nt50_tl_bs$fit + (nt50_tl_bs$se * 1.96)

pt50_ta_bs <- predict(tt4_s, newdata = nt50_ta, se.fit = T, type = 'response')
nt50_ta_bs <- nt50_ta
nt50_ta_bs$fit <- pt50_ta_bs$fit
nt50_ta_bs$se <- pt50_ta_bs$se.fit
nt50_ta_bs$lci <- nt50_ta_bs$fit - (nt50_ta_bs$se * 1.96)
nt50_ta_bs$uci <- nt50_ta_bs$fit + (nt50_ta_bs$se * 1.96)

pt50_th_bs <- predict(tt4_s, newdata = nt50_th, se.fit = T, type = 'response')
nt50_th_bs <- nt50_th
nt50_th_bs$fit <- pt50_th_bs$fit
nt50_th_bs$se <- pt50_th_bs$se.fit
nt50_th_bs$lci <- nt50_th_bs$fit - (nt50_th_bs$se * 1.96)
nt50_th_bs$uci <- nt50_th_bs$fit + (nt50_th_bs$se * 1.96)


# ad
pt50_tl_ad <- predict(tt4_s, newdata = nt50_tl, se.fit = T, type = 'response')
nt50_tl_ad <- nt50_tl
nt50_tl_ad$fit <- pt50_tl_ad$fit
nt50_tl_ad$se <- pt50_tl_ad$se.fit
nt50_tl_ad$lci <- nt50_tl_ad$fit - (nt50_tl_ad$se * 1.96)
nt50_tl_ad$uci <- nt50_tl_ad$fit + (nt50_tl_ad$se * 1.96)

pt50_ta_ad <- predict(tt4_s, newdata = nt50_ta, se.fit = T, type = 'response')
nt50_ta_ad <- nt50_ta
nt50_ta_ad$fit <- pt50_ta_ad$fit
nt50_ta_ad$se <- pt50_ta_ad$se.fit
nt50_ta_ad$lci <- nt50_ta_ad$fit - (nt50_ta_ad$se * 1.96)
nt50_ta_ad$uci <- nt50_ta_ad$fit + (nt50_ta_ad$se * 1.96)

pt50_th_ad <- predict(tt4_s, newdata = nt50_th, se.fit = T, type = 'response')
nt50_th_ad <- nt50_th
nt50_th_ad$fit <- pt50_th_ad$fit
nt50_th_ad$se <- pt50_th_ad$se.fit
nt50_th_ad$lci <- nt50_th_ad$fit - (nt50_th_ad$se * 1.96)
nt50_th_ad$uci <- nt50_th_ad$fit + (nt50_th_ad$se * 1.96)



# adc
pt50_tl_adc <- predict(tt4_s, newdata = nt50_tl, se.fit = T, type = 'response')
nt50_tl_adc <- nt50_tl
nt50_tl_adc$fit <- pt50_tl_adc$fit
nt50_tl_adc$se <- pt50_tl_adc$se.fit
nt50_tl_adc$lci <- nt50_tl_adc$fit - (nt50_tl_adc$se * 1.96)
nt50_tl_adc$uci <- nt50_tl_adc$fit + (nt50_tl_adc$se * 1.96)

pt50_ta_adc <- predict(tt4_s, newdata = nt50_ta, se.fit = T, type = 'response')
nt50_ta_adc <- nt50_ta
nt50_ta_adc$fit <- pt50_ta_adc$fit
nt50_ta_adc$se <- pt50_ta_adc$se.fit
nt50_ta_adc$lci <- nt50_ta_adc$fit - (nt50_ta_adc$se * 1.96)
nt50_ta_adc$uci <- nt50_ta_adc$fit + (nt50_ta_adc$se * 1.96)

pt50_th_adc <- predict(tt4_s, newdata = nt50_th, se.fit = T, type = 'response')
nt50_th_adc <- nt50_th
nt50_th_adc$fit <- pt50_th_adc$fit
nt50_th_adc$se <- pt50_th_adc$se.fit
nt50_th_adc$lci <- nt50_th_adc$fit - (nt50_th_adc$se * 1.96)
nt50_th_adc$uci <- nt50_th_adc$fit + (nt50_th_adc$se * 1.96)




# Plots
dev.new(width = 20, height = 12, noRStudioGD = T, dpi = 300)
par(mfrow = c(2,3), mar = c(8,6,3,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))


plot.default(nt50_taso$Treatment, nt50_taso$fit, type = 'p', pch = 19, xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(5, 23), yaxt = "n")
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
axis(side = 2, at = seq(5,23,2), labels = seq(5,23,2), cex.axis = 1.4, las = 1)
axis(side = 2, at = seq(5,23,1), labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.4, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_taso$lci, c(1:6), nt50_taso$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(tt4_so), 3), sep = ""), line = -1.5, cex = 1.2)
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlso$fit, pch = 19, col = "steelblue2")
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlso$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlso$uci, length = 0.05, code = 3, angle = 90, col = 'steelblue2')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thso$fit, pch = 19, col = 'red')
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thso$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thso$uci, length = 0.05, code = 3, angle = 90, col = 'red')
mtext(expression(bold("(a) GAM - s(bs = 'tp') optimised")), cex = 2)



plot.default(nt50_tasco$Treatment, nt50_tasco$fit, type = 'p', pch = 19, xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(5, 23), yaxt = "n")
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
axis(side = 2, at = seq(5,23,2), labels = seq(5,23,2), cex.axis = 1.4, las = 1)
axis(side = 2, at = seq(5,23,1), labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.4, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_tasco$lci, c(1:6), nt50_tasco$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(tt4_sco), 3), sep = ""), line = -1.5, cex = 1.2)
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlsco$fit, pch = 19, col = "steelblue2")
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlsco$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tlsco$uci, length = 0.05, code = 3, angle = 90, col = 'steelblue2')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thsco$fit, pch = 19, col = 'red')
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thsco$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_thsco$uci, length = 0.05, code = 3, angle = 90, col = 'red')
mtext(expression(bold("(b) GAM - ti(bs = 'cr') optimised")), cex = 2)



plot.default(nt50_ta_ts$Treatment, nt50_ta_ts$fit, type = 'p', pch = 19, xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(5, 23), yaxt = "n")
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
axis(side = 2, at = seq(5,23,2), labels = seq(5,23,2), cex.axis = 1.4, las = 1)
axis(side = 2, at = seq(5,23,1), labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.4, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_ta_ts$lci, c(1:6), nt50_ta_ts$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(tt4_ts), 3), sep = ""), line = -1.5, cex = 1.2)
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_ts$fit, pch = 19, col = "steelblue2")
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_ts$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_ts$uci, length = 0.05, code = 3, angle = 90, col = 'steelblue2')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_ts$fit, pch = 19, col = 'red')
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_ts$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_ts$uci, length = 0.05, code = 3, angle = 90, col = 'red')
mtext(expression(bold("(c) GAM - bs = 'ts'")), cex = 2)

par(xpd = NA)
legend(x = 6.2, y = 23, legend = c("0 fires", "3 fires", "6 fires"), col = c("steelblue2", "black", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)



plot.default(nt50_ta_ps$Treatment, nt50_ta_ps$fit, type = 'p', pch = 19, xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(5, 23), yaxt = "n")
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
axis(side = 2, at = seq(5,23,2), labels = seq(5,23,2), cex.axis = 1.4, las = 1)
axis(side = 2, at = seq(5,23,1), labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.4, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_ta_ps$lci, c(1:6), nt50_ta_ps$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(tt4_ps), 3), sep = ""), line = -1.5, cex = 1.2)
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_ps$fit, pch = 19, col = "steelblue2")
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_ps$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_ps$uci, length = 0.05, code = 3, angle = 90, col = 'steelblue2')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_ps$fit, pch = 19, col = 'red')
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_ps$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_ps$uci, length = 0.05, code = 3, angle = 90, col = 'red')
mtext(expression(bold("(d) GAM - bs = 'ps'")), cex = 2)


plot.default(nt50_ta_bs$Treatment, nt50_ta_bs$fit, type = 'p', pch = 19, xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(5, 23), yaxt = "n")
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
axis(side = 2, at = seq(5,23,2), labels = seq(5,23,2), cex.axis = 1.4, las = 1)
axis(side = 2, at = seq(5,23,1), labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.4, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_ta_bs$lci, c(1:6), nt50_ta_bs$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(tt4_bs), 3), sep = ""), line = -1.5, cex = 1.2)
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_bs$fit, pch = 19, col = "steelblue2")
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_bs$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_bs$uci, length = 0.05, code = 3, angle = 90, col = 'steelblue2')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_bs$fit, pch = 19, col = 'red')
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_bs$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_bs$uci, length = 0.05, code = 3, angle = 90, col = 'red')
mtext(expression(bold("(e) GAM - bs = 'bs'")), cex = 2)



plot.default(nt50_ta_adc$Treatment, nt50_ta_adc$fit, type = 'p', pch = 19, xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 1.4, ylim = c(5, 23), yaxt = "n")
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
axis(side = 2, at = seq(5,23,2), labels = seq(5,23,2), cex.axis = 1.4, las = 1)
axis(side = 2, at = seq(5,23,1), labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.4, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(side = 1, expression(bold("Treatment")), line = 6, cex = 1.5)
mtext(side = 2, expression(bold("Time to 50% germination (days)")), line = 3, cex = 1.5)
arrows(c(1:6), nt50_ta_adc$lci, c(1:6), nt50_ta_adc$uci, length = 0.05, code = 3, angle = 90)
mtext(paste("AICc = ", round(AICc(tt4_adc), 3), sep = ""), line = -1.5, cex = 1.2)
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_adc$fit, pch = 19, col = "steelblue2")
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_adc$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_tl_adc$uci, length = 0.05, code = 3, angle = 90, col = 'steelblue2')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_adc$fit, pch = 19, col = 'red')
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_adc$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_th_adc$uci, length = 0.05, code = 3, angle = 90, col = 'red')
mtext(expression(bold("(f) GAM - bs = 'ad' and ti('cr')")), cex = 2)





# 5. QUESTION 2 :How does recent fire activity influence population age structure and female fecundity? ----
# Proportions of seedling, saplings, recruits and number of cones as response
# TSF
# TSF * latitude
# TSF * FPC
# TSF * precipitation seasonality
# TSF * temperature seasonality
  
  
### For fecundity only - Maybe we are using the tree level data, then we nest individualin transect in location? Otherwise we can't really investigate height 
# TSF + height 
# TSF + height * latitude
# TSF * latitude + height
# TSF * FPC + height
# TSF + height * FPC
# TSF + height * precipitation seasonality
# TSF * precipitation seasonality + height

# Rescale variables for GLMER
tor_transects$r_TSF <- scale(tor_transects$TSF)
lit_transects$r_TSF <- scale(lit_transects$TSF)
tor_transects$r_Latitude <- scale(tor_transects$Latitude)
lit_transects$r_Latitude <- scale(lit_transects$Latitude)
tor_transects$r_FPC <- scale(tor_transects$FPC)
lit_transects$r_FPC <- scale(lit_transects$FPC)
tor_transects$r_Precip <- scale(tor_transects$Precip)
lit_transects$r_Precip <- scale(lit_transects$Precip)
tor_transects$r_Temp <- scale(tor_transects$Temp)
lit_transects$r_Temp <- scale(lit_transects$Temp)

# 5.1 Population structure ----
# 5.1.1 Proportion seedlings ---- 
rect_null <- glmer(Proportion_seedlings ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
rect_gnull <- gam(Proportion_seedlings ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')

recl_null <- glmer(Proportion_seedlings ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)
recl_gnull <- gam(Proportion_seedlings ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')



rec_m1t <- glmer(Proportion_seedlings ~ r_TSF + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(rec_m1t)
rec_g1t <- gam(Proportion_seedlings ~ s(TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g1t)
plot(rec_g1t)
rec_g1.1t <- gam(Proportion_seedlings ~ s(TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g1.1t)
plot(rec_g1.1t)

rec_m1l <- glmer(Proportion_seedlings ~ r_TSF + (1 | Location/Transect), family = binomial, data = lit_transects)
summary(rec_m1l)
rec_g1l <- gam(Proportion_seedlings ~ s(TSF, k = 4), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g1l)
plot(rec_g1l)
rec_g1.1l <- gam(Proportion_seedlings ~ s(TSF, bs = 'cr', k = 4), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g1.1l)
plot(rec_g1.1l)


rec_m2t <- glmer(Proportion_seedlings ~ r_TSF * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rec_m2t)
rec_g2t <- gam(Proportion_seedlings ~ s(TSF) + s(Latitude) + ti(Latitude, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g2t)
plot(rec_g2t)
summary(rec_g2t)
rec_g2.1t <- gam(Proportion_seedlings ~ s(TSF) + s(Latitude) + ti(Latitude, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g2.1t)
plot(rec_g2.1t)
rec_g2.2t <- gam(Proportion_seedlings ~ s(TSF, bs = 'cr', k = 5) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g2.1t)
plot(rec_g2.2t)


#rec_m2l <- glmer(Proportion_seedlings ~ TSF * Latitude + (1|Location/Transect), family = binomial, data = lit_transects) # This model fails due to the low number of points available to model an interaction. We can only model an additive effect of TSF and latitude
rec_m2.1l <- glmer(Proportion_seedlings ~ r_TSF + r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rec_m2.1l)
rec_g2l <- gam(Proportion_seedlings ~ s(TSF, k = 4) + s(Latitude, k = 4) + ti(Latitude, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g2l)
plot(rec_g2l)
rec_g2.1l <- gam(Proportion_seedlings ~ s(TSF, k = 4) + s(Latitude, k = 4) + ti(Latitude, by = TSF, bs = 'tp', k = 8), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g2.1l)
plot(rec_g2.1l)
rec_g2.2l <- gam(Proportion_seedlings ~ s(TSF, bs = 'cr', k = 4) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g2.1l)
plot(rec_g2.2l)



rec_m3t <- glmer(Proportion_seedlings ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
rec_g3t <- gam(Proportion_seedlings ~ s(TSF) + s(FPC) + ti(FPC, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g3t)
plot(rec_g3t)
rec_g3.1t <- gam(Proportion_seedlings ~ s(TSF) + s(FPC) + ti(FPC, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g3.1t)
plot(rec_g3.1t)
rec_g3.2t <- gam(Proportion_seedlings ~ s(TSF, bs = 'cr', k = 5) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g3.1t)
plot(rec_g3.2t)


rec_m3l <- glmer(Proportion_seedlings ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rec_m3l)
rec_g3l <- gam(Proportion_seedlings ~ s(TSF, k = 4) + s(FPC, k = 4) + ti(FPC, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g3l)
plot(rec_g3l)
rec_g3.1l <- gam(Proportion_seedlings ~ s(TSF, k = 4) + s(FPC, k = 4) + ti(FPC, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g4.1l)
plot(rec_g4.1l)
rec_g3.2l <- gam(Proportion_seedlings ~ s(TSF, bs = 'cr', k = 4) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g3.2l)
plot(rec_g3.2l)


rec_m4t <- glmer(Proportion_seedlings ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
rec_g4t <- gam(Proportion_seedlings ~ s(TSF) + s(Precip) + ti(Precip, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g4t)
plot(rec_g4t)
rec_g4.1t <- gam(Proportion_seedlings ~ s(TSF) + s(Precip) + ti(Precip, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g4.1t)
plot(rec_g4.1t)
rec_g4.2t <- gam(Proportion_seedlings ~ s(TSF, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g4.1t)
plot(rec_g4.2t)


rec_m4l <- glmer(Proportion_seedlings ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rec_m4l)
rec_g4l <- gam(Proportion_seedlings ~ s(TSF, k = 4) + s(Precip, k = 4) + ti(Precip, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g4l)
plot(rec_g4l)
rec_g4.1l <- gam(Proportion_seedlings ~ s(TSF, k = 4) + s(Precip, k = 4) + ti(Precip, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g4.1l)
plot(rec_g4.1l)
rec_g4.2l <- gam(Proportion_seedlings ~ s(TSF, bs = 'cr', k = 4) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g4.2l)
plot(rec_g4.2l)

rec_m5t <- glmer(Proportion_seedlings ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
rec_g5t <- gam(Proportion_seedlings ~ s(TSF) + s(Temp) + ti(Temp, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g5t)
plot(rec_g5t)
rec_g5.1t <- gam(Proportion_seedlings ~ s(TSF) + s(Temp) + ti(Temp, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g5.1t)
plot(rec_g5.1t)
rec_g5.2t <- gam(Proportion_seedlings ~ s(TSF, bs = 'cr', k = 5) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g5.1t)
plot(rec_g5.2t)


rec_m5l <- glmer(Proportion_seedlings ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=40)))
summary(rec_m5l)
rec_g5l <- gam(Proportion_seedlings ~ s(TSF, k = 4) + s(Temp, k = 4) + ti(Temp, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g5l)
plot(rec_g5l)
rec_g5.1l <- gam(Proportion_seedlings ~ s(TSF, k = 4) + s(Temp, k = 4) + ti(Temp, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g5.1l)
plot(rec_g5.1l)
rec_g5.2l <- gam(Proportion_seedlings ~ s(TSF, bs = 'cr', k = 4) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rec_g5.2l)
plot(rec_g5.2l)


# Compare AICs for each modelling type
# GLMER
rec_aic_t <- list(rect_null, rec_m1t, rec_m2t, rec_m3t, rec_m4t, rec_m5t)
aictab(rec_aic_t) # Null model is best

rec_aic_l <- list(recl_null, rec_m1l, rec_m2.1l, rec_m3l, rec_m4l, rec_m5l)
aictab(rec_aic_l) # Null model is best

# GAM defaults
rec_td_aic <- as.data.frame(1:6)
rec_td_aic$AICc <- "NA"
rec_td_aic$Model <- "NA"
rec_td_aic$LL <- "NA"
rec_td_aic$AICc[1] <- AICc(rect_gnull)
rec_td_aic$Model[1] <- "Null"
rec_td_aic$LL[1] <- logLik(rect_gnull)
rec_td_aic$AICc[2] <- AICc(rec_g1t)
rec_td_aic$Model[2] <- "m1"
rec_td_aic$LL[2] <- logLik(rec_g1t)
rec_td_aic$AICc[3] <- AICc(rec_g2t)
rec_td_aic$Model[3] <- "m2"
rec_td_aic$LL[3] <- logLik(rec_g2t)
rec_td_aic$AICc[4] <- AICc(rec_g3t)
rec_td_aic$Model[4] <- "m3"
rec_td_aic$LL[4] <- logLik(rec_g3t)
rec_td_aic$AICc[5] <- AICc(rec_g4t)
rec_td_aic$Model[5] <- "m4"
rec_td_aic$LL[5] <- logLik(rec_g4t)
rec_td_aic$AICc[6] <- AICc(rec_g5t)
rec_td_aic$Model[6] <- "m5"
rec_td_aic$LL[6] <- logLik(rec_g5t)
rec_td_aic$AICc <- as.numeric(rec_td_aic$AICc)
rec_td_aic$LL <- as.numeric(rec_td_aic$LL)
rec_td_aic <- rec_td_aic[order(rec_td_aic$AICc), ]
rec_td_aic$Delta_AICc <- "0.00"
rec_td_aic$Delta_AICc[2] <- round(rec_td_aic$AICc[1]-rec_td_aic$AICc[2], 2)
rec_td_aic$Delta_AICc[3] <- round(rec_td_aic$AICc[1]-rec_td_aic$AICc[3], 2)
rec_td_aic$Delta_AICc[4] <- round(rec_td_aic$AICc[1]-rec_td_aic$AICc[4], 2)
rec_td_aic$Delta_AICc[5] <- round(rec_td_aic$AICc[1]-rec_td_aic$AICc[5], 2)
rec_td_aic$Delta_AICc[6] <- round(rec_td_aic$AICc[1]-rec_td_aic$AICc[6], 2)
rec_td_aic  #Null model is best

rec_ld_aic <- as.data.frame(1:6)
rec_ld_aic$AICc <- "NA"
rec_ld_aic$Model <- "NA"
rec_ld_aic$LL <- "NA"
rec_ld_aic$AICc[1] <- AICc(recl_gnull)
rec_ld_aic$Model[1] <- "Null"
rec_ld_aic$LL[1] <- logLik(recl_gnull)
rec_ld_aic$AICc[2] <- AICc(rec_g1l)
rec_ld_aic$Model[2] <- "m1"
rec_ld_aic$LL[2] <- logLik(rec_g1l)
rec_ld_aic$AICc[3] <- AICc(rec_g2l)
rec_ld_aic$Model[3] <- "m2"
rec_ld_aic$LL[3] <- logLik(rec_g2l)
rec_ld_aic$AICc[4] <- AICc(rec_g3l)
rec_ld_aic$Model[4] <- "m3"
rec_ld_aic$LL[4] <- logLik(rec_g3l)
rec_ld_aic$AICc[5] <- AICc(rec_g4l)
rec_ld_aic$Model[5] <- "m4"
rec_ld_aic$LL[5] <- logLik(rec_g4l)
rec_ld_aic$AICc[6] <- AICc(rec_g5l)
rec_ld_aic$Model[6] <- "m5"
rec_ld_aic$LL[6] <- logLik(rec_g5l)
rec_ld_aic$AICc <- as.numeric(rec_ld_aic$AICc)
rec_ld_aic$LL <- as.numeric(rec_ld_aic$LL)
rec_ld_aic <- rec_ld_aic[order(rec_ld_aic$AICc), ]
rec_ld_aic$Delta_AICc <- "0.00"
rec_ld_aic$Delta_AICc[2] <- round(rec_ld_aic$AICc[1]-rec_ld_aic$AICc[2], 2)
rec_ld_aic$Delta_AICc[3] <- round(rec_ld_aic$AICc[1]-rec_ld_aic$AICc[3], 2)
rec_ld_aic$Delta_AICc[4] <- round(rec_ld_aic$AICc[1]-rec_ld_aic$AICc[4], 2)
rec_ld_aic$Delta_AICc[5] <- round(rec_ld_aic$AICc[1]-rec_ld_aic$AICc[5], 2)
rec_ld_aic$Delta_AICc[6] <- round(rec_ld_aic$AICc[1]-rec_ld_aic$AICc[6], 2)
rec_ld_aic  #Null model is best



# GAM s(bs = 'tp')
rec_ts_aic <- as.data.frame(1:6)
rec_ts_aic$AICc <- "NA"
rec_ts_aic$Model <- "NA"
rec_ts_aic$LL <- "NA"
rec_ts_aic$AICc[1] <- AICc(rect_gnull)
rec_ts_aic$Model[1] <- "Null"
rec_ts_aic$LL[1] <- logLik(rect_gnull)
rec_ts_aic$AICc[2] <- AICc(rec_g1.1t)
rec_ts_aic$Model[2] <- "m1"
rec_ts_aic$LL[2] <- logLik(rec_g1.1t)
rec_ts_aic$AICc[3] <- AICc(rec_g2.1t)
rec_ts_aic$Model[3] <- "m2"
rec_ts_aic$LL[3] <- logLik(rec_g2.1t)
rec_ts_aic$AICc[4] <- AICc(rec_g3.1t)
rec_ts_aic$Model[4] <- "m3"
rec_ts_aic$LL[4] <- logLik(rec_g3.1t)
rec_ts_aic$AICc[5] <- AICc(rec_g4.1t)
rec_ts_aic$Model[5] <- "m4"
rec_ts_aic$LL[5] <- logLik(rec_g4.1t)
rec_ts_aic$AICc[6] <- AICc(rec_g5.1t)
rec_ts_aic$Model[6] <- "m5"
rec_ts_aic$LL[6] <- logLik(rec_g5.1t)
rec_ts_aic$AICc <- as.numeric(rec_ts_aic$AICc)
rec_ts_aic$LL <- as.numeric(rec_ts_aic$LL)
rec_ts_aic <- rec_ts_aic[order(rec_ts_aic$AICc), ]
rec_ts_aic$Delta_AICc <- "0.00"
rec_ts_aic$Delta_AICc[2] <- round(rec_ts_aic$AICc[1]-rec_ts_aic$AICc[2], 2)
rec_ts_aic$Delta_AICc[3] <- round(rec_ts_aic$AICc[1]-rec_ts_aic$AICc[3], 2)
rec_ts_aic$Delta_AICc[4] <- round(rec_ts_aic$AICc[1]-rec_ts_aic$AICc[4], 2)
rec_ts_aic$Delta_AICc[5] <- round(rec_ts_aic$AICc[1]-rec_ts_aic$AICc[5], 2)
rec_ts_aic$Delta_AICc[6] <- round(rec_ts_aic$AICc[1]-rec_ts_aic$AICc[6], 2)
rec_ts_aic  #Null model is best

rec_ls_aic <- as.data.frame(1:6)
rec_ls_aic$AICc <- "NA"
rec_ls_aic$Model <- "NA"
rec_ls_aic$LL <- "NA"
rec_ls_aic$AICc[1] <- AICc(recl_gnull)
rec_ls_aic$Model[1] <- "Null"
rec_ls_aic$LL[1] <- logLik(recl_gnull)
rec_ls_aic$AICc[2] <- AICc(rec_g1.1l)
rec_ls_aic$Model[2] <- "m1"
rec_ls_aic$LL[2] <- logLik(rec_g1.1l)
rec_ls_aic$AICc[3] <- AICc(rec_g2.1l)
rec_ls_aic$Model[3] <- "m2"
rec_ls_aic$LL[3] <- logLik(rec_g2.1l)
rec_ls_aic$AICc[4] <- AICc(rec_g3.1l)
rec_ls_aic$Model[4] <- "m3"
rec_ls_aic$LL[4] <- logLik(rec_g3.1l)
rec_ls_aic$AICc[5] <- AICc(rec_g4.1l)
rec_ls_aic$Model[5] <- "m4"
rec_ls_aic$LL[5] <- logLik(rec_g4.1l)
rec_ls_aic$AICc[6] <- AICc(rec_g5.1l)
rec_ls_aic$Model[6] <- "m5"
rec_ls_aic$LL[6] <- logLik(rec_g5.1l)
rec_ls_aic$AICc <- as.numeric(rec_ls_aic$AICc)
rec_ls_aic$LL <- as.numeric(rec_ls_aic$LL)
rec_ls_aic <- rec_ls_aic[order(rec_ls_aic$AICc), ]
rec_ls_aic$Delta_AICc <- "0.00"
rec_ls_aic$Delta_AICc[2] <- round(rec_ls_aic$AICc[1]-rec_ls_aic$AICc[2], 2)
rec_ls_aic$Delta_AICc[3] <- round(rec_ls_aic$AICc[1]-rec_ls_aic$AICc[3], 2)
rec_ls_aic$Delta_AICc[4] <- round(rec_ls_aic$AICc[1]-rec_ls_aic$AICc[4], 2)
rec_ls_aic$Delta_AICc[5] <- round(rec_ls_aic$AICc[1]-rec_ls_aic$AICc[5], 2)
rec_ls_aic$Delta_AICc[6] <- round(rec_ls_aic$AICc[1]-rec_ls_aic$AICc[6], 2)
rec_ls_aic  #Null model is best


# GAM ti(bs='cr')
rec_tsc_aic <- as.data.frame(1:6)
rec_tsc_aic$AICc <- "NA"
rec_tsc_aic$Model <- "NA"
rec_tsc_aic$LL <- "NA"
rec_tsc_aic$AICc[1] <- AICc(rect_gnull)
rec_tsc_aic$Model[1] <- "Null"
rec_tsc_aic$LL[1] <- logLik(rect_gnull)
rec_tsc_aic$AICc[2] <- AICc(rec_g1t)
rec_tsc_aic$Model[2] <- "m1"
rec_tsc_aic$LL[2] <- logLik(rec_g1t)
rec_tsc_aic$AICc[3] <- AICc(rec_g2.2t)
rec_tsc_aic$Model[3] <- "m2"
rec_tsc_aic$LL[3] <- logLik(rec_g2.2t)
rec_tsc_aic$AICc[4] <- AICc(rec_g3.2t)
rec_tsc_aic$Model[4] <- "m3"
rec_tsc_aic$LL[4] <- logLik(rec_g3.2t)
rec_tsc_aic$AICc[5] <- AICc(rec_g4.2t)
rec_tsc_aic$Model[5] <- "m4"
rec_tsc_aic$LL[5] <- logLik(rec_g4.2t)
rec_tsc_aic$AICc[6] <- AICc(rec_g5.2t)
rec_tsc_aic$Model[6] <- "m5"
rec_tsc_aic$LL[6] <- logLik(rec_g5.2t)
rec_tsc_aic$AICc <- as.numeric(rec_tsc_aic$AICc)
rec_tsc_aic$LL <- as.numeric(rec_tsc_aic$LL)
rec_tsc_aic <- rec_tsc_aic[order(rec_tsc_aic$AICc), ]
rec_tsc_aic$Delta_AICc <- "0.00"
rec_tsc_aic$Delta_AICc[2] <- round(rec_tsc_aic$AICc[1]-rec_tsc_aic$AICc[2], 2)
rec_tsc_aic$Delta_AICc[3] <- round(rec_tsc_aic$AICc[1]-rec_tsc_aic$AICc[3], 2)
rec_tsc_aic$Delta_AICc[4] <- round(rec_tsc_aic$AICc[1]-rec_tsc_aic$AICc[4], 2)
rec_tsc_aic$Delta_AICc[5] <- round(rec_tsc_aic$AICc[1]-rec_tsc_aic$AICc[5], 2)
rec_tsc_aic$Delta_AICc[6] <- round(rec_tsc_aic$AICc[1]-rec_tsc_aic$AICc[6], 2)
rec_tsc_aic  #Null model is best

rec_lsc_aic <- as.data.frame(1:6)
rec_lsc_aic$AICc <- "NA"
rec_lsc_aic$Model <- "NA"
rec_lsc_aic$LL <- "NA"
rec_lsc_aic$AICc[1] <- AICc(recl_gnull)
rec_lsc_aic$Model[1] <- "Null"
rec_lsc_aic$LL[1] <- logLik(recl_gnull)
rec_lsc_aic$AICc[2] <- AICc(rec_g1l)
rec_lsc_aic$Model[2] <- "m1"
rec_lsc_aic$LL[2] <- logLik(rec_g1l)
rec_lsc_aic$AICc[3] <- AICc(rec_g2.2l)
rec_lsc_aic$Model[3] <- "m2"
rec_lsc_aic$LL[3] <- logLik(rec_g2.2l)
rec_lsc_aic$AICc[4] <- AICc(rec_g3.2l)
rec_lsc_aic$Model[4] <- "m3"
rec_lsc_aic$LL[4] <- logLik(rec_g3.2l)
rec_lsc_aic$AICc[5] <- AICc(rec_g4.2l)
rec_lsc_aic$Model[5] <- "m4"
rec_lsc_aic$LL[5] <- logLik(rec_g4.2l)
rec_lsc_aic$AICc[6] <- AICc(rec_g5.2l)
rec_lsc_aic$Model[6] <- "m5"
rec_lsc_aic$LL[6] <- logLik(rec_g5.2l)
rec_lsc_aic$AICc <- as.numeric(rec_lsc_aic$AICc)
rec_lsc_aic$LL <- as.numeric(rec_lsc_aic$LL)
rec_lsc_aic <- rec_lsc_aic[order(rec_lsc_aic$AICc), ]
rec_lsc_aic$Delta_AICc <- "0.00"
rec_lsc_aic$Delta_AICc[2] <- round(rec_lsc_aic$AICc[1]-rec_lsc_aic$AICc[2], 2)
rec_lsc_aic$Delta_AICc[3] <- round(rec_lsc_aic$AICc[1]-rec_lsc_aic$AICc[3], 2)
rec_lsc_aic$Delta_AICc[4] <- round(rec_lsc_aic$AICc[1]-rec_lsc_aic$AICc[4], 2)
rec_lsc_aic$Delta_AICc[5] <- round(rec_lsc_aic$AICc[1]-rec_lsc_aic$AICc[5], 2)
rec_lsc_aic$Delta_AICc[6] <- round(rec_lsc_aic$AICc[1]-rec_lsc_aic$AICc[6], 2)
rec_lsc_aic  #Null model is best



# 5.1.2 Proportion saplings ----
rsap_tnull <- glmer(Proportion_saplings ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
rsap_gt_null <- gam(Proportion_saplings ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')

rsap_lnull <- glmer(Proportion_saplings ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)
rsap_gl_null <- gam(Proportion_saplings ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')



rsap_m1t <- glmer(Proportion_saplings ~ r_TSF + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m1t)
rsap_g1t <- gam(Proportion_saplings ~ s(TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g1t)
plot(rsap_g1t)
rsap_g1.1t <- gam(Proportion_saplings ~ s(TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g1.1t)
plot(rsap_g1.1t)

rsap_m1l <- glmer(Proportion_saplings ~ r_TSF + (1 | Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(rsap_m1l) # Not a very good model
rsap_g1l <- gam(Proportion_saplings ~ s(TSF, k = 4), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g1l)
plot(rsap_g1l)
rsap_g1.1l <- gam(Proportion_saplings ~ s(TSF, bs = 'cr', k = 4), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g1.1l)
plot(rsap_g1.1l)


rsap_m2t <- glmer(Proportion_saplings ~ r_TSF * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m2t)
rsap_g2t <- gam(Proportion_saplings ~ s(TSF) + s(Latitude) + ti(Latitude, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g2t)
plot(rsap_g2t)
rsap_g2.1t <- gam(Proportion_saplings ~ s(TSF) + s(Latitude) + ti(Latitude, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g2.1t)
plot(rsap_g2.1t)
rsap_g2.2t <- gam(Proportion_saplings ~ s(TSF, bs = 'cr', k = 5) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g2.1t)
plot(rsap_g2.2t)


#rsap_m2l <- glmer(Proportion_saplings ~ TSF * Latitude + (1|Location/Transect), family = binomial, data = lit_transects) # This model fails due to the low number of points available to model an interaction. We can only model an additive effect of TSF and latitude
rsap_m2.1l <- glmer(Proportion_saplings ~ r_TSF + r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(rsap_m2.1l) # Not very good model
rsap_g2l <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(Latitude, k = 4) + ti(Latitude, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g2l)
plot(rsap_g2l)
rsap_g2.1l <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(Latitude, k = 4) + ti(Latitude, by = TSF, bs = 'tp', k = 8), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g2.1l)
plot(rsap_g2.1l)
rsap_g2.2l <- gam(Proportion_saplings ~ s(TSF, bs = 'cr', k = 4) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g2.1l)
plot(rsap_g2.2l)



rsap_m3t <- glmer(Proportion_saplings ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
rsap_g3t <- gam(Proportion_saplings ~ s(TSF) + s(FPC) + ti(FPC, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g3t)
plot(rsap_g3t)
rsap_g3.1t <- gam(Proportion_saplings ~ s(TSF) + s(FPC) + ti(FPC, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g3.1t)
plot(rsap_g3.1t)
rsap_g3.2t <- gam(Proportion_saplings ~ s(TSF, bs = 'cr', k = 5) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g3.1t)
plot(rsap_g3.2t)


rsap_m3l <- glmer(Proportion_saplings ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=200)))
summary(rsap_m3l) # Not a very good model
rsap_g3l <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(FPC, k = 4) + ti(FPC, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g3l)
plot(rsap_g3l)
rsap_g3.1l <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(FPC, k = 4) + ti(FPC, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g3.1t)
plot(rsap_g2.1t)
rsap_g3.2l <- gam(Proportion_saplings ~ s(TSF, bs = 'cr', k = 4) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g3.2t)
plot(rsap_g3.2l)



rsap_m4t <-  glmer(Proportion_saplings ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m4t) # Not a very good model
rsap_g4t <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(Precip, k = 4) + ti(Precip, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g4t)
plot(rsap_g4t)
rsap_g4.1t <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(Precip, k = 4) + ti(Precip, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g4.1t)
plot(rsap_g4.1t)
rsap_g4.2t <- gam(Proportion_saplings ~ s(TSF, bs = 'cr', k = 4) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g4.2t)
plot(rsap_g4.2t)



rsap_m4l <-  glmer(Proportion_saplings ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=130)))
summary(rsap_m4l) # Not a very good model
rsap_g4l <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(Precip, k = 4) + ti(Precip, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g4l)
plot(rsap_g4l)
rsap_g4.1l <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(Precip, k = 4) + ti(Precip, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g4.1l)
plot(rsap_g4.1l)
rsap_g4.2l <- gam(Proportion_saplings ~ s(TSF, bs = 'cr', k = 4) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g4.2l)
plot(rsap_g4.2l)


rsap_m5t <-  glmer(Proportion_saplings ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m5t) # Not a very good model
rsap_g5t <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(Temp, k = 4) + ti(Temp, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g5t)
plot(rsap_g5t)
rsap_g5.1t <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(Temp, k = 4) + ti(Temp, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g5.1t)
plot(rsap_g5.1t)
rsap_g5.2t <- gam(Proportion_saplings ~ s(TSF, bs = 'cr', k = 4) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g5.2t)
plot(rsap_g5.2t)
rsap_g5.3t <- gam(Proportion_saplings ~ s(TSF, bs = 'cr', k = 3) + s(Temp, bs = 'cc', k = 3) + ti(Temp, by = TSF, bs = 'cc'), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g5.3t)
plot(rsap_g5.3t)

rsap_m5l <-  glmer(Proportion_saplings ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rsap_m5l) # Not a very good model
rsap_g5l <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(Temp, k = 4) + ti(Temp, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g5l)
plot(rsap_g5l)
rsap_g5.1l <- gam(Proportion_saplings ~ s(TSF, k = 4) + s(Temp, k = 4) + ti(Temp, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g5.1l)
plot(rsap_g5.1l)
rsap_g5.2l <- gam(Proportion_saplings ~ s(TSF, bs = 'cr', k = 4) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(rsap_g5.2l)
plot(rsap_g5.2l)

# Compare AICs for each modelling type
# GLMER
rsap_aic_t <- list(rsap_tnull, rsap_m1t, rsap_m2t, rsap_m3t, rsap_m4t, rsap_m5t)
aictab(rsap_aic_t) # Model 5 and 6 are the best models - those including precipitation and temperature seasonality

rsap_aic_l <- list(rsap_lnull, rsap_m1l, rsap_m2.1l, rsap_m3l, rsap_m4l, rsap_m5l)
aictab(rsap_aic_l) # Null model is best

# GAM defaults
rsap_td_aic <- as.data.frame(1:6)
rsap_td_aic$AICc <- "NA"
rsap_td_aic$Model <- "NA"
rsap_td_aic$LL <- "NA"
rsap_td_aic$AICc[1] <- AICc(rsap_gt_null)
rsap_td_aic$Model[1] <- "Null"
rsap_td_aic$LL[1] <- logLik(rsap_gt_null)
rsap_td_aic$AICc[2] <- AICc(rsap_g1t)
rsap_td_aic$Model[2] <- "m1"
rsap_td_aic$LL[2] <- logLik(rsap_g1t)
rsap_td_aic$AICc[3] <- AICc(rsap_g2t)
rsap_td_aic$Model[3] <- "m2"
rsap_td_aic$LL[3] <- logLik(rsap_g2t)
rsap_td_aic$AICc[4] <- AICc(rsap_g3t)
rsap_td_aic$Model[4] <- "m3"
rsap_td_aic$LL[4] <- logLik(rsap_g3t)
rsap_td_aic$AICc[5] <- AICc(rsap_g4t)
rsap_td_aic$Model[5] <- "m4"
rsap_td_aic$LL[5] <- logLik(rsap_g4t)
rsap_td_aic$AICc[6] <- AICc(rsap_g5t)
rsap_td_aic$Model[6] <- "m5"
rsap_td_aic$LL[6] <- logLik(rsap_g5t)
rsap_td_aic$AICc <- as.numeric(rsap_td_aic$AICc)
rsap_td_aic$LL <- as.numeric(rsap_td_aic$LL)
rsap_td_aic <- rsap_td_aic[order(rsap_td_aic$AICc), ]
rsap_td_aic$Delta_AICc <- "0.00"
rsap_td_aic$Delta_AICc[2] <- round(rsap_td_aic$AICc[1]-rsap_td_aic$AICc[2], 2)
rsap_td_aic$Delta_AICc[3] <- round(rsap_td_aic$AICc[1]-rsap_td_aic$AICc[3], 2)
rsap_td_aic$Delta_AICc[4] <- round(rsap_td_aic$AICc[1]-rsap_td_aic$AICc[4], 2)
rsap_td_aic$Delta_AICc[5] <- round(rsap_td_aic$AICc[1]-rsap_td_aic$AICc[5], 2)
rsap_td_aic$Delta_AICc[6] <- round(rsap_td_aic$AICc[1]-rsap_td_aic$AICc[6], 2)
rsap_td_aic  # Model 5 is best - temperature only. No other models ranked within delta AICc <2.  

rsap_ld_aic <- as.data.frame(1:6)
rsap_ld_aic$AICc <- "NA"
rsap_ld_aic$Model <- "NA"
rsap_ld_aic$LL <- "NA"
rsap_ld_aic$AICc[1] <- AICc(rsap_gl_null)
rsap_ld_aic$Model[1] <- "Null"
rsap_ld_aic$LL[1] <- logLik(rsap_gl_null)
rsap_ld_aic$AICc[2] <- AICc(rsap_g1l)
rsap_ld_aic$Model[2] <- "m1"
rsap_ld_aic$LL[2] <- logLik(rsap_g1l)
rsap_ld_aic$AICc[3] <- AICc(rsap_g2l)
rsap_ld_aic$Model[3] <- "m2"
rsap_ld_aic$LL[3] <- logLik(rsap_g2l)
rsap_ld_aic$AICc[4] <- AICc(rsap_g3l)
rsap_ld_aic$Model[4] <- "m3"
rsap_ld_aic$LL[4] <- logLik(rsap_g3l)
rsap_ld_aic$AICc[5] <- AICc(rsap_g4l)
rsap_ld_aic$Model[5] <- "m4"
rsap_ld_aic$LL[5] <- logLik(rsap_g4l)
rsap_ld_aic$AICc[6] <- AICc(rsap_g5l)
rsap_ld_aic$Model[6] <- "m5"
rsap_ld_aic$LL[6] <- logLik(rsap_g5l)
rsap_ld_aic$AICc <- as.numeric(rsap_ld_aic$AICc)
rsap_ld_aic$LL <- as.numeric(rsap_ld_aic$LL)
rsap_ld_aic <- rsap_ld_aic[order(rsap_ld_aic$AICc), ]
rsap_ld_aic$Delta_AICc <- "0.00"
rsap_ld_aic$Delta_AICc[2] <- round(rsap_ld_aic$AICc[1]-rsap_ld_aic$AICc[2], 2)
rsap_ld_aic$Delta_AICc[3] <- round(rsap_ld_aic$AICc[1]-rsap_ld_aic$AICc[3], 2)
rsap_ld_aic$Delta_AICc[4] <- round(rsap_ld_aic$AICc[1]-rsap_ld_aic$AICc[4], 2)
rsap_ld_aic$Delta_AICc[5] <- round(rsap_ld_aic$AICc[1]-rsap_ld_aic$AICc[5], 2)
rsap_ld_aic$Delta_AICc[6] <- round(rsap_ld_aic$AICc[1]-rsap_ld_aic$AICc[6], 2)
rsap_ld_aic  #Null model is best



# GAM s(bs = 'tp')
rsap_ts_aic <- as.data.frame(1:6)
rsap_ts_aic$AICc <- "NA"
rsap_ts_aic$Model <- "NA"
rsap_ts_aic$LL <- "NA"
rsap_ts_aic$AICc[1] <- AICc(rsap_gt_null)
rsap_ts_aic$Model[1] <- "Null"
rsap_ts_aic$LL[1] <- logLik(rsap_gt_null)
rsap_ts_aic$AICc[2] <- AICc(rsap_g1.1t)
rsap_ts_aic$Model[2] <- "m1"
rsap_ts_aic$LL[2] <- logLik(rsap_g1.1t)
rsap_ts_aic$AICc[3] <- AICc(rsap_g2.1t)
rsap_ts_aic$Model[3] <- "m2"
rsap_ts_aic$LL[3] <- logLik(rsap_g2.1t)
rsap_ts_aic$AICc[4] <- AICc(rsap_g3.1t)
rsap_ts_aic$Model[4] <- "m3"
rsap_ts_aic$LL[4] <- logLik(rsap_g3.1t)
rsap_ts_aic$AICc[5] <- AICc(rsap_g4.1t)
rsap_ts_aic$Model[5] <- "m4"
rsap_ts_aic$LL[5] <- logLik(rsap_g4.1t)
rsap_ts_aic$AICc[6] <- AICc(rsap_g5.1t)
rsap_ts_aic$Model[6] <- "m5"
rsap_ts_aic$LL[6] <- logLik(rsap_g5.1t)
rsap_ts_aic$AICc <- as.numeric(rsap_ts_aic$AICc)
rsap_ts_aic$LL <- as.numeric(rsap_ts_aic$LL)
rsap_ts_aic <- rsap_ts_aic[order(rsap_ts_aic$AICc), ]
rsap_ts_aic$Delta_AICc <- "0.00"
rsap_ts_aic$Delta_AICc[2] <- round(rsap_ts_aic$AICc[1]-rsap_ts_aic$AICc[2], 2)
rsap_ts_aic$Delta_AICc[3] <- round(rsap_ts_aic$AICc[1]-rsap_ts_aic$AICc[3], 2)
rsap_ts_aic$Delta_AICc[4] <- round(rsap_ts_aic$AICc[1]-rsap_ts_aic$AICc[4], 2)
rsap_ts_aic$Delta_AICc[5] <- round(rsap_ts_aic$AICc[1]-rsap_ts_aic$AICc[5], 2)
rsap_ts_aic$Delta_AICc[6] <- round(rsap_ts_aic$AICc[1]-rsap_ts_aic$AICc[6], 2)
rsap_ts_aic  #Model 5 is best.

rsap_ls_aic <- as.data.frame(1:6)
rsap_ls_aic$AICc <- "NA"
rsap_ls_aic$Model <- "NA"
rsap_ls_aic$LL <- "NA"
rsap_ls_aic$AICc[1] <- AICc(rsap_gl_null)
rsap_ls_aic$Model[1] <- "Null"
rsap_ls_aic$LL[1] <- logLik(rsap_gl_null)
rsap_ls_aic$AICc[2] <- AICc(rsap_g1.1l)
rsap_ls_aic$Model[2] <- "m1"
rsap_ls_aic$LL[2] <- logLik(rsap_g1.1l)
rsap_ls_aic$AICc[3] <- AICc(rsap_g2.1l)
rsap_ls_aic$Model[3] <- "m2"
rsap_ls_aic$LL[3] <- logLik(rsap_g2.1l)
rsap_ls_aic$AICc[4] <- AICc(rsap_g3.1l)
rsap_ls_aic$Model[4] <- "m3"
rsap_ls_aic$LL[4] <- logLik(rsap_g3.1l)
rsap_ls_aic$AICc[5] <- AICc(rsap_g4.1l)
rsap_ls_aic$Model[5] <- "m4"
rsap_ls_aic$LL[5] <- logLik(rsap_g4.1l)
rsap_ls_aic$AICc[6] <- AICc(rsap_g5.1l)
rsap_ls_aic$Model[6] <- "m5"
rsap_ls_aic$LL[6] <- logLik(rsap_g5.1l)
rsap_ls_aic$AICc <- as.numeric(rsap_ls_aic$AICc)
rsap_ls_aic$LL <- as.numeric(rsap_ls_aic$LL)
rsap_ls_aic <- rsap_ls_aic[order(rsap_ls_aic$AICc), ]
rsap_ls_aic$Delta_AICc <- "0.00"
rsap_ls_aic$Delta_AICc[2] <- round(rsap_ls_aic$AICc[1]-rsap_ls_aic$AICc[2], 2)
rsap_ls_aic$Delta_AICc[3] <- round(rsap_ls_aic$AICc[1]-rsap_ls_aic$AICc[3], 2)
rsap_ls_aic$Delta_AICc[4] <- round(rsap_ls_aic$AICc[1]-rsap_ls_aic$AICc[4], 2)
rsap_ls_aic$Delta_AICc[5] <- round(rsap_ls_aic$AICc[1]-rsap_ls_aic$AICc[5], 2)
rsap_ls_aic$Delta_AICc[6] <- round(rsap_ls_aic$AICc[1]-rsap_ls_aic$AICc[6], 2)
rsap_ls_aic  #Null model is best


# GAM ti(bs='cr')
rsap_tsc_aic <- as.data.frame(1:6)
rsap_tsc_aic$AICc <- "NA"
rsap_tsc_aic$Model <- "NA"
rsap_tsc_aic$LL <- "NA"
rsap_tsc_aic$AICc[1] <- AICc(rsap_gt_null)
rsap_tsc_aic$Model[1] <- "Null"
rsap_tsc_aic$LL[1] <- logLik(rsap_gt_null)
rsap_tsc_aic$AICc[2] <- AICc(rsap_g1t)
rsap_tsc_aic$Model[2] <- "m1"
rsap_tsc_aic$LL[2] <- logLik(rsap_g1t)
rsap_tsc_aic$AICc[3] <- AICc(rsap_g2.2t)
rsap_tsc_aic$Model[3] <- "m2"
rsap_tsc_aic$LL[3] <- logLik(rsap_g2.2t)
rsap_tsc_aic$AICc[4] <- AICc(rsap_g3.2t)
rsap_tsc_aic$Model[4] <- "m3"
rsap_tsc_aic$LL[4] <- logLik(rsap_g3.2t)
rsap_tsc_aic$AICc[5] <- AICc(rsap_g4.2t)
rsap_tsc_aic$Model[5] <- "m4"
rsap_tsc_aic$LL[5] <- logLik(rsap_g4.2t)
rsap_tsc_aic$AICc[6] <- AICc(rsap_g5.2t)
rsap_tsc_aic$Model[6] <- "m5"
rsap_tsc_aic$LL[6] <- logLik(rsap_g5.2t)
rsap_tsc_aic$AICc <- as.numeric(rsap_tsc_aic$AICc)
rsap_tsc_aic$LL <- as.numeric(rsap_tsc_aic$LL)
rsap_tsc_aic <- rsap_tsc_aic[order(rsap_tsc_aic$AICc), ]
rsap_tsc_aic$Delta_AICc <- "0.00"
rsap_tsc_aic$Delta_AICc[2] <- round(rsap_tsc_aic$AICc[1]-rsap_tsc_aic$AICc[2], 2)
rsap_tsc_aic$Delta_AICc[3] <- round(rsap_tsc_aic$AICc[1]-rsap_tsc_aic$AICc[3], 2)
rsap_tsc_aic$Delta_AICc[4] <- round(rsap_tsc_aic$AICc[1]-rsap_tsc_aic$AICc[4], 2)
rsap_tsc_aic$Delta_AICc[5] <- round(rsap_tsc_aic$AICc[1]-rsap_tsc_aic$AICc[5], 2)
rsap_tsc_aic$Delta_AICc[6] <- round(rsap_tsc_aic$AICc[1]-rsap_tsc_aic$AICc[6], 2)
rsap_tsc_aic  #Model 5 is best.

rsap_lsc_aic <- as.data.frame(1:6)
rsap_lsc_aic$AICc <- "NA"
rsap_lsc_aic$Model <- "NA"
rsap_lsc_aic$LL <- "NA"
rsap_lsc_aic$AICc[1] <- AICc(rsap_gl_null)
rsap_lsc_aic$Model[1] <- "Null"
rsap_lsc_aic$LL[1] <- logLik(rsap_gl_null)
rsap_lsc_aic$AICc[2] <- AICc(rsap_g1l)
rsap_lsc_aic$Model[2] <- "m1"
rsap_lsc_aic$LL[2] <- logLik(rsap_g1l)
rsap_lsc_aic$AICc[3] <- AICc(rsap_g2.2l)
rsap_lsc_aic$Model[3] <- "m2"
rsap_lsc_aic$LL[3] <- logLik(rsap_g2.2l)
rsap_lsc_aic$AICc[4] <- AICc(rsap_g3.2l)
rsap_lsc_aic$Model[4] <- "m3"
rsap_lsc_aic$LL[4] <- logLik(rsap_g3.2l)
rsap_lsc_aic$AICc[5] <- AICc(rsap_g4.2l)
rsap_lsc_aic$Model[5] <- "m4"
rsap_lsc_aic$LL[5] <- logLik(rsap_g4.2l)
rsap_lsc_aic$AICc[6] <- AICc(rsap_g5.2l)
rsap_lsc_aic$Model[6] <- "m5"
rsap_lsc_aic$LL[6] <- logLik(rsap_g5.2l)
rsap_lsc_aic$AICc <- as.numeric(rsap_lsc_aic$AICc)
rsap_lsc_aic$LL <- as.numeric(rsap_lsc_aic$LL)
rsap_lsc_aic <- rsap_lsc_aic[order(rsap_lsc_aic$AICc), ]
rsap_lsc_aic$Delta_AICc <- "0.00"
rsap_lsc_aic$Delta_AICc[2] <- round(rsap_lsc_aic$AICc[1]-rsap_lsc_aic$AICc[2], 2)
rsap_lsc_aic$Delta_AICc[3] <- round(rsap_lsc_aic$AICc[1]-rsap_lsc_aic$AICc[3], 2)
rsap_lsc_aic$Delta_AICc[4] <- round(rsap_lsc_aic$AICc[1]-rsap_lsc_aic$AICc[4], 2)
rsap_lsc_aic$Delta_AICc[5] <- round(rsap_lsc_aic$AICc[1]-rsap_lsc_aic$AICc[5], 2)
rsap_lsc_aic$Delta_AICc[6] <- round(rsap_lsc_aic$AICc[1]-rsap_lsc_aic$AICc[6], 2)
rsap_lsc_aic  #Null model is best



# Predict for proportion saplings data ----
nsap_m5l <- expand.grid(TSF = min(tor_transects$TSF),
                         Temp = seq(min(tor_transects$Temp), max(tor_transects$Temp), length = 50))
nsap_m5a <- expand.grid(TSF = mean(tor_transects$TSF),
                        Temp = seq(min(tor_transects$Temp), max(tor_transects$Temp), length = 50))
nsap_m5h <- expand.grid(TSF = max(tor_transects$TSF),
                        Temp = seq(min(tor_transects$Temp), max(tor_transects$Temp), length = 50))



# GLMER 
nsap_m5l_lm <-  expand.grid(r_TSF = min(tor_transects$r_TSF),
                            r_Temp = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 50))
nsap_m5a_lm <- expand.grid(r_TSF = mean(tor_transects$r_TSF),
                           r_Temp = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 50))
nsap_m5h_lm <- expand.grid(r_TSF = max(tor_transects$r_TSF),
                           r_Temp = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 50))


nsap_m4l_lm <-  expand.grid(r_TSF = min(tor_transects$r_TSF),
                            r_Precip = seq(min(tor_transects$r_Precip), max(tor_transects$r_Precip), length = 50))
nsap_m4a_lm <- expand.grid(r_TSF = mean(tor_transects$r_TSF),
                           r_Precip = seq(min(tor_transects$r_Precip), max(tor_transects$r_Precip), length = 50))
nsap_m4h_lm <- expand.grid(r_TSF = max(tor_transects$r_TSF),
                           r_Precip = seq(min(tor_transects$r_Precip), max(tor_transects$r_Precip), length = 50))




psap_m5l_lm <- predictSE(rsap_m5t, newdata = nsap_m5l_lm, se.fit = T, type = 'link')
nsap_m5l_lm$fit.link <- psap_m5l_lm$fit
nsap_m5l_lm$se.link <- psap_m5l_lm$se.fit
nsap_m5l_lm$lci.link <- nsap_m5l_lm$fit.link - (nsap_m5l_lm$se.link * 1.96)
nsap_m5l_lm$uci.link <- nsap_m5l_lm$fit.link + (nsap_m5l_lm$se.link * 1.96)
nsap_m5l_lm$fit <- invlogit(nsap_m5l_lm$fit.link)
nsap_m5l_lm$se <- invlogit(nsap_m5l_lm$se.link)
nsap_m5l_lm$lci <- invlogit(nsap_m5l_lm$lci.link)
nsap_m5l_lm$uci <- invlogit(nsap_m5l_lm$uci.link)

psap_m5a_lm <- predictSE(rsap_m5t, newdata = nsap_m5a_lm, se.fit = T, type = 'link')
nsap_m5a_lm$fit.link <- psap_m5a_lm$fit
nsap_m5a_lm$se.link <- psap_m5a_lm$se.fit
nsap_m5a_lm$lci.link <- nsap_m5a_lm$fit.link - (nsap_m5a_lm$se.link * 1.96)
nsap_m5a_lm$uci.link <- nsap_m5a_lm$fit.link + (nsap_m5a_lm$se.link * 1.96)
nsap_m5a_lm$fit <- invlogit(nsap_m5a_lm$fit.link)
nsap_m5a_lm$se <- invlogit(nsap_m5a_lm$se.link)
nsap_m5a_lm$lci <- invlogit(nsap_m5a_lm$lci.link)
nsap_m5a_lm$uci <- invlogit(nsap_m5a_lm$uci.link)

psap_m5h_lm <- predictSE(rsap_m5t, newdata = nsap_m5h_lm, se.fit = T, type = 'link')
nsap_m5h_lm$fit.link <- psap_m5h_lm$fit
nsap_m5h_lm$se.link <- psap_m5h_lm$se.fit
nsap_m5h_lm$lci.link <- nsap_m5h_lm$fit.link - (nsap_m5h_lm$se.link * 1.96)
nsap_m5h_lm$uci.link <- nsap_m5h_lm$fit.link + (nsap_m5h_lm$se.link * 1.96)
nsap_m5h_lm$fit <- invlogit(nsap_m5h_lm$fit.link)
nsap_m5h_lm$se <- invlogit(nsap_m5h_lm$se.link)
nsap_m5h_lm$lci <- invlogit(nsap_m5h_lm$lci.link)
nsap_m5h_lm$uci <- invlogit(nsap_m5h_lm$uci.link)




psap_m4l_lm <- predictSE(rsap_m4t, newdata = nsap_m4l_lm, se.fit = T, type = 'link')
nsap_m4l_lm$fit.link <- psap_m4l_lm$fit
nsap_m4l_lm$se.link <- psap_m4l_lm$se.fit
nsap_m4l_lm$lci.link <- nsap_m4l_lm$fit.link - (nsap_m4l_lm$se.link * 1.96)
nsap_m4l_lm$uci.link <- nsap_m4l_lm$fit.link + (nsap_m4l_lm$se.link * 1.96)
nsap_m4l_lm$fit <- invlogit(nsap_m4l_lm$fit.link)
nsap_m4l_lm$se <- invlogit(nsap_m4l_lm$se.link)
nsap_m4l_lm$lci <- invlogit(nsap_m4l_lm$lci.link)
nsap_m4l_lm$uci <- invlogit(nsap_m4l_lm$uci.link)

psap_m4a_lm <- predictSE(rsap_m4t, newdata = nsap_m4a_lm, se.fit = T, type = 'link')
nsap_m4a_lm$fit.link <- psap_m4a_lm$fit
nsap_m4a_lm$se.link <- psap_m4a_lm$se.fit
nsap_m4a_lm$lci.link <- nsap_m4a_lm$fit.link - (nsap_m4a_lm$se.link * 1.96)
nsap_m4a_lm$uci.link <- nsap_m4a_lm$fit.link + (nsap_m4a_lm$se.link * 1.96)
nsap_m4a_lm$fit <- invlogit(nsap_m4a_lm$fit.link)
nsap_m4a_lm$se <- invlogit(nsap_m4a_lm$se.link)
nsap_m4a_lm$lci <- invlogit(nsap_m4a_lm$lci.link)
nsap_m4a_lm$uci <- invlogit(nsap_m4a_lm$uci.link)

psap_m4h_lm <- predictSE(rsap_m4t, newdata = nsap_m4h_lm, se.fit = T, type = 'link')
nsap_m4h_lm$fit.link <- psap_m4h_lm$fit
nsap_m4h_lm$se.link <- psap_m4h_lm$se.fit
nsap_m4h_lm$lci.link <- nsap_m4h_lm$fit.link - (nsap_m4h_lm$se.link * 1.96)
nsap_m4h_lm$uci.link <- nsap_m4h_lm$fit.link + (nsap_m4h_lm$se.link * 1.96)
nsap_m4h_lm$fit <- invlogit(nsap_m4h_lm$fit.link)
nsap_m4h_lm$se <- invlogit(nsap_m4h_lm$se.link)
nsap_m4h_lm$lci <- invlogit(nsap_m4h_lm$lci.link)
nsap_m4h_lm$uci <- invlogit(nsap_m4h_lm$uci.link)



# GAM default
nsap_m5l_d <- nsap_m5l
nsap_m5a_d <- nsap_m5a
nsap_m5h_d <- nsap_m5h

psap_m5l_d <- predict(rsap_g5t, newdata = nsap_m5l_d, se.fit = T, type = 'response')
nsap_m5l_d$fit <- psap_m5l_d$fit
nsap_m5l_d$se <- psap_m5l_d$se.fit
nsap_m5l_d$lci <- nsap_m5l_d$fit - (nsap_m5l_d$se * 1.96)
nsap_m5l_d$uci <- nsap_m5l_d$fit + (nsap_m5l_d$se * 1.96)


psap_m5a_d <- predict(rsap_g5t, newdata = nsap_m5a_d, se.fit = T, type = 'response')
nsap_m5a_d$fit <- psap_m5a_d$fit
nsap_m5a_d$se <- psap_m5a_d$se.fit
nsap_m5a_d$lci <- nsap_m5a_d$fit - (nsap_m5a_d$se * 1.96)
nsap_m5a_d$uci <- nsap_m5a_d$fit + (nsap_m5a_d$se * 1.96)

psap_m5h_d <- predict(rsap_g5t, newdata = nsap_m5h_d, se.fit = T, type = 'response')
nsap_m5h_d$fit <- psap_m5h_d$fit
nsap_m5h_d$se <- psap_m5h_d$se.fit
nsap_m5h_d$lci <- nsap_m5h_d$fit - (nsap_m5h_d$se * 1.96)
nsap_m5h_d$uci <- nsap_m5h_d$fit + (nsap_m5h_d$se * 1.96)


# GAM bs = 'tp'
nsap_m5l_t <- nsap_m5l
nsap_m5a_t <- nsap_m5a
nsap_m5h_t <- nsap_m5h

psap_m5l_t <- predict(rsap_g5.1t, newdata = nsap_m5l_t, se.fit = T, type = 'response')
nsap_m5l_t$fit <- psap_m5l_t$fit
nsap_m5l_t$se <- psap_m5l_t$se.fit
nsap_m5l_t$lci <- nsap_m5l_t$fit - (nsap_m5l_t$se * 1.96)
nsap_m5l_t$uci <- nsap_m5l_t$fit + (nsap_m5l_t$se * 1.96)

psap_m5a_t <- predict(rsap_g5.1t, newdata = nsap_m5a_t, se.fit = T, type = 'response')
nsap_m5a_t$fit <- psap_m5a_t$fit
nsap_m5a_t$se <- psap_m5a_t$se.fit
nsap_m5a_t$lci <- nsap_m5a_t$fit - (nsap_m5a_t$se * 1.96)
nsap_m5a_t$uci <- nsap_m5a_t$fit + (nsap_m5a_t$se * 1.96)

psap_m5h_t <- predict(rsap_g5.1t, newdata = nsap_m5h_t, se.fit = T, type = 'response')
nsap_m5h_t$fit <- psap_m5h_t$fit
nsap_m5h_t$se <- psap_m5h_t$se.fit
nsap_m5h_t$lci <- nsap_m5h_t$fit - (nsap_m5h_t$se * 1.96)
nsap_m5h_t$uci <- nsap_m5h_t$fit + (nsap_m5h_t$se * 1.96)


# GAM bs = 'cr'
nsap_m5l_c <- nsap_m5l
nsap_m5a_c <- nsap_m5a
nsap_m5h_c <- nsap_m5h

psap_m5l_c <- predict(rsap_g5.2t, newdata = nsap_m5l_c, se.fit = T, type = 'response')
nsap_m5l_c$fit <- psap_m5l_c$fit
nsap_m5l_c$se <- psap_m5l_c$se.fit
nsap_m5l_c$lci <- nsap_m5l_c$fit - (nsap_m5l_c$se * 1.96)
nsap_m5l_c$uci <- nsap_m5l_c$fit + (nsap_m5l_c$se * 1.96)

psap_m5a_c <- predict(rsap_g5.2t, newdata = nsap_m5a_c, se.fit = T, type = 'response')
nsap_m5a_c$fit <- psap_m5a_c$fit
nsap_m5a_c$se <- psap_m5a_c$se.fit
nsap_m5a_c$lci <- nsap_m5a_c$fit - (nsap_m5a_c$se * 1.96)
nsap_m5a_c$uci <- nsap_m5a_c$fit + (nsap_m5a_c$se * 1.96)

psap_m5h_c <- predict(rsap_g5.2t, newdata = nsap_m5h_c, se.fit = T, type = 'response')
nsap_m5h_c$fit <- psap_m5h_c$fit
nsap_m5h_c$se <- psap_m5h_c$se.fit
nsap_m5h_c$lci <- nsap_m5h_c$fit - (nsap_m5h_c$se * 1.96)
nsap_m5h_c$uci <- nsap_m5h_c$fit + (nsap_m5h_c$se * 1.96)

# GAM bs = 'cc' for Temp.
nsap_m5l_cc <- nsap_m5l
nsap_m5a_cc <- nsap_m5a
nsap_m5h_cc <- nsap_m5h

psap_m5l_cc <- predict(rsap_g5.3t, newdata = nsap_m5l_cc, se.fit = T, type = 'response')
nsap_m5l_cc$fit <- psap_m5l_cc$fit
nsap_m5l_cc$se <- psap_m5l_cc$se.fit
nsap_m5l_cc$lci <- nsap_m5l_cc$fit - (nsap_m5l_cc$se * 1.96)
nsap_m5l_cc$uci <- nsap_m5l_cc$fit + (nsap_m5l_cc$se * 1.96)

psap_m5a_cc <- predict(rsap_g5.3t, newdata = nsap_m5a_cc, se.fit = T, type = 'response')
nsap_m5a_cc$fit <- psap_m5a_cc$fit
nsap_m5a_cc$se <- psap_m5a_cc$se.fit
nsap_m5a_cc$lci <- nsap_m5a_cc$fit - (nsap_m5a_cc$se * 1.96)
nsap_m5a_cc$uci <- nsap_m5a_cc$fit + (nsap_m5a_cc$se * 1.96)

psap_m5h_cc <- predict(rsap_g5.3t, newdata = nsap_m5h_cc, se.fit = T, type = 'response')
nsap_m5h_cc$fit <- psap_m5h_cc$fit
nsap_m5h_cc$se <- psap_m5h_cc$se.fit
nsap_m5h_cc$lci <- nsap_m5h_cc$fit - (nsap_m5h_cc$se * 1.96)
nsap_m5h_cc$uci <- nsap_m5h_cc$fit + (nsap_m5h_cc$se * 1.96)

# Plot predictions proportion saplings -----
# GLMER
plot(nsap_m5l_lm$r_Temp, nsap_m5l_lm$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, xlim = c(-1.7, 1.5), xaxt = "n")
axis(side = 1, at = seq(-1.7, 1.5, 0.4), labels = seq(385, 425, 5))
mtext(side = 1, expression(bold("Temperature")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(rsap_m5t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'r_Temp', data = "nsap_m5l_lm", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5a_lm$r_Temp, nsap_m5a_lm$fit, type = "l", col = "black")
pg.ci(x = 'r_Temp', data = "nsap_m5a_lm", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5h_lm$r_Temp, nsap_m5h_lm$fit, type = "l", col = "red")
pg.ci(x = 'r_Temp', data = "nsap_m5h_lm", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(nsap_m4l_lm$r_Precip, nsap_m4l_lm$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, xlim = c(-1.7, 1.5), xaxt = "n")
axis(side = 1, at = seq(-1.7, 1.5, 0.4), labels = seq(385, 425, 5))
mtext(side = 1, expression(bold("Preciperature")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(rsap_m4t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'r_Precip', data = "nsap_m4l_lm", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m4a_lm$r_Precip, nsap_m4a_lm$fit, type = "l", col = "black")
pg.ci(x = 'r_Precip', data = "nsap_m4a_lm", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m4h_lm$r_Precip, nsap_m4h_lm$fit, type = "l", col = "red")
pg.ci(x = 'r_Precip', data = "nsap_m4h_lm", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')




# GLMER and GAM comparisons
dev.new(width = 16, height = 12, noRStudioGD = T, dpi = 300)
par(mfrow = c(2,2), mar = c(6,6,2,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))

plot(nsap_m5l_lm$r_Temp, nsap_m5l_lm$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(-1.7, 1.5), xaxt = "n")
axis(side = 1, at = seq(-1.7, 1.5, 0.4), labels = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(rsap_m5t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'r_Temp', data = "nsap_m5l_lm", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5a_lm$r_Temp, nsap_m5a_lm$fit, type = "l", col = "black")
pg.ci(x = 'r_Temp', data = "nsap_m5a_lm", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5h_lm$r_Temp, nsap_m5h_lm$fit, type = "l", col = "red")
pg.ci(x = 'r_Temp', data = "nsap_m5h_lm", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) GLMER")), cex = 2)



par(xpd = NA)
legend(x = 5.9, y = 1, legend = c("1 year", "9.5 years", "25 years"), col = c("blue", "black", 'red'), title = expression(bold("Time since \n fire")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)




plot(nsap_m5l_d$Temp, nsap_m5l_d$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(rsap_g5t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nsap_m5l_d", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5a_d$Temp, nsap_m5a_d$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nsap_m5a_d", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5h_d$Temp, nsap_m5h_d$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nsap_m5h_d", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b) GAM defaults")), cex = 2)


plot(nsap_m5l_t$Temp, nsap_m5l_t$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4,  col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(rsap_g5.1t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nsap_m5l_t", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5a_t$Temp, nsap_m5a_t$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nsap_m5a_t", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5h_t$Temp, nsap_m5h_t$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nsap_m5h_t", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) GAM bs = 'tp'")), cex = 2)


plot(nsap_m5l_c$Temp, nsap_m5l_c$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(rsap_g5.2t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nsap_m5l_c", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5a_c$Temp, nsap_m5a_c$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nsap_m5a_c", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5h_c$Temp, nsap_m5h_c$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nsap_m5h_c", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(d) GAM bs = 'cr'")), cex = 2)


# Modelling as a cyclic cubic regression not useful for time since fire, just have flat responses. Unable to optimise model fit as too few data. 
plot(nsap_m5l_cc$Temp, nsap_m5l_cc$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(rsap_g5.3t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nsap_m5l_cc", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5a_cc$Temp, nsap_m5a_cc$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nsap_m5a_cc", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsap_m5h_cc$Temp, nsap_m5h_cc$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nsap_m5h_cc", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(d) GAM bs = 'cc'")), cex = 2)



# 5.1.3 Proportion recruits ----
recruit_tnull <- glmer(Proportion_recruits ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
recruit_gt_null <- gam(Proportion_recruits ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')

recruit_lnull <- glmer(Proportion_recruits ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)
recruit_gl_null <- gam(Proportion_recruits ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')



recruit_m1t <- glmer(Proportion_recruits ~ r_TSF + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m1t)
recruit_g1t <- gam(Proportion_recruits ~ s(TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g1t)
plot(recruit_g1t)
recruit_g1.1t <- gam(Proportion_recruits ~ s(TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g1.1t)
plot(recruit_g1.1t)

recruit_m1l <- glmer(Proportion_recruits ~ r_TSF + (1 | Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(recruit_m1l) # Not a very good model
recruit_g1l <- gam(Proportion_recruits ~ s(TSF, k = 4), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g1l)
plot(recruit_g1l)
recruit_g1.1l <- gam(Proportion_recruits ~ s(TSF, bs = 'cr', k = 4), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g1.1l)
plot(recruit_g1.1l)


recruit_m2t <- glmer(Proportion_recruits ~ r_TSF * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m2t)
recruit_g2t <- gam(Proportion_recruits ~ s(TSF) + s(Latitude) + ti(Latitude, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g2t)
plot(recruit_g2t)
recruit_g2.1t <- gam(Proportion_recruits ~ s(TSF) + s(Latitude) + ti(Latitude, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g2.1t)
plot(recruit_g2.1t)
recruit_g2.2t <- gam(Proportion_recruits ~ s(TSF, bs = 'cr', k = 5) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g2.1t)
plot(recruit_g2.2t)


#recruit_m2l <- glmer(Proportion_recruits ~ TSF * Latitude + (1|Location/Transect), family = binomial, data = lit_transects) # This model fails due to the low number of points available to model an interaction. We can only model an additive effect of TSF and latitude
recruit_m2.1l <- glmer(Proportion_recruits ~ r_TSF + r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(recruit_m2.1l) # Not very good model
recruit_g2l <- gam(Proportion_recruits ~ s(TSF, k = 4) + s(Latitude, k = 4) + ti(Latitude, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g2l)
plot(recruit_g2l)
recruit_g2.1l <- gam(Proportion_recruits ~ s(TSF, k = 4) + s(Latitude, k = 4) + ti(Latitude, by = TSF, bs = 'tp', k = 8), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g2.1l)
plot(recruit_g2.1l)
recruit_g2.2l <- gam(Proportion_recruits ~ s(TSF, bs = 'cr', k = 4) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g2.1l)
plot(recruit_g2.2l)



recruit_m3t <- glmer(Proportion_recruits ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
recruit_g3t <- gam(Proportion_recruits ~ s(TSF) + s(FPC) + ti(FPC, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g3t)
plot(recruit_g3t)
recruit_g3.1t <- gam(Proportion_recruits ~ s(TSF) + s(FPC) + ti(FPC, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g3.1t)
plot(recruit_g3.1t)
recruit_g3.2t <- gam(Proportion_recruits ~ s(TSF, bs = 'cr', k = 5) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g3.1t)
plot(recruit_g3.2t)


recruit_m3l <- glmer(Proportion_recruits ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(recruit_m3l) # Not a very good model
recruit_g3l <- gam(Proportion_recruits ~ s(TSF, k = 4) + s(FPC, k = 4) + ti(FPC, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g3l)
plot(recruit_g3l)
recruit_g3.1l <- gam(Proportion_recruits ~ s(TSF, k = 4) + s(FPC, k = 4) + ti(FPC, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g3.1l)
plot(recruit_g3.1l)
recruit_g3.2l <- gam(Proportion_recruits ~ s(TSF, bs = 'cr', k = 4) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g3.2l)
plot(recruit_g3.2l)


recruit_m4t <- glmer(Proportion_recruits ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
recruit_g4t <- gam(Proportion_recruits ~ s(TSF) + s(Precip) + ti(Precip, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g4t)
plot(recruit_g4t)
recruit_g4.1t <- gam(Proportion_recruits ~ s(TSF) + s(Precip) + ti(Precip, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g4.1t)
plot(recruit_g4.1t)
recruit_g4.2t <- gam(Proportion_recruits ~ s(TSF, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g4.1t)
plot(recruit_g4.2t)


recruit_m4l <- glmer(Proportion_recruits ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=250)))
summary(recruit_m4l)
recruit_g4l <- gam(Proportion_recruits ~ s(TSF, k = 4) + s(Precip, k = 4) + ti(Precip, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g4l)
plot(recruit_g4l)
recruit_g4.1l <- gam(Proportion_recruits ~ s(TSF, k = 4) + s(Precip, k = 4) + ti(Precip, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g4.1l)
plot(recruit_g4.1l)
recruit_g4.2l <- gam(Proportion_recruits ~ s(TSF, bs = 'cr', k = 4) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g4.2l)
plot(recruit_g4.2l)



recruit_m5t <- glmer(Proportion_recruits ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
recruit_g5t <- gam(Proportion_recruits ~ s(TSF) + s(Temp) + ti(Temp, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g5t)
plot(recruit_g5t)
recruit_g5.1t <- gam(Proportion_recruits ~ s(TSF) + s(Temp) + ti(Temp, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g5.1t)
plot(recruit_g5.1t)
recruit_g5.2t <- gam(Proportion_recruits ~ s(TSF, bs = 'cr', k = 5) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g5.1t)
plot(recruit_g5.2t)
recruit_g5.3t <- gam(Proportion_recruits ~ s(TSF, bs = 'cr', k = 5) + s(Temp, bs = 'cc', k = 30) + ti(Temp, by = TSF, bs = 'cc', k = 20), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g5.1t)
plot(recruit_g5.3t) # Not a good model, cannot optimise model fit by increasing or decreasing k for those terms including the cyclic cubic regression which is appropriate for climatic variables


recruit_m5l <- glmer(Proportion_recruits ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=250)))
summary(recruit_m5l) 
recruit_g5l <- gam(Proportion_recruits ~ s(TSF, k = 4) + s(Temp, k = 4) + ti(Temp, by = TSF), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g5l)
plot(recruit_g5l)
recruit_g5.1l <- gam(Proportion_recruits ~ s(TSF, k = 4) + s(Temp, k = 4) + ti(Temp, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g5.1l)
plot(recruit_g5.1l)
recruit_g5.2l <- gam(Proportion_recruits ~ s(TSF, bs = 'cr', k = 4) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recruit_g5.2l)
plot(recruit_g5.2l)

# Compare AICs for each modelling type
# GLMER
recruit_aic_t <- list(recruit_tnull, recruit_m1t, recruit_m2t, recruit_m3t, recruit_m4t, recruit_m5t)
aictab(recruit_aic_t) # Model 5 is best with no model ranked within AICc <2. 

recruit_aic_l <- list(recruit_lnull, recruit_m1l, recruit_m2.1l, recruit_m3l, recruit_m4l, recruit_m5t)
aictab(recruit_aic_l) # Model 2 is best but null model is ranked within delta AICc <2. 

# GAM defaults
recruit_td_aic <- as.data.frame(1:6)
recruit_td_aic$AICc <- "NA"
recruit_td_aic$Model <- "NA"
recruit_td_aic$LL <- "NA"
recruit_td_aic$AICc[1] <- AICc(recruit_gt_null)
recruit_td_aic$Model[1] <- "Null"
recruit_td_aic$LL[1] <- logLik(recruit_gt_null)
recruit_td_aic$AICc[2] <- AICc(recruit_g1t)
recruit_td_aic$Model[2] <- "m1"
recruit_td_aic$LL[2] <- logLik(recruit_g1t)
recruit_td_aic$AICc[3] <- AICc(recruit_g2t)
recruit_td_aic$Model[3] <- "m2"
recruit_td_aic$LL[3] <- logLik(recruit_g2t)
recruit_td_aic$AICc[4] <- AICc(recruit_g3t)
recruit_td_aic$Model[4] <- "m3"
recruit_td_aic$LL[4] <- logLik(recruit_g3t)
recruit_td_aic$AICc[5] <- AICc(recruit_g4t)
recruit_td_aic$Model[5] <- "m4"
recruit_td_aic$LL[5] <- logLik(recruit_g4t)
recruit_td_aic$AICc[6] <- AICc(recruit_g5t)
recruit_td_aic$Model[6] <- "m5"
recruit_td_aic$LL[6] <- logLik(recruit_g5t)
recruit_td_aic$AICc <- as.numeric(recruit_td_aic$AICc)
recruit_td_aic$LL <- as.numeric(recruit_td_aic$LL)
recruit_td_aic <- recruit_td_aic[order(recruit_td_aic$AICc), ]
recruit_td_aic$Delta_AICc <- "0.00"
recruit_td_aic$Delta_AICc[2] <- round(recruit_td_aic$AICc[1]-recruit_td_aic$AICc[2], 2)
recruit_td_aic$Delta_AICc[3] <- round(recruit_td_aic$AICc[1]-recruit_td_aic$AICc[3], 2)
recruit_td_aic$Delta_AICc[4] <- round(recruit_td_aic$AICc[1]-recruit_td_aic$AICc[4], 2)
recruit_td_aic$Delta_AICc[5] <- round(recruit_td_aic$AICc[1]-recruit_td_aic$AICc[5], 2)
recruit_td_aic$Delta_AICc[6] <- round(recruit_td_aic$AICc[1]-recruit_td_aic$AICc[6], 2)
recruit_td_aic  # Model 5 is best with no model ranked within delta AICc <2. 

recruit_ld_aic <- as.data.frame(1:6)
recruit_ld_aic$AICc <- "NA"
recruit_ld_aic$Model <- "NA"
recruit_ld_aic$LL <- "NA"
recruit_ld_aic$AICc[1] <- AICc(recruit_gl_null)
recruit_ld_aic$Model[1] <- "Null"
recruit_ld_aic$LL[1] <- logLik(recruit_gl_null)
recruit_ld_aic$AICc[2] <- AICc(recruit_g1l)
recruit_ld_aic$Model[2] <- "m1"
recruit_ld_aic$LL[2] <- logLik(recruit_g1l)
recruit_ld_aic$AICc[3] <- AICc(recruit_g2l)
recruit_ld_aic$Model[3] <- "m2"
recruit_ld_aic$LL[3] <- logLik(recruit_g2l)
recruit_ld_aic$AICc[4] <- AICc(recruit_g3l)
recruit_ld_aic$Model[4] <- "m3"
recruit_ld_aic$LL[4] <- logLik(recruit_g3l)
recruit_ld_aic$AICc[5] <- AICc(recruit_g4l)
recruit_ld_aic$Model[5] <- "m4"
recruit_ld_aic$LL[5] <- logLik(recruit_g4l)
recruit_ld_aic$AICc[6] <- AICc(recruit_g5l)
recruit_ld_aic$Model[6] <- "m5"
recruit_ld_aic$LL[6] <- logLik(recruit_g5l)
recruit_ld_aic$AICc <- as.numeric(recruit_ld_aic$AICc)
recruit_ld_aic$LL <- as.numeric(recruit_ld_aic$LL)
recruit_ld_aic <- recruit_ld_aic[order(recruit_ld_aic$AICc), ]
recruit_ld_aic$Delta_AICc <- "0.00"
recruit_ld_aic$Delta_AICc[2] <- round(recruit_ld_aic$AICc[1]-recruit_ld_aic$AICc[2], 2)
recruit_ld_aic$Delta_AICc[3] <- round(recruit_ld_aic$AICc[1]-recruit_ld_aic$AICc[3], 2)
recruit_ld_aic$Delta_AICc[4] <- round(recruit_ld_aic$AICc[1]-recruit_ld_aic$AICc[4], 2)
recruit_ld_aic$Delta_AICc[5] <- round(recruit_ld_aic$AICc[1]-recruit_ld_aic$AICc[5], 2)
recruit_ld_aic$Delta_AICc[6] <- round(recruit_ld_aic$AICc[1]-recruit_ld_aic$AICc[6], 2)
recruit_ld_aic  # Model 5 is best with no model ranked within delta AICc <2.



# GAM s(bs = 'tp')
recruit_ts_aic <- as.data.frame(1:6)
recruit_ts_aic$AICc <- "NA"
recruit_ts_aic$Model <- "NA"
recruit_ts_aic$LL <- "NA"
recruit_ts_aic$AICc[1] <- AICc(recruit_gt_null)
recruit_ts_aic$Model[1] <- "Null"
recruit_ts_aic$LL[1] <- logLik(recruit_gt_null)
recruit_ts_aic$AICc[2] <- AICc(recruit_g1t)
recruit_ts_aic$Model[2] <- "m1"
recruit_ts_aic$LL[2] <- logLik(recruit_g1t)
recruit_ts_aic$AICc[3] <- AICc(recruit_g2.1t)
recruit_ts_aic$Model[3] <- "m2"
recruit_ts_aic$LL[3] <- logLik(recruit_g2.1t)
recruit_ts_aic$AICc[4] <- AICc(recruit_g3.1t)
recruit_ts_aic$Model[4] <- "m3"
recruit_ts_aic$LL[4] <- logLik(recruit_g3.1t)
recruit_ts_aic$AICc[5] <- AICc(recruit_g4.1t)
recruit_ts_aic$Model[5] <- "m4"
recruit_ts_aic$LL[5] <- logLik(recruit_g4.1t)
recruit_ts_aic$AICc[6] <- AICc(recruit_g5.1t)
recruit_ts_aic$Model[6] <- "m5"
recruit_ts_aic$LL[6] <- logLik(recruit_g5.1t)
recruit_ts_aic$AICc <- as.numeric(recruit_ts_aic$AICc)
recruit_ts_aic$LL <- as.numeric(recruit_ts_aic$LL)
recruit_ts_aic <- recruit_ts_aic[order(recruit_ts_aic$AICc), ]
recruit_ts_aic$Delta_AICc <- "0.00"
recruit_ts_aic$Delta_AICc[2] <- round(recruit_ts_aic$AICc[1]-recruit_ts_aic$AICc[2], 2)
recruit_ts_aic$Delta_AICc[3] <- round(recruit_ts_aic$AICc[1]-recruit_ts_aic$AICc[3], 2)
recruit_ts_aic$Delta_AICc[4] <- round(recruit_ts_aic$AICc[1]-recruit_ts_aic$AICc[4], 2)
recruit_ts_aic$Delta_AICc[5] <- round(recruit_ts_aic$AICc[1]-recruit_ts_aic$AICc[5], 2)
recruit_ts_aic$Delta_AICc[6] <- round(recruit_ts_aic$AICc[1]-recruit_ts_aic$AICc[6], 2)
recruit_ts_aic  # Model 5 is best with no model ranked within delta AICc <2.

recruit_ls_aic <- as.data.frame(1:6)
recruit_ls_aic$AICc <- "NA"
recruit_ls_aic$Model <- "NA"
recruit_ls_aic$LL <- "NA"
recruit_ls_aic$AICc[1] <- AICc(recruit_gl_null)
recruit_ls_aic$Model[1] <- "Null"
recruit_ls_aic$LL[1] <- logLik(recruit_gl_null)
recruit_ls_aic$AICc[2] <- AICc(recruit_g1l)
recruit_ls_aic$Model[2] <- "m1"
recruit_ls_aic$LL[2] <- logLik(recruit_g1l)
recruit_ls_aic$AICc[3] <- AICc(recruit_g2.1l)
recruit_ls_aic$Model[3] <- "m2"
recruit_ls_aic$LL[3] <- logLik(recruit_g2.1l)
recruit_ls_aic$AICc[4] <- AICc(recruit_g3.1l)
recruit_ls_aic$Model[4] <- "m3"
recruit_ls_aic$LL[4] <- logLik(recruit_g3.1l)
recruit_ls_aic$AICc[5] <- AICc(recruit_g4.1l)
recruit_ls_aic$Model[5] <- "m4"
recruit_ls_aic$LL[5] <- logLik(recruit_g4.1l)
recruit_ls_aic$AICc[6] <- AICc(recruit_g5.1l)
recruit_ls_aic$Model[6] <- "m5"
recruit_ls_aic$LL[6] <- logLik(recruit_g5.1l)
recruit_ls_aic$AICc <- as.numeric(recruit_ls_aic$AICc)
recruit_ls_aic$LL <- as.numeric(recruit_ls_aic$LL)
recruit_ls_aic <- recruit_ls_aic[order(recruit_ls_aic$AICc), ]
recruit_ls_aic$Delta_AICc <- "0.00"
recruit_ls_aic$Delta_AICc[2] <- round(recruit_ls_aic$AICc[1]-recruit_ls_aic$AICc[2], 2)
recruit_ls_aic$Delta_AICc[3] <- round(recruit_ls_aic$AICc[1]-recruit_ls_aic$AICc[3], 2)
recruit_ls_aic$Delta_AICc[4] <- round(recruit_ls_aic$AICc[1]-recruit_ls_aic$AICc[4], 2)
recruit_ls_aic$Delta_AICc[5] <- round(recruit_ls_aic$AICc[1]-recruit_ls_aic$AICc[5], 2)
recruit_ls_aic$Delta_AICc[6] <- round(recruit_ls_aic$AICc[1]-recruit_ls_aic$AICc[6], 2)
recruit_ls_aic  #Null model is best but m1 is ranked within delta AICc <2.


# GAM ti(bs='cr')
recruit_tsc_aic <- as.data.frame(1:6)
recruit_tsc_aic$AICc <- "NA"
recruit_tsc_aic$Model <- "NA"
recruit_tsc_aic$LL <- "NA"
recruit_tsc_aic$AICc[1] <- AICc(recruit_gt_null)
recruit_tsc_aic$Model[1] <- "Null"
recruit_tsc_aic$LL[1] <- logLik(recruit_gt_null)
recruit_tsc_aic$AICc[2] <- AICc(recruit_g1.1t)
recruit_tsc_aic$Model[2] <- "m1"
recruit_tsc_aic$LL[2] <- logLik(recruit_g1.1t)
recruit_tsc_aic$AICc[3] <- AICc(recruit_g2.2t)
recruit_tsc_aic$Model[3] <- "m2"
recruit_tsc_aic$LL[3] <- logLik(recruit_g2.2t)
recruit_tsc_aic$AICc[4] <- AICc(recruit_g3.2t)
recruit_tsc_aic$Model[4] <- "m3"
recruit_tsc_aic$LL[4] <- logLik(recruit_g3.2t)
recruit_tsc_aic$AICc[5] <- AICc(recruit_g4.2t)
recruit_tsc_aic$Model[5] <- "m4"
recruit_tsc_aic$LL[5] <- logLik(recruit_g4.2t)
recruit_tsc_aic$AICc[6] <- AICc(recruit_g5.2t)
recruit_tsc_aic$Model[6] <- "m5"
recruit_tsc_aic$LL[6] <- logLik(recruit_g5.2t)
recruit_tsc_aic$AICc <- as.numeric(recruit_tsc_aic$AICc)
recruit_tsc_aic$LL <- as.numeric(recruit_tsc_aic$LL)
recruit_tsc_aic <- recruit_tsc_aic[order(recruit_tsc_aic$AICc), ]
recruit_tsc_aic$Delta_AICc <- "0.00"
recruit_tsc_aic$Delta_AICc[2] <- round(recruit_tsc_aic$AICc[1]-recruit_tsc_aic$AICc[2], 2)
recruit_tsc_aic$Delta_AICc[3] <- round(recruit_tsc_aic$AICc[1]-recruit_tsc_aic$AICc[3], 2)
recruit_tsc_aic$Delta_AICc[4] <- round(recruit_tsc_aic$AICc[1]-recruit_tsc_aic$AICc[4], 2)
recruit_tsc_aic$Delta_AICc[5] <- round(recruit_tsc_aic$AICc[1]-recruit_tsc_aic$AICc[5], 2)
recruit_tsc_aic$Delta_AICc[6] <- round(recruit_tsc_aic$AICc[1]-recruit_tsc_aic$AICc[6], 2)
recruit_tsc_aic  # Model 5 is best with no model ranked within delta AICc <2.

recruit_lsc_aic <- as.data.frame(1:6)
recruit_lsc_aic$AICc <- "NA"
recruit_lsc_aic$Model <- "NA"
recruit_lsc_aic$LL <- "NA"
recruit_lsc_aic$AICc[1] <- AICc(recruit_gl_null)
recruit_lsc_aic$Model[1] <- "Null"
recruit_lsc_aic$LL[1] <- logLik(recruit_gl_null)
recruit_lsc_aic$AICc[2] <- AICc(recruit_g1l)
recruit_lsc_aic$Model[2] <- "m1"
recruit_lsc_aic$LL[2] <- logLik(recruit_g1l)
recruit_lsc_aic$AICc[3] <- AICc(recruit_g2.2l)
recruit_lsc_aic$Model[3] <- "m2"
recruit_lsc_aic$LL[3] <- logLik(recruit_g2.2l)
recruit_lsc_aic$AICc[4] <- AICc(recruit_g3.2l)
recruit_lsc_aic$Model[4] <- "m3"
recruit_lsc_aic$LL[4] <- logLik(recruit_g3.2l)
recruit_lsc_aic$AICc[5] <- AICc(recruit_g4.2l)
recruit_lsc_aic$Model[5] <- "m4"
recruit_lsc_aic$LL[5] <- logLik(recruit_g4.2l)
recruit_lsc_aic$AICc[6] <- AICc(recruit_g5.2l)
recruit_lsc_aic$Model[6] <- "m5"
recruit_lsc_aic$LL[6] <- logLik(recruit_g5.2l)
recruit_lsc_aic$AICc <- as.numeric(recruit_lsc_aic$AICc)
recruit_lsc_aic$LL <- as.numeric(recruit_lsc_aic$LL)
recruit_lsc_aic <- recruit_lsc_aic[order(recruit_lsc_aic$AICc), ]
recruit_lsc_aic$Delta_AICc <- "0.00"
recruit_lsc_aic$Delta_AICc[2] <- round(recruit_lsc_aic$AICc[1]-recruit_lsc_aic$AICc[2], 2)
recruit_lsc_aic$Delta_AICc[3] <- round(recruit_lsc_aic$AICc[1]-recruit_lsc_aic$AICc[3], 2)
recruit_lsc_aic$Delta_AICc[4] <- round(recruit_lsc_aic$AICc[1]-recruit_lsc_aic$AICc[4], 2)
recruit_lsc_aic$Delta_AICc[5] <- round(recruit_lsc_aic$AICc[1]-recruit_lsc_aic$AICc[5], 2)
recruit_lsc_aic$Delta_AICc[6] <- round(recruit_lsc_aic$AICc[1]-recruit_lsc_aic$AICc[6], 2)
recruit_lsc_aic  #Null model is best but m1 is rakned within delta AICc <2.

# Predict for proportion recruits ----
# We can use the same dataframes that were produced previously for the GLMER and GAM model 

# GLMER
nrecruit_m5l_lm <- nsap_m5l_lm[, 1:2]
nrecruit_m5a_lm <- nsap_m5a_lm[, 1:2]
nrecruit_m5h_lm <- nsap_m5h_lm[, 1:2]

precruit_m5l_lm <- predictSE(recruit_m5t, newdata = nsap_m5l_lm, se.fit = T, type = 'link')
nrecruit_m5l_lm$fit.link <- precruit_m5l_lm$fit
nrecruit_m5l_lm$se.link <- precruit_m5l_lm$se.fit
nrecruit_m5l_lm$lci.link <- nrecruit_m5l_lm$fit.link - (nrecruit_m5l_lm$se.link * 1.96)
nrecruit_m5l_lm$uci.link <- nrecruit_m5l_lm$fit.link + (nrecruit_m5l_lm$se.link * 1.96)
nrecruit_m5l_lm$fit <- invlogit(nrecruit_m5l_lm$fit.link)
nrecruit_m5l_lm$se <- invlogit(nrecruit_m5l_lm$se.link)
nrecruit_m5l_lm$lci <- invlogit(nrecruit_m5l_lm$lci.link)
nrecruit_m5l_lm$uci <- invlogit(nrecruit_m5l_lm$uci.link)

precruit_m5a_lm <- predictSE(recruit_m5t, newdata = nsap_m5a_lm, se.fit = T, type = 'link')
nrecruit_m5a_lm$fit.link <- precruit_m5a_lm$fit
nrecruit_m5a_lm$se.link <- precruit_m5a_lm$se.fit
nrecruit_m5a_lm$lci.link <- nrecruit_m5a_lm$fit.link - (nrecruit_m5a_lm$se.link * 1.96)
nrecruit_m5a_lm$uci.link <- nrecruit_m5a_lm$fit.link + (nrecruit_m5a_lm$se.link * 1.96)
nrecruit_m5a_lm$fit <- invlogit(nrecruit_m5a_lm$fit.link)
nrecruit_m5a_lm$se <- invlogit(nrecruit_m5a_lm$se.link)
nrecruit_m5a_lm$lci <- invlogit(nrecruit_m5a_lm$lci.link)
nrecruit_m5a_lm$uci <- invlogit(nrecruit_m5a_lm$uci.link)

precruit_m5h_lm <- predictSE(recruit_m5t, newdata = nsap_m5h_lm, se.fit = T, type = 'link')
nrecruit_m5h_lm$fit.link <- precruit_m5h_lm$fit
nrecruit_m5h_lm$se.link <- precruit_m5h_lm$se.fit
nrecruit_m5h_lm$lci.link <- nrecruit_m5h_lm$fit.link - (nrecruit_m5h_lm$se.link * 1.96)
nrecruit_m5h_lm$uci.link <- nrecruit_m5h_lm$fit.link + (nrecruit_m5h_lm$se.link * 1.96)
nrecruit_m5h_lm$fit <- invlogit(nrecruit_m5h_lm$fit.link)
nrecruit_m5h_lm$se <- invlogit(nrecruit_m5h_lm$se.link)
nrecruit_m5h_lm$lci <- invlogit(nrecruit_m5h_lm$lci.link)
nrecruit_m5h_lm$uci <- invlogit(nrecruit_m5h_lm$uci.link)






# GAM default
nrecruit_m5l_d <- nsap_m5l
nrecruit_m5a_d <- nsap_m5a
nrecruit_m5h_d <- nsap_m5h

precruit_m5l_d <- predict(recruit_g5t, newdata = nsap_m5l_d, se.fit = T, type = 'response')
nrecruit_m5l_d$fit <- precruit_m5l_d$fit
nrecruit_m5l_d$se <- precruit_m5l_d$se.fit
nrecruit_m5l_d$lci <- nrecruit_m5l_d$fit - (nrecruit_m5l_d$se * 1.96)
nrecruit_m5l_d$uci <- nrecruit_m5l_d$fit + (nrecruit_m5l_d$se * 1.96)


precruit_m5a_d <- predict(recruit_g5t, newdata = nsap_m5a_d, se.fit = T, type = 'response')
nrecruit_m5a_d$fit <- precruit_m5a_d$fit
nrecruit_m5a_d$se <- precruit_m5a_d$se.fit
nrecruit_m5a_d$lci <- nrecruit_m5a_d$fit - (nrecruit_m5a_d$se * 1.96)
nrecruit_m5a_d$uci <- nrecruit_m5a_d$fit + (nrecruit_m5a_d$se * 1.96)

precruit_m5h_d <- predict(recruit_g5t, newdata = nsap_m5h_d, se.fit = T, type = 'response')
nrecruit_m5h_d$fit <- precruit_m5h_d$fit
nrecruit_m5h_d$se <- precruit_m5h_d$se.fit
nrecruit_m5h_d$lci <- nrecruit_m5h_d$fit - (nrecruit_m5h_d$se * 1.96)
nrecruit_m5h_d$uci <- nrecruit_m5h_d$fit + (nrecruit_m5h_d$se * 1.96)


# GAM bs = 'tp'
nrecruit_m5l_t <- nsap_m5l
nrecruit_m5a_t <- nsap_m5a
nrecruit_m5h_t <- nsap_m5h

precruit_m5l_t <- predict(recruit_g5.1t, newdata = nsap_m5l_t, se.fit = T, type = 'response')
nrecruit_m5l_t$fit <- precruit_m5l_t$fit
nrecruit_m5l_t$se <- precruit_m5l_t$se.fit
nrecruit_m5l_t$lci <- nrecruit_m5l_t$fit - (nrecruit_m5l_t$se * 1.96)
nrecruit_m5l_t$uci <- nrecruit_m5l_t$fit + (nrecruit_m5l_t$se * 1.96)

precruit_m5a_t <- predict(recruit_g5.1t, newdata = nsap_m5a_t, se.fit = T, type = 'response')
nrecruit_m5a_t$fit <- precruit_m5a_t$fit
nrecruit_m5a_t$se <- precruit_m5a_t$se.fit
nrecruit_m5a_t$lci <- nrecruit_m5a_t$fit - (nrecruit_m5a_t$se * 1.96)
nrecruit_m5a_t$uci <- nrecruit_m5a_t$fit + (nrecruit_m5a_t$se * 1.96)

precruit_m5h_t <- predict(recruit_g5.1t, newdata = nsap_m5h_t, se.fit = T, type = 'response')
nrecruit_m5h_t$fit <- precruit_m5h_t$fit
nrecruit_m5h_t$se <- precruit_m5h_t$se.fit
nrecruit_m5h_t$lci <- nrecruit_m5h_t$fit - (nrecruit_m5h_t$se * 1.96)
nrecruit_m5h_t$uci <- nrecruit_m5h_t$fit + (nrecruit_m5h_t$se * 1.96)


# GAM bs = 'cr'
nrecruit_m5l_c <- nsap_m5l
nrecruit_m5a_c <- nsap_m5a
nrecruit_m5h_c <- nsap_m5h

precruit_m5l_c <- predict(recruit_g5.2t, newdata = nsap_m5l_c, se.fit = T, type = 'response')
nrecruit_m5l_c$fit <- precruit_m5l_c$fit
nrecruit_m5l_c$se <- precruit_m5l_c$se.fit
nrecruit_m5l_c$lci <- nrecruit_m5l_c$fit - (nrecruit_m5l_c$se * 1.96)
nrecruit_m5l_c$uci <- nrecruit_m5l_c$fit + (nrecruit_m5l_c$se * 1.96)

precruit_m5a_c <- predict(recruit_g5.2t, newdata = nsap_m5a_c, se.fit = T, type = 'response')
nrecruit_m5a_c$fit <- precruit_m5a_c$fit
nrecruit_m5a_c$se <- precruit_m5a_c$se.fit
nrecruit_m5a_c$lci <- nrecruit_m5a_c$fit - (nrecruit_m5a_c$se * 1.96)
nrecruit_m5a_c$uci <- nrecruit_m5a_c$fit + (nrecruit_m5a_c$se * 1.96)

precruit_m5h_c <- predict(recruit_g5.2t, newdata = nsap_m5h_c, se.fit = T, type = 'response')
nrecruit_m5h_c$fit <- precruit_m5h_c$fit
nrecruit_m5h_c$se <- precruit_m5h_c$se.fit
nrecruit_m5h_c$lci <- nrecruit_m5h_c$fit - (nrecruit_m5h_c$se * 1.96)
nrecruit_m5h_c$uci <- nrecruit_m5h_c$fit + (nrecruit_m5h_c$se * 1.96)



# Plot predictions proportion recruits -----
dev.new(width = 16, height = 12, noRStudioGD = T, dpi = 300)
par(mfrow = c(2,2), mar = c(6,6,2,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))

plot(nrecruit_m5l_lm$r_Temp, nrecruit_m5l_lm$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(-1.7, 1.5), xaxt = "n")
axis(side = 1, at = seq(-1.7, 1.5, 0.4), labels = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(recruit_m5t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'r_Temp', data = "nrecruit_m5l_lm", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nrecruit_m5a_lm$r_Temp, nrecruit_m5a_lm$fit, type = "l", col = "black")
pg.ci(x = 'r_Temp', data = "nrecruit_m5a_lm", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nrecruit_m5h_lm$r_Temp, nrecruit_m5h_lm$fit, type = "l", col = "red")
pg.ci(x = 'r_Temp', data = "nrecruit_m5h_lm", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) GLMER")), cex = 2)



par(xpd = NA)
legend(x = 5.9, y = 1, legend = c("1 year", "9.5 years", "25 years"), col = c("blue", "black", 'red'), title = expression(bold("Time since \n fire")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)




plot(nrecruit_m5l_d$Temp, nrecruit_m5l_d$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(recruit_g5t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nrecruit_m5l_d", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nrecruit_m5a_d$Temp, nrecruit_m5a_d$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nrecruit_m5a_d", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nrecruit_m5h_d$Temp, nrecruit_m5h_d$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nrecruit_m5h_d", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b) GAM defaults")), cex = 2)


plot(nrecruit_m5l_t$Temp, nrecruit_m5l_t$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4,  col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(recruit_g5.1t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nrecruit_m5l_t", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nrecruit_m5a_t$Temp, nrecruit_m5a_t$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nrecruit_m5a_t", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nrecruit_m5h_t$Temp, nrecruit_m5h_t$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nrecruit_m5h_t", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) GAM bs = 'tp'")), cex = 2)


plot(nrecruit_m5l_c$Temp, nrecruit_m5l_c$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(recruit_g5.2t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nrecruit_m5l_c", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nrecruit_m5a_c$Temp, nrecruit_m5a_c$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nrecruit_m5a_c", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nrecruit_m5h_c$Temp, nrecruit_m5h_c$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nrecruit_m5h_c", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(d) GAM bs = 'cr'")), cex = 2)






# 5.2 Female fecundity ----
# This should use the tor_tree and lit_tree data as this will give us more information on how height variation influences responses
### For fecundity only - Maybe we are using the tree level data, then we nest individual in transect in location? Otherwise we can't really investigate height 
# TSF 
# TSF + height 
# TSF * latitude
# TSF * FPC
# TSF * precipitation seasonality
# TSF * temperature seasonality

tor_fecundity <- tor_tree[which(tor_tree$Cone_presence == "Y"),]
lit_fecundity <- lit_tree[which(lit_tree$Cone_presence == "Y"),] # Note there are some NAs present as the Gatton transects did not record the number of cones if they were present.

tor_fecundity$Height_cm <- tor_fecundity$Height_.m.*100
lit_fecundity$Height_cm <- lit_fecundity$Height_.m.*100

str(tor_fecundity)

tor_fecundity$r_TSF <- scale(tor_fecundity$TSF)
tor_fecundity$r_Latitude <- scale(tor_fecundity$Latitude)
tor_fecundity$r_height <- scale(tor_fecundity$Height_cm)
tor_fecundity$r_FPC <- scale(tor_fecundity$FPC)
tor_fecundity$r_Precip <- scale(tor_fecundity$Precip)
tor_fecundity$r_Temp <- scale(tor_fecundity$Temp)

lit_fecundity$r_TSF <- scale(lit_fecundity$TSF)
lit_fecundity$r_Latitude <- scale(lit_fecundity$Latitude)
lit_fecundity$r_height <- scale(lit_fecundity$Height_cm)
lit_fecundity$r_FPC <- scale(lit_fecundity$FPC)
lit_fecundity$r_Precip <- scale(lit_fecundity$Precip)
lit_fecundity$r_Temp <- scale(lit_fecundity$Temp)


fecundity_tnull <- glmer(Cone_number ~ 1 + (1 | Location/Transect), family = poisson, data = tor_fecundity)
fecundity_gt_null <- gam(Cone_number ~ 1, random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')

fecundity_lnull <- glmer(Cone_number ~ 1 + (1 | Location/Transect), family = poisson, data = lit_fecundity)
fecundity_gl_null <- gam(Cone_number ~ 1, random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')



fecundity_m1t <- glmer(Cone_number ~ r_TSF + (1 | Location/Transect), family = poisson, data = tor_fecundity)
summary(fecundity_m1t)
fecundity_g1t <- gam(Cone_number ~ s(TSF, k = 9), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g1t)
plot(fecundity_g1t)
fecundity_g1.1t <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g1.1t)
plot(fecundity_g1.1t)
fecundity_g1ot <- gam(Cone_number ~ s(TSF, k = 4), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g1ot) # k' and edf are always close but this gives us a nice smooth curve
plot(fecundity_g1ot)
fecundity_g1.1ot <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 4), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g1.1ot)
plot(fecundity_g1.1ot)


fecundity_m1l <- glmer(Cone_number ~ r_TSF + (1 | Location/Transect), family = poisson, data = lit_fecundity)
summary(fecundity_m1l) 
fecundity_g1l <- gam(Cone_number ~ s(TSF, k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g1l)
plot(fecundity_g1l)
fecundity_g1.1l <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g1.1l)
plot(fecundity_g1.1l)
fecundity_g1ol <- gam(Cone_number ~ s(TSF, k = 4), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g1ol) # k' and edf are always close but this gives us a nice smooth curve
plot(fecundity_g1ol)
fecundity_g1.1ol <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 4), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g1.1ol)
plot(fecundity_g1.1ol)



fecundity_m2t <- glmer(Cone_number ~ r_TSF * r_Latitude + (1|Location/Transect), family = poisson, data = tor_fecundity)
summary(fecundity_m2t)
fecundity_g2t <- gam(Cone_number ~ s(TSF, k = 9) + s(Latitude) + ti(Latitude, by = TSF), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g2t)
plot(fecundity_g2t)
fecundity_g2.1t <- gam(Cone_number ~ s(TSF, k = 9) + s(Latitude) + ti(Latitude, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g2.1t)
plot(fecundity_g2.1t)
fecundity_g2.2t <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g2.2t)
plot(fecundity_g2.2t)
fecundity_g2.1ot <- gam(Cone_number ~ s(TSF, k = 5) + s(Latitude, k = 5) + ti(Latitude, by = TSF, bs = 'tp', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g2.1ot)
plot(fecundity_g2.1ot)
fecundity_g2.2ot <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 4) + s(Latitude, bs = 'cr', k = 4) + ti(Latitude, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g2.2ot)
plot(fecundity_g2.2ot)

fecundity_m2l <- glmer(Cone_number ~ TSF * Latitude + (1|Location/Transect), family = poisson, data = lit_fecundity) 
summary(fecundity_m2l)
fecundity_g2l <- gam(Cone_number ~ s(TSF, k = 5) + s(Latitude) + ti(Latitude, by = TSF), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g2l)
plot(fecundity_g2l)
fecundity_g2.1l <- gam(Cone_number ~ s(TSF, k = 5) + s(Latitude) + ti(Latitude, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g2.1l)
plot(fecundity_g2.1l)
fecundity_g2.2l <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g2.2l)
plot(fecundity_g2.2l)



fecundity_m3t <- glmer(Cone_number ~ r_TSF * r_FPC + (1|Location/Transect), family = poisson, data = tor_fecundity)
fecundity_g3t <- gam(Cone_number ~ s(TSF, k = 9) + s(FPC) + ti(FPC, by = TSF), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g3t)
plot(fecundity_g3t)
fecundity_g3.1t <- gam(Cone_number ~ s(TSF, k = 9) + s(FPC) + ti(FPC, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g3.1t)
plot(fecundity_g3.1t)
fecundity_g3.2t <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g3.2t)
plot(fecundity_g3.2t)
fecundity_g3.1ot <- gam(Cone_number ~ s(TSF, k = 4) + s(FPC, k = 4) + ti(FPC, by = TSF, bs = 'tp', k = 13), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g3.1ot)
plot(fecundity_g3.1ot)
fecundity_g3.2ot <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 4) + s(FPC, bs = 'cr', k = 4) + ti(FPC, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g3.2ot)
plot(fecundity_g3.2ot)


fecundity_m3l <- glmer(Cone_number ~ r_TSF * r_FPC + (1|Location/Transect), family = poisson, data = lit_fecundity)
summary(fecundity_m3l) 
fecundity_g3l <- gam(Cone_number ~ s(TSF, k = 4) + s(FPC) + ti(FPC, by = TSF), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g3l)
plot(fecundity_g3l)
fecundity_g3.1l <- gam(Cone_number ~ s(TSF, k = 4) + s(FPC) + ti(FPC, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g3.1l)
plot(fecundity_g3.1l)
fecundity_g3.2l <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 4) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g3.2l)
plot(fecundity_g3.2l)
fecundity_g3.1ol <- gam(Cone_number ~ s(TSF, k = 4) + s(FPC, k = 8) + ti(FPC, by = TSF, bs = 'tp', k = 4), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g3.1ol)
plot(fecundity_g3.1ol)
fecundity_g3.2ol <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 4) + s(FPC, bs = 'cr', k = 4) + ti(FPC, by = TSF, bs = 'cr', k = 4), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g3.2ol)
plot(fecundity_g3.2ol)

fecundity_m4t <- glmer(Cone_number ~ r_TSF * r_Precip + (1|Location/Transect), family = poisson, data = tor_fecundity)
fecundity_g4t <- gam(Cone_number ~ s(TSF, k = 9) + s(Precip) + ti(Precip, by = TSF), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g4t)
plot(fecundity_g4t)
fecundity_g4.1t <- gam(Cone_number ~ s(TSF, k = 9) + s(Precip) + ti(Precip, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g4.1t)
plot(fecundity_g4.1t)
fecundity_g4.2t <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g4.2t)
plot(fecundity_g4.2t)
fecundity_g4.1ot <- gam(Cone_number ~ s(TSF, k = 5) + s(Precip, k = 4) + ti(Precip, by = TSF, bs = 'tp', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g4.1ot)
plot(fecundity_g4.1ot)
fecundity_g4.2ot <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 4) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g4.2ot)
plot(fecundity_g4.2ot)


fecundity_m4l <- glmer(Cone_number ~ r_TSF * r_Precip + (1|Location/Transect), family = poisson, data = lit_fecundity)
summary(fecundity_m4l)
fecundity_g4l <- gam(Cone_number ~ s(TSF, k = 5) + s(Precip, k = 5) + ti(Precip, by = TSF), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g4l)
plot(fecundity_g4l)
fecundity_g4.1l <- gam(Cone_number ~ s(TSF, k = 5) + s(Precip, k = 5) + ti(Precip, by = TSF, bs = 'tp', k = 6), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g4.1l)
plot(fecundity_g4.1l)
fecundity_g4.2l <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g4.2l)
plot(fecundity_g4.2l)
fecundity_g4.1ol <- gam(Cone_number ~ s(TSF, k = 3) + s(Precip, k = 3) + ti(Precip, by = TSF, bs = 'tp', k = 3), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g4.1ol)
plot(fecundity_g4.1ol)
fecundity_g4.2ol <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 4) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g4.2ol)
plot(fecundity_g4.2ol)


fecundity_m5t <- glmer(Cone_number ~ r_TSF * r_Precip + (1|Location/Transect), family = poisson, data = tor_fecundity)
fecundity_g5t <- gam(Cone_number ~ s(TSF, k = 9) + s(Precip) + ti(Precip, by = TSF), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g5t)
plot(fecundity_g5t)
fecundity_g5.1t <- gam(Cone_number ~ s(TSF, k = 9) + s(Precip) + ti(Precip, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g5.1t)
plot(fecundity_g5.1t)
fecundity_g5.2t <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g5.2t)
plot(fecundity_g5.2t)
fecundity_g5.1ot <- gam(Cone_number ~ s(TSF, k = 4) + s(Precip, k = 4) + ti(Precip, by = TSF, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g5.1ot)
plot(fecundity_g5.1ot)
fecundity_g5.2ot <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 3) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g5.2ot)
plot(fecundity_g5.2ot)


fecundity_m5l <- glmer(Cone_number ~ r_TSF * r_Precip + (1|Location/Transect), family = poisson, data = lit_fecundity)
summary(fecundity_m5l)
fecundity_g5l <- gam(Cone_number ~ s(TSF, k = 5) + s(Precip, k = 5) + ti(Precip, by = TSF), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g5l)
plot(fecundity_g5l)
fecundity_g5.1l <- gam(Cone_number ~ s(TSF, k = 5) + s(Precip, k = 5) + ti(Precip, by = TSF, bs = 'tp', k = 6), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g5.1l)
plot(fecundity_g5.1l)
fecundity_g5.2l <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g5.2l)
plot(fecundity_g5.2l)
fecundity_g5.1ol <- gam(Cone_number ~ s(TSF, k = 5) + s(Precip, k = 4) + ti(Precip, by = TSF, bs = 'tp', k = 4), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g5.1ol)
plot(fecundity_g5.1ol)
fecundity_g5.2ol <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = TSF, bs = 'cr', k = 3), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g5.2ol)
plot(fecundity_g5.2ol)



fecundity_m6t <- glmer(Cone_number ~ r_TSF + r_height + (1|Location/Transect), family = poisson, data = tor_fecundity)
fecundity_g6t <- gam(Cone_number ~ s(TSF, k = 9) + s(Height_cm), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g6t)
plot(fecundity_g6t)
fecundity_g6.1t <- gam(Cone_number ~ s(TSF, k = 9) + s(Height_cm) + ti(Height_cm, by = TSF, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g6.1t)
plot(fecundity_g6.1t)
fecundity_g6.2t <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5) + s(Height_cm, bs = 'cr', k = 5) + ti(Height_cm, by = TSF, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g6.2t)
plot(fecundity_g6.2t)
fecundity_g6.1ot <- gam(Cone_number ~ s(TSF, k = 5) + s(Height_cm, k = 5) + ti(Height_cm, by = TSF, bs = 'tp', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g6.1ot)
plot(fecundity_g6.1ot)



fecundity_m6l <- glmer(Cone_number ~ r_TSF + r_height + (1|Location/Transect), family = poisson, data = lit_fecundity)
summary(fecundity_m6l)
fecundity_g6l <- gam(Cone_number ~ s(TSF, k = 5) + s(Height_cm), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g6l)
plot(fecundity_g6l)
fecundity_g6.1l <- gam(Cone_number ~ s(TSF, k = 5) + s(Height_cm), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g6.1l)
plot(fecundity_g6.1l)
fecundity_g6.2l <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 5) + s(Height_cm, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g6.2l)
plot(fecundity_g6.2l)
fecundity_g6.1ol <- gam(Cone_number ~ s(TSF, k = 5) + s(Height_cm, k = 4), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g6.1ol)
plot(fecundity_g6.1ol)
fecundity_g6.2ol <- gam(Cone_number ~ s(TSF, bs = 'cr', k = 4) + s(Height_cm, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundity_g6.2ol)
plot(fecundity_g6.2ol)



# Compare AICs for each modelling type
# GLMER
fecundity_aic_t <- list(fecundity_tnull, fecundity_m1t, fecundity_m2t, fecundity_m3t, fecundity_m4t, fecundity_m5t)
aictab(fecundity_aic_t) # Model 3 is best with TSF and FPC

fecundity_aic_l <- list(fecundity_lnull, fecundity_m1l, fecundity_m2l, fecundity_m3l, fecundity_m4l, fecundity_m5l)
aictab(fecundity_aic_l) # Model 3 is best with TSF and FPC

# GAM defaults
fecundity_td_aic <- as.data.frame(1:7)
fecundity_td_aic$AICc <- "NA"
fecundity_td_aic$Model <- "NA"
fecundity_td_aic$LL <- "NA"
fecundity_td_aic$AICc[1] <- AICc(fecundity_gt_null)
fecundity_td_aic$Model[1] <- "Null"
fecundity_td_aic$LL[1] <- logLik(fecundity_gt_null)
fecundity_td_aic$AICc[2] <- AICc(fecundity_g1t)
fecundity_td_aic$Model[2] <- "m1"
fecundity_td_aic$LL[2] <- logLik(fecundity_g1t)
fecundity_td_aic$AICc[3] <- AICc(fecundity_g2t)
fecundity_td_aic$Model[3] <- "m2"
fecundity_td_aic$LL[3] <- logLik(fecundity_g2t)
fecundity_td_aic$AICc[4] <- AICc(fecundity_g3t)
fecundity_td_aic$Model[4] <- "m3"
fecundity_td_aic$LL[4] <- logLik(fecundity_g3t)
fecundity_td_aic$AICc[5] <- AICc(fecundity_g4t)
fecundity_td_aic$Model[5] <- "m4"
fecundity_td_aic$LL[5] <- logLik(fecundity_g4t)
fecundity_td_aic$AICc[6] <- AICc(fecundity_g5t)
fecundity_td_aic$Model[6] <- "m5"
fecundity_td_aic$LL[6] <- logLik(fecundity_g5t)
fecundity_td_aic$AICc[7] <- AICc(fecundity_g6t)
fecundity_td_aic$Model[7] <- "m6"
fecundity_td_aic$LL[7] <- logLik(fecundity_g6t)
fecundity_td_aic$AICc <- as.numeric(fecundity_td_aic$AICc)
fecundity_td_aic$LL <- as.numeric(fecundity_td_aic$LL)
fecundity_td_aic <- fecundity_td_aic[order(fecundity_td_aic$AICc), ]
fecundity_td_aic$Delta_AICc <- "0.00"
fecundity_td_aic$Delta_AICc[2] <- round(fecundity_td_aic$AICc[1]-fecundity_td_aic$AICc[2], 2)
fecundity_td_aic$Delta_AICc[3] <- round(fecundity_td_aic$AICc[1]-fecundity_td_aic$AICc[3], 2)
fecundity_td_aic$Delta_AICc[4] <- round(fecundity_td_aic$AICc[1]-fecundity_td_aic$AICc[4], 2)
fecundity_td_aic$Delta_AICc[5] <- round(fecundity_td_aic$AICc[1]-fecundity_td_aic$AICc[5], 2)
fecundity_td_aic$Delta_AICc[6] <- round(fecundity_td_aic$AICc[1]-fecundity_td_aic$AICc[6], 2)
fecundity_td_aic$Delta_AICc[7] <- round(fecundity_td_aic$AICc[1]-fecundity_td_aic$AICc[7], 2)
fecundity_td_aic  # Model 3 is best.

fecundity_ld_aic <- as.data.frame(1:7)
fecundity_ld_aic$AICc <- "NA"
fecundity_ld_aic$Model <- "NA"
fecundity_ld_aic$LL <- "NA"
fecundity_ld_aic$AICc[1] <- AICc(fecundity_gl_null)
fecundity_ld_aic$Model[1] <- "Null"
fecundity_ld_aic$LL[1] <- logLik(fecundity_gl_null)
fecundity_ld_aic$AICc[2] <- AICc(fecundity_g1l)
fecundity_ld_aic$Model[2] <- "m1"
fecundity_ld_aic$LL[2] <- logLik(fecundity_g1l)
fecundity_ld_aic$AICc[3] <- AICc(fecundity_g2l)
fecundity_ld_aic$Model[3] <- "m2"
fecundity_ld_aic$LL[3] <- logLik(fecundity_g2l)
fecundity_ld_aic$AICc[4] <- AICc(fecundity_g3l)
fecundity_ld_aic$Model[4] <- "m3"
fecundity_ld_aic$LL[4] <- logLik(fecundity_g3l)
fecundity_ld_aic$AICc[5] <- AICc(fecundity_g4l)
fecundity_ld_aic$Model[5] <- "m4"
fecundity_ld_aic$LL[5] <- logLik(fecundity_g4l)
fecundity_ld_aic$AICc[6] <- AICc(fecundity_g5l)
fecundity_ld_aic$Model[6] <- "m5"
fecundity_ld_aic$LL[6] <- logLik(fecundity_g5l)
fecundity_ld_aic$AICc[7] <- AICc(fecundity_g6l)
fecundity_ld_aic$Model[7] <- "m6"
fecundity_ld_aic$LL[7] <- logLik(fecundity_g6l)
fecundity_ld_aic$AICc <- as.numeric(fecundity_ld_aic$AICc)
fecundity_ld_aic$LL <- as.numeric(fecundity_ld_aic$LL)
fecundity_ld_aic <- fecundity_ld_aic[order(fecundity_ld_aic$AICc), ]
fecundity_ld_aic$Delta_AICc <- "0.00"
fecundity_ld_aic$Delta_AICc[2] <- round(fecundity_ld_aic$AICc[1]-fecundity_ld_aic$AICc[2], 2)
fecundity_ld_aic$Delta_AICc[3] <- round(fecundity_ld_aic$AICc[1]-fecundity_ld_aic$AICc[3], 2)
fecundity_ld_aic$Delta_AICc[4] <- round(fecundity_ld_aic$AICc[1]-fecundity_ld_aic$AICc[4], 2)
fecundity_ld_aic$Delta_AICc[5] <- round(fecundity_ld_aic$AICc[1]-fecundity_ld_aic$AICc[5], 2)
fecundity_ld_aic$Delta_AICc[6] <- round(fecundity_ld_aic$AICc[1]-fecundity_ld_aic$AICc[6], 2)
fecundity_ld_aic$Delta_AICc[7] <- round(fecundity_ld_aic$AICc[1]-fecundity_ld_aic$AICc[7], 2)
fecundity_ld_aic  #Model 6 is best



# GAM s(bs = 'tp')
fecundity_ts_aic <- as.data.frame(1:7)
fecundity_ts_aic$AICc <- "NA"
fecundity_ts_aic$Model <- "NA"
fecundity_ts_aic$LL <- "NA"
fecundity_ts_aic$AICc[1] <- AICc(fecundity_gt_null)
fecundity_ts_aic$Model[1] <- "Null"
fecundity_ts_aic$LL[1] <- logLik(fecundity_gt_null)
fecundity_ts_aic$AICc[2] <- AICc(fecundity_g1t)
fecundity_ts_aic$Model[2] <- "m1"
fecundity_ts_aic$LL[2] <- logLik(fecundity_g1t)
fecundity_ts_aic$AICc[3] <- AICc(fecundity_g2.1t)
fecundity_ts_aic$Model[3] <- "m2"
fecundity_ts_aic$LL[3] <- logLik(fecundity_g2.1t)
fecundity_ts_aic$AICc[4] <- AICc(fecundity_g3.1t)
fecundity_ts_aic$Model[4] <- "m3"
fecundity_ts_aic$LL[4] <- logLik(fecundity_g3.1t)
fecundity_ts_aic$AICc[5] <- AICc(fecundity_g4.1t)
fecundity_ts_aic$Model[5] <- "m4"
fecundity_ts_aic$LL[5] <- logLik(fecundity_g4.1t)
fecundity_ts_aic$AICc[6] <- AICc(fecundity_g5.1t)
fecundity_ts_aic$Model[6] <- "m5"
fecundity_ts_aic$LL[6] <- logLik(fecundity_g5.1t)
fecundity_ts_aic$AICc[7] <- AICc(fecundity_g6.1t)
fecundity_ts_aic$Model[7] <- "m6"
fecundity_ts_aic$LL[7] <- logLik(fecundity_g6.1t)
fecundity_ts_aic$AICc <- as.numeric(fecundity_ts_aic$AICc)
fecundity_ts_aic$LL <- as.numeric(fecundity_ts_aic$LL)
fecundity_ts_aic <- fecundity_ts_aic[order(fecundity_ts_aic$AICc), ]
fecundity_ts_aic$Delta_AICc <- "0.00"
fecundity_ts_aic$Delta_AICc[2] <- round(fecundity_ts_aic$AICc[1]-fecundity_ts_aic$AICc[2], 2)
fecundity_ts_aic$Delta_AICc[3] <- round(fecundity_ts_aic$AICc[1]-fecundity_ts_aic$AICc[3], 2)
fecundity_ts_aic$Delta_AICc[4] <- round(fecundity_ts_aic$AICc[1]-fecundity_ts_aic$AICc[4], 2)
fecundity_ts_aic$Delta_AICc[5] <- round(fecundity_ts_aic$AICc[1]-fecundity_ts_aic$AICc[5], 2)
fecundity_ts_aic$Delta_AICc[6] <- round(fecundity_ts_aic$AICc[1]-fecundity_ts_aic$AICc[6], 2)
fecundity_ts_aic$Delta_AICc[7] <- round(fecundity_ts_aic$AICc[1]-fecundity_ts_aic$AICc[7], 2)
fecundity_ts_aic  #Model 3 is best.

fecundity_ls_aic <- as.data.frame(1:7)
fecundity_ls_aic$AICc <- "NA"
fecundity_ls_aic$Model <- "NA"
fecundity_ls_aic$LL <- "NA"
fecundity_ls_aic$AICc[1] <- AICc(fecundity_gl_null)
fecundity_ls_aic$Model[1] <- "Null"
fecundity_ls_aic$LL[1] <- logLik(fecundity_gl_null)
fecundity_ls_aic$AICc[2] <- AICc(fecundity_g1l)
fecundity_ls_aic$Model[2] <- "m1"
fecundity_ls_aic$LL[2] <- logLik(fecundity_g1l)
fecundity_ls_aic$AICc[3] <- AICc(fecundity_g2.1l)
fecundity_ls_aic$Model[3] <- "m2"
fecundity_ls_aic$LL[3] <- logLik(fecundity_g2.1l)
fecundity_ls_aic$AICc[4] <- AICc(fecundity_g3.1l)
fecundity_ls_aic$Model[4] <- "m3"
fecundity_ls_aic$LL[4] <- logLik(fecundity_g3.1l)
fecundity_ls_aic$AICc[5] <- AICc(fecundity_g4.1l)
fecundity_ls_aic$Model[5] <- "m4"
fecundity_ls_aic$LL[5] <- logLik(fecundity_g4.1l)
fecundity_ls_aic$AICc[6] <- AICc(fecundity_g5.1l)
fecundity_ls_aic$Model[6] <- "m5"
fecundity_ls_aic$LL[6] <- logLik(fecundity_g5.1l)
fecundity_ls_aic$AICc[7] <- AICc(fecundity_g6.1l)
fecundity_ls_aic$Model[7] <- "m6"
fecundity_ls_aic$LL[7] <- logLik(fecundity_g6.1l)
fecundity_ls_aic$AICc <- as.numeric(fecundity_ls_aic$AICc)
fecundity_ls_aic$LL <- as.numeric(fecundity_ls_aic$LL)
fecundity_ls_aic <- fecundity_ls_aic[order(fecundity_ls_aic$AICc), ]
fecundity_ls_aic$Delta_AICc <- "0.00"
fecundity_ls_aic$Delta_AICc[2] <- round(fecundity_ls_aic$AICc[1]-fecundity_ls_aic$AICc[2], 2)
fecundity_ls_aic$Delta_AICc[3] <- round(fecundity_ls_aic$AICc[1]-fecundity_ls_aic$AICc[3], 2)
fecundity_ls_aic$Delta_AICc[4] <- round(fecundity_ls_aic$AICc[1]-fecundity_ls_aic$AICc[4], 2)
fecundity_ls_aic$Delta_AICc[5] <- round(fecundity_ls_aic$AICc[1]-fecundity_ls_aic$AICc[5], 2)
fecundity_ls_aic$Delta_AICc[6] <- round(fecundity_ls_aic$AICc[1]-fecundity_ls_aic$AICc[6], 2)
fecundity_ls_aic$Delta_AICc[7] <- round(fecundity_ls_aic$AICc[1]-fecundity_ls_aic$AICc[7], 2)
fecundity_ls_aic  #Model 6 is best

# GAM s(bs = 'tp') optimised
fecundity_tso_aic <- as.data.frame(1:7)
fecundity_tso_aic$AICc <- "NA"
fecundity_tso_aic$Model <- "NA"
fecundity_tso_aic$LL <- "NA"
fecundity_tso_aic$AICc[1] <- AICc(fecundity_gt_null)
fecundity_tso_aic$Model[1] <- "Null"
fecundity_tso_aic$LL[1] <- logLik(fecundity_gt_null)
fecundity_tso_aic$AICc[2] <- AICc(fecundity_g1t)
fecundity_tso_aic$Model[2] <- "m1"
fecundity_tso_aic$LL[2] <- logLik(fecundity_g1t)
fecundity_tso_aic$AICc[3] <- AICc(fecundity_g2.1ot)
fecundity_tso_aic$Model[3] <- "m2"
fecundity_tso_aic$LL[3] <- logLik(fecundity_g2.1ot)
fecundity_tso_aic$AICc[4] <- AICc(fecundity_g3.1ot)
fecundity_tso_aic$Model[4] <- "m3"
fecundity_tso_aic$LL[4] <- logLik(fecundity_g3.1ot)
fecundity_tso_aic$AICc[5] <- AICc(fecundity_g4.1ot)
fecundity_tso_aic$Model[5] <- "m4"
fecundity_tso_aic$LL[5] <- logLik(fecundity_g4.1ot)
fecundity_tso_aic$AICc[6] <- AICc(fecundity_g5.1ot)
fecundity_tso_aic$Model[6] <- "m5"
fecundity_tso_aic$LL[6] <- logLik(fecundity_g5.1ot)
fecundity_tso_aic$AICc[7] <- AICc(fecundity_g6.1ot)
fecundity_tso_aic$Model[7] <- "m6"
fecundity_tso_aic$LL[7] <- logLik(fecundity_g6.1ot)
fecundity_tso_aic$AICc <- as.numeric(fecundity_tso_aic$AICc)
fecundity_tso_aic$LL <- as.numeric(fecundity_tso_aic$LL)
fecundity_tso_aic <- fecundity_tso_aic[order(fecundity_tso_aic$AICc), ]
fecundity_tso_aic$Delta_AICc <- "0.00"
fecundity_tso_aic$Delta_AICc[2] <- round(fecundity_tso_aic$AICc[1]-fecundity_tso_aic$AICc[2], 2)
fecundity_tso_aic$Delta_AICc[3] <- round(fecundity_tso_aic$AICc[1]-fecundity_tso_aic$AICc[3], 2)
fecundity_tso_aic$Delta_AICc[4] <- round(fecundity_tso_aic$AICc[1]-fecundity_tso_aic$AICc[4], 2)
fecundity_tso_aic$Delta_AICc[5] <- round(fecundity_tso_aic$AICc[1]-fecundity_tso_aic$AICc[5], 2)
fecundity_tso_aic$Delta_AICc[6] <- round(fecundity_tso_aic$AICc[1]-fecundity_tso_aic$AICc[6], 2)
fecundity_tso_aic$Delta_AICc[7] <- round(fecundity_tso_aic$AICc[1]-fecundity_tso_aic$AICc[7], 2)
fecundity_tso_aic  #Model 3 is best.

fecundity_lso_aic <- as.data.frame(1:7)
fecundity_lso_aic$AICc <- "NA"
fecundity_lso_aic$Model <- "NA"
fecundity_lso_aic$LL <- "NA"
fecundity_lso_aic$AICc[1] <- AICc(fecundity_gl_null)
fecundity_lso_aic$Model[1] <- "Null"
fecundity_lso_aic$LL[1] <- logLik(fecundity_gl_null)
fecundity_lso_aic$AICc[2] <- AICc(fecundity_g1l)
fecundity_lso_aic$Model[2] <- "m1"
fecundity_lso_aic$LL[2] <- logLik(fecundity_g1l)
fecundity_lso_aic$AICc[3] <- AICc(fecundity_g2.1ol)
fecundity_lso_aic$Model[3] <- "m2"
fecundity_lso_aic$LL[3] <- logLik(fecundity_g2.1ol)
fecundity_lso_aic$AICc[4] <- AICc(fecundity_g3.1ol)
fecundity_lso_aic$Model[4] <- "m3"
fecundity_lso_aic$LL[4] <- logLik(fecundity_g3.1ol)
fecundity_lso_aic$AICc[5] <- AICc(fecundity_g4.1ol)
fecundity_lso_aic$Model[5] <- "m4"
fecundity_lso_aic$LL[5] <- logLik(fecundity_g4.1ol)
fecundity_lso_aic$AICc[6] <- AICc(fecundity_g5.1ol)
fecundity_lso_aic$Model[6] <- "m5"
fecundity_lso_aic$LL[6] <- logLik(fecundity_g5.1ol)
fecundity_lso_aic$AICc[7] <- AICc(fecundity_g6.1ol)
fecundity_lso_aic$Model[7] <- "m6"
fecundity_lso_aic$LL[7] <- logLik(fecundity_g6.1ol)
fecundity_lso_aic$AICc <- as.numeric(fecundity_lso_aic$AICc)
fecundity_lso_aic$LL <- as.numeric(fecundity_lso_aic$LL)
fecundity_lso_aic <- fecundity_lso_aic[order(fecundity_lso_aic$AICc), ]
fecundity_lso_aic$Delta_AICc <- "0.00"
fecundity_lso_aic$Delta_AICc[2] <- round(fecundity_lso_aic$AICc[1]-fecundity_lso_aic$AICc[2], 2)
fecundity_lso_aic$Delta_AICc[3] <- round(fecundity_lso_aic$AICc[1]-fecundity_lso_aic$AICc[3], 2)
fecundity_lso_aic$Delta_AICc[4] <- round(fecundity_lso_aic$AICc[1]-fecundity_lso_aic$AICc[4], 2)
fecundity_lso_aic$Delta_AICc[5] <- round(fecundity_lso_aic$AICc[1]-fecundity_lso_aic$AICc[5], 2)
fecundity_lso_aic$Delta_AICc[6] <- round(fecundity_lso_aic$AICc[1]-fecundity_lso_aic$AICc[6], 2)
fecundity_lso_aic$Delta_AICc[7] <- round(fecundity_lso_aic$AICc[1]-fecundity_lso_aic$AICc[7], 2)
fecundity_lso_aic  #Model 3 is best


# GAM ti(bs='cr')
fecundity_tsc_aic <- as.data.frame(1:7)
fecundity_tsc_aic$AICc <- "NA"
fecundity_tsc_aic$Model <- "NA"
fecundity_tsc_aic$LL <- "NA"
fecundity_tsc_aic$AICc[1] <- AICc(fecundity_gt_null)
fecundity_tsc_aic$Model[1] <- "Null"
fecundity_tsc_aic$LL[1] <- logLik(fecundity_gt_null)
fecundity_tsc_aic$AICc[2] <- AICc(fecundity_g1.1t)
fecundity_tsc_aic$Model[2] <- "m1"
fecundity_tsc_aic$LL[2] <- logLik(fecundity_g1.1t)
fecundity_tsc_aic$AICc[3] <- AICc(fecundity_g2.2t)
fecundity_tsc_aic$Model[3] <- "m2"
fecundity_tsc_aic$LL[3] <- logLik(fecundity_g2.2t)
fecundity_tsc_aic$AICc[4] <- AICc(fecundity_g3.2t)
fecundity_tsc_aic$Model[4] <- "m3"
fecundity_tsc_aic$LL[4] <- logLik(fecundity_g3.2t)
fecundity_tsc_aic$AICc[5] <- AICc(fecundity_g4.2t)
fecundity_tsc_aic$Model[5] <- "m4"
fecundity_tsc_aic$LL[5] <- logLik(fecundity_g4.2t)
fecundity_tsc_aic$AICc[6] <- AICc(fecundity_g5.2t)
fecundity_tsc_aic$Model[6] <- "m5"
fecundity_tsc_aic$LL[6] <- logLik(fecundity_g5.2t)
fecundity_tsc_aic$AICc[7] <- AICc(fecundity_g6.2t)
fecundity_tsc_aic$Model[7] <- "m6"
fecundity_tsc_aic$LL[7] <- logLik(fecundity_g6.2t)
fecundity_tsc_aic$AICc <- as.numeric(fecundity_tsc_aic$AICc)
fecundity_tsc_aic$LL <- as.numeric(fecundity_tsc_aic$LL)
fecundity_tsc_aic <- fecundity_tsc_aic[order(fecundity_tsc_aic$AICc), ]
fecundity_tsc_aic$Delta_AICc <- "0.00"
fecundity_tsc_aic$Delta_AICc[2] <- round(fecundity_tsc_aic$AICc[1]-fecundity_tsc_aic$AICc[2], 2)
fecundity_tsc_aic$Delta_AICc[3] <- round(fecundity_tsc_aic$AICc[1]-fecundity_tsc_aic$AICc[3], 2)
fecundity_tsc_aic$Delta_AICc[4] <- round(fecundity_tsc_aic$AICc[1]-fecundity_tsc_aic$AICc[4], 2)
fecundity_tsc_aic$Delta_AICc[5] <- round(fecundity_tsc_aic$AICc[1]-fecundity_tsc_aic$AICc[5], 2)
fecundity_tsc_aic$Delta_AICc[6] <- round(fecundity_tsc_aic$AICc[1]-fecundity_tsc_aic$AICc[6], 2)
fecundity_tsc_aic$Delta_AICc[7] <- round(fecundity_tsc_aic$AICc[1]-fecundity_tsc_aic$AICc[7], 2)
fecundity_tsc_aic  #Model 6 is best

fecundity_lsc_aic <- as.data.frame(1:7)
fecundity_lsc_aic$AICc <- "NA"
fecundity_lsc_aic$Model <- "NA"
fecundity_lsc_aic$LL <- "NA"
fecundity_lsc_aic$AICc[1] <- AICc(fecundity_gl_null)
fecundity_lsc_aic$Model[1] <- "Null"
fecundity_lsc_aic$LL[1] <- logLik(fecundity_gl_null)
fecundity_lsc_aic$AICc[2] <- AICc(fecundity_g1.1l)
fecundity_lsc_aic$Model[2] <- "m1"
fecundity_lsc_aic$LL[2] <- logLik(fecundity_g1.1l)
fecundity_lsc_aic$AICc[3] <- AICc(fecundity_g2.2l)
fecundity_lsc_aic$Model[3] <- "m2"
fecundity_lsc_aic$LL[3] <- logLik(fecundity_g2.2l)
fecundity_lsc_aic$AICc[4] <- AICc(fecundity_g3.2l)
fecundity_lsc_aic$Model[4] <- "m3"
fecundity_lsc_aic$LL[4] <- logLik(fecundity_g3.2l)
fecundity_lsc_aic$AICc[5] <- AICc(fecundity_g4.2l)
fecundity_lsc_aic$Model[5] <- "m4"
fecundity_lsc_aic$LL[5] <- logLik(fecundity_g4.2l)
fecundity_lsc_aic$AICc[6] <- AICc(fecundity_g5.2l)
fecundity_lsc_aic$Model[6] <- "m5"
fecundity_lsc_aic$LL[6] <- logLik(fecundity_g5.2l)
fecundity_lsc_aic$AICc[7] <- AICc(fecundity_g6.2l)
fecundity_lsc_aic$Model[7] <- "m6"
fecundity_lsc_aic$LL[7] <- logLik(fecundity_g6.2l)
fecundity_lsc_aic$AICc <- as.numeric(fecundity_lsc_aic$AICc)
fecundity_lsc_aic$LL <- as.numeric(fecundity_lsc_aic$LL)
fecundity_lsc_aic <- fecundity_lsc_aic[order(fecundity_lsc_aic$AICc), ]
fecundity_lsc_aic$Delta_AICc <- "0.00"
fecundity_lsc_aic$Delta_AICc[2] <- round(fecundity_lsc_aic$AICc[1]-fecundity_lsc_aic$AICc[2], 2)
fecundity_lsc_aic$Delta_AICc[3] <- round(fecundity_lsc_aic$AICc[1]-fecundity_lsc_aic$AICc[3], 2)
fecundity_lsc_aic$Delta_AICc[4] <- round(fecundity_lsc_aic$AICc[1]-fecundity_lsc_aic$AICc[4], 2)
fecundity_lsc_aic$Delta_AICc[5] <- round(fecundity_lsc_aic$AICc[1]-fecundity_lsc_aic$AICc[5], 2)
fecundity_lsc_aic$Delta_AICc[6] <- round(fecundity_lsc_aic$AICc[1]-fecundity_lsc_aic$AICc[6], 2)
fecundity_lsc_aic$Delta_AICc[7] <- round(fecundity_lsc_aic$AICc[1]-fecundity_lsc_aic$AICc[7], 2)
fecundity_lsc_aic  # Model 6 is best


# GAM (bs = 'cr') optimised
fecundity_tsco_aic <- as.data.frame(1:7)
fecundity_tsco_aic$AICc <- "NA"
fecundity_tsco_aic$Model <- "NA"
fecundity_tsco_aic$LL <- "NA"
fecundity_tsco_aic$AICc[1] <- AICc(fecundity_gt_null)
fecundity_tsco_aic$Model[1] <- "Null"
fecundity_tsco_aic$LL[1] <- logLik(fecundity_gt_null)
fecundity_tsco_aic$AICc[2] <- AICc(fecundity_g1.1ot)
fecundity_tsco_aic$Model[2] <- "m1"
fecundity_tsco_aic$LL[2] <- logLik(fecundity_g1.1ot)
fecundity_tsco_aic$AICc[3] <- AICc(fecundity_g2.2ot)
fecundity_tsco_aic$Model[3] <- "m2"
fecundity_tsco_aic$LL[3] <- logLik(fecundity_g2.2ot)
fecundity_tsco_aic$AICc[4] <- AICc(fecundity_g3.2ot)
fecundity_tsco_aic$Model[4] <- "m3"
fecundity_tsco_aic$LL[4] <- logLik(fecundity_g3.2ot)
fecundity_tsco_aic$AICc[5] <- AICc(fecundity_g4.2ot)
fecundity_tsco_aic$Model[5] <- "m4"
fecundity_tsco_aic$LL[5] <- logLik(fecundity_g4.2ot)
fecundity_tsco_aic$AICc[6] <- AICc(fecundity_g5.2ot)
fecundity_tsco_aic$Model[6] <- "m5"
fecundity_tsco_aic$LL[6] <- logLik(fecundity_g5.2ot)
fecundity_tsco_aic$AICc[7] <- AICc(fecundity_g6.2t)
fecundity_tsco_aic$Model[7] <- "m6"
fecundity_tsco_aic$LL[7] <- logLik(fecundity_g6.2t)
fecundity_tsco_aic$AICc <- as.numeric(fecundity_tsco_aic$AICc)
fecundity_tsco_aic$LL <- as.numeric(fecundity_tsco_aic$LL)
fecundity_tsco_aic <- fecundity_tsco_aic[order(fecundity_tsco_aic$AICc), ]
fecundity_tsco_aic$Delta_AICc <- "0.00"
fecundity_tsco_aic$Delta_AICc[2] <- round(fecundity_tsco_aic$AICc[1]-fecundity_tsco_aic$AICc[2], 2)
fecundity_tsco_aic$Delta_AICc[3] <- round(fecundity_tsco_aic$AICc[1]-fecundity_tsco_aic$AICc[3], 2)
fecundity_tsco_aic$Delta_AICc[4] <- round(fecundity_tsco_aic$AICc[1]-fecundity_tsco_aic$AICc[4], 2)
fecundity_tsco_aic$Delta_AICc[5] <- round(fecundity_tsco_aic$AICc[1]-fecundity_tsco_aic$AICc[5], 2)
fecundity_tsco_aic$Delta_AICc[6] <- round(fecundity_tsco_aic$AICc[1]-fecundity_tsco_aic$AICc[6], 2)
fecundity_tsco_aic$Delta_AICc[7] <- round(fecundity_tsco_aic$AICc[1]-fecundity_tsco_aic$AICc[7], 2)
fecundity_tsco_aic  #Model 6 is best.

fecundity_lsco_aic <- as.data.frame(1:7)
fecundity_lsco_aic$AICc <- "NA"
fecundity_lsco_aic$Model <- "NA"
fecundity_lsco_aic$LL <- "NA"
fecundity_lsco_aic$AICc[1] <- AICc(fecundity_gl_null)
fecundity_lsco_aic$Model[1] <- "Null"
fecundity_lsco_aic$LL[1] <- logLik(fecundity_gl_null)
fecundity_lsco_aic$AICc[2] <- AICc(fecundity_g1.10l)
fecundity_lsco_aic$Model[2] <- "m1"
fecundity_lsco_aic$LL[2] <- logLik(fecundity_g1.1ol)
fecundity_lsco_aic$AICc[3] <- AICc(fecundity_g2.2l)
fecundity_lsco_aic$Model[3] <- "m2"
fecundity_lsco_aic$LL[3] <- logLik(fecundity_g2.2l)
fecundity_lsco_aic$AICc[4] <- AICc(fecundity_g3.2l)
fecundity_lsco_aic$Model[4] <- "m3"
fecundity_lsco_aic$LL[4] <- logLik(fecundity_g3.2l)
fecundity_lsco_aic$AICc[5] <- AICc(fecundity_g4.2l)
fecundity_lsco_aic$Model[5] <- "m4"
fecundity_lsco_aic$LL[5] <- logLik(fecundity_g4.2l)
fecundity_lsco_aic$AICc[6] <- AICc(fecundity_g5.2l)
fecundity_lsco_aic$Model[6] <- "m5"
fecundity_lsco_aic$LL[6] <- logLik(fecundity_g5.2l)
fecundity_lsco_aic$AICc[7] <- AICc(fecundity_g6.2l)
fecundity_lsco_aic$Model[7] <- "m6"
fecundity_lsco_aic$LL[7] <- logLik(fecundity_g6.2l)
fecundity_lsco_aic$AICc <- as.numeric(fecundity_lsco_aic$AICc)
fecundity_lsco_aic$LL <- as.numeric(fecundity_lsco_aic$LL)
fecundity_lsco_aic <- fecundity_lsco_aic[order(fecundity_lsco_aic$AICc), ]
fecundity_lsco_aic$Delta_AICc <- "0.00"
fecundity_lsco_aic$Delta_AICc[2] <- round(fecundity_lsco_aic$AICc[1]-fecundity_lsco_aic$AICc[2], 2)
fecundity_lsco_aic$Delta_AICc[3] <- round(fecundity_lsco_aic$AICc[1]-fecundity_lsco_aic$AICc[3], 2)
fecundity_lsco_aic$Delta_AICc[4] <- round(fecundity_lsco_aic$AICc[1]-fecundity_lsco_aic$AICc[4], 2)
fecundity_lsco_aic$Delta_AICc[5] <- round(fecundity_lsco_aic$AICc[1]-fecundity_lsco_aic$AICc[5], 2)
fecundity_lsco_aic$Delta_AICc[6] <- round(fecundity_lsco_aic$AICc[1]-fecundity_lsco_aic$AICc[6], 2)
fecundity_lsco_aic$Delta_AICc[7] <- round(fecundity_lsco_aic$AICc[1]-fecundity_lsco_aic$AICc[7], 2)
fecundity_lsco_aic  #Model 6 is best



# Predict for fecundity in response to TSF ----
# For GLMER m3 
# For littoralis GAM default, tp, cr, and cr optimised m6
# For littoralis GAM tp optimised m3
# For torulosa GAM default and tp m3
# For torulosa GAM cr m6

fec_m3tl <-  data.frame(TSF = min(tor_fecundity$TSF),
                       FPC = seq(min(tor_fecundity$FPC, na.rm = T), max(tor_fecundity$FPC, na.rm = T), length = 50))
fec_m3ta <-  data.frame(TSF = mean(tor_fecundity$TSF),
                       FPC = seq(min(tor_fecundity$FPC, na.rm = T), max(tor_fecundity$FPC, na.rm = T), length = 50))
fec_m3th <-  data.frame(TSF = max(tor_fecundity$TSF),
                       FPC = seq(min(tor_fecundity$FPC, na.rm = T), max(tor_fecundity$FPC, na.rm = T), length = 50))

fec_m3tl_lm <-  data.frame(r_TSF = min(tor_fecundity$r_TSF),
                           r_FPC = seq(min(tor_fecundity$r_FPC, na.rm = T), max(tor_fecundity$r_FPC, na.rm = T), length = 50))
fec_m3ta_lm <-  data.frame(r_TSF = mean(tor_fecundity$r_TSF),
                           r_FPC = seq(min(tor_fecundity$r_FPC, na.rm = T), max(tor_fecundity$r_FPC, na.rm = T), length = 50))
fec_m3th_lm <-  data.frame(r_TSF = max(tor_fecundity$r_TSF),
                           r_FPC = seq(min(tor_fecundity$r_FPC, na.rm = T), max(tor_fecundity$r_FPC, na.rm = T), length = 50))



fec_m6tl <-  data.frame(TSF = min(tor_fecundity$TSF),
                        Height_cm = seq(min(tor_fecundity$Height_cm, na.rm = T), max(tor_fecundity$Height_cm, na.rm = T), length = 50))
fec_m6ta <-  data.frame(TSF = mean(tor_fecundity$TSF),
                        Height_cm = seq(min(tor_fecundity$Height_cm, na.rm = T), max(tor_fecundity$Height_cm, na.rm = T), length = 50))
fec_m6th <-  data.frame(TSF = max(tor_fecundity$TSF),
                        Height_cm = seq(min(tor_fecundity$Height_cm, na.rm = T), max(tor_fecundity$Height_cm, na.rm = T), length = 50))




fec_m3ll_lm <-  data.frame(r_TSF = min(lit_fecundity$r_TSF),
                           r_FPC = seq(min(lit_fecundity$r_FPC, na.rm = T), max(lit_fecundity$r_FPC, na.rm = T), length = 50))
fec_m3la_lm <-  data.frame(r_TSF = mean(lit_fecundity$r_TSF),
                           r_FPC = seq(min(lit_fecundity$r_FPC, na.rm = T), max(lit_fecundity$r_FPC, na.rm = T), length = 50))
fec_m3lh_lm <-  data.frame(r_TSF = max(lit_fecundity$r_TSF),
                           r_FPC = seq(min(lit_fecundity$r_FPC, na.rm = T), max(lit_fecundity$r_FPC, na.rm = T), length = 50))





fec_m3ll <-  data.frame(TSF = min(lit_fecundity$TSF),
                           FPC = seq(min(lit_fecundity$FPC, na.rm = T), max(lit_fecundity$FPC, na.rm = T), length = 50))
fec_m3la <-  data.frame(TSF = mean(lit_fecundity$TSF),
                           FPC = seq(min(lit_fecundity$FPC, na.rm = T), max(lit_fecundity$FPC, na.rm = T), length = 50))
fec_m3lh <-  data.frame(TSF = max(lit_fecundity$TSF),
                           FPC = seq(min(lit_fecundity$FPC, na.rm = T), max(lit_fecundity$FPC, na.rm = T), length = 50))


fec_m6ll <-  data.frame(TSF = min(lit_fecundity$TSF),
                        Height_cm = seq(min(lit_fecundity$Height_cm, na.rm = T), max(lit_fecundity$Height_cm, na.rm = T), length = 50))
fec_m6la <-  data.frame(TSF = mean(lit_fecundity$TSF),
                        Height_cm = seq(min(lit_fecundity$Height_cm, na.rm = T), max(lit_fecundity$Height_cm, na.rm = T), length = 50))
fec_m6lh <-  data.frame(TSF = max(lit_fecundity$TSF),
                        Height_cm = seq(min(lit_fecundity$Height_cm, na.rm = T), max(lit_fecundity$Height_cm, na.rm = T), length = 50))

# GLMER
pfec_m3tl_lm <- predictSE(fecundity_m3t, newdata = fec_m3tl_lm, se.fit = T, type = 'link')
fec_m3tl_lm$fit.link <- pfec_m3tl_lm$fit
fec_m3tl_lm$se.link <- pfec_m3tl_lm$se.fit
fec_m3tl_lm$lci.link <- fec_m3tl_lm$fit.link - (fec_m3tl_lm$se.link * 1.96)
fec_m3tl_lm$uci.link <- fec_m3tl_lm$fit.link + (fec_m3tl_lm$se.link * 1.96)
fec_m3tl_lm$fit <- exp(fec_m3tl_lm$fit.link)
fec_m3tl_lm$se <- exp(fec_m3tl_lm$se.link)
fec_m3tl_lm$lci <- exp(fec_m3tl_lm$lci.link)
fec_m3tl_lm$uci <- exp(fec_m3tl_lm$uci.link)
fec_m3tl_lm 


pfec_m3ta_lm <- predictSE(fecundity_m3t, newdata = fec_m3ta_lm, se.fit = T, type = 'link')
fec_m3ta_lm$fit.link <- pfec_m3ta_lm$fit
fec_m3ta_lm$se.link <- pfec_m3ta_lm$se.fit
fec_m3ta_lm$lci.link <- fec_m3ta_lm$fit.link - (fec_m3ta_lm$se.link * 1.96)
fec_m3ta_lm$uci.link <- fec_m3ta_lm$fit.link + (fec_m3ta_lm$se.link * 1.96)
fec_m3ta_lm$fit <- exp(fec_m3ta_lm$fit.link)
fec_m3ta_lm$se <- exp(fec_m3ta_lm$se.link)
fec_m3ta_lm$lci <- exp(fec_m3ta_lm$lci.link)
fec_m3ta_lm$uci <- exp(fec_m3ta_lm$uci.link)



pfec_m3th_lm <- predictSE(fecundity_m3t, newdata = fec_m3th_lm, se.fit = T, type = 'link')
fec_m3th_lm$fit.link <- pfec_m3th_lm$fit
fec_m3th_lm$se.link <- pfec_m3th_lm$se.fit
fec_m3th_lm$lci.link <- fec_m3th_lm$fit.link - (fec_m3th_lm$se.link * 1.96)
fec_m3th_lm$uci.link <- fec_m3th_lm$fit.link + (fec_m3th_lm$se.link * 1.96)
fec_m3th_lm$fit <- exp(fec_m3th_lm$fit.link)
fec_m3th_lm$se <- exp(fec_m3th_lm$se.link)
fec_m3th_lm$lci <- exp(fec_m3th_lm$lci.link)
fec_m3th_lm$uci <- exp(fec_m3th_lm$uci.link)



pfec_m3ll_lm <- predictSE(fecundity_m3l, newdata = fec_m3ll_lm, se.fit = T, type = 'link')
fec_m3ll_lm$fit.link <- pfec_m3ll_lm$fit
fec_m3ll_lm$se.link <- pfec_m3ll_lm$se.fit
fec_m3ll_lm$lci.link <- fec_m3ll_lm$fit.link - (fec_m3ll_lm$se.link * 1.96)
fec_m3ll_lm$uci.link <- fec_m3ll_lm$fit.link + (fec_m3ll_lm$se.link * 1.96)
fec_m3ll_lm$fit <- exp(fec_m3ll_lm$fit.link)
fec_m3ll_lm$se <- exp(fec_m3ll_lm$se.link)
fec_m3ll_lm$lci <- exp(fec_m3ll_lm$lci.link)
fec_m3ll_lm$uci <- exp(fec_m3ll_lm$uci.link)

pfec_m3la_lm <- predictSE(fecundity_m3l, newdata = fec_m3la_lm, se.fit = T, type = 'link')
fec_m3la_lm$fit.link <- pfec_m3la_lm$fit
fec_m3la_lm$se.link <- pfec_m3la_lm$se.fit
fec_m3la_lm$lci.link <- fec_m3la_lm$fit.link - (fec_m3la_lm$se.link * 1.96)
fec_m3la_lm$uci.link <- fec_m3la_lm$fit.link + (fec_m3la_lm$se.link * 1.96)
fec_m3la_lm$fit <- exp(fec_m3la_lm$fit.link)
fec_m3la_lm$se <- exp(fec_m3la_lm$se.link)
fec_m3la_lm$lci <- exp(fec_m3la_lm$lci.link)
fec_m3la_lm$uci <- exp(fec_m3la_lm$uci.link)

pfec_m3lh_lm <- predictSE(fecundity_m3l, newdata = fec_m3lh_lm, se.fit = T, type = 'link')
fec_m3lh_lm$fit.link <- pfec_m3lh_lm$fit
fec_m3lh_lm$se.link <- pfec_m3lh_lm$se.fit
fec_m3lh_lm$lci.link <- fec_m3lh_lm$fit.link - (fec_m3lh_lm$se.link * 1.96)
fec_m3lh_lm$uci.link <- fec_m3lh_lm$fit.link + (fec_m3lh_lm$se.link * 1.96)
fec_m3lh_lm$fit <- exp(fec_m3lh_lm$fit.link)
fec_m3lh_lm$se <- exp(fec_m3lh_lm$se.link)
fec_m3lh_lm$lci <- exp(fec_m3lh_lm$lci.link)
fec_m3lh_lm$uci <- exp(fec_m3lh_lm$uci.link)



# GAM default
nfec_g3tl <- fec_m3tl
nfec_g3ta <- fec_m3ta
nfec_g3th <- fec_m3th

nfec_g6ll <- fec_m6ll
nfec_g6la <- fec_m6la
nfec_g6lh <- fec_m6lh



pfec_g3tl <- predict(fecundity_g3t, newdata = nfec_g3tl, se.fit = T, type = 'link')
nfec_g3tl$fit <- pfec_g3tl$fit
nfec_g3tl$se <- pfec_g3tl$se.fit
nfec_g3tl$lci <- nfec_g3tl$fit - (nfec_g3tl$se * 1.96)
nfec_g3tl$uci <- nfec_g3tl$fit + (nfec_g3tl$se * 1.96)
nfec_g3tl # Not a good model, we get negative values. Think this needs a different modelling approach.


pfec_g3ta <- predict(fecundity_g3t, newdata = nfec_g3tl, se.fit = T, type = 'response')
nfec_g3ta$fit <- pfec_g3ta$fit
nfec_g3ta$se <- pfec_g3ta$se.fit
nfec_g3ta$lci <- nfec_g3ta$fit - (nfec_g3ta$se * 1.96)
nfec_g3ta$uci <- nfec_g3ta$fit + (nfec_g3ta$se * 1.96)

pfec_g3th <- predict(fecundity_g3t, newdata = nfec_g3tl, se.fit = T, type = 'response')
nfec_g3th$fit <- pfec_g3th$fit
nfec_g3th$se <- pfec_g3th$se.fit
nfec_g3th$lci <- nfec_g3th$fit - (nfec_g3th$se * 1.96)
nfec_g3th$uci <- nfec_g3th$fit + (nfec_g3th$se * 1.96)



pfec_g6ll <- predict(fecundity_g6l, newdata = nfec_g6ll, se.fit = T, type = 'response')
nfec_g6ll$fit <- pfec_g6ll$fit
nfec_g6ll$se <- pfec_g6ll$se.fit
nfec_g6ll$lci <- nfec_g6ll$fit - (nfec_g6ll$se * 1.96)
nfec_g6ll$uci <- nfec_g6ll$fit + (nfec_g6ll$se * 1.96)
nfec_g6ll # Somewhat reasonable but does not perform well for plants below 1m in height.


pfec_g6la <- predict(fecundity_g6l, newdata = nfec_g6la, se.fit = T, type = 'response')
nfec_g6la$fit <- pfec_g6la$fit
nfec_g6la$se <- pfec_g6la$se.fit
nfec_g6la$lci <- nfec_g6la$fit - (nfec_g6la$se * 1.96)
nfec_g6la$uci <- nfec_g6la$fit + (nfec_g6la$se * 1.96)

pfec_g6lh <- predict(fecundity_g6l, newdata = nfec_g6lh, se.fit = T, type = 'response')
nfec_g6lh$fit <- pfec_g6lh$fit
nfec_g6lh$se <- pfec_g6lh$se.fit
nfec_g6lh$lci <- nfec_g6lh$fit - (nfec_g6lh$se * 1.96)
nfec_g6lh$uci <- nfec_g6lh$fit + (nfec_g6lh$se * 1.96)


# GAM tp
nfec_g3.1tl <- fec_m3tl
nfec_g3.1ta <- fec_m3ta
nfec_g3.1th <- fec_m3th

nfec_g6.1ll <- fec_m6ll
nfec_g6.1la <- fec_m6la
nfec_g6.1lh <- fec_m6lh

pfec_g3.1tl <- predict(fecundity_g3.1t, newdata = nfec_g3.1tl, se.fit = T, type = 'link')
nfec_g3.1tl$fit <- pfec_g3.1tl$fit
nfec_g3.1tl$se <- pfec_g3.1tl$se.fit
nfec_g3.1tl$lci <- nfec_g3.1tl$fit - (nfec_g3.1tl$se * 1.96)
nfec_g3.1tl$uci <- nfec_g3.1tl$fit + (nfec_g3.1tl$se * 1.96)

pfec_g3.1ta <- predict(fecundity_g3.1t, newdata = nfec_g3.1tl, se.fit = T, type = 'response')
nfec_g3.1ta$fit <- pfec_g3.1ta$fit
nfec_g3.1ta$se <- pfec_g3.1ta$se.fit
nfec_g3.1ta$lci <- nfec_g3.1ta$fit - (nfec_g3.1ta$se * 1.96)
nfec_g3.1ta$uci <- nfec_g3.1ta$fit + (nfec_g3.1ta$se * 1.96)

pfec_g3.1th <- predict(fecundity_g3.1t, newdata = nfec_g3.1tl, se.fit = T, type = 'response')
nfec_g3.1th$fit <- pfec_g3.1th$fit
nfec_g3.1th$se <- pfec_g3.1th$se.fit
nfec_g3.1th$lci <- nfec_g3.1th$fit - (nfec_g3.1th$se * 1.96)
nfec_g3.1th$uci <- nfec_g3.1th$fit + (nfec_g3.1th$se * 1.96)



pfec_g6.1ll <- predict(fecundity_g6.1l, newdata = nfec_g6.1ll, se.fit = T, type = 'response')
nfec_g6.1ll$fit <- pfec_g6.1ll$fit
nfec_g6.1ll$se <- pfec_g6.1ll$se.fit
nfec_g6.1ll$lci <- nfec_g6.1ll$fit - (nfec_g6.1ll$se * 1.96)
nfec_g6.1ll$uci <- nfec_g6.1ll$fit + (nfec_g6.1ll$se * 1.96)

pfec_g6.1la <- predict(fecundity_g6.1l, newdata = nfec_g6.1la, se.fit = T, type = 'response')
nfec_g6.1la$fit <- pfec_g6.1la$fit
nfec_g6.1la$se <- pfec_g6.1la$se.fit
nfec_g6.1la$lci <- nfec_g6.1la$fit - (nfec_g6.1la$se * 1.96)
nfec_g6.1la$uci <- nfec_g6.1la$fit + (nfec_g6.1la$se * 1.96)

pfec_g6.1lh <- predict(fecundity_g6.1l, newdata = nfec_g6.1lh, se.fit = T, type = 'response')
nfec_g6.1lh$fit <- pfec_g6.1lh$fit
nfec_g6.1lh$se <- pfec_g6.1lh$se.fit
nfec_g6.1lh$lci <- nfec_g6.1lh$fit - (nfec_g6.1lh$se * 1.96)
nfec_g6.1lh$uci <- nfec_g6.1lh$fit + (nfec_g6.1lh$se * 1.96)

# GAM tp opt
nfec_g3.1otl <- fec_m3tl
nfec_g3.1ota <- fec_m3ta
nfec_g3.1oth <- fec_m3th

nfec_g3.1oll <- fec_m3ll
nfec_g3.1ola <- fec_m3la
nfec_g3.1olh <- fec_m3lh

pfec_g3.1otl <- predict(fecundity_g3.1ot, newdata = nfec_g3.1otl, se.fit = T, type = 'link')
nfec_g3.1otl$fit <- pfec_g3.1otl$fit
nfec_g3.1otl$se <- pfec_g3.1otl$se.fit
nfec_g3.1otl$lci <- nfec_g3.1otl$fit - (nfec_g3.1otl$se * 1.96)
nfec_g3.1otl$uci <- nfec_g3.1otl$fit + (nfec_g3.1otl$se * 1.96)

pfec_g3.1ota <- predict(fecundity_g3.1ot, newdata = nfec_g3.1otl, se.fit = T, type = 'response')
nfec_g3.1ota$fit <- pfec_g3.1ota$fit
nfec_g3.1ota$se <- pfec_g3.1ota$se.fit
nfec_g3.1ota$lci <- nfec_g3.1ota$fit - (nfec_g3.1ota$se * 1.96)
nfec_g3.1ota$uci <- nfec_g3.1ota$fit + (nfec_g3.1ota$se * 1.96)

pfec_g3.1oth <- predict(fecundity_g3.1ot, newdata = nfec_g3.1otl, se.fit = T, type = 'response')
nfec_g3.1oth$fit <- pfec_g3.1oth$fit
nfec_g3.1oth$se <- pfec_g3.1oth$se.fit
nfec_g3.1oth$lci <- nfec_g3.1oth$fit - (nfec_g3.1oth$se * 1.96)
nfec_g3.1oth$uci <- nfec_g3.1oth$fit + (nfec_g3.1oth$se * 1.96)



pfec_g3.1oll <- predict(fecundity_g3.1ol, newdata = nfec_g3.1oll, se.fit = T, type = 'response')
nfec_g3.1oll$fit <- pfec_g3.1oll$fit
nfec_g3.1oll$se <- pfec_g3.1oll$se.fit
nfec_g3.1oll$lci <- nfec_g3.1oll$fit - (nfec_g3.1oll$se * 1.96)
nfec_g3.1oll$uci <- nfec_g3.1oll$fit + (nfec_g3.1oll$se * 1.96)

pfec_g3.1ola <- predict(fecundity_g3.1ol, newdata = nfec_g3.1ola, se.fit = T, type = 'response')
nfec_g3.1ola$fit <- pfec_g3.1ola$fit
nfec_g3.1ola$se <- pfec_g3.1ola$se.fit
nfec_g3.1ola$lci <- nfec_g3.1ola$fit - (nfec_g3.1ola$se * 1.96)
nfec_g3.1ola$uci <- nfec_g3.1ola$fit + (nfec_g3.1ola$se * 1.96)

pfec_g3.1olh <- predict(fecundity_g3.1ol, newdata = nfec_g3.1olh, se.fit = T, type = 'response')
nfec_g3.1olh$fit <- pfec_g3.1olh$fit
nfec_g3.1olh$se <- pfec_g3.1olh$se.fit
nfec_g3.1olh$lci <- nfec_g3.1olh$fit - (nfec_g3.1olh$se * 1.96)
nfec_g3.1olh$uci <- nfec_g3.1olh$fit + (nfec_g3.1olh$se * 1.96)

# GAM cr
nfec_g6.2tl <- fec_m6tl
nfec_g6.2ta <- fec_m6ta
nfec_g6.2th <- fec_m6th

nfec_g6.2ll <- fec_m6ll
nfec_g6.2la <- fec_m6la
nfec_g6.2lh <- fec_m6lh

pfec_g6.2tl <- predict(fecundity_g6.2t, newdata = nfec_g6.2tl, se.fit = T, type = 'link')
nfec_g6.2tl$fit <- pfec_g6.2tl$fit
nfec_g6.2tl$se <- pfec_g6.2tl$se.fit
nfec_g6.2tl$lci <- nfec_g6.2tl$fit - (nfec_g6.2tl$se * 1.96)
nfec_g6.2tl$uci <- nfec_g6.2tl$fit + (nfec_g6.2tl$se * 1.96)

pfec_g6.2ta <- predict(fecundity_g6.2t, newdata = nfec_g6.2tl, se.fit = T, type = 'response')
nfec_g6.2ta$fit <- pfec_g6.2ta$fit
nfec_g6.2ta$se <- pfec_g6.2ta$se.fit
nfec_g6.2ta$lci <- nfec_g6.2ta$fit - (nfec_g6.2ta$se * 1.96)
nfec_g6.2ta$uci <- nfec_g6.2ta$fit + (nfec_g6.2ta$se * 1.96)

pfec_g6.2th <- predict(fecundity_g6.2t, newdata = nfec_g6.2tl, se.fit = T, type = 'response')
nfec_g6.2th$fit <- pfec_g6.2th$fit
nfec_g6.2th$se <- pfec_g6.2th$se.fit
nfec_g6.2th$lci <- nfec_g6.2th$fit - (nfec_g6.2th$se * 1.96)
nfec_g6.2th$uci <- nfec_g6.2th$fit + (nfec_g6.2th$se * 1.96)



pfec_g6.2ll <- predict(fecundity_g6.2l, newdata = nfec_g6.2ll, se.fit = T, type = 'response')
nfec_g6.2ll$fit <- pfec_g6.2ll$fit
nfec_g6.2ll$se <- pfec_g6.2ll$se.fit
nfec_g6.2ll$lci <- nfec_g6.2ll$fit - (nfec_g6.2ll$se * 1.96)
nfec_g6.2ll$uci <- nfec_g6.2ll$fit + (nfec_g6.2ll$se * 1.96)

pfec_g6.2la <- predict(fecundity_g6.2l, newdata = nfec_g6.2la, se.fit = T, type = 'response')
nfec_g6.2la$fit <- pfec_g6.2la$fit
nfec_g6.2la$se <- pfec_g6.2la$se.fit
nfec_g6.2la$lci <- nfec_g6.2la$fit - (nfec_g6.2la$se * 1.96)
nfec_g6.2la$uci <- nfec_g6.2la$fit + (nfec_g6.2la$se * 1.96)

pfec_g6.2lh <- predict(fecundity_g6.2l, newdata = nfec_g6.2lh, se.fit = T, type = 'response')
nfec_g6.2lh$fit <- pfec_g6.2lh$fit
nfec_g6.2lh$se <- pfec_g6.2lh$se.fit
nfec_g6.2lh$lci <- nfec_g6.2lh$fit - (nfec_g6.2lh$se * 1.96)
nfec_g6.2lh$uci <- nfec_g6.2lh$fit + (nfec_g6.2lh$se * 1.96)

# GAM cr opt
nfec_g6.2tl <- fec_m6tl
nfec_g6.2ta <- fec_m6ta
nfec_g6.2th <- fec_m6th

nfec_g6.2oll <- fec_m6ll
nfec_g6.2ola <- fec_m6la
nfec_g6.2olh <- fec_m6lh

pfec_g6.2tl <- predict(fecundity_g6.2t, newdata = nfec_g6.2tl, se.fit = T, type = 'link')
nfec_g6.2tl$fit <- pfec_g6.2tl$fit
nfec_g6.2tl$se <- pfec_g6.2tl$se.fit
nfec_g6.2tl$lci <- nfec_g6.2tl$fit - (nfec_g6.2tl$se * 1.96)
nfec_g6.2tl$uci <- nfec_g6.2tl$fit + (nfec_g6.2tl$se * 1.96)

pfec_g6.2ta <- predict(fecundity_g6.2t, newdata = nfec_g6.2tl, se.fit = T, type = 'response')
nfec_g6.2ta$fit <- pfec_g6.2ta$fit
nfec_g6.2ta$se <- pfec_g6.2ta$se.fit
nfec_g6.2ta$lci <- nfec_g6.2ta$fit - (nfec_g6.2ta$se * 1.96)
nfec_g6.2ta$uci <- nfec_g6.2ta$fit + (nfec_g6.2ta$se * 1.96)

pfec_g6.2th <- predict(fecundity_g6.2t, newdata = nfec_g6.2tl, se.fit = T, type = 'response')
nfec_g6.2th$fit <- pfec_g6.2th$fit
nfec_g6.2th$se <- pfec_g6.2th$se.fit
nfec_g6.2th$lci <- nfec_g6.2th$fit - (nfec_g6.2th$se * 1.96)
nfec_g6.2th$uci <- nfec_g6.2th$fit + (nfec_g6.2th$se * 1.96)



pfec_g6.2oll <- predict(fecundity_g6.2ol, newdata = nfec_g6.2oll, se.fit = T, type = 'response')
nfec_g6.2oll$fit <- pfec_g6.2oll$fit
nfec_g6.2oll$se <- pfec_g6.2oll$se.fit
nfec_g6.2oll$lci <- nfec_g6.2oll$fit - (nfec_g6.2oll$se * 1.96)
nfec_g6.2oll$uci <- nfec_g6.2oll$fit + (nfec_g6.2oll$se * 1.96)

pfec_g6.2ola <- predict(fecundity_g6.2ol, newdata = nfec_g6.2ola, se.fit = T, type = 'response')
nfec_g6.2ola$fit <- pfec_g6.2ola$fit
nfec_g6.2ola$se <- pfec_g6.2ola$se.fit
nfec_g6.2ola$lci <- nfec_g6.2ola$fit - (nfec_g6.2ola$se * 1.96)
nfec_g6.2ola$uci <- nfec_g6.2ola$fit + (nfec_g6.2ola$se * 1.96)

pfec_g6.2olh <- predict(fecundity_g6.2ol, newdata = nfec_g6.2olh, se.fit = T, type = 'response')
nfec_g6.2olh$fit <- pfec_g6.2olh$fit
nfec_g6.2olh$se <- pfec_g6.2olh$se.fit
nfec_g6.2olh$lci <- nfec_g6.2olh$fit - (nfec_g6.2olh$se * 1.96)
nfec_g6.2olh$uci <- nfec_g6.2olh$fit + (nfec_g6.2olh$se * 1.96)



# Plot GLMERs for Fecundity ~ TSF ----
# We know looking at the data that GLMER is the only one producing reasonable estimates as the GAMs tend to predict high cone numbers at low heights or for FPC may give negative values
# NOTE we limit the y axis of the plot as cone number increases exponentially
dev.new(width = 12, height = 10, noRStudioGD = T, dpi = 300)
par(mfrow = c(1,2), mar = c(8,6,3,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))

plot(fec_m3ll_lm$r_FPC, fec_m3ll_lm$fit, type = 'l', col = 'blue', ylab = "", xlab = "", las = 1, ylim = c(0, 1000), xaxt = "n", cex.axis = 1.4)
axis(side = 1, at = seq(-2.4, 2.6, 0.5), labels = seq(36, 46, 1), cex.axis = 1.4)
mtext(side = 1, expression(bold("Foliage Projective Cover (%)")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Female fecundity (no. cones)")), line = 3.5, cex = 1.5)
mtext(paste("AICc = ", round(AICc(fecundity_m3l), 3), sep = ""), line = -1.5, cex = 1.2)
mtext(expression(bold("(a) Allocasuarina littoralis TSF")), cex = 2)
pg.ci(x = 'r_FPC', data = 'fec_m3ll_lm', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(fec_m3la_lm$r_FPC, fec_m3la_lm$fit, type = 'l', col = 'black')
pg.ci(x = 'r_FPC', data = 'fec_m3la_lm', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(fec_m3lh_lm$r_FPC, fec_m3lh_lm$fit, type = 'l', col = 'red')
pg.ci(x = 'r_FPC', data = 'fec_m3lh_lm', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(fec_m3tl_lm$r_FPC, fec_m3tl_lm$fit, type = 'l', col = 'blue', ylab = "", xlab = "", las = 1, ylim = c(0, 1000))
pg.ci(x = 'r_FPC', data = 'fec_m3tl_lm', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(fec_m3ta_lm$r_FPC, fec_m3ta_lm$fit, type = 'l', col = 'black')
pg.ci(x = 'r_FPC', data = 'fec_m3ta_lm', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(fec_m3th_lm$r_FPC, fec_m3th_lm$fit, type = 'l', col = 'red')
pg.ci(x = 'r_FPC', data = 'fec_m3th_lm', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')






# 6. QUESTION 3: How does contemporary fire history (i.e., fire frequency) and environmental attributes influence reproductive traits? ----
# Proportions of seedling, saplings, recruits and number of cones, seed size as response
# Fire frequency 
# Fire frequency * latitude
# Fire frequency * FPC



### For fecundity only
# TSF + height 
# TSF + height 
# TSF * latitude
# TSF * FPC
# TSF * precipitation seasonality
# TSF * temperature seasonality

# May include nested effects of location and transect

# Rescale fire frequency for GLMER
tor_transects$r_Fire_freq <- scale(tor_transects$Fire_freq)
lit_transects$r_Fire_freq <- scale(lit_transects$Fire_freq)

s# 6.1 Population structure ----
# 6.1.1 Proportion seedlings ---- 
seed_fft_null <- glmer(Proportion_seedlings ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
seed_fft_gnull <- gam(Proportion_seedlings ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')

seed_ffl_null <- glmer(Proportion_seedlings ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)
seed_ffl_gnull <- gam(Proportion_seedlings ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')



seed_ff_m1t <- glmer(Proportion_seedlings ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m1t)
seed_ff_g1t <- gam(Proportion_seedlings ~ s(Fire_freq, k = 8), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g1t)
plot(seed_ff_g1t)
seed_ff_g1.1t <- gam(Proportion_seedlings ~ s(Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g1.1t)
plot(seed_ff_g1.1t)

seed_ff_m1l <- glmer(Proportion_seedlings ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m1l)
seed_ff_g1l <- gam(Proportion_seedlings ~ s(Fire_freq, k = 3), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g1l)
plot(seed_ff_g1l)
seed_ff_g1.1l <- gam(Proportion_seedlings ~ s(Fire_freq, bs = 'cr', k = 3), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g1.1l)
plot(seed_ff_g1.1l)


seed_ff_m2t <- glmer(Proportion_seedlings ~ r_Fire_freq * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m2t)
seed_ff_g2t <- gam(Proportion_seedlings ~ s(Fire_freq, k = 8) + s(Latitude) + ti(Latitude, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g2t)
plot(seed_ff_g2t)
summary(seed_ff_g2t)
seed_ff_g2.1t <- gam(Proportion_seedlings ~ s(Fire_freq, k = 8) + s(Latitude) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g2.1t)
plot(seed_ff_g2.1t)
seed_ff_g2.2t <- gam(Proportion_seedlings ~ s(Fire_freq, bs = 'cr', k = 5) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g2.1t)
plot(seed_ff_g2.2t)


seed_ff_m2l <- glmer(Proportion_seedlings ~ Fire_freq * Latitude + (1|Location/Transect), family = binomial, data = lit_transects) # This model fails due to the low number of points available to model an interaction. We can only model an additive effect of Fire_freq and latitude
summary(seed_ff_m2l)
seed_ff_g2l <- gam(Proportion_seedlings ~ s(Fire_freq, k = 3) + s(Latitude, k = 8) + ti(Latitude, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g2l)
plot(seed_ff_g2l)
seed_ff_g2.1l <- gam(Proportion_seedlings ~ s(Fire_freq, k = 3) + s(Latitude, k = 8) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 8), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g2.1l)
plot(seed_ff_g2.1l)
seed_ff_g2.2l <- gam(Proportion_seedlings ~ s(Fire_freq, bs = 'cr', k = 3) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g2.1tl)
plot(seed_ff_g2.2l)



seed_ff_m3t <- glmer(Proportion_seedlings ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m3t)
seed_ff_g3t <- gam(Proportion_seedlings ~ s(Fire_freq, k = 8) + s(FPC) + ti(FPC, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g3t)
plot(seed_ff_g3t)
seed_ff_g3.1t <- gam(Proportion_seedlings ~ s(Fire_freq, k = 8) + s(FPC) + ti(FPC, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g3.1t)
plot(seed_ff_g3.1t)
seed_ff_g3.2t <- gam(Proportion_seedlings ~ s(Fire_freq, bs = 'cr', k = 5) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g3.1t)
plot(seed_ff_g3.2t)


seed_ff_m3l <- glmer(Proportion_seedlings ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m3l)
seed_ff_g3l <- gam(Proportion_seedlings ~ s(Fire_freq, k = 3) + s(FPC, k = 7) + ti(FPC, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g3l)
plot(seed_ff_g3l)
seed_ff_g3.1l <- gam(Proportion_seedlings ~ s(Fire_freq, k = 3) + s(FPC, k = 7) + ti(FPC, by = Fire_freq, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g3.1l)
plot(seed_ff_g3.1l)
seed_ff_g3.2l <- gam(Proportion_seedlings ~ s(Fire_freq, bs = 'cr', k = 3) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g3.2l)
plot(seed_ff_g3.2l)


seed_ff_m4t <- glmer(Proportion_seedlings ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m4t)
seed_ff_g4t <- gam(Proportion_seedlings ~ s(Fire_freq, k = 8) + s(Precip) + ti(Precip, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g4t)
plot(seed_ff_g4t)
seed_ff_g4.1t <- gam(Proportion_seedlings ~ s(Fire_freq, k = 8) + s(Precip) + ti(Precip, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g4.1t)
plot(seed_ff_g4.1t)
seed_ff_g4.2t <- gam(Proportion_seedlings ~ s(Fire_freq, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g4.1t)
plot(seed_ff_g4.2t)


seed_ff_m4l <- glmer(Proportion_seedlings ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m4l)
seed_ff_g4l <- gam(Proportion_seedlings ~ s(Fire_freq, k = 3) + s(Precip, k = 7) + ti(Precip, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g4l)
plot(seed_ff_g4l)
seed_ff_g4.1l <- gam(Proportion_seedlings ~ s(Fire_freq, k = 3) + s(Precip, k = 7) + ti(Precip, by = Fire_freq, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g4.1l)
plot(seed_ff_g4.1l)
seed_ff_g4.2l <- gam(Proportion_seedlings ~ s(Fire_freq, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g4.2l)
plot(seed_ff_g4.2l)



seed_ff_m5t <- glmer(Proportion_seedlings ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m5t)
seed_ff_g5t <- gam(Proportion_seedlings ~ s(Fire_freq, k = 8) + s(Temp) + ti(Temp, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g5t)
plot(seed_ff_g5t)
seed_ff_g5.1t <- gam(Proportion_seedlings ~ s(Fire_freq, k = 8) + s(Temp) + ti(Temp, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g5.1t)
plot(seed_ff_g5.1t)
seed_ff_g5.2t <- gam(Proportion_seedlings ~ s(Fire_freq, bs = 'cr', k = 5) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g5.1t)
plot(seed_ff_g5.2t)


seed_ff_m5l <- glmer(Proportion_seedlings ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m5l)
seed_ff_g5l <- gam(Proportion_seedlings ~ s(Fire_freq, k = 3) + s(Temp, k = 7) + ti(Temp, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g5l)
plot(seed_ff_g5l)
seed_ff_g5.1l <- gam(Proportion_seedlings ~ s(Fire_freq, k = 3) + s(Temp, k = 7) + ti(Temp, by = Fire_freq, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g5.1l)
plot(seed_ff_g5.1l)
seed_ff_g5.2l <- gam(Proportion_seedlings ~ s(Fire_freq, bs = 'cr', k = 3) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(seed_ff_g5.2l)
plot(seed_ff_g5.2l)


# Compare AICs for each modelling type
# GLMER
seed_ff_aic_t <- list(seed_fft_null, seed_ff_m1t, seed_ff_m2t, seed_ff_m3t, seed_ff_m4t, seed_ff_m5t)
aictab(seed_ff_aic_t) # Null model is best but model 2 is ranked within delta AICc <2.

seed_ff_aic_l <- list(seed_ffl_null, seed_ff_m1l, seed_ff_m2l, seed_ff_m3l, seed_ff-m4l, seed_ff_m5l)
aictab(seed_ff_aic_l) # Null model is best

# GAM defaults
seed_ff_td_aic <- as.data.frame(1:6)
seed_ff_td_aic$AICc <- "NA"
seed_ff_td_aic$Model <- "NA"
seed_ff_td_aic$LL <- "NA"
seed_ff_td_aic$AICc[1] <- AICc(seed_fft_gnull)
seed_ff_td_aic$Model[1] <- "Null"
seed_ff_td_aic$LL[1] <- logLik(seed_fft_gnull)
seed_ff_td_aic$AICc[2] <- AICc(seed_ff_g1t)
seed_ff_td_aic$Model[2] <- "m1"
seed_ff_td_aic$LL[2] <- logLik(seed_ff_g1t)
seed_ff_td_aic$AICc[3] <- AICc(seed_ff_g2t)
seed_ff_td_aic$Model[3] <- "m2"
seed_ff_td_aic$LL[3] <- logLik(seed_ff_g2t)
seed_ff_td_aic$AICc[4] <- AICc(seed_ff_g3t)
seed_ff_td_aic$Model[4] <- "m3"
seed_ff_td_aic$LL[4] <- logLik(seed_ff_g3t)
seed_ff_td_aic$AICc[5] <- AICc(seed_ff_g4t)
seed_ff_td_aic$Model[5] <- "m4"
seed_ff_td_aic$LL[5] <- logLik(seed_ff_g4t)
seed_ff_td_aic$AICc[6] <- AICc(seed_ff_g5t)
seed_ff_td_aic$Model[6] <- "m5"
seed_ff_td_aic$LL[6] <- logLik(seed_ff_g5t)
seed_ff_td_aic$AICc <- as.numeric(seed_ff_td_aic$AICc)
seed_ff_td_aic$LL <- as.numeric(seed_ff_td_aic$LL)
seed_ff_td_aic <- seed_ff_td_aic[order(seed_ff_td_aic$AICc), ]
seed_ff_td_aic$Delta_AICc <- "0.00"
seed_ff_td_aic$Delta_AICc[2] <- round(seed_ff_td_aic$AICc[1]-seed_ff_td_aic$AICc[2], 2)
seed_ff_td_aic$Delta_AICc[3] <- round(seed_ff_td_aic$AICc[1]-seed_ff_td_aic$AICc[3], 2)
seed_ff_td_aic$Delta_AICc[4] <- round(seed_ff_td_aic$AICc[1]-seed_ff_td_aic$AICc[4], 2)
seed_ff_td_aic$Delta_AICc[5] <- round(seed_ff_td_aic$AICc[1]-seed_ff_td_aic$AICc[5], 2)
seed_ff_td_aic$Delta_AICc[6] <- round(seed_ff_td_aic$AICc[1]-seed_ff_td_aic$AICc[6], 2)
seed_ff_td_aic  #Null model is best

seed_ff_ld_aic <- as.data.frame(1:6)
seed_ff_ld_aic$AICc <- "NA"
seed_ff_ld_aic$Model <- "NA"
seed_ff_ld_aic$LL <- "NA"
seed_ff_ld_aic$AICc[1] <- AICc(seed_ffl_gnull)
seed_ff_ld_aic$Model[1] <- "Null"
seed_ff_ld_aic$LL[1] <- logLik(seed_ffl_gnull)
seed_ff_ld_aic$AICc[2] <- AICc(seed_ff_g1l)
seed_ff_ld_aic$Model[2] <- "m1"
seed_ff_ld_aic$LL[2] <- logLik(seed_ff_g1l)
seed_ff_ld_aic$AICc[3] <- AICc(seed_ff_g2l)
seed_ff_ld_aic$Model[3] <- "m2"
seed_ff_ld_aic$LL[3] <- logLik(seed_ff_g2l)
seed_ff_ld_aic$AICc[4] <- AICc(seed_ff_g3l)
seed_ff_ld_aic$Model[4] <- "m3"
seed_ff_ld_aic$LL[4] <- logLik(seed_ff_g3l)
seed_ff_ld_aic$AICc[5] <- AICc(seed_ff_g4l)
seed_ff_ld_aic$Model[5] <- "m4"
seed_ff_ld_aic$LL[5] <- logLik(seed_ff_g4l)
seed_ff_ld_aic$AICc[6] <- AICc(seed_ff_g5l)
seed_ff_ld_aic$Model[6] <- "m5"
seed_ff_ld_aic$LL[6] <- logLik(seed_ff_g5l)
seed_ff_ld_aic$AICc <- as.numeric(seed_ff_ld_aic$AICc)
seed_ff_ld_aic$LL <- as.numeric(seed_ff_ld_aic$LL)
seed_ff_ld_aic <- seed_ff_ld_aic[order(seed_ff_ld_aic$AICc), ]
seed_ff_ld_aic$Delta_AICc <- "0.00"
seed_ff_ld_aic$Delta_AICc[2] <- round(seed_ff_ld_aic$AICc[1]-seed_ff_ld_aic$AICc[2], 2)
seed_ff_ld_aic$Delta_AICc[3] <- round(seed_ff_ld_aic$AICc[1]-seed_ff_ld_aic$AICc[3], 2)
seed_ff_ld_aic$Delta_AICc[4] <- round(seed_ff_ld_aic$AICc[1]-seed_ff_ld_aic$AICc[4], 2)
seed_ff_ld_aic$Delta_AICc[5] <- round(seed_ff_ld_aic$AICc[1]-seed_ff_ld_aic$AICc[5], 2)
seed_ff_ld_aic$Delta_AICc[6] <- round(seed_ff_ld_aic$AICc[1]-seed_ff_ld_aic$AICc[6], 2)
seed_ff_ld_aic  #Null model is best



# GAM s(bs = 'tp')
seed_ff_ts_aic <- as.data.frame(1:6)
seed_ff_ts_aic$AICc <- "NA"
seed_ff_ts_aic$Model <- "NA"
seed_ff_ts_aic$LL <- "NA"
seed_ff_ts_aic$AICc[1] <- AICc(seed_fft_gnull)
seed_ff_ts_aic$Model[1] <- "Null"
seed_ff_ts_aic$LL[1] <- logLik(seed_fft_gnull)
seed_ff_ts_aic$AICc[2] <- AICc(seed_ff_g1.1t)
seed_ff_ts_aic$Model[2] <- "m1"
seed_ff_ts_aic$LL[2] <- logLik(seed_ff_g1.1t)
seed_ff_ts_aic$AICc[3] <- AICc(seed_ff_g2.1t)
seed_ff_ts_aic$Model[3] <- "m2"
seed_ff_ts_aic$LL[3] <- logLik(seed_ff_g2.1t)
seed_ff_ts_aic$AICc[4] <- AICc(seed_ff_g3.1t)
seed_ff_ts_aic$Model[4] <- "m3"
seed_ff_ts_aic$LL[4] <- logLik(seed_ff_g3.1t)
seed_ff_ts_aic$AICc[5] <- AICc(seed_ff_g4.1t)
seed_ff_ts_aic$Model[5] <- "m4"
seed_ff_ts_aic$LL[5] <- logLik(seed_ff_g4.1t)
seed_ff_ts_aic$AICc[6] <- AICc(seed_ff_g5.1t)
seed_ff_ts_aic$Model[6] <- "m5"
seed_ff_ts_aic$LL[6] <- logLik(seed_ff_g5.1t)
seed_ff_ts_aic$AICc <- as.numeric(seed_ff_ts_aic$AICc)
seed_ff_ts_aic$LL <- as.numeric(seed_ff_ts_aic$LL)
seed_ff_ts_aic <- seed_ff_ts_aic[order(seed_ff_ts_aic$AICc), ]
seed_ff_ts_aic$Delta_AICc <- "0.00"
seed_ff_ts_aic$Delta_AICc[2] <- round(seed_ff_ts_aic$AICc[1]-seed_ff_ts_aic$AICc[2], 2)
seed_ff_ts_aic$Delta_AICc[3] <- round(seed_ff_ts_aic$AICc[1]-seed_ff_ts_aic$AICc[3], 2)
seed_ff_ts_aic$Delta_AICc[4] <- round(seed_ff_ts_aic$AICc[1]-seed_ff_ts_aic$AICc[4], 2)
seed_ff_ts_aic$Delta_AICc[5] <- round(seed_ff_ts_aic$AICc[1]-seed_ff_ts_aic$AICc[5], 2)
seed_ff_ts_aic$Delta_AICc[6] <- round(seed_ff_ts_aic$AICc[1]-seed_ff_ts_aic$AICc[6], 2)
seed_ff_ts_aic  #Null model is best

seed_ff_ls_aic <- as.data.frame(1:6)
seed_ff_ls_aic$AICc <- "NA"
seed_ff_ls_aic$Model <- "NA"
seed_ff_ls_aic$LL <- "NA"
seed_ff_ls_aic$AICc[1] <- AICc(seed_ffl_gnull)
seed_ff_ls_aic$Model[1] <- "Null"
seed_ff_ls_aic$LL[1] <- logLik(seed_ffl_gnull)
seed_ff_ls_aic$AICc[2] <- AICc(seed_ff_g1.1l)
seed_ff_ls_aic$Model[2] <- "m1"
seed_ff_ls_aic$LL[2] <- logLik(seed_ff_g1.1l)
seed_ff_ls_aic$AICc[3] <- AICc(seed_ff_g2.1l)
seed_ff_ls_aic$Model[3] <- "m2"
seed_ff_ls_aic$LL[3] <- logLik(seed_ff_g2.1l)
seed_ff_ls_aic$AICc[4] <- AICc(seed_ff_g3.1l)
seed_ff_ls_aic$Model[4] <- "m3"
seed_ff_ls_aic$LL[4] <- logLik(seed_ff_g3.1l)
seed_ff_ls_aic$AICc[5] <- AICc(seed_ff_g4.1l)
seed_ff_ls_aic$Model[5] <- "m4"
seed_ff_ls_aic$LL[5] <- logLik(seed_ff_g4.1l)
seed_ff_ls_aic$AICc[6] <- AICc(seed_ff_g5.1l)
seed_ff_ls_aic$Model[6] <- "m5"
seed_ff_ls_aic$LL[6] <- logLik(seed_ff_g5.1l)
seed_ff_ls_aic$AICc <- as.numeric(seed_ff_ls_aic$AICc)
seed_ff_ls_aic$LL <- as.numeric(seed_ff_ls_aic$LL)
seed_ff_ls_aic <- seed_ff_ls_aic[order(seed_ff_ls_aic$AICc), ]
seed_ff_ls_aic$Delta_AICc <- "0.00"
seed_ff_ls_aic$Delta_AICc[2] <- round(seed_ff_ls_aic$AICc[1]-seed_ff_ls_aic$AICc[2], 2)
seed_ff_ls_aic$Delta_AICc[3] <- round(seed_ff_ls_aic$AICc[1]-seed_ff_ls_aic$AICc[3], 2)
seed_ff_ls_aic$Delta_AICc[4] <- round(seed_ff_ls_aic$AICc[1]-seed_ff_ls_aic$AICc[4], 2)
seed_ff_ls_aic$Delta_AICc[5] <- round(seed_ff_ls_aic$AICc[1]-seed_ff_ls_aic$AICc[5], 2)
seed_ff_ls_aic$Delta_AICc[6] <- round(seed_ff_ls_aic$AICc[1]-seed_ff_ls_aic$AICc[6], 2)
seed_ff_ls_aic  #Null model is best


# GAM ti(bs='cr')
seed_ff_tsc_aic <- as.data.frame(1:6)
seed_ff_tsc_aic$AICc <- "NA"
seed_ff_tsc_aic$Model <- "NA"
seed_ff_tsc_aic$LL <- "NA"
seed_ff_tsc_aic$AICc[1] <- AICc(seed_fft_gnull)
seed_ff_tsc_aic$Model[1] <- "Null"
seed_ff_tsc_aic$LL[1] <- logLik(seed_fft_gnull)
seed_ff_tsc_aic$AICc[2] <- AICc(seed_ff_g1t)
seed_ff_tsc_aic$Model[2] <- "m1"
seed_ff_tsc_aic$LL[2] <- logLik(seed_ff_g1t)
seed_ff_tsc_aic$AICc[3] <- AICc(seed_ff_g2.2t)
seed_ff_tsc_aic$Model[3] <- "m2"
seed_ff_tsc_aic$LL[3] <- logLik(seed_ff_g2.2t)
seed_ff_tsc_aic$AICc[4] <- AICc(seed_ff_g3.2t)
seed_ff_tsc_aic$Model[4] <- "m3"
seed_ff_tsc_aic$LL[4] <- logLik(seed_ff_g3.2t)
seed_ff_tsc_aic$AICc[5] <- AICc(seed_ff_g4.2t)
seed_ff_tsc_aic$Model[5] <- "m4"
seed_ff_tsc_aic$LL[5] <- logLik(seed_ff_g4.2t)
seed_ff_tsc_aic$AICc[6] <- AICc(seed_ff_g5.2t)
seed_ff_tsc_aic$Model[6] <- "m5"
seed_ff_tsc_aic$LL[6] <- logLik(seed_ff_g5.2t)
seed_ff_tsc_aic$AICc <- as.numeric(seed_ff_tsc_aic$AICc)
seed_ff_tsc_aic$LL <- as.numeric(seed_ff_tsc_aic$LL)
seed_ff_tsc_aic <- seed_ff_tsc_aic[order(seed_ff_tsc_aic$AICc), ]
seed_ff_tsc_aic$Delta_AICc <- "0.00"
seed_ff_tsc_aic$Delta_AICc[2] <- round(seed_ff_tsc_aic$AICc[1]-seed_ff_tsc_aic$AICc[2], 2)
seed_ff_tsc_aic$Delta_AICc[3] <- round(seed_ff_tsc_aic$AICc[1]-seed_ff_tsc_aic$AICc[3], 2)
seed_ff_tsc_aic$Delta_AICc[4] <- round(seed_ff_tsc_aic$AICc[1]-seed_ff_tsc_aic$AICc[4], 2)
seed_ff_tsc_aic$Delta_AICc[5] <- round(seed_ff_tsc_aic$AICc[1]-seed_ff_tsc_aic$AICc[5], 2)
seed_ff_tsc_aic$Delta_AICc[6] <- round(seed_ff_tsc_aic$AICc[1]-seed_ff_tsc_aic$AICc[6], 2)
seed_ff_tsc_aic  #Null model is best

seed_ff_lsc_aic <- as.data.frame(1:6)
seed_ff_lsc_aic$AICc <- "NA"
seed_ff_lsc_aic$Model <- "NA"
seed_ff_lsc_aic$LL <- "NA"
seed_ff_lsc_aic$AICc[1] <- AICc(seed_ffl_gnull)
seed_ff_lsc_aic$Model[1] <- "Null"
seed_ff_lsc_aic$LL[1] <- logLik(seed_ffl_gnull)
seed_ff_lsc_aic$AICc[2] <- AICc(seed_ff_g1l)
seed_ff_lsc_aic$Model[2] <- "m1"
seed_ff_lsc_aic$LL[2] <- logLik(seed_ff_g1l)
seed_ff_lsc_aic$AICc[3] <- AICc(seed_ff_g2.2l)
seed_ff_lsc_aic$Model[3] <- "m2"
seed_ff_lsc_aic$LL[3] <- logLik(seed_ff_g2.2l)
seed_ff_lsc_aic$AICc[4] <- AICc(seed_ff_g3.2l)
seed_ff_lsc_aic$Model[4] <- "m3"
seed_ff_lsc_aic$LL[4] <- logLik(seed_ff_g3.2l)
seed_ff_lsc_aic$AICc[5] <- AICc(seed_ff_g4.2l)
seed_ff_lsc_aic$Model[5] <- "m4"
seed_ff_lsc_aic$LL[5] <- logLik(seed_ff_g4.2l)
seed_ff_lsc_aic$AICc[6] <- AICc(seed_ff_g5.2l)
seed_ff_lsc_aic$Model[6] <- "m5"
seed_ff_lsc_aic$LL[6] <- logLik(seed_ff_g5.2l)
seed_ff_lsc_aic$AICc <- as.numeric(seed_ff_lsc_aic$AICc)
seed_ff_lsc_aic$LL <- as.numeric(seed_ff_lsc_aic$LL)
seed_ff_lsc_aic <- seed_ff_lsc_aic[order(seed_ff_lsc_aic$AICc), ]
seed_ff_lsc_aic$Delta_AICc <- "0.00"
seed_ff_lsc_aic$Delta_AICc[2] <- round(seed_ff_lsc_aic$AICc[1]-seed_ff_lsc_aic$AICc[2], 2)
seed_ff_lsc_aic$Delta_AICc[3] <- round(seed_ff_lsc_aic$AICc[1]-seed_ff_lsc_aic$AICc[3], 2)
seed_ff_lsc_aic$Delta_AICc[4] <- round(seed_ff_lsc_aic$AICc[1]-seed_ff_lsc_aic$AICc[4], 2)
seed_ff_lsc_aic$Delta_AICc[5] <- round(seed_ff_lsc_aic$AICc[1]-seed_ff_lsc_aic$AICc[5], 2)
seed_ff_lsc_aic$Delta_AICc[6] <- round(seed_ff_lsc_aic$AICc[1]-seed_ff_lsc_aic$AICc[6], 2)
seed_ff_lsc_aic  #Null model is best



# 6.1.2 Proportion saplings ----
sap_ff_tnull <- glmer(Proportion_saplings ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
sap_ff_gt_null <- gam(Proportion_saplings ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')

sap_ff_lnull <- glmer(Proportion_saplings ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)
sap_ff_gl_null <- gam(Proportion_saplings ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')



sap_ff_m1t <- glmer(Proportion_saplings ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m1t)
sap_ff_g1t <- gam(Proportion_saplings ~ s(Fire_freq, k = 8), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g1t)
plot(sap_ff_g1t)
sap_ff_g1.1t <- gam(Proportion_saplings ~ s(Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g1.1t)
plot(sap_ff_g1.1t)

sap_ff_m1l <- glmer(Proportion_saplings ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(sap_ff_m1l) # Not a very good model
sap_ff_g1l <- gam(Proportion_saplings ~ s(Fire_freq, k = 3), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g1l)
plot(sap_ff_g1l)
sap_ff_g1.1l <- gam(Proportion_saplings ~ s(Fire_freq, bs = 'cr', k = 3), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g1.1l)
plot(sap_ff_g1.1l)


sap_ff_m2t <- glmer(Proportion_saplings ~ r_Fire_freq * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m2t)
sap_ff_g2t <- gam(Proportion_saplings ~ s(Fire_freq, k = 8) + s(Latitude) + ti(Latitude, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g2t)
plot(sap_ff_g2t)
sap_ff_g2.1t <- gam(Proportion_saplings ~ s(Fire_freq, k = 8) + s(Latitude) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g2.1t)
plot(sap_ff_g2.1t)
sap_ff_g2.2t <- gam(Proportion_saplings ~ s(Fire_freq, bs = 'cr', k = 5) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g2.1t)
plot(sap_ff_g2.2t)


#sap_ff_m2l <- glmer(Proportion_saplings ~ Fire_freq * Latitude + (1|Location/Transect), family = binomial, data = lit_transects) # This model fails due to the low number of points available to model an interaction. We can only model an additive effect of Fire_freq and latitude
sap_ff_m2.1l <- glmer(Proportion_saplings ~ r_Fire_freq + r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=80)))
summary(sap_ff_m2.1l) # Not very good model
sap_ff_g2l <- gam(Proportion_saplings ~ s(Fire_freq, k = 3) + s(Latitude, k = 8) + ti(Latitude, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g2l)
plot(sap_ff_g2l)
sap_ff_g2.1l <- gam(Proportion_saplings ~ s(Fire_freq, k = 3) + s(Latitude, k = 8) + ti(Latitude, by = Fire_freq, bs = 'tp'), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g2.1l)
plot(sap_ff_g2.1l)
sap_ff_g2.2l <- gam(Proportion_saplings ~ s(Fire_freq, bs = 'cr', k = 3) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g2.1l)
plot(sap_ff_g2.2l)



sap_ff_m3t <- glmer(Proportion_saplings ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m3t)
sap_ff_g3t <- gam(Proportion_saplings ~ s(Fire_freq, k = 8) + s(FPC) + ti(FPC, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g3t)
plot(sap_ff_g3t)
sap_ff_g3.1t <- gam(Proportion_saplings ~ s(Fire_freq, k = 8) + s(FPC) + ti(FPC, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g3.1t)
plot(sap_ff_g3.1t)
sap_ff_g3.2t <- gam(Proportion_saplings ~ s(Fire_freq, bs = 'cr', k = 5) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g3.1t)
plot(sap_ff_g3.2t)


sap_ff_m3l <- glmer(Proportion_saplings ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects)
summary(sap_ff_m3l)
sap_ff_g3l <- gam(Proportion_saplings ~ s(Fire_freq, k = 3) + s(FPC, k = 7) + ti(FPC, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g3l)
plot(sap_ff_g3l)
sap_ff_g3.1l <- gam(Proportion_saplings ~ s(Fire_freq, k = 3) + s(FPC, k = 7) + ti(FPC, by = Fire_freq, bs = 'tp', k = =7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g3.1l)
plot(sap_ff_g3.1l)
sap_ff_g3.2l <- gam(Proportion_saplings ~ s(Fire_freq, bs = 'cr', k = 3) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g3.2l)
plot(sap_ff_g3.2l)


sap_ff_m4t <- glmer(Proportion_saplings ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m4t)
sap_ff_g4t <- gam(Proportion_saplings ~ s(Fire_freq, k = 8) + s(Precip) + ti(Precip, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g4t)
plot(sap_ff_g4t)
sap_ff_g4.1t <- gam(Proportion_saplings ~ s(Fire_freq, k = 8) + s(Precip) + ti(Precip, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g4.1t)
plot(sap_ff_g4.1t)
sap_ff_g4.2t <- gam(Proportion_saplings ~ s(Fire_freq, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g4.1t)
plot(sap_ff_g4.2t)


sap_ff_m4l <- glmer(Proportion_saplings ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(sap_ff_m4l)
sap_ff_g4l <- gam(Proportion_saplings ~ s(Fire_freq, k = 3) + s(Precip, k = 7) + ti(Precip, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g4l)
plot(sap_ff_g4l)
sap_ff_g4.1l <- gam(Proportion_saplings ~ s(Fire_freq, k = 3) + s(Precip, k = 7) + ti(Precip, by = Fire_freq, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g4.1l)
plot(sap_ff_g4.1l)
sap_ff_g4.2l <- gam(Proportion_saplings ~ s(Fire_freq, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g4.2l)
plot(sap_ff_g4.2l)


sap_ff_m5t <- glmer(Proportion_saplings ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m5t)

sap_ff_m5ta <- glmer(Proportion_saplings ~ r_Fire_freq + r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
sap_ff_m5t2 <- glmer(Proportion_saplings ~ r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)

sap_ff_g5t <- gam(Proportion_saplings ~ s(Fire_freq, k = 8) + s(Temp) + ti(Temp, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g5t)
plot(sap_ff_g5t)
sap_ff_g5.1t <- gam(Proportion_saplings ~ s(Fire_freq, k = 8) + s(Temp) + ti(Temp, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g5.1t)
plot(sap_ff_g5.1t)
sap_ff_g5.2t <- gam(Proportion_saplings ~ s(Fire_freq, bs = 'cr', k = 5) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g5.2t)
plot(sap_ff_g5.2t)
sap_ff_g5.3t <- gam(Proportion_saplings ~ s(Fire_freq, bs = 'cr', k = 5) + s(Temp, bs = 'cc', k = 4) + ti(Temp, by = Fire_freq, bs = 'cc', k = 4), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g5.3t)
plot(sap_ff_g5.3t)

sap_ff_m5l <- glmer(Proportion_saplings ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects)
summary(sap_ff_m5l)
sap_ff_g5l <- gam(Proportion_saplings ~ s(Fire_freq, k = 3) + s(Temp, k = 7) + ti(Temp, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g5l)
plot(sap_ff_g5l)
sap_ff_g5.1l <- gam(Proportion_saplings ~ s(Fire_freq, k = 3) + s(Temp, k = 7) + ti(Temp, by = Fire_freq, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g5.1l)
plot(sap_ff_g5.1l)
sap_ff_g5.2l <- gam(Proportion_saplings ~ s(Fire_freq, bs = 'cr', k = 3) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(sap_ff_g5.2l)
plot(sap_ff_g5.2l)

# Compare AICs for each modelling type
# GLMER
sap_ff_aic_t <- list(sap_ff_tnull, sap_ff_m1t, sap_ff_m2t, sap_ff_m3t, sap_ff_m4t, sap_ff_m5t)
aictab(sap_ff_aic_t) # Model 5 is best but null model is ranked within delta AICc <2. 

sap_ff_aic_l <- list(sap_ff_lnull, sap_ff_m1l, sap_ff_m2.1l, sap_ff_m3l, sap_ff_m4l, sap_ff_m5l)
aictab(sap_ff_aic_l) # Null model is best

# GAM defaults
sap_ff_td_aic <- as.data.frame(1:6)
sap_ff_td_aic$AICc <- "NA"
sap_ff_td_aic$Model <- "NA"
sap_ff_td_aic$LL <- "NA"
sap_ff_td_aic$AICc[1] <- AICc(sap_ff_gt_null)
sap_ff_td_aic$Model[1] <- "Null"
sap_ff_td_aic$LL[1] <- logLik(sap_ff_gt_null)
sap_ff_td_aic$AICc[2] <- AICc(sap_ff_g1t)
sap_ff_td_aic$Model[2] <- "m1"
sap_ff_td_aic$LL[2] <- logLik(sap_ff_g1t)
sap_ff_td_aic$AICc[3] <- AICc(sap_ff_g2t)
sap_ff_td_aic$Model[3] <- "m2"
sap_ff_td_aic$LL[3] <- logLik(sap_ff_g2t)
sap_ff_td_aic$AICc[4] <- AICc(sap_ff_g3t)
sap_ff_td_aic$Model[4] <- "m3"
sap_ff_td_aic$LL[4] <- logLik(sap_ff_g3t)
sap_ff_td_aic$AICc[5] <- AICc(sap_ff_g4t)
sap_ff_td_aic$Model[5] <- "m4"
sap_ff_td_aic$LL[5] <- logLik(sap_ff_g4t)
sap_ff_td_aic$AICc[6] <- AICc(sap_ff_g5t)
sap_ff_td_aic$Model[6] <- "m5"
sap_ff_td_aic$LL[6] <- logLik(sap_ff_g5t)
sap_ff_td_aic$AICc <- as.numeric(sap_ff_td_aic$AICc)
sap_ff_td_aic$LL <- as.numeric(sap_ff_td_aic$LL)
sap_ff_td_aic <- sap_ff_td_aic[order(sap_ff_td_aic$AICc), ]
sap_ff_td_aic$Delta_AICc <- "0.00"
sap_ff_td_aic$Delta_AICc[2] <- round(sap_ff_td_aic$AICc[1]-sap_ff_td_aic$AICc[2], 2)
sap_ff_td_aic$Delta_AICc[3] <- round(sap_ff_td_aic$AICc[1]-sap_ff_td_aic$AICc[3], 2)
sap_ff_td_aic$Delta_AICc[4] <- round(sap_ff_td_aic$AICc[1]-sap_ff_td_aic$AICc[4], 2)
sap_ff_td_aic$Delta_AICc[5] <- round(sap_ff_td_aic$AICc[1]-sap_ff_td_aic$AICc[5], 2)
sap_ff_td_aic$Delta_AICc[6] <- round(sap_ff_td_aic$AICc[1]-sap_ff_td_aic$AICc[6], 2)
sap_ff_td_aic  # Model 5 is best but null model is ranked within delta AICc <2. 

sap_ff_ld_aic <- as.data.frame(1:6)
sap_ff_ld_aic$AICc <- "NA"
sap_ff_ld_aic$Model <- "NA"
sap_ff_ld_aic$LL <- "NA"
sap_ff_ld_aic$AICc[1] <- AICc(sap_ff_gl_null)
sap_ff_ld_aic$Model[1] <- "Null"
sap_ff_ld_aic$LL[1] <- logLik(sap_ff_gl_null)
sap_ff_ld_aic$AICc[2] <- AICc(sap_ff_g1l)
sap_ff_ld_aic$Model[2] <- "m1"
sap_ff_ld_aic$LL[2] <- logLik(sap_ff_g1l)
sap_ff_ld_aic$AICc[3] <- AICc(sap_ff_g2l)
sap_ff_ld_aic$Model[3] <- "m2"
sap_ff_ld_aic$LL[3] <- logLik(sap_ff_g2l)
sap_ff_ld_aic$AICc[4] <- AICc(sap_ff_g3l)
sap_ff_ld_aic$Model[4] <- "m3"
sap_ff_ld_aic$LL[4] <- logLik(sap_ff_g3l)
sap_ff_ld_aic$AICc[5] <- AICc(sap_ff_g4l)
sap_ff_ld_aic$Model[5] <- "m4"
sap_ff_ld_aic$LL[5] <- logLik(sap_ff_g4l)
sap_ff_ld_aic$AICc[6] <- AICc(sap_ff_g5l)
sap_ff_ld_aic$Model[6] <- "m5"
sap_ff_ld_aic$LL[6] <- logLik(sap_ff_g5l)
sap_ff_ld_aic$AICc <- as.numeric(sap_ff_ld_aic$AICc)
sap_ff_ld_aic$LL <- as.numeric(sap_ff_ld_aic$LL)
sap_ff_ld_aic <- sap_ff_ld_aic[order(sap_ff_ld_aic$AICc), ]
sap_ff_ld_aic$Delta_AICc <- "0.00"
sap_ff_ld_aic$Delta_AICc[2] <- round(sap_ff_ld_aic$AICc[1]-sap_ff_ld_aic$AICc[2], 2)
sap_ff_ld_aic$Delta_AICc[3] <- round(sap_ff_ld_aic$AICc[1]-sap_ff_ld_aic$AICc[3], 2)
sap_ff_ld_aic$Delta_AICc[4] <- round(sap_ff_ld_aic$AICc[1]-sap_ff_ld_aic$AICc[4], 2)
sap_ff_ld_aic$Delta_AICc[5] <- round(sap_ff_ld_aic$AICc[1]-sap_ff_ld_aic$AICc[5], 2)
sap_ff_ld_aic$Delta_AICc[6] <- round(sap_ff_ld_aic$AICc[1]-sap_ff_ld_aic$AICc[6], 2)
sap_ff_ld_aic  #Null model is best



# GAM s(bs = 'tp')
sap_ff_ts_aic <- as.data.frame(1:6)
sap_ff_ts_aic$AICc <- "NA"
sap_ff_ts_aic$Model <- "NA"
sap_ff_ts_aic$LL <- "NA"
sap_ff_ts_aic$AICc[1] <- AICc(sap_ff_gt_null)
sap_ff_ts_aic$Model[1] <- "Null"
sap_ff_ts_aic$LL[1] <- logLik(sap_ff_gt_null)
sap_ff_ts_aic$AICc[2] <- AICc(sap_ff_g1.1t)
sap_ff_ts_aic$Model[2] <- "m1"
sap_ff_ts_aic$LL[2] <- logLik(sap_ff_g1.1t)
sap_ff_ts_aic$AICc[3] <- AICc(sap_ff_g2.1t)
sap_ff_ts_aic$Model[3] <- "m2"
sap_ff_ts_aic$LL[3] <- logLik(sap_ff_g2.1t)
sap_ff_ts_aic$AICc[4] <- AICc(sap_ff_g3.1t)
sap_ff_ts_aic$Model[4] <- "m3"
sap_ff_ts_aic$LL[4] <- logLik(sap_ff_g3.1t)
sap_ff_ts_aic$AICc[5] <- AICc(sap_ff_g4.1t)
sap_ff_ts_aic$Model[5] <- "m4"
sap_ff_ts_aic$LL[5] <- logLik(sap_ff_g4.1t)
sap_ff_ts_aic$AICc[6] <- AICc(sap_ff_g5.1t)
sap_ff_ts_aic$Model[6] <- "m5"
sap_ff_ts_aic$LL[6] <- logLik(sap_ff_g5.1t)
sap_ff_ts_aic$AICc <- as.numeric(sap_ff_ts_aic$AICc)
sap_ff_ts_aic$LL <- as.numeric(sap_ff_ts_aic$LL)
sap_ff_ts_aic <- sap_ff_ts_aic[order(sap_ff_ts_aic$AICc), ]
sap_ff_ts_aic$Delta_AICc <- "0.00"
sap_ff_ts_aic$Delta_AICc[2] <- round(sap_ff_ts_aic$AICc[1]-sap_ff_ts_aic$AICc[2], 2)
sap_ff_ts_aic$Delta_AICc[3] <- round(sap_ff_ts_aic$AICc[1]-sap_ff_ts_aic$AICc[3], 2)
sap_ff_ts_aic$Delta_AICc[4] <- round(sap_ff_ts_aic$AICc[1]-sap_ff_ts_aic$AICc[4], 2)
sap_ff_ts_aic$Delta_AICc[5] <- round(sap_ff_ts_aic$AICc[1]-sap_ff_ts_aic$AICc[5], 2)
sap_ff_ts_aic$Delta_AICc[6] <- round(sap_ff_ts_aic$AICc[1]-sap_ff_ts_aic$AICc[6], 2)
sap_ff_ts_aic  # Model 5 is best but null model is ranked within delta AICc <2. 

sap_ff_ls_aic <- as.data.frame(1:6)
sap_ff_ls_aic$AICc <- "NA"
sap_ff_ls_aic$Model <- "NA"
sap_ff_ls_aic$LL <- "NA"
sap_ff_ls_aic$AICc[1] <- AICc(sap_ff_gl_null)
sap_ff_ls_aic$Model[1] <- "Null"
sap_ff_ls_aic$LL[1] <- logLik(sap_ff_gl_null)
sap_ff_ls_aic$AICc[2] <- AICc(sap_ff_g1.1l)
sap_ff_ls_aic$Model[2] <- "m1"
sap_ff_ls_aic$LL[2] <- logLik(sap_ff_g1.1l)
sap_ff_ls_aic$AICc[3] <- AICc(sap_ff_g2.1l)
sap_ff_ls_aic$Model[3] <- "m2"
sap_ff_ls_aic$LL[3] <- logLik(sap_ff_g2.1l)
sap_ff_ls_aic$AICc[4] <- AICc(sap_ff_g3.1l)
sap_ff_ls_aic$Model[4] <- "m3"
sap_ff_ls_aic$LL[4] <- logLik(sap_ff_g3.1l)
sap_ff_ls_aic$AICc[5] <- AICc(sap_ff_g4.1l)
sap_ff_ls_aic$Model[5] <- "m4"
sap_ff_ls_aic$LL[5] <- logLik(sap_ff_g4.1l)
sap_ff_ls_aic$AICc[6] <- AICc(sap_ff_g5.1l)
sap_ff_ls_aic$Model[6] <- "m5"
sap_ff_ls_aic$LL[6] <- logLik(sap_ff_g5.1l)
sap_ff_ls_aic$AICc <- as.numeric(sap_ff_ls_aic$AICc)
sap_ff_ls_aic$LL <- as.numeric(sap_ff_ls_aic$LL)
sap_ff_ls_aic <- sap_ff_ls_aic[order(sap_ff_ls_aic$AICc), ]
sap_ff_ls_aic$Delta_AICc <- "0.00"
sap_ff_ls_aic$Delta_AICc[2] <- round(sap_ff_ls_aic$AICc[1]-sap_ff_ls_aic$AICc[2], 2)
sap_ff_ls_aic$Delta_AICc[3] <- round(sap_ff_ls_aic$AICc[1]-sap_ff_ls_aic$AICc[3], 2)
sap_ff_ls_aic$Delta_AICc[4] <- round(sap_ff_ls_aic$AICc[1]-sap_ff_ls_aic$AICc[4], 2)
sap_ff_ls_aic$Delta_AICc[5] <- round(sap_ff_ls_aic$AICc[1]-sap_ff_ls_aic$AICc[5], 2)
sap_ff_ls_aic$Delta_AICc[6] <- round(sap_ff_ls_aic$AICc[1]-sap_ff_ls_aic$AICc[6], 2)
sap_ff_ls_aic  #Null model is best


# GAM ti(bs='cr')
sap_ff_tsc_aic <- as.data.frame(1:6)
sap_ff_tsc_aic$AICc <- "NA"
sap_ff_tsc_aic$Model <- "NA"
sap_ff_tsc_aic$LL <- "NA"
sap_ff_tsc_aic$AICc[1] <- AICc(sap_ff_gt_null)
sap_ff_tsc_aic$Model[1] <- "Null"
sap_ff_tsc_aic$LL[1] <- logLik(sap_ff_gt_null)
sap_ff_tsc_aic$AICc[2] <- AICc(sap_ff_g1t)
sap_ff_tsc_aic$Model[2] <- "m1"
sap_ff_tsc_aic$LL[2] <- logLik(sap_ff_g1t)
sap_ff_tsc_aic$AICc[3] <- AICc(sap_ff_g2.2t)
sap_ff_tsc_aic$Model[3] <- "m2"
sap_ff_tsc_aic$LL[3] <- logLik(sap_ff_g2.2t)
sap_ff_tsc_aic$AICc[4] <- AICc(sap_ff_g3.2t)
sap_ff_tsc_aic$Model[4] <- "m3"
sap_ff_tsc_aic$LL[4] <- logLik(sap_ff_g3.2t)
sap_ff_tsc_aic$AICc[5] <- AICc(sap_ff_g4.2t)
sap_ff_tsc_aic$Model[5] <- "m4"
sap_ff_tsc_aic$LL[5] <- logLik(sap_ff_g4.2t)
sap_ff_tsc_aic$AICc[6] <- AICc(sap_ff_g5.2t)
sap_ff_tsc_aic$Model[6] <- "m5"
sap_ff_tsc_aic$LL[6] <- logLik(sap_ff_g5.2t)
sap_ff_tsc_aic$AICc <- as.numeric(sap_ff_tsc_aic$AICc)
sap_ff_tsc_aic$LL <- as.numeric(sap_ff_tsc_aic$LL)
sap_ff_tsc_aic <- sap_ff_tsc_aic[order(sap_ff_tsc_aic$AICc), ]
sap_ff_tsc_aic$Delta_AICc <- "0.00"
sap_ff_tsc_aic$Delta_AICc[2] <- round(sap_ff_tsc_aic$AICc[1]-sap_ff_tsc_aic$AICc[2], 2)
sap_ff_tsc_aic$Delta_AICc[3] <- round(sap_ff_tsc_aic$AICc[1]-sap_ff_tsc_aic$AICc[3], 2)
sap_ff_tsc_aic$Delta_AICc[4] <- round(sap_ff_tsc_aic$AICc[1]-sap_ff_tsc_aic$AICc[4], 2)
sap_ff_tsc_aic$Delta_AICc[5] <- round(sap_ff_tsc_aic$AICc[1]-sap_ff_tsc_aic$AICc[5], 2)
sap_ff_tsc_aic$Delta_AICc[6] <- round(sap_ff_tsc_aic$AICc[1]-sap_ff_tsc_aic$AICc[6], 2)
sap_ff_tsc_aic  # Model 5 is best but null model is ranked within delta AICc <2. 

sap_ff_lsc_aic <- as.data.frame(1:6)
sap_ff_lsc_aic$AICc <- "NA"
sap_ff_lsc_aic$Model <- "NA"
sap_ff_lsc_aic$LL <- "NA"
sap_ff_lsc_aic$AICc[1] <- AICc(sap_ff_gl_null)
sap_ff_lsc_aic$Model[1] <- "Null"
sap_ff_lsc_aic$LL[1] <- logLik(sap_ff_gl_null)
sap_ff_lsc_aic$AICc[2] <- AICc(sap_ff_g1l)
sap_ff_lsc_aic$Model[2] <- "m1"
sap_ff_lsc_aic$LL[2] <- logLik(sap_ff_g1l)
sap_ff_lsc_aic$AICc[3] <- AICc(sap_ff_g2.2l)
sap_ff_lsc_aic$Model[3] <- "m2"
sap_ff_lsc_aic$LL[3] <- logLik(sap_ff_g2.2l)
sap_ff_lsc_aic$AICc[4] <- AICc(sap_ff_g3.2l)
sap_ff_lsc_aic$Model[4] <- "m3"
sap_ff_lsc_aic$LL[4] <- logLik(sap_ff_g3.2l)
sap_ff_lsc_aic$AICc[5] <- AICc(sap_ff_g4.2l)
sap_ff_lsc_aic$Model[5] <- "m4"
sap_ff_lsc_aic$LL[5] <- logLik(sap_ff_g4.2l)
sap_ff_lsc_aic$AICc[6] <- AICc(sap_ff_g5.2l)
sap_ff_lsc_aic$Model[6] <- "m5"
sap_ff_lsc_aic$LL[6] <- logLik(sap_ff_g5.2l)
sap_ff_lsc_aic$AICc <- as.numeric(sap_ff_lsc_aic$AICc)
sap_ff_lsc_aic$LL <- as.numeric(sap_ff_lsc_aic$LL)
sap_ff_lsc_aic <- sap_ff_lsc_aic[order(sap_ff_lsc_aic$AICc), ]
sap_ff_lsc_aic$Delta_AICc <- "0.00"
sap_ff_lsc_aic$Delta_AICc[2] <- round(sap_ff_lsc_aic$AICc[1]-sap_ff_lsc_aic$AICc[2], 2)
sap_ff_lsc_aic$Delta_AICc[3] <- round(sap_ff_lsc_aic$AICc[1]-sap_ff_lsc_aic$AICc[3], 2)
sap_ff_lsc_aic$Delta_AICc[4] <- round(sap_ff_lsc_aic$AICc[1]-sap_ff_lsc_aic$AICc[4], 2)
sap_ff_lsc_aic$Delta_AICc[5] <- round(sap_ff_lsc_aic$AICc[1]-sap_ff_lsc_aic$AICc[5], 2)
sap_ff_lsc_aic$Delta_AICc[6] <- round(sap_ff_lsc_aic$AICc[1]-sap_ff_lsc_aic$AICc[6], 2)
sap_ff_lsc_aic  #Null model is best




# Predict for proportion saplings data ----
new_ff_m5l <-  expand.grid(Fire_freq = min(tor_transects$Fire_freq),
                        Temp = seq(min(tor_transects$Temp), max(tor_transects$Temp), length = 50))
new_ff_m5a <- expand.grid(Fire_freq = mean(tor_transects$Fire_freq),
                        Temp = seq(min(tor_transects$Temp), max(tor_transects$Temp), length = 50))
new_ff_m5h <- expand.grid(Fire_freq = max(tor_transects$Fire_freq),
                        Temp = seq(min(tor_transects$Temp), max(tor_transects$Temp), length = 50))



# GLMER 
new_ff_m5l_lm <-  expand.grid(r_Fire_freq = min(tor_transects$r_Fire_freq),
                           r_Temp = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 50))
new_ff_m5a_lm <- expand.grid(r_Fire_freq = mean(tor_transects$Fire_freq),
                          r_Temp = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 50))
new_ff_m5h_lm <- expand.grid(r_Fire_freq = max(tor_transects$Fire_freq),
                          r_Temp = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 50))



nsapff_m5l_lm <- new_ff_m5l_lm
nsapff_m5a_lm <- new_ff_m5a_lm
nsapff_m5h_lm <- new_ff_m5h_lm

psapff_m5l_lm <- predictSE(sap_ff_m5t, newdata = nsapff_m5l_lm, se.fit = T, type = 'link')
nsapff_m5l_lm$fit.link <- psapff_m5l_lm$fit
nsapff_m5l_lm$se.link <- psapff_m5l_lm$se.fit
nsapff_m5l_lm$lci.link <- nsapff_m5l_lm$fit.link - (nsapff_m5l_lm$se.link * 1.96)
nsapff_m5l_lm$uci.link <- nsapff_m5l_lm$fit.link + (nsapff_m5l_lm$se.link * 1.96)
nsapff_m5l_lm$fit <- invlogit(nsapff_m5l_lm$fit.link)
nsapff_m5l_lm$se <- invlogit(nsapff_m5l_lm$se.link)
nsapff_m5l_lm$lci <- invlogit(nsapff_m5l_lm$lci.link)
nsapff_m5l_lm$uci <- invlogit(nsapff_m5l_lm$uci.link)

psapff_m5a_lm <- predictSE(sap_ff_m5t, newdata = nsapff_m5a_lm, se.fit = T, type = 'link')
nsapff_m5a_lm$fit.link <- psapff_m5a_lm$fit
nsapff_m5a_lm$se.link <- psapff_m5a_lm$se.fit
nsapff_m5a_lm$lci.link <- nsapff_m5a_lm$fit.link - (nsapff_m5a_lm$se.link * 1.96)
nsapff_m5a_lm$uci.link <- nsapff_m5a_lm$fit.link + (nsapff_m5a_lm$se.link * 1.96)
nsapff_m5a_lm$fit <- invlogit(nsapff_m5a_lm$fit.link)
nsapff_m5a_lm$se <- invlogit(nsapff_m5a_lm$se.link)
nsapff_m5a_lm$lci <- invlogit(nsapff_m5a_lm$lci.link)
nsapff_m5a_lm$uci <- invlogit(nsapff_m5a_lm$uci.link)

psapff_m5h_lm <- predictSE(sap_ff_m5t, newdata = nsapff_m5h_lm, se.fit = T, type = 'link')
nsapff_m5h_lm$fit.link <- psapff_m5h_lm$fit
nsapff_m5h_lm$se.link <- psapff_m5h_lm$se.fit
nsapff_m5h_lm$lci.link <- nsapff_m5h_lm$fit.link - (nsapff_m5h_lm$se.link * 1.96)
nsapff_m5h_lm$uci.link <- nsapff_m5h_lm$fit.link + (nsapff_m5h_lm$se.link * 1.96)
nsapff_m5h_lm$fit <- invlogit(nsapff_m5h_lm$fit.link)
nsapff_m5h_lm$se <- invlogit(nsapff_m5h_lm$se.link)
nsapff_m5h_lm$lci <- invlogit(nsapff_m5h_lm$lci.link)
nsapff_m5h_lm$uci <- invlogit(nsapff_m5h_lm$uci.link)




# GAM default
nsapff_m5l_d <- new_ff_m5l
nsapff_m5a_d <- new_ff_m5a
nsapff_m5h_d <- new_ff_m5h

psapff_m5l_d <- predict(sap_ff_g5t, newdata = nsapff_m5l_d, se.fit = T, type = 'response')
nsapff_m5l_d$fit <- psapff_m5l_d$fit
nsapff_m5l_d$se <- psapff_m5l_d$se.fit
nsapff_m5l_d$lci <- nsapff_m5l_d$fit - (nsapff_m5l_d$se * 1.96)
nsapff_m5l_d$uci <- nsapff_m5l_d$fit + (nsapff_m5l_d$se * 1.96)


psapff_m5a_d <- predict(sap_ff_g5t, newdata = nsapff_m5a_d, se.fit = T, type = 'response')
nsapff_m5a_d$fit <- psapff_m5a_d$fit
nsapff_m5a_d$se <- psapff_m5a_d$se.fit
nsapff_m5a_d$lci <- nsapff_m5a_d$fit - (nsapff_m5a_d$se * 1.96)
nsapff_m5a_d$uci <- nsapff_m5a_d$fit + (nsapff_m5a_d$se * 1.96)

psapff_m5h_d <- predict(sap_ff_g5t, newdata = nsapff_m5h_d, se.fit = T, type = 'response')
nsapff_m5h_d$fit <- psapff_m5h_d$fit
nsapff_m5h_d$se <- psapff_m5h_d$se.fit
nsapff_m5h_d$lci <- nsapff_m5h_d$fit - (nsapff_m5h_d$se * 1.96)
nsapff_m5h_d$uci <- nsapff_m5h_d$fit + (nsapff_m5h_d$se * 1.96)


# GAM bs = 'tp'
nsapff_m5l_t <- new_ff_m5l
nsapff_m5a_t <- new_ff_m5a
nsapff_m5h_t <- new_ff_m5h

psapff_m5l_t <- predict(sap_ff_g5.1t, newdata = nsapff_m5l_t, se.fit = T, type = 'response')
nsapff_m5l_t$fit <- psapff_m5l_t$fit
nsapff_m5l_t$se <- psapff_m5l_t$se.fit
nsapff_m5l_t$lci <- nsapff_m5l_t$fit - (nsapff_m5l_t$se * 1.96)
nsapff_m5l_t$uci <- nsapff_m5l_t$fit + (nsapff_m5l_t$se * 1.96)

psapff_m5a_t <- predict(sap_ff_g5.1t, newdata = nsapff_m5a_t, se.fit = T, type = 'response')
nsapff_m5a_t$fit <- psapff_m5a_t$fit
nsapff_m5a_t$se <- psapff_m5a_t$se.fit
nsapff_m5a_t$lci <- nsapff_m5a_t$fit - (nsapff_m5a_t$se * 1.96)
nsapff_m5a_t$uci <- nsapff_m5a_t$fit + (nsapff_m5a_t$se * 1.96)

psapff_m5h_t <- predict(sap_ff_g5.1t, newdata = nsapff_m5h_t, se.fit = T, type = 'response')
nsapff_m5h_t$fit <- psapff_m5h_t$fit
nsapff_m5h_t$se <- psapff_m5h_t$se.fit
nsapff_m5h_t$lci <- nsapff_m5h_t$fit - (nsapff_m5h_t$se * 1.96)
nsapff_m5h_t$uci <- nsapff_m5h_t$fit + (nsapff_m5h_t$se * 1.96)


# GAM bs = 'cr'
nsapff_m5l_c <- new_ff_m5l
nsapff_m5a_c <- new_ff_m5a
nsapff_m5h_c <- new_ff_m5h

psapff_m5l_c <- predict(sap_ff_g5.2t, newdata = nsapff_m5l_c, se.fit = T, type = 'response')
nsapff_m5l_c$fit <- psapff_m5l_c$fit
nsapff_m5l_c$se <- psapff_m5l_c$se.fit
nsapff_m5l_c$lci <- nsapff_m5l_c$fit - (nsapff_m5l_c$se * 1.96)
nsapff_m5l_c$uci <- nsapff_m5l_c$fit + (nsapff_m5l_c$se * 1.96)

psapff_m5a_c <- predict(sap_ff_g5.2t, newdata = nsapff_m5a_c, se.fit = T, type = 'response')
nsapff_m5a_c$fit <- psapff_m5a_c$fit
nsapff_m5a_c$se <- psapff_m5a_c$se.fit
nsapff_m5a_c$lci <- nsapff_m5a_c$fit - (nsapff_m5a_c$se * 1.96)
nsapff_m5a_c$uci <- nsapff_m5a_c$fit + (nsapff_m5a_c$se * 1.96)

psapff_m5h_c <- predict(sap_ff_g5.2t, newdata = nsapff_m5h_c, se.fit = T, type = 'response')
nsapff_m5h_c$fit <- psapff_m5h_c$fit
nsapff_m5h_c$se <- psapff_m5h_c$se.fit
nsapff_m5h_c$lci <- nsapff_m5h_c$fit - (nsapff_m5h_c$se * 1.96)
nsapff_m5h_c$uci <- nsapff_m5h_c$fit + (nsapff_m5h_c$se * 1.96)

# GAM bs = 'cc'
nsapff_m5l_cc <- new_ff_m5l
nsapff_m5a_cc <- new_ff_m5a
nsapff_m5h_cc <- new_ff_m5h

psapff_m5l_cc <- predict(sap_ff_g5.3t, newdata = nsapff_m5l_cc, se.fit = T, type = 'response')
nsapff_m5l_cc$fit <- psapff_m5l_cc$fit
nsapff_m5l_cc$se <- psapff_m5l_cc$se.fit
nsapff_m5l_cc$lci <- nsapff_m5l_cc$fit - (nsapff_m5l_cc$se * 1.96)
nsapff_m5l_cc$uci <- nsapff_m5l_cc$fit + (nsapff_m5l_cc$se * 1.96)

psapff_m5a_cc <- predict(sap_ff_g5.3t, newdata = nsapff_m5a_cc, se.fit = T, type = 'response')
nsapff_m5a_cc$fit <- psapff_m5a_cc$fit
nsapff_m5a_cc$se <- psapff_m5a_cc$se.fit
nsapff_m5a_cc$lci <- nsapff_m5a_cc$fit - (nsapff_m5a_cc$se * 1.96)
nsapff_m5a_cc$uci <- nsapff_m5a_cc$fit + (nsapff_m5a_cc$se * 1.96)

psapff_m5h_cc <- predict(sap_ff_g5.3t, newdata = nsapff_m5h_cc, se.fit = T, type = 'response')
nsapff_m5h_cc$fit <- psapff_m5h_cc$fit
nsapff_m5h_cc$se <- psapff_m5h_cc$se.fit
nsapff_m5h_cc$lci <- nsapff_m5h_cc$fit - (nsapff_m5h_cc$se * 1.96)
nsapff_m5h_cc$uci <- nsapff_m5h_cc$fit + (nsapff_m5h_cc$se * 1.96)

# Plot predictions proportion saplings -----
# GLMER and GAM comparisons
dev.new(width = 16, height = 12, noRStudioGD = T, dpi = 300)
par(mfrow = c(2,2), mar = c(6,6,2,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))

plot(nsapff_m5l_lm$r_Temp, nsapff_m5l_lm$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(-1.7, 1.5), xaxt = "n")
axis(side = 1, at = seq(-1.7, 1.5, 0.4), labels = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(sap_ff_m5t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'r_Temp', data = "nsapff_m5l_lm", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsapff_m5a_lm$r_Temp, nsapff_m5a_lm$fit, type = "l", col = "black")
pg.ci(x = 'r_Temp', data = "nsapff_m5a_lm", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsapff_m5h_lm$r_Temp, nsapff_m5h_lm$fit, type = "l", col = "red")
pg.ci(x = 'r_Temp', data = "nsapff_m5h_lm", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) GLMER")), cex = 2)



par(xpd = NA)
legend(x = 6.3, y = 1, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", "black", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)




plot(nsapff_m5l_d$Temp, nsapff_m5l_d$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(sap_ff_g5t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nsapff_m5l_d", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsapff_m5a_d$Temp, nsapff_m5a_d$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nsapff_m5a_d", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsapff_m5h_d$Temp, nsapff_m5h_d$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nsapff_m5h_d", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b) GAM defaults")), cex = 2)


plot(nsapff_m5l_t$Temp, nsapff_m5l_t$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4,  col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(sap_ff_g5.1t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nsapff_m5l_t", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsapff_m5a_t$Temp, nsapff_m5a_t$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nsapff_m5a_t", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsapff_m5h_t$Temp, nsapff_m5h_t$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nsapff_m5h_t", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) GAM bs = 'tp'")), cex = 2)


plot(nsapff_m5l_c$Temp, nsapff_m5l_c$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(sap_ff_g5.2t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nsapff_m5l_c", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsapff_m5a_c$Temp, nsapff_m5a_c$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nsapff_m5a_c", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsapff_m5h_c$Temp, nsapff_m5h_c$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nsapff_m5h_c", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(d) GAM bs = 'cr'")), cex = 2)



plot(nsapff_m5l_cc$Temp, nsapff_m5l_cc$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(sap_ff_g5.2t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nsapff_m5l_cc", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nsapff_m5a_cc$Temp, nsapff_m5a_cc$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nsapff_m5a_cc", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nsapff_m5h_cc$Temp, nsapff_m5h_cc$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nsapff_m5h_cc", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(e) GAM bs = 'cc'")), cex = 2)

# Same issue as with TSF, reduces to flat lines. 

# 6.1.3 Proportion recruits ----
recff_tnull <- glmer(Proportion_recruits ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
recff_gt_null <- gam(Proportion_recruits ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')

recff_lnull <- glmer(Proportion_recruits ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)
recff_gl_null <- gam(Proportion_recruits ~ 1, random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')



recff_m1t <- glmer(Proportion_recruits ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(recff_m1t)
recff_g1t <- gam(Proportion_recruits ~ s(Fire_freq, k = 8), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g1t)
plot(recff_g1t)
recff_g1.1t <- gam(Proportion_recruits ~ s(Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g1.1t)
plot(recff_g1.1t)

recff_m1l <- glmer(Proportion_recruits ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(recff_m1l) # Not a very good model
recff_g1l <- gam(Proportion_recruits ~ s(Fire_freq, k = 3), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g1l)
plot(recff_g1l)
recff_g1.1l <- gam(Proportion_recruits ~ s(Fire_freq, bs = 'cr', k = 3), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g1.1l)
plot(recff_g1.1l)


recff_m2t <- glmer(Proportion_recruits ~ r_Fire_freq * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recff_m2t)
recff_g2t <- gam(Proportion_recruits ~ s(Fire_freq, k = 8) + s(Latitude) + ti(Latitude, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g2t)
plot(recff_g2t)
recff_g2.1t <- gam(Proportion_recruits ~ s(Fire_freq, k = 8) + s(Latitude) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g2.1t)
plot(recff_g2.1t)
recff_g2.2t <- gam(Proportion_recruits ~ s(Fire_freq, bs = 'cr', k = 5) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g2.1t)
plot(recff_g2.2t)


recff_m2l <- glmer(Proportion_recruits ~ Fire_freq * Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30))) 
summary(recff_m2l)
recff_g2l <- gam(Proportion_recruits ~ s(Fire_freq, k = 3) + s(Latitude, k = 8) + ti(Latitude, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g2l)
plot(recff_g2l)
recff_g2.1l <- gam(Proportion_recruits ~ s(Fire_freq, k = 3) + s(Latitude, k = 8) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 8), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g2.1l)
plot(recff_g2.1l)
recff_g2.2l <- gam(Proportion_recruits ~ s(Fire_freq, bs = 'cr', k = 3) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g2.1l)
plot(recff_g2.2l)



recff_m3t <- glmer(Proportion_recruits ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
recff_g3t <- gam(Proportion_recruits ~ s(Fire_freq, k = 8) + s(FPC) + ti(FPC, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g3t)
plot(recff_g3t)
recff_g3.1t <- gam(Proportion_recruits ~ s(Fire_freq, k= 8) + s(FPC) + ti(FPC, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g3.1t)
plot(recff_g3.1t)
recff_g3.2t <- gam(Proportion_recruits ~ s(Fire_freq, bs = 'cr', k = 5) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g3.1t)
plot(recff_g3.2t)


recff_m3l <- glmer(Proportion_recruits ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(recff_m3l) # Not a very good model
recff_g3l <- gam(Proportion_recruits ~ s(Fire_freq, k = 3) + s(FPC, k = 4) + ti(FPC, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g3l)
plot(recff_g3l)
recff_g3.1l <- gam(Proportion_recruits ~ s(Fire_freq, k = 3) + s(FPC, k = 4) + ti(FPC, by = Fire_freq, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g3.1l)
plot(recff_g3.1l)
recff_g3.2l <- gam(Proportion_recruits ~ s(Fire_freq, bs = 'cr', k = 3) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g3.2l)
plot(recff_g3.2l)



recff_m4t <- glmer(Proportion_recruits ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects) # Warning, model fails to converge
recff_g4t <- gam(Proportion_recruits ~ s(Fire_freq, k = 8) + s(Precip) + ti(Precip, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g4t)
plot(recff_g4t)
recff_g4.1t <- gam(Proportion_recruits ~ s(Fire_freq, k= 8) + s(Precip) + ti(Precip, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g4.1t)
plot(recff_g4.1t)
recff_g4.2t <- gam(Proportion_recruits ~ s(Fire_freq, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g4.1t)
plot(recff_g4.2t)


recff_m4l <- glmer(Proportion_recruits ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=80)))
summary(recff_m4l) # Not a very good model
recff_g4l <- gam(Proportion_recruits ~ s(Fire_freq, k = 3) + s(Precip, k = 4) + ti(Precip, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g4l)
plot(recff_g4l)
recff_g4.1l <- gam(Proportion_recruits ~ s(Fire_freq, k = 3) + s(Precip, k = 4) + ti(Precip, by = Fire_freq, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g4.1l)
plot(recff_g4.1l)
recff_g4.2l <- gam(Proportion_recruits ~ s(Fire_freq, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g4.2l)
plot(recff_g4.2l)



recff_m5t <- glmer(Proportion_recruits ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
recff_g5t <- gam(Proportion_recruits ~ s(Fire_freq, k = 8) + s(Temp) + ti(Temp, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g5t)
plot(recff_g5t)
recff_g5.1t <- gam(Proportion_recruits ~ s(Fire_freq, k= 8) + s(Temp) + ti(Temp, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g5.1t)
plot(recff_g5.1t)
recff_g5.2t <- gam(Proportion_recruits ~ s(Fire_freq, bs = 'cr', k = 5) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = tor_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g5.1t)
plot(recff_g5.2t)

recff_m5l <- glmer(Proportion_recruits ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects)
summary(recff_m5l) # Not a very good model
recff_g5l <- gam(Proportion_recruits ~ s(Fire_freq, k = 3) + s(Temp, k = 4) + ti(Temp, by = Fire_freq), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g5l)
plot(recff_g5l)
recff_g5.1l <- gam(Proportion_recruits ~ s(Fire_freq, k = 3) + s(Temp, k = 4) + ti(Temp, by = Fire_freq, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g5.1l)
plot(recff_g5.1l)
recff_g5.2l <- gam(Proportion_recruits ~ s(Fire_freq, bs = 'cr', k = 3) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = binomial, data = lit_transects, method = 'ML')
par(mfrow = c(2,2)); gam.check(recff_g5.2l)
plot(recff_g5.2l)

# Compare AICs for each modelling type
# GLMER
recff_aic_t <- list(recff_tnull, recff_m1t, recff_m2t, recff_m3t, recff_m4t, recff_m5t)
aictab(recff_aic_t) # Model 5 is best

recff_aic_l <- list(recff_lnull, recff_m1l, recff_m2.1l, recff_m3l, recff_m4l, recff_m5l)
aictab(recff_aic_l) # Null model is best

# GAM defaults
recff_td_aic <- as.data.frame(1:6)
recff_td_aic$AICc <- "NA"
recff_td_aic$Model <- "NA"
recff_td_aic$LL <- "NA"
recff_td_aic$AICc[1] <- AICc(recff_gt_null)
recff_td_aic$Model[1] <- "Null"
recff_td_aic$LL[1] <- logLik(recff_gt_null)
recff_td_aic$AICc[2] <- AICc(recff_g1t)
recff_td_aic$Model[2] <- "m1"
recff_td_aic$LL[2] <- logLik(recff_g1t)
recff_td_aic$AICc[3] <- AICc(recff_g2t)
recff_td_aic$Model[3] <- "m2"
recff_td_aic$LL[3] <- logLik(recff_g2t)
recff_td_aic$AICc[4] <- AICc(recff_g3t)
recff_td_aic$Model[4] <- "m3"
recff_td_aic$LL[4] <- logLik(recff_g3t)
recff_td_aic$AICc[5] <- AICc(recff_g4t)
recff_td_aic$Model[5] <- "m4"
recff_td_aic$LL[5] <- logLik(recff_g4t)
recff_td_aic$AICc[6] <- AICc(recff_g5t)
recff_td_aic$Model[6] <- "m5"
recff_td_aic$LL[6] <- logLik(recff_g5t)
recff_td_aic$AICc <- as.numeric(recff_td_aic$AICc)
recff_td_aic$LL <- as.numeric(recff_td_aic$LL)
recff_td_aic <- recff_td_aic[order(recff_td_aic$AICc), ]
recff_td_aic$Delta_AICc <- "0.00"
recff_td_aic$Delta_AICc[2] <- round(recff_td_aic$AICc[1]-recff_td_aic$AICc[2], 2)
recff_td_aic$Delta_AICc[3] <- round(recff_td_aic$AICc[1]-recff_td_aic$AICc[3], 2)
recff_td_aic$Delta_AICc[4] <- round(recff_td_aic$AICc[1]-recff_td_aic$AICc[4], 2)
recff_td_aic$Delta_AICc[5] <- round(recff_td_aic$AICc[1]-recff_td_aic$AICc[5], 2)
recff_td_aic$Delta_AICc[6] <- round(recff_td_aic$AICc[1]-recff_td_aic$AICc[6], 2)
recff_td_aic  # Model 5 is best

recff_ld_aic <- as.data.frame(1:6)
recff_ld_aic$AICc <- "NA"
recff_ld_aic$Model <- "NA"
recff_ld_aic$LL <- "NA"
recff_ld_aic$AICc[1] <- AICc(recff_gl_null)
recff_ld_aic$Model[1] <- "Null"
recff_ld_aic$LL[1] <- logLik(recff_gl_null)
recff_ld_aic$AICc[2] <- AICc(recff_g1l)
recff_ld_aic$Model[2] <- "m1"
recff_ld_aic$LL[2] <- logLik(recff_g1l)
recff_ld_aic$AICc[3] <- AICc(recff_g2l)
recff_ld_aic$Model[3] <- "m2"
recff_ld_aic$LL[3] <- logLik(recff_g2l)
recff_ld_aic$AICc[4] <- AICc(recff_g3l)
recff_ld_aic$Model[4] <- "m3"
recff_ld_aic$LL[4] <- logLik(recff_g3l)
recff_ld_aic$AICc[5] <- AICc(recff_g4l)
recff_ld_aic$Model[5] <- "m4"
recff_ld_aic$LL[5] <- logLik(recff_g4l)
recff_ld_aic$AICc[6] <- AICc(recff_g5l)
recff_ld_aic$Model[6] <- "m5"
recff_ld_aic$LL[6] <- logLik(recff_g5l)
recff_ld_aic$AICc <- as.numeric(recff_ld_aic$AICc)
recff_ld_aic$LL <- as.numeric(recff_ld_aic$LL)
recff_ld_aic <- recff_ld_aic[order(recff_ld_aic$AICc), ]
recff_ld_aic$Delta_AICc <- "0.00"
recff_ld_aic$Delta_AICc[2] <- round(recff_ld_aic$AICc[1]-recff_ld_aic$AICc[2], 2)
recff_ld_aic$Delta_AICc[3] <- round(recff_ld_aic$AICc[1]-recff_ld_aic$AICc[3], 2)
recff_ld_aic$Delta_AICc[4] <- round(recff_ld_aic$AICc[1]-recff_ld_aic$AICc[4], 2)
recff_ld_aic$Delta_AICc[5] <- round(recff_ld_aic$AICc[1]-recff_ld_aic$AICc[5], 2)
recff_ld_aic$Delta_AICc[6] <- round(recff_ld_aic$AICc[1]-recff_ld_aic$AICc[6], 2)
recff_ld_aic  #Null model is best but model 1 is ranked within AICc <2.



# GAM s(bs = 'tp')
recff_ts_aic <- as.data.frame(1:6)
recff_ts_aic$AICc <- "NA"
recff_ts_aic$Model <- "NA"
recff_ts_aic$LL <- "NA"
recff_ts_aic$AICc[1] <- AICc(recff_gt_null)
recff_ts_aic$Model[1] <- "Null"
recff_ts_aic$LL[1] <- logLik(recff_gt_null)
recff_ts_aic$AICc[2] <- AICc(recff_g1t)
recff_ts_aic$Model[2] <- "m1"
recff_ts_aic$LL[2] <- logLik(recff_g1t)
recff_ts_aic$AICc[3] <- AICc(recff_g2.1t)
recff_ts_aic$Model[3] <- "m2"
recff_ts_aic$LL[3] <- logLik(recff_g2.1t)
recff_ts_aic$AICc[4] <- AICc(recff_g3.1t)
recff_ts_aic$Model[4] <- "m3"
recff_ts_aic$LL[4] <- logLik(recff_g3.1t)
recff_ts_aic$AICc[5] <- AICc(recff_g4.1t)
recff_ts_aic$Model[5] <- "m4"
recff_ts_aic$LL[5] <- logLik(recff_g4.1t)
recff_ts_aic$AICc[6] <- AICc(recff_g5.1t)
recff_ts_aic$Model[6] <- "m5"
recff_ts_aic$LL[6] <- logLik(recff_g5.1t)
recff_ts_aic$AICc <- as.numeric(recff_ts_aic$AICc)
recff_ts_aic$LL <- as.numeric(recff_ts_aic$LL)
recff_ts_aic <- recff_ts_aic[order(recff_ts_aic$AICc), ]
recff_ts_aic$Delta_AICc <- "0.00"
recff_ts_aic$Delta_AICc[2] <- round(recff_ts_aic$AICc[1]-recff_ts_aic$AICc[2], 2)
recff_ts_aic$Delta_AICc[3] <- round(recff_ts_aic$AICc[1]-recff_ts_aic$AICc[3], 2)
recff_ts_aic$Delta_AICc[4] <- round(recff_ts_aic$AICc[1]-recff_ts_aic$AICc[4], 2)
recff_ts_aic$Delta_AICc[5] <- round(recff_ts_aic$AICc[1]-recff_ts_aic$AICc[5], 2)
recff_ts_aic$Delta_AICc[6] <- round(recff_ts_aic$AICc[1]-recff_ts_aic$AICc[6], 2)
recff_ts_aic  # Model 5 is best.

recff_ls_aic <- as.data.frame(1:6)
recff_ls_aic$AICc <- "NA"
recff_ls_aic$Model <- "NA"
recff_ls_aic$LL <- "NA"
recff_ls_aic$AICc[1] <- AICc(recff_gl_null)
recff_ls_aic$Model[1] <- "Null"
recff_ls_aic$LL[1] <- logLik(recff_gl_null)
recff_ls_aic$AICc[2] <- AICc(recff_g1l)
recff_ls_aic$Model[2] <- "m1"
recff_ls_aic$LL[2] <- logLik(recff_g1l)
recff_ls_aic$AICc[3] <- AICc(recff_g2.1l)
recff_ls_aic$Model[3] <- "m2"
recff_ls_aic$LL[3] <- logLik(recff_g2.1l)
recff_ls_aic$AICc[4] <- AICc(recff_g3.1l)
recff_ls_aic$Model[4] <- "m3"
recff_ls_aic$LL[4] <- logLik(recff_g3.1l)
recff_ls_aic$AICc[5] <- AICc(recff_g4.1l)
recff_ls_aic$Model[5] <- "m4"
recff_ls_aic$LL[5] <- logLik(recff_g4.1l)
recff_ls_aic$AICc[6] <- AICc(recff_g5.1l)
recff_ls_aic$Model[6] <- "m5"
recff_ls_aic$LL[6] <- logLik(recff_g5.1l)
recff_ls_aic$AICc <- as.numeric(recff_ls_aic$AICc)
recff_ls_aic$LL <- as.numeric(recff_ls_aic$LL)
recff_ls_aic <- recff_ls_aic[order(recff_ls_aic$AICc), ]
recff_ls_aic$Delta_AICc <- "0.00"
recff_ls_aic$Delta_AICc[2] <- round(recff_ls_aic$AICc[1]-recff_ls_aic$AICc[2], 2)
recff_ls_aic$Delta_AICc[3] <- round(recff_ls_aic$AICc[1]-recff_ls_aic$AICc[3], 2)
recff_ls_aic$Delta_AICc[4] <- round(recff_ls_aic$AICc[1]-recff_ls_aic$AICc[4], 2)
recff_ls_aic$Delta_AICc[5] <- round(recff_ls_aic$AICc[1]-recff_ls_aic$AICc[5], 2)
recff_ls_aic$Delta_AICc[6] <- round(recff_ls_aic$AICc[1]-recff_ls_aic$AICc[6], 2)
recff_ls_aic  #Null model is best but model 1 is ranked within AICc <2.


# GAM ti(bs='cr')
recff_tsc_aic <- as.data.frame(1:6)
recff_tsc_aic$AICc <- "NA"
recff_tsc_aic$Model <- "NA"
recff_tsc_aic$LL <- "NA"
recff_tsc_aic$AICc[1] <- AICc(recff_gt_null)
recff_tsc_aic$Model[1] <- "Null"
recff_tsc_aic$LL[1] <- logLik(recff_gt_null)
recff_tsc_aic$AICc[2] <- AICc(recff_g1.1t)
recff_tsc_aic$Model[2] <- "m1"
recff_tsc_aic$LL[2] <- logLik(recff_g1.1t)
recff_tsc_aic$AICc[3] <- AICc(recff_g2.2t)
recff_tsc_aic$Model[3] <- "m2"
recff_tsc_aic$LL[3] <- logLik(recff_g2.2t)
recff_tsc_aic$AICc[4] <- AICc(recff_g3.2t)
recff_tsc_aic$Model[4] <- "m3"
recff_tsc_aic$LL[4] <- logLik(recff_g3.2t)
recff_tsc_aic$AICc[5] <- AICc(recff_g4.2t)
recff_tsc_aic$Model[5] <- "m4"
recff_tsc_aic$LL[5] <- logLik(recff_g4.2t)
recff_tsc_aic$AICc[6] <- AICc(recff_g5.2t)
recff_tsc_aic$Model[6] <- "m5"
recff_tsc_aic$LL[6] <- logLik(recff_g5.2t)
recff_tsc_aic$AICc <- as.numeric(recff_tsc_aic$AICc)
recff_tsc_aic$LL <- as.numeric(recff_tsc_aic$LL)
recff_tsc_aic <- recff_tsc_aic[order(recff_tsc_aic$AICc), ]
recff_tsc_aic$Delta_AICc <- "0.00"
recff_tsc_aic$Delta_AICc[2] <- round(recff_tsc_aic$AICc[1]-recff_tsc_aic$AICc[2], 2)
recff_tsc_aic$Delta_AICc[3] <- round(recff_tsc_aic$AICc[1]-recff_tsc_aic$AICc[3], 2)
recff_tsc_aic$Delta_AICc[4] <- round(recff_tsc_aic$AICc[1]-recff_tsc_aic$AICc[4], 2)
recff_tsc_aic$Delta_AICc[5] <- round(recff_tsc_aic$AICc[1]-recff_tsc_aic$AICc[5], 2)
recff_tsc_aic$Delta_AICc[6] <- round(recff_tsc_aic$AICc[1]-recff_tsc_aic$AICc[6], 2)
recff_tsc_aic  # Model 5 is best

recff_lsc_aic <- as.data.frame(1:6)
recff_lsc_aic$AICc <- "NA"
recff_lsc_aic$Model <- "NA"
recff_lsc_aic$LL <- "NA"
recff_lsc_aic$AICc[1] <- AICc(recff_gl_null)
recff_lsc_aic$Model[1] <- "Null"
recff_lsc_aic$LL[1] <- logLik(recff_gl_null)
recff_lsc_aic$AICc[2] <- AICc(recff_g1l)
recff_lsc_aic$Model[2] <- "m1"
recff_lsc_aic$LL[2] <- logLik(recff_g1l)
recff_lsc_aic$AICc[3] <- AICc(recff_g2.2l)
recff_lsc_aic$Model[3] <- "m2"
recff_lsc_aic$LL[3] <- logLik(recff_g2.2l)
recff_lsc_aic$AICc[4] <- AICc(recff_g3.2l)
recff_lsc_aic$Model[4] <- "m3"
recff_lsc_aic$LL[4] <- logLik(recff_g3.2l)
recff_lsc_aic$AICc[5] <- AICc(recff_g4.2l)
recff_lsc_aic$Model[5] <- "m4"
recff_lsc_aic$LL[5] <- logLik(recff_g4.2l)
recff_lsc_aic$AICc[6] <- AICc(recff_g5.2l)
recff_lsc_aic$Model[6] <- "m5"
recff_lsc_aic$LL[6] <- logLik(recff_g5.2l)
recff_lsc_aic$AICc <- as.numeric(recff_lsc_aic$AICc)
recff_lsc_aic$LL <- as.numeric(recff_lsc_aic$LL)
recff_lsc_aic <- recff_lsc_aic[order(recff_lsc_aic$AICc), ]
recff_lsc_aic$Delta_AICc <- "0.00"
recff_lsc_aic$Delta_AICc[2] <- round(recff_lsc_aic$AICc[1]-recff_lsc_aic$AICc[2], 2)
recff_lsc_aic$Delta_AICc[3] <- round(recff_lsc_aic$AICc[1]-recff_lsc_aic$AICc[3], 2)
recff_lsc_aic$Delta_AICc[4] <- round(recff_lsc_aic$AICc[1]-recff_lsc_aic$AICc[4], 2)
recff_lsc_aic$Delta_AICc[5] <- round(recff_lsc_aic$AICc[1]-recff_lsc_aic$AICc[5], 2)
recff_lsc_aic$Delta_AICc[6] <- round(recff_lsc_aic$AICc[1]-recff_lsc_aic$AICc[6], 2)
recff_lsc_aic  #Null model is best but model 1 is ranked within AICc <2. 

# Predict for proportion recruits data ----
# GLMER
nrecruitff_m5l_lm <- new_ff_m5l_lm
nrecruitff_m5a_lm <- new_ff_m5a_lm
nrecruitff_m5h_lm <- new_ff_m5h_lm

precruitff_m5l_lm <- predictSE(recff_m5t, newdata = nrecruitff_m5l_lm, se.fit = T, type = 'link')
nrecruitff_m5l_lm$fit.link <- precruitff_m5l_lm$fit
nrecruitff_m5l_lm$se.link <- precruitff_m5l_lm$se.fit
nrecruitff_m5l_lm$lci.link <- nrecruitff_m5l_lm$fit.link - (nrecruitff_m5l_lm$se.link * 1.96)
nrecruitff_m5l_lm$uci.link <- nrecruitff_m5l_lm$fit.link + (nrecruitff_m5l_lm$se.link * 1.96)
nrecruitff_m5l_lm$fit <- invlogit(nrecruitff_m5l_lm$fit.link)
nrecruitff_m5l_lm$se <- invlogit(nrecruitff_m5l_lm$se.link)
nrecruitff_m5l_lm$lci <- invlogit(nrecruitff_m5l_lm$lci.link)
nrecruitff_m5l_lm$uci <- invlogit(nrecruitff_m5l_lm$uci.link)

precruitff_m5a_lm <- predictSE(recff_m5t, newdata = nrecruitff_m5a_lm, se.fit = T, type = 'link')
nrecruitff_m5a_lm$fit.link <- precruitff_m5a_lm$fit
nrecruitff_m5a_lm$se.link <- precruitff_m5a_lm$se.fit
nrecruitff_m5a_lm$lci.link <- nrecruitff_m5a_lm$fit.link - (nrecruitff_m5a_lm$se.link * 1.96)
nrecruitff_m5a_lm$uci.link <- nrecruitff_m5a_lm$fit.link + (nrecruitff_m5a_lm$se.link * 1.96)
nrecruitff_m5a_lm$fit <- invlogit(nrecruitff_m5a_lm$fit.link)
nrecruitff_m5a_lm$se <- invlogit(nrecruitff_m5a_lm$se.link)
nrecruitff_m5a_lm$lci <- invlogit(nrecruitff_m5a_lm$lci.link)
nrecruitff_m5a_lm$uci <- invlogit(nrecruitff_m5a_lm$uci.link)

precruitff_m5h_lm <- predictSE(recff_m5t, newdata = nrecruitff_m5h_lm, se.fit = T, type = 'link')
nrecruitff_m5h_lm$fit.link <- precruitff_m5h_lm$fit
nrecruitff_m5h_lm$se.link <- precruitff_m5h_lm$se.fit
nrecruitff_m5h_lm$lci.link <- nrecruitff_m5h_lm$fit.link - (nrecruitff_m5h_lm$se.link * 1.96)
nrecruitff_m5h_lm$uci.link <- nrecruitff_m5h_lm$fit.link + (nrecruitff_m5h_lm$se.link * 1.96)
nrecruitff_m5h_lm$fit <- invlogit(nrecruitff_m5h_lm$fit.link)
nrecruitff_m5h_lm$se <- invlogit(nrecruitff_m5h_lm$se.link)
nrecruitff_m5h_lm$lci <- invlogit(nrecruitff_m5h_lm$lci.link)
nrecruitff_m5h_lm$uci <- invlogit(nrecruitff_m5h_lm$uci.link)




# GAM default
nrecruitff_m5l_d <- new_ff_m5l
nrecruitff_m5a_d <- new_ff_m5a
nrecruitff_m5h_d <- new_ff_m5h

precruitff_m5l_d <- predict(recff_g5t, newdata = nrecruitff_m5l_d, se.fit = T, type = 'response')
nrecruitff_m5l_d$fit <- precruitff_m5l_d$fit
nrecruitff_m5l_d$se <- precruitff_m5l_d$se.fit
nrecruitff_m5l_d$lci <- nrecruitff_m5l_d$fit - (nrecruitff_m5l_d$se * 1.96)
nrecruitff_m5l_d$uci <- nrecruitff_m5l_d$fit + (nrecruitff_m5l_d$se * 1.96)


precruitff_m5a_d <- predict(recff_g5t, newdata = nrecruitff_m5a_d, se.fit = T, type = 'response')
nrecruitff_m5a_d$fit <- precruitff_m5a_d$fit
nrecruitff_m5a_d$se <- precruitff_m5a_d$se.fit
nrecruitff_m5a_d$lci <- nrecruitff_m5a_d$fit - (nrecruitff_m5a_d$se * 1.96)
nrecruitff_m5a_d$uci <- nrecruitff_m5a_d$fit + (nrecruitff_m5a_d$se * 1.96)

precruitff_m5h_d <- predict(recff_g5t, newdata = nrecruitff_m5h_d, se.fit = T, type = 'response')
nrecruitff_m5h_d$fit <- precruitff_m5h_d$fit
nrecruitff_m5h_d$se <- precruitff_m5h_d$se.fit
nrecruitff_m5h_d$lci <- nrecruitff_m5h_d$fit - (nrecruitff_m5h_d$se * 1.96)
nrecruitff_m5h_d$uci <- nrecruitff_m5h_d$fit + (nrecruitff_m5h_d$se * 1.96)


# GAM bs = 'tp'
nrecruitff_m5l_t <- new_ff_m5l
nrecruitff_m5a_t <- new_ff_m5a
nrecruitff_m5h_t <- new_ff_m5h

precruitff_m5l_t <- predict(recff_g5.1t, newdata = nrecruitff_m5l_t, se.fit = T, type = 'response')
nrecruitff_m5l_t$fit <- precruitff_m5l_t$fit
nrecruitff_m5l_t$se <- precruitff_m5l_t$se.fit
nrecruitff_m5l_t$lci <- nrecruitff_m5l_t$fit - (nrecruitff_m5l_t$se * 1.96)
nrecruitff_m5l_t$uci <- nrecruitff_m5l_t$fit + (nrecruitff_m5l_t$se * 1.96)

precruitff_m5a_t <- predict(recff_g5.1t, newdata = nrecruitff_m5a_t, se.fit = T, type = 'response')
nrecruitff_m5a_t$fit <- precruitff_m5a_t$fit
nrecruitff_m5a_t$se <- precruitff_m5a_t$se.fit
nrecruitff_m5a_t$lci <- nrecruitff_m5a_t$fit - (nrecruitff_m5a_t$se * 1.96)
nrecruitff_m5a_t$uci <- nrecruitff_m5a_t$fit + (nrecruitff_m5a_t$se * 1.96)

precruitff_m5h_t <- predict(recff_g5.1t, newdata = nrecruitff_m5h_t, se.fit = T, type = 'response')
nrecruitff_m5h_t$fit <- precruitff_m5h_t$fit
nrecruitff_m5h_t$se <- precruitff_m5h_t$se.fit
nrecruitff_m5h_t$lci <- nrecruitff_m5h_t$fit - (nrecruitff_m5h_t$se * 1.96)
nrecruitff_m5h_t$uci <- nrecruitff_m5h_t$fit + (nrecruitff_m5h_t$se * 1.96)


# GAM bs = 'cr'
nrecruitff_m5l_c <- new_ff_m5l
nrecruitff_m5a_c <- new_ff_m5a
nrecruitff_m5h_c <- new_ff_m5h

precruitff_m5l_c <- predict(recff_g5.2t, newdata = nrecruitff_m5l_c, se.fit = T, type = 'response')
nrecruitff_m5l_c$fit <- precruitff_m5l_c$fit
nrecruitff_m5l_c$se <- precruitff_m5l_c$se.fit
nrecruitff_m5l_c$lci <- nrecruitff_m5l_c$fit - (nrecruitff_m5l_c$se * 1.96)
nrecruitff_m5l_c$uci <- nrecruitff_m5l_c$fit + (nrecruitff_m5l_c$se * 1.96)

precruitff_m5a_c <- predict(recff_g5.2t, newdata = nrecruitff_m5a_c, se.fit = T, type = 'response')
nrecruitff_m5a_c$fit <- precruitff_m5a_c$fit
nrecruitff_m5a_c$se <- precruitff_m5a_c$se.fit
nrecruitff_m5a_c$lci <- nrecruitff_m5a_c$fit - (nrecruitff_m5a_c$se * 1.96)
nrecruitff_m5a_c$uci <- nrecruitff_m5a_c$fit + (nrecruitff_m5a_c$se * 1.96)

precruitff_m5h_c <- predict(recff_g5.2t, newdata = nrecruitff_m5h_c, se.fit = T, type = 'response')
nrecruitff_m5h_c$fit <- precruitff_m5h_c$fit
nrecruitff_m5h_c$se <- precruitff_m5h_c$se.fit
nrecruitff_m5h_c$lci <- nrecruitff_m5h_c$fit - (nrecruitff_m5h_c$se * 1.96)
nrecruitff_m5h_c$uci <- nrecruitff_m5h_c$fit + (nrecruitff_m5h_c$se * 1.96)



# Plot predictions proportion recruits -----
# GLMER and GAM comparisons
dev.new(width = 16, height = 12, noRStudioGD = T, dpi = 300)
par(mfrow = c(2,2), mar = c(6,6,2,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))

plot(nrecruitff_m5l_lm$r_Temp, nrecruitff_m5l_lm$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(-1.7, 1.5), xaxt = "n")
axis(side = 1, at = seq(-1.7, 1.5, 0.4), labels = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(recff_m5t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'r_Temp', data = "nrecruitff_m5l_lm", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nrecruitff_m5a_lm$r_Temp, nrecruitff_m5a_lm$fit, type = "l", col = "black")
pg.ci(x = 'r_Temp', data = "nrecruitff_m5a_lm", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nrecruitff_m5h_lm$r_Temp, nrecruitff_m5h_lm$fit, type = "l", col = "red")
pg.ci(x = 'r_Temp', data = "nrecruitff_m5h_lm", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a) GLMER")), cex = 2)



par(xpd = NA)
legend(x = 6.3, y = 1, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", "black", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)




plot(nrecruitff_m5l_d$Temp, nrecruitff_m5l_d$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(recff_g5t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nrecruitff_m5l_d", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nrecruitff_m5a_d$Temp, nrecruitff_m5a_d$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nrecruitff_m5a_d", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nrecruitff_m5h_d$Temp, nrecruitff_m5h_d$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nrecruitff_m5h_d", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b) GAM defaults")), cex = 2)


plot(nrecruitff_m5l_t$Temp, nrecruitff_m5l_t$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4,  col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(recff_g5.1t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nrecruitff_m5l_t", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nrecruitff_m5a_t$Temp, nrecruitff_m5a_t$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nrecruitff_m5a_t", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nrecruitff_m5h_t$Temp, nrecruitff_m5h_t$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nrecruitff_m5h_t", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(c) GAM bs = 'tp'")), cex = 2)


plot(nrecruitff_m5l_c$Temp, nrecruitff_m5l_c$fit, type = 'l', las = 1, ylim = c(0,1), xlab = "", ylab = "", cex.axis = 1.4, col = 'blue', xlim = c(385, 425), xaxt = "n")
axis(side = 1, at = seq(385, 425, 5), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 4, cex = 1.5)
mtext(side = 2, expression(bold("Proportion saplings")), line = 3, cex = 1.5)
mtext(paste("AICc = ", round(AICc(recff_g5.2t ),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Temp', data = "nrecruitff_m5l_c", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(nrecruitff_m5a_c$Temp, nrecruitff_m5a_c$fit, type = "l", col = "black")
pg.ci(x = 'Temp', data = "nrecruitff_m5a_c", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(nrecruitff_m5h_c$Temp, nrecruitff_m5h_c$fit, type = "l", col = "red")
pg.ci(x = 'Temp', data = "nrecruitff_m5h_c", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(d) GAM bs = 'cr'")), cex = 2)






# 6.2 Female fecundity  ----
# FF 
# FF + height 
# FF * latitude
# FF * FPC
# FF * precipitation seasonality
# FF * temperature seasonality

tor_fecundity$r_Fire_freq <- scale(tor_fecundity$Fire_freq)
lit_fecundity$r_Fire_freq <- scale(lit_fecundity$Fire_freq)

fecundityff_tnull <- glmer(Cone_number ~ 1 + (1 | Location/Transect), family = poisson, data = tor_fecundity)
fecundityff_gt_null <- gam(Cone_number ~ 1, random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')

fecundityff_lnull <- glmer(Cone_number ~ 1 + (1 | Location/Transect), family = poisson, data = lit_fecundity)
fecundityff_gl_null <- gam(Cone_number ~ 1, random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')



fecundityff_m1t <- glmer(Cone_number ~ r_Fire_freq + (1 | Location/Transect), family = poisson, data = tor_fecundity)
summary(fecundityff_m1t)
fecundityff_g1t <- gam(Cone_number ~ s(Fire_freq, k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g1t)
plot(fecundityff_g1t)
fecundityff_g1.1t <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g1.1t)
plot(fecundityff_g1.1t)
fecundityff_g1ot <- gam(Cone_number ~ s(Fire_freq, k = 4), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g1ot) # k' and edf are always close but this gives us a nice smooth curve
plot(fecundityff_g1ot)
fecundityff_g1.1ot <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 4), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g1.1ot)
plot(fecundityff_g1.1ot)


fecundityff_m1l <- glmer(Cone_number ~ r_Fire_freq + (1 | Location/Transect), family = poisson, data = lit_fecundity)
summary(fecundityff_m1l) 
fecundityff_g1l <- gam(Cone_number ~ s(Fire_freq, k = 3), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g1l)
plot(fecundityff_g1l)
fecundityff_g1.1l <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 3), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g1.1l)
plot(fecundityff_g1.1l)




fecundityff_m2t <- glmer(Cone_number ~ r_Fire_freq * r_Latitude + (1|Location/Transect), family = poisson, data = tor_fecundity)
summary(fecundityff_m2t)
fecundityff_g2t <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(Latitude) + ti(Latitude, by = Fire_freq), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g2t)
plot(fecundityff_g2t)
fecundityff_g2.1t <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(Latitude) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g2.1t)
plot(fecundityff_g2.1t)
fecundityff_g2.2t <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 5) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g2.2t)
plot(fecundityff_g2.2t)
fecundityff_g2.1ot <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(Latitude, k = 5) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g2.1ot)
plot(fecundityff_g2.1ot)
fecundityff_g2.2ot <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 4) + s(Latitude, bs = 'cr', k = 4) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g2.2ot)
plot(fecundityff_g2.2ot)

fecundityff_m2l <- glmer(Cone_number ~ r_Fire_freq * Latitude + (1|Location/Transect), family = poisson, data = lit_fecundity) 
summary(fecundityff_m2l)
fecundityff_g2l <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(Latitude) + ti(Latitude, by = Fire_freq), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g2l)
plot(fecundityff_g2l)
fecundityff_g2.1l <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(Latitude) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g2.1l)
plot(fecundityff_g2.1l)
fecundityff_g2.2l <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 3) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g2.2l)
plot(fecundityff_g2.2l)



fecundityff_m3t <- glmer(Cone_number ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = poisson, data = tor_fecundity)
fecundityff_g3t <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(FPC) + ti(FPC, by = Fire_freq), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g3t)
plot(fecundityff_g3t)
fecundityff_g3.1t <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(FPC) + ti(FPC, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g3.1t)
plot(fecundityff_g3.1t)
fecundityff_g3.2t <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 5) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g3.2t)
plot(fecundityff_g3.2t)
fecundityff_g3.1ot <- gam(Cone_number ~ s(Fire_freq, k = 4) + s(FPC, k = 4) + ti(FPC, by = Fire_freq, bs = 'tp', k = 13), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g3.1ot)
plot(fecundityff_g3.1ot)
fecundityff_g3.2ot <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 4) + s(FPC, bs = 'cr', k = 4) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g3.2ot)
plot(fecundityff_g3.2ot)


fecundityff_m3l <- glmer(Cone_number ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = poisson, data = lit_fecundity)
summary(fecundityff_m3l) 
fecundityff_g3l <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(FPC) + ti(FPC, by = Fire_freq), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g3l)
plot(fecundityff_g3l)
fecundityff_g3.1l <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(FPC) + ti(FPC, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g3.1l)
plot(fecundityff_g3.1l)
fecundityff_g3.2l <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 3) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g3.2l)
plot(fecundityff_g3.2l)
fecundityff_g3.1ol <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(FPC, k = 6) + ti(FPC, by = Fire_freq, bs = 'tp', k = 6), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g3.1ol)
plot(fecundityff_g3.1ol)
fecundityff_g3.2ol <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 3) + s(FPC, bs = 'cr', k = 4) + ti(FPC, by = Fire_freq, bs = 'cr', k = 4), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g3.2ol)
plot(fecundityff_g3.2ol)

fecundityff_m4t <- glmer(Cone_number ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = poisson, data = tor_fecundity)
fecundityff_g4t <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(Precip) + ti(Precip, by = Fire_freq), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g4t)
plot(fecundityff_g4t)
fecundityff_g4.1t <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(Precip) + ti(Precip, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g4.1t)
plot(fecundityff_g4.1t)
fecundityff_g4.2t <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g4.2t)
plot(fecundityff_g4.2t)
fecundityff_g4.1ot <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(Precip, k = 4) + ti(Precip, by = Fire_freq, bs = 'tp', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g4.1ot)
plot(fecundityff_g4.1ot)
fecundityff_g4.2ot <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 4) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g4.2ot)
plot(fecundityff_g4.2ot)


fecundityff_m4l <- glmer(Cone_number ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = poisson, data = lit_fecundity)
summary(fecundityff_m4l)
fecundityff_g4l <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(Precip, k = 5) + ti(Precip, by = Fire_freq), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g4l)
plot(fecundityff_g4l)
fecundityff_g4.1l <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(Precip, k = 5) + ti(Precip, by = Fire_freq, bs = 'tp', k = 6), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g4.1l)
plot(fecundityff_g4.1l)
fecundityff_g4.2l <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g4.2l)
plot(fecundityff_g4.2l)
fecundityff_g4.1ol <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(Precip, k = 3) + ti(Precip, by = Fire_freq, bs = 'tp', k = 3), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g4.1ol)
plot(fecundityff_g4.1ol)
fecundityff_g4.2ol <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g4.2ol)
plot(fecundityff_g4.2ol)


fecundityff_m5t <- glmer(Cone_number ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = poisson, data = tor_fecundity)
fecundityff_g5t <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(Precip) + ti(Precip, by = Fire_freq), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g5t)
plot(fecundityff_g5t)
fecundityff_g5.1t <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(Precip) + ti(Precip, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g5.1t)
plot(fecundityff_g5.1t)
fecundityff_g5.2t <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g5.2t)
plot(fecundityff_g5.2t)
fecundityff_g5.1ot <- gam(Cone_number ~ s(Fire_freq, k = 4) + s(Precip, k = 4) + ti(Precip, by = Fire_freq, bs = 'tp', k = 7), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g5.1ot)
plot(fecundityff_g5.1ot)
fecundityff_g5.2ot <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 3) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g5.2ot)
plot(fecundityff_g5.2ot)


fecundityff_m5l <- glmer(Cone_number ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = poisson, data = lit_fecundity)
summary(fecundityff_m5l)
fecundityff_g5l <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(Precip, k = 5) + ti(Precip, by = Fire_freq), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g5l)
plot(fecundityff_g5l)
fecundityff_g5.1l <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(Precip, k = 5) + ti(Precip, by = Fire_freq, bs = 'tp', k = 6), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g5.1l)
plot(fecundityff_g5.1l)
fecundityff_g5.2l <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g5.2l)
plot(fecundityff_g5.2l)
fecundityff_g5.1ol <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(Precip, k = 4) + ti(Precip, by = Fire_freq, bs = 'tp', k = 4), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g5.1ol)
plot(fecundityff_g5.1ol)
fecundityff_g5.2ol <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 3), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g5.2ol)
plot(fecundityff_g5.2ol)



fecundityff_m6t <- glmer(Cone_number ~ r_Fire_freq + r_height + (1|Location/Transect), family = poisson, data = tor_fecundity)
fecundityff_g6t <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(Height_cm), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g6t)
plot(fecundityff_g6t)
fecundityff_g6.1t <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(Height_cm) + ti(Height_cm, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g6.1t)
plot(fecundityff_g6.1t)
fecundityff_g6.2t <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 5) + s(Height_cm, bs = 'cr', k = 5) + ti(Height_cm, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g6.2t)
plot(fecundityff_g6.2t)
fecundityff_g6.1ot <- gam(Cone_number ~ s(Fire_freq, k = 5) + s(Height_cm, k = 5) + ti(Height_cm, by = Fire_freq, bs = 'tp', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g6.1ot)
plot(fecundityff_g6.1ot)
fecundityff_g6.2ot <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 4) + s(Height_cm, bs = 'cr', k = 4) + ti(Height_cm, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = tor_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g6.2ot)
plot(fecundityff_g6.2ot)


fecundityff_m6l <- glmer(Cone_number ~ r_Fire_freq + r_height + (1|Location/Transect), family = poisson, data = lit_fecundity)
summary(fecundityff_m6l)
fecundityff_g6l <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(Height_cm), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g6l)
plot(fecundityff_g6l)
fecundityff_g6.1l <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(Height_cm), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g6.1l)
plot(fecundityff_g6.1l)
fecundityff_g6.2l <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 3) + s(Height_cm, bs = 'cr', k = 5), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g6.2l)
plot(fecundityff_g6.2l)
fecundityff_g6.1ol <- gam(Cone_number ~ s(Fire_freq, k = 3) + s(Height_cm, k = 4), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g6.1ol)
plot(fecundityff_g6.1ol)
fecundityff_g6.2ol <- gam(Cone_number ~ s(Fire_freq, bs = 'cr', k = 3) + s(Height_cm, bs = 'cr', k = 3), random = ~ (1 | Transect/Location), family = poisson, data = lit_fecundity, method = 'ML')
par(mfrow = c(2,2)); gam.check(fecundityff_g6.2ol)
plot(fecundityff_g6.2ol)



# Compare AICs for each modelling type
# GLMER
fecundityff_aic_t <- list(fecundityff_tnull, fecundityff_m1t, fecundityff_m2t, fecundityff_m3t, fecundityff_m4t, fecundityff_m5t)
aictab(fecundityff_aic_t) # Model 3 is best with FF and FPC

fecundityff_aic_l <- list(fecundityff_lnull, fecundityff_m1l, fecundityff_m2l, fecundityff_m3l, fecundityff_m4l, fecundityff_m5l)
aictab(fecundityff_aic_l) # Model 3 is best with FF and FPC

# GAM defaults
fecundityff_td_aic <- as.data.frame(1:7)
fecundityff_td_aic$AICc <- "NA"
fecundityff_td_aic$Model <- "NA"
fecundityff_td_aic$LL <- "NA"
fecundityff_td_aic$AICc[1] <- AICc(fecundityff_gt_null)
fecundityff_td_aic$Model[1] <- "Null"
fecundityff_td_aic$LL[1] <- logLik(fecundityff_gt_null)
fecundityff_td_aic$AICc[2] <- AICc(fecundityff_g1t)
fecundityff_td_aic$Model[2] <- "m1"
fecundityff_td_aic$LL[2] <- logLik(fecundityff_g1t)
fecundityff_td_aic$AICc[3] <- AICc(fecundityff_g2t)
fecundityff_td_aic$Model[3] <- "m2"
fecundityff_td_aic$LL[3] <- logLik(fecundityff_g2t)
fecundityff_td_aic$AICc[4] <- AICc(fecundityff_g3t)
fecundityff_td_aic$Model[4] <- "m3"
fecundityff_td_aic$LL[4] <- logLik(fecundityff_g3t)
fecundityff_td_aic$AICc[5] <- AICc(fecundityff_g4t)
fecundityff_td_aic$Model[5] <- "m4"
fecundityff_td_aic$LL[5] <- logLik(fecundityff_g4t)
fecundityff_td_aic$AICc[6] <- AICc(fecundityff_g5t)
fecundityff_td_aic$Model[6] <- "m5"
fecundityff_td_aic$LL[6] <- logLik(fecundityff_g5t)
fecundityff_td_aic$AICc[7] <- AICc(fecundityff_g6t)
fecundityff_td_aic$Model[7] <- "m6"
fecundityff_td_aic$LL[7] <- logLik(fecundityff_g6t)
fecundityff_td_aic$AICc <- as.numeric(fecundityff_td_aic$AICc)
fecundityff_td_aic$LL <- as.numeric(fecundityff_td_aic$LL)
fecundityff_td_aic <- fecundityff_td_aic[order(fecundityff_td_aic$AICc), ]
fecundityff_td_aic$Delta_AICc <- "0.00"
fecundityff_td_aic$Delta_AICc[2] <- round(fecundityff_td_aic$AICc[1]-fecundityff_td_aic$AICc[2], 2)
fecundityff_td_aic$Delta_AICc[3] <- round(fecundityff_td_aic$AICc[1]-fecundityff_td_aic$AICc[3], 2)
fecundityff_td_aic$Delta_AICc[4] <- round(fecundityff_td_aic$AICc[1]-fecundityff_td_aic$AICc[4], 2)
fecundityff_td_aic$Delta_AICc[5] <- round(fecundityff_td_aic$AICc[1]-fecundityff_td_aic$AICc[5], 2)
fecundityff_td_aic$Delta_AICc[6] <- round(fecundityff_td_aic$AICc[1]-fecundityff_td_aic$AICc[6], 2)
fecundityff_td_aic$Delta_AICc[7] <- round(fecundityff_td_aic$AICc[1]-fecundityff_td_aic$AICc[7], 2)
fecundityff_td_aic  # Model 6 is best.

fecundityff_ld_aic <- as.data.frame(1:7)
fecundityff_ld_aic$AICc <- "NA"
fecundityff_ld_aic$Model <- "NA"
fecundityff_ld_aic$LL <- "NA"
fecundityff_ld_aic$AICc[1] <- AICc(fecundityff_gl_null)
fecundityff_ld_aic$Model[1] <- "Null"
fecundityff_ld_aic$LL[1] <- logLik(fecundityff_gl_null)
fecundityff_ld_aic$AICc[2] <- AICc(fecundityff_g1l)
fecundityff_ld_aic$Model[2] <- "m1"
fecundityff_ld_aic$LL[2] <- logLik(fecundityff_g1l)
fecundityff_ld_aic$AICc[3] <- AICc(fecundityff_g2l)
fecundityff_ld_aic$Model[3] <- "m2"
fecundityff_ld_aic$LL[3] <- logLik(fecundityff_g2l)
fecundityff_ld_aic$AICc[4] <- AICc(fecundityff_g3l)
fecundityff_ld_aic$Model[4] <- "m3"
fecundityff_ld_aic$LL[4] <- logLik(fecundityff_g3l)
fecundityff_ld_aic$AICc[5] <- AICc(fecundityff_g4l)
fecundityff_ld_aic$Model[5] <- "m4"
fecundityff_ld_aic$LL[5] <- logLik(fecundityff_g4l)
fecundityff_ld_aic$AICc[6] <- AICc(fecundityff_g5l)
fecundityff_ld_aic$Model[6] <- "m5"
fecundityff_ld_aic$LL[6] <- logLik(fecundityff_g5l)
fecundityff_ld_aic$AICc[7] <- AICc(fecundityff_g6l)
fecundityff_ld_aic$Model[7] <- "m6"
fecundityff_ld_aic$LL[7] <- logLik(fecundityff_g6l)
fecundityff_ld_aic$AICc <- as.numeric(fecundityff_ld_aic$AICc)
fecundityff_ld_aic$LL <- as.numeric(fecundityff_ld_aic$LL)
fecundityff_ld_aic <- fecundityff_ld_aic[order(fecundityff_ld_aic$AICc), ]
fecundityff_ld_aic$Delta_AICc <- "0.00"
fecundityff_ld_aic$Delta_AICc[2] <- round(fecundityff_ld_aic$AICc[1]-fecundityff_ld_aic$AICc[2], 2)
fecundityff_ld_aic$Delta_AICc[3] <- round(fecundityff_ld_aic$AICc[1]-fecundityff_ld_aic$AICc[3], 2)
fecundityff_ld_aic$Delta_AICc[4] <- round(fecundityff_ld_aic$AICc[1]-fecundityff_ld_aic$AICc[4], 2)
fecundityff_ld_aic$Delta_AICc[5] <- round(fecundityff_ld_aic$AICc[1]-fecundityff_ld_aic$AICc[5], 2)
fecundityff_ld_aic$Delta_AICc[6] <- round(fecundityff_ld_aic$AICc[1]-fecundityff_ld_aic$AICc[6], 2)
fecundityff_ld_aic$Delta_AICc[7] <- round(fecundityff_ld_aic$AICc[1]-fecundityff_ld_aic$AICc[7], 2)
fecundityff_ld_aic  #Model 6 is best



# GAM s(bs = 'tp')
fecundityff_ts_aic <- as.data.frame(1:7)
fecundityff_ts_aic$AICc <- "NA"
fecundityff_ts_aic$Model <- "NA"
fecundityff_ts_aic$LL <- "NA"
fecundityff_ts_aic$AICc[1] <- AICc(fecundityff_gt_null)
fecundityff_ts_aic$Model[1] <- "Null"
fecundityff_ts_aic$LL[1] <- logLik(fecundityff_gt_null)
fecundityff_ts_aic$AICc[2] <- AICc(fecundityff_g1t)
fecundityff_ts_aic$Model[2] <- "m1"
fecundityff_ts_aic$LL[2] <- logLik(fecundityff_g1t)
fecundityff_ts_aic$AICc[3] <- AICc(fecundityff_g2.1t)
fecundityff_ts_aic$Model[3] <- "m2"
fecundityff_ts_aic$LL[3] <- logLik(fecundityff_g2.1t)
fecundityff_ts_aic$AICc[4] <- AICc(fecundityff_g3.1t)
fecundityff_ts_aic$Model[4] <- "m3"
fecundityff_ts_aic$LL[4] <- logLik(fecundityff_g3.1t)
fecundityff_ts_aic$AICc[5] <- AICc(fecundityff_g4.1t)
fecundityff_ts_aic$Model[5] <- "m4"
fecundityff_ts_aic$LL[5] <- logLik(fecundityff_g4.1t)
fecundityff_ts_aic$AICc[6] <- AICc(fecundityff_g5.1t)
fecundityff_ts_aic$Model[6] <- "m5"
fecundityff_ts_aic$LL[6] <- logLik(fecundityff_g5.1t)
fecundityff_ts_aic$AICc[7] <- AICc(fecundityff_g6.1t)
fecundityff_ts_aic$Model[7] <- "m6"
fecundityff_ts_aic$LL[7] <- logLik(fecundityff_g6.1t)
fecundityff_ts_aic$AICc <- as.numeric(fecundityff_ts_aic$AICc)
fecundityff_ts_aic$LL <- as.numeric(fecundityff_ts_aic$LL)
fecundityff_ts_aic <- fecundityff_ts_aic[order(fecundityff_ts_aic$AICc), ]
fecundityff_ts_aic$Delta_AICc <- "0.00"
fecundityff_ts_aic$Delta_AICc[2] <- round(fecundityff_ts_aic$AICc[1]-fecundityff_ts_aic$AICc[2], 2)
fecundityff_ts_aic$Delta_AICc[3] <- round(fecundityff_ts_aic$AICc[1]-fecundityff_ts_aic$AICc[3], 2)
fecundityff_ts_aic$Delta_AICc[4] <- round(fecundityff_ts_aic$AICc[1]-fecundityff_ts_aic$AICc[4], 2)
fecundityff_ts_aic$Delta_AICc[5] <- round(fecundityff_ts_aic$AICc[1]-fecundityff_ts_aic$AICc[5], 2)
fecundityff_ts_aic$Delta_AICc[6] <- round(fecundityff_ts_aic$AICc[1]-fecundityff_ts_aic$AICc[6], 2)
fecundityff_ts_aic$Delta_AICc[7] <- round(fecundityff_ts_aic$AICc[1]-fecundityff_ts_aic$AICc[7], 2)
fecundityff_ts_aic  #Model 2 is best.

fecundityff_ls_aic <- as.data.frame(1:7)
fecundityff_ls_aic$AICc <- "NA"
fecundityff_ls_aic$Model <- "NA"
fecundityff_ls_aic$LL <- "NA"
fecundityff_ls_aic$AICc[1] <- AICc(fecundityff_gl_null)
fecundityff_ls_aic$Model[1] <- "Null"
fecundityff_ls_aic$LL[1] <- logLik(fecundityff_gl_null)
fecundityff_ls_aic$AICc[2] <- AICc(fecundityff_g1l)
fecundityff_ls_aic$Model[2] <- "m1"
fecundityff_ls_aic$LL[2] <- logLik(fecundityff_g1l)
fecundityff_ls_aic$AICc[3] <- AICc(fecundityff_g2.1l)
fecundityff_ls_aic$Model[3] <- "m2"
fecundityff_ls_aic$LL[3] <- logLik(fecundityff_g2.1l)
fecundityff_ls_aic$AICc[4] <- AICc(fecundityff_g3.1l)
fecundityff_ls_aic$Model[4] <- "m3"
fecundityff_ls_aic$LL[4] <- logLik(fecundityff_g3.1l)
fecundityff_ls_aic$AICc[5] <- AICc(fecundityff_g4.1l)
fecundityff_ls_aic$Model[5] <- "m4"
fecundityff_ls_aic$LL[5] <- logLik(fecundityff_g4.1l)
fecundityff_ls_aic$AICc[6] <- AICc(fecundityff_g5.1l)
fecundityff_ls_aic$Model[6] <- "m5"
fecundityff_ls_aic$LL[6] <- logLik(fecundityff_g5.1l)
fecundityff_ls_aic$AICc[7] <- AICc(fecundityff_g6.1l)
fecundityff_ls_aic$Model[7] <- "m6"
fecundityff_ls_aic$LL[7] <- logLik(fecundityff_g6.1l)
fecundityff_ls_aic$AICc <- as.numeric(fecundityff_ls_aic$AICc)
fecundityff_ls_aic$LL <- as.numeric(fecundityff_ls_aic$LL)
fecundityff_ls_aic <- fecundityff_ls_aic[order(fecundityff_ls_aic$AICc), ]
fecundityff_ls_aic$Delta_AICc <- "0.00"
fecundityff_ls_aic$Delta_AICc[2] <- round(fecundityff_ls_aic$AICc[1]-fecundityff_ls_aic$AICc[2], 2)
fecundityff_ls_aic$Delta_AICc[3] <- round(fecundityff_ls_aic$AICc[1]-fecundityff_ls_aic$AICc[3], 2)
fecundityff_ls_aic$Delta_AICc[4] <- round(fecundityff_ls_aic$AICc[1]-fecundityff_ls_aic$AICc[4], 2)
fecundityff_ls_aic$Delta_AICc[5] <- round(fecundityff_ls_aic$AICc[1]-fecundityff_ls_aic$AICc[5], 2)
fecundityff_ls_aic$Delta_AICc[6] <- round(fecundityff_ls_aic$AICc[1]-fecundityff_ls_aic$AICc[6], 2)
fecundityff_ls_aic$Delta_AICc[7] <- round(fecundityff_ls_aic$AICc[1]-fecundityff_ls_aic$AICc[7], 2)
fecundityff_ls_aic  #Model 3 is best

# GAM s(bs = 'tp') optimised
fecundityff_tso_aic <- as.data.frame(1:7)
fecundityff_tso_aic$AICc <- "NA"
fecundityff_tso_aic$Model <- "NA"
fecundityff_tso_aic$LL <- "NA"
fecundityff_tso_aic$AICc[1] <- AICc(fecundityff_gt_null)
fecundityff_tso_aic$Model[1] <- "Null"
fecundityff_tso_aic$LL[1] <- logLik(fecundityff_gt_null)
fecundityff_tso_aic$AICc[2] <- AICc(fecundityff_g1t)
fecundityff_tso_aic$Model[2] <- "m1"
fecundityff_tso_aic$LL[2] <- logLik(fecundityff_g1t)
fecundityff_tso_aic$AICc[3] <- AICc(fecundityff_g2.1ot)
fecundityff_tso_aic$Model[3] <- "m2"
fecundityff_tso_aic$LL[3] <- logLik(fecundityff_g2.1ot)
fecundityff_tso_aic$AICc[4] <- AICc(fecundityff_g3.1ot)
fecundityff_tso_aic$Model[4] <- "m3"
fecundityff_tso_aic$LL[4] <- logLik(fecundityff_g3.1ot)
fecundityff_tso_aic$AICc[5] <- AICc(fecundityff_g4.1ot)
fecundityff_tso_aic$Model[5] <- "m4"
fecundityff_tso_aic$LL[5] <- logLik(fecundityff_g4.1ot)
fecundityff_tso_aic$AICc[6] <- AICc(fecundityff_g5.1ot)
fecundityff_tso_aic$Model[6] <- "m5"
fecundityff_tso_aic$LL[6] <- logLik(fecundityff_g5.1ot)
fecundityff_tso_aic$AICc[7] <- AICc(fecundityff_g6.1ot)
fecundityff_tso_aic$Model[7] <- "m6"
fecundityff_tso_aic$LL[7] <- logLik(fecundityff_g6.1ot)
fecundityff_tso_aic$AICc <- as.numeric(fecundityff_tso_aic$AICc)
fecundityff_tso_aic$LL <- as.numeric(fecundityff_tso_aic$LL)
fecundityff_tso_aic <- fecundityff_tso_aic[order(fecundityff_tso_aic$AICc), ]
fecundityff_tso_aic$Delta_AICc <- "0.00"
fecundityff_tso_aic$Delta_AICc[2] <- round(fecundityff_tso_aic$AICc[1]-fecundityff_tso_aic$AICc[2], 2)
fecundityff_tso_aic$Delta_AICc[3] <- round(fecundityff_tso_aic$AICc[1]-fecundityff_tso_aic$AICc[3], 2)
fecundityff_tso_aic$Delta_AICc[4] <- round(fecundityff_tso_aic$AICc[1]-fecundityff_tso_aic$AICc[4], 2)
fecundityff_tso_aic$Delta_AICc[5] <- round(fecundityff_tso_aic$AICc[1]-fecundityff_tso_aic$AICc[5], 2)
fecundityff_tso_aic$Delta_AICc[6] <- round(fecundityff_tso_aic$AICc[1]-fecundityff_tso_aic$AICc[6], 2)
fecundityff_tso_aic$Delta_AICc[7] <- round(fecundityff_tso_aic$AICc[1]-fecundityff_tso_aic$AICc[7], 2)
fecundityff_tso_aic  #Model 6 is best.

fecundityff_lso_aic <- as.data.frame(1:7)
fecundityff_lso_aic$AICc <- "NA"
fecundityff_lso_aic$Model <- "NA"
fecundityff_lso_aic$LL <- "NA"
fecundityff_lso_aic$AICc[1] <- AICc(fecundityff_gl_null)
fecundityff_lso_aic$Model[1] <- "Null"
fecundityff_lso_aic$LL[1] <- logLik(fecundityff_gl_null)
fecundityff_lso_aic$AICc[2] <- AICc(fecundityff_g1l)
fecundityff_lso_aic$Model[2] <- "m1"
fecundityff_lso_aic$LL[2] <- logLik(fecundityff_g1l)
fecundityff_lso_aic$AICc[3] <- AICc(fecundityff_g2.1l)
fecundityff_lso_aic$Model[3] <- "m2"
fecundityff_lso_aic$LL[3] <- logLik(fecundityff_g2.1l)
fecundityff_lso_aic$AICc[4] <- AICc(fecundityff_g3.1ol)
fecundityff_lso_aic$Model[4] <- "m3"
fecundityff_lso_aic$LL[4] <- logLik(fecundityff_g3.1ol)
fecundityff_lso_aic$AICc[5] <- AICc(fecundityff_g4.1ol)
fecundityff_lso_aic$Model[5] <- "m4"
fecundityff_lso_aic$LL[5] <- logLik(fecundityff_g4.1ol)
fecundityff_lso_aic$AICc[6] <- AICc(fecundityff_g5.1ol)
fecundityff_lso_aic$Model[6] <- "m5"
fecundityff_lso_aic$LL[6] <- logLik(fecundityff_g5.1ol)
fecundityff_lso_aic$AICc[7] <- AICc(fecundityff_g6.1ol)
fecundityff_lso_aic$Model[7] <- "m6"
fecundityff_lso_aic$LL[7] <- logLik(fecundityff_g6.1ol)
fecundityff_lso_aic$AICc <- as.numeric(fecundityff_lso_aic$AICc)
fecundityff_lso_aic$LL <- as.numeric(fecundityff_lso_aic$LL)
fecundityff_lso_aic <- fecundityff_lso_aic[order(fecundityff_lso_aic$AICc), ]
fecundityff_lso_aic$Delta_AICc <- "0.00"
fecundityff_lso_aic$Delta_AICc[2] <- round(fecundityff_lso_aic$AICc[1]-fecundityff_lso_aic$AICc[2], 2)
fecundityff_lso_aic$Delta_AICc[3] <- round(fecundityff_lso_aic$AICc[1]-fecundityff_lso_aic$AICc[3], 2)
fecundityff_lso_aic$Delta_AICc[4] <- round(fecundityff_lso_aic$AICc[1]-fecundityff_lso_aic$AICc[4], 2)
fecundityff_lso_aic$Delta_AICc[5] <- round(fecundityff_lso_aic$AICc[1]-fecundityff_lso_aic$AICc[5], 2)
fecundityff_lso_aic$Delta_AICc[6] <- round(fecundityff_lso_aic$AICc[1]-fecundityff_lso_aic$AICc[6], 2)
fecundityff_lso_aic$Delta_AICc[7] <- round(fecundityff_lso_aic$AICc[1]-fecundityff_lso_aic$AICc[7], 2)
fecundityff_lso_aic  #Model 3 is best


# GAM ti(bs='cr')
fecundityff_tsc_aic <- as.data.frame(1:7)
fecundityff_tsc_aic$AICc <- "NA"
fecundityff_tsc_aic$Model <- "NA"
fecundityff_tsc_aic$LL <- "NA"
fecundityff_tsc_aic$AICc[1] <- AICc(fecundityff_gt_null)
fecundityff_tsc_aic$Model[1] <- "Null"
fecundityff_tsc_aic$LL[1] <- logLik(fecundityff_gt_null)
fecundityff_tsc_aic$AICc[2] <- AICc(fecundityff_g1.1t)
fecundityff_tsc_aic$Model[2] <- "m1"
fecundityff_tsc_aic$LL[2] <- logLik(fecundityff_g1.1t)
fecundityff_tsc_aic$AICc[3] <- AICc(fecundityff_g2.2t)
fecundityff_tsc_aic$Model[3] <- "m2"
fecundityff_tsc_aic$LL[3] <- logLik(fecundityff_g2.2t)
fecundityff_tsc_aic$AICc[4] <- AICc(fecundityff_g3.2t)
fecundityff_tsc_aic$Model[4] <- "m3"
fecundityff_tsc_aic$LL[4] <- logLik(fecundityff_g3.2t)
fecundityff_tsc_aic$AICc[5] <- AICc(fecundityff_g4.2t)
fecundityff_tsc_aic$Model[5] <- "m4"
fecundityff_tsc_aic$LL[5] <- logLik(fecundityff_g4.2t)
fecundityff_tsc_aic$AICc[6] <- AICc(fecundityff_g5.2t)
fecundityff_tsc_aic$Model[6] <- "m5"
fecundityff_tsc_aic$LL[6] <- logLik(fecundityff_g5.2t)
fecundityff_tsc_aic$AICc[7] <- AICc(fecundityff_g6.2t)
fecundityff_tsc_aic$Model[7] <- "m6"
fecundityff_tsc_aic$LL[7] <- logLik(fecundityff_g6.2t)
fecundityff_tsc_aic$AICc <- as.numeric(fecundityff_tsc_aic$AICc)
fecundityff_tsc_aic$LL <- as.numeric(fecundityff_tsc_aic$LL)
fecundityff_tsc_aic <- fecundityff_tsc_aic[order(fecundityff_tsc_aic$AICc), ]
fecundityff_tsc_aic$Delta_AICc <- "0.00"
fecundityff_tsc_aic$Delta_AICc[2] <- round(fecundityff_tsc_aic$AICc[1]-fecundityff_tsc_aic$AICc[2], 2)
fecundityff_tsc_aic$Delta_AICc[3] <- round(fecundityff_tsc_aic$AICc[1]-fecundityff_tsc_aic$AICc[3], 2)
fecundityff_tsc_aic$Delta_AICc[4] <- round(fecundityff_tsc_aic$AICc[1]-fecundityff_tsc_aic$AICc[4], 2)
fecundityff_tsc_aic$Delta_AICc[5] <- round(fecundityff_tsc_aic$AICc[1]-fecundityff_tsc_aic$AICc[5], 2)
fecundityff_tsc_aic$Delta_AICc[6] <- round(fecundityff_tsc_aic$AICc[1]-fecundityff_tsc_aic$AICc[6], 2)
fecundityff_tsc_aic$Delta_AICc[7] <- round(fecundityff_tsc_aic$AICc[1]-fecundityff_tsc_aic$AICc[7], 2)
fecundityff_tsc_aic  #Model 6 is best

fecundityff_lsc_aic <- as.data.frame(1:7)
fecundityff_lsc_aic$AICc <- "NA"
fecundityff_lsc_aic$Model <- "NA"
fecundityff_lsc_aic$LL <- "NA"
fecundityff_lsc_aic$AICc[1] <- AICc(fecundityff_gl_null)
fecundityff_lsc_aic$Model[1] <- "Null"
fecundityff_lsc_aic$LL[1] <- logLik(fecundityff_gl_null)
fecundityff_lsc_aic$AICc[2] <- AICc(fecundityff_g1.1l)
fecundityff_lsc_aic$Model[2] <- "m1"
fecundityff_lsc_aic$LL[2] <- logLik(fecundityff_g1.1l)
fecundityff_lsc_aic$AICc[3] <- AICc(fecundityff_g2.2l)
fecundityff_lsc_aic$Model[3] <- "m2"
fecundityff_lsc_aic$LL[3] <- logLik(fecundityff_g2.2l)
fecundityff_lsc_aic$AICc[4] <- AICc(fecundityff_g3.2l)
fecundityff_lsc_aic$Model[4] <- "m3"
fecundityff_lsc_aic$LL[4] <- logLik(fecundityff_g3.2l)
fecundityff_lsc_aic$AICc[5] <- AICc(fecundityff_g4.2l)
fecundityff_lsc_aic$Model[5] <- "m4"
fecundityff_lsc_aic$LL[5] <- logLik(fecundityff_g4.2l)
fecundityff_lsc_aic$AICc[6] <- AICc(fecundityff_g5.2l)
fecundityff_lsc_aic$Model[6] <- "m5"
fecundityff_lsc_aic$LL[6] <- logLik(fecundityff_g5.2l)
fecundityff_lsc_aic$AICc[7] <- AICc(fecundityff_g6.2l)
fecundityff_lsc_aic$Model[7] <- "m6"
fecundityff_lsc_aic$LL[7] <- logLik(fecundityff_g6.2l)
fecundityff_lsc_aic$AICc <- as.numeric(fecundityff_lsc_aic$AICc)
fecundityff_lsc_aic$LL <- as.numeric(fecundityff_lsc_aic$LL)
fecundityff_lsc_aic <- fecundityff_lsc_aic[order(fecundityff_lsc_aic$AICc), ]
fecundityff_lsc_aic$Delta_AICc <- "0.00"
fecundityff_lsc_aic$Delta_AICc[2] <- round(fecundityff_lsc_aic$AICc[1]-fecundityff_lsc_aic$AICc[2], 2)
fecundityff_lsc_aic$Delta_AICc[3] <- round(fecundityff_lsc_aic$AICc[1]-fecundityff_lsc_aic$AICc[3], 2)
fecundityff_lsc_aic$Delta_AICc[4] <- round(fecundityff_lsc_aic$AICc[1]-fecundityff_lsc_aic$AICc[4], 2)
fecundityff_lsc_aic$Delta_AICc[5] <- round(fecundityff_lsc_aic$AICc[1]-fecundityff_lsc_aic$AICc[5], 2)
fecundityff_lsc_aic$Delta_AICc[6] <- round(fecundityff_lsc_aic$AICc[1]-fecundityff_lsc_aic$AICc[6], 2)
fecundityff_lsc_aic$Delta_AICc[7] <- round(fecundityff_lsc_aic$AICc[1]-fecundityff_lsc_aic$AICc[7], 2)
fecundityff_lsc_aic  # Model 3 is best


# GAM (bs = 'cr') optimised
fecundityff_tsco_aic <- as.data.frame(1:7)
fecundityff_tsco_aic$AICc <- "NA"
fecundityff_tsco_aic$Model <- "NA"
fecundityff_tsco_aic$LL <- "NA"
fecundityff_tsco_aic$AICc[1] <- AICc(fecundityff_gt_null)
fecundityff_tsco_aic$Model[1] <- "Null"
fecundityff_tsco_aic$LL[1] <- logLik(fecundityff_gt_null)
fecundityff_tsco_aic$AICc[2] <- AICc(fecundityff_g1.1t)
fecundityff_tsco_aic$Model[2] <- "m1"
fecundityff_tsco_aic$LL[2] <- logLik(fecundityff_g1.1t)
fecundityff_tsco_aic$AICc[3] <- AICc(fecundityff_g2.2ot)
fecundityff_tsco_aic$Model[3] <- "m2"
fecundityff_tsco_aic$LL[3] <- logLik(fecundityff_g2.2ot)
fecundityff_tsco_aic$AICc[4] <- AICc(fecundityff_g3.2ot)
fecundityff_tsco_aic$Model[4] <- "m3"
fecundityff_tsco_aic$LL[4] <- logLik(fecundityff_g3.2ot)
fecundityff_tsco_aic$AICc[5] <- AICc(fecundityff_g4.2ot)
fecundityff_tsco_aic$Model[5] <- "m4"
fecundityff_tsco_aic$LL[5] <- logLik(fecundityff_g4.2ot)
fecundityff_tsco_aic$AICc[6] <- AICc(fecundityff_g5.2ot)
fecundityff_tsco_aic$Model[6] <- "m5"
fecundityff_tsco_aic$LL[6] <- logLik(fecundityff_g5.2ot)
fecundityff_tsco_aic$AICc[7] <- AICc(fecundityff_g6.2t)
fecundityff_tsco_aic$Model[7] <- "m6"
fecundityff_tsco_aic$LL[7] <- logLik(fecundityff_g6.2t)
fecundityff_tsco_aic$AICc <- as.numeric(fecundityff_tsco_aic$AICc)
fecundityff_tsco_aic$LL <- as.numeric(fecundityff_tsco_aic$LL)
fecundityff_tsco_aic <- fecundityff_tsco_aic[order(fecundityff_tsco_aic$AICc), ]
fecundityff_tsco_aic$Delta_AICc <- "0.00"
fecundityff_tsco_aic$Delta_AICc[2] <- round(fecundityff_tsco_aic$AICc[1]-fecundityff_tsco_aic$AICc[2], 2)
fecundityff_tsco_aic$Delta_AICc[3] <- round(fecundityff_tsco_aic$AICc[1]-fecundityff_tsco_aic$AICc[3], 2)
fecundityff_tsco_aic$Delta_AICc[4] <- round(fecundityff_tsco_aic$AICc[1]-fecundityff_tsco_aic$AICc[4], 2)
fecundityff_tsco_aic$Delta_AICc[5] <- round(fecundityff_tsco_aic$AICc[1]-fecundityff_tsco_aic$AICc[5], 2)
fecundityff_tsco_aic$Delta_AICc[6] <- round(fecundityff_tsco_aic$AICc[1]-fecundityff_tsco_aic$AICc[6], 2)
fecundityff_tsco_aic$Delta_AICc[7] <- round(fecundityff_tsco_aic$AICc[1]-fecundityff_tsco_aic$AICc[7], 2)
fecundityff_tsco_aic  #Model 6 is best.

fecundityff_lsco_aic <- as.data.frame(1:7)
fecundityff_lsco_aic$AICc <- "NA"
fecundityff_lsco_aic$Model <- "NA"
fecundityff_lsco_aic$LL <- "NA"
fecundityff_lsco_aic$AICc[1] <- AICc(fecundityff_gl_null)
fecundityff_lsco_aic$Model[1] <- "Null"
fecundityff_lsco_aic$LL[1] <- logLik(fecundityff_gl_null)
fecundityff_lsco_aic$AICc[2] <- AICc(fecundityff_g1.1l)
fecundityff_lsco_aic$Model[2] <- "m1"
fecundityff_lsco_aic$LL[2] <- logLik(fecundityff_g1.1l)
fecundityff_lsco_aic$AICc[3] <- AICc(fecundityff_g2.2l)
fecundityff_lsco_aic$Model[3] <- "m2"
fecundityff_lsco_aic$LL[3] <- logLik(fecundityff_g2.2l)
fecundityff_lsco_aic$AICc[4] <- AICc(fecundityff_g3.2l)
fecundityff_lsco_aic$Model[4] <- "m3"
fecundityff_lsco_aic$LL[4] <- logLik(fecundityff_g3.2l)
fecundityff_lsco_aic$AICc[5] <- AICc(fecundityff_g4.2l)
fecundityff_lsco_aic$Model[5] <- "m4"
fecundityff_lsco_aic$LL[5] <- logLik(fecundityff_g4.2l)
fecundityff_lsco_aic$AICc[6] <- AICc(fecundityff_g5.2l)
fecundityff_lsco_aic$Model[6] <- "m5"
fecundityff_lsco_aic$LL[6] <- logLik(fecundityff_g5.2l)
fecundityff_lsco_aic$AICc[7] <- AICc(fecundityff_g6.2l)
fecundityff_lsco_aic$Model[7] <- "m6"
fecundityff_lsco_aic$LL[7] <- logLik(fecundityff_g6.2l)
fecundityff_lsco_aic$AICc <- as.numeric(fecundityff_lsco_aic$AICc)
fecundityff_lsco_aic$LL <- as.numeric(fecundityff_lsco_aic$LL)
fecundityff_lsco_aic <- fecundityff_lsco_aic[order(fecundityff_lsco_aic$AICc), ]
fecundityff_lsco_aic$Delta_AICc <- "0.00"
fecundityff_lsco_aic$Delta_AICc[2] <- round(fecundityff_lsco_aic$AICc[1]-fecundityff_lsco_aic$AICc[2], 2)
fecundityff_lsco_aic$Delta_AICc[3] <- round(fecundityff_lsco_aic$AICc[1]-fecundityff_lsco_aic$AICc[3], 2)
fecundityff_lsco_aic$Delta_AICc[4] <- round(fecundityff_lsco_aic$AICc[1]-fecundityff_lsco_aic$AICc[4], 2)
fecundityff_lsco_aic$Delta_AICc[5] <- round(fecundityff_lsco_aic$AICc[1]-fecundityff_lsco_aic$AICc[5], 2)
fecundityff_lsco_aic$Delta_AICc[6] <- round(fecundityff_lsco_aic$AICc[1]-fecundityff_lsco_aic$AICc[6], 2)
fecundityff_lsco_aic$Delta_AICc[7] <- round(fecundityff_lsco_aic$AICc[1]-fecundityff_lsco_aic$AICc[7], 2)
fecundityff_lsco_aic  #Model 3 is best


# Predict for female fecundity and fire frequency ----
# For torulosa GLMER m3
# For torulosa GAM default, tp opt, cr and cr opt m6
# For torulosa GAM tp m2

# For littoralis GLMER and GAM all but default m3
# For littoralis GAM default m6


lffec_m3_l <- data.frame(Fire_freq = min(lit_fecundity$Fire_freq),
                       FPC = seq(min(lit_fecundity$FPC, na.rm = T), max(lit_fecundity$FPC, na.rm = T), length = 50))
lffec_m3_a <- data.frame(Fire_freq = mean(lit_fecundity$Fire_freq),
                         FPC = seq(min(lit_fecundity$FPC, na.rm = T), max(lit_fecundity$FPC, na.rm = T), length = 50))
lffec_m3_h <- data.frame(Fire_freq = max(lit_fecundity$Fire_freq),
                         FPC = seq(min(lit_fecundity$FPC, na.rm = T), max(lit_fecundity$FPC, na.rm = T), length = 50))

lffec_m3lm_l <- data.frame(r_Fire_freq = min(lit_fecundity$r_Fire_freq),
                         r_FPC = seq(min(lit_fecundity$r_FPC, na.rm = T), max(lit_fecundity$r_FPC, na.rm = T), length = 50))
lffec_m3lm_a <- data.frame(r_Fire_freq = mean(lit_fecundity$Fire_freq),
                         r_FPC = seq(min(lit_fecundity$r_FPC, na.rm = T), max(lit_fecundity$r_FPC, na.rm = T), length = 50))
lffec_m3lm_h <- data.frame(r_Fire_freq = max(lit_fecundity$Fire_freq),
                         r_FPC = seq(min(lit_fecundity$r_FPC, na.rm = T), max(lit_fecundity$r_FPC, na.rm = T), length = 50))



tffec_m3lm_l <- data.frame(r_Fire_freq = min(tor_fecundity$r_Fire_freq),
                           r_FPC = seq(min(tor_fecundity$r_FPC, na.rm = T), max(tor_fecundity$r_FPC, na.rm = T), length = 50))
tffec_m3lm_a <- data.frame(r_Fire_freq = mean(tor_fecundity$Fire_freq),
                           r_FPC = seq(min(tor_fecundity$r_FPC, na.rm = T), max(tor_fecundity$r_FPC, na.rm = T), length = 50))
tffec_m3lm_h <- data.frame(r_Fire_freq = max(tor_fecundity$Fire_freq),
                           r_FPC = seq(min(tor_fecundity$r_FPC, na.rm = T), max(tor_fecundity$r_FPC, na.rm = T), length = 50))



# We've already decided that GLMER is best for littoralis so let's have a quick look at one prediction from a GAM for torulosa to see if the fit is reasonable
tffec_m6_l <- data.frame(Fire_freq = min(lit_fecundity$Fire_freq),
                         Height_cm = seq(min(lit_fecundity$Height_cm, na.rm = T), max(lit_fecundity$Height_cm, na.rm = T), length = 50))
tffec_m6_a <- data.frame(Fire_freq = mean(lit_fecundity$Fire_freq),
                         Height_cm = seq(min(lit_fecundity$Height_cm, na.rm = T), max(lit_fecundity$Height_cm, na.rm = T), length = 50))
tffec_m6_h <- data.frame(Fire_freq = max(lit_fecundity$Fire_freq),
                         Height_cm = seq(min(lit_fecundity$Height_cm, na.rm = T), max(lit_fecundity$Height_cm, na.rm = T), length = 50))

ptfec_m6_l <- predict(fecundityff_g6t, newdata = tffec_m6_l, se.fit = T, type = 'response')
tffec_m6_l$fit <- ptfec_m6_l$fit
tffec_m6_l$se <- ptfec_m6_l$se.fit
tffec_m6_l$lci <- tffec_m6_l$fit - (tffec_m6_l$se * 1.96)
tffec_m6_l$uci <- tffec_m6_l$fit + (tffec_m6_l$se * 1.96)
# Not reasonable, we get ridiculously high numbers for fit when height is small. This is not what we would expect. Only modelling GLMER for fecundity





# GLMER
pffec_m3lm_l <- predictSE(fecundityff_m3l, newdata = lffec_m3lm_l, se.fit = T, type = 'link')
lffec_m3lm_l$fit.link <- pffec_m3lm_l$fit
lffec_m3lm_l$se.link <- pffec_m3lm_l$se.fit
lffec_m3lm_l$lci.link <- lffec_m3lm_l$fit.link - (lffec_m3lm_l$se.link * 1.96)
lffec_m3lm_l$uci.link <- lffec_m3lm_l$fit.link + (lffec_m3lm_l$se.link * 1.96)
lffec_m3lm_l$fit <- exp(lffec_m3lm_l$fit.link)
lffec_m3lm_l$se <- exp(lffec_m3lm_l$se.link)
lffec_m3lm_l$lci <- exp(lffec_m3lm_l$lci.link)
lffec_m3lm_l$uci <- exp(lffec_m3lm_l$uci.link)
lffec_m3lm_l # Fit does not go into the negatives. 

pffec_m3lm_a <- predictSE(fecundityff_m3l, newdata = lffec_m3lm_a, se.fit = T, type = 'link')
lffec_m3lm_a$fit.link <- pffec_m3lm_a$fit
lffec_m3lm_a$se.link <- pffec_m3lm_a$se.fit
lffec_m3lm_a$lci.link <- lffec_m3lm_a$fit.link - (lffec_m3lm_a$se.link * 1.96)
lffec_m3lm_a$uci.link <- lffec_m3lm_a$fit.link + (lffec_m3lm_a$se.link * 1.96)
lffec_m3lm_a$fit <- exp(lffec_m3lm_a$fit.link)
lffec_m3lm_a$se <- exp(lffec_m3lm_a$se.link)
lffec_m3lm_a$lci <- exp(lffec_m3lm_a$lci.link)
lffec_m3lm_a$uci <- exp(lffec_m3lm_a$uci.link)

pffec_m3lm_h <- predictSE(fecundityff_m3l, newdata = lffec_m3lm_h, se.fit = T, type = 'link')
lffec_m3lm_h$fit.link <- pffec_m3lm_h$fit
lffec_m3lm_h$se.link <- pffec_m3lm_h$se.fit
lffec_m3lm_h$lci.link <- lffec_m3lm_h$fit.link - (lffec_m3lm_h$se.link * 1.96)
lffec_m3lm_h$uci.link <- lffec_m3lm_h$fit.link + (lffec_m3lm_h$se.link * 1.96)
lffec_m3lm_h$fit <- exp(lffec_m3lm_h$fit.link)
lffec_m3lm_h$se <- exp(lffec_m3lm_h$se.link)
lffec_m3lm_h$lci <- exp(lffec_m3lm_h$lci.link)
lffec_m3lm_h$uci <- exp(lffec_m3lm_h$uci.link)



ptffec_m3lm_l <- predictSE(fecundityff_m3l, newdata = tffec_m3lm_l, se.fit = T, type = 'link')
tffec_m3lm_l$fit.link <- ptffec_m3lm_l$fit
tffec_m3lm_l$se.link <- ptffec_m3lm_l$se.fit
tffec_m3lm_l$lci.link <- tffec_m3lm_l$fit.link - (tffec_m3lm_l$se.link * 1.96)
tffec_m3lm_l$uci.link <- tffec_m3lm_l$fit.link + (tffec_m3lm_l$se.link * 1.96)
tffec_m3lm_l$fit <- exp(tffec_m3lm_l$fit.link)
tffec_m3lm_l$se <- exp(tffec_m3lm_l$se.link)
tffec_m3lm_l$lci <- exp(tffec_m3lm_l$lci.link)
tffec_m3lm_l$uci <- exp(tffec_m3lm_l$uci.link)
tffec_m3lm_l # Fit does not go into the negatives. 

ptffec_m3lm_a <- predictSE(fecundityff_m3l, newdata = tffec_m3lm_a, se.fit = T, type = 'link')
tffec_m3lm_a$fit.link <- ptffec_m3lm_a$fit
tffec_m3lm_a$se.link <- ptffec_m3lm_a$se.fit
tffec_m3lm_a$lci.link <- tffec_m3lm_a$fit.link - (tffec_m3lm_a$se.link * 1.96)
tffec_m3lm_a$uci.link <- tffec_m3lm_a$fit.link + (tffec_m3lm_a$se.link * 1.96)
tffec_m3lm_a$fit <- exp(tffec_m3lm_a$fit.link)
tffec_m3lm_a$se <- exp(tffec_m3lm_a$se.link)
tffec_m3lm_a$lci <- exp(tffec_m3lm_a$lci.link)
tffec_m3lm_a$uci <- exp(tffec_m3lm_a$uci.link)

ptffec_m3lm_h <- predictSE(fecundityff_m3l, newdata = tffec_m3lm_h, se.fit = T, type = 'link')
tffec_m3lm_h$fit.link <- ptffec_m3lm_h$fit
tffec_m3lm_h$se.link <- ptffec_m3lm_h$se.fit
tffec_m3lm_h$lci.link <- tffec_m3lm_h$fit.link - (tffec_m3lm_h$se.link * 1.96)
tffec_m3lm_h$uci.link <- tffec_m3lm_h$fit.link + (tffec_m3lm_h$se.link * 1.96)
tffec_m3lm_h$fit <- exp(tffec_m3lm_h$fit.link)
tffec_m3lm_h$se <- exp(tffec_m3lm_h$se.link)
tffec_m3lm_h$lci <- exp(tffec_m3lm_h$lci.link)
tffec_m3lm_h$uci <- exp(tffec_m3lm_h$uci.link)

# GAM tp
lffec_m3t_l <- lffec_m3_l
lffec_m3t_a <- lffec_m3_a
lffec_m3t_h <- lffec_m3_h


pffec_m3t_l <- predict(fecundityff_g3.1l, newdata = lffec_m3t_l, se.fit = T, type = 'link')
lffec_m3t_l$fit <- pffec_m3t_l$fit
lffec_m3t_l$se <- pffec_m3t_l$se.fit
lffec_m3t_l$lci <- lffec_m3t_l$fit - (lffec_m3t_l$se * 1.96)
lffec_m3t_l$uci <- lffec_m3t_l$fit + (lffec_m3t_l$se * 1.96)
lffec_m3t_l # Fit goes into negatives, this model is not good.

pffec_m3t_a <- predict(fecundityff_g3.1l, newdata = lffec_m3t_a, se.fit = T, type = 'link')
lffec_m3t_a$fit <- pffec_m3t_a$fit
lffec_m3t_a$se <- pffec_m3t_a$se.fit
lffec_m3t_a$lci <- lffec_m3t_a$fit - (lffec_m3t_a$se * 1.96)
lffec_m3t_a$uci <- lffec_m3t_a$fit + (lffec_m3t_a$se * 1.96)


pffec_m3t_h <- predict(fecundityff_g3.1l, newdata = lffec_m3t_h, se.fit = T, type = 'link')
lffec_m3t_h$fit <- pffec_m3t_h$fit
lffec_m3t_h$se <- pffec_m3t_h$se.fit
lffec_m3t_h$lci <- lffec_m3t_h$fit - (lffec_m3t_h$se * 1.96)
lffec_m3t_h$uci <- lffec_m3t_h$fit + (lffec_m3t_h$se * 1.96)



# GAM tp opt
lffec_m3to_l <- lffec_m3_l
lffec_m3to_a <- lffec_m3_a
lffec_m3to_h <- lffec_m3_h

pffec_m3to_l <- predict(fecundityff_g3.1ol, newdata = lffec_m3to_l, se.fit = T, type = 'link')
lffec_m3to_l$fit <- pffec_m3to_l$fit
lffec_m3to_l$se <- pffec_m3to_l$se.fit
lffec_m3to_l$lci <- lffec_m3to_l$fit - (lffec_m3to_l$se * 1.96)
lffec_m3to_l$uci <- lffec_m3to_l$fit + (lffec_m3to_l$se * 1.96)
lffec_m3to_l # Fit goes into negatives, this model is not good


pffec_m3to_a <- predict(fecundityff_g3.1ol, newdata = lffec_m3to_a, se.fit = T, type = 'link')
lffec_m3to_a$fit <- pffec_m3to_a$fit
lffec_m3to_a$se <- pffec_m3to_a$se.fit
lffec_m3to_a$lci <- lffec_m3to_a$fit - (lffec_m3to_a$se * 1.96)
lffec_m3to_a$uci <- lffec_m3to_a$fit + (lffec_m3to_a$se * 1.96)
lffec_m3to_a$fit <- exp(lffec_m3to_a$fit)
lffec_m3to_a$se <- exp(lffec_m3to_a$se)
lffec_m3to_a$lci <- exp(lffec_m3to_a$lci)
lffec_m3to_a$uci <- exp(lffec_m3to_a$uci)

pffec_m3to_h <- predict(fecundityff_g3.1ol, newdata = lffec_m3to_h, se.fit = T, type = 'link')
lffec_m3to_h$fit <- pffec_m3to_h$fit
lffec_m3to_h$se <- pffec_m3to_h$se.fit
lffec_m3to_h$lci <- lffec_m3to_h$fit - (lffec_m3to_h$se * 1.96)
lffec_m3to_h$uci <- lffec_m3to_h$fit + (lffec_m3to_h$se * 1.96)



# GAM cr
lffec_m3c_l <- lffec_m3_l
lffec_m3c_a <- lffec_m3_a
lffec_m3c_h <- lffec_m3_h

pffec_m3c_l <- predict(fecundityff_g3.2l, newdata = lffec_m3c_l, se.fit = T, type = 'link')
lffec_m3c_l$fit <- pffec_m3c_l$fit
lffec_m3c_l$se <- pffec_m3c_l$se.fit
lffec_m3c_l$lci <- lffec_m3c_l$fit - (lffec_m3c_l$se * 1.96)
lffec_m3c_l$uci <- lffec_m3c_l$fit + (lffec_m3c_l$se * 1.96)
lffec_m3c_l #Fit does not go into negatives

pffec_m3c_a <- predict(fecundityff_g3.2l, newdata = lffec_m3c_a, se.fit = T, type = 'link')
lffec_m3c_a$fit <- pffec_m3c_a$fit
lffec_m3c_a$se <- pffec_m3c_a$se.fit
lffec_m3c_a$lci <- lffec_m3c_a$fit - (lffec_m3c_a$se * 1.96)
lffec_m3c_a$uci <- lffec_m3c_a$fit + (lffec_m3c_a$se * 1.96)


pffec_m3c_h <- predict(fecundityff_g3.2l, newdata = lffec_m3c_h, se.fit = T, type = 'link')
lffec_m3c_h$fit <- pffec_m3c_h$fit
lffec_m3c_h$se <- pffec_m3c_h$se.fit
lffec_m3c_h$lci <- lffec_m3c_h$fit - (lffec_m3c_h$se * 1.96)
lffec_m3c_h$uci <- lffec_m3c_h$fit + (lffec_m3c_h$se * 1.96)

# GAM cr opt
lffec_m3co_l <- lffec_m3_l
lffec_m3co_a <- lffec_m3_a
lffec_m3co_h <- lffec_m3_h

pffec_m3co_l <- predict(fecundityff_g3.2ol, newdata = lffec_m3co_l, se.fit = T, type = 'link')
lffec_m3co_l$fit <- pffec_m3co_l$fit
lffec_m3co_l$se <- pffec_m3co_l$se.fit
lffec_m3co_l$lci <- lffec_m3co_l$fit - (lffec_m3co_l$se * 1.96)
lffec_m3co_l$uci <- lffec_m3co_l$fit + (lffec_m3co_l$se * 1.96)
lffec_m3co_l # Fit does not go into negatives

pffec_m3co_a <- predict(fecundityff_g3.2ol, newdata = lffec_m3co_a, se.fit = T, type = 'link')
lffec_m3co_a$fit <- pffec_m3co_a$fit
lffec_m3co_a$se <- pffec_m3co_a$se.fit
lffec_m3co_a$lci <- lffec_m3co_a$fit - (lffec_m3co_a$se * 1.96)
lffec_m3co_a$uci <- lffec_m3co_a$fit + (lffec_m3co_a$se * 1.96)
# Fit goes into negatives

pffec_m3co_h <- predict(fecundityff_g3.2ol, newdata = lffec_m3co_h, se.fit = T, type = 'link')
lffec_m3co_h$fit <- pffec_m3co_h$fit
lffec_m3co_h$se <- pffec_m3co_h$se.fit
lffec_m3co_h$lci <- lffec_m3co_h$fit - (lffec_m3co_h$se * 1.96)
lffec_m3co_h$uci <- lffec_m3co_h$fit + (lffec_m3co_h$se * 1.96)
# fit is negative

# Plot predictions for female fecundity and FF ----
# NOTE for littoralis we are not plotting GAMs with bs = tp as the fit went into the negatives which we do not want. GLMs produce the more reasonable plots here. Predictions 

dev.new(width = 20, height = 12, noRStudioGD = T, dpi = 300)
par(mfrow = c(2,3), mar = c(8,6,3,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))

plot(lffec_m3c_l$FPC, lffec_m3c_l$fit, type = 'l', col = 'blue', ylab = "", xlab = "", las = 1, ylim = c(1, 6))
pg.ci(x = 'FPC', data = 'lffec_m3c_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lffec_m3c_a$FPC, lffec_m3c_a$fit, type = 'l', col = 'black')
lines(lffec_m3c_h$FPC, lffec_m3c_h$fit, type = 'l', col = 'red')

plot(lffec_m3lm_l$r_FPC, lffec_m3lm_l$fit, type = 'l', col = 'blue', ylab = "", xlab = "", las = 1, ylim = c(0, 500))
pg.ci(x = 'r_FPC', data = 'lffec_m3lm_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lffec_m3lm_a$r_FPC, lffec_m3lm_a$fit, type = 'l', col = 'black')
pg.ci(x = 'r_FPC', data = 'lffec_m3lm_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lffec_m3lm_h$r_FPC, lffec_m3lm_h$fit, type = 'l', col = 'red')
pg.ci(x = 'r_FPC', data = 'lffec_m3lm_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')




dev.new(width = 12, height = 10, noRStudioGD = T, dpi = 300)
par(mfrow = c(1,2), mar = c(8,6,3,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))

plot(lffec_m3lm_l$r_FPC, lffec_m3lm_l$fit, type = 'l', col = 'blue', ylab = "", xlab = "", las = 1, ylim = c(0, 500))
pg.ci(x = 'r_FPC', data = 'lffec_m3lm_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lffec_m3lm_a$r_FPC, lffec_m3lm_a$fit, type = 'l', col = 'black')
pg.ci(x = 'r_FPC', data = 'lffec_m3lm_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lffec_m3lm_h$r_FPC, lffec_m3lm_h$fit, type = 'l', col = 'red')
pg.ci(x = 'r_FPC', data = 'lffec_m3lm_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(tffec_m3lm_l$r_FPC, tffec_m3lm_l$fit, type = 'l', col = 'blue', ylab = "", xlab = "", las = 1, ylim = c(0, 500))
pg.ci(x = 'r_FPC', data = 'tffec_m3lm_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(tffec_m3lm_a$r_FPC, tffec_m3lm_a$fit, type = 'l', col = 'black')
pg.ci(x = 'r_FPC', data = 'tffec_m3lm_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(tffec_m3lm_h$r_FPC, tffec_m3lm_h$fit, type = 'l', col = 'red')
pg.ci(x = 'r_FPC', data = 'tffec_m3lm_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')



# 6.3 Seed size ----
tor_cum.prop$r_Latitude <- scale(tor_cum.prop$Latitude)
lit_cum.prop$r_Latitude <- scale(lit_cum.prop$Latitude)
tor_cum.prop$r_FPC <- scale(tor_cum.prop$FPC)
lit_cum.prop$r_FPC <- scale(lit_cum.prop$FPC)
tor_cum.prop$r_Temp <- scale(tor_cum.prop$Temp)
lit_cum.prop$r_Temp <- scale(lit_cum.prop$Temp)
tor_cum.prop$r_Precip <- scale(tor_cum.prop$Temp)
lit_cum.prop$r_Precip <- scale(lit_cum.prop$Precip)

seedwt_t_null <- glmer(seed_wt_mg ~ 1 + (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
seedwt_t_gnull <- gam(seed_wt_mg ~ 1, random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')

seedwt_l_null <- glmer(seed_wt_mg ~ 1 + (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
seedwt_l_gnull <- gam(seed_wt_mg ~ 1, random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')



seedwt_m1t <- glmer(seed_wt_mg ~ r_fire_freq + (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m1t)
seedwt_g1t <- gam(seed_wt_mg ~ s(Fire_freq, k = 7), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g1t)
plot(seedwt_g1t)
seedwt_g1.1t <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g1.1t)
plot(seedwt_g1.1t)
seedwt_g1.1ot <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 3), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g1.1ot)
plot(seedwt_g1.1ot)

seedwt_m1l <- glmer(seed_wt_mg ~ r_fire_freq + (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m1l)
seedwt_g1l <- gam(seed_wt_mg ~ s(Fire_freq, k = 3), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g1l)
plot(seedwt_g1l)
seedwt_g1.1l <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 3), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g1.1l)
plot(seedwt_g1.1l)



seedwt_m2t <- glmer(seed_wt_mg ~ r_fire_freq * r_Latitude + (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m2t)
seedwt_g2t <- gam(seed_wt_mg ~ s(Fire_freq, k = 7) + s(Latitude) + ti(Latitude, by = Fire_freq), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g2t)
plot(seedwt_g2t)
summary(seedwt_g2t)
seedwt_g2.1t <- gam(seed_wt_mg ~ s(Fire_freq, k = 7) + s(Latitude) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g2.1t)
plot(seedwt_g2.1t)
seedwt_g2.2t <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 5) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g2.2t)
plot(seedwt_g2.2t)
seedwt_g2.1ot <- gam(seed_wt_mg ~ s(Fire_freq, k = 6) + s(Latitude, k = 5) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g2.1ot)
plot(seedwt_g2.1ot)
seedwt_g2.2ot <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 4) + s(Latitude, bs = 'cr', k = 3) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 4), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g2.2ot)
plot(seedwt_g2.2ot)


seedwt_m2l <- glmer(seed_wt_mg ~ r_fire_freq * r_Latitude + (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop) 
summary(seedwt_m2l)
seedwt_g2l <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(Latitude) + ti(Latitude, by = Fire_freq), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g2l)
plot(seedwt_g2l)
seedwt_g2.1l <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(Latitude) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g2.1l)
plot(seedwt_g2.1l)
seedwt_g2.2l <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 3) + s(Latitude, bs = 'cr', k = 5) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g2.1l)
plot(seedwt_g2.2l)
seedwt_g2.1ol <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(Latitude, k = 3) + ti(Latitude, by = Fire_freq, bs = 'tp', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g2.1ol)
plot(seedwt_g2.1ol)
seedwt_g2.2ol <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 3) + s(Latitude, bs = 'cr', k = 4) + ti(Latitude, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g2.1ol)
plot(seedwt_g2.2ol)


seedwt_m3t <- glmer(seed_wt_mg ~ r_fire_freq * r_FPC + (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m3t)
seedwt_g3t <- gam(seed_wt_mg ~ s(Fire_freq, k = 7) + s(FPC) + ti(FPC, by = Fire_freq), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g3t)
plot(seedwt_g3t)
seedwt_g3.1t <- gam(seed_wt_mg ~ s(Fire_freq, k = 7) + s(FPC) + ti(FPC, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g3.1t)
plot(seedwt_g3.1t)
seedwt_g3.2t <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 5) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g3.2t)
plot(seedwt_g3.2t)
seedwt_g3.1ot <- gam(seed_wt_mg ~ s(Fire_freq, k = 6) + s(FPC, k= 7) + ti(FPC, by = Fire_freq, bs = 'tp', k = 8), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g3.1ot)
plot(seedwt_g3.1ot)
seedwt_g3.2ot <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 5) + s(FPC, bs = 'cr', k = 4) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g3.2ot)
plot(seedwt_g3.2ot)

seedwt_m3l <- glmer(seed_wt_mg ~ r_fire_freq * r_FPC + (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m3l)
seedwt_g3l <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(FPC) + ti(FPC, by = Fire_freq), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g3l)
plot(seedwt_g3l)
seedwt_g3.1l <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(FPC) + ti(FPC, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g3.1l)
plot(seedwt_g3.1l)
seedwt_g3.2l <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 3) + s(FPC, bs = 'cr', k = 5) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g3.2l)
plot(seedwt_g3.2l)
seedwt_g3.1ol <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(FPC,k = 7) + ti(FPC, by = Fire_freq, bs = 'tp', k = 8), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g3.1ol)
plot(seedwt_g3.1ol)
seedwt_g3.2ol <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 3) + s(FPC, bs = 'cr', k = 3) + ti(FPC, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g3.2ol)
plot(seedwt_g3.2ol)


seedwt_m4t <- glmer(seed_wt_mg ~ r_fire_freq * r_Precip + (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m4t)
seedwt_g4t <- gam(seed_wt_mg ~ s(Fire_freq, k = 7) + s(Precip) + ti(Precip, by = Fire_freq), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g4t)
plot(seedwt_g4t)
seedwt_g4.1t <- gam(seed_wt_mg ~ s(Fire_freq, k = 7) + s(Precip) + ti(Precip, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g4.1t)
plot(seedwt_g4.1t)
seedwt_g4.2t <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 5) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g4.2t)
plot(seedwt_g4.2t)
seedwt_g4.1ot <- gam(seed_wt_mg ~ s(Fire_freq, k = 4) + s(Precip, k = 5) + ti(Precip, by = Fire_freq, bs = 'tp', k = 8), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g4.1ot)
plot(seedwt_g4.1ot)
seedwt_g4.2ot <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 4) + s(Precip, bs = 'cr', k = 4) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g4.2ot)
plot(seedwt_g4.2ot)


seedwt_m4l <- glmer(seed_wt_mg ~ r_fire_freq * r_Precip + (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m4l)
seedwt_g4l <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(Precip, k = 9) + ti(Precip, by = Fire_freq), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g4l)
plot(seedwt_g4l)
seedwt_g4.1l <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(Precip, k = 9) + ti(Precip, by = Fire_freq, bs = 'tp', k = 9), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g4.1l)
plot(seedwt_g4.1l)
seedwt_g4.2l <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g4.2l)
plot(seedwt_g4.2l)
seedwt_g4.1ol <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(Precip, k = 5) + ti(Precip, by = Fire_freq, bs = 'tp', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g4.1ol)
plot(seedwt_g4.1ol)
seedwt_g4.2ol <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 3) + s(Precip, bs = 'cr', k = 5) + ti(Precip, by = Fire_freq, bs = 'cr', k = 3), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g4.2ol)
plot(seedwt_g4.2ol)
seedwt_g4.3l <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 3) + s(Precip, bs = 'cc', k = 5) + ti(Precip, by = Fire_freq, bs = 'cc', k = 4), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g4.3l)
plot(seedwt_g4.3l)

seedwt_m5t <- glmer(seed_wt_mg ~ r_fire_freq * r_Temp + (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m5t)
seedwt_g5t <- gam(seed_wt_mg ~ s(Fire_freq, k = 7) + s(Temp) + ti(Temp, by = Fire_freq), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g5t)
plot(seedwt_g5t)
seedwt_g5.1t <- gam(seed_wt_mg ~ s(Fire_freq, k = 7) + s(Temp) + ti(Temp, by = Fire_freq, bs = 'tp', k = 10), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g5.1t)
plot(seedwt_g5.1t)
seedwt_g5.2t <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 5) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g5.2t)
plot(seedwt_g5.2t)
seedwt_g5.1ot <- gam(seed_wt_mg ~ s(Fire_freq, k = 5) + s(Temp, k= 4) + ti(Temp, by = Fire_freq, bs = 'tp', k = 7), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g5.1ot)
plot(seedwt_g5.1ot)
seedwt_g5.2ot <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 4) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = Fire_freq, bs = 'cr', k = 4), random = ~ (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g5.2ot)
plot(seedwt_g5.2ot)


seedwt_m5l <- glmer(seed_wt_mg ~ r_fire_freq * r_Temp + (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m5l)
seedwt_g5l <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(Temp, k = 9) + ti(Temp, by = Fire_freq), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g5l)
plot(seedwt_g5l)
seedwt_g5.1l <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(Temp, k = 9) + ti(Temp, by = Fire_freq, bs = 'tp', k = 9), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g5.1l)
plot(seedwt_g5.1l)
seedwt_g5.2l <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 3) + s(Temp, bs = 'cr', k = 5) + ti(Temp, by = Fire_freq, bs = 'cr', k = 5), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g5.2l)
plot(seedwt_g5.2l)
seedwt_g5.1ol <- gam(seed_wt_mg ~ s(Fire_freq, k = 3) + s(Temp, k = 5) + ti(Temp, by = Fire_freq, bs = 'tp', k = 4), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g5.1ol)
plot(seedwt_g5.1ol)
seedwt_g5.2ol <- gam(seed_wt_mg ~ s(Fire_freq, bs = 'cr', k = 3) + s(Temp, bs = 'cr', k = 3) + ti(Temp, by = Fire_freq, bs = 'cr', k = 3), random = ~ (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop, method = 'ML')
par(mfrow = c(2,2)); gam.check(seedwt_g5.2ol)
plot(seedwt_g5.2ol)

# Compare AICs for each modelling type
# GLMER
seedwt_aic_t <- list(seedwt_t_null, seedwt_m1t, seedwt_m2t, seedwt_m3t, seedwt_m4t, seedwt_m5t)
aictab(seedwt_aic_t) # Model 2 is best

seedwt_aic_l <- list(seedwt_l_null, seedwt_m1l, seedwt_m2l, seedwt_m3l, seedwt_m4l, seedwt_m5l)
aictab(seedwt_aic_l) # Model 5 is best

# GAM defaults
seedwt_td_aic <- as.data.frame(1:6)
seedwt_td_aic$AICc <- "NA"
seedwt_td_aic$Model <- "NA"
seedwt_td_aic$LL <- "NA"
seedwt_td_aic$AICc[1] <- AICc(seedwt_t_gnull)
seedwt_td_aic$Model[1] <- "Null"
seedwt_td_aic$LL[1] <- logLik(seedwt_t_gnull)
seedwt_td_aic$AICc[2] <- AICc(seedwt_g1t)
seedwt_td_aic$Model[2] <- "m1"
seedwt_td_aic$LL[2] <- logLik(seedwt_g1t)
seedwt_td_aic$AICc[3] <- AICc(seedwt_g2t)
seedwt_td_aic$Model[3] <- "m2"
seedwt_td_aic$LL[3] <- logLik(seedwt_g2t)
seedwt_td_aic$AICc[4] <- AICc(seedwt_g3t)
seedwt_td_aic$Model[4] <- "m3"
seedwt_td_aic$LL[4] <- logLik(seedwt_g3t)
seedwt_td_aic$AICc[5] <- AICc(seedwt_g4t)
seedwt_td_aic$Model[5] <- "m4"
seedwt_td_aic$LL[5] <- logLik(seedwt_g4t)
seedwt_td_aic$AICc[6] <- AICc(seedwt_g5t)
seedwt_td_aic$Model[6] <- "m5"
seedwt_td_aic$LL[6] <- logLik(seedwt_g5t)
seedwt_td_aic$AICc <- as.numeric(seedwt_td_aic$AICc)
seedwt_td_aic$LL <- as.numeric(seedwt_td_aic$LL)
seedwt_td_aic <- seedwt_td_aic[order(seedwt_td_aic$AICc), ]
seedwt_td_aic$Delta_AICc <- "0.00"
seedwt_td_aic$Delta_AICc[2] <- round(seedwt_td_aic$AICc[1]-seedwt_td_aic$AICc[2], 2)
seedwt_td_aic$Delta_AICc[3] <- round(seedwt_td_aic$AICc[1]-seedwt_td_aic$AICc[3], 2)
seedwt_td_aic$Delta_AICc[4] <- round(seedwt_td_aic$AICc[1]-seedwt_td_aic$AICc[4], 2)
seedwt_td_aic$Delta_AICc[5] <- round(seedwt_td_aic$AICc[1]-seedwt_td_aic$AICc[5], 2)
seedwt_td_aic$Delta_AICc[6] <- round(seedwt_td_aic$AICc[1]-seedwt_td_aic$AICc[6], 2)
seedwt_td_aic  # Model 3 is best

seedwt_ld_aic <- as.data.frame(1:6)
seedwt_ld_aic$AICc <- "NA"
seedwt_ld_aic$Model <- "NA"
seedwt_ld_aic$LL <- "NA"
seedwt_ld_aic$AICc[1] <- AICc(seedwt_l_gnull)
seedwt_ld_aic$Model[1] <- "Null"
seedwt_ld_aic$LL[1] <- logLik(seedwt_l_gnull)
seedwt_ld_aic$AICc[2] <- AICc(seedwt_g1l)
seedwt_ld_aic$Model[2] <- "m1"
seedwt_ld_aic$LL[2] <- logLik(seedwt_g1l)
seedwt_ld_aic$AICc[3] <- AICc(seedwt_g2l)
seedwt_ld_aic$Model[3] <- "m2"
seedwt_ld_aic$LL[3] <- logLik(seedwt_g2l)
seedwt_ld_aic$AICc[4] <- AICc(seedwt_g3l)
seedwt_ld_aic$Model[4] <- "m3"
seedwt_ld_aic$LL[4] <- logLik(seedwt_g3l)
seedwt_ld_aic$AICc[5] <- AICc(seedwt_g4l)
seedwt_ld_aic$Model[5] <- "m4"
seedwt_ld_aic$LL[5] <- logLik(seedwt_g4l)
seedwt_ld_aic$AICc[6] <- AICc(seedwt_g5l)
seedwt_ld_aic$Model[6] <- "m5"
seedwt_ld_aic$LL[6] <- logLik(seedwt_g5l)
seedwt_ld_aic$AICc <- as.numeric(seedwt_ld_aic$AICc)
seedwt_ld_aic$LL <- as.numeric(seedwt_ld_aic$LL)
seedwt_ld_aic <- seedwt_ld_aic[order(seedwt_ld_aic$AICc), ]
seedwt_ld_aic$Delta_AICc <- "0.00"
seedwt_ld_aic$Delta_AICc[2] <- round(seedwt_ld_aic$AICc[1]-seedwt_ld_aic$AICc[2], 2)
seedwt_ld_aic$Delta_AICc[3] <- round(seedwt_ld_aic$AICc[1]-seedwt_ld_aic$AICc[3], 2)
seedwt_ld_aic$Delta_AICc[4] <- round(seedwt_ld_aic$AICc[1]-seedwt_ld_aic$AICc[4], 2)
seedwt_ld_aic$Delta_AICc[5] <- round(seedwt_ld_aic$AICc[1]-seedwt_ld_aic$AICc[5], 2)
seedwt_ld_aic$Delta_AICc[6] <- round(seedwt_ld_aic$AICc[1]-seedwt_ld_aic$AICc[6], 2)
seedwt_ld_aic  # Model 2 is best



# GAM s(bs = 'tp')
seedwt_ts_aic <- as.data.frame(1:6)
seedwt_ts_aic$AICc <- "NA"
seedwt_ts_aic$Model <- "NA"
seedwt_ts_aic$LL <- "NA"
seedwt_ts_aic$AICc[1] <- AICc(seedwt_t_gnull)
seedwt_ts_aic$Model[1] <- "Null"
seedwt_ts_aic$LL[1] <- logLik(seedwt_t_gnull)
seedwt_ts_aic$AICc[2] <- AICc(seedwt_g1t)
seedwt_ts_aic$Model[2] <- "m1"
seedwt_ts_aic$LL[2] <- logLik(seedwt_g1t)
seedwt_ts_aic$AICc[3] <- AICc(seedwt_g2.1t)
seedwt_ts_aic$Model[3] <- "m2"
seedwt_ts_aic$LL[3] <- logLik(seedwt_g2.1t)
seedwt_ts_aic$AICc[4] <- AICc(seedwt_g3.1t)
seedwt_ts_aic$Model[4] <- "m3"
seedwt_ts_aic$LL[4] <- logLik(seedwt_g3.1t)
seedwt_ts_aic$AICc[5] <- AICc(seedwt_g4.1t)
seedwt_ts_aic$Model[5] <- "m4"
seedwt_ts_aic$LL[5] <- logLik(seedwt_g4.1t)
seedwt_ts_aic$AICc[6] <- AICc(seedwt_g5.1t)
seedwt_ts_aic$Model[6] <- "m5"
seedwt_ts_aic$LL[6] <- logLik(seedwt_g5.1t)
seedwt_ts_aic$AICc <- as.numeric(seedwt_ts_aic$AICc)
seedwt_ts_aic$LL <- as.numeric(seedwt_ts_aic$LL)
seedwt_ts_aic <- seedwt_ts_aic[order(seedwt_ts_aic$AICc), ]
seedwt_ts_aic$Delta_AICc <- "0.00"
seedwt_ts_aic$Delta_AICc[2] <- round(seedwt_ts_aic$AICc[1]-seedwt_ts_aic$AICc[2], 2)
seedwt_ts_aic$Delta_AICc[3] <- round(seedwt_ts_aic$AICc[1]-seedwt_ts_aic$AICc[3], 2)
seedwt_ts_aic$Delta_AICc[4] <- round(seedwt_ts_aic$AICc[1]-seedwt_ts_aic$AICc[4], 2)
seedwt_ts_aic$Delta_AICc[5] <- round(seedwt_ts_aic$AICc[1]-seedwt_ts_aic$AICc[5], 2)
seedwt_ts_aic$Delta_AICc[6] <- round(seedwt_ts_aic$AICc[1]-seedwt_ts_aic$AICc[6], 2)
seedwt_ts_aic  # Model 3 is best

seedwt_ls_aic <- as.data.frame(1:6)
seedwt_ls_aic$AICc <- "NA"
seedwt_ls_aic$Model <- "NA"
seedwt_ls_aic$LL <- "NA"
seedwt_ls_aic$AICc[1] <- AICc(seedwt_l_gnull)
seedwt_ls_aic$Model[1] <- "Null"
seedwt_ls_aic$LL[1] <- logLik(seedwt_l_gnull)
seedwt_ls_aic$AICc[2] <- AICc(seedwt_g1l)
seedwt_ls_aic$Model[2] <- "m1"
seedwt_ls_aic$LL[2] <- logLik(seedwt_g1l)
seedwt_ls_aic$AICc[3] <- AICc(seedwt_g2.1l)
seedwt_ls_aic$Model[3] <- "m2"
seedwt_ls_aic$LL[3] <- logLik(seedwt_g2.1l)
seedwt_ls_aic$AICc[4] <- AICc(seedwt_g3.1l)
seedwt_ls_aic$Model[4] <- "m3"
seedwt_ls_aic$LL[4] <- logLik(seedwt_g3.1l)
seedwt_ls_aic$AICc[5] <- AICc(seedwt_g4.1l)
seedwt_ls_aic$Model[5] <- "m4"
seedwt_ls_aic$LL[5] <- logLik(seedwt_g4.1l)
seedwt_ls_aic$AICc[6] <- AICc(seedwt_g5.1l)
seedwt_ls_aic$Model[6] <- "m5"
seedwt_ls_aic$LL[6] <- logLik(seedwt_g5.1l)
seedwt_ls_aic$AICc <- as.numeric(seedwt_ls_aic$AICc)
seedwt_ls_aic$LL <- as.numeric(seedwt_ls_aic$LL)
seedwt_ls_aic <- seedwt_ls_aic[order(seedwt_ls_aic$AICc), ]
seedwt_ls_aic$Delta_AICc <- "0.00"
seedwt_ls_aic$Delta_AICc[2] <- round(seedwt_ls_aic$AICc[1]-seedwt_ls_aic$AICc[2], 2)
seedwt_ls_aic$Delta_AICc[3] <- round(seedwt_ls_aic$AICc[1]-seedwt_ls_aic$AICc[3], 2)
seedwt_ls_aic$Delta_AICc[4] <- round(seedwt_ls_aic$AICc[1]-seedwt_ls_aic$AICc[4], 2)
seedwt_ls_aic$Delta_AICc[5] <- round(seedwt_ls_aic$AICc[1]-seedwt_ls_aic$AICc[5], 2)
seedwt_ls_aic$Delta_AICc[6] <- round(seedwt_ls_aic$AICc[1]-seedwt_ls_aic$AICc[6], 2)
seedwt_ls_aic  # Model 2 is best

# GAM s(bs = 'tp') optimised
seedwt_tso_aic <- as.data.frame(1:6)
seedwt_tso_aic$AICc <- "NA"
seedwt_tso_aic$Model <- "NA"
seedwt_tso_aic$LL <- "NA"
seedwt_tso_aic$AICc[1] <- AICc(seedwt_t_gnull)
seedwt_tso_aic$Model[1] <- "Null"
seedwt_tso_aic$LL[1] <- logLik(seedwt_t_gnull)
seedwt_tso_aic$AICc[2] <- AICc(seedwt_g1t)
seedwt_tso_aic$Model[2] <- "m1"
seedwt_tso_aic$LL[2] <- logLik(seedwt_g1t)
seedwt_tso_aic$AICc[3] <- AICc(seedwt_g2.1ot)
seedwt_tso_aic$Model[3] <- "m2"
seedwt_tso_aic$LL[3] <- logLik(seedwt_g2.1ot)
seedwt_tso_aic$AICc[4] <- AICc(seedwt_g3.1ot)
seedwt_tso_aic$Model[4] <- "m3"
seedwt_tso_aic$LL[4] <- logLik(seedwt_g3.1ot)
seedwt_tso_aic$AICc[5] <- AICc(seedwt_g4.1ot)
seedwt_tso_aic$Model[5] <- "m4"
seedwt_tso_aic$LL[5] <- logLik(seedwt_g4.1ot)
seedwt_tso_aic$AICc[6] <- AICc(seedwt_g5.1ot)
seedwt_tso_aic$Model[6] <- "m5"
seedwt_tso_aic$LL[6] <- logLik(seedwt_g5.1ot)
seedwt_tso_aic$AICc <- as.numeric(seedwt_tso_aic$AICc)
seedwt_tso_aic$LL <- as.numeric(seedwt_tso_aic$LL)
seedwt_tso_aic <- seedwt_tso_aic[order(seedwt_tso_aic$AICc), ]
seedwt_tso_aic$Delta_AICc <- "0.00"
seedwt_tso_aic$Delta_AICc[2] <- round(seedwt_tso_aic$AICc[1]-seedwt_tso_aic$AICc[2], 2)
seedwt_tso_aic$Delta_AICc[3] <- round(seedwt_tso_aic$AICc[1]-seedwt_tso_aic$AICc[3], 2)
seedwt_tso_aic$Delta_AICc[4] <- round(seedwt_tso_aic$AICc[1]-seedwt_tso_aic$AICc[4], 2)
seedwt_tso_aic$Delta_AICc[5] <- round(seedwt_tso_aic$AICc[1]-seedwt_tso_aic$AICc[5], 2)
seedwt_tso_aic$Delta_AICc[6] <- round(seedwt_tso_aic$AICc[1]-seedwt_tso_aic$AICc[6], 2)
seedwt_tso_aic  # Model 3 is best

seedwt_lso_aic <- as.data.frame(1:6)
seedwt_lso_aic$AICc <- "NA"
seedwt_lso_aic$Model <- "NA"
seedwt_lso_aic$LL <- "NA"
seedwt_lso_aic$AICc[1] <- AICc(seedwt_l_gnull)
seedwt_lso_aic$Model[1] <- "Null"
seedwt_lso_aic$LL[1] <- logLik(seedwt_l_gnull)
seedwt_lso_aic$AICc[2] <- AICc(seedwt_g1l)
seedwt_lso_aic$Model[2] <- "m1"
seedwt_lso_aic$LL[2] <- logLik(seedwt_g1l)
seedwt_lso_aic$AICc[3] <- AICc(seedwt_g2.1ol)
seedwt_lso_aic$Model[3] <- "m2"
seedwt_lso_aic$LL[3] <- logLik(seedwt_g2.1ol)
seedwt_lso_aic$AICc[4] <- AICc(seedwt_g3.1ol)
seedwt_lso_aic$Model[4] <- "m3"
seedwt_lso_aic$LL[4] <- logLik(seedwt_g3.1ol)
seedwt_lso_aic$AICc[5] <- AICc(seedwt_g4.1ol)
seedwt_lso_aic$Model[5] <- "m4"
seedwt_lso_aic$LL[5] <- logLik(seedwt_g4.1ol)
seedwt_lso_aic$AICc[6] <- AICc(seedwt_g5.1ol)
seedwt_lso_aic$Model[6] <- "m5"
seedwt_lso_aic$LL[6] <- logLik(seedwt_g5.1ol)
seedwt_lso_aic$AICc <- as.numeric(seedwt_lso_aic$AICc)
seedwt_lso_aic$LL <- as.numeric(seedwt_lso_aic$LL)
seedwt_lso_aic <- seedwt_lso_aic[order(seedwt_lso_aic$AICc), ]
seedwt_lso_aic$Delta_AICc <- "0.00"
seedwt_lso_aic$Delta_AICc[2] <- round(seedwt_lso_aic$AICc[1]-seedwt_lso_aic$AICc[2], 2)
seedwt_lso_aic$Delta_AICc[3] <- round(seedwt_lso_aic$AICc[1]-seedwt_lso_aic$AICc[3], 2)
seedwt_lso_aic$Delta_AICc[4] <- round(seedwt_lso_aic$AICc[1]-seedwt_lso_aic$AICc[4], 2)
seedwt_lso_aic$Delta_AICc[5] <- round(seedwt_lso_aic$AICc[1]-seedwt_lso_aic$AICc[5], 2)
seedwt_lso_aic$Delta_AICc[6] <- round(seedwt_lso_aic$AICc[1]-seedwt_lso_aic$AICc[6], 2)
seedwt_lso_aic  # Model 4 is best



# GAM ti(bs='cr')
seedwt_tsc_aic <- as.data.frame(1:6)
seedwt_tsc_aic$AICc <- "NA"
seedwt_tsc_aic$Model <- "NA"
seedwt_tsc_aic$LL <- "NA"
seedwt_tsc_aic$AICc[1] <- AICc(seedwt_t_gnull)
seedwt_tsc_aic$Model[1] <- "Null"
seedwt_tsc_aic$LL[1] <- logLik(seedwt_t_gnull)
seedwt_tsc_aic$AICc[2] <- AICc(seedwt_g1.1t)
seedwt_tsc_aic$Model[2] <- "m1"
seedwt_tsc_aic$LL[2] <- logLik(seedwt_g1.1t)
seedwt_tsc_aic$AICc[3] <- AICc(seedwt_g2.2t)
seedwt_tsc_aic$Model[3] <- "m2"
seedwt_tsc_aic$LL[3] <- logLik(seedwt_g2.2t)
seedwt_tsc_aic$AICc[4] <- AICc(seedwt_g3.2t)
seedwt_tsc_aic$Model[4] <- "m3"
seedwt_tsc_aic$LL[4] <- logLik(seedwt_g3.2t)
seedwt_tsc_aic$AICc[5] <- AICc(seedwt_g4.2t)
seedwt_tsc_aic$Model[6] <- "m4"
seedwt_tsc_aic$LL[5] <- logLik(seedwt_g4.2t)
seedwt_tsc_aic$AICc[6] <- AICc(seedwt_g5.2t)
seedwt_tsc_aic$Model[6] <- "m5"
seedwt_tsc_aic$LL[6] <- logLik(seedwt_g5.2t)
seedwt_tsc_aic$AICc <- as.numeric(seedwt_tsc_aic$AICc)
seedwt_tsc_aic$LL <- as.numeric(seedwt_tsc_aic$LL)
seedwt_tsc_aic <- seedwt_tsc_aic[order(seedwt_tsc_aic$AICc), ]
seedwt_tsc_aic$Delta_AICc <- "0.00"
seedwt_tsc_aic$Delta_AICc[2] <- round(seedwt_tsc_aic$AICc[1]-seedwt_tsc_aic$AICc[2], 2)
seedwt_tsc_aic$Delta_AICc[3] <- round(seedwt_tsc_aic$AICc[1]-seedwt_tsc_aic$AICc[3], 2)
seedwt_tsc_aic$Delta_AICc[4] <- round(seedwt_tsc_aic$AICc[1]-seedwt_tsc_aic$AICc[4], 2)
seedwt_tsc_aic$Delta_AICc[5] <- round(seedwt_tsc_aic$AICc[1]-seedwt_tsc_aic$AICc[5], 2)
seedwt_tsc_aic$Delta_AICc[6] <- round(seedwt_tsc_aic$AICc[1]-seedwt_tsc_aic$AICc[6], 2)
seedwt_tsc_aic  # Model 3 is best

seedwt_lsc_aic <- as.data.frame(1:6)
seedwt_lsc_aic$AICc <- "NA"
seedwt_lsc_aic$Model <- "NA"
seedwt_lsc_aic$LL <- "NA"
seedwt_lsc_aic$AICc[1] <- AICc(seedwt_l_gnull)
seedwt_lsc_aic$Model[1] <- "Null"
seedwt_lsc_aic$LL[1] <- logLik(seedwt_l_gnull)
seedwt_lsc_aic$AICc[2] <- AICc(seedwt_g1.1l)
seedwt_lsc_aic$Model[2] <- "m1"
seedwt_lsc_aic$LL[2] <- logLik(seedwt_g1.1l)
seedwt_lsc_aic$AICc[3] <- AICc(seedwt_g2.2l)
seedwt_lsc_aic$Model[3] <- "m2"
seedwt_lsc_aic$LL[3] <- logLik(seedwt_g2.2l)
seedwt_lsc_aic$AICc[4] <- AICc(seedwt_g3.2l)
seedwt_lsc_aic$Model[4] <- "m3"
seedwt_lsc_aic$LL[4] <- logLik(seedwt_g3.2l)
seedwt_lsc_aic$AICc[5] <- AICc(seedwt_g4.2l)
seedwt_lsc_aic$Model[5] <- "m4"
seedwt_lsc_aic$LL[5] <- logLik(seedwt_g4.2l)
seedwt_lsc_aic$AICc[6] <- AICc(seedwt_g5.2l)
seedwt_lsc_aic$Model[6] <- "m5"
seedwt_lsc_aic$LL[6] <- logLik(seedwt_g5.2l)
seedwt_lsc_aic$AICc <- as.numeric(seedwt_lsc_aic$AICc)
seedwt_lsc_aic$LL <- as.numeric(seedwt_lsc_aic$LL)
seedwt_lsc_aic <- seedwt_lsc_aic[order(seedwt_lsc_aic$AICc), ]
seedwt_lsc_aic$Delta_AICc <- "0.00"
seedwt_lsc_aic$Delta_AICc[2] <- round(seedwt_lsc_aic$AICc[1]-seedwt_lsc_aic$AICc[2], 2)
seedwt_lsc_aic$Delta_AICc[3] <- round(seedwt_lsc_aic$AICc[1]-seedwt_lsc_aic$AICc[3], 2)
seedwt_lsc_aic$Delta_AICc[4] <- round(seedwt_lsc_aic$AICc[1]-seedwt_lsc_aic$AICc[4], 2)
seedwt_lsc_aic$Delta_AICc[5] <- round(seedwt_lsc_aic$AICc[1]-seedwt_lsc_aic$AICc[5], 2)
seedwt_lsc_aic$Delta_AICc[6] <- round(seedwt_lsc_aic$AICc[1]-seedwt_lsc_aic$AICc[6], 2)
seedwt_lsc_aic  # Model 4 is best


# GAM ti(bs='cr') optimised
seedwt_tsco_aic <- as.data.frame(1:6)
seedwt_tsco_aic$AICc <- "NA"
seedwt_tsco_aic$Model <- "NA"
seedwt_tsco_aic$LL <- "NA"
seedwt_tsco_aic$AICc[1] <- AICc(seedwt_t_gnull)
seedwt_tsco_aic$Model[1] <- "Null"
seedwt_tsco_aic$LL[1] <- logLik(seedwt_t_gnull)
seedwt_tsco_aic$AICc[2] <- AICc(seedwt_g1.1t)
seedwt_tsco_aic$Model[2] <- "m1"
seedwt_tsco_aic$LL[2] <- logLik(seedwt_g1.1ot)
seedwt_tsco_aic$AICc[3] <- AICc(seedwt_g2.2ot)
seedwt_tsco_aic$Model[3] <- "m2"
seedwt_tsco_aic$LL[3] <- logLik(seedwt_g2.2ot)
seedwt_tsco_aic$AICc[4] <- AICc(seedwt_g3.2ot)
seedwt_tsco_aic$Model[4] <- "m3"
seedwt_tsco_aic$LL[4] <- logLik(seedwt_g3.2ot)
seedwt_tsco_aic$AICc[5] <- AICc(seedwt_g4.2ot)
seedwt_tsco_aic$Model[5] <- "m4"
seedwt_tsco_aic$LL[5] <- logLik(seedwt_g4.2ot)
seedwt_tsco_aic$AICc[6] <- AICc(seedwt_g5.2ot)
seedwt_tsco_aic$Model[6] <- "m5"
seedwt_tsco_aic$LL[6] <- logLik(seedwt_g5.2ot)
seedwt_tsco_aic$AICc <- as.numeric(seedwt_tsco_aic$AICc)
seedwt_tsco_aic$LL <- as.numeric(seedwt_tsco_aic$LL)
seedwt_tsco_aic <- seedwt_tsco_aic[order(seedwt_tsco_aic$AICc), ]
seedwt_tsco_aic$Delta_AICc <- "0.00"
seedwt_tsco_aic$Delta_AICc[2] <- round(seedwt_tsco_aic$AICc[1]-seedwt_tsco_aic$AICc[2], 2)
seedwt_tsco_aic$Delta_AICc[3] <- round(seedwt_tsco_aic$AICc[1]-seedwt_tsco_aic$AICc[3], 2)
seedwt_tsco_aic$Delta_AICc[4] <- round(seedwt_tsco_aic$AICc[1]-seedwt_tsco_aic$AICc[4], 2)
seedwt_tsco_aic$Delta_AICc[5] <- round(seedwt_tsco_aic$AICc[1]-seedwt_tsco_aic$AICc[5], 2)
seedwt_tsco_aic$Delta_AICc[6] <- round(seedwt_tsco_aic$AICc[1]-seedwt_tsco_aic$AICc[6], 2)
seedwt_tsco_aic  # Model 3 is best

seedwt_lsco_aic <- as.data.frame(1:6)
seedwt_lsco_aic$AICc <- "NA"
seedwt_lsco_aic$Model <- "NA"
seedwt_lsco_aic$LL <- "NA"
seedwt_lsco_aic$AICc[1] <- AICc(seedwt_l_gnull)
seedwt_lsco_aic$Model[1] <- "Null"
seedwt_lsco_aic$LL[1] <- logLik(seedwt_l_gnull)
seedwt_lsco_aic$AICc[2] <- AICc(seedwt_g1.1l)
seedwt_lsco_aic$Model[2] <- "m1"
seedwt_lsco_aic$LL[2] <- logLik(seedwt_g1.1l)
seedwt_lsco_aic$AICc[3] <- AICc(seedwt_g2.2ol)
seedwt_lsco_aic$Model[3] <- "m2"
seedwt_lsco_aic$LL[3] <- logLik(seedwt_g2.2ol)
seedwt_lsco_aic$AICc[4] <- AICc(seedwt_g3.2ol)
seedwt_lsco_aic$Model[4] <- "m3"
seedwt_lsco_aic$LL[4] <- logLik(seedwt_g3.2ol)
seedwt_lsco_aic$AICc[5] <- AICc(seedwt_g4.2ol)
seedwt_lsco_aic$Model[5] <- "m4"
seedwt_lsco_aic$LL[5] <- logLik(seedwt_g4.2ol)
seedwt_lsco_aic$AICc[6] <- AICc(seedwt_g5.2ol)
seedwt_lsco_aic$Model[6] <- "m5"
seedwt_lsco_aic$LL[6] <- logLik(seedwt_g5.2ol)
seedwt_lsco_aic$AICc <- as.numeric(seedwt_lsco_aic$AICc)
seedwt_lsco_aic$LL <- as.numeric(seedwt_lsco_aic$LL)
seedwt_lsco_aic <- seedwt_lsco_aic[order(seedwt_lsco_aic$AICc), ]
seedwt_lsco_aic$Delta_AICc <- "0.00"
seedwt_lsco_aic$Delta_AICc[2] <- round(seedwt_lsco_aic$AICc[1]-seedwt_lsco_aic$AICc[2], 2)
seedwt_lsco_aic$Delta_AICc[3] <- round(seedwt_lsco_aic$AICc[1]-seedwt_lsco_aic$AICc[3], 2)
seedwt_lsco_aic$Delta_AICc[4] <- round(seedwt_lsco_aic$AICc[1]-seedwt_lsco_aic$AICc[4], 2)
seedwt_lsco_aic$Delta_AICc[5] <- round(seedwt_lsco_aic$AICc[1]-seedwt_lsco_aic$AICc[5], 2)
seedwt_lsco_aic$Delta_AICc[6] <- round(seedwt_lsco_aic$AICc[1]-seedwt_lsco_aic$AICc[6], 2)
seedwt_lsco_aic  # Model 4 is best




# Predict to new data for seed size ----
# GLMER torulosa m2 littoralis m5
# GAM torulosa m3
# GAM littoralis m2 but m4 for bs = tp optimised, cr and cr optimised.

tseed_m3_l <- data.frame(Fire_freq = min(tor_cum.prop$Fire_freq),
                         FPC = seq(min(tor_cum.prop$FPC, na.rm = T), max(tor_cum.prop$FPC, na.rm = T), length = 50))
tseed_m3_a <- data.frame(Fire_freq = mean(tor_cum.prop$Fire_freq),
                         FPC = seq(min(tor_cum.prop$FPC, na.rm = T), max(tor_cum.prop$FPC, na.rm = T), length = 50))
tseed_m3_h <- data.frame(Fire_freq = max(tor_cum.prop$Fire_freq),
                         FPC = seq(min(tor_cum.prop$FPC, na.rm = T), max(tor_cum.prop$FPC, na.rm = T), length = 50))

lseed_m2_l <- data.frame(Fire_freq = min(lit_cum.prop$Fire_freq),
                         Latitude = seq(min(lit_cum.prop$Latitude), max(lit_cum.prop$Latitude), length = 50))
lseed_m2_a <- data.frame(Fire_freq = mean(lit_cum.prop$Fire_freq),
                         Latitude = seq(min(lit_cum.prop$Latitude), max(lit_cum.prop$Latitude), length = 50))
lseed_m2_h <- data.frame(Fire_freq = max(lit_cum.prop$Fire_freq),
                         Latitude = seq(min(lit_cum.prop$Latitude), max(lit_cum.prop$Latitude), length = 50))


lseed_m4_l <- data.frame(Fire_freq = min(lit_cum.prop$Fire_freq),
                         Precip = seq(min(lit_cum.prop$Precip), max(lit_cum.prop$Precip), length = 50))
lseed_m4_a <- data.frame(Fire_freq = mean(lit_cum.prop$Fire_freq),
                         Precip = seq(min(lit_cum.prop$Precip), max(lit_cum.prop$Precip), length = 50))
lseed_m4_h <- data.frame(Fire_freq = max(lit_cum.prop$Fire_freq),
                         Precip = seq(min(lit_cum.prop$Precip), max(lit_cum.prop$Precip), length = 50))


# GLMER
tseed_m2lm_l <- data.frame(r_fire_freq = min(tor_cum.prop$r_fire_freq),
                           r_Latitude = seq(min(tor_cum.prop$r_Latitude), max(tor_cum.prop$r_Latitude), length = 50))
tseed_m2lm_a <- data.frame(r_fire_freq = mean(tor_cum.prop$r_fire_freq),
                           r_Latitude = seq(min(tor_cum.prop$r_Latitude), max(tor_cum.prop$r_Latitude), length = 50))
tseed_m2lm_h <- data.frame(r_fire_freq = max(tor_cum.prop$r_fire_freq),
                           r_Latitude = seq(min(tor_cum.prop$r_Latitude), max(tor_cum.prop$r_Latitude), length = 50))

lseed_m5lm_l <- data.frame(r_fire_freq = min(lit_cum.prop$r_fire_freq),
                         r_Temp = seq(min(lit_cum.prop$r_Temp), max(lit_cum.prop$r_Temp), length = 50))
lseed_m5lm_a <- data.frame(r_fire_freq = mean(lit_cum.prop$r_fire_freq),
                         r_Temp = seq(min(lit_cum.prop$r_Temp), max(lit_cum.prop$r_Temp), length = 50))
lseed_m5lm_h <- data.frame(r_fire_freq = max(lit_cum.prop$r_fire_freq),
                         r_Temp = seq(min(lit_cum.prop$r_Temp), max(lit_cum.prop$r_Temp), length = 50))


ptseed_m2lm_l <- predictSE(seedwt_m2t, newdata = tseed_m2lm_l, se.fit = T, type = 'link')
tseed_m2lm_l$fit.link <- ptseed_m2lm_l$fit
tseed_m2lm_l$se.link <- ptseed_m2lm_l$se.fit
tseed_m2lm_l$lci.link <- tseed_m2lm_l$fit.link - (tseed_m2lm_l$se.link * 1.96)
tseed_m2lm_l$uci.link <- tseed_m2lm_l$fit.link + (tseed_m2lm_l$se.link * 1.96)
tseed_m2lm_l$fit <- exp(tseed_m2lm_l$fit.link)
tseed_m2lm_l$se <- exp(tseed_m2lm_l$se.link)
tseed_m2lm_l$lci <- exp(tseed_m2lm_l$lci.link)
tseed_m2lm_l$uci <- exp(tseed_m2lm_l$uci.link)

ptseed_m2lm_a <- predictSE(seedwt_m2t, newdata = tseed_m2lm_a, se.fit = T, type = 'link')
tseed_m2lm_a$fit.link <- ptseed_m2lm_a$fit
tseed_m2lm_a$se.link <- ptseed_m2lm_a$se.fit
tseed_m2lm_a$lci.link <- tseed_m2lm_a$fit.link - (tseed_m2lm_a$se.link * 1.96)
tseed_m2lm_a$uci.link <- tseed_m2lm_a$fit.link + (tseed_m2lm_a$se.link * 1.96)
tseed_m2lm_a$fit <- exp(tseed_m2lm_a$fit.link)
tseed_m2lm_a$se <- exp(tseed_m2lm_a$se.link)
tseed_m2lm_a$lci <- exp(tseed_m2lm_a$lci.link)
tseed_m2lm_a$uci <- exp(tseed_m2lm_a$uci.link)

ptseed_m2lm_h <- predictSE(seedwt_m2t, newdata = tseed_m2lm_h, se.fit = T, type = 'link')
tseed_m2lm_h$fit.link <- ptseed_m2lm_h$fit
tseed_m2lm_h$se.link <- ptseed_m2lm_h$se.fit
tseed_m2lm_h$lci.link <- tseed_m2lm_h$fit.link - (tseed_m2lm_h$se.link * 1.96)
tseed_m2lm_h$uci.link <- tseed_m2lm_h$fit.link + (tseed_m2lm_h$se.link * 1.96)
tseed_m2lm_h$fit <- exp(tseed_m2lm_h$fit.link)
tseed_m2lm_h$se <- exp(tseed_m2lm_h$se.link)
tseed_m2lm_h$lci <- exp(tseed_m2lm_h$lci.link)
tseed_m2lm_h$uci <- exp(tseed_m2lm_h$uci.link)



plseed_m5lm_l <- predictSE(seedwt_m5t, newdata = lseed_m5lm_l, se.fit = T, type = 'link')
lseed_m5lm_l$fit.link <- plseed_m5lm_l$fit
lseed_m5lm_l$se.link <- plseed_m5lm_l$se.fit
lseed_m5lm_l$lci.link <- lseed_m5lm_l$fit.link - (lseed_m5lm_l$se.link * 1.96)
lseed_m5lm_l$uci.link <- lseed_m5lm_l$fit.link + (lseed_m5lm_l$se.link * 1.96)
lseed_m5lm_l$fit <- exp(lseed_m5lm_l$fit.link)
lseed_m5lm_l$se <- exp(lseed_m5lm_l$se.link)
lseed_m5lm_l$lci <- exp(lseed_m5lm_l$lci.link)
lseed_m5lm_l$uci <- exp(lseed_m5lm_l$uci.link)

plseed_m5lm_a <- predictSE(seedwt_m5t, newdata = lseed_m5lm_a, se.fit = T, type = 'link')
lseed_m5lm_a$fit.link <- plseed_m5lm_a$fit
lseed_m5lm_a$se.link <- plseed_m5lm_a$se.fit
lseed_m5lm_a$lci.link <- lseed_m5lm_a$fit.link - (lseed_m5lm_a$se.link * 1.96)
lseed_m5lm_a$uci.link <- lseed_m5lm_a$fit.link + (lseed_m5lm_a$se.link * 1.96)
lseed_m5lm_a$fit <- exp(lseed_m5lm_a$fit.link)
lseed_m5lm_a$se <- exp(lseed_m5lm_a$se.link)
lseed_m5lm_a$lci <- exp(lseed_m5lm_a$lci.link)
lseed_m5lm_a$uci <- exp(lseed_m5lm_a$uci.link)

plseed_m5lm_h <- predictSE(seedwt_m5t, newdata = lseed_m5lm_h, se.fit = T, type = 'link')
lseed_m5lm_h$fit.link <- plseed_m5lm_h$fit
lseed_m5lm_h$se.link <- plseed_m5lm_h$se.fit
lseed_m5lm_h$lci.link <- lseed_m5lm_h$fit.link - (lseed_m5lm_h$se.link * 1.96)
lseed_m5lm_h$uci.link <- lseed_m5lm_h$fit.link + (lseed_m5lm_h$se.link * 1.96)
lseed_m5lm_h$fit <- exp(lseed_m5lm_h$fit.link)
lseed_m5lm_h$se <- exp(lseed_m5lm_h$se.link)
lseed_m5lm_h$lci <- exp(lseed_m5lm_h$lci.link)
lseed_m5lm_h$uci <- exp(lseed_m5lm_h$uci.link)



# GAM default
tseed_m3d_l <- tseed_m3_l
tseed_m3d_a <- tseed_m3_a
tseed_m3d_h <- tseed_m3_h

lseed_m2d_l <- lseed_m2_l
lseed_m2d_a <- lseed_m2_a
lseed_m2d_h <- lseed_m2_h


ptseed_m3d_l <- predict(seedwt_g3t, newdata = tseed_m3d_l, se.fit = T, type = 'response')
tseed_m3d_l$fit <- ptseed_m3d_l$fit
tseed_m3d_l$se <- ptseed_m3d_l$se.fit
tseed_m3d_l$lci <- tseed_m3d_l$fit - (tseed_m3d_l$se* 1.96)
tseed_m3d_l$uci <- tseed_m3d_l$fit + (tseed_m3d_l$se * 1.96)

ptseed_m3d_a <- predict(seedwt_g3t, newdata = tseed_m3d_a, se.fit = T, type = 'response')
tseed_m3d_a$fit <- ptseed_m3d_a$fit
tseed_m3d_a$se <- ptseed_m3d_a$se.fit
tseed_m3d_a$lci <- tseed_m3d_a$fit - (tseed_m3d_a$se* 1.96)
tseed_m3d_a$uci <- tseed_m3d_a$fit + (tseed_m3d_a$se * 1.96)

ptseed_m3d_h <- predict(seedwt_g3t, newdata = tseed_m3d_h, se.fit = T, type = 'response')
tseed_m3d_h$fit <- ptseed_m3d_h$fit
tseed_m3d_h$se <- ptseed_m3d_h$se.fit
tseed_m3d_h$lci <- tseed_m3d_h$fit - (tseed_m3d_h$se* 1.96)
tseed_m3d_h$uci <- tseed_m3d_h$fit + (tseed_m3d_h$se * 1.96)



plseed_m2d_l <- predict(seedwt_g2l, newdata = lseed_m2d_l, se.fit = T, type = 'response')
lseed_m2d_l$fit <- plseed_m2d_l$fit
lseed_m2d_l$se <- plseed_m2d_l$se.fit
lseed_m2d_l$lci <- lseed_m2d_l$fit - (lseed_m2d_l$se * 1.96)
lseed_m2d_l$uci <- lseed_m2d_l$fit + (lseed_m2d_l$se * 1.96)

plseed_m2d_a <- predict(seedwt_g2l, newdata = lseed_m2d_a, se.fit = T, type = 'response')
lseed_m2d_a$fit <- plseed_m2d_a$fit
lseed_m2d_a$se <- plseed_m2d_a$se.fit
lseed_m2d_a$lci <- lseed_m2d_a$fit - (lseed_m2d_a$se * 1.96)
lseed_m2d_a$uci <- lseed_m2d_a$fit + (lseed_m2d_a$se * 1.96)

plseed_m2d_h <- predict(seedwt_g2l, newdata = lseed_m2d_h, se.fit = T, type = 'response')
lseed_m2d_h$fit <- plseed_m2d_h$fit
lseed_m2d_h$se <- plseed_m2d_h$se.fit
lseed_m2d_h$lci <- lseed_m2d_h$fit - (lseed_m2d_h$se * 1.96)
lseed_m2d_h$uci <- lseed_m2d_h$fit + (lseed_m2d_h$se * 1.96)


# GAM bs = tp
tseed_m3t_l <- tseed_m3_l
tseed_m3t_a <- tseed_m3_a
tseed_m3t_h <- tseed_m3_h

lseed_m2t_l <- lseed_m2_l
lseed_m2t_a <- lseed_m2_a
lseed_m2t_h <- lseed_m2_h


ptseed_m3t_l <- predict(seedwt_g3.1t, newdata = tseed_m3t_l, se.fit = T, type = 'response')
tseed_m3t_l$fit <- ptseed_m3t_l$fit
tseed_m3t_l$se <- ptseed_m3t_l$se.fit
tseed_m3t_l$lci <- tseed_m3t_l$fit - (tseed_m3t_l$se* 1.96)
tseed_m3t_l$uci <- tseed_m3t_l$fit + (tseed_m3t_l$se * 1.96)

ptseed_m3t_a <- predict(seedwt_g3.1t, newdata = tseed_m3t_a, se.fit = T, type = 'response')
tseed_m3t_a$fit <- ptseed_m3t_a$fit
tseed_m3t_a$se <- ptseed_m3t_a$se.fit
tseed_m3t_a$lci <- tseed_m3t_a$fit - (tseed_m3t_a$se* 1.96)
tseed_m3t_a$uci <- tseed_m3t_a$fit + (tseed_m3t_a$se * 1.96)

ptseed_m3t_h <- predict(seedwt_g3.1t, newdata = tseed_m3t_h, se.fit = T, type = 'response')
tseed_m3t_h$fit <- ptseed_m3t_h$fit
tseed_m3t_h$se <- ptseed_m3t_h$se.fit
tseed_m3t_h$lci <- tseed_m3t_h$fit - (tseed_m3t_h$se* 1.96)
tseed_m3t_h$uci <- tseed_m3t_h$fit + (tseed_m3t_h$se * 1.96)


plseed_m2t_l <- predict(seedwt_g2.1l, newdata = lseed_m2t_l, se.fit = T, type = 'response')
lseed_m2t_l$fit <- plseed_m2t_l$fit
lseed_m2t_l$se <- plseed_m2t_l$se.fit
lseed_m2t_l$lci <- lseed_m2t_l$fit - (lseed_m2t_l$se * 1.96)
lseed_m2t_l$uci <- lseed_m2t_l$fit + (lseed_m2t_l$se * 1.96)

plseed_m2t_a <- predict(seedwt_g2.1l, newdata = lseed_m2t_a, se.fit = T, type = 'response')
lseed_m2t_a$fit <- plseed_m2t_a$fit
lseed_m2t_a$se <- plseed_m2t_a$se.fit
lseed_m2t_a$lci <- lseed_m2t_a$fit - (lseed_m2t_a$se * 1.96)
lseed_m2t_a$uci <- lseed_m2t_a$fit + (lseed_m2t_a$se * 1.96)

plseed_m2t_h <- predict(seedwt_g2.1l, newdata = lseed_m2t_h, se.fit = T, type = 'response')
lseed_m2t_h$fit <- plseed_m2t_h$fit
lseed_m2t_h$se <- plseed_m2t_h$se.fit
lseed_m2t_h$lci <- lseed_m2t_h$fit - (lseed_m2t_h$se * 1.96)
lseed_m2t_h$uci <- lseed_m2t_h$fit + (lseed_m2t_h$se * 1.96)


# GAM bs = tp optimised
tseed_m3to_l <- tseed_m3_l
tseed_m3to_a <- tseed_m3_a
tseed_m3to_h <- tseed_m3_h


lseed_m4to_l <- lseed_m4_l
lseed_m4to_a <- lseed_m4_a
lseed_m4to_h <- lseed_m4_h

ptseed_m3to_l <- predict(seedwt_g3.1ot, newdata = tseed_m3to_l, se.fit = T, type = 'response')
tseed_m3to_l$fit <- ptseed_m3to_l$fit
tseed_m3to_l$se <- ptseed_m3to_l$se.fit
tseed_m3to_l$lci <- tseed_m3to_l$fit - (tseed_m3to_l$se* 1.96)
tseed_m3to_l$uci <- tseed_m3to_l$fit + (tseed_m3to_l$se * 1.96)

ptseed_m3to_a <- predict(seedwt_g3.1ot, newdata = tseed_m3to_a, se.fit = T, type = 'response')
tseed_m3to_a$fit <- ptseed_m3to_a$fit
tseed_m3to_a$se <- ptseed_m3to_a$se.fit
tseed_m3to_a$lci <- tseed_m3to_a$fit - (tseed_m3to_a$se* 1.96)
tseed_m3to_a$uci <- tseed_m3to_a$fit + (tseed_m3to_a$se * 1.96)

ptseed_m3to_h <- predict(seedwt_g3.1ot, newdata = tseed_m3to_h, se.fit = T, type = 'response')
tseed_m3to_h$fit <- ptseed_m3to_h$fit
tseed_m3to_h$se <- ptseed_m3to_h$se.fit
tseed_m3to_h$lci <- tseed_m3to_h$fit - (tseed_m3to_h$se* 1.96)
tseed_m3to_h$uci <- tseed_m3to_h$fit + (tseed_m3to_h$se * 1.96)



plseed_m4to_l <- predict(seedwt_g4.1ol, newdata = lseed_m4to_l, se.fit = T, type = 'response')
lseed_m4to_l$fit <- plseed_m4to_l$fit
lseed_m4to_l$se <- plseed_m4to_l$se.fit
lseed_m4to_l$lci <- lseed_m4to_l$fit - (lseed_m4to_l$se * 1.96)
lseed_m4to_l$uci <- lseed_m4to_l$fit + (lseed_m4to_l$se * 1.96)

plseed_m4to_a <- predict(seedwt_g4.1ol, newdata = lseed_m4to_a, se.fit = T, type = 'response')
lseed_m4to_a$fit <- plseed_m4to_a$fit
lseed_m4to_a$se <- plseed_m4to_a$se.fit
lseed_m4to_a$lci <- lseed_m4to_a$fit - (lseed_m4to_a$se * 1.96)
lseed_m4to_a$uci <- lseed_m4to_a$fit + (lseed_m4to_a$se * 1.96)

plseed_m4to_h <- predict(seedwt_g4.1ol, newdata = lseed_m4to_h, se.fit = T, type = 'response')
lseed_m4to_h$fit <- plseed_m4to_h$fit
lseed_m4to_h$se <- plseed_m4to_h$se.fit
lseed_m4to_h$lci <- lseed_m4to_h$fit - (lseed_m4to_h$se * 1.96)
lseed_m4to_h$uci <- lseed_m4to_h$fit + (lseed_m4to_h$se * 1.96)

# GAM bs = cr
tseed_m3c_l <- tseed_m3_l
tseed_m3c_a <- tseed_m3_a
tseed_m3c_h <- tseed_m3_h

lseed_m4c_l <- lseed_m4_l
lseed_m4c_a <- lseed_m4_a
lseed_m4c_h <- lseed_m4_h

ptseed_m3c_l <- predict(seedwt_g3.2t, newdata = tseed_m3c_l, se.fit = T, type = 'response')
tseed_m3c_l$fit <- ptseed_m3c_l$fit
tseed_m3c_l$se <- ptseed_m3c_l$se.fit
tseed_m3c_l$lci <- tseed_m3c_l$fit - (tseed_m3c_l$se* 1.96)
tseed_m3c_l$uci <- tseed_m3c_l$fit + (tseed_m3c_l$se * 1.96)

ptseed_m3c_a <- predict(seedwt_g3.2t, newdata = tseed_m3c_a, se.fit = T, type = 'response')
tseed_m3c_a$fit <- ptseed_m3c_a$fit
tseed_m3c_a$se <- ptseed_m3c_a$se.fit
tseed_m3c_a$lci <- tseed_m3c_a$fit - (tseed_m3c_a$se* 1.96)
tseed_m3c_a$uci <- tseed_m3c_a$fit + (tseed_m3c_a$se * 1.96)

ptseed_m3c_h <- predict(seedwt_g3.2t, newdata = tseed_m3c_h, se.fit = T, type = 'response')
tseed_m3c_h$fit <- ptseed_m3c_h$fit
tseed_m3c_h$se <- ptseed_m3c_h$se.fit
tseed_m3c_h$lci <- tseed_m3c_h$fit - (tseed_m3c_h$se* 1.96)
tseed_m3c_h$uci <- tseed_m3c_h$fit + (tseed_m3c_h$se * 1.96)




plseed_m4c_l <- predict(seedwt_g4.2l, newdata = lseed_m4c_l, se.fit = T, type = 'response')
lseed_m4c_l$fit <- plseed_m4c_l$fit
lseed_m4c_l$se <- plseed_m4c_l$se.fit
lseed_m4c_l$lci <- lseed_m4c_l$fit - (lseed_m4c_l$se * 1.96)
lseed_m4c_l$uci <- lseed_m4c_l$fit + (lseed_m4c_l$se * 1.96)

plseed_m4c_a <- predict(seedwt_g4.2l, newdata = lseed_m4c_a, se.fit = T, type = 'response')
lseed_m4c_a$fit <- plseed_m4c_a$fit
lseed_m4c_a$se <- plseed_m4c_a$se.fit
lseed_m4c_a$lci <- lseed_m4c_a$fit - (lseed_m4c_a$se * 1.96)
lseed_m4c_a$uci <- lseed_m4c_a$fit + (lseed_m4c_a$se * 1.96)

plseed_m4c_h <- predict(seedwt_g4.2l, newdata = lseed_m4c_h, se.fit = T, type = 'response')
lseed_m4c_h$fit <- plseed_m4c_h$fit
lseed_m4c_h$se <- plseed_m4c_h$se.fit
lseed_m4c_h$lci <- lseed_m4c_h$fit - (lseed_m4c_h$se * 1.96)
lseed_m4c_h$uci <- lseed_m4c_h$fit + (lseed_m4c_h$se * 1.96)


# GAM bs = cr optimised
tseed_m3co_l <- tseed_m3_l
tseed_m3co_a <- tseed_m3_a
tseed_m3co_h <- tseed_m3_h

lseed_m4co_l <- lseed_m4_l
lseed_m4co_a <- lseed_m4_a
lseed_m4co_h <- lseed_m4_h


ptseed_m3co_l <- predict(seedwt_g3.2ot, newdata = tseed_m3co_l, se.fit = T, type = 'response')
tseed_m3co_l$fit <- ptseed_m3co_l$fit
tseed_m3co_l$se <- ptseed_m3co_l$se.fit
tseed_m3co_l$lci <- tseed_m3co_l$fit - (tseed_m3co_l$se* 1.96)
tseed_m3co_l$uci <- tseed_m3co_l$fit + (tseed_m3co_l$se * 1.96)

ptseed_m3co_a <- predict(seedwt_g3.2ot, newdata = tseed_m3co_a, se.fit = T, type = 'response')
tseed_m3co_a$fit <- ptseed_m3co_a$fit
tseed_m3co_a$se <- ptseed_m3co_a$se.fit
tseed_m3co_a$lci <- tseed_m3co_a$fit - (tseed_m3co_a$se* 1.96)
tseed_m3co_a$uci <- tseed_m3co_a$fit + (tseed_m3co_a$se * 1.96)

ptseed_m3co_h <- predict(seedwt_g3.2ot, newdata = tseed_m3co_h, se.fit = T, type = 'response')
tseed_m3co_h$fit <- ptseed_m3co_h$fit
tseed_m3co_h$se <- ptseed_m3co_h$se.fit
tseed_m3co_h$lci <- tseed_m3co_h$fit - (tseed_m3co_h$se* 1.96)
tseed_m3co_h$uci <- tseed_m3co_h$fit + (tseed_m3co_h$se * 1.96)


plseed_m4co_l <- predict(seedwt_g4.2ol, newdata = lseed_m4co_l, se.fit = T, type = 'response')
lseed_m4co_l$fit <- plseed_m4co_l$fit
lseed_m4co_l$se <- plseed_m4co_l$se.fit
lseed_m4co_l$lci <- lseed_m4co_l$fit - (lseed_m4co_l$se * 1.96)
lseed_m4co_l$uci <- lseed_m4co_l$fit + (lseed_m4co_l$se * 1.96)

plseed_m4co_a <- predict(seedwt_g4.2ol, newdata = lseed_m4co_a, se.fit = T, type = 'response')
lseed_m4co_a$fit <- plseed_m4co_a$fit
lseed_m4co_a$se <- plseed_m4co_a$se.fit
lseed_m4co_a$lci <- lseed_m4co_a$fit - (lseed_m4co_a$se * 1.96)
lseed_m4co_a$uci <- lseed_m4co_a$fit + (lseed_m4co_a$se * 1.96)

plseed_m4co_h <- predict(seedwt_g4.2ol, newdata = lseed_m4co_h, se.fit = T, type = 'response')
lseed_m4co_h$fit <- plseed_m4co_h$fit
lseed_m4co_h$se <- plseed_m4co_h$se.fit
lseed_m4co_h$lci <- lseed_m4co_h$fit - (lseed_m4co_h$se * 1.96)
lseed_m4co_h$uci <- lseed_m4co_h$fit + (lseed_m4co_h$se * 1.96)


# GAM bs = 'cc'
lseed_m4cc_l <- lseed_m4_l
lseed_m4cc_a <- lseed_m4_a
lseed_m4cc_h <- lseed_m4_h


plseed_m4cc_l <- predict(seedwt_g4.3l, newdata = lseed_m4cc_l, se.fit = T, type = 'response')
lseed_m4cc_l$fit <- plseed_m4cc_l$fit
lseed_m4cc_l$se <- plseed_m4cc_l$se.fit
lseed_m4cc_l$lci <- lseed_m4cc_l$fit - (lseed_m4cc_l$se * 1.96)
lseed_m4cc_l$uci <- lseed_m4cc_l$fit + (lseed_m4cc_l$se * 1.96)

plseed_m4cc_a <- predict(seedwt_g4.3l, newdata = lseed_m4cc_a, se.fit = T, type = 'response')
lseed_m4cc_a$fit <- plseed_m4cc_a$fit
lseed_m4cc_a$se <- plseed_m4cc_a$se.fit
lseed_m4cc_a$lci <- lseed_m4cc_a$fit - (lseed_m4cc_a$se * 1.96)
lseed_m4cc_a$uci <- lseed_m4cc_a$fit + (lseed_m4cc_a$se * 1.96)

plseed_m4cc_h <- predict(seedwt_g4.3l, newdata = lseed_m4cc_h, se.fit = T, type = 'response')
lseed_m4cc_h$fit <- plseed_m4cc_h$fit
lseed_m4cc_h$se <- plseed_m4cc_h$se.fit
lseed_m4cc_h$lci <- lseed_m4cc_h$fit - (lseed_m4cc_h$se * 1.96)
lseed_m4cc_h$uci <- lseed_m4cc_h$fit + (lseed_m4cc_h$se * 1.96)


# Plot predictions for seed size ----
# Torulosa
dev.new(width = 20, height = 12, noRStudioGD = T, dpi = 300)
par(mfrow = c(2,3), mar = c(8,6,3,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))

plot(tseed_m3d_l$FPC, tseed_m3d_l$fit, type = 'l', xlim = c(35, 65), xlab = "", ylab = "", las = 1, ylim =c (3,8), col = 'blue', cex.axis = 1.4)
mtext(side = 1, expression(bold("Foliage Projective Cover (%)")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 2.5, cex = 1.5)
mtext(expression(bold("(a) GAM defaults")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g3t), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'FPC', data = "tseed_m3d_l", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m3d_a$FPC, tseed_m3d_a$fit, type = 'l', col = 'black')
pg.ci(x = 'FPC', data = "tseed_m3d_a", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m3d_h$FPC, tseed_m3d_h$fit, type = 'l', col = 'red')
pg.ci(x = 'FPC', data = "tseed_m3d_h", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(tseed_m3t_l$FPC, tseed_m3t_l$fit, type = 'l', xlim = c(35, 65), xlab = "", ylab = "", las = 1, ylim =c (3,8), col = 'blue', cex.axis = 1.4)
mtext(side = 1, expression(bold("Foliage Projective Cover (%)")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 2.5, cex = 1.5)
mtext(expression(bold("(b) GAM bs = 'tp'")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g3.1t), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'FPC', data = "tseed_m3t_l", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m3t_a$FPC, tseed_m3t_a$fit, type = 'l', col = 'black')
pg.ci(x = 'FPC', data = "tseed_m3t_a", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m3t_h$FPC, tseed_m3t_h$fit, type = 'l', col = 'red')
pg.ci(x = 'FPC', data = "tseed_m3t_h", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(tseed_m3to_l$FPC, tseed_m3to_l$fit, type = 'l', xlim = c(35, 65), xlab = "", ylab = "", las = 1, ylim =c (3,8), col = 'blue', cex.axis = 1.4)
mtext(side = 1, expression(bold("Foliage Projective Cover (%)")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 2.5, cex = 1.5)
mtext(expression(bold("(c) GAM bs = 'tp' optimised")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g3.1ot), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'FPC', data = "tseed_m3to_l", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m3to_a$FPC, tseed_m3to_a$fit, type = 'l', col = 'black')
pg.ci(x = 'FPC', data = "tseed_m3to_a", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m3to_h$FPC, tseed_m3to_h$fit, type = 'l', col = 'red')
pg.ci(x = 'FPC', data = "tseed_m3to_h", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


par(xpd = NA)
legend(x = 67, y = 8, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", "black", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)


plot(tseed_m3c_l$FPC, tseed_m3c_l$fit, type = 'l', xlim = c(35, 65), xlab = "", ylab = "", las = 1, ylim =c (3,8), col = 'blue', cex.axis = 1.4)
mtext(side = 1, expression(bold("Foliage Projective Cover (%)")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 2.5, cex = 1.5)
mtext(expression(bold("(d) GAM bs = 'cr'")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g3.2t), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'FPC', data = "tseed_m3c_l", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m3c_a$FPC, tseed_m3c_a$fit, type = 'l', col = 'black')
pg.ci(x = 'FPC', data = "tseed_m3c_a", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m3c_h$FPC, tseed_m3c_h$fit, type = 'l', col = 'red')
pg.ci(x = 'FPC', data = "tseed_m3c_h", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(tseed_m3co_l$FPC, tseed_m3co_l$fit, type = 'l', xlim = c(35, 65), xlab = "", ylab = "", las = 1, ylim =c (3,8), col = 'blue', cex.axis = 1.4)
mtext(side = 1, expression(bold("Foliage Projective Cover (%)")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 2.5, cex = 1.5)
mtext(expression(bold("(e) GAM bs = 'cr' optimised")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g3.2ot), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'FPC', data = "tseed_m3co_l", colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m3co_a$FPC, tseed_m3co_a$fit, type = 'l', col = 'black')
pg.ci(x = 'FPC', data = "tseed_m3co_a", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m3co_h$FPC, tseed_m3co_h$fit, type = 'l', col = 'red')
pg.ci(x = 'FPC', data = "tseed_m3co_h", colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')



plot(tseed_m2lm_l$r_Latitude, tseed_m2lm_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4, xaxt = "n", ylim =c (3,8),)
axis(side = 1, at = seq(-0.9, 1.3, 0.2), labels = seq(-28.28, -27.45, 0.07), cex.axis = 1.4)
mtext(side = 1, expression(bold("Latitude")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 2.5, cex = 1.5)
mtext(expression(bold("(f) GLMER")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_m2t), 3), sep = ""), cex =1.2, line = -1.5)
pg.ci(x = 'r_Latitude', data = 'tseed_m2lm_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m2lm_a$r_Latitude, tseed_m2lm_a$fit, type = 'l', col = 'black')
pg.ci(x = 'r_Latitude', data = "tseed_m2lm_a", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(tseed_m2lm_h$r_Latitude, tseed_m2lm_h$fit, type = 'l', col = 'red')
pg.ci(x = 'r_Latitude', data = 'tseed_m2lm_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')



# Littoralis
dev.new(width = 20, height = 12, noRStudioGD = T, dpi = 300)
par(mfrow = c(2,3), mar = c(8,6,3,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))

plot(lseed_m2d_l$Latitude, lseed_m2d_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4) 
mtext(side = 1, expression(bold("Latitude")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 2.5, cex = 1.5)
mtext(expression(bold("(a) GAM default")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g2l), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Latitude', data = 'lseed_m2d_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m2d_a$Latitude, lseed_m2d_a$fit, type = 'l', col = 'black')
pg.ci(x = 'Latitude', data = 'lseed_m2d_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m2d_h$Latitude, lseed_m2d_h$fit, type = 'l', col = 'red')
pg.ci(x = 'Latitude', data = 'lseed_m2d_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')

plot(lseed_m2t_l$Latitude, lseed_m2t_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4) 
mtext(side = 1, expression(bold("Latitude")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 2.5, cex = 1.5)
mtext(expression(bold("(b) GAM bs = 'tp'")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g2.1l), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Latitude', data = 'lseed_m2t_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m2t_a$Latitude, lseed_m2t_a$fit, type = 'l', col = 'black')
pg.ci(x = 'Latitude', data = 'lseed_m2t_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m2t_h$Latitude, lseed_m2t_h$fit, type = 'l', col = 'red')
pg.ci(x = 'Latitude', data = 'lseed_m2t_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(lseed_m4to_l$Precip, lseed_m4to_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4, ylim = c(1.5, 4)) 
mtext(side = 1, expression(bold("Precipitation seasonality")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 3, cex = 1.5)
mtext(expression(bold("(c) GAM bs = 'tp' optimised")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g2.1ol), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Precip', data = 'lseed_m4to_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4to_a$Precip, lseed_m4to_a$fit, type = 'l', col = 'black')
pg.ci(x = 'Precip', data = 'lseed_m4to_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4to_h$Precip, lseed_m4to_h$fit, type = 'l', col = 'red')
pg.ci(x = 'Precip', data = 'lseed_m4to_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


par(xpd = NA)
legend(x = 49.8, y = 4, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", "black", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)


plot(lseed_m4c_l$Precip, lseed_m4c_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4, ylim = c(1.5, 4)) 
mtext(side = 1, expression(bold("Precipitation seasonality")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 3, cex = 1.5)
mtext(expression(bold("(d) GAM bs = 'cr'")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g2.2l), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Precip', data = 'lseed_m4c_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4c_a$Precip, lseed_m4c_a$fit, type = 'l', col = 'black')
pg.ci(x = 'Precip', data = 'lseed_m4c_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4c_h$Precip, lseed_m4c_h$fit, type = 'l', col = 'red')
pg.ci(x = 'Precip', data = 'lseed_m4c_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(lseed_m4co_l$Precip, lseed_m4co_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4, ylim = c(1.5, 4)) 
mtext(side = 1, expression(bold("Precipitation seasonality")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 3, cex = 1.5)
mtext(expression(bold("(e) GAM bs = 'cr' optimised")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g2.2ol),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Precip', data = 'lseed_m4co_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4co_a$Precip, lseed_m4co_a$fit, type = 'l', col = 'black')
pg.ci(x = 'Precip', data = 'lseed_m4co_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4co_h$Precip, lseed_m4co_h$fit, type = 'l', col = 'red')
pg.ci(x = 'Precip', data = 'lseed_m4co_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(lseed_m5lm_l$r_Temp, lseed_m5lm_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4, xaxt = "n", xlim = c(-1.1, 1.9), ylim = c(4.4,5.2))
axis(side = 1, at = seq(-1.1, 1.9, 0.3), labels = seq(392, 432, 4), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 3, cex = 1.5)
mtext(expression(bold("(f) GLMER")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_m5l), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'r_Temp', data = 'lseed_m5lm_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m5lm_a$r_Temp, lseed_m5lm_a$fit, type = 'l', col = 'black')
pg.ci(x = 'r_Temp', data = 'lseed_m5lm_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m5lm_h$r_Temp, lseed_m5lm_h$fit, type = 'l', col = 'red')
pg.ci(x = 'r_Temp', data = 'lseed_m5lm_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
# CIs so small can't be seen?


# Littoralis compare cubic regression splines
dev.new(width = 20, height = 12, noRStudioGD = T, dpi = 300)
par(mfrow = c(2,3), mar = c(8,6,3,2), mgp = c(2.7,1,0), oma = c(0,0,0,10))


plot(lseed_m4to_l$Precip, lseed_m4to_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4, ylim = c(1.5, 4)) 
mtext(side = 1, expression(bold("Precipitation seasonality")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 3, cex = 1.5)
mtext(expression(bold("(a) GAM bs = 'tp' optimised")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g2.1ol), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Precip', data = 'lseed_m4to_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4to_a$Precip, lseed_m4to_a$fit, type = 'l', col = 'black')
pg.ci(x = 'Precip', data = 'lseed_m4to_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4to_h$Precip, lseed_m4to_h$fit, type = 'l', col = 'red')
pg.ci(x = 'Precip', data = 'lseed_m4to_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(lseed_m4c_l$Precip, lseed_m4c_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4, ylim = c(1.5, 4)) 
mtext(side = 1, expression(bold("Precipitation seasonality")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 3, cex = 1.5)
mtext(expression(bold("(b) GAM bs = 'cr'")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g2.2l), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Precip', data = 'lseed_m4c_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4c_a$Precip, lseed_m4c_a$fit, type = 'l', col = 'black')
pg.ci(x = 'Precip', data = 'lseed_m4c_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4c_h$Precip, lseed_m4c_h$fit, type = 'l', col = 'red')
pg.ci(x = 'Precip', data = 'lseed_m4c_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(lseed_m4co_l$Precip, lseed_m4co_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4, ylim = c(1.5, 4)) 
mtext(side = 1, expression(bold("Precipitation seasonality")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 3, cex = 1.5)
mtext(expression(bold("(c) GAM bs = 'cr' optimised")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g2.2ol),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Precip', data = 'lseed_m4co_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4co_a$Precip, lseed_m4co_a$fit, type = 'l', col = 'black')
pg.ci(x = 'Precip', data = 'lseed_m4co_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4co_h$Precip, lseed_m4co_h$fit, type = 'l', col = 'red')
pg.ci(x = 'Precip', data = 'lseed_m4co_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')



par(xpd = NA)
legend(x = 49.8, y = 4, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", "black", 'red'), title = expression(bold("Fire frequency")), lty = 1, lwd = 2, cex = 1.8, bty = "n")
par(xpd = F)



plot(lseed_m4cc_l$Precip, lseed_m4cc_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4, ylim = c(1.5, 4)) 
mtext(side = 1, expression(bold("Precipitation seasonality")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 3, cex = 1.5)
mtext(expression(bold("(d) GAM bs = 'cc' optimised")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_g2.2ol),3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'Precip', data = 'lseed_m4cc_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4cc_a$Precip, lseed_m4cc_a$fit, type = 'l', col = 'black')
pg.ci(x = 'Precip', data = 'lseed_m4cc_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m4cc_h$Precip, lseed_m4cc_h$fit, type = 'l', col = 'red')
pg.ci(x = 'Precip', data = 'lseed_m4cc_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')


plot(lseed_m5lm_l$r_Temp, lseed_m5lm_l$fit, type = 'l', xlab = "", ylab = "", las = 1, col = 'blue', cex.axis = 1.4, xaxt = "n", xlim = c(-1.1, 1.9), ylim = c(4.4,5.2))
axis(side = 1, at = seq(-1.1, 1.9, 0.3), labels = seq(392, 432, 4), cex.axis = 1.4)
mtext(side = 1, expression(bold("Temperature seasonality")), line = 3, cex = 1.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), line = 3, cex = 1.5)
mtext(expression(bold("(e) GLMER")), cex = 2)
mtext(paste("AICc = ", round(AICc(seedwt_m5l), 3), sep = ""), line = -1.5, cex = 1.2)
pg.ci(x = 'r_Temp', data = 'lseed_m5lm_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m5lm_a$r_Temp, lseed_m5lm_a$fit, type = 'l', col = 'black')
pg.ci(x = 'r_Temp', data = 'lseed_m5lm_a', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m5lm_h$r_Temp, lseed_m5lm_h$fit, type = 'l', col = 'red')
pg.ci(x = 'r_Temp', data = 'lseed_m5lm_h', colour = rgb(1,0,0,0.1), lower = 'lci', upper = 'uci')
# CIs so small can't be seen?


# Notes on analyses: ----
# 1. We use a binomial regression family as proportion germination is bounded by 0 and 1 and beta regression was not suitable as we have one-inflated data. 


# 2. We use Akaike's information criterion corrected for small sample sizes to rank models with the best model chosen as the model with the lowest AICc, which improves model fit over the null model by a change in AICc > 2. 


# 3. Smooth terms in gam are can be many things, some of which are 's' or 'te/t2' or ' ti'. When fitting interactions we want to use tensor product smooths (t2/te/ti) as this accounts for variables being on different scales (i.e., seed weight = mg, germination rates = % or days), thus, not requiring any variable rescaling. However, for fitting models with the main effect and interaction terms, it is most suitable to use 'ti' for interaction terms and 's' for the main effects. However, if the main effect is a factor variable, we do not fit this with a smoother (see #7 for exceptions).
# https://stats.stackexchange.com/questions/619212/in-r-mgcv-whats-the-difference-between-sx-sy-tix-y-and-tix-ti
# https://stats.stackexchange.com/questions/622071/mgcv-use-of-s-or-te-with-interactions-in-gams
# https://stats.stackexchange.com/questions/395557/difference-between-s-and-ti-terms-in-mgcv-package-when-applied-to-one-variab/409193#409193


# 4. Models with only a main effect of factor variables are basically a linear model  but if we were to fit this as a linear model and compare to the rest of the models fit as a gam() this changes assumptions for best model selection.
# https://stackoverflow.com/questions/63726043/gam-with-only-categorical-logical


# 5. We use mgcv::gam() as ti(), most suitable for modelling interactions, does not work with gamm4() and it is unclear whether t2 functions in the same manner. However, AIC does not work well with mgcv, so we manually calculate these values for model comparisons.


# 6. The marginal basis function (the basis for the smooth term) can equal a few things; e.g. fs (random factor smooth interactions), re (random effects)  tp (default numeric), cr (cubic regression spline) or cc (cyclic regression spline). When we fit interactions with smooth terms, rather than using fs, we specify that the interaction should be by the factor variable so that separate smooths are fit for interactions with each level of the factor variable.
# https://stats.stackexchange.com/questions/608112/predicting-with-gam-mgcv-and-categorical-factor-covariate-in-r
# https://stackoverflow.com/questions/63726043/gam-with-only-categorical-logical


# 7. Where the main effect of a factor variable is included in models with interactive terms we fit the model with s(bs = 're'), this ensures intercepts are calculated for the factor main effect.
# Page 17 https://peerj.com/articles/6876/


# 8. Generally for ecological datasets it is recommended to specify method = 'REML/ML' to account for small sample sizes. This will reduce the likelihood of over fitting resulting in splines that are too wiggly by the default GCV method. However, method = 'REML' is not appropriate when comparing models with different fixed effects or non-fully penalized smooths (i.e., most basis functions except 'fs' or 're') so we instead use ML.
# https://stats.stackexchange.com/questions/301364/gam-optimization-methods-in-mgcv-r-package-which-to-choose




# For Proportion_seedlings ~ TSF choose modelling method being used as the null model is best no matter if we consider GLMER or GAM. 


save.image('./02_Workspaces/002_full_analysis_exploration.RData')



