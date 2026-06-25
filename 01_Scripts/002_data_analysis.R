# Written by Felicity Charles
# Caveat emptor
# Date: 12th March 2025

# Analysis of full factorial germination experiment and population age structure surveys. 


#1. Load packages ----
library(AICcmodavg)
library(mgcv)
library(arm)
library(lme4)
library(ggplot2)
library(cowplot)
library(emmeans)
library(multcomp)
library(multcompView)

# 1.1 Load custom functions ----
invisible(lapply(paste("./04_Functions/", dir ("04_Functions"), sep = ""), function(x) source (x)))


# 2. Read in data ----
# For analyses, we want to analyse the species separately as we expect them to have different responses. 
# Germination data with seed characteristics and environmental data
dat_cum.prop <- read.csv('./00_Data/Seeds_data/Full_experiment/Full_experiment_cumulative_germ.csv', header = T, stringsAsFactors = T)
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
#tor_tree[90, 10] <- "sapling" # fix incorrect age category label
tor_tree$Age <- factor(tor_tree$Age, levels = c('seedling', 'sapling', 'mature'))
lit_tree <- dat_tree[dat_tree$Species == "littoralis", ]
# For GB1, this transect was not sampled as comprehensively as other transects as this was the first transect surveyed so much tree level data is missing. It has been limited to 50m, where it was originally measured out to 100m.
lit_tree$Age <- factor(lit_tree$Age, levels = c('seedling', 'sapling', 'mature'))

# Extract data for investigating cone number
tor_cones <- tor_tree[which(tor_tree$Cone_presence == "Y"),]
lit_cones <- lit_tree[which(lit_tree$Cone_presence == "Y"),] 
unique(tor_cones$Cone_number)
unique(lit_cones$Cone_number) # Note there are some NAs present as the Gatton transects did not record the number of cones if they were present. We will remove these rows.
lit_cones <- lit_cones[!is.na(lit_cones$Cone_number),]

# Calculate height in cm for cones data
tor_cones$Height_cm <- tor_cones$Height_.m.*100
lit_cones$Height_cm <- lit_cones$Height_.m.*100




# Look at what size do these species have cones, plot DBH vs height histograms

plot(tor_cones$Height_.m., tor_cones$Cone_number, pch = 20, xlim = c(2, 17), xaxt = 'n', xlab = "Height (m)", ylab = "Cone number")
axis(side = 1, at = seq(2,17, 1))
# All torulosa individuals bearing cones are over 1m tall, so reproductive maturity is reached when they are greater than 1m. Need to look at a dataset that includes all individuals

# If cones are present what is the distribution

dev.new(width = 18, height = 8, dpi = 80, pointsize = 16, noRStudioGD = T)
par(mfrow = c(1,2), mar = c(4,4,4,2))

plot(lit_cones$Height_.m., lit_cones$Cone_number, pch = 20, col = lit_cones$Age, xlim = c(0,10), las = 1, cex.axis = 1.3, xaxt = 'n', xlab = '', ylab = '', cex = 1.8)
axis(side = 1, at = seq(0,10,1), cex.axis = 1.3) # One individual which which was tiny with many cones. Only outlier, all others are reproductively mature.
axis(side = 2, at = seq(0,200,10), labels = F, lwd = 0.3)
mtext(side = 1, expression(bold("Height (m)")), line = 2.7, cex = 1.75)
mtext(side = 2, expression(bold("Number of cones")), line = 2.7, cex = 1.75)
mtext("(a)", cex = 1.75, at = 0, line = 0.2)
mtext(~italic(Allocasuarina~littoralis), cex = 1.75, at = 10, line = 2)

plot(lit_cones$DBH, lit_cones$Cone_number, pch = 20, col = lit_cones$Age, las = 1, cex.axis = 1.3, xaxt = 'n', xlab = '', ylab = '', cex = 1.8, xlim = c(0,12))
axis(side = 1, at = seq(0,12,1), cex.axis = 1.3)
axis(side = 2, at = seq(0, 500, 10), labels = F, lwd = 0.3)
mtext(side = 1, expression(bold("Diameter at breast height (cm)")), line = 2.7, cex = 1.75)
mtext(side = 2, expression(bold("Number of cones")), line = 2.7, cex = 1.75)
mtext("(b)", cex = 1.75, at = 0, line = 0.2)


dev.new(width = 18, height = 8, dpi = 80, pointsize = 16, noRStudioGD = T)
par(mfrow = c(1,2), mar = c(4,4,4,2))

plot(tor_cones$Height_.m., tor_cones$Cone_number, pch = 20, col = tor_cones$Age, cex.axis = 1.3, xaxt = 'n', xlab = '', ylab = '', cex = 1.8, xlim = c(0,18), las = 1)
axis(side = 1, at = seq(0,18,1), cex.axis = 1.3)
axis(side = 2, at = seq(0, 500, 10), labels = F, lwd = 0.3)
mtext(side = 1, expression(bold("Height (m)")), line = 2.7, cex = 1.75)
mtext(side = 2, expression(bold("Number of cones")), line = 2.7, cex = 1.75)
mtext("(a)", cex = 1.75, at = 0, line = 0.2)
mtext(~italic(Allocasuarina~torulosa), cex = 1.75, at = 18, line = 2)

plot(tor_cones$DBH, tor_cones$Cone_number, pch = 20, col = tor_cones$Age, cex.axis = 1.3, xaxt = 'n', xlab = '', ylab = '', cex = 1.8, xlim = c(0, 35), las = 1)
axis(side = 1, at = seq(0,36,2), cex.axis = 1.3)
axis(side = 2, at = seq(0, 500, 10), labels = F, lwd = 0.3)
mtext(side = 1, expression(bold("Diameter at breast height (cm)")), line = 2.7, cex = 1.75)
mtext(side = 2, expression(bold("Number of cones")), line = 2.7, cex = 1.75)
mtext("(b)", cex = 1.75, at = 0, line = 0.2)






dev.new(width = 20, height = 8, dpi = 80, pointsize = 12, noRStudioGD = T)
par(mfrow = c(1,2), mar = c(4,6,4,2))


plot(lit_tree$DBH, lit_tree$Height_.m., pch = 20, xlab = "", ylab = '', xaxt = 'n', yaxt = 'n', xlim = c(0,10), ylim = c(0,16), cex = 1.8)
# Generally proportional increase in DBH with height.
axis(side = 1, at = seq(0,10, 2), cex.axis = 1.3)
axis(side = 2, at = seq(0,16, 2), las = 1, cex.axis = 1.3)
arrows(x0 = -2, y0 = 1, x1 = 55, y1 = 1, code = 0) # Seedlings
text(expression(bold("Seedlings")), x = 8, y = 0.5, cex = 1.3)
arrows(x0 = 3, y0 = 1, x1 = 3, y1 = 17, code = 0)
text(expression(bold("Saplings")), x = 1, y = 16, cex = 1.3)
text(expression(bold("Adults")), x = 8, y = 16, cex = 1.3)
mtext(side = 1, expression(bold("Diameter at breast height (cm)")), line = 2.7, cex = 1.75)
mtext(side = 2, expression(bold("Height (m)")), line = 2.5, cex = 1.75)
mtext("(a) "~italic(Allocasuarina~littoralis), cex = 1.75, adj = 0.001)


plot(tor_tree$DBH, tor_tree$Height_.m., pch = 20, xlab = '', ylab = ' ', xaxt = 'n', yaxt = 'n', xlim = c(0,49), cex = 1.8)
# Generally proportional increase in DBH with height.
axis(side = 1, at = seq(0,48, 2), cex.axis = 1.3)
axis(side = 2, at = seq(0,16, 2), las = 1, cex.axis = 1.3)
arrows(x0 = -2, y0 = 1, x1 = 55, y1 = 1, code = 0) # Seedlings
text(expression(bold("Seedlings")), x = 26, y = 0.5, cex = 1.3)
arrows(x0 = 3, y0 = 1, x1 = 3, y1 = 17, code = 0)
text(expression(bold("Saplings")), x = 1, y = 16, cex = 1.3)
text(expression(bold("Adults")), x = 42, y = 16, cex = 1.3)
mtext(side = 1, expression(bold("Diameter at breast height (cm)")), line = 2.7, cex = 1.75)
mtext(side = 2, expression(bold("Height (m)")), line = 2.5, cex = 1.75)
mtext("(b) "~italic(Allocasuarina~torulosa), cex = 1.75, adj = 0.001)


# <1m = seedling, >1m but DBH <3cm = sapling, adult = >1m and DBH >3cm






# 4. Check for correlations and rescale variables ----
# We consider Spearmans correlation to be highly correlated if >0.7.
lit_germc <- ggstatsplot::ggcorrmat(lit_cum.prop,
                                    type = "non-parametric",
                                    label = T,
                                    cor.vars = c("Latitude", "seed_wt_mg", "Fire_freq", "TSF", "FPC", "Precip", "Temp", "TWI"),
                                    size = 2)
lit_germc # Temperature and latitude are highly correlated, this is to be expected as temperature increases with decreasing latitude. These variables will not be included in the same analyses so we do not need to exclude either. TSF and Fire frequency are highly correlated but these variables are not included in the same analyses.  TSF and temperature are also highly correlated so some caution needs to be taken if using a model including these terms for Allocasuarina littoralis germination. 

tor_germc <- ggstatsplot::ggcorrmat(tor_cum.prop,
                                    type = "non-parametric",
                                    label = T,
                                    cor.vars = c("Latitude", "seed_wt_mg", "Fire_freq", "TSF", "FPC", "Precip", "Temp", "TWI"),
                                    size = 2)
tor_germc # No variables were highly correlated.


lit_transect_cor <- ggstatsplot::ggcorrmat(lit_transects,
                                           type = "non-parametric",
                                           label = T,
                                           cor.vars = c("Location", "Latitude", "Fire_freq", "TSF", "FPC", "Precip", "Temp", "TWI"),
                                           size = 2)
lit_transect_cor # As with the germination data temperature and latitude and time since fire and fire frequency are highly correlated, but these variables will not be included in the same analyses. We also expect that these variables will be highly correlated

tor_transect_cor <- ggstatsplot::ggcorrmat(tor_transects,
                                           type = "non-parametric",
                                           label = T,
                                           cor.vars = c("Location", "Latitude", "Fire_freq", "TSF", "FPC", "Precip", "Temp", "TWI"),
                                           size = 2)
tor_transect_cor # No variables were highly correlated


lit_tree_cor <- ggstatsplot::ggcorrmat(lit_tree,
                                           type = "non-parametric",
                                           label = T,
                                           cor.vars = c("Location", "Latitude", "Fire_freq", "TSF", "FPC", "Precip", "Temp", "Height_.m.", "TWI"),
                                           size = 2)
lit_tree_cor # Temperature and latitude are highly correlated, will not be included in the same analyses. 

tor_tree_cor <- ggstatsplot::ggcorrmat(tor_tree,
                                       type = "non-parametric",
                                       label = T,
                                       cor.vars = c("Location", "Latitude", "Fire_freq", "TSF", "FPC", "Precip", "Temp", "Height_.m.", "TWI"),
                                       size = 2)
tor_tree_cor # TSF and fire frequency are highly correlated but will not be included in the same analyses.


# We need to rescale predictor variables as our response variable is proportional data. We need to wrap scale in c() for its dimension-stripping properties to ensure rescaling works correctly. Before rescaling latitude, we need to remove the '-', we can add this back on to axis labels later if we need to but if we keep this, we cannot produce plots correctly (i.e., more southern latitudes will be at the start of the plot.)
# 4.1 Rescale germination experiment data
tor_cum.prop$r_seed_wt <- c(scale(tor_cum.prop$seed_wt_mg))
tor_cum.prop$r_fire_freq <- c(scale(tor_cum.prop$Fire_freq))
tor_cum.prop$Latitude <- as.numeric(gsub("-", "", tor_cum.prop$Latitude))
tor_cum.prop$r_Latitude <- c(scale(tor_cum.prop$Latitude))
tor_cum.prop$r_FPC <- c(scale(tor_cum.prop$FPC))
tor_cum.prop$r_Temp <- c(scale(tor_cum.prop$Temp))
tor_cum.prop$r_Precip <- c(scale(tor_cum.prop$Temp))
tor_cum.prop$r_TWI <- c(scale(tor_cum.prop$TWI))

lit_cum.prop$r_seed_wt <- c(scale(lit_cum.prop$seed_wt_mg))
lit_cum.prop$r_fire_freq <- c(scale(lit_cum.prop$Fire_freq))
lit_cum.prop$Latitude <- as.numeric(gsub("-", "", lit_cum.prop$Latitude))
lit_cum.prop$r_Latitude <- c(scale(lit_cum.prop$Latitude))
lit_cum.prop$r_FPC <- c(scale(lit_cum.prop$FPC))
lit_cum.prop$r_Temp <- c(scale(lit_cum.prop$Temp))
lit_cum.prop$r_Precip <- c(scale(lit_cum.prop$Precip))
lit_cum.prop$r_TWI <- c(scale(lit_cum.prop$TWI))

# 4.2 Rescale transect level data
tor_transects$r_TSF <- c(scale(tor_transects$TSF))
tor_transects$Latitude <- as.numeric(gsub("-", "", tor_transects$Latitude))
tor_transects$r_Latitude <- c(scale(tor_transects$Latitude))
tor_transects$r_FPC <- c(scale(tor_transects$FPC))
tor_transects$r_Precip <- c(scale(tor_transects$Precip))
tor_transects$r_Temp <- c(scale(tor_transects$Temp))
tor_transects$r_Fire_freq <- c(scale(tor_transects$Fire_freq))
tor_transects$r_TWI <- c(scale(tor_transects$TWI))

lit_transects$r_TSF <- c(scale(lit_transects$TSF))
lit_transects$Latitude <- as.numeric(gsub("-", "", lit_transects$Latitude))
lit_transects$r_Latitude <- c(scale(lit_transects$Latitude))
lit_transects$r_FPC <- c(scale(lit_transects$FPC))
lit_transects$r_Precip <- c(scale(lit_transects$Precip))
lit_transects$r_Temp <- c(scale(lit_transects$Temp))
lit_transects$r_Fire_freq <- c(scale(lit_transects$Fire_freq))
lit_transects$r_TWI <- c(scale(lit_transects$TWI))

# 4.3 Rescale cone number data
tor_cones$r_TSF <- c(scale(tor_cones$TSF))
tor_cones$r_Fire_freq <- c(scale(tor_cones$Fire_freq))
tor_cones$r_Height <- c(scale(tor_cones$Height_cm))
tor_cones$Latitude <- as.numeric(gsub("-", "", tor_cones$Latitude))
tor_cones$r_Latitude <- c(scale(tor_cones$Latitude))
tor_cones$r_FPC <- c(scale(tor_cones$FPC))
tor_cones$r_Precip <- c(scale(tor_cones$Precip))
tor_cones$r_Temp <- c(scale(tor_cones$Temp))
tor_cones$r_TWI <- c(scale(tor_cones$TWI))

lit_cones$r_TSF <- c(scale(lit_cones$TSF))
lit_cones$r_Fire_freq <- c(scale(lit_cones$Fire_freq))
lit_cones$r_Height <- c(scale(lit_cones$Height_cm))
lit_cones$Latitude <- as.numeric(gsub("-", "", lit_cones$Latitude))
lit_cones$r_Latitude <- c(scale(lit_cones$Latitude))
lit_cones$r_FPC <- c(scale(lit_cones$FPC))
lit_cones$r_Precip <- c(scale(lit_cones$Precip))
lit_cones$r_Temp <- c(scale(lit_cones$Temp))
lit_cones$r_TWI <- c(scale(lit_cones$TWI))

# Notes on analyses: ----
# 1. We use a binomial regression family as proportion germination is bounded by 0 and 1 and beta regression was not suitable as we have one-inflated data. 

# 2. We use Akaike's information criterion corrected for small sample sizes to rank models with the best model chosen as the model with the lowest AICc, which improves model fit over the null model by a change in AICc > 2. 




# QUESTION 1: How are germination rates influenced by seed treatment, seed attributes and/or fire frequency? ----
# 5.1 Proportion germination ----
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
aictab(tp_aic) # Model 6 is the best model, no model is ranked within delta AICc <2.

lp_aic <- list(mnull_l, m1_l, m2_l, m3_l, m4_l, m5_l, m6_l)
aictab(lp_aic) # Model 2 is the best model, no model is ranked within delta AICc <2. 



# 5.1.1 Predict to new data for proportion germination ----
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


# 5.1.2 Plot predictions for proportion germination ----
dev.new(width = 20, height = 8, dpi = 300, noRStudioGD = T)
par(mfrow = c(1,2), mar = c(4,6,2.5,2))

plot(lit_cum.prop$r_seed_wt, lit_cum.prop$Proportion_germ, pch = 20, cex = 1, ylim = c(0,1), ylab = "", las = 1, xlab = "", xlim = c(-3, 2.6), xaxt = "n", cex.axis = 2)
axis(side = 1, at = c(-3, 0, 2.6), labels = seq(1, 3, 1), cex.axis = 2)
lines(new_lp$r_seed_wt, new_lp$fit)
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3, cex = 2.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.7, cex = 2.5)
pg.ci(x = "r_seed_wt", data = "new_lp", colour = rgb(0/255, 0/255, 0/255, 0.1), lower = "lci", upper = "uci")
mtext("(a) "~italic(Allocasuarina~littoralis), cex = 2.5, adj = 0.001)


plot(tor_cum.prop$r_seed_wt[tor_cum.prop$Fire_freq == 0], tor_cum.prop$Proportion_germ[tor_cum.prop$Fire_freq == 0], pch = 20, cex = 1, ylim = c(0,1), ylab = "", las = 1, xlab = "", xlim = c(-2.4, 2.6), xaxt = "n", cex.axis = 2, col = "blue")
points(tor_cum.prop$r_seed_wt[tor_cum.prop$Fire_freq == 3], tor_cum.prop$Proportion_germ[tor_cum.prop$Fire_freq == 3], col = 'gray36', pch = 20, cex = 1)
points(tor_cum.prop$r_seed_wt[tor_cum.prop$Fire_freq == 6], tor_cum.prop$Proportion_germ[tor_cum.prop$Fire_freq == 6], col = 'red', pch = 20, cex = 1)
axis(side = 1, at = seq(-2.4, 2.6, 1), labels = round(seq(min(tor_cum.prop$seed_wt_mg), max(tor_cum.prop$seed_wt_mg), length = 6)), cex.axis = 2)
lines(new_tpl$r_seed_wt, new_tpl$fit, col = "blue")
mtext(side = 1, expression(bold("Seed weight (mg)")), line = 3, cex = 2.5)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.7, cex = 2.5)
pg.ci(x = "r_seed_wt", data = "new_tpl", colour = rgb(0, 0, 1, 0.1), lower = "lci", upper = "uci")
lines(new_tpa$r_seed_wt, new_tpa$fit)
pg.ci(x = "r_seed_wt", data = "new_tpa", colour = rgb(92/255, 92/255, 92/255, 0.1), lower = 'lci', upper = 'uci')
lines(new_tph$r_seed_wt, new_tph$fit, col = 'red')
pg.ci(x = 'r_seed_wt', data = 'new_tph', colour = rgb(255/255, 0/255, 0/255, 0.1), lower = 'lci', upper = 'uci')
mtext("(b) "~italic(Allocasuarina~torulosa), cex = 2.5, adj = 0.001)
legend(x = 0.9, y = 0.33, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", 'gray36', 'red'), title = expression(bold("Fire frequency")), lty = 1, cex = 2, bty = "n")



# 5.1.3 Create a plot to compare proportion germination between treatments for SI ----
lit_prop <- data.frame(Treatment = as.factor(c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke")))
litprop_p <- predictSE(m1_l, newdata = lit_prop, se.fit = T, type = 'response')
lit_prop$fit <- litprop_p$fit
lit_prop$se <- litprop_p$se.fit
lit_prop$lci <- lit_prop$fit - (lit_prop$se * 1.96)
lit_prop$uci <- lit_prop$fit + (lit_prop$se * 1.96)

# Torulosa
tor_prop <- data.frame(Treatment = as.factor(c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke")))
torprop_p <- predictSE(m1_t, newdata = tor_prop, se.fit = T, type = 'response')
tor_prop$fit <- torprop_p$fit
tor_prop$se <- torprop_p$se.fit
tor_prop$lci <- tor_prop$fit - (tor_prop$se * 1.96)
tor_prop$uci <- tor_prop$fit + (tor_prop$se * 1.96)

dev.new(width = 17, height = 8, dpi = 300, noRStudioGD = T)
par(mfrow = c(1,2), mar = c(8,5.5,3.5,2))

plot.default(c(1:6), lit_prop$fit, type = 'p', pch = 20, xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 2, ylim = c(0,1.05), xlim = c(0.5, 6.5), cex = 2.5)
axis(side = 1, at = c(1:6), labels = c("Control", "80°C", "95°C", "Smoke", "", ""), cex.axis = 2)
axis(side = 1, at = c(5), labels = "80°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
axis(side = 1, at = c(6), labels = "95°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
axis(side = 1, at = 2, labels = "80°C", cex.axis = 2, tick = F)
axis(side = 1, at = 4, labels = "Smoke", cex.axis = 2, tick = F)
mtext(side = 1, expression(bold("Treatment")), line = 4.5, cex = 2.2)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 2.2)
arrows(c(1:6), lit_prop$lci, c(1:6), lit_prop$uci, length = 0.05, code = 3, angle = 90)
mtext("(a) "~italic(Allocasuarina~littoralis), cex = 2.5, adj = 0.001)


# Torulosa
plot.default(c(1:6), tor_prop$fit, type = 'p', pch = 20, xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 2, ylim = c(0,1.05), xlim = c(0.5, 6.5), cex = 2.5)
axis(side = 1, at = c(1:6), labels = c("Control", "80°C", "95°C", "Smoke", "", ""), cex.axis = 2)
axis(side = 1, at = c(5), labels = "80°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
axis(side = 1, at = c(6), labels = "95°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
axis(side = 1, at = 2, labels = "80°C", cex.axis = 2, tick = F)
axis(side = 1, at = 4, labels = "Smoke", cex.axis = 2, tick = F)
mtext(side = 1, expression(bold("Treatment")), line = 4.5, cex = 2.2)
mtext(side = 2, expression(bold("Proportion germination")), line = 3.5, cex = 2.2)
arrows(c(1:6), tor_prop$lci, c(1:6), tor_prop$uci, length = 0.05, code = 3, angle = 90)
mtext("(b) "~italic(Allocasuarina~torulosa), cex = 2.5, adj = 0.001)




# 5.2 Time to 50% Germination ----
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
aictab(tt_aic, modnames = c("Null", "Treatment", "Seed weight", "FF", "Treatment*Seed weight", "Treatment*FF", "Seed weight*FF")) # Models 1, 4, and 5 were ranked within delta AICc <2 and are better than the null

lt_aic <- list(mn_lt50, m1_lt50, m2_lt50, m3_lt50, m4_lt50, m5_lt50, m6_lt50)
aictab(lt_aic, modnames = c("Null", "Treatment", "Seed weight", "FF", "Treatment*Seed weight", "Treatment*FF", "Seed weight*FF")) # Model 1 is better than the null butthree the null is within delta AICc <2. 


# 5.2.1 Predict to new data for time to 50% germination ---- 
# Littoralis
nt50_m1_l <- data.frame(Treatment = as.factor(c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke")))
pt50_m1_l <- predictSE(m1_lt50, newdata = nt50_m1_l, se.fit = T, type = 'response')
nt50_m1_l$fit <- pt50_m1_l$fit
nt50_m1_l$se <- pt50_m1_l$se.fit
nt50_m1_l$lci <- nt50_m1_l$fit - (nt50_m1_l$se * 1.96)
nt50_m1_l$uci <- nt50_m1_l$fit + (nt50_m1_l$se * 1.96)

# Torulosa
# Model 1 
nt50_m1_t <- data.frame(Treatment = as.factor(c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke")))
pt50_m1_t <- predictSE(m1_tt50, newdata = nt50_m1_t, se.fit = T, type = 'response')
nt50_m1_t$fit <- pt50_m1_t$fit
nt50_m1_t$se <- pt50_m1_t$se.fit
nt50_m1_t$lci <- nt50_m1_t$fit - (nt50_m1_t$se * 1.96)
nt50_m1_t$uci <- nt50_m1_t$fit + (nt50_m1_t$se * 1.96)



# Model 4
nt50_m4_tl <- expand.grid(r_seed_wt = min(tor_cum.prop$r_seed_wt),
                          Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
pt50_m4_tl <- predictSE(m4_tt50, newdata = nt50_m4_tl, se.fit = T, type = 'response')
nt50_m4_tl$fit <- pt50_m4_tl$fit
nt50_m4_tl$se <- pt50_m4_tl$se.fit
nt50_m4_tl$lci <- nt50_m4_tl$fit - (nt50_m4_tl$se * 1.96)
nt50_m4_tl$uci <- nt50_m4_tl$fit + (nt50_m4_tl$se * 1.96)


nt50_m4_ta <- expand.grid(r_seed_wt = mean(tor_cum.prop$r_seed_wt),
                          Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
pt50_m4_ta <- predictSE(m4_tt50, newdata = nt50_m4_ta, se.fit = T, type = 'response')
nt50_m4_ta$fit <- pt50_m4_ta$fit
nt50_m4_ta$se <- pt50_m4_ta$se.fit
nt50_m4_ta$lci <- nt50_m4_ta$fit - (nt50_m4_ta$se * 1.96)
nt50_m4_ta$uci <- nt50_m4_ta$fit + (nt50_m4_ta$se * 1.96)



nt50_m4_th <- expand.grid(r_seed_wt = max(tor_cum.prop$r_seed_wt),
                          Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
pt50_m4_th <- predictSE(m4_tt50, newdata = nt50_m4_th, se.fit = T, type = 'response')
nt50_m4_th$fit <- pt50_m4_th$fit
nt50_m4_th$se <- pt50_m4_th$se.fit
nt50_m4_th$lci <- nt50_m4_th$fit - (nt50_m4_th$se * 1.96)
nt50_m4_th$uci <- nt50_m4_th$fit + (nt50_m4_th$se * 1.96)



# Model 5
nt50_m5_tl <- expand.grid(r_fire_freq = min(tor_cum.prop$r_fire_freq),
                          Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
pt50_m5_tl <- predictSE(m5_tt50, newdata = nt50_m5_tl, se.fit = T, type = 'response')
nt50_m5_tl$fit <- pt50_m5_tl$fit
nt50_m5_tl$se <- pt50_m5_tl$se.fit
nt50_m5_tl$lci <- nt50_m5_tl$fit - (nt50_m5_tl$se * 1.96)
nt50_m5_tl$uci <- nt50_m5_tl$fit + (nt50_m5_tl$se * 1.96)



nt50_m5_ta <- expand.grid(r_fire_freq = mean(tor_cum.prop$r_fire_freq),
                          Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
pt50_m5_ta <- predictSE(m5_tt50, newdata = nt50_m5_ta, se.fit = T, type = 'response')
nt50_m5_ta$fit <- pt50_m5_ta$fit
nt50_m5_ta$se <- pt50_m5_ta$se.fit
nt50_m5_ta$lci <- nt50_m5_ta$fit - (nt50_m5_ta$se * 1.96)
nt50_m5_ta$uci <- nt50_m5_ta$fit + (nt50_m5_ta$se * 1.96)



nt50_m5_th <- expand.grid(r_fire_freq = max(tor_cum.prop$r_fire_freq),
                          Treatment = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
pt50_m5_th <- predictSE(m5_tt50, newdata = nt50_m5_th, se.fit = T, type = 'response')
nt50_m5_th$fit <- pt50_m5_th$fit
nt50_m5_th$se <- pt50_m5_th$se.fit
nt50_m5_th$lci <- nt50_m5_th$fit - (nt50_m5_th$se * 1.96)
nt50_m5_th$uci <- nt50_m5_th$fit + (nt50_m5_th$se * 1.96)



# Plot predictions for time to 50% germination ----
dev.new(width = 24, height = 24, dpi = 300, noRStudioGD = T)
par(mfrow = c(2,2), mar = c(8,5.5,3.5,2))



# Torulosa
plot.default(c(1:6), nt50_m1_t$fit, type = 'p', pch = 20, xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 2, ylim = c(4,18), cex = 2)
axis(side = 1, at = c(1:6), labels = c("Control", "80°C", "95°C", "Smoke", "", ""), cex.axis = 2)
axis(side = 1, at = c(5), labels = "80°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
axis(side = 1, at = c(6), labels = "95°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
axis(side = 1, at = c(2), labels = "80°C", cex.axis = 2, tick = F)
axis(side = 1, at = c(4), labels = "Smoke", cex.axis = 2, tick = F)
mtext(side = 1, expression(bold("Treatment")), line = 4.5, cex = 1.8)
mtext(side = 2, expression(bold("Days to 50% germination")), line = 3, cex = 1.8)
arrows(c(1:6), nt50_m1_t$lci, c(1:6), nt50_m1_t$uci, length = 0.05, code = 3, angle = 90)
mtext("(a) ", cex = 2, adj = 0.001, line = 0.3)


plot.default(nt50_m4_ta$Treatment, nt50_m4_ta$fit, pch = 20, type = 'p', xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 2, ylim = c(4, 18), col = "#D94801", cex = 2)
axis(side = 1, at = c(1:6), labels = c("Control", "80°C", "95°C", "Smoke", "", ""), cex.axis = 2)
axis(side = 1, at = c(5), labels = "80°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
axis(side = 1, at = c(2), labels = "80°C", cex.axis = 2, tick = F)
axis(side = 1, at = c(4), labels = "Smoke", cex.axis = 2, tick = F)
axis(side = 1, at = c(6), labels = "95°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
mtext(side = 1, expression(bold("Treatment")), line = 4.5, cex = 1.8)
mtext(side = 2, expression(bold("Days to 50% germination")), line = 3, cex = 1.8)
arrows(c(1:6), nt50_m4_ta$lci, c(1:6), nt50_m4_ta$uci, length = 0.05, code = 3, angle = 90, col = "#D94801")
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_m4_tl$fit, pch = 20, col = "#FD8D3C", cex = 2)
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_m4_tl$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_m4_tl$uci, length = 0.05, code = 3, angle = 90, col = "#FD8D3C")
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_m4_th$fit, pch = 20, col = "#7F2704", cex = 2)
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_m4_th$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_m4_th$uci, length = 0.05, code = 3, angle = 90, col = "#7F2704")
legend(x = 0.9, y = 18.5, legend = c("2.3 mg", "4.5 mg", "7 mg"), col = c("#FD8D3C", "#D94801", "#7F2704"), title = expression(bold("Seed weight")), lty = 1, pch = 19, cex = 1.8, bty = "n")
mtext("(b) ", cex = 2, adj = 0.001, line = 0.3)


plot.default(nt50_m5_ta$Treatment, nt50_m5_ta$fit, pch = 20, type = 'p', xlab  = "", ylab = "", xaxt = "n", las = 1, cex.axis = 2, ylim = c(4, 18), col = 'gray36', cex = 2)
axis(side = 1, at = c(1:6), labels = c("Control", "80°C", "95°C", "Smoke", "", ""), cex.axis = 2)
axis(side = 1, at = c(5), labels = "80°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
axis(side = 1, at = c(6), labels = "95°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
axis(side = 1, at = c(2), labels = "80°C", cex.axis = 2, tick = F)
axis(side = 1, at = c(4), labels = "Smoke", cex.axis = 2, tick = F)
mtext(side = 1, expression(bold("Treatment")), line = 4.5, cex = 1.8)
mtext(side = 2, expression(bold("Days to 50% germination")), line = 3, cex = 1.8)
arrows(c(1:6), nt50_m5_ta$lci, c(1:6), nt50_m5_ta$uci, length = 0.05, code = 3, angle = 90, col = 'gray36')
points(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_m5_tl$fit, pch = 20, col = "blue", cex = 2)
arrows(c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_m5_tl$lci,c(0.9, 1.9, 2.9, 3.9, 4.9, 5.9), nt50_m5_tl$uci, length = 0.05, code = 3, angle = 90, col = 'blue')
points(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_m5_th$fit, pch = 20, col = 'red', cex = 2)
arrows(c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_m5_th$lci, c(1.1, 2.1, 3.1, 4.1, 5.1, 6.1), nt50_m5_th$uci, length = 0.05, code = 3, angle = 90, col = 'red')
legend(x = 0.9, y = 18.5, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", 'gray36', 'red'), title = expression(bold("Fire frequency")), lty = 1, pch = 19, cex = 1.8, bty = "n")
mtext("(c) ", cex = 2, adj = 0.001, line = 0.3)



# Littoralis
plot.default(c(1:6), nt50_m1_l$fit, type = 'p', pch = 20, xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 2, ylim = c(4,18), cex = 2)
axis(side = 1, at = c(1:6), labels = c("Control", "80°C", "95°C", "Smoke", "", ""), cex.axis = 2)
axis(side = 1, at = c(5), labels = "80°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
axis(side = 1, at = c(6), labels = "95°C\n &\n smoke", cex.axis = 2, line = 4.5, tick = F)
axis(side = 1, at = c(2), labels = "80°C", cex.axis = 2, tick = F)
axis(side = 1, at = c(4), labels = "Smoke", cex.axis = 2, tick = F)
mtext(side = 1, expression(bold("Treatment")), line = 4.5, cex = 1.8)
mtext(side = 2, expression(bold("Days to 50% germination")), line = 3, cex = 1.8)
arrows(c(1:6), nt50_m1_l$lci, c(1:6), nt50_m1_l$uci, length = 0.05, code = 3, angle = 90)
mtext("(a) "~italic(Allocasuarina~littoralis), cex = 2, adj = 0.001)


# Plot effect sizes for t50 ----

# Create coefficient tables for plotting
# Treatment only m1
treat_coefs <- data.frame(
  Estimate = c(summary(m1_tt50)$coefficients[,1]),
  SE = c(summary(m1_tt50)$coefficients[,2]),
  Term = c("Intercept", "80°C", '95°C', "Smoke", "80°C+smoke", "90°C+smoke")
)
treat_coefs$lci <- treat_coefs$Estimate - (treat_coefs$SE * 1.96)
treat_coefs$uci <- treat_coefs$Estimate + (treat_coefs$SE * 1.96)
rownames(treat_coefs) <- treat_coefs$Term


# Treatment and seed weight m4
swt_coefs <- data.frame(
  Estimate = c(summary(m4_tt50)$coefficients[,1]),
  SE = c(summary(m4_tt50)$coefficients[,2]),
  Term = c("Intercept", "80°C", '95°C', "Smoke", "80°C+smoke", "95°C+smoke", "seed_wt", "80°C:seed_wt", "95°C:seed_wt", "Smoke:seed_wt", "80°C+smoke:seed_wt", "95°C+smoke:seed_wt")
)
swt_coefs$lci <- swt_coefs$Estimate - (swt_coefs$SE * 1.96)
swt_coefs$uci <- swt_coefs$Estimate + (swt_coefs$SE * 1.96)
rownames(swt_coefs) <- swt_coefs$Term


# Treatment and fire frequency m5
ff_coefs <- data.frame(
  Estimate = c(summary(m5_tt50)$coefficients[,1]),
  SE = c(summary(m5_tt50)$coefficients[,2]),
  Term = c("Intercept", "80°C", '95°C', "Smoke", "80°C+smoke", "95°C+smoke", "Fire_freq", "80°C:Fire_freq", "95°C:Fire_freq", "Smoke:Fire_freq", "80°C+smoke:Fire_freq", "95°C+smoke:Fire_freq")
)
ff_coefs$lci <- ff_coefs$Estimate - (ff_coefs$SE * 1.96)
ff_coefs$uci <- ff_coefs$Estimate + (ff_coefs$SE * 1.96)
rownames(ff_coefs) <- ff_coefs$Term


dev.new(height=10,width=12,dpi=80,pointsize=12,noRStudioGD = T)
par(mfrow=c(2,2),mar=c(5,12,2,3),mgp=c(2.7,1,0), cex = 1, cex.axis = 1.3, cex.lab = 1.75)

# Treatment only
plot(treat_coefs$Estimate, rev(1:nrow(treat_coefs)), xlim = c(min(treat_coefs$lci), max(treat_coefs$uci)), las = 1, cex = 1.8, ylab = "", xlab = expression(bold("Effect size")), pch = 20, yaxt = 'n', col = "black")
axis(side = 2, at = rev(1:nrow(treat_coefs)), labels = rownames(treat_coefs), las = 1)
arrows(treat_coefs$uci, rev(1:nrow(treat_coefs)), treat_coefs$lci, rev(1:nrow(treat_coefs)), code = 0, lwd = 0.8)
arrows(0,0,0, 7, code = 0, lwd = 0.8)
mtext('(a)', adj = 0000, line = 0.2, cex = 1.75)

plot.new()
plot(swt_coefs$Estimate, rev(1:nrow(swt_coefs)), xlim = c(min(swt_coefs$lci), max(swt_coefs$uci)), las = 1, cex = 1.8, ylab = "", xlab =  expression(bold("Effect size")), pch = 20, yaxt = 'n', col = "black")
axis(side = 2, at = rev(1:nrow(swt_coefs)), labels = rownames(swt_coefs), las = 1)
arrows(swt_coefs$uci, rev(1:nrow(swt_coefs)), swt_coefs$lci, rev(1:nrow(swt_coefs)), code = 0, lwd = 0.8)
arrows(0,0,0, 13, code = 0, lwd = 0.8)
mtext('(b)', adj = 0000, line = 0.2, cex = 1.75)


plot(ff_coefs$Estimate, rev(1:nrow(ff_coefs)), xlim = c(min(ff_coefs$lci), max(ff_coefs$uci)), las = 1, cex = 1.8, ylab = "", xlab =  expression(bold("Effect size")), pch = 20, yaxt = 'n', col = "black")
axis(side = 2, at = rev(1:nrow(ff_coefs)), labels = rownames(ff_coefs), las = 1)
arrows(ff_coefs$uci, rev(1:nrow(ff_coefs)), ff_coefs$lci, rev(1:nrow(ff_coefs)), code = 0, lwd = 0.8)
arrows(0,0,0, 13, code = 0, lwd = 0.8)
mtext('(c)', adj = 0000, line = 0.2, cex = 1.75)

# This would indicate that seed weight and fire frequency are uninformative, with the interactive effects of these terms not providing any useful information. 

# Investigate least-square means -----
sw_ls <- emmeans(m4_tt50, c("Treatment", 'r_seed_wt'), at = list(r_seed_wt = c(-2.14103, 0, 2.52076)), infer = c(T, T))
sw_ls
sw_ls_minmn <- emmeans(m4_tt50, c("Treatment", 'r_seed_wt'), at = list(r_seed_wt = c(-2.14103, 0)), infer = c(T, T))
cont_sw_minmn <- contrast(sw_ls_minmn, 'trt.vs.ctrl', infer = c(T,T))
sw_ls_minmx <- emmeans(m4_tt50, c("Treatment", 'r_seed_wt'), at = list(r_seed_wt = c(2.52076, -2.14103)), infer = c(T, T))
cont_sw_minmx <- contrast(sw_ls_minmx, 'trt.vs.ctrl', infer = c(T,T))
sw_ls_mnmx <- emmeans(m4_tt50, c("Treatment", 'r_seed_wt'), at = list(r_seed_wt = c(0, 2.52076)), infer = c(T, T))
cont_sw_mnmx <- contrast(sw_ls_mnmx, 'trt.vs.ctrl', infer = c(T,T))


# Create new dataframes 
contsw_minmn <- as.data.frame(1:11)
contsw_minmn$`1:11` <- c("80°C Min seed wt - Control Min seed wt", "95°C Min seed wt - Control Min seed wt", "Smoke Min seed wt - Control Min seed wt", "80°C+smoke Min seed wt - Control Min seed wt", "95°C+smoke Min seed wt - Control Min seed wt",
                         "Control Mean seed wt - Control Min seed wt", "80°C Mean seed wt - Control Min seed wt", "95°C Mean seed wt - Control Min seed wt", "Smoke Mean seed wt - Control Min seed wt", "80°C+smoke Mean seed wt - Control Min seed wt", "95°C+smoke Mean seed wt - Control Min seed wt")
colnames(contsw_minmn) <- "comp"
contsw_minmn$estimate <- summary(cont_sw_minmn)$estimate
contsw_minmn$se <- summary(cont_sw_minmn)$SE
contsw_minmn$lci <- contsw_minmn$estimate - (contsw_minmn$se * 1.96)
contsw_minmn$uci <- contsw_minmn$estimate + (contsw_minmn$se * 1.96)
contsw_minmn$lab <- c("80°C:Control 2.36 mg", "95°C:Control 2.36 mg", "Smoke:Control 2.36 mg", "80°C+smoke:Control 2.36 mg", "95°C+smoke:Control 2.36 mg",
                         "Control:Control 2.36:4.51 mg", "80°C:Control 2.36:4.51 mg", "95°C:Control 2.36:4.51 mg", "Smoke:Control 2.36:4.51 mg", "80°C+smoke:Control 2.36:4.51 mg", "95°C+smoke:Control 2.36:4.51 mg")

contsw_minmx <- as.data.frame(1:11)
contsw_minmx$`1:11` <- c("80°C Max seed wt - Control Max seed wt", "95°C Max seed wt - Control Max seed wt", "Smoke Max seed wt - Control Max seed wt", "80°C+smoke Max seed wt - Control Max seed wt", "95°C+smoke Max seed wt - Control Max seed wt",
                         "Control Max seed wt - Control ", "80°C Max seed wt - Control ", "95°C Max seed wt - Control ", "Smoke Max seed wt - Control ", "80°C+smoke Max seed wt - Control ", "95°C+smoke Max seed wt - Control ")
colnames(contsw_minmx) <- "comp"
contsw_minmx$estimate <- summary(cont_sw_minmx)$estimate
contsw_minmx$se <- summary(cont_sw_minmx)$SE
contsw_minmx$lci <- contsw_minmx$estimate - (contsw_minmx$se * 1.96)
contsw_minmx$uci <- contsw_minmx$estimate + (contsw_minmx$se * 1.96)
contsw_minmx$lab <- c("80°C:Control 7.04 mg", "95°C:Control 7.04 mg", "Smoke:Control 7.04 mg", "80°C+smoke:Control 7.04 mg", "95°C+smoke:Control 7.04 mg",
                      "Control:Control 2.36:7.04 mg", "80°C:Control 2.36:7.04 mg", "95°C:Control 2.36:7.04 mg", "Smoke:Control 2.36:7.04 mg", "80°C+smoke:Control 2.36:7.04 mg", "95°C+smoke:Control 2.36:7.04 mg")

contsw_mnmx <- as.data.frame(1:11)
contsw_mnmx$`1:11` <- c("80°C Mean seed wt - Control Mean seed wt", "95°C Mean seed wt - Control Mean seed wt", "Smoke Mean seed wt - Control Mean seed wt", "80°C+smoke Mean seed wt - Control Mean seed wt", "95°C+smoke Mean seed wt - Control Mean seed wt",
                         "Control Max seed wt - Control Mean seed wt", "80°C Max seed wt - Control Mean seed wt", "95°C Max seed wt - Control Mean seed wt", "Smoke Max seed wt - Control Mean seed wt", "80°C+smoke Max seed wt - Control Mean seed wt", "95°C+smoke Max seed wt - Control Mean seed wt")
colnames(contsw_mnmx) <- "comp"
contsw_mnmx$estimate <- summary(cont_sw_mnmx)$estimate
contsw_mnmx$se <- summary(cont_sw_mnmx)$SE
contsw_mnmx$lci <- contsw_mnmx$estimate - (contsw_mnmx$se * 1.96)
contsw_mnmx$uci <- contsw_mnmx$estimate + (contsw_mnmx$se * 1.96)
contsw_mnmx$lab <- c("80°C:Control 4.51 mg", "95°C:Control 4.51 mg", "Smoke:Control 4.51 mg", "80°C+smoke:Control 4.51 mg", "95°C+smoke:Control 4.51:7.04 mg",
                      "Control:Control 4.51:7.04 mg", "80°C:Control 4.51:7.04 mg", "95°C:Control 4.51:7.04 mg", "Smoke:Control 4.51:7.04 mg", "80°C+smoke:Control 4.51:7.04 mg", "95°C+smoke:Control 4.51:7.04 mg")


m5_tt50_ls_minmn<- emmeans(m5_tt50, c("Treatment", 'r_fire_freq'), at = list(r_fire_freq = c(-1.91386, 0)))# Specify at to override default reference grid with fire frequency only kept at the mean. Use by to ensure we get separate sets of predictions for each fire frequency category
m5_tt50_ls_minmn
cont_m5_minmn <- contrast(m5_tt50_ls_minmn, "trt.vs.ctrl", infer = c(T,T))
m5_tt50_ls_minmx <- emmeans(m5_tt50, c("Treatment", 'r_fire_freq'), at = list(r_fire_freq = c(2.09967, -1.91386))) 
m5_tt50_ls_minmx
cont_m5_minmx <- contrast(m5_tt50_ls_minmx, "trt.vs.ctrl", infer = c(T,T))
m5_tt50_ls_mnmx <- emmeans(m5_tt50, c("Treatment", 'r_fire_freq'), at = list(r_fire_freq = c(0,  2.09967))) 
m5_tt50_ls_mnmx
cont_m5_mnmx <- contrast(m5_tt50_ls_mnmx, 'trt.vs.ctrl', infer = c(T,T))


# Create new dataframe containing the comparisons we want
cont_m5_tt50 <- as.data.frame(1:11)
cont_m5_tt50$`1:11` <- c("80°C 0 fires - Control 0 fires", "95°C 0 fires - Control 0 fires", "Smoke 0 fires - Control 0 fires", "80°C+smoke 0 fires - Control 0 fires", "95°C+smoke 0 fires - Control 0 fires",
                                "Control 3 fires - Control 0 fires", "80°C 3 fires - Control 0 fires", "95°C 3 fires - Control 0 fires", "Smoke 3 fires - Control 0 fires", "80°C+smoke 3 fires - Control 0 fires", "95°C+smoke 3 fires - Control 0 fires")
colnames(cont_m5_tt50) <- "comparison"
cont_m5_tt50$estimate <- summary(cont_m5_minmn)$estimate
cont_m5_tt50$se <- summary(cont_m5_minmn)$SE
cont_m5_tt50$lci <- cont_m5_tt50$estimate - (cont_m5_tt50$se * 1.96)
cont_m5_tt50$uci <- cont_m5_tt50$estimate + (cont_m5_tt50$se * 1.96)
cont_m5_tt50$lab <- c("80°C:Control 0 fires", "95°C:Control 0 fires", "Smoke:Control 0 fires", "80°C+smoke:Control 0 fires", "95°C+smoke:Control 0 fires",
                     "Control:Control 0:3 fires", "80°C:Control 0:3 fires", "95°C:Control 0:3 fires", "Smoke:Control 0:3 fires", "80°C+smoke:Control 0:3 fires", "95°C+smoke:Control 0:3 fires")



cont_m5_tt50_1 <- as.data.frame(1:11)
cont_m5_tt50_1$`1:11` <- c("80°C 6 fires - Control 6 fires", "95°C 6 fires - Control 6 fires", "Smoke 6 fires - Control 6 fires", "80°C+smoke 6 fires - Control 6 fires", "95°C+smoke 6 fires - Control 6 fires",
                         "Control 0 fires - Control 6 fires", "80°C 0 fires - Control 6 fires", "95°C 0 fires - Control 6 fires", "Smoke 0 fires - Control 6 fires", "80°C+smoke 0 fires - Control 6 fires", "95°C+smoke 0 fires - Control 6 fires")
colnames(cont_m5_tt50_1) <- "comparison"
cont_m5_tt50_1$estimate <- summary(cont_m5_minmx)$estimate
cont_m5_tt50_1$se <- summary(cont_m5_minmx)$SE
cont_m5_tt50_1$lci <- cont_m5_tt50_1$estimate - (cont_m5_tt50_1$se * 1.96)
cont_m5_tt50_1$uci <- cont_m5_tt50_1$estimate + (cont_m5_tt50_1$se * 1.96)
cont_m5_tt50_1$lab <- c("80°C:Control 6 fires", "95°C:Control 6 fires", "Smoke:Control 6 fires", "80°C+smoke:Control 6 fires", "95°C+smoke:Control 6 fires",
                      "Control:Control 0:6 fires", "80°C:Control 0:6 fires", "95°C:Control 0:6 fires", "Smoke:Control 0:6 fires", "80°C+smoke:Control 0:6 fires", "95°C+smoke:Control 0:6 fires")

cont_m5_tt50_2 <- as.data.frame(1:11)
cont_m5_tt50_2$`1:11` <- c("80°C 3 fires - Control 3 fires", "95°C 3 fires - Control 3 fires", "Smoke 3 fires - Control 3 fires", "80°C+smoke 3 fires - Control 3 fires", "95°C+smoke 3 fires - Control 3 fires",
                           "Control 6 fires - Control 3 fires", "80°C 6 fires - Control 3 fires", "95°C 6 fires - Control 3 fires", "Smoke 6 fires - Control 3 fires", "80°C+smoke 6 fires - Control 3 fires", "95°C+smoke 6 fires - Control 3 fires")
colnames(cont_m5_tt50_2) <- "comparison"
cont_m5_tt50_2$estimate <- summary(cont_m5_mnmx)$estimate
cont_m5_tt50_2$se <- summary(cont_m5_mnmx)$SE
cont_m5_tt50_2$lci <- cont_m5_tt50_2$estimate - (cont_m5_tt50_2$se * 1.96)
cont_m5_tt50_2$uci <- cont_m5_tt50_2$estimate + (cont_m5_tt50_2$se * 1.96)
cont_m5_tt50_2$lab <- c("80°C:Control 3 fires", "95°C:Control 3 fires", "Smoke:Control 3 fires", "80°C+smoke:Control 3 fires", "95°C+smoke:Control 3 fires",
                      "Control:Control 3:6 fires", "80°C:Control 3:6 fires", "95°C:Control 3:6 fires", "Smoke:Control 3:6 fires", "80°C+smoke:Control 3:6 fires", "95°C+smoke:Control 3:6 fires")




# Plot least square means ----
dev.new(height=15,width=15,dpi=80,pointsize=12,noRStudioGD = T)
par(mfrow=c(3,2),mar=c(4,12,5,2),mgp=c(2.7,1,0), cex = 1, cex.axis = 1, cex.lab = 1.5, oma = c(0,2,0,0))

plot(contsw_minmn$estimate, rev(1:nrow(contsw_minmn)), xlim = c(-0.4,0.9), las = 1, cex = 1.8, ylab = "", xlab = expression(bold("Least squares mean difference")), pch = 20, yaxt = "n", col = "black")
axis(side = 2, at = rev(1:nrow(contsw_minmn)), labels = contsw_minmn$lab, las = 1)
arrows(contsw_minmn$uci, rev(1:nrow(contsw_minmn)), contsw_minmn$lci, rev(1:nrow(contsw_minmn)), code = 0, lwd = 0.8)
arrows(0,0,0, 12, code = 0, lwd = 0.8)
mtext("(a) 2.36 mg compared to 4.51 mg", adj = 1, line = 0.2, cex = 1.5)
mtext(expression(bold("Seed weight")), cex = 2, line = 2, adj = 0.5)


plot(cont_m5_tt50$estimate, rev(1:nrow(cont_m5_tt50)), xlim = c(-0.4,0.9), las = 1, cex = 1.8, ylab = "", xlab = expression(bold("Least squares mean difference")), pch = 20, yaxt = "n", col = "black")
axis(side = 2, at = rev(1:nrow(cont_m5_tt50)), labels = cont_m5_tt50$lab, las = 1)
arrows(cont_m5_tt50$uci, rev(1:nrow(cont_m5_tt50)), cont_m5_tt50$lci, rev(1:nrow(cont_m5_tt50)), code = 0, lwd = 0.8)
arrows(0,0,0,12, code = 0, lwd = 0.8)
mtext('(d) 0 fires compared to 3 fires', adj = -1, line = 0.2, cex = 1.5)
mtext(expression(bold("Fire frequency")), cex = 2, line = 2, adj = 0.5)




plot(contsw_minmx$estimate, rev(1:nrow(contsw_minmx)), xlim = c(-0.4,0.9), las = 1, cex = 1.8, ylab = "", xlab = expression(bold("Least squares mean difference")), pch = 20, yaxt = "n", col = 'black')
axis(side = 2, at = rev(1:nrow(contsw_minmx)), labels = contsw_minmx$lab, las = 1)
arrows(contsw_minmx$uci, rev(1:nrow(contsw_minmx)), contsw_minmx$lci, rev(1:nrow(contsw_minmx)), code = 0, lwd = 0.8)
arrows(0,0,0, 12, code = 0, lwd = 0.8)
mtext("(b) 2.36 mg compared to 7.04 mg", adj = 1, line = 0.2, cex = 1.5)


plot(cont_m5_tt50_1$estimate, rev(1:nrow(cont_m5_tt50_1)), xlim = c(-0.4,0.9), las = 1, cex = 1.8, ylab = "", xlab = expression(bold("Least squares mean difference")), pch = 20, yaxt = "n", col = "black")
axis(side = 2, at = rev(1:nrow(cont_m5_tt50_1)), labels = cont_m5_tt50_1$lab, las = 1)
arrows(cont_m5_tt50_1$uci, rev(1:nrow(cont_m5_tt50_1)), cont_m5_tt50_1$lci, rev(1:nrow(cont_m5_tt50_1)), code = 0, lwd = 0.8)
arrows(0,0,0,12, code = 0, lwd = 0.8)
mtext('(e) 0 fires compared to 6 fires', adj = -1, line = 0.2, cex = 1.5)





plot(contsw_mnmx$estimate, rev(1:nrow(contsw_mnmx)), xlim = c(-0.4,0.9), las = 1, cex = 1.8, ylab = "", xlab = expression(bold("Least squares mean difference")), pch = 20, yaxt = "n", col = 'black')
axis(side = 2, at = rev(1:nrow(contsw_mnmx)), labels = contsw_mnmx$lab, las = 1)
arrows(contsw_mnmx$uci, rev(1:nrow(contsw_mnmx)), contsw_mnmx$lci, rev(1:nrow(contsw_mnmx)), code = 0, lwd = 0.8)
arrows(0,0,0, 12, code = 0, lwd = 0.8)
mtext("(c) 4.51 mg compared to 7.04 mg", adj = 1, line = 0.2, cex = 1.5)



plot(cont_m5_tt50_2$estimate, rev(1:nrow(cont_m5_tt50_2)), xlim = c(-0.4,0.9), las = 1, cex = 1.8, ylab = "", xlab = expression(bold("Least squares mean difference")), pch = 20, yaxt = "n", col = "black")
axis(side = 2, at = rev(1:nrow(cont_m5_tt50_2)), labels = cont_m5_tt50_2$lab, las = 1)
arrows(cont_m5_tt50_2$uci, rev(1:nrow(cont_m5_tt50_2)), cont_m5_tt50_2$lci, rev(1:nrow(cont_m5_tt50_2)), code = 0, lwd = 0.8)
arrows(0,0,0,12, code = 0, lwd = 0.8)
mtext('(f) 3 fires compared to 6 fires', adj = -1, line = 0.2, cex = 1.5)








# QUESTION 2 : How does recent fire activity influence population age structure and number of cones? ----
# 6.1 Proportion of seedlings to adults ----
rect_null <- glmer(Proportion_seedlings ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
recl_null <- glmer(Proportion_seedlings ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)


rec_m1t <- glmer(Proportion_seedlings ~ r_TSF + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(rec_m1t)
rec_m1l <- glmer(Proportion_seedlings ~ r_TSF + (1 | Location/Transect), family = binomial, data = lit_transects)
summary(rec_m1l)


rec_m2t <- glmer(Proportion_seedlings ~ r_TSF * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rec_m2t)
rec_m2l <- glmer(Proportion_seedlings ~ r_TSF * r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rec_m2l) 


rec_m3t <- glmer(Proportion_seedlings ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rec_m3t)
rec_m3l <- glmer(Proportion_seedlings ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rec_m3l)


rec_m4t <- glmer(Proportion_seedlings ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rec_m4t)
rec_m4l <- glmer(Proportion_seedlings ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rec_m4l)


rec_m5t <- glmer(Proportion_seedlings ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rec_m5t)
rec_m5l <- glmer(Proportion_seedlings ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=40)))
summary(rec_m5l)


rec_m6t <- glmer(Proportion_seedlings ~ r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rec_m6t)
rec_m6l <- glmer(Proportion_seedlings ~ r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rec_m6l)


rec_m7t <- glmer(Proportion_seedlings ~ r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rec_m7t)
rec_m7l <- glmer(Proportion_seedlings ~ r_Temp + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rec_m7l)


rec_m8t <- glmer(Proportion_seedlings ~ r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rec_m8t)
rec_m8l <- glmer(Proportion_seedlings ~ r_FPC + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rec_m8l)


rec_m9t <- glmer(Proportion_seedlings ~ r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rec_m9t)
rec_m9l <- glmer(Proportion_seedlings ~ r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=140)))
summary(rec_m9l)


rec_m10t <- glmer(Proportion_seedlings ~ r_TSF * r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rec_m10t)
rec_m10l <- glmer(Proportion_seedlings ~ r_TSF * r_TWI + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rec_m10l)


rec_m11t <- glmer(Proportion_seedlings ~ r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rec_m11t)
rec_m11l <- glmer(Proportion_seedlings ~ r_TWI + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rec_m11l)


# Compare models
aictab(cand.set = c(rect_null, rec_m1t, rec_m2t, rec_m3t, rec_m4t, rec_m5t, rec_m6t, rec_m7t, rec_m8t, rec_m9t, rec_m10t, rec_m11t), modnames = c("Null", "TSF", "TSF*Lat", "TSF*FPC", "TSF*Precip", "TSF*Temp", "Precip", "Temp", "FPC", "Lat", "TSF*TWI", "TWI"))  
# M8 is best. Proportion of seedlings is significantly influenced by foliage projective cover and climatic controls rather than time since fire, which is ranked lower than these univariate environmental models.

aictab(cand.set = c(recl_null, rec_m1l, rec_m2l, rec_m3l, rec_m4l, rec_m5l, rec_m6l, rec_m7l, rec_m8l, rec_m9l, rec_m10l, rec_m11l), modnames = c("Null", "TSF", "TSF*Lat", "TSF*FPC", "TSF*Precip", "TSF*Temp", "Precip", "Temp", "FPC", "Lat", "TSF*TWI", "TWI")) 
# Null model is best but temperature seasonality has a stronger influence on proportion of seedlings for torulosa than time since fire (Mod 2). 


# Predict to new data for proportion seedlings 
tor_seedtsf <- data.frame(r_FPC = seq(min(tor_transects$r_FPC, na.rm = T), max(tor_transects$r_FPC, na.rm = T), length = 50))
pt_seedtsf <- predictSE(rec_m8t, tor_seedtsf, se.fit = T, type = 'link')
tor_seedtsf$fit.link <- pt_seedtsf$fit
tor_seedtsf$se.link <- pt_seedtsf$se.fit
tor_seedtsf$lci.link <- tor_seedtsf$fit.link - (tor_seedtsf$se.link * 1.96)
tor_seedtsf$uci.link <- tor_seedtsf$fit.link + (tor_seedtsf$se.link * 1.96)
tor_seedtsf$fit <- invlogit(tor_seedtsf$fit.link)
tor_seedtsf$se <- invlogit(tor_seedtsf$se.link)
tor_seedtsf$lci <- invlogit(tor_seedtsf$lci.link)
tor_seedtsf$uci <- invlogit(tor_seedtsf$uci.link)


plot(tor_seedtsf$r_FPC, tor_seedtsf$fit, type = 'l', ylim = c(0,1))
pg.ci(x = 'r_FPC', data = 'tor_seedtsf', colour = rgb(0, 0, 0, 0.1), lower = 'lci', upper = 'uci')

# 6.2 Proportion of saplings to adults ----
rsap_tnull <- glmer(Proportion_saplings ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
rsap_lnull <- glmer(Proportion_saplings ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)


rsap_m1t <- glmer(Proportion_saplings ~ r_TSF + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m1t)
rsap_m1l <- glmer(Proportion_saplings ~ r_TSF + (1 | Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))# Not a very good model
summary(rsap_m1l) 


rsap_m2t <- glmer(Proportion_saplings ~ r_TSF * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m2t)
rsap_m2l <- glmer(Proportion_saplings ~ r_TSF * r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=40))) # Not very good model, as before the low number of samples means we can only model this additively
summary(rsap_m2l) 


rsap_m3t <- glmer(Proportion_saplings ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m3t)
rsap_m3l <- glmer(Proportion_saplings ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=200)))
summary(rsap_m3l)


rsap_m4t <-  glmer(Proportion_saplings ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m4t)
rsap_m4l <-  glmer(Proportion_saplings ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=130)))
summary(rsap_m4l)


rsap_m5t <-  glmer(Proportion_saplings ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m5t)
rsap_m5l <-  glmer(Proportion_saplings ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rsap_m5l)


rsap_m6t <- glmer(Proportion_saplings ~ r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m6t)
rsap_m6l <- glmer(Proportion_saplings ~ r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rsap_m6l)


rsap_m7t <- glmer(Proportion_saplings ~ r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m7t)
rsap_m7l <- glmer(Proportion_saplings ~ r_Temp + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl = list(maxfun = 30)))
summary(rsap_m7l)


rsap_m8t <- glmer(Proportion_saplings ~ r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m8t)
rsap_m8l <- glmer(Proportion_saplings ~ r_FPC + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl = list(maxfun = 30)))
summary(rsap_m8l)


rsap_m9t <- glmer(Proportion_saplings ~ r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m9t)
rsap_m9l <- glmer(Proportion_saplings ~ r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl = list(maxfun = 30)))
summary(rsap_m9l)


rsap_m10t <- glmer(Proportion_saplings ~ r_TSF * r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m10t)
rsap_m10l <- glmer(Proportion_saplings ~ r_TSF * r_TWI + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl = list(maxfun = 50)))
summary(rsap_m10l)


rsap_m11t <- glmer(Proportion_saplings ~ r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m11t)
rsap_m11l <- glmer(Proportion_saplings ~ r_TWI + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rsap_m11l)

# Compare models
aictab(cand.set = c(rsap_tnull, rsap_m1t, rsap_m2t, rsap_m3t, rsap_m4t, rsap_m5t, rsap_m6t, rsap_m7t, rsap_m8t, rsap_m9t, rsap_m10t, rsap_m11t), modnames = c("Null", "TSF", "TSF*Lat", "TSF*FPC", "TSF*Precip", "TSF*Temp", "Precip", "Temp", "FPC", "Lat", "TSF*TWI", "TWI")) 
# Model 7 is best, so temperature seasonality has a stronger influence on proportion of saplings than time since fire

aictab(cand.set = c(rsap_lnull, rsap_m1l, rsap_m2l, rsap_m3l, rsap_m4l, rsap_m5l, rsap_m6l, rsap_m7l, rsap_m8l, rsap_m9l, rsap_m10l, rsap_m11l), modnames = c("Null", "TSF", "TSF*Lat", "TSF*FPC", "TSF*Precip", "TSF*Temp", "Precip", "Temp", "FPC", "Lat", "TSF*TWI", "TWI")) 
# Null model is best, but temperature seasonality and latityde have a stronger influence on proportion of saplings than time since fire.


# 6.3 Proportion of recruits ----
recruit_tnull <- glmer(Proportion_recruits ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
recruit_lnull <- glmer(Proportion_recruits ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)


recruit_m1t <- glmer(Proportion_recruits ~ r_TSF + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m1t)
recruit_m1l <- glmer(Proportion_recruits ~ r_TSF + (1 | Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30))) # Not a very good model
summary(recruit_m1l)


recruit_m2t <- glmer(Proportion_recruits ~ r_TSF * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m2t)
summary(recruit_m2t)
recruit_m2l <- glmer(Proportion_recruits ~ r_TSF * r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=40))) # Not a very good model 
summary(recruit_m2l) 

recruit_m3t <- glmer(Proportion_recruits ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m3t)
recruit_m3l <- glmer(Proportion_recruits ~ r_TSF * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30))) # Not a very good model
summary(recruit_m3l)


recruit_m4t <- glmer(Proportion_recruits ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m4t)
recruit_m4l <- glmer(Proportion_recruits ~ r_TSF * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=250)))
summary(recruit_m4l)


recruit_m5t <- glmer(Proportion_recruits ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m5t)
recruit_m5l <- glmer(Proportion_recruits ~ r_TSF * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=250)))
summary(recruit_m5l) 


recruit_m6t <- glmer(Proportion_recruits ~ r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m6t)
recruit_m6l <- glmer(Proportion_recruits ~ r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(recruit_m6l)


recruit_m7t <- glmer(Proportion_recruits ~ r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m7t)
recruit_m7l <- glmer(Proportion_recruits ~ r_Temp + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl = list(maxfun = 30)))
summary(recruit_m7l)


recruit_m8t <- glmer(Proportion_recruits ~ r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m8t)
recruit_m8l <- glmer(Proportion_recruits ~ r_FPC + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl = list(maxfun = 20)))
summary(recruit_m8l)


recruit_m9t <- glmer(Proportion_recruits ~ r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m9t)
recruit_m9l <- glmer(Proportion_recruits ~ r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl = list(maxfun = 30)))
summary(recruit_m9l)


recruit_m10t <- glmer(Proportion_recruits ~ r_TSF * r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m10t)
recruit_m10l <- glmer(Proportion_recruits ~ r_TSF * r_TWI + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=50)))
summary(recruit_m10l)


recruit_m11t <- glmer(Proportion_recruits ~ r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recruit_m11t)
recruit_m11l <- glmer(Proportion_recruits ~ r_TWI + (1|Location/Transect), family = binomial, data = lit_transects)
summary(recruit_m11l)

# Compare models
aictab(cand.set = c(recruit_tnull, recruit_m1t, recruit_m2t, recruit_m3t, recruit_m4t, recruit_m5t, recruit_m6t, recruit_m7t, recruit_m8t, recruit_m9t, recruit_m10t, recruit_m11t), modnames = c("Null", "TSF", "TSF*Lat", "TSF*FPC", "TSF*Precip", "TSF*Temp", "Precip", "Temp", "FPC", "Lat", "TSF*TWI", "TWI")) 
# Model 7 is best, while model 5 is ranked within delta AICc <2, this includes TSF and temperature seasonality but temperature seasonality alone is a better predictor of proportion recruitment. 


aictab(cand.set = c(recruit_lnull, recruit_m1l, recruit_m2l, recruit_m3l, recruit_m4l, recruit_m5l, recruit_m6l, recruit_m7l, recruit_m8l, recruit_m9l, recruit_m10l, recruit_m11l), modnames = c("Null", "TSF", "TSF*Lat", "TSF*FPC", "TSF*Precip", "TSF*Temp", "Precip", "Temp", "FPC", "Lat", "TSF*TWI", "TWI")) 
# Null model is best but temperature seasonality is also a better predictor of proportion of recruits than time since fire. 


# Predict to new data for torulosa proportion recruits ~ TSF ----
rect_m7 <- data.frame(r_Temp = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 50))

rect_m5_l <- data.frame(r_Temp = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 50),
                        r_TSF = min(tor_transects$r_TSF))
rect_m5_a <- data.frame(r_Temp = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 50),
                        r_TSF = mean(tor_transects$r_TSF))
rect_m5_h <- data.frame(r_Temp = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 50),
                        r_TSF = max(tor_transects$r_TSF))


prtm7 <- predictSE(recruit_m7t, newdata = rect_m7, se.fit = T, type = 'link')
rect_m7$fit.link <- prtm7$fit
rect_m7$se.link <- prtm7$se.fit
rect_m7$lci.link <- rect_m7$fit.link - (rect_m7$se.link * 1.96)
rect_m7$uci.link <- rect_m7$fit.link + (rect_m7$se.link * 1.96)
rect_m7$lci85.link <- rect_m7$fit.link - (rect_m7$se.link * 1.440)
rect_m7$uci85.link <- rect_m7$fit.link + (rect_m7$se.link * 1.440)
rect_m7$fit <- invlogit(rect_m7$fit.link)
rect_m7$se <- invlogit(rect_m7$se.link)
rect_m7$lci <- invlogit(rect_m7$lci.link)
rect_m7$uci <- invlogit(rect_m7$uci.link)
rect_m7$lci85 <- invlogit(rect_m7$lci85.link)
rect_m7$uci85 <- invlogit(rect_m7$uci85.link)


prtm5_l <- predictSE(recruit_m5t, newdata = rect_m5_l, se.fit = T, type = 'link')
rect_m5_l$fit.link <- prtm5_l$fit
rect_m5_l$se.link <- prtm5_l$se.fit
rect_m5_l$lci.link <- rect_m5_l$fit.link - (rect_m5_l$se.link * 1.96)
rect_m5_l$uci.link <- rect_m5_l$fit.link + (rect_m5_l$se.link * 1.96)
rect_m5_l$lci85.link <- rect_m5_l$fit.link - (rect_m5_l$se.link * 1.440)
rect_m5_l$uci85.link <- rect_m5_l$fit.link + (rect_m5_l$se.link * 1.440)
rect_m5_l$fit <- invlogit(rect_m5_l$fit.link)
rect_m5_l$se <- invlogit(rect_m5_l$se.link)
rect_m5_l$lci <- invlogit(rect_m5_l$lci.link)
rect_m5_l$uci <- invlogit(rect_m5_l$uci.link)
rect_m5_l$lci85 <- invlogit(rect_m5_l$lci85.link)
rect_m5_l$uci85 <- invlogit(rect_m5_l$uci85.link)


prtm5_a <- predictSE(recruit_m5t, newdata = rect_m5_a, se.fit = T, type = 'link')
rect_m5_a$fit.link <- prtm5_a$fit
rect_m5_a$se.link <- prtm5_a$se.fit
rect_m5_a$lci.link <- rect_m5_a$fit.link - (rect_m5_a$se.link * 1.96)
rect_m5_a$uci.link <- rect_m5_a$fit.link + (rect_m5_a$se.link * 1.96)
rect_m5_a$lci85.link <- rect_m5_a$fit.link - (rect_m5_a$se.link * 1.440)
rect_m5_a$uci85.link <- rect_m5_a$fit.link + (rect_m5_a$se.link * 1.440)
rect_m5_a$fit <- invlogit(rect_m5_a$fit.link)
rect_m5_a$se <- invlogit(rect_m5_a$se.link)
rect_m5_a$lci <- invlogit(rect_m5_a$lci.link)
rect_m5_a$uci <- invlogit(rect_m5_a$uci.link)
rect_m5_a$lci85 <- invlogit(rect_m5_a$lci85.link)
rect_m5_a$uci85 <- invlogit(rect_m5_a$uci85.link)


prtm5_h <- predictSE(recruit_m5t, newdata = rect_m5_h, se.fit = T, type = 'link')
rect_m5_h$fit.link <- prtm5_h$fit
rect_m5_h$se.link <- prtm5_h$se.fit
rect_m5_h$lci.link <- rect_m5_h$fit.link - (rect_m5_h$se.link * 1.96)
rect_m5_h$uci.link <- rect_m5_h$fit.link + (rect_m5_h$se.link * 1.96)
rect_m5_h$lci85.link <- rect_m5_h$fit.link - (rect_m5_h$se.link * 1.440)
rect_m5_h$uci85.link <- rect_m5_h$fit.link + (rect_m5_h$se.link * 1.440)
rect_m5_h$fit <- invlogit(rect_m5_h$fit.link)
rect_m5_h$se <- invlogit(rect_m5_h$se.link)
rect_m5_h$lci <- invlogit(rect_m5_h$lci.link)
rect_m5_h$uci <- invlogit(rect_m5_h$uci.link)
rect_m5_h$lci85 <- invlogit(rect_m5_h$lci85.link)
rect_m5_h$uci85 <- invlogit(rect_m5_h$uci85.link)

# Plot predictions for torulosa proportion recruits ~ TSF ----
dev.new(width = 22, height = 8, dpi = 300, noRStudioGD = T)
par(mfrow = c(1,2), mar = c(4,5.5,2,2), oma = c(0, 0, 0, 10))

plot(tor_transects$r_Temp, tor_transects$Proportion_recruits, pch = 19, cex = 0.5, xlab = "", ylab = "", las = 1, cex.axis = 1.8, xaxt = "n", col = 'gray36', ylim = c(0,1))
axis(side = 1, at = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 10), labels = round(seq(min(tor_transects$Temp), max(tor_transects$Temp), length = 10), 0), cex.axis = 1.7)
axis(side = 1, at = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 10), labels = c("", "391", "", "399", "", "408", "", "416", "", "424"), cex.axis = 1.7)
lines(rect_m7$r_Temp, rect_m7$fit, col = 'gray36')
pg.ci(x = 'r_Temp', data = 'rect_m7', colour = rgb(92/255, 92/255, 92/255, 0.1), lower = 'lci', upper = 'uci')
mtext(side = 1, expression(bold("Temperature seasonality")), line = 3, cex = 2)
mtext(side = 2, expression(bold("Proportion recruits to adults")), cex = 2, line = 3.5)
mtext("(a) ", cex = 2, adj = 0.001, line = 0.2)


plot(tor_transects$r_Temp, tor_transects$Proportion_recruits, pch = 19, cex = 0.5, xlab = "", ylab = "", las = 1, cex.axis = 1.8, xaxt = "n", col = 'gray36', ylim = c(0,1))
axis(side = 1, at = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 10), labels = round(seq(min(tor_transects$Temp), max(tor_transects$Temp), length = 10), 0), cex.axis = 1.7)
axis(side = 1, at = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 10), labels = c("", "391", "", "399", "", "408", "", "416", "", "424"), cex.axis = 1.7)
lines(rect_m5_h$r_Temp, rect_m5_h$fit, col = 'firebrick')
pg.ci(x = 'r_Temp', data = 'rect_m5_h', colour = rgb(178/255, 34/255, 34/255, 0.1), lower = 'lci', upper = 'uci')
lines(rect_m5_a$r_Temp, rect_m5_a$fit, col = 'gray52')
pg.ci(x = 'r_Temp', data = 'rect_m5_a', colour = rgb(133/255, 133/255, 133/255, 0.2), lower = 'lci', upper = 'uci')
lines(rect_m5_l$r_Temp, rect_m5_l$fit, col = 'deepskyblue3')
pg.ci(x = 'r_Temp', data = 'rect_m5_l', colour = rgb(0, 154/255, 205/255, 0.1), lower = 'lci', upper = 'uci')
mtext(side = 1, expression(bold("Temperature seasonality")), line = 3, cex = 2)
mtext(side = 2, expression(bold("Proportion recruits to adults")), cex = 2, line = 3.5)
mtext("(b) ", cex = 2, adj = 0.001, line = 0.2)

par(xpd = NA)
legend(x = 1.4, y = 1, legend = c("1 year", "9 years", "25 years"), col = c("deepskyblue3", 'gray52', 'firebrick'), title = expression(bold("Time since \n fire")), lty = 1, cex = 2, bty = "n")
par(xpd = F)



# 85% CI
dev.new(width = 22, height = 8, dpi = 300, noRStudioGD = T)
par(mfrow = c(1,2), mar = c(4,5.5,2,2), oma = c(0, 0, 0, 10))

plot(tor_transects$r_Temp, tor_transects$Proportion_recruits, pch = 19, cex = 0.5, xlab = "", ylab = "", las = 1, cex.axis = 1.8, xaxt = "n", col = 'gray36', ylim = c(0,1))
axis(side = 1, at = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 10), labels = round(seq(min(tor_transects$Temp), max(tor_transects$Temp), length = 10), 0), cex.axis = 1.7)
axis(side = 1, at = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 10), labels = c("", "391", "", "399", "", "408", "", "416", "", "424"), cex.axis = 1.7)
lines(rect_m7$r_Temp, rect_m7$fit, col = 'gray36')
pg.ci(x = 'r_Temp', data = 'rect_m7', colour = rgb(92/255, 92/255, 92/255, 0.1), lower = 'lci85', upper = 'uci85')
mtext(side = 1, expression(bold("Temperature seasonality")), line = 3, cex = 2)
mtext(side = 2, expression(bold("Proportion recruits to adults")), cex = 2, line = 3.5)
mtext("(a) ", cex = 2, adj = 0.001, line = 0.2)


plot(tor_transects$r_Temp, tor_transects$Proportion_recruits, pch = 19, cex = 0.5, xlab = "", ylab = "", las = 1, cex.axis = 1.8, xaxt = "n", col = 'gray36', ylim = c(0,1))
axis(side = 1, at = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 10), labels = round(seq(min(tor_transects$Temp), max(tor_transects$Temp), length = 10), 0), cex.axis = 1.7)
axis(side = 1, at = seq(min(tor_transects$r_Temp), max(tor_transects$r_Temp), length = 10), labels = c("", "391", "", "399", "", "408", "", "416", "", "424"), cex.axis = 1.7)
lines(rect_m5_h$r_Temp, rect_m5_h$fit, col = 'firebrick')
pg.ci(x = 'r_Temp', data = 'rect_m5_h', colour = rgb(178/255, 34/255, 34/255, 0.1), lower = 'lci85', upper = 'uci85')
lines(rect_m5_a$r_Temp, rect_m5_a$fit, col = 'gray52')
pg.ci(x = 'r_Temp', data = 'rect_m5_a', colour = rgb(133/255, 133/255, 133/255, 0.2), lower = 'lci85', upper = 'uci85')
lines(rect_m5_l$r_Temp, rect_m5_l$fit, col = 'deepskyblue3')
pg.ci(x = 'r_Temp', data = 'rect_m5_l', colour = rgb(0, 154/255, 205/255, 0.1), lower = 'lci85', upper = 'uci85')
mtext(side = 1, expression(bold("Temperature seasonality")), line = 3, cex = 2)
mtext(side = 2, expression(bold("Proportion recruits to adults")), cex = 2, line = 3.5)
mtext("(b) ", cex = 2, adj = 0.001, line = 0.2)

par(xpd = NA)
legend(x = 1.4, y = 1, legend = c("1 year", "9 years", "25 years"), col = c("deepskyblue3", 'gray52', 'firebrick'), title = expression(bold("Time since \n fire")), lty = 1, cex = 2, bty = "n")
par(xpd = F)


# Plot effect sizes for time since fire and proportion recruits ----

# Create coefficient tables for plotting
recruit_m5t
recruit_m7t

# Temperature only m7
temp_coefs <- data.frame(
  Estimate = c(summary(recruit_m7t)$coefficients[,1]),
  SE = c(summary(recruit_m7t)$coefficients[,2]),
  Term = c("Intercept", "Temp")
)
temp_coefs$lci <- temp_coefs$Estimate - (temp_coefs$SE * 1.96)
temp_coefs$uci <- temp_coefs$Estimate + (temp_coefs$SE * 1.96)
rownames(temp_coefs) <- temp_coefs$Term


# Temperature and fire frequency m5
tff_coefs <- data.frame(
  Estimate = c(summary(recruit_m5t)$coefficients[,1]),
  SE = c(summary(recruit_m5t)$coefficients[,2]),
  Term = c("Intercept", "TSF", "Temp", "TSF:Temp")
)
tff_coefs$lci <- tff_coefs$Estimate - (tff_coefs$SE * 1.96)
tff_coefs$uci <- tff_coefs$Estimate + (tff_coefs$SE * 1.96)
rownames(tff_coefs) <- tff_coefs$Term



dev.new(height=5,width=10, dpi=80,pointsize=12,noRStudioGD = T)
par(mfrow=c(1,2),mar=c(5,6, 2,1),mgp=c(2.7,1,0), cex = 1, cex.axis = 1, cex.lab = 1.5)

plot(temp_coefs$Estimate, rev(1:nrow(temp_coefs)), xlim = c(min(temp_coefs$lci), max(temp_coefs$uci)), las = 1, cex = 1.8, ylab = "", xlab = expression(bold("Effect size")), pch = 20, yaxt = 'n', col = "black")
axis(side = 2, at = rev(1:nrow(temp_coefs)), labels = rownames(temp_coefs), las = 1)
arrows(temp_coefs$uci, rev(1:nrow(temp_coefs)), temp_coefs$lci, rev(1:nrow(temp_coefs)), code = 0, lwd = 0.8)
arrows(0,0,0, 7, code = 0, lwd = 0.8)
mtext('(a)', adj = 0000, line = 0.2, cex = 1.5)

plot(tff_coefs$Estimate, rev(1:nrow(tff_coefs)), xlim = c(min(tff_coefs$lci), max(tff_coefs$uci)), las = 1, cex = 1.8, ylab = "", xlab =  expression(bold("Effect size")), pch = 20, yaxt = 'n', col = "black")
axis(side = 2, at = rev(1:nrow(tff_coefs)), labels = rownames(tff_coefs), las = 1)
arrows(tff_coefs$uci, rev(1:nrow(tff_coefs)), tff_coefs$lci, rev(1:nrow(tff_coefs)), code = 0, lwd = 0.8)
arrows(0,0,0, 13, code = 0, lwd = 0.8)
mtext('(b)', adj = 0000, line = 0.2, cex = 1.5)


# This would indicate that time since fire is an uninformative variable and the interaction is not providing any new information.



# 6.4 Number of cones ----
cone_tnull <- glmer(Cone_number ~ 1 + (1 | Location/Transect), family = poisson, data = tor_cones)
cone_lnull <- glmer(Cone_number ~ 1 + (1 | Location/Transect), family = poisson, data = lit_cones)

cone_m1t <- glmer(Cone_number ~ r_TSF + (1 | Location/Transect), family = poisson, data = tor_cones)
summary(cone_m1t)
cone_m1l <- glmer(Cone_number ~ r_TSF + (1 | Location/Transect), family = poisson, data = lit_cones)
summary(cone_m1l) 

# Compare models
aictab(cand.set = c(cone_tnull, cone_m1t)) # Null model is best
aictab(cand.set = c(cone_lnull, cone_m1l)) # Null model is best

# 6.5 Relationship between number of cones and tree size ----
# Ella Plumanns-Pouton paper analyses number of cones ~ TSF + height
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2435.14619
# Similar to above but additive and interactive DBH https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2745.13480
# Also interesting but not quite the analyses we want https://www-publish-csiro-au.ap1.proxy.openathens.net/wf/Fulltext/WF14045#materials|methods

# We can produce models like those in these papers but we have so few samples for cone number and also littoralis TSF that we may be overparameterising these models.




# QUESTION 3 :  How does contemporary fire history (i.e., fire frequency) and environmental attributes influence reproductive traits? ----
# 7.1 Proportion of seedlings to adults ----
seed_fft_null <- glmer(Proportion_seedlings ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
seed_ffl_null <- glmer(Proportion_seedlings ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)


seed_ff_m1t <- glmer(Proportion_seedlings ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m1t)
seed_ff_m1l <- glmer(Proportion_seedlings ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m1l)


seed_ff_m2t <- glmer(Proportion_seedlings ~ r_Fire_freq * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m2t)
seed_ff_m2l <- glmer(Proportion_seedlings ~ Fire_freq * r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=20)))
summary(seed_ff_m2l)


seed_ff_m3t <- glmer(Proportion_seedlings ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m3t)
seed_ff_m3l <- glmer(Proportion_seedlings ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m3l)


seed_ff_m4t <- glmer(Proportion_seedlings ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m4t)
seed_ff_m4l <- glmer(Proportion_seedlings ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m4l)


seed_ff_m5t <- glmer(Proportion_seedlings ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m5t)
seed_ff_m5l <- glmer(Proportion_seedlings ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m5l)


seed_ff_m6t <- glmer(Proportion_seedlings ~ r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m6t)
seed_ff_m6l <- glmer(Proportion_seedlings ~ r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m6l)


seed_ff_m7t <- glmer(Proportion_seedlings ~ r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m7t)
seed_ff_m7l <- glmer(Proportion_seedlings ~ r_Temp + (1|Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m7l)


seed_ff_m8t <- glmer(Proportion_seedlings ~ r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m8t)
seed_ff_m8l <- glmer(Proportion_seedlings ~ r_FPC + (1|Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m8l)


seed_ff_m9t <- glmer(Proportion_seedlings ~ r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m9t)
seed_ff_m9l <- glmer(Proportion_seedlings ~ r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=140)))
summary(seed_ff_m9l)


seed_ff_m10t <- glmer(Proportion_seedlings ~ r_Fire_freq * r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m10t)
seed_ff_m10l <- glmer(Proportion_seedlings ~ r_Fire_freq * r_TWI + (1|Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m10l)


seed_ff_m11t <- glmer(Proportion_seedlings ~ r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(seed_ff_m11t)
seed_ff_m11l <- glmer(Proportion_seedlings ~ r_TWI + (1|Location/Transect), family = binomial, data = lit_transects)
summary(seed_ff_m11l)

# Compare models
aictab(cand.set = c(seed_fft_null, seed_ff_m1t, seed_ff_m2t, seed_ff_m3t, seed_ff_m4t, seed_ff_m5t, seed_ff_m6t, seed_ff_m7t, seed_ff_m8t, seed_ff_m9t, seed_ff_m10t, seed_ff_m11t), modnames = c("Null", "FF", "FF*Lat", "FF*FPC", "FF*Precip", "FF*Temp", "Precip", "Temp", "FPC", "Lat", "FF*TWI", "TWI")) 
# Model 8 is best, so FPC has a stronger influence on proportion of seedlings than fire frequency. The univariate precipitation seasonality model is also ranked higher than the univariate fire frequency model, suggesting that both FPC and precipitation seasonality have stronger influence on proportion of seedlings than fire frequency.


aictab(cand.set = c(seed_ffl_null, seed_ff_m1l, seed_ff_m2l, seed_ff_m3l, seed_ff_m4l, seed_ff_m5l, seed_ff_m6l, seed_ff_m7l, seed_ff_m8l, seed_ff_m9l, seed_ff_m10l, seed_ff_m11l), modnames = c("Null", "FF", "FF*Lat", "FF*FPC", "FF*Precip", "FF*Temp", "Precip", "Temp", "FPC", "Lat", "FF*TWI", "TWI"))
# The null model is best, but the univariate temperature seasonality and latitude models are ranked higher than the univariate fire frequency model, suggesting that temperature seasonalityand latitude have a stronger influence on proportion of seedlings than fire frequency.


# 7.2 Proportion of saplings to adults ----
sap_ff_tnull <- glmer(Proportion_saplings ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
sap_ff_lnull <- glmer(Proportion_saplings ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)


sap_ff_m1t <- glmer(Proportion_saplings ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m1t)
sap_ff_m1l <- glmer(Proportion_saplings ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30))) # Not a very good model
summary(sap_ff_m1l)


sap_ff_m2t <- glmer(Proportion_saplings ~ r_Fire_freq * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m2t)
sap_ff_m2l <- glmer(Proportion_saplings ~ r_Fire_freq * r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=250)))
summary(sap_ff_m2l) 


sap_ff_m3t <- glmer(Proportion_saplings ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m3t)
sap_ff_m3l <- glmer(Proportion_saplings ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects)
summary(sap_ff_m3l)


sap_ff_m4t <- glmer(Proportion_saplings ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m4t)
sap_ff_m4l <- glmer(Proportion_saplings ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(sap_ff_m4l)


sap_ff_m5t <- glmer(Proportion_saplings ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m5t)
sap_ff_m5l <- glmer(Proportion_saplings ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects)
summary(sap_ff_m5l)


sap_ff_m6t <- glmer(Proportion_saplings ~ r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m6t)
sap_ff_m6l <- glmer(Proportion_saplings ~ r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(sap_ff_m6l)


sap_ff_m7t <- glmer(Proportion_saplings ~ r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m7t)
sap_ff_m7l <- glmer(Proportion_saplings ~ r_Temp + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(sap_ff_m7l)


sap_ff_m8t <- glmer(Proportion_saplings ~ r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m8t)
sap_ff_m8l <- glmer(Proportion_saplings ~ r_FPC + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(sap_ff_m8l)


sap_ff_m9t <- glmer(Proportion_saplings ~ r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(sap_ff_m9t)
sap_ff_m9l <- glmer(Proportion_saplings ~ r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(sap_ff_m9l)


sap_ff_m10t <- glmer(Proportion_saplings ~ r_Fire_freq * r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m10t)
sap_ff_m10l <- glmer(Proportion_saplings ~ r_Fire_freq * r_TWI + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=50)))
summary(rsap_m10l)


sap_ff_m11t <- glmer(Proportion_saplings ~ r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(rsap_m11t)
sap_ff_m11l <- glmer(Proportion_saplings ~ r_TWI + (1|Location/Transect), family = binomial, data = lit_transects)
summary(rsap_m11l)

# Compare models
aictab(cand.set = c(sap_ff_tnull, sap_ff_m1t, sap_ff_m2t, sap_ff_m3t, sap_ff_m4t, sap_ff_m5t, sap_ff_m6t, sap_ff_m7t, sap_ff_m8t, sap_ff_m9t, sap_ff_m10t, sap_ff_m11t), modnames = c("Null", "FF", "FF*Lat", "FF*FPC", "FF*Precip", "FF*Temp", "Precip", "Temp", "FPC", "Lat", "FF*TWI", "TWI")) 
# The univariate temperature model is best. So temperature seasonality alone has a stronger influence on proportion of saplings than  fire frequency.

aictab(cand.set = c(sap_ff_lnull, sap_ff_m1l, sap_ff_m2l, sap_ff_m3l, sap_ff_m4l, sap_ff_m5l, sap_ff_m6l, sap_ff_m7l, sap_ff_m8l, sap_ff_m9l, sap_ff_m10l, sap_ff_m11l), modnames = c("Null", "FF", "FF*Lat", "FF*FPC", "FF*Precip", "FF*Temp", "Precip", "Temp", "FPC", "Lat", "FF*TWI", "TWI"))
# The null model is best, but fire frequency had a stronger influence on proportion of saplings than environmental covariates.


# 7.3 Proportion of recruits to adults ----
recff_tnull <- glmer(Proportion_recruits ~ 1 + (1 | Location/Transect), family = binomial, data = tor_transects)
recff_lnull <- glmer(Proportion_recruits ~ 1 + (1 | Location/Transect), family = binomial, data = lit_transects)


recff_m1t <- glmer(Proportion_recruits ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = tor_transects)
summary(recff_m1t)
recff_m1l <- glmer(Proportion_recruits ~ r_Fire_freq + (1 | Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30))) # Not a very good model
summary(recff_m1l) 


recff_m2t <- glmer(Proportion_recruits ~ r_Fire_freq * r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recff_m2t)
recff_m2l <- glmer(Proportion_recruits ~ Fire_freq * Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30))) 
summary(recff_m2l)


recff_m3t <- glmer(Proportion_recruits ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recff_m3t)
recff_m3l <- glmer(Proportion_recruits ~ r_Fire_freq * r_FPC + (1|Location/Transect), family = binomial, data = lit_transects)
summary(recff_m3l)


recff_m4t <- glmer(Proportion_recruits ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recff_m4t)
recff_m4l <- glmer(Proportion_recruits ~ r_Fire_freq * r_Precip + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=80)))# Not a very good model
summary(recff_m4l) 


recff_m5t <- glmer(Proportion_recruits ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recff_m5t)
recff_m5l <- glmer(Proportion_recruits ~ r_Fire_freq * r_Temp + (1|Location/Transect), family = binomial, data = lit_transects)
summary(recff_m5l) 


recff_m6t <- glmer(Proportion_recruits ~ r_Precip + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recff_m6t)
recff_m6l <- glmer(Proportion_recruits ~ r_Precip + (1|Location/Transect), family = binomial, data = lit_transects)
summary(recff_m6l) 


recff_m7t <- glmer(Proportion_recruits ~ r_Temp + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recff_m7t)
recff_m7l <- glmer(Proportion_recruits ~ r_Temp + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(recff_m7l) 


recff_m8t <- glmer(Proportion_recruits ~ r_FPC + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recff_m8t)
recff_m8l <- glmer(Proportion_recruits ~ r_FPC + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=20)))
summary(recff_m8l) 


recff_m9t <- glmer(Proportion_recruits ~ r_Latitude + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recff_m9t)
recff_m9l <- glmer(Proportion_recruits ~ r_Latitude + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=30)))
summary(recff_m9l) 



recff_m10t <- glmer(Proportion_recruits ~ r_Fire_freq * r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recff_m10t)
recff_m10l <- glmer(Proportion_recruits ~ r_Fire_freq * r_TWI + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=140)))
summary(recff_m10l)


recff_m11t <- glmer(Proportion_recruits ~ r_TWI + (1|Location/Transect), family = binomial, data = tor_transects)
summary(recff_m11t)
recff_m11l <- glmer(Proportion_recruits ~ r_TWI + (1|Location/Transect), family = binomial, data = lit_transects, control = glmerControl(optCtrl=list(maxfun=140)))
summary(recff_m11l)


# Compare models
aictab(cand.set = c(recff_tnull, recff_m1t, recff_m2t, recff_m3t, recff_m4t, recff_m5t, recff_m6t, recff_m7t, recff_m8t, recff_m9t, recff_m10t, recff_m11t), modnames = c("Null", "FF", "FF*Lat", "FF*FPC", "FF*Precip", "FF*Temp", "Precip", "Temp", "FPC", "Lat", "FF*TWI", "TWI")) 
# The univariate temperature model is the best model. So temperature seasonality along with all environmental covariates have a stronger influence on proportion of recruits than fire frequency. 

aictab(cand.set = c(recff_lnull, recff_m1l, recff_m2l, recff_m3l, recff_m4l, recff_m5l, recff_m6l, recff_m7l, recff_m8l, recff_m9l, recff_m10l, recff_m11l), modnames = c("Null", "FF", "FF*Lat", "FF*FPC", "FF*Precip", "FF*Temp", "Precip", "Temp", "FPC", "Lat", "FF*TWI", "TWI"))
# The null model is best but fire frequency has a stronger influence on proportion of recruits than any other environmental covariates. 




# Plot effect sizes for top models of population structure ----

# Create coefficient tables for plotting

# Torulosa
seedt <- data.frame(
  Estimate = c(summary(seed_ff_m8t)$coefficients[,1]),
  SE = c(summary(seed_ff_m8t)$coefficients[,2]),
  Term = c("Intercept", "FPC")
)
seedt$lci <-seedt$Estimate - (seedt$SE * 1.96)
seedt$uci <- seedt$Estimate + (seedt$SE * 1.96)
rownames(seedt) <- seedt$Term

sapt <- data.frame(
  Estimate = c(summary(sap_ff_m7t)$coefficients[,1]),
  SE = c(summary(sap_ff_m7t)$coefficients[,2]),
  Term = c("Intercept", "Temperature \n seasonality")
)
sapt$lci <-sapt$Estimate - (sapt$SE * 1.96)
sapt$uci <- sapt$Estimate + (sapt$SE * 1.96)
rownames(sapt) <- sapt$Term

rect <- data.frame(
  Estimate = c(summary(recff_m7t)$coefficients[,1]),
  SE = c(summary(recff_m7t)$coefficients[,2]),
  Term = c("Intercept", "Temperature \n seasonality")
)
rect$lci <-rect$Estimate - (rect$SE * 1.96)
rect$uci <- rect$Estimate + (rect$SE * 1.96)
rownames(rect) <- rect$Term




dev.new(height=5,width=15, dpi=80,pointsize=12,noRStudioGD = T)
par(mfrow=c(1,3),mar=c(4,6, 2,1),mgp=c(2.7,1,0), cex = 1, cex.axis = 1, cex.lab = 1.5)


plot(seedt$Estimate, rev(1:nrow(seedt)), xlim = c(-6,6), las = 1, cex = 1.2, ylab = "", xlab = expression(bold("Effect size")), pch = 20, yaxt = 'n', col = "black")
axis(side = 2, at = rev(1:nrow(seedt)), labels = rownames(seedt), las = 1)
arrows(seedt$uci, rev(1:nrow(seedt)), seedt$lci, rev(1:nrow(seedt)), code = 0, lwd = 0.8)
arrows(0,0,0, 7, code = 0, lwd = 0.8)
mtext('(a) Proportion seedlings', adj = 0000, line = 0.2, cex = 1.5)

plot(sapt$Estimate, rev(1:nrow(sapt)), xlim = c(-6,6), las = 1, cex = 1.2, ylab = "", xlab =  expression(bold("Effect size")), pch = 20, yaxt = 'n', col = "black")
axis(side = 2, at = rev(1:nrow(sapt)), labels = rownames(sapt), las = 1)
arrows(sapt$uci, rev(1:nrow(sapt)), sapt$lci, rev(1:nrow(sapt)), code = 0, lwd = 0.8)
arrows(0,0,0, 13, code = 0, lwd = 0.8)
mtext('(b) Proportion saplings', adj = 0000, line = 0.2, cex = 1.5)

plot(rect$Estimate, rev(1:nrow(rect)), xlim = c(-6,6), las = 1, cex = 1.2, ylab = "", xlab =  expression(bold("Effect size")), pch = 20, yaxt = 'n', col = "black")
axis(side = 2, at = rev(1:nrow(rect)), labels = rownames(rect), las = 1)
arrows(rect$uci, rev(1:nrow(rect)), rect$lci, rev(1:nrow(rect)), code = 0, lwd = 0.8)
arrows(0,0,0, 13, code = 0, lwd = 0.8)
mtext('(c) Proportion recruits', adj = 0000, line = 0.2, cex = 1.5)


# 7.4 Number of cones ----
cone_m2t <- glmer(Cone_number ~ r_Fire_freq + (1 | Location/Transect), family = poisson, data = tor_cones)
summary(cone_m1t)
cone_m2l <- glmer(Cone_number ~ r_Fire_freq + (1 | Location/Transect), family = poisson, data = lit_cones)
summary(cone_m1l) 

# Compare models
aictab(cand.set = c(cone_tnull, cone_m1t, cone_m2t), modnames = c("Null", "TSF", "FF"))
# Fire frequency is best
aictab(cand.set = c(cone_lnull, cone_m1l, cone_m2l), modnames = c("Null", "TSF", "FF")) 
# Fire frequency is best

# Determine whether any environmental attributes are more influential on cone number than fire frequency
cone_m3t <- glmer(Cone_number ~ r_Latitude + (1 | Location/Transect), family = poisson, data = tor_cones)
summary(cone_m3t)
cone_m3l <- glmer(Cone_number ~ r_Latitude + (1 | Location/Transect), family = poisson, data = lit_cones)
summary(cone_m3l)

cone_m4t <- glmer(Cone_number ~ r_FPC + (1 | Location/Transect), family = poisson, data = tor_cones)
summary(cone_m4t)
cone_m4l <- glmer(Cone_number ~ r_FPC + (1 | Location/Transect), family = poisson, data = lit_cones)
summary(cone_m4l)

cone_m5t <- glmer(Cone_number ~ r_Precip + (1 | Location/Transect), family = poisson, data = tor_cones)
summary(cone_m5t)
cone_m5l <- glmer(Cone_number ~ r_Precip + (1 | Location/Transect), family = poisson, data = lit_cones)
summary(cone_m5l)

cone_m6t <- glmer(Cone_number ~ r_Temp + (1 | Location/Transect), family = poisson, data = tor_cones)
summary(cone_m6t)
cone_m6l <- glmer(Cone_number ~ r_Temp + (1 | Location/Transect), family = poisson, data = lit_cones)
summary(cone_m6l)

cone_m7t <- glmer(Cone_number ~ r_TWI  + (1 | Location/Transect), family = poisson, data = tor_cones)
summary(cone_m7t)
cone_m7l <- glmer(Cone_number ~ r_TWI + (1 | Location/Transect), family = poisson, data = lit_cones)
summary(cone_m7l)




# Compare models
aictab(cand.set = c(cone_tnull, cone_m1t, cone_m2t, cone_m3t, cone_m4t, cone_m5t, cone_m6t, cone_m7t), modnames = c("Null", "TSF", "FF", "Lat", "FPC", "Precip", "Temp", "TWI")) # The FPC only model is better than fire frequency only

aictab(cand.set = c(cone_lnull, cone_m1l, cone_m2l, cone_m3l, cone_m4l, cone_m5l, cone_m6l, cone_m7l), modnames = c("Null", "TSF", "FF", "Lat", "FPC", "Precip", "Temp", "TWI"))
# The FPC only model is best and the FPC and TWI only models are better than fire frequency only. 




# Predict to new data for number of cones ~ FPC ----
cff_t <- data.frame(r_FPC = seq(min(tor_cones$r_FPC, na.rm = T), max(tor_cones$r_FPC, na.rm = T), length = 50))
cff_l <- data.frame(r_FPC = seq(min(lit_cones$r_FPC, na.rm = T), max(lit_cones$r_FPC, na.rm = T), length = 50))

pcfft <- predictSE(cone_m4t, newdata = cff_t, se.fit = T, type = 'response')
cff_t$fit <- pcfft$fit
cff_t$se <- pcfft$se.fit
cff_t$lci <- cff_t$fit - (cff_t$se * 1.96)
cff_t$uci <- cff_t$fit + (cff_t$se * 1.96)


pcffl <- predictSE(cone_m4l, newdata = cff_l, se.fit = T, type = 'response')
cff_l$fit <- pcffl$fit
cff_l$se <- pcffl$se.fit
cff_l$lci <- cff_l$fit - (cff_l$se * 1.96)
cff_l$uci <- cff_l$fit + (cff_l$se * 1.96)


# Plot predictions for number of cones ~ FPC ----
dev.new(width = 20, height = 8, dpi = 300, noRStudioGD = T)
par(mfrow = c(1,2), mar = c(4,5,2,2))

plot(cff_l$r_FPC, cff_l$fit, type = 'l', las = 1, cex.axis = 1.4, ylab = "", xlab = "", ylim = c(0, 500), xlim = c(-3.3, 3.1), xaxt = "n")
axis(side = 1, at = seq(min(lit_cones$r_FPC, na.rm = T), max(lit_cones$r_FPC, na.rm = T), length = 10), labels = round(seq(min(lit_cones$FPC, na.rm = T), max(lit_cones$FPC, na.rm = T), length = 10)), cex.axis = 1.4)
mtext(side = 2, expression(bold("Number of cones")), cex = 1.5, line = 3.5)
mtext(side = 1, expression(bold("Foliage projective cover")), cex = 1.5, line = 3)
pg.ci(x = "r_FPC", data = "cff_l", colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(a)")), cex = 2, adj = 0.001)
points(lit_cones$r_FPC, lit_cones$Cone_number, pch = 19, cex = 0.4, col = 'gray36')

plot(cff_t$r_FPC, cff_t$fit, type = 'l', las = 1, cex.axis = 1.4, ylab = "", xlab = "", ylim = c(0, 500), xlim = c(-1.5, 2.3), xaxt = "n")
axis(side = 1, at = seq(min(tor_cones$r_FPC, na.rm = T), max(tor_cones$r_FPC, na.rm = T), length = 14), labels = round(seq(min(tor_cones$FPC, na.rm = T), max(tor_cones$FPC, na.rm = T), length = 14)), cex.axis = 1.4)
mtext(side = 2, expression(bold("Number of cones")), cex = 1.5, line = 3.5)
mtext(side = 1, expression(bold("Foliage projective cover")), cex = 1.5, line = 3)
pg.ci(x = 'r_FPC', data = 'cff_t', colour = rgb(0,0,0,0.1), lower = 'lci', upper = 'uci')
mtext(expression(bold("(b)")), cex = 2, adj = 0.001)
points(tor_cones$r_FPC, tor_cones$Cone_number, pch = 19, cex = 0.4, col = 'gray36')

# 7.5 Relationship between number of cones and tree size ----
dim(tor_cones) # 35 samples 
unique(tor_cones$Fire_freq) # 5 fire frequencies
dim(lit_cones) # 45 samples
unique(lit_cones$Fire_freq) # 3 fire frequencies

# I think we have too few variables to model cone number with height. Not finding anything in the literature that has proportion of cones depending on height.



# Plot effect sizes for top model of reproductive output ----

# Create coefficient tables for plotting
cone_m4t
cone_m4l

# Littoralis
litcon_coef <- data.frame(
  Estimate = c(summary(cone_m4l)$coefficients[,1]),
  SE = c(summary(cone_m4l)$coefficients[,2]),
  Term = c("Intercept", "FPC")
)
litcon_coef$lci <-litcon_coef$Estimate - (litcon_coef$SE * 1.96)
litcon_coef$uci <- litcon_coef$Estimate + (litcon_coef$SE * 1.96)
rownames(litcon_coef) <- litcon_coef$Term


# Torulosa
torcon_coef <- data.frame(
  Estimate = c(summary(cone_m4t)$coefficients[,1]),
  SE = c(summary(cone_m4t)$coefficients[,2]),
  Term = c("Intercept", "FPC")
)
torcon_coef$lci <- torcon_coef$Estimate - (torcon_coef$SE * 1.96)
torcon_coef$uci <- torcon_coef$Estimate + (torcon_coef$SE * 1.96)
rownames(torcon_coef) <- torcon_coef$Term



dev.new(height=5,width=10, dpi=80,pointsize=12,noRStudioGD = T)
par(mfrow=c(1,2),mar=c(5,6, 2,1),mgp=c(2.7,1,0), cex = 1, cex.axis = 1, cex.lab = 1.5)

plot(litcon_coef$Estimate, rev(1:nrow(litcon_coef)), xlim = c(-2,6), las = 1, cex = 1.2, ylab = "", xlab = expression(bold("Effect size")), pch = 20, yaxt = 'n', col = "black")
axis(side = 2, at = rev(1:nrow(litcon_coef)), labels = rownames(litcon_coef), las = 1)
arrows(litcon_coef$uci, rev(1:nrow(litcon_coef)), litcon_coef$lci, rev(1:nrow(litcon_coef)), code = 0, lwd = 0.8)
arrows(0,0,0, 7, code = 0, lwd = 0.8)
mtext('(a)'~italic(Allocasuarina~littoralis), adj = 0000, line = 0.2, cex = 1.5)

plot(torcon_coef$Estimate, rev(1:nrow(torcon_coef)), xlim = c(min(torcon_coef$lci), max(torcon_coef$uci)), las = 1, cex = 1.2, ylab = "", xlab =  expression(bold("Effect size")), pch = 20, yaxt = 'n', col = "black")
axis(side = 2, at = rev(1:nrow(torcon_coef)), labels = rownames(torcon_coef), las = 1)
arrows(torcon_coef$uci, rev(1:nrow(torcon_coef)), torcon_coef$lci, rev(1:nrow(torcon_coef)), code = 0, lwd = 0.8)
arrows(0,0,0, 13, code = 0, lwd = 0.8)
mtext('(b)'~italic(Allocasuarina~torulosa), adj = 0000, line = 0.2, cex = 1.5)





# 7.6 Seed size ----
# Problem - unable to predict to new data with default Gamma as the inverse link function is not supported by AICcmodavg::predictSE. This is what the canonical link function error actually means. So we either have to use a different family or change the link function to link = log
# Solutions:
  # 1. Run MASS::glmmPQL instead but we cannot compare AIC as this glmmPQL uses a quasi-loglikelihood meaning log likelihood is not calculated.
  # 2. Fit glmer with Gamma(link = 'log'), apparently uncommon to fit Gamma with the natural or default link (i.e., inverse link)


# In reality the link = log is appropriate for this data as it constrains the value to be positive >=0.

seedwt_t_null <- glmer(seed_wt_mg ~ 1 + (1|Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
seedwt_l_null <- glmer(seed_wt_mg ~ 1 + (1|Individual), family = Gamma(link = 'log'), data = lit_cum.prop)


seedwt_m1t <- glmer(seed_wt_mg ~ r_fire_freq + (1 | Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m1t)
seedwt_m1l <- glmer(seed_wt_mg ~ r_fire_freq + (1 | Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m1l)


seedwt_m2t <- glmer(seed_wt_mg ~ r_fire_freq * r_Latitude + (1 | Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m2t)
seedwt_m2l <- glmer(seed_wt_mg ~ r_fire_freq * r_Latitude + (1 | Individual), family = Gamma(link = 'log'), data = lit_cum.prop) 
summary(seedwt_m2l)


seedwt_m3t <- glmer(seed_wt_mg ~ r_fire_freq * r_FPC + (1 | Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m3t)
seedwt_m3l <- glmer(seed_wt_mg ~ r_fire_freq * r_FPC + (1 | Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m3l)


seedwt_m4t <- glmer(seed_wt_mg ~ r_fire_freq * r_Precip + (1 | Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m4t)
seedwt_m4l <- glmer(seed_wt_mg ~ r_fire_freq * r_Precip + (1 | Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m4l)


seedwt_m5t <- glmer(seed_wt_mg ~ r_fire_freq * r_Temp + (1 | Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m5t)
seedwt_m5l <- glmer(seed_wt_mg ~ r_fire_freq * r_Temp + (1 | Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m5l)


seedwt_m6t <- glmer(seed_wt_mg ~ r_Precip + (1 | Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m6t)
seedwt_m6l <- glmer(seed_wt_mg ~ r_Precip + (1 | Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m6l)


seedwt_m7t <- glmer(seed_wt_mg ~ r_Temp + (1 | Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m7t)
seedwt_m7l <- glmer(seed_wt_mg ~ r_Temp + (1 | Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m7l)


seedwt_m8t <- glmer(seed_wt_mg ~ r_FPC + (1 | Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m8t)
seedwt_m8l <- glmer(seed_wt_mg ~ r_FPC + (1 | Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m8l)


seedwt_m9t <- glmer(seed_wt_mg ~ r_Latitude + (1 | Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m9t)
seedwt_m9l <- glmer(seed_wt_mg ~ r_Latitude + (1 | Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m9l)


seedwt_m10t <- glmer(seed_wt_mg ~ r_fire_freq * r_TWI + (1 | Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m10t)
seedwt_m10l <- glmer(seed_wt_mg ~ r_fire_freq * r_TWI + (1 | Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m10l)


seedwt_m11t <- glmer(seed_wt_mg ~ r_TWI + (1 | Individual), family = Gamma(link = 'log'), data = tor_cum.prop)
summary(seedwt_m11t)
seedwt_m11l <- glmer(seed_wt_mg ~ r_TWI + (1 | Individual), family = Gamma(link = 'log'), data = lit_cum.prop)
summary(seedwt_m11l)

# Compare models
aictab(cand.set = c(seedwt_t_null, seedwt_m1t, seedwt_m2t, seedwt_m3t, seedwt_m4t, seedwt_m5t, seedwt_m6t, seedwt_m7t, seedwt_m8t, seedwt_m9t, seedwt_m10t, seedwt_m11t), modnames = c("Null", "FF", "FF*Lat", "FF*FPC", "FF*Precip", "FF*Temp", "Precip", "Temp", "FPC", "Lat", "FF*TWI", "TWI"))
# The model with fire frequency and Latitude is the best model and the model with fire frequency and TWI is ranked within delta AICc <2. 


aictab(cand.set = c(seedwt_l_null, seedwt_m1l, seedwt_m2l, seedwt_m3l, seedwt_m4l, seedwt_m5l, seedwt_m6l, seedwt_m7l, seedwt_m8l, seedwt_m9l, seedwt_m10l, seedwt_m11l), modnames = c("Null", "FF", "FF*Lat", "FF*FPC", "FF*Precip", "FF*Temp", "Precip", "Temp", "FPC", "Lat", "FF*TWI", "TWI"))
# Model with fire frequency and temperature is the best model.


save.image('./02_Workspaces/002_data_analysis.RData')



# Predict to new data for seed size ----
# Littoralis
lseed_m5_l <- data.frame(r_fire_freq = min(lit_cum.prop$r_fire_freq),
                         r_Temp = seq(min(lit_cum.prop$r_Temp), max(lit_cum.prop$r_Temp), length = 50))
lseed_m5_a <- data.frame(r_fire_freq = mean(lit_cum.prop$r_fire_freq),
                         r_Temp = seq(min(lit_cum.prop$r_Temp), max(lit_cum.prop$r_Temp), length = 50))
lseed_m5_h <- data.frame(r_fire_freq = max(lit_cum.prop$r_fire_freq),
                         r_Temp = seq(min(lit_cum.prop$r_Temp), max(lit_cum.prop$r_Temp), length = 50))

plseed_m5_l <- predictSE(seedwt_m5l, newdata = lseed_m5_l, se.fit = T, type = 'response')
lseed_m5_l$fit <- plseed_m5_l$fit
lseed_m5_l$se <- plseed_m5_l$se.fit
lseed_m5_l$lci <- lseed_m5_l$fit - (lseed_m5_l$se * 1.96)
lseed_m5_l$uci <- lseed_m5_l$fit + (lseed_m5_l$se * 1.96)


plseed_m5_a <- predictSE(seedwt_m5l, newdata = lseed_m5_a, se.fit = T, type = 'response')
lseed_m5_a$fit <- plseed_m5_a$fit
lseed_m5_a$se <- plseed_m5_a$se.fit
lseed_m5_a$lci <- lseed_m5_a$fit - (lseed_m5_a$se * 1.96)
lseed_m5_a$uci <- lseed_m5_a$fit + (lseed_m5_a$se * 1.96)


plseed_m5_h <- predictSE(seedwt_m5l, newdata = lseed_m5_h, se.fit = T, type = 'response')
lseed_m5_h$fit <- plseed_m5_h$fit
lseed_m5_h$se <- plseed_m5_h$se.fit
lseed_m5_h$lci <- lseed_m5_h$fit - (lseed_m5_h$se * 1.96)
lseed_m5_h$uci <- lseed_m5_h$fit + (lseed_m5_h$se * 1.96)




# Torulosa
tseed_m2_l <- data.frame(r_fire_freq = min(tor_cum.prop$r_fire_freq),
                          r_Latitude = seq(min(tor_cum.prop$r_Latitude), max(tor_cum.prop$r_Latitude), length = 50))
tseed_m2_a <- data.frame(r_fire_freq = mean(tor_cum.prop$r_fire_freq),
                          r_Latitude = seq(min(tor_cum.prop$r_Latitude), max(tor_cum.prop$r_Latitude), length = 50))
tseed_m2_h <- data.frame(r_fire_freq = max(tor_cum.prop$r_fire_freq),
                          r_Latitude = seq(min(tor_cum.prop$r_Latitude), max(tor_cum.prop$r_Latitude), length = 50))

ptseed_m2_l <- predictSE(seedwt_m2t, newdata = tseed_m2_l, se.fit = T, type = 'response')
tseed_m2_l$fit <- ptseed_m2_l$fit
tseed_m2_l$se <- ptseed_m2_l$se.fit
tseed_m2_l$lci <- tseed_m2_l$fit - (tseed_m2_l$se * 1.96)
tseed_m2_l$uci <- tseed_m2_l$fit + (tseed_m2_l$se * 1.96)


ptseed_m2_a <- predictSE(seedwt_m2t, newdata = tseed_m2_a, se.fit = T, type = 'response')
tseed_m2_a$fit <- ptseed_m2_a$fit
tseed_m2_a$se <- ptseed_m2_a$se.fit
tseed_m2_a$lci <- tseed_m2_a$fit - (tseed_m2_a$se * 1.96)
tseed_m2_a$uci <- tseed_m2_a$fit + (tseed_m2_a$se * 1.96)


ptseed_m2_h <- predictSE(seedwt_m2t, newdata = tseed_m2_h, se.fit = T, type = 'response')
tseed_m2_h$fit <- ptseed_m2_h$fit
tseed_m2_h$se <- ptseed_m2_h$se.fit
tseed_m2_h$lci <- tseed_m2_h$fit - (tseed_m2_h$se * 1.96)
tseed_m2_h$uci <- tseed_m2_h$fit + (tseed_m2_h$se * 1.96)





# Plot prediction for seed size ----
dev.new(width = 16, height = 8, dpi = 300, noRStudioGD = T)
par(mfrow = c(1,2), mar = c(5,5,3,2))


plot(lit_cum.prop$r_Temp, lit_cum.prop$seed_wt_mg, pch = 19, cex = 0.5, xlab = "", ylab = "", las = 1, cex.axis = 2, xaxt = "n", col = 'gray36', ylim = c(1,8))
axis(side = 1, at = seq(min(lit_cum.prop$r_Temp), max(lit_cum.prop$r_Temp), length = 8), labels = round(seq(min(lit_cum.prop$Temp), max(lit_cum.prop$Temp), length = 8), 0), cex.axis = 2)
lines(lseed_m5_l$r_Temp, lseed_m5_l$fit, col = 'blue')
pg.ci(x = 'r_Temp', data = 'lseed_m5_l', colour = rgb(0,0,1,0.1), lower = 'lci', upper = 'uci')
lines(lseed_m5_a$r_Temp, lseed_m5_a$fit, col = 'gray36')
pg.ci(x = 'r_Temp', data = 'lseed_m5_a', colour = rgb(92/255, 92/255, 92/255, 0.1), lower = 'lci', upper = 'uci')
lines(lseed_m5_h$r_Temp, lseed_m5_h$fit, col = 'red')
pg.ci(x = 'r_Temp', data = 'lseed_m5_h', colour = rgb(1, 0, 0, 0.1), lower = 'lci', upper = 'uci')
mtext(side = 1, expression(bold("Temperature seasonality")), line = 3.5, cex = 2.5)
mtext(side = 2, expression(bold("Seed weight (mg)")), cex = 2.5, line = 2.5)
legend(x = -1.1, y = 8.2, legend = c("0 fires", "1 fires", "4 fires"), col = c("blue", 'gray36', 'red'), title = expression(bold("Fire frequency")), lty = 1, cex = 2, bty = "n")
mtext("(a) "~italic(Allocasuarina~littoralis), cex = 2.5, adj = 0.001)



plot(tor_cum.prop$r_Latitude, tor_cum.prop$seed_wt_mg, pch = 19, cex = 0.5, xlab  = "", ylab = "", las = 1, cex.axis = 2, col = 'gray36', xaxt = "n", ylim = c(1, 8))
axis(side = 1, at = seq(min(tor_cum.prop$r_Latitude), max(tor_cum.prop$r_Latitude), length.out = 9), labels = paste("-", round(seq(min(tor_cum.prop$Latitude), max(tor_cum.prop$Latitude), length.out = 9), 1), sep = ""), cex.axis = 1.7)
axis(side = 1, at = seq(min(tor_cum.prop$r_Latitude), max(tor_cum.prop$r_Latitude), length.out = 9), labels = c("", "-27.6", "", "-27.8", "", "-28", "", "28.2", ""), cex.axis = 1.7)
lines(tseed_m2_l$r_Latitude, tseed_m2_l$fit, col = 'blue')
lines(tseed_m2_a$r_Latitude, tseed_m2_a$fit, col = 'gray36')
lines(tseed_m2_h$r_Latitude, tseed_m2_h$fit, col = 'red')
pg.ci(x = 'r_Latitude', data = 'tseed_m2_l', colour = rgb(0, 0, 1, 0.1), lower = 'lci', upper = 'uci')
pg.ci(x = 'r_Latitude', data = 'tseed_m2_a', colour = rgb(92/255, 92/255, 92/255, 0.1), lower = 'lci', upper = 'uci')
pg.ci(x = 'r_Latitude', data = 'tseed_m2_h', colour = rgb(1, 0, 0, 0.1), lower = 'lci', upper = 'uci')
mtext(side = 2, expression(bold("Seed weight (mg)")), cex = 2.5, line = 2.5)
mtext(side = 1, expression(bold("Latitude")), cex = 2.5, line = 3.5)
legend(x = -1.3, y = 8.2, legend = c("0 fires", "3 fires", "6 fires"), col = c("blue", 'gray36', 'red'), title = expression(bold("Fire frequency")), lty = 1, cex = 2, bty = "n")
mtext("(b) "~italic(Allocasuarina~torulosa), cex = 2.5, adj = 0.001)
mtext("North", side = 1, line = 3.3, adj = 0.001, cex = 1.8)
mtext("South", side = 1, line = 3.3, adj = 1, cex = 1.8)

par(xpd = NA)
arrows(x0 = c(-0.56470963, 0.17167138), x1 = c(-1.05563030, 0.66259206), y0 = c(-0.2, -0.2), y1  = c(-0.2, -0.2), code = 2, length = 0.2)
par(xpd = F)

# 8. Plot raw data plot of number of plants in each size class ~ fire frequency ----
# Reformat the transects data
tran_dat <- dat_transects[, c(2:4, 13)]
tran_dat <- rbind(tran_dat, tran_dat, tran_dat)
tran_dat$Number <- c(dat_transects[,6], dat_transects[,7], dat_transects[,8])
tran_dat$Size <- "NA"
tran_dat[1:40, 6] <- "Seedling"
tran_dat[41:80, 6] <- "Sapling"
tran_dat[81:120, 6] <- "Mature"
head(tran_dat)
tail(tran_dat)
tran_dat$Size <- as.factor(tran_dat$Size)


mypallette <- c("#BDBDBD", "#525252", "#000000")

# Create new headings for the facets
to_string <- as_labeller(c('littoralis' = "Allocasuarina littoralis", 'torulosa' = "Allocasuarina torulosa"))


# Produce the plot
tran_dat_l <- tran_dat[]


p1 <- ggplot(data = tran_dat[tran_dat$Species == "littoralis",],
       aes(x = Fire_freq, y = Number))+
  geom_col(aes(fill = Size))+
  scale_x_continuous(name = "Fire frequency", breaks = seq(-1,7, 1), limits = c(-0.5,7,1))+
  scale_y_continuous(name = "Number of plants", limits = c(0,400))+
  scale_fill_manual(values = mypallette, labels = c("Seedling", "Sapling", "Mature"), name = "Age class")+
  theme_bw()+
  theme(strip.text = element_text(face = "bold.italic"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "none")

p2 <- ggplot(data = tran_dat[tran_dat$Species == "torulosa",],
       aes(x = Fire_freq, y = Number))+
  geom_col(aes(fill = Size))+
  scale_x_continuous(name = "Fire frequency", breaks = seq(from = min(tran_dat$Fire_freq), to = max(tran_dat$Fire_freq), by = 1))+
  scale_y_continuous(name = "Number of plants", limits = c(0,400))+
  scale_fill_manual(values = mypallette, labels = c("Seedling", "Sapling", "Mature"), name = "Age class")+
  theme_bw()+
  theme(strip.text = element_text(face = "bold.italic"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

dev.new(width = 12, height = 4, dpi = 300, noRStudioGD = T)
par(mfrow = c(1,2), mar = c(5,5,3,2))

plot_grid(p1, p2, labels = c("(a)", "(b)"), label_fontface = "plain", label_size = 20, hjust = 0, rel_widths = c(0.8, 1))


save.image('./02_Workspaces/002_data_analysis.RData')
