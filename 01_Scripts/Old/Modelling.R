
# Question 1: Hows does treatment affect the proportion of seed germinated? Is there any affect of environmental or seed characteristics on proportion germination? -----
# Does seed pre-treatments, fire frequency and/or seed characteristics influence germination of A. torulosa and A. littoralis?

prop_null <- glmer(Proportion_germ ~ 1 + (1 | Individual) + (1 | Set), family = binomial, data = dat_cum.prop)


m1 <- glmer(Proportion_germ ~ Treatment * Species + (1 | Individual) + (1 | Set), family = binomial, data = dat_cum.prop)
summary(m1) # No significant effects

m2 <- glmer(Proportion_germ ~ Fire_freq * Species + (1 | Individual) + (1 | Set), family = binomial, data = dat_cum.prop)
summary(m2)# No significant effects

m3 <- glmer(Proportion_germ ~ seed_weight * Species + (1 | Individual) + (1 | Set), family = binomial, data = dat_cum.prop, control = glmerControl(optCtrl=list(maxfun=250)))
summary(m3) # Significant effect of species

m4 <- glmer(Proportion_germ ~ Treatment * Species * Fire_freq + (1 | Individual) + (1 | Set), family = binomial, data = dat_cum.prop)
summary(m4) # No significant effects

m5 <- glmer(Proportion_germ ~ Treatment * Species + seed_weight + (1 | Individual) + (1 | Set), family = binomial, data = dat_cum.prop)
summary(m5) # Significant effect of seed weight


m6 <- glmer(Proportion_germ ~ Fire_freq * Species + seed_weight + (1 | Individual) + (1 | Set), family = binomial, data = dat_cum.prop)
summary(m6) # Significant effect of seed weight



cand_set <- list(prop_null, m1, m2, m3, m4, m5, m6)
aictab(cand_set)

# The best model is model 6 with the interactive effect of fire frequency and species and additive effect of seed weight, with no model ranked within AICc <2. This model improved upon the null by AICc 19.84


# The proportion of seed germination was not affected by seed pre-treatment, however, it was affected by environmental (e.g., fire frequency) and seed (e.g., seed weight) characteristics. 


# Question 1: Generate predictions ----
m6 <- glmer(Proportion_germ ~ Fire_freq * Species + seed_weight + (1 | Individual) + (1 | Set), family = binomial, data = dat_cum.prop)
summary(m6) # Significant effect of seed weight


# We treat interactions simply as effects, such that the effect of x1 on y holding x2 constant. Treating these the same as we would for any other model predictions
# Produce new data from which to generate predictions

new_FF <- data.frame(Fire_freq = seq(min(dat_cum.prop$Fire_freq), max(dat_cum.prop$Fire_freq), length.out = 50),
                     Species = as.factor(c(rep("littoralis", 50), rep("torulosa", 50))),
                     seed_weight = mean(dat_cum.prop$seed_weight))



new_seed <- data.frame(Fire_freq = mean(dat_cum.prop$Fire_freq),
                       Species = as.factor(c(rep("littoralis", 50), rep("torulosa", 50))),
                       seed_weight = seq(min(dat_cum.prop$seed_weight), max(dat_cum.prop$seed_weight), length.out = 50))


# Generate predictions from newdata and the model
FF.pred <- predictSE(mod = m6, newdata = new_FF, se.fit = T, type = 'response')
new_FF$fit <- FF.pred$fit
new_FF$se <- FF.pred$se.fit
new_FF$lci <- new_FF$fit - (new_FF$se * 1.96)
new_FF$uci <- new_FF$fit + (new_FF$se * 1.96)
head(new_FF)



seed.pred <- predictSE(mod = m6, newdata = new_seed, se.fit = T, type = 'response')
new_seed$fit <- seed.pred$fit
new_seed$se <- seed.pred$se.fit
new_seed$lci <- new_seed$fit - (new_seed$se * 1.96)
new_seed$uci <- new_seed$fit + (new_seed$se * 1.96)
head(new_seed)


# Plot predictions ----

# We want to see how the interaction between fire frequency and species influenced proportion germination.
# Seed weight having an additive effect can be plotted separately
dev.new(width=22, height=8, res = 300, dpi=80, pointsize=18, noRStudioGD = T)
par(mfrow = c(1,2), mar = c(4,5,2,0), oma = c(0,0,0,17))

plot(new_FF$Fire_freq[new_FF$Species == "torulosa"], new_FF$fit[new_FF$Species == "torulosa"], type = 'l', ylim = c(min(new_FF$lci - 0.1), max(new_FF$uci)), col = 'red', las = 1, ylab = expression(bold("Proportion germination")), xlab = expression(bold("Fire frequency")), yaxt = "n", cex.lab = 1.5, cex.axis = 1.5)
points(new_FF$Fire_freq[new_FF$Species == "littoralis"], new_FF$fit[new_FF$Species == "littoralis"], type = 'l')
axis(side = 2, at = seq(0.2, 1, 0.2), las = 1, cex.axis = 1.5)
axis(side = 2, at = seq(0.2, 1, 0.1), labels = F, las = 1, cex.axis = 1.5)

lines(new_FF$Fire_freq[new_FF$Species == "torulosa"], new_FF$lci[new_FF$Species == "torulosa"], lty = 2, col = 'red')
lines(new_FF$Fire_freq[new_FF$Species == "torulosa"], new_FF$uci[new_FF$Species == "torulosa"], lty = 2, col = 'red')
lines(new_FF$Fire_freq[new_FF$Species == "littoralis"], new_FF$lci[new_FF$Species == "littoralis"], lty = 2)
lines(new_FF$Fire_freq[new_FF$Species == "littoralis"], new_FF$uci[new_FF$Species == "littoralis"], lty = 2)

par(xpd = NA)
legend(x = 10.9, y = 1, legend = c(expression(italic("Allocasuarina littoralis")), expression(italic("Allocasuarina torulosa"))), col = c("black", "red"), lty = 1, lwd = 2, bty = "n", cex = 1.5)
par(xpd = F)


plot(new_seed$seed_weight[new_seed$Species == "littoralis"], new_seed$fit[new_seed$Species == "littoralis"], type = 'l', ylim = c(min(new_FF$lci - 0.1), max(new_FF$uci)), las = 1, ylab = "", xlab = expression(bold("Seed weight (g)")), yaxt = "n", cex.lab = 1.5, cex.axis = 1.5)
points(new_seed$seed_weight[new_seed$Species == "torulosa"], new_seed$fit[new_seed$Species == "torulosa"], type = 'l', col = 'red')
axis(side = 2, at = seq(0.2, 1, 0.2), las = 1, cex.axis = 1.5)
axis(side = 2, at = seq(0.2, 1, 0.1), labels = F, las = 1, cex.axis = 1.5)
axis(side = 1, at = seq(0.01, 0.02, 0.01), las = 1, cex.axis = 1.5)

lines(new_seed$seed_weight[new_seed$Species == "littoralis"], new_seed$lci[new_seed$Species == "littoralis"], lty = 2)
lines(new_seed$seed_weight[new_seed$Species == "littoralis"], new_seed$uci[new_seed$Species == "littoralis"], lty = 2)
lines(new_seed$seed_weight[new_seed$Species == "torulosa"], new_seed$lci[new_seed$Species == "torulosa"], lty = 2, col = 'red')
lines(new_seed$seed_weight[new_seed$Species == "torulosa"], new_seed$uci[new_seed$Species == "torulosa"], lty = 2, col = 'red')




# Question 2: Does post-fire reproductive mode influence tolerances? Is this influenced by fire frequency. -----
# Does post-fire reproductive mode and fire frequency influence species tolerances to fire-related germination cues? 

# We ran models as part of Question 1 for selecting the most appropriate model by dropping model terms. As a reminder these are the models 

m1 <- glmer(Proportion_germ ~ Treatment * Species + (1 | Individual) + (1 | Set), family = binomial, data = dat_cum.prop)
summary(m1) # No significant effects



m4 <- glmer(Proportion_germ ~ Treatment * Species * Fire_freq + (1 | Individual) + (1 | Set), family = binomial, data = dat_cum.prop)
summary(m4) # No significant effects




# Compare only the relevant models for this question
cand_set2 <- list(prop_null, m1, m4)
aictab(cand_set2)
# The null model was the best model with no model ranked within AICc <2. 

# Post-fire reproductive mode did not have a significant influence on seed tolerances and was unaffected by the fire frequency which individuals were exposed. 






# Question 3: Do fire-related germination cues promote the speed seed germination? -----
# Consider time to 50% germination. Treatment must be included as a model term here as we are investigating the influence of fire-related germination cues.

t50_null <- glmer(t50 ~ 1 + (1 | Individual) + (1 | Set), family = poisson, data = dat_cum.prop)

m8 <- glmer(t50 ~ Treatment * Species + (1 | Individual) + (1 | Set), family = poisson, data = dat_cum.prop)
summary(m8) # Significant effects of some treatments and their interactions with species


cand_set4 <- list(t50_null, m8)
aictab(cand_set4) # The model with treatment was ranked higher than the null model, with a change in AICc of 120.92. So time to 50% germination was significantly influence by the treatment seeds had been exposed between the species. 

# Question 3: Generate predictions ----
m8 <- glmer(t50 ~ Treatment * Species + (1 | Individual) + (1 | Set), family = poisson, data = dat_cum.prop)
summary(m8) # Significant effects of some treatments and their interactions with species



new_treat_l <- data.frame(Treatment = as.factor(c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke")),
                          Species = as.factor("littoralis"))
new_treat_t <- data.frame(Treatment = as.factor(c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke")),
                          Species = as.factor("torulosa"))

new_treat <- rbind(new_treat_l, new_treat_t)

pred_treat <- predictSE(mod = m8, newdata = new_treat, se.fit = T, type = 'response')
new_treat$fit <- pred_treat$fit
new_treat$se <- pred_treat$se.fit
new_treat$lci <- new_treat$fit - (new_treat$se * 1.96)
new_treat$uci <- new_treat$fit + (new_treat$se * 1.96)
str(new_treat$Treatment)
new_treat$Treatment <- factor(new_treat$Treatment, levels = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))


# Question 3: Plot predictions -----

# Treatment x species

dev.new(width=12, height=8, res = 300, dpi=80, pointsize=18, noRStudioGD = T)
par(mfrow = c(1,1), mar = c(7,5,2,2), oma = c(0,0,0,10.3))


plot.default(new_treat$Treatment[new_treat$Species == "torulosa"], new_treat$fit[new_treat$Species == "torulosa"], type = 'p', pch = 19, ylim = c(7,15), ylab = expression(bold("Time to 50% germination")), xlab = "", xaxt = "n", col = 'red', cex.lab = 1.5, yaxt = "n")
axis(side = 1, at = c(1:6), cex.axis = 1.5, labels = F)
text(x = c(1:6), y = par("usr")[3] - 1.6, srt = 45, labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), xpd = T, cex = 1.5)
mtext(expression(bold("Treatment")), side = 1, cex = 1.5, line = 6)
arrows(c(1:6), new_treat$lci[new_treat$Species == "torulosa"], c(1:6), new_treat$uci[new_treat$Species == "torulosa"], length = 0.05, code = 3, angle = 90, col = "red")
axis(side = 2, at = seq(7,15, 1), las = 1, cex.axis = 1.5)

points(new_treat$Treatment[new_treat$Species == "littoralis"], new_treat$fit[new_treat$Species == "littoralis"], pch = 19)
arrows(c(1:6), new_treat$lci[new_treat$Species == "littoralis"], c(1:6), new_treat$uci[new_treat$Species == "littoralis"], length = 0.05, code = 3, angle = 90)


par(xpd = NA)
legend(x = 6.3, y = 15, legend = c(expression(italic("Allocasuarina littoralis")), expression(italic("Allocasuarina torulosa"))), col = c("black", "red"), lty = 1, lwd = 2, bty = "n")
par(xpd = F)
