# Written by Felicity Charles
# Caveat emptor
# Date: 14th October 2024

# Full experiment preliminary analysis


#1. Load packages
# Having a problem with this germinationmetrics::germination.indices(). A warning message was updated but does not seem to be showing the update message when I run the function. Cannot install from github.

library(germinationmetrics)
library(dplyr)




# 2. Read in data and combine by species grouping
lit1 <- read.csv('./00_Data/Full_experiment/Set1/Littoralis1.csv', header = T, stringsAsFactors = T, )
lit2 <- read.csv('./00_Data/Full_experiment/Set2/littoralis2.csv', header = T, stringsAsFactors = T)
lit3 <- read.csv('./00_Data/Full_experiment/Set3/littoralis3.csv', header = T, stringsAsFactors = T)
littoralis <- rbind(lit1, lit2, lit3)
head(littoralis); dim(littoralis)
str(littoralis)
littoralis$Treatment <- factor(littoralis$Treatment, levels = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
str(littoralis)
head(littoralis); dim(littoralis)
unique(littoralis$Treatment)


torlow1 <- read.csv('./00_Data/Full_experiment/Set1/torlow1.csv', header = T, stringsAsFactors = T)
torlow2 <- read.csv('./00_Data/Full_experiment/Set2/torlow2.csv', header = T, stringsAsFactors = T)
torlow3 <- read.csv('./00_Data/Full_experiment/Set3/torlow3.csv', header = T, stringsAsFactors = T)
torlow <- rbind(torlow1, torlow2, torlow3)
head(torlow); dim(torlow)
str(torlow)
torlow$Treatment <- factor(torlow$Treatment, levels = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
unique(torlow$Treatment)


torhigh1 <- read.csv('./00_Data/Full_experiment/Set1/torhigh1.csv', header = T, stringsAsFactors = T)
torhigh2 <- read.csv('./00_Data/Full_experiment/Set2/torhigh2.csv', header = T, stringsAsFactors = T)
torhigh3 <- read.csv('./00_Data/Full_experiment/Set3/torhigh3.csv', header = T, stringsAsFactors = T)
torhigh <- rbind(torhigh1, torhigh2, torhigh3)
head(torhigh); dim(torhigh)
str(torhigh)
torhigh$Treatment <- factor(torhigh$Treatment, levels = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
unique(torhigh$Treatment)



# 3. Calculate the cumulative proportions of germination 

# We want to calculate the adjusted proportion germinated as this will account for the fact that some seeds were not vaible, it does not appear from the preliminary investigations of the correlations between the two viability measurements that any one treatment had a significant effect on post-germination seed viability so there should not be any issue that those that were considered unviable were exposed to a treatment that necessarily killed the seeds 
lit_cum.sum <- data.frame(t(apply(littoralis[, 12:ncol(littoralis)], 1, FUN = function(x) cumsum(unlist(x)))))
littoralis$Total_adjusted <- littoralis$Total_germination + littoralis$TTC_num_viable # Calculate the number of seeds that were actually viable as we want to calculate the adjusted proportion germination for the number of viable seeds

lit_cum.prop <- lit_cum.sum/littoralis$Total_seeds
lit_cum.prop <- cbind(littoralis[, 1:11], lit_cum.prop)
head(lit_cum.prop); dim(lit_cum.prop)


lit_cum.prop_adjusted <- lit_cum.sum/littoralis$Total_adjusted
lit_cum.prop_adjusted <- cbind(littoralis[, 1:11], lit_cum.prop_adjusted)
head(lit_cum.prop_adjusted)


torl_cum.sum <- data.frame(t(apply(torlow[, 12:ncol(torlow)], 1, FUN = function(x) cumsum(unlist(x)))))
torl_cum.prop <- torl_cum.sum/torlow$Total_seeds
torl_cum.prop <- cbind(torlow[, 1:11], torl_cum.prop)
head(torl_cum.prop); dim(torl_cum.prop)



torh_cum.sum <- data.frame(t(apply(torhigh[, 12:ncol(torhigh)], 1, FUN = function(x) cumsum(unlist(x)))))
torh_cum.prop <- torh_cum.sum/torhigh$Total_seeds
torh_cum.prop <- cbind(torhigh[, 1:11], torh_cum.prop)
head(torh_cum.prop); dim(torh_cum.prop)



# 4. Boxplots of proportion germinated by treatment

### Littoralis
dev.new(height = 18, width = 40, dpi = 80, pointsize = 18, noRStudioGD = T)
par(mar = c(6,4,2,1), mgp = c(3,1,0), mfrow = c(1,3), oma = c(1, 0.5, 3, 0))

boxplot(TTC_proportion_viable ~ Treatment, data = littoralis[which(littoralis$Set == 1),], xlab = "", ylab = "Proportion germinated", las = 2)
title('Rep. 1', font.main = 1)
title(xlab = "Treatment", line = 5)


boxplot(TTC_proportion_viable ~ Treatment, data = littoralis[which(littoralis$Set == 2),], xlab = "", ylab = "Proportion germinated", las = 2)
title('Rep. 2', font.main = 1)
title(xlab = "Treatment", line = 5)
mtext(expression(bolditalic(Allocasuarina~littoralis)), side = 3, line = 3)

boxplot(TTC_proportion_viable ~ Treatment, data = littoralis[which(littoralis$Set == 3),], xlab = "", ylab = "Proportion germinated", las = 2)
title('Rep. 3', font.main = 1)
title(xlab = "Treatment", line = 5)



boxplot(TTC_proportion_viable ~ Treatment, data = littoralis, xlab = "Treatment", ylab = "Proportion germinated")



### Torulosa low fire

dev.new(height = 18, width = 40, dpi = 80, pointsize = 18, noRStudioGD = T)
par(mar = c(6,4,2,1), mgp = c(3,1,0), mfrow = c(1,3), oma = c(1, 0.5, 3, 0))


boxplot(TTC_proportion_viable ~ Treatment, data = torlow[which(torlow$Set == 1),], xlab = "", ylab = "Proportion germinated", las = 2)
title('Rep. 1', font.main = 1)
title(xlab = "Treatment", line = 5)


boxplot(TTC_proportion_viable ~ Treatment, data = torlow[which(torlow$Set == 2),], xlab = "", ylab = "Proportion germinated", las = 2)
title('Rep. 2', font.main = 1)
title(xlab = "Treatment", line = 5)
mtext(expression(bold(bolditalic(Allocasuarina~torulosa)*' low fire')), side = 3, line = 3)


boxplot(TTC_proportion_viable ~ Treatment, data = torlow[which(torlow$Set == 3),], xlab = "", ylab = "Proportion germinated", las = 2)
title('Rep. 3', font.main = 1)
title(xlab = "Treatment", line = 5)



boxplot(TTC_proportion_viable ~ Treatment, data = torlow, xlab = "Treatment", ylab = "Proportion germinated")


### Torulosa high

dev.new(height = 18, width = 40, dpi = 80, pointsize = 18, noRStudioGD = T)
par(mar = c(6,4,2,1), mgp = c(3,1,0), mfrow = c(1,3), oma = c(1, 0.5, 3, 0))


boxplot(TTC_proportion_viable ~ Treatment, data = torhigh[which(torhigh$Set == 1),], xlab = "", ylab = "Proportion germinated", las = 2)
title('Rep. 1', font.main = 1)
title(xlab = "Treatment", line = 5)


boxplot(TTC_proportion_viable ~ Treatment, data = torhigh[which(torhigh$Set == 2),], xlab = "", ylab = "Proportion germinated", las = 2)
title('Rep. 2', font.main = 1)
title(xlab = "Treatment", line = 5)
mtext(expression(bold(bolditalic(Allocasuarina~torulosa)*' high fire')), side = 3, line = 3)


boxplot(TTC_proportion_viable ~ Treatment, data = torhigh[which(torhigh$Set == 3),], xlab = "", ylab = "Proportion germinated", las = 2)
title('Rep. 3', font.main = 1)
title(xlab = "Treatment", line = 5)



boxplot(TTC_proportion_viable ~ Treatment, data = torhigh, xlab = "Treatment", ylab = "Proportion germinated")








dev.new(height = 110, width = 90, dpi = 80, pointsize = 18, noRStudioGD = T)
par(mar = c(8,4,2,1), mgp = c(3,1,0), mfrow = c(3,3), oma = c(1, 2, 3, 0))

boxplot(TTC_proportion_viable ~ Treatment, data = littoralis[which(littoralis$Set == 1),], xlab = "", ylab = expression(bold("Proportion germinated")), xaxt = "n", ylim = c(0,1), cex.lab = 1.1, las = 2)
title('Rep. 1', font.main = 1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


boxplot(TTC_proportion_viable ~ Treatment, data = littoralis[which(littoralis$Set == 2),], xlab = "", ylab = "", xaxt = "n", ylim = c(0,1), cex.lab = 1.1, las = 2)
title('Rep. 2', font.main = 1)
mtext(expression(bolditalic(Allocasuarina~littoralis)), side = 3, line = 2)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


boxplot(TTC_proportion_viable ~ Treatment, data = littoralis[which(littoralis$Set == 3),], xlab = "", ylab = "", xaxt = "n", ylim = c(0,1), cex.lab = 1.1, las = 2)
title('Rep. 3', font.main = 1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)



### Torulosa low fire
boxplot(TTC_proportion_viable ~ Treatment, data = torlow[which(torlow$Set == 1),], xlab = "", ylab = expression(bold("Proportion germinated")), xaxt = "n", ylim = c(0,1), cex.lab = 1.1, las = 2)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


boxplot(TTC_proportion_viable ~ Treatment, data = torlow[which(torlow$Set == 2),], xlab = "", ylab = "", xaxt = "n", ylim = c(0,1), cex.lab = 1.1, las = 2)
mtext(expression(bold(bolditalic(Allocasuarina~torulosa)*' low fire')), side = 3, line = 0.75)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


boxplot(TTC_proportion_viable ~ Treatment, data = torlow[which(torlow$Set == 3),], xlab = "", ylab = "", xaxt = "n", ylim = c(0,1), cex.lab = 1.1, las = 2)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


### Torulosa high
boxplot(TTC_proportion_viable ~ Treatment, data = torhigh[which(torhigh$Set == 1),], xlab = "", ylab = expression(bold("Proportion germinated")), xaxt = "n", ylim = c(0,1), cex.lab = 1.1, las = 2)
title(xlab = expression(bold("Treatment")), line = 6, cex.lab = 1.1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


boxplot(TTC_proportion_viable ~ Treatment, data = torhigh[which(torhigh$Set == 2),], xlab = "", ylab = "", xaxt = "n", ylim = c(0,1), cex.lab = 1.1, las = 2)
title(xlab = expression(bold("Treatment")), line = 6, cex.lab = 1.1)
mtext(expression(bold(bolditalic(Allocasuarina~torulosa)*' high fire')), side = 3, line = 0.75)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


boxplot(TTC_proportion_viable ~ Treatment, data = torhigh[which(torhigh$Set == 3),], xlab = "", ylab = "", xaxt = "n", ylim = c(0,1), cex.lab = 1.1, las = 2)
title(xlab = expression(bold("Treatment")), line = 6, cex.lab = 1.1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)







# Germination metrics ----

# Clean the datasets to remove unneccessary columns
littoralis_cumulative <- littoralis[, c(1, 2, 4,8)]
colnames(littoralis_cumulative) <- c("Treatment", "Rep", "Individual", "Total_seeds")
littoralis_cumulative <- cbind(littoralis_cumulative, lit_cum.sum)
littoralis_cumulative[5:65] <- lapply(littoralis_cumulative[5:65], as.numeric)
head(littoralis_cumulative); dim(littoralis_cumulative) # We can add back in the information on treatment and any thing else we need after we have the germination metrics calculated.

torlow_cumulative <- torlow[, c(1,2,4,8)]
colnames(torlow_cumulative) <- c("Treatment", "Rep", "Individual", "Total_seeds")
torlow_cumulative <- cbind(torlow_cumulative, torl_cum.sum)
torlow_cumulative[5:65] <- lapply(torlow_cumulative[5:65], as.numeric)
head(torlow_cumulative); dim(torlow_cumulative)



torhigh_cumulative <- torhigh[, c(1,2,4,8)]
colnames(torhigh_cumulative) <- c("Treatment", "Rep", "Individual", "Total_seeds")
torhigh_cumulative <- cbind(torhigh_cumulative, torh_cum.sum)
torhigh_cumulative[5:65] <- lapply(torhigh_cumulative[5:65], as.numeric)
head(torhigh_cumulative); dim(torhigh_cumulative)


# We cannot get the germinationmetrics::germination.indeices functionality to work as it either fails because of some incorrect error thrown when t50 is included, or it decides that the resulting vector would be too long, despite working for the a control only group. Will have to do it for each row separate then rbind them back together, tried creating a for loop but it fails. 

int <- 1:61

# This works but do not want to repeat this process 270 times.
l1 <- littoralis_cumulative[1,] %>%
  mutate(time_50 = t50(germ.counts = as.numeric(littoralis_cumulative[1,5:65]),intervals = int, partial = F),
         peak_germ_percent = PeakGermPercent(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), total.seeds = littoralis_cumulative[1, 4], intervals = int, partial = F),
         timespread = TimeSpreadGerm(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F),
         mean_germ_time = MeanGermTime(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F),
         var_germ_time = VarGermTime(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F),
         se_germ_time = SEGermTime(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F),
         cv_germ_time = CVGermTime(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F),
         mean_germ_rate = MeanGermRate(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F),
         coefficient_velocity_germ = CVG(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F),
         variance_germ_rate = VarGermRate(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F),
         se_germ_rate = SEGermRate(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F),
         germ_speed = GermSpeed(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F, percent = F, total.seeds = littoralis_cumulative[1, 4]),
         weighted_germ_percent = WeightGermPercent(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative[1, 4]),
         mean_germ_percent = MeanGermPercent(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative[1, 4]),
         coefficient_uniformity_germ = CUGerm(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F),
         germ_synchrony = GermSynchrony(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F),
         germ_uncertainty = GermUncertainty(germ.counts = as.numeric(littoralis_cumulative[1, 5:65]), intervals = int, partial = F))


# Throws an error message but ignored this as the expected results are returned. However, we need to split the data into blocks of 50 to make sure this is working correctly
dim(littoralis_cumulative)
littoralis_cumulative1 <- littoralis_cumulative[1:50,]
dim(littoralis_cumulative1)

littoralis_cumulative2 <- littoralis_cumulative[51:100,]
littoralis_cumulative3 <- littoralis_cumulative[101:150,]
littoralis_cumulative4 <- littoralis_cumulative[151:200,]
littoralis_cumulative5 <- littoralis_cumulative[201:250,]
littoralis_cumulative6 <- littoralis_cumulative[251:nrow(littoralis_cumulative),]


# Littoralis block 1
for(i in 1:nrow(littoralis_cumulative1)){
  
  littoralis_cumulative1$t50[i] <- t50(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative1[i, 4])
  littoralis_cumulative1$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$germ_end[i] <- LastGermTime(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F, percent = F, total.seeds = littoralis_cumulative1[i, 4])
  littoralis_cumulative1$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative1[i, 4])
  littoralis_cumulative1$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative1[i, 4])
  littoralis_cumulative1$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative1$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(littoralis_cumulative1[i, 5:65]), intervals = int, partial = F)
  
  
}

littoralis_cumulative1

# Littoralis block 2 
for(i in 1:nrow(littoralis_cumulative2)){
  
  littoralis_cumulative2$t50[i] <- t50(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative2[i, 4])
  littoralis_cumulative2$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$germ_end[i] <- LastGermTime(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F, percent = F, total.seeds = littoralis_cumulative2[i, 4])
  littoralis_cumulative2$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative2[i, 4])
  littoralis_cumulative2$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative2[i, 4])
  littoralis_cumulative2$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative2$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(littoralis_cumulative2[i, 5:65]), intervals = int, partial = F)
  
  
}

littoralis_cumulative2


# Littoralis block 3
for(i in 1:nrow(littoralis_cumulative3)){
  
  littoralis_cumulative3$t50[i] <- t50(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative3[i, 4])
  littoralis_cumulative3$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$germ_end[i] <- LastGermTime(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F, percent = F, total.seeds = littoralis_cumulative3[i, 4])
  littoralis_cumulative3$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative3[i, 4])
  littoralis_cumulative3$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative3[i, 4])
  littoralis_cumulative3$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative3$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(littoralis_cumulative3[i, 5:65]), intervals = int, partial = F)
  
  
}

littoralis_cumulative3

# Littoralis block 4
for(i in 1:nrow(littoralis_cumulative4)){
  
  littoralis_cumulative4$t50[i] <- t50(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative4[i, 4])
  littoralis_cumulative4$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$germ_end[i] <- LastGermTime(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F, percent = F, total.seeds = littoralis_cumulative4[i, 4])
  littoralis_cumulative4$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative4[i, 4])
  littoralis_cumulative4$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative4[i, 4])
  littoralis_cumulative4$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative4$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(littoralis_cumulative4[i, 5:65]), intervals = int, partial = F)
  
  
}

littoralis_cumulative4

# Littoralis block 5
for(i in 1:nrow(littoralis_cumulative5)){
  
  littoralis_cumulative5$t50[i] <- t50(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative5[i, 4])
  littoralis_cumulative5$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$germ_end[i] <- LastGermTime(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F, percent = F, total.seeds = littoralis_cumulative5[i, 4])
  littoralis_cumulative5$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative5[i, 4])
  littoralis_cumulative5$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative5[i, 4])
  littoralis_cumulative5$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative5$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(littoralis_cumulative5[i, 5:65]), intervals = int, partial = F)
  
  
}

littoralis_cumulative5

# Littoralis block 6
for(i in 1:nrow(littoralis_cumulative6)){
  
  littoralis_cumulative6$t50[i] <- t50(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative6[i, 4])
  littoralis_cumulative6$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$germ_end[i] <- LastGermTime(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F, percent = F, total.seeds = littoralis_cumulative6[i, 4])
  littoralis_cumulative6$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative6[i, 4])
  littoralis_cumulative6$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative6[i, 4])
  littoralis_cumulative6$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative6$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(littoralis_cumulative6[i, 5:65]), intervals = int, partial = F)
  
  
}

littoralis_cumulative6


# Combine littoralis blocks back into one set 
littoralis_cumulative <- rbind(littoralis_cumulative1, littoralis_cumulative2, littoralis_cumulative3, littoralis_cumulative4, littoralis_cumulative5, littoralis_cumulative6)
head(littoralis_cumulative); tail(littoralis_cumulative); dim(littoralis_cumulative)

# Check if rows match between the for loop version or a row by row calculation ----
littoralis_cumulative[1,]
l1


l15 <- littoralis_cumulative[15,] %>%
  mutate(t50 = t50(germ.counts = as.numeric(littoralis_cumulative[15,5:65]),intervals = int, partial = F),
         peak_germ_percent = PeakGermPercent(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), total.seeds = littoralis_cumulative[15, 4], intervals = int, partial = F),
         timespread = TimeSpreadGerm(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F),
         mean_germ_time = MeanGermTime(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F),
         var_germ_time = VarGermTime(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F),
         se_germ_time = SEGermTime(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F),
         cv_germ_time = CVGermTime(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F),
         mean_germ_rate = MeanGermRate(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F),
         coefficient_velocity_germ = CVG(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F),
         variance_germ_rate = VarGermRate(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F),
         se_germ_rate = SEGermRate(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F),
         germ_speed = GermSpeed(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F, percent = F, total.seeds = littoralis_cumulative[15, 4]),
         weighted_germ_percent = WeightGermPercent(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative[15, 4]),
         mean_germ_percent = MeanGermPercent(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative[15, 4]),
         coefficient_uniformity_germ = CUGerm(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F),
         germ_synchrony = GermSynchrony(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F),
         germ_uncertainty = GermUncertainty(germ.counts = as.numeric(littoralis_cumulative[15, 5:65]), intervals = int, partial = F))

l15
littoralis_cumulative[15,]



l72 <- littoralis_cumulative[72,] %>%
  mutate(t50 = t50(germ.counts = as.numeric(littoralis_cumulative[72,5:65]),intervals = int, partial = F),
         peak_germ_percent = PeakGermPercent(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), total.seeds = littoralis_cumulative[72, 4], intervals = int, partial = F),
         timespread = TimeSpreadGerm(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F),
         mean_germ_time = MeanGermTime(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F),
         var_germ_time = VarGermTime(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F),
         se_germ_time = SEGermTime(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F),
         cv_germ_time = CVGermTime(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F),
         mean_germ_rate = MeanGermRate(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F),
         coefficient_velocity_germ = CVG(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F),
         variance_germ_rate = VarGermRate(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F),
         se_germ_rate = SEGermRate(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F),
         germ_speed = GermSpeed(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F, percent = F, total.seeds = littoralis_cumulative[72, 4]),
         weighted_germ_percent = WeightGermPercent(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative[72, 4]),
         mean_germ_percent = MeanGermPercent(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative[72, 4]),
         coefficient_uniformity_germ = CUGerm(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F),
         germ_synchrony = GermSynchrony(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F),
         germ_uncertainty = GermUncertainty(germ.counts = as.numeric(littoralis_cumulative[72, 5:65]), intervals = int, partial = F))

l72
littoralis_cumulative[72,] # Now that we have done the calculation in blocks of 50, rather than as a whole large set for the for loop, the values match between the dplyr implementation and the for loop implementation




# Torulosa low fire ----

# Split data into blocks

torlow_cumulative1 <- torlow_cumulative[1:50,]
torlow_cumulative2 <- torlow_cumulative[51:100,]
torlow_cumulative3 <- torlow_cumulative[101:150,]
torlow_cumulative4 <- torlow_cumulative[151:200,]
torlow_cumulative5 <- torlow_cumulative[201:250,]
torlow_cumulative6 <- torlow_cumulative[251:nrow(torlow_cumulative),]

# Torulosa low block 1
for(i in 1:nrow(torlow_cumulative1)){
  torlow_cumulative1$t50[i] <- t50(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                  intervals = int,
                                  partial = F)
  torlow_cumulative1$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative1[i, 4])
  torlow_cumulative1$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torlow_cumulative1$germ_end[i] <- LastGermTime(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                intervals = int,
                                                partial = F)
  torlow_cumulative1$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative1$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative1$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative1$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative1$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative1$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative1$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torlow_cumulative1$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torlow_cumulative1$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                               intervals = int,
                                               partial = F)
  torlow_cumulative1$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                            intervals = int,
                                            partial = F,
                                            total.seeds = torlow_cumulative1[i, 4])
  torlow_cumulative1$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                                  intervals = int,
                                                                  partial = F,
                                                                  total.seeds = torlow_cumulative1[i, 4])
  torlow_cumulative1$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative1[i, 4])
  torlow_cumulative1$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                          intervals = int,
                                                          partial = F)
  torlow_cumulative1$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torlow_cumulative1$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torlow_cumulative1[i, 5:65]),
                                                           intervals = int,
                                                           partial = F)
  
}
torlow_cumulative1


# Torulosa low block 2
for(i in 1:nrow(torlow_cumulative2)){
  torlow_cumulative2$t50[i] <- t50(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                  intervals = int,
                                  partial = F)
  torlow_cumulative2$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative2[i, 4])
  torlow_cumulative2$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torlow_cumulative2$germ_end[i] <- LastGermTime(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                intervals = int,
                                                partial = F)
  torlow_cumulative2$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative2$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative2$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative2$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative2$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative2$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative2$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torlow_cumulative2$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torlow_cumulative2$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative2$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                               intervals = int,
                                               partial = F,
                                               total.seeds = torlow_cumulative2[i, 4])
  torlow_cumulative2$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                                  intervals = int,
                                                                  partial = F,
                                                                  total.seeds = torlow_cumulative2[i, 4])
  torlow_cumulative2$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative2[i, 4])
  torlow_cumulative2$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                             intervals = int,
                                                             partial = F)
  torlow_cumulative2$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torlow_cumulative2$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torlow_cumulative2[i, 5:65]),
                                                           intervals = int,
                                                           partial = F)
  
}
torlow_cumulative2


# Torulosa low block 3
for(i in 1:nrow(torlow_cumulative3)){
  torlow_cumulative3$t50[i] <- t50(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                  intervals = int,
                                  partial = F)
  torlow_cumulative3$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative3[i, 4])
  torlow_cumulative3$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torlow_cumulative3$germ_end[i] <- LastGermTime(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                intervals = int,
                                                partial = F)
  torlow_cumulative3$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative3$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative3$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative3$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative3$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative3$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative3$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torlow_cumulative3$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torlow_cumulative3$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative3$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                               intervals = int,
                                               partial = F,
                                               total.seeds = torlow_cumulative3[i, 4])
  torlow_cumulative3$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                                  intervals = int,
                                                                  partial = F,
                                                                  total.seeds = torlow_cumulative3[i, 4])
  torlow_cumulative3$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative3[i, 4])
  torlow_cumulative3$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                             intervals = int,
                                                             partial = F)
  torlow_cumulative3$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torlow_cumulative3$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torlow_cumulative3[i, 5:65]),
                                                           intervals = int,
                                                           partial = F)
  
}
torlow_cumulative3

# Torulosa low block 4
for(i in 1:nrow(torlow_cumulative4)){
  torlow_cumulative4$t50[i] <- t50(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                  intervals = int,
                                  partial = F)
  torlow_cumulative4$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative4[i, 4])
  torlow_cumulative4$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torlow_cumulative4$germ_end[i] <- LastGermTime(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                intervals = int,
                                                partial = F)
  torlow_cumulative4$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative4$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative4$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative4$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative4$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative4$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative4$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torlow_cumulative4$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torlow_cumulative4$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative4$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                               intervals = int,
                                               partial = F,
                                               total.seeds = torlow_cumulative4[i, 4])
  torlow_cumulative4$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                                  intervals = int,
                                                                  partial = F,
                                                                  total.seeds = torlow_cumulative4[i, 4])
  torlow_cumulative4$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative4[i, 4])
  torlow_cumulative4$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                             intervals = int,
                                                             partial = F)
  torlow_cumulative4$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torlow_cumulative4$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torlow_cumulative4[i, 5:65]),
                                                           intervals = int,
                                                           partial = F)
  
}
torlow_cumulative4

# Torulosa low block 5
for(i in 1:nrow(torlow_cumulative5)){
  torlow_cumulative5$t50[i] <- t50(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                  intervals = int,
                                  partial = F)
  torlow_cumulative5$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative5[i, 4])
  torlow_cumulative5$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torlow_cumulative5$germ_end[i] <- LastGermTime(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                intervals = int,
                                                partial = F)
  torlow_cumulative5$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative5$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative5$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative5$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative5$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative5$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative5$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torlow_cumulative5$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torlow_cumulative5$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative5$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                               intervals = int,
                                               partial = F,
                                               total.seeds = torlow_cumulative5[i, 4])
  torlow_cumulative5$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                                  intervals = int,
                                                                  partial = F,
                                                                  total.seeds = torlow_cumulative5[i, 4])
  torlow_cumulative5$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative5[i, 4])
  torlow_cumulative5$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                             intervals = int,
                                                             partial = F)
  torlow_cumulative5$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torlow_cumulative5$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torlow_cumulative5[i, 5:65]),
                                                           intervals = int,
                                                           partial = F)
  
}
torlow_cumulative5


# Torulosa low block 6
for(i in 1:nrow(torlow_cumulative6)){
  torlow_cumulative6$t50[i] <- t50(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                  intervals = int,
                                  partial = F)
  torlow_cumulative6$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative6[i, 4])
  torlow_cumulative6$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torlow_cumulative6$germ_end[i] <- LastGermTime(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                intervals = int,
                                                partial = F)
  torlow_cumulative6$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative6$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative6$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative6$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative6$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative6$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative6$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torlow_cumulative6$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torlow_cumulative6$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative6$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                               intervals = int,
                                               partial = F,
                                               total.seeds = torlow_cumulative6[i, 4])
  torlow_cumulative6$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                                  intervals = int,
                                                                  partial = F,
                                                                  total.seeds = torlow_cumulative6[i, 4])
  torlow_cumulative6$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative6[i, 4])
  torlow_cumulative6$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                             intervals = int,
                                                             partial = F)
  torlow_cumulative6$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torlow_cumulative6$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torlow_cumulative6[i, 5:65]),
                                                           intervals = int,
                                                           partial = F)
  
}
torlow_cumulative6


# Combine torulosa low blocks

torlow_cumulative <- rbind(torlow_cumulative1, torlow_cumulative2, torlow_cumulative3, torlow_cumulative4, torlow_cumulative5, torlow_cumulative6)
head(torlow_cumulative); tail(torlow_cumulative); dim(torlow_cumulative)




# Torulosa high fire -----
# Split torulosa high fire into blocks
torhigh_cumulative1 <- torhigh_cumulative[1:50,]
torhigh_cumulative2 <- torhigh_cumulative[51:100,]
torhigh_cumulative3 <- torhigh_cumulative[101:150,]
torhigh_cumulative4 <- torhigh_cumulative[151:200,]
torhigh_cumulative5 <- torhigh_cumulative[201:250,]
torhigh_cumulative6 <- torhigh_cumulative[251:nrow(torhigh_cumulative),]

# Torulosa high block 1
for(i in 1:nrow(torhigh_cumulative1)){
  torhigh_cumulative1$t50[i] <- t50(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                  intervals = int,
                                  partial = F)
  torhigh_cumulative1$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative1[i, 4])
  torhigh_cumulative1$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative1$germ_finish[i] <- LastGermTime(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative1$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative1$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative1$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative1$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torhigh_cumulative1$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative1$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative1$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torhigh_cumulative1$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                          intervals = int,
                                                          partial = F)
  torhigh_cumulative1$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative1$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                intervals = int,
                                                partial = F,
                                                total.seeds = torhigh_cumulative1[i, 4])
  torhigh_cumulative1$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                                   intervals = int,
                                                                   partial = F,
                                                                   total.seeds = torhigh_cumulative1[i, 4])
  torhigh_cumulative1$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative1[i, 4])
  torhigh_cumulative1$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                              intervals = int,
                                                              partial = F)
  torhigh_cumulative1$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torhigh_cumulative1$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torhigh_cumulative1[i, 5:65]),
                                                            intervals = int,
                                                            partial = F)
}
torhigh_cumulative1

# Torulosa high block 2
for(i in 1:nrow(torhigh_cumulative2)){
  torhigh_cumulative2$t50[i] <- t50(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                   intervals = int,
                                   partial = F)
  torhigh_cumulative2$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative2[i, 4])
  torhigh_cumulative2$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative2$germ_finish[i] <- LastGermTime(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative2$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative2$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative2$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative2$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative2$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative2$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative2$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torhigh_cumulative2$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                          intervals = int,
                                                          partial = F)
  torhigh_cumulative2$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative2$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                intervals = int,
                                                partial = F,
                                                total.seeds = torhigh_cumulative2[i, 4])
  torhigh_cumulative2$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                                   intervals = int,
                                                                   partial = F,
                                                                   total.seeds = torhigh_cumulative2[i, 4])
  torhigh_cumulative2$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative2[i, 4])
  torhigh_cumulative2$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                              intervals = int,
                                                              partial = F)
  torhigh_cumulative2$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torhigh_cumulative2$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torhigh_cumulative2[i, 5:65]),
                                                            intervals = int,
                                                            partial = F)
}
torhigh_cumulative2

# Torulosa high block 3
for(i in 1:nrow(torhigh_cumulative3)){
  torhigh_cumulative3$t50[i] <- t50(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                   intervals = int,
                                   partial = F)
  torhigh_cumulative3$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative3[i, 4])
  torhigh_cumulative3$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative3$germ_finish[i] <- LastGermTime(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative3$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative3$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative3$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative3$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative3$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative3$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative3$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torhigh_cumulative3$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                          intervals = int,
                                                          partial = F)
  torhigh_cumulative3$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative3$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                intervals = int,
                                                partial = F,
                                                total.seeds = torhigh_cumulative3[i, 4])
  torhigh_cumulative3$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                                   intervals = int,
                                                                   partial = F,
                                                                   total.seeds = torhigh_cumulative3[i, 4])
  torhigh_cumulative3$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative3[i, 4])
  torhigh_cumulative3$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                              intervals = int,
                                                              partial = F)
  torhigh_cumulative3$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torhigh_cumulative3$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torhigh_cumulative3[i, 5:65]),
                                                            intervals = int,
                                                            partial = F)
}
torhigh_cumulative3


# Torulosa high block 4
for(i in 1:nrow(torhigh_cumulative4)){
  torhigh_cumulative4$t50[i] <- t50(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                   intervals = int,
                                   partial = F)
  torhigh_cumulative4$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative4[i, 4])
  torhigh_cumulative4$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative4$germ_finish[i] <- LastGermTime(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative4$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative4$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative4$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative4$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative4$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative4$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative4$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torhigh_cumulative4$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                          intervals = int,
                                                          partial = F)
  torhigh_cumulative4$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative4$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                intervals = int,
                                                partial = F,
                                                total.seeds = torhigh_cumulative4[i, 4])
  torhigh_cumulative4$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                                   intervals = int,
                                                                   partial = F,
                                                                   total.seeds = torhigh_cumulative4[i, 4])
  torhigh_cumulative4$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative4[i, 4])
  torhigh_cumulative4$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                              intervals = int,
                                                              partial = F)
  torhigh_cumulative4$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torhigh_cumulative4$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torhigh_cumulative4[i, 5:65]),
                                                            intervals = int,
                                                            partial = F)
}
torhigh_cumulative4

# Torulosa high block 5
for(i in 1:nrow(torhigh_cumulative5)){
  torhigh_cumulative5$t50[i] <- t50(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                   intervals = int,
                                   partial = F)
  torhigh_cumulative5$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative5[i, 4])
  torhigh_cumulative5$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative5$germ_finish[i] <- LastGermTime(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative5$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative5$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative5$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative5$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative5$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative5$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative5$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torhigh_cumulative5$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                          intervals = int,
                                                          partial = F)
  torhigh_cumulative5$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative5$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                intervals = int,
                                                partial = F,
                                                total.seeds = torhigh_cumulative5[i, 4])
  torhigh_cumulative5$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                                   intervals = int,
                                                                   partial = F,
                                                                   total.seeds = torhigh_cumulative5[i, 4])
  torhigh_cumulative5$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative5[i, 4])
  torhigh_cumulative5$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                              intervals = int,
                                                              partial = F)
  torhigh_cumulative5$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torhigh_cumulative5$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torhigh_cumulative5[i, 5:65]),
                                                            intervals = int,
                                                            partial = F)
}
torhigh_cumulative5

# Torulosa high block 6
for(i in 1:nrow(torhigh_cumulative6)){
  torhigh_cumulative6$t50[i] <- t50(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                   intervals = int,
                                   partial = F)
  torhigh_cumulative6$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative6[i, 4])
  torhigh_cumulative6$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative6$germ_finish[i] <- LastGermTime(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative6$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative6$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative6$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative6$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative6$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative6$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative6$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torhigh_cumulative6$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                          intervals = int,
                                                          partial = F)
  torhigh_cumulative6$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative6$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                intervals = int,
                                                partial = F,
                                                total.seeds = torhigh_cumulative6[i, 4])
  torhigh_cumulative6$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                                   intervals = int,
                                                                   partial = F,
                                                                   total.seeds = torhigh_cumulative6[i, 4])
  torhigh_cumulative6$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative6[i, 4])
  torhigh_cumulative6$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                              intervals = int,
                                                              partial = F)
  torhigh_cumulative6$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torhigh_cumulative6$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torhigh_cumulative6[i, 5:65]),
                                                            intervals = int,
                                                            partial = F)
}
torhigh_cumulative6



# Combine torulosa high fire blocks
torhigh_cumulative <- rbind(torhigh_cumulative1, torhigh_cumulative2, torhigh_cumulative3, torhigh_cumulative4, torhigh_cumulative5, torhigh_cumulative6)
head(torhigh_cumulative); tail(torhigh_cumulative); dim(torhigh_cumulative)

# Make boxplots of the time to 50% germination per treatment group

# Want to do two plots, top plot showing time to 50% germination and then a second plot showing something about the speed of germination, possibly time spread but thinking germination speed would be better.


#### Littoralis
dev.new(height = 20, width = 25, dpi = 80, pointsize = 18, noRStudioGD = T)
par(mar = c(7.5,5,2,1), mgp = c(3,1,0), mfrow = c(2,3), oma = c(2, 2, 3, 0))

# Time to 50%
boxplot(t50 ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 1),], xaxt = "n", xlab = "", ylab = expression(bold("Time to 50% germination")), ylim = c(0,20), cex.lab = 1.25, las = 2)
title('Rep. 1', font.main = 1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)

boxplot(t50 ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 2),], xaxt = "n", xlab = "", ylab = "", las = 2, ylim = c(0,20), cex.lab = 1.25)
title('Rep. 2', font.main = 1)
mtext(expression(bolditalic(Allocasuarina~littoralis)), side = 3, line = 2)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)

boxplot(t50 ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 3),], xaxt = "n", xlab = "", ylab = "", las = 2, ylim = c(0,20), cex.lab = 1.25)
title('Rep. 3', font.main = 1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


# Germination speed
boxplot(germ_speed ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 1),], xaxt = "n", xlab = "", ylab = expression(bold("Germination speed")), ylim = c(0,3), cex.lab = 1.25, las = 2)
title(xlab = expression(bold("Treatment")), line = 6.5, cex.lab = 1.25)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)

boxplot(germ_speed ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 2),], xaxt = "n", xlab = "", ylab = "", ylim = c(0,3), cex.lab = 1.25, las = 2)
title(xlab = expression(bold("Treatment")), line = 6.5, cex.lab = 1.25)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)

boxplot(germ_speed ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 3),], xaxt = "n", xlab = "", ylab = "", ylim = c(0,3), cex.lab = 1.25, las = 2)
title(xlab = expression(bold("Treatment")), line = 6.5, cex.lab = 1.25)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)





### Torulosa low fire
dev.new(height = 20, width = 25, dpi = 80, pointsize = 18, noRStudioGD = T)
par(mar = c(7.5,5,2,1), mgp = c(3,1,0), mfrow = c(2,3), oma = c(2, 2, 3, 0))

# Time to 50%
boxplot(t50 ~ Treatment, data = torlow_cumulative[which(torlow_cumulative$Rep == 1),], xaxt = "n", xlab = "", ylab = expression(bold("Time to 50% germination")), ylim = c(0,20), cex.lab = 1.25, las = 2)
title('Rep. 1', font.main = 1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


boxplot(t50 ~ Treatment, data = torlow_cumulative[which(torlow_cumulative$Rep == 2),], xaxt = "n", xlab = "", ylab = "", ylim = c(0,20), cex.lab = 1.25, las = 2)
mtext(expression(bolditalic(Allocasuarina~torulosa)*bold(' low fire')), side = 3, line = 2)
title('Rep. 2', font.main = 1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


boxplot(t50 ~ Treatment, data = torlow_cumulative[which(torlow_cumulative$Rep == 3),], xaxt = "n", xlab = "", ylab = "", ylim = c(0,20), cex.lab = 1.25, las = 2)
title('Rep. 3', font.main = 1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)



# Germination speed
boxplot(germ_speed ~ Treatment, data = torlow_cumulative[which(torlow_cumulative$Rep == 1),], xaxt = "n", xlab = "", ylab = expression(bold("Germination speed")), ylim = c(0, 3), cex.lab = 1.25, las = 2)
title(xlab = expression(bold("Treatment")), line = 6.5, cex.lab = 1.25)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)

boxplot(germ_speed ~ Treatment, data = torlow_cumulative[which(torlow_cumulative$Rep == 2),],xaxt = "n", xlab = "", ylab = "", las = 2, ylim = c(0,3), cex.lab = 1.25)
title(xlab = expression(bold("Treatment")), line = 6.5, cex.lab = 1.25)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)

boxplot(germ_speed ~ Treatment, data = torlow_cumulative[which(torlow_cumulative$Rep == 3),], xaxt = "n", xlab = "", ylab = "", las = 2, ylim = c(0,3), cex.lab = 1.25)
title(xlab = expression(bold("Treatment")), line = 6.5, cex.lab = 1.25)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)






### Torulosa high fire
dev.new(height = 20, width = 25, dpi = 80, pointsize = 18, noRStudioGD = T)
par(mar = c(7.5,5,2,1), mgp = c(3,1,0), mfrow = c(2,3), oma = c(2, 2, 3, 0))

# Time to 50%                   
boxplot(t50 ~ Treatment, data = torhigh_cumulative[which(torhigh_cumulative$Rep == 1),], xaxt = "n", xlab = "", ylab = expression(bold("Time to 50% germination")), ylim = c(0,20), cex.lab = 1.25, las = 2)
title('Rep. 1', font.main = 1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)

boxplot(t50 ~ Treatment, data = torhigh_cumulative[which(torhigh_cumulative$Rep == 2),], xaxt = "n", xlab = "", ylab = "", ylim = c(0,20), cex.lab = 1.25, las = 2)
mtext(expression(bolditalic(Allocasuarina~torulosa)*bold(' high fire')), side = 3, line = 2)
title('Rep. 2', font.main = 1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)

boxplot(t50 ~ Treatment, data = torhigh_cumulative[which(torhigh_cumulative$Rep == 3),], xaxt = "n", xlab = "", ylab = "", ylim = c(0,20), cex.lab = 1.25, las = 2)
title('Rep. 3', font.main = 1)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


# Germination speed
boxplot(germ_speed ~ Treatment, data = torhigh_cumulative[which(torhigh_cumulative$Rep == 1),], xaxt = "n", xlab = "", ylab = expression(bold("Germination speed")), ylim = c(0,3), cex.lab = 1.25, las = 2)
title(xlab = expression(bold("Treatment")), line = 6.5, cex.lab = 1.25)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


boxplot(germ_speed ~ Treatment, data = torhigh_cumulative[which(torhigh_cumulative$Rep == 2),], xaxt = "n", xlab = "", ylab = "", ylim = c(0,3), cex.lab = 1.25, las = 2)
title(xlab = expression(bold("Treatment")), line = 6.5, cex.lab = 1.25)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


boxplot(germ_speed ~ Treatment, data = torhigh_cumulative[which(torhigh_cumulative$Rep == 3),], xaxt = "n", xlab = "", ylab = "", ylim = c(0,3), cex.lab = 1.25, las = 2)
title(xlab = expression(bold("Treatment")), line = 6.5, cex.lab = 1.25)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)

save.image('./02_Workspaces/Full_experiment_prelim_analysis.RData')
