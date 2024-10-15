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


# Throws an error message but ignored this as the expected results are returned
for(i in 1:nrow(littoralis_cumulative)){
  
  littoralis_cumulative$t50[i] <- t50(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative[i, 4])
  littoralis_cumulative$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$germ_end[i] <- LastGermTime(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F, percent = F, total.seeds = littoralis_cumulative[i, 4])
  littoralis_cumulative$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative[i, 4])
  littoralis_cumulative$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F, total.seeds = littoralis_cumulative[i, 4])
  littoralis_cumulative$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  littoralis_cumulative$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(littoralis_cumulative[i, 5:65]), intervals = int, partial = F)
  
  
}

littoralis_cumulative 


# Check if a couple of the rows match between the for loop version or a row by row calculation where we do not get an error message
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
littoralis_cumulative[72,] # Some values do change slightly, most we will probably use did not change. Could be related to the multiple peaks in some cases.



# We are going to say that this for loop calculation works how it should, if we think any values are out of range of what we expect we can always recalculate for the specific rows that are problematic.

# Torulosa low fire ----
for(i in 1:nrow(torlow_cumulative)){
  torlow_cumulative$t50[i] <- t50(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                  intervals = int,
                                  partial = F)
  torlow_cumulative$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative[i, 4])
  torlow_cumulative$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torlow_cumulative$germ_end[i] <- LastGermTime(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                intervals = int,
                                                partial = F)
  torlow_cumulative$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torlow_cumulative$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torlow_cumulative$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                      intervals = int,
                                                      partial = F)
  torlow_cumulative$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torlow_cumulative$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torlow_cumulative$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                               intervals = int,
                                               partial = F)
  torlow_cumulative$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                            intervals = int,
                                            partial = F,
                                            total.seeds = torlow_cumulative[i, 4])
  torlow_cumulative$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                                  intervals = int,
                                                                  partial = F,
                                                                  total.seeds = torlow_cumulative[i, 4])
  torlow_cumulative$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                            intervals = int,
                                                            partial = F,
                                                            total.seeds = torlow_cumulative[i, 4])
  torlow_cumulative$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                          intervals = int,
                                                          partial = F)
  torlow_cumulative$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torlow_cumulative$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torlow_cumulative[i, 5:65]),
                                                           intervals = int,
                                                           partial = F)
  
}
torlow_cumulative



# Torulosa high fire -----
for(i in 1:nrow(torhigh_cumulative)){
  torhigh_cumulative$t50[i] <- t50(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                  intervals = int,
                                  partial = F)
  torhigh_cumulative$peak_germ_percent[i] <- PeakGermPercent(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative[i, 4])
  torhigh_cumulative$germ_start[i] <- FirstGermTime(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative$germ_finish[i] <- LastGermTime(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                    intervals = int,
                                                    partial = F)
  torhigh_cumulative$timespread[i] <- TimeSpreadGerm(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative$mean_germ_time[i] <- MeanGermTime(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative$var_germ_time[i] <- VarGermTime(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                     intervals = int,
                                                     partial = F)
  torhigh_cumulative$se_germ_time[i] <- SEGermTime(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                  intervals = int,
                                                  partial = F)
  torhigh_cumulative$cv_germ_time[i] <- CVGermTime(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative$mean_germ_rate[i] <- MeanGermRate(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                       intervals = int,
                                                       partial = F)
  torhigh_cumulative$coefficient_velocity_germ[i] <- CVG(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                         intervals = int,
                                                         partial = F)
  torhigh_cumulative$variance_germ_rate[i] <- VarGermRate(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                          intervals = int,
                                                          partial = F)
  torhigh_cumulative$se_germ_rate[i] <- SEGermRate(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                   intervals = int,
                                                   partial = F)
  torhigh_cumulative$germ_speed[i] <- GermSpeed(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                intervals = int,
                                                partial = F,
                                                total.seeds = torhigh_cumulative[i, 4])
  torhigh_cumulative$weighted_germ_percent[i] <- WeightGermPercent(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                                   intervals = int,
                                                                   partial = F,
                                                                   total.seeds = torhigh_cumulative[i, 4])
  torhigh_cumulative$mean_germ_percent[i] <- MeanGermPercent(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                             intervals = int,
                                                             partial = F,
                                                             total.seeds = torhigh_cumulative[i, 4])
  torhigh_cumulative$coefficient_uniformity_germ[i] <- CUGerm(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                              intervals = int,
                                                              partial = F)
  torhigh_cumulative$germ_synchrony[i] <- GermSynchrony(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                        intervals = int,
                                                        partial = F)
  torhigh_cumulative$germ_uncertainty[i] <- GermUncertainty(germ.counts = as.numeric(torhigh_cumulative[i, 5:65]),
                                                            intervals = int,
                                                            partial = F)
}

torhigh_cumulative




# Make boxplots of the time to 50% germination per treatment group

# Want to do two plots, top plot showing time to 50% germination and then a second plot showing the timespread of germination for each replicate and each species grouping

dev.new(height = 20, width = 35, dpi = 80, pointsize = 18, noRStudioGD = T)
par(mar = c(1,4,2,1), mgp = c(3,1,0), mfrow = c(2,3), oma = c(6, 0.5, 3, 0))

boxplot(t50 ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 1),], xaxt = "n", xlab = "", ylab = expression(bold("Time to 50% germination")), las = 2)
title('Rep. 1', font.main = 1)

boxplot(t50 ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 2),], xaxt = "n", xlab = "", ylab = "", las = 2, ylim = c(0,20))
title('Rep. 2', font.main = 1)

boxplot(t50 ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 3),], xaxt = "n", xlab = "", ylab = "", las = 2, ylim = c(0,20))
title('Rep. 3', font.main = 1)




boxplot(timespread ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 1),], xlab = "", ylab = expression(bold("Time spread")), las = 2, ylim = c(0,20))
title(xlab = "Treatment", line = 5)

boxplot(timespread ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 2),], xlab = "", ylab = "", las = 2, ylim = c(0,20))
title(xlab = "Treatment", line = 5)

boxplot(timespread ~ Treatment, data = littoralis_cumulative[which(littoralis_cumulative$Rep == 3),], xlab = "", ylab = "", las = 2, ylim = c(0,20))
title(xlab = "Treatment", line = 5)


save.image('./02_Workspaces/Full_experiment_prelim_analysis.RData')
