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
littoralis$Group <- "NA"
str(littoralis)
littoralis <- littoralis[, c(3, 73, 1, 2, 4:ncol(littoralis)-1)]
littoralis <- littoralis[, c(1:4, 6:ncol(littoralis))]
head(littoralis); dim(littoralis)
unique(littoralis$Treatment)

torlow1 <- read.csv('./00_Data/Full_experiment/Set1/torlow1.csv', header = T, stringsAsFactors = T)
torlow2 <- read.csv('./00_Data/Full_experiment/Set2/torlow2.csv', header = T, stringsAsFactors = T)
torlow3 <- read.csv('./00_Data/Full_experiment/Set3/torlow3.csv', header = T, stringsAsFactors = T)
torlow <- rbind(torlow1, torlow2, torlow3)
head(torlow); dim(torlow)
torlow$Species <- "torulosa"
torlow$Group <- "lowfi"
torlow <- torlow[, c(3, 73, 1, 2, 4:ncol(torlow)-1)]
head(torlow); dim(torlow)
torlow <- torlow[, c(1:4, 6:ncol(torlow))]
str(torlow)
torlow$Treatment <- factor(torlow$Treatment, levels = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
unique(torlow$Treatment)


torhigh1 <- read.csv('./00_Data/Full_experiment/Set1/torhigh1.csv', header = T, stringsAsFactors = T)
torhigh2 <- read.csv('./00_Data/Full_experiment/Set2/torhigh2.csv', header = T, stringsAsFactors = T)
torhigh3 <- read.csv('./00_Data/Full_experiment/Set3/torhigh3.csv', header = T, stringsAsFactors = T)
torhigh <- rbind(torhigh1, torhigh2, torhigh3)
head(torhigh); dim(torhigh)
torhigh$Species <- "torulosa"
torhigh$Group <- "hifi"
torhigh <- torhigh[, c(3, 73, 1, 2, 4:ncol(torhigh)-1)]
head(torhigh); dim(torhigh)
torhigh <- torhigh[, c(1:4, 6:ncol(torhigh))]
str(torhigh)
torhigh$Treatment <- factor(torhigh$Treatment, levels = c("Control", "80", "95", "Smoke", "80+smoke", "95+smoke"))
unique(torhigh$Treatment)


# Combine into one dataset 
dat <- rbind(littoralis, torlow, torhigh)
head(dat); tail(dat); dim(dat)

# 3. Calculate the cumulative proportions of germination 
dat_cum.sum <- data.frame(t(apply(dat[ 13:ncol(dat)], 1, FUN = function (x) cumsum(unlist(x)))))
dat_cum.prop <- dat_cum.sum/dat$Total_seeds
dat_cum.prop <- cbind(dat[, 1:12], dat_cum.prop)
dat_cumulative <- cbind(dat[, 1:12], lapply(dat_cum.sum, as.numeric))
head(dat_cum.prop)
head(dat_cumulative)
str(dat_cumulative)


# 4. Boxplot of proportion germinated by treatment and species


dev.new(height = 15, width = 30, dpi = 80, pointsize = 18, noRStudioGD = T)
par(mar = c(10,4.5,2.5,1), mgp = c(3,1,0), mfrow = c(1,3), oma = c(1, 1, 0, 0), cex = 1.1, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

### Littoralis
boxplot(TTC_proportion_viable ~ Treatment, data = dat[which(dat$Species == "littoralis"),], xlab = "", xaxt = "n", ylab = expression(bold("Proportion germinated")), las = 2, pch = 20)
title('(a)', adj = 0, line = 0.5)
title(xlab = expression(bold("Treatment")), line = 9)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)


### Torulosa low fire
boxplot(TTC_proportion_viable ~ Treatment, data = dat[which(dat$Species == "torulosa" & dat$Group == "lowfi"),], xlab = "", xaxt = "n", ylab = "", las = 2, pch = 20)
title('(b)', adj = 0, line = 0.5)
title(xlab = expression(bold("Treatment")), line = 9)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)

### Torulosa high
boxplot(TTC_proportion_viable ~ Treatment, data = dat[which(dat$Species == "torulosa" & dat$Group == "hifi"),], xlab = "", xaxt = "n", ylab = "", las = 2, pch = 20)
title('(c)', adj = 0, line = 0.5)
title(xlab = expression(bold("Treatment")), line = 9)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)







# Germination metrics ----
head(dat_cumulative)

# We cannot get the germinationmetrics::germination.indices functionality to work as it either fails because of some incorrect error thrown when t50 is included, or it decides that the resulting vector would be too long, despite working for the a control only group. We will instead create our own functions to run the calculations as even doing a for loop with the germinationmetrics standalone functions was not working correctly in some instances such as t50. 

# Specify the number of intervals
int <- 1:61

# Calculate 50% germination 
dat_cumulative$Perc50 <- dat_cumulative$Total_germination/2


# Run calculations for each germination metric for each row of data
for(i in 1:nrow(dat_cumulative)){
 dat_cumulative$t50[i] <- paste(colnames(dat_cumulative[i, 13:73])[which(dat_cumulative[i, 13:73] >= dat_cumulative$Perc50)])[1]
 dat_cumulative$t50 <- as.numeric(sub("Day", "", dat_cumulative$t50)) 
 dat_cumulative$firstgerm[i] <- int[min(which(dat_cumulative[i, 13:73] != 0))]
 dat_cumulative$lastgerm[i] <- int[min(which(dat_cumulative[i, 13:73] == dat_cumulative$Total_germination))]
 dat_cumulative$timespread[i] <-  int[min(which(dat_cumulative[i, 13:73] == dat_cumulative$Total_germination))] - int[min(which(dat_cumulative[i, 13:73] != 0))]
 dat_cumulative$germspeed[i] <- sum((dat_cumulative[i, 13:73]/dat_cumulative$Total_seeds)/int)
 dat_cumulative[is.na(dat_cumulative)] <- 0
} 

head(dat_cumulative); tail(dat_cumulative); dim(dat_cumulative)

str(dat_cumulative)

# Create plots for time to 50% germination and germination speed 

dev.new(height = 12, width = 30, dpi = 80, pointsize = 14, noRStudioGD = T)
par(mar = c(14,7,3,1), mgp = c(3,1,0), mfrow = c(1,3), oma = c(0, 0, 0, 0), cex.lab = 3, cex.axis = 2.5, cex.main = 3)

###### TIME TO 50% GERMINATION
### Littoralis
boxplot(t50 ~ Treatment, data = dat_cumulative[dat_cumulative$Species == "littoralis",], xaxt = "n", xlab = "", yaxt = "n", ylab = "", ylim = c(0,24), pch = 19, cex = 1.5)
title('(a)', adj = 0, line = 0.5)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)
title(ylab = expression(bold("Time to 50%")), line = 4)
axis(side = 2, at = seq(from = 0, to = 30, by = 2), las = 1)
title(xlab = expression(bold("Treatment")), line = 13, adj = 0.45)


### Torulosa low fire
boxplot(t50 ~ Treatment, data = dat_cumulative[dat_cumulative$Species == "torulosa" & dat_cumulative$Group == "lowfi",], xaxt = "n", xlab = "", yaxt = "n", ylab = "", ylim = c(0,24), pch = 19, cex = 1.5)
title('(b)', adj = 0, line = 0.5)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)
axis(side = 2, at = seq(from = 0, to = 30, by = 2), las = 1)
title(xlab = expression(bold("Treatment")), line = 13, adj = 0.45)


### Torulosa high fire
boxplot(t50 ~ Treatment, data = dat_cumulative[dat_cumulative$Species == "torulosa" & dat_cumulative$Group == "hifi",], xaxt = "n", xlab = "", yaxt = "n", ylab = "", ylim = c(0,24), pch = 19, cex = 1.5)
title('(c)', adj = 0, line = 0.5)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)
axis(side = 2, at = seq(from = 0, to = 30, by = 2), las = 1)
title(xlab = expression(bold("Treatment")), line = 13, adj = 0.45)





###### GERMINATION SPEED
dev.new(height = 12, width = 30, dpi = 80, pointsize = 14, noRStudioGD = T)
par(mar = c(14,7,3,1), mgp = c(3,1,0), mfrow = c(1,3), oma = c(0, 0, 0, 0), cex.lab = 3, cex.axis = 2.5, cex.main = 3)


### Littoralis
boxplot(germspeed ~ Treatment, data = dat_cumulative[dat_cumulative$Species == "littoralis",], xaxt = "n", xlab = "", ylab = "", ylim = c(0, 2.5), las = 1, pch = 19, cex = 1.5)
title('(a)', adj = 0, line = 0.5)
title(xlab = expression(bold("Treatment")), line = 12, adj = 0.45)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)
title(ylab = expression(bold("Germination speed")), line = 4)
rug(x = c(0:2.5 + 0.25, 0:2.5 + 0.75), ticksize = -0.025, side = 2)


### Torulosa low fire
boxplot(germspeed ~ Treatment, data = dat_cumulative[dat_cumulative$Species == "torulosa" & dat_cumulative$Group == "lowfi",], xaxt = "n", xlab = "", ylab = "", ylim = c(0, 2.5), las = 1, pch = 19, cex = 1.5)
title('(b)', adj = 0, line = 0.5)
title(xlab = expression(bold("Treatment")), line = 12, adj = 0.45)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)
rug(x = c(0:2.5 + 0.25, 0:2.5 + 0.75), ticksize = -0.025, side = 2)


### Torulosa high fire
boxplot(germspeed ~ Treatment, data = dat_cumulative[dat_cumulative$Species == "torulosa" & dat_cumulative$Group == "hifi",], xaxt = "n", xlab = "", ylab = "", ylim = c(0, 2.5), las = 1, pch = 19, cex = 1.5)
title('(c)', adj = 0, line = 0.5)
title(xlab = expression(bold("Treatment")), line = 12, adj = 0.45)
axis(side = 1, at = c(1,2,3,4,5,6), labels = c("Control", "80°C", "95°C", "Smoke", "80°C+smoke", "95°C+smoke"), las = 2)
rug(x = c(0:2.5 + 2.5, 0:2.5 + 0.75), ticksize = -0.025, side = 2)


save.image('./02_Workspaces/Full_experiment_prelim_analysis.RData')



