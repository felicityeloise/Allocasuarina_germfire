# Written by Felicity Charles
# Caveat emptor
# Date: 12th June 2024

# Optimisation testing - analysis of the test treatments for seed germination heat shock temperatures and smoke application
library(germinationmetrics)

# 1. Heat shock temperature optimisation
heat <- read.table('./00_Data/Optimisation_tests/heat_tests.txt', header = T)
head(heat); dim(heat)


heat$Time_spread <- heat$Time_to_finish - heat$Time_to_germ

heat.seed <- heat[,which(colnames(heat) == "X2"): which(colnames(heat) == "X39")]

heat_cum.sum <- data.frame(t(apply(heat.seed, 1, FUN = function(x) cumsum(unlist(x)))))
heat_cum.prop <- heat_cum.sum/heat$Total_seeds_adjusted
heat_cum.prop2<- cbind(heat[,c(1:4, 25)], heat_cum.prop)
head(heat_cum.prop2); dim(heat_cum.prop2)
str(heat_cum.prop2)

# Split the data by species for plotting
lit_heat <- subset(heat_cum.prop2, heat_cum.prop2$Species == "littoralis")
rownames(lit_heat) <- 1:nrow(lit_heat)
head(lit_heat); dim(lit_heat)

tor_heat <- subset(heat_cum.prop2, heat_cum.prop2$Species == "torulosa")
rownames(tor_heat) <- 1:nrow(tor_heat)
head(tor_heat);dim(tor_heat)


# Set colour palette to be used for plotting
# We have 5 temperatures and the control so six lines total
# Would like control to be plotted in black, with the others ranging from the coolest temp in blue to hottest temp in red

library(RColorBrewer)
brewer.pal(5, "RdYlBu")
pal <- c("black", '#ABD9E9', '#2C7BB6', '#FFFFBF', '#FDAE61', "#D7191C")
duration.col <- data.frame(col = pal, duration= unique(heat$Duration))

# Extract the data just relating to the proportion germination
prop.heat_lit <- lit_heat[,which(colnames(lit_heat) == "X2"): which(colnames(lit_heat) == "X39")]
prop.heat_tor <- tor_heat[, which(colnames(tor_heat) == "X2"): which(colnames(tor_heat) == "X39")]

# Subset the data further by heat treatment
# For the purposes of plotting we can ignore the temperatures of 125 and 150 degrees as we either had no germination at these temperatures or only 1 or 2 seeds which germinated so we won't be pursuing these further. 

lit_heat1 <- subset(lit_heat, lit_heat$Temperature == '0' | lit_heat$Temperature == '80')
lit_heat2 <- subset(lit_heat, lit_heat$Temperature == '0' | lit_heat$Temperature == '95')
lit_heat3 <- subset(lit_heat, lit_heat$Temperature == '0' | lit_heat$Temperature == '110')

tor_heat1 <- subset(tor_heat, tor_heat$Temperature == '0' | tor_heat$Temperature == '80')
tor_heat2 <- subset(tor_heat, tor_heat$Temperature == '0' | tor_heat$Temperature == '95')
tor_heat3 <- subset(tor_heat, tor_heat$Temperature == '0' | tor_heat$Temperature == '110')



# Create the plots for heat shocks

dev.new(height = 10, width = 18, dpi = 80, pointsize = 18, noRStudioGD = T)
par(mar = c(6,4,2,1), mgp = c(3,1,0), mfrow = c(2,3), oma = c(0,0, 3, 8))

plot(1:ncol(prop.heat_lit), 1:ncol(prop.heat_lit), ylim = c(0,1), type = "n", las = 1, ylab = "Proportion germinated", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.heat_lit), labels = c(2,4,7,9,11,14,16,18,21,23,25,28,30,32,35,37,39), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(main = "80°C", font.main = 1)
for(i in 1:nrow(lit_heat1)){
  dat.thisrun <- lit_heat1[i,6:ncol(lit_heat1)]
  duration.thisrun <- lit_heat1$Duration[i]
  col.thisrun <- duration.col$col[duration.col$duration == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}


plot(1:ncol(prop.heat_lit), 1:ncol(prop.heat_lit), ylim = c(0,1), type = "n", las = 1, ylab = "", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.heat_lit), labels = c(2,4,7,9,11,14,16,18,21,23,25,28,30,32,35,37,39), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(main = "95°C", font.main = 1)
for(i in 1:nrow(lit_heat2)){
  dat.thisrun <- lit_heat2[i,6:ncol(lit_heat2)]
  duration.thisrun <- lit_heat2$Duration[i]
  col.thisrun <- duration.col$col[duration.col$duration == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}


plot(1:ncol(prop.heat_lit), 1:ncol(prop.heat_lit), ylim = c(0,1), type = "n", las = 1, ylab = "", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.heat_lit), labels = c(2,4,7,9,11,14,16,18,21,23,25,28,30,32,35,37,39), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(main = "110°C", font.main = 1)
for(i in 1:nrow(lit_heat3)){
  dat.thisrun <- lit_heat3[i,6:ncol(lit_heat3)]
  duration.thisrun <- lit_heat2$Duration[i]
  col.thisrun <- duration.col$col[duration.col$duration == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}
mtext(expression(bold('(a) ')~bolditalic(Allocasuarina~littoralis)), side = 3, outer = T, font = )

mtext(expression(bold('(b) ')~bolditalic(Allocasuarina~torulosa)), side = 3, line = -24, outer = T, font = 2)


plot(1:ncol(prop.heat_tor), 1:ncol(prop.heat_tor), ylim = c(0,1), type = "n", las = 1, ylab = "Proportion germinated", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.heat_tor), labels = c(2,4,7,9,11,14,16,18,21,23,25,28,30,32,35,37,39), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(xlab = "Days", mgp = c(2.2, 1, 0))
title(main = "80°C", font.main = 1)
for(i in 1:nrow(tor_heat1)){
  dat.thisrun <- tor_heat1[i,6:ncol(tor_heat1)]
  duration.thisrun <- tor_heat1$Duration[i]
  col.thisrun <- duration.col$col[duration.col$duration == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}


plot(1:ncol(prop.heat_tor), 1:ncol(prop.heat_tor), ylim = c(0,1), type = "n", las = 1, ylab = '', xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.heat_tor), labels = c(2,4,7,9,11,14,16,18,21,23,25,28,30,32,35,37,39), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(xlab = "Days", mgp = c(2.2, 1, 0))
title(main = "95°C", font.main = 1)
for(i in 1:nrow(tor_heat2)){
  dat.thisrun <- tor_heat2[i, 6:ncol(tor_heat2)]
  duration.thisrun <- tor_heat2$Duration[i]
  col.thisrun <- duration.col$col[duration.col$duration == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}

plot(1:ncol(prop.heat_tor), 1:ncol(prop.heat_tor), ylim = c(0,1), type = "n", las = 1, ylab = "", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.heat_tor), labels = c(2,4,7,9,11,14,16,18,21,23,25,28,30,32,35,37,39), cex.axis = 0.8, mgp = c(2.2, 0.7,0))
title(xlab = "Days", mgp = c(2.2, 1, 0))
title(main = "110°C", font.main = 1)
for(i in 1:nrow(tor_heat3)){
  dat.thisrun <- tor_heat3[i, 6:ncol(tor_heat3)]
  duration.thisrun  <- tor_heat3$Duration[i]
  col.thisrun <- duration.col$col[duration.col$duration == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}

par(xpd = NA)
legend(x = 18, y = 2.5, legend = c('Control', '30 sec', '1 min', '2 min', '5 min', '10 min'), col = duration.col$col, lty = 1, lwd = 4, cex = 0.95, text.width = 0.2, title = 'Exposure duration', bty = "n")
par(xpd = F)

# Maybe we want to re plot the data as we decide against particular durations to make it easier to see what is going on and make further decisions. Also want to think about delving into this further, looking at Time to 50% germination and those other germination metrics to back up an decisions.



# Look at some summary stats ----
heat_long <- read.table('./00_Data/Optimisation_tests/heat_long.txt', header = T)
head(heat_long);dim(heat_long)

# Split the data for analysing
# Controls
T1_con <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "Control")
T2_con <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "Control")
WPS73_con <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "Control")

WPS89_con <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "Control")
WPS107_con <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "Control")
WPS154_con <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "Control")



# 80 deg 30 sec
T1_80_30 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "80" & heat_long$Duration == '30')
T2_80_30 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "80" & heat_long$Duration == '30')
WPS73_80_30 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "80" & heat_long$Duration == '30')

WPS89_80_30 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "80" & heat_long$Duration == '30')
WPS107_80_30 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "80" & heat_long$Duration == '30')
WPS154_80_30 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "80" & heat_long$Duration == '30')



# 80 deg 1 min
T1_80_1 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "80" & heat_long$Duration == '1')
T2_80_1 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "80" & heat_long$Duration == '1')
WPS73_80_1 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "80" & heat_long$Duration == '1')


WPS89_80_1 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "80" & heat_long$Duration == '1')
WPS107_80_1 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "80" & heat_long$Duration == '1')
WPS154_80_1 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "80" & heat_long$Duration == '1')


# 80 deg 2 min
T1_80_2 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "80" & heat_long$Duration == '2')
T2_80_2 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "80" & heat_long$Duration == '2')
WPS73_80_2 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "80" & heat_long$Duration == '2')


WPS89_80_2 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "80" & heat_long$Duration == '2')
WPS107_80_2 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "80" & heat_long$Duration == '2')
WPS154_80_2 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "80" & heat_long$Duration == '2')



# 80 deg 5 min
T1_80_5 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "80" & heat_long$Duration == '30')
T2_80_5 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "80" & heat_long$Duration == '30')
WPS73_80_5 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "80" & heat_long$Duration == '30')


WPS89_80_5 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "80" & heat_long$Duration == '30')
WPS107_80_5 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "80" & heat_long$Duration == '30')
WPS154_80_5 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "80" & heat_long$Duration == '30')



# 80 deg 10 min
T1_80_10 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "80" & heat_long$Duration == '10')
T2_80_10 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "80" & heat_long$Duration == '10')
WPS73_80_10 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "80" & heat_long$Duration == '10')


WPS89_80_10 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "80" & heat_long$Duration == '10')
WPS107_80_10 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "80" & heat_long$Duration == '10')
WPS154_80_10 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "80" & heat_long$Duration == '10')




# 95 deg 30 sec
T1_95_30 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "95" & heat_long$Duration == '30')
T2_95_30 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "95" & heat_long$Duration == '30')
WPS73_95_30 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "95" & heat_long$Duration == '30')

WPS89_95_30 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "95" & heat_long$Duration == '30')
WPS107_95_30 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "95" & heat_long$Duration == '30')
WPS154_95_30 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "95" & heat_long$Duration == '30')



# 95 deg 1 min
T1_95_1 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "95" & heat_long$Duration == '1')
T2_95_1 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "95" & heat_long$Duration == '1')
WPS73_95_1 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "95" & heat_long$Duration == '1')


WPS89_95_1 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "95" & heat_long$Duration == '1')
WPS107_95_1 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "95" & heat_long$Duration == '1')
WPS154_95_1 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "95" & heat_long$Duration == '1')


# 95 deg 2 min
T1_95_2 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "95" & heat_long$Duration == '2')
T2_95_2 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "95" & heat_long$Duration == '2')
WPS73_95_2 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "95" & heat_long$Duration == '2')


WPS89_95_2 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "95" & heat_long$Duration == '2')
WPS107_95_2 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "95" & heat_long$Duration == '2')
WPS154_95_2 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "95" & heat_long$Duration == '2')



# 95 deg 5 min
T1_95_5 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "95" & heat_long$Duration == '30')
T2_95_5 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "95" & heat_long$Duration == '30')
WPS73_95_5 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "95" & heat_long$Duration == '30')


WPS89_95_5 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "95" & heat_long$Duration == '30')
WPS107_95_5 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "95" & heat_long$Duration == '30')
WPS154_95_5 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "95" & heat_long$Duration == '30')



# 95 deg 10 min
T1_95_10 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "95" & heat_long$Duration == '10')
T2_95_10 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "95" & heat_long$Duration == '10')
WPS73_95_10 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "95" & heat_long$Duration == '10')


WPS89_95_10 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "95" & heat_long$Duration == '10')
WPS107_95_10 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "95" & heat_long$Duration == '10')
WPS154_95_10 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "95" & heat_long$Duration == '10')




# 110 deg 30 sec
T1_110_30 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "110" & heat_long$Duration == '30')
T2_110_30 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "110" & heat_long$Duration == '30')
WPS73_110_30 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "110" & heat_long$Duration == '30')

WPS89_110_30 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "110" & heat_long$Duration == '30')
WPS107_110_30 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "110" & heat_long$Duration == '30')
WPS154_110_30 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "110" & heat_long$Duration == '30')



# 110 deg 1 min
T1_110_1 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "110" & heat_long$Duration == '1')
T2_110_1 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "110" & heat_long$Duration == '1')
WPS73_110_1 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "110" & heat_long$Duration == '1')


WPS89_110_1 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "110" & heat_long$Duration == '1')
WPS107_110_1 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "110" & heat_long$Duration == '1')
WPS154_110_1 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "110" & heat_long$Duration == '1')


# 110 deg 2 min
T1_110_2 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "110" & heat_long$Duration == '2')
T2_110_2 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "110" & heat_long$Duration == '2')
WPS73_110_2 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "110" & heat_long$Duration == '2')


WPS89_110_2 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "110" & heat_long$Duration == '2')
WPS107_110_2 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "110" & heat_long$Duration == '2')
WPS154_110_2 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "110" & heat_long$Duration == '2')



# 110 deg 5 min
T1_110_5 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "110" & heat_long$Duration == '30')
T2_110_5 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "110" & heat_long$Duration == '30')
WPS73_110_5 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "110" & heat_long$Duration == '30')


WPS89_110_5 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "110" & heat_long$Duration == '30')
WPS107_110_5 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "110" & heat_long$Duration == '30')
WPS154_110_5 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "110" & heat_long$Duration == '30')



# 110 deg 10 min
T1_110_10 <- subset(heat_long, heat_long$Seed_lot == "T1" & heat_long$Temperature == "110" & heat_long$Duration == '10')
T2_110_10 <- subset(heat_long, heat_long$Seed_lot == "T2" & heat_long$Temperature == "110" & heat_long$Duration == '10')
WPS73_110_10 <- subset(heat_long, heat_long$Seed_lot == "WPS73" & heat_long$Temperature == "110" & heat_long$Duration == '10')


WPS89_110_10 <- subset(heat_long, heat_long$Seed_lot == "WPS89" & heat_long$Temperature == "110" & heat_long$Duration == '10')
WPS107_110_10 <- subset(heat_long, heat_long$Seed_lot == "WPS107" & heat_long$Temperature == "110" & heat_long$Duration == '10')
WPS154_110_10 <- subset(heat_long, heat_long$Seed_lot == "WPS154" & heat_long$Temperature == "110" & heat_long$Duration == '10')


# Analyse the data
head(T1_con)
GermPercent(germ.counts = T1_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_con$Count, intervals = T1_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_con$Count, intervals = T1_con$Interval, partial = T)
t50(germ.counts = T1_con$Count, intervals = T1_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_con$Count, intervals = T1_con$Interval, partial = T)
WeightGermPercent(germ.counts = T1_con$Count, total.seeds = 20, intervals = T1_con$Interval, partial = T)

head(T2_con)
GermPercent(germ.counts = T2_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_con$Count, intervals = T2_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_con$Count, intervals = T2_con$Interval, partial = T)
t50(germ.counts = T2_con$Count, intervals = T2_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_con$Count, intervals = T2_con$Interval, partial = T)
WeightGermPercent(germ.counts = T2_con$Count, total.seeds = 20, intervals = T2_con$Interval, partial = T)


head(WPS73_con)
GermPercent(germ.counts = WPS73_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_con$Count, intervals = WPS73_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_con$Count, intervals = WPS73_con$Interval, partial = T)
t50(germ.counts = WPS73_con$Count, intervals = WPS73_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_con$Count, intervals = WPS73_con$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_con$Count, total.seeds = 20, intervals = WPS73_con$Interval, partial = T)




head(T1_80_30)
GermPercent(germ.counts = T1_80_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_80_30$Count, intervals = T1_80_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_80_30$Count, intervals = T1_80_30$Interval, partial = T)
t50(germ.counts = T1_80_30$Count, intervals = T1_80_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_80_30$Count, intervals = T1_80_30$Interval, partial = T)
WeightGermPercent(germ.counts = T1_80_30$Count, total.seeds = 20, intervals = T1_80_30$Interval, partial = T)

head(T2_80_30)
GermPercent(germ.counts = T2_80_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_80_30$Count, intervals = T2_80_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_80_30$Count, intervals = T2_80_30$Interval, partial = T)
t50(germ.counts = T2_80_30$Count, intervals = T2_80_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_80_30$Count, intervals = T2_80_30$Interval, partial = T)
WeightGermPercent(germ.counts = T2_80_30$Count, total.seeds = 20, intervals = T2_80_30$Interval, partial = T)


head(WPS73_80_30)
GermPercent(germ.counts = WPS73_80_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_80_30$Count, intervals = WPS73_80_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_80_30$Count, intervals = WPS73_80_30$Interval, partial = T)
t50(germ.counts = WPS73_80_30$Count, intervals = WPS73_80_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_80_30$Count, intervals = WPS73_80_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_80_30$Count, total.seeds = 20, intervals = WPS73_80_30$Interval, partial = T)



head(T1_80_1)
GermPercent(germ.counts = T1_80_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_80_1$Count, intervals = T1_80_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_80_1$Count, intervals = T1_80_1$Interval, partial = T)
t50(germ.counts = T1_80_1$Count, intervals = T1_80_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_80_1$Count, intervals = T1_80_1$Interval, partial = T)
WeightGermPercent(germ.counts = T1_80_1$Count, total.seeds = 20, intervals = T1_80_1$Interval, partial = T)

head(T2_80_1)
GermPercent(germ.counts = T2_80_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_80_1$Count, intervals = T2_80_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_80_1$Count, intervals = T2_80_1$Interval, partial = T)
t50(germ.counts = T2_80_1$Count, intervals = T2_80_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_80_1$Count, intervals = T2_80_1$Interval, partial = T)
WeightGermPercent(germ.counts = T2_80_1$Count, total.seeds = 20, intervals = T2_80_1$Interval, partial = T)


head(WPS73_80_1)
GermPercent(germ.counts = WPS73_80_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_80_1$Count, intervals = WPS73_80_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_80_1$Count, intervals = WPS73_80_1$Interval, partial = T)
t50(germ.counts = WPS73_80_1$Count, intervals = WPS73_80_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_80_1$Count, intervals = WPS73_80_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_80_1$Count, total.seeds = 20, intervals = WPS73_80_1$Interval, partial = T)



head(T1_80_2)
GermPercent(germ.counts = T1_80_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_80_2$Count, intervals = T1_80_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_80_2$Count, intervals = T1_80_2$Interval, partial = T)
t50(germ.counts = T1_80_2$Count, intervals = T1_80_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_80_2$Count, intervals = T1_80_2$Interval, partial = T)
WeightGermPercent(germ.counts = T1_80_2$Count, total.seeds = 20, intervals = T1_80_2$Interval, partial = T)

head(T2_80_2)
GermPercent(germ.counts = T2_80_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_80_2$Count, intervals = T2_80_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_80_2$Count, intervals = T2_80_2$Interval, partial = T)
t50(germ.counts = T2_80_2$Count, intervals = T2_80_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_80_2$Count, intervals = T2_80_2$Interval, partial = T)
WeightGermPercent(germ.counts = T2_80_2$Count, total.seeds = 20, intervals = T2_80_2$Interval, partial = T)


head(WPS73_80_2)
GermPercent(germ.counts = WPS73_80_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_80_2$Count, intervals = WPS73_80_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_80_2$Count, intervals = WPS73_80_2$Interval, partial = T)
t50(germ.counts = WPS73_80_2$Count, intervals = WPS73_80_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_80_2$Count, intervals = WPS73_80_2$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_80_2$Count, total.seeds = 20, intervals = WPS73_80_2$Interval, partial = T)



head(T1_80_5)
GermPercent(germ.counts = T1_80_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_80_5$Count, intervals = T1_80_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_80_5$Count, intervals = T1_80_5$Interval, partial = T)
t50(germ.counts = T1_80_5$Count, intervals = T1_80_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_80_5$Count, intervals = T1_80_5$Interval, partial = T)
WeightGermPercent(germ.counts = T1_80_5$Count, total.seeds = 20, intervals = T1_80_5$Interval, partial = T)

head(T2_80_5)
GermPercent(germ.counts = T2_80_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_80_5$Count, intervals = T2_80_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_80_5$Count, intervals = T2_80_5$Interval, partial = T)
t50(germ.counts = T2_80_5$Count, intervals = T2_80_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_80_5$Count, intervals = T2_80_5$Interval, partial = T)
WeightGermPercent(germ.counts = T2_80_5$Count, total.seeds = 20, intervals = T2_80_5$Interval, partial = T)


head(WPS73_80_5)
GermPercent(germ.counts = WPS73_80_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_80_5$Count, intervals = WPS73_80_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_80_5$Count, intervals = WPS73_80_5$Interval, partial = T)
t50(germ.counts = WPS73_80_5$Count, intervals = WPS73_80_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_80_5$Count, intervals = WPS73_80_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_80_5$Count, total.seeds = 20, intervals = WPS73_80_5$Interval, partial = T)



head(T1_80_10)
GermPercent(germ.counts = T1_80_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_80_10$Count, intervals = T1_80_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_80_10$Count, intervals = T1_80_10$Interval, partial = T)
t50(germ.counts = T1_80_10$Count, intervals = T1_80_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_80_10$Count, intervals = T1_80_10$Interval, partial = T)
WeightGermPercent(germ.counts = T1_80_10$Count, total.seeds = 20, intervals = T1_80_10$Interval, partial = T)

head(T2_80_10)
GermPercent(germ.counts = T2_80_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_80_10$Count, intervals = T2_80_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_80_10$Count, intervals = T2_80_10$Interval, partial = T)
t50(germ.counts = T2_80_10$Count, intervals = T2_80_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_80_10$Count, intervals = T2_80_10$Interval, partial = T)
WeightGermPercent(germ.counts = T2_80_10$Count, total.seeds = 20, intervals = T2_80_10$Interval, partial = T)


head(WPS73_80_10)
GermPercent(germ.counts = WPS73_80_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_80_10$Count, intervals = WPS73_80_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_80_10$Count, intervals = WPS73_80_10$Interval, partial = T)
t50(germ.counts = WPS73_80_10$Count, intervals = WPS73_80_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_80_10$Count, intervals = WPS73_80_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_80_10$Count, total.seeds = 20, intervals = WPS73_80_10$Interval, partial = T)





head(T1_95_30)
GermPercent(germ.counts = T1_95_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_95_30$Count, intervals = T1_95_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_95_30$Count, intervals = T1_95_30$Interval, partial = T)
t50(germ.counts = T1_95_30$Count, intervals = T1_95_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_95_30$Count, intervals = T1_95_30$Interval, partial = T)
WeightGermPercent(germ.counts = T1_95_30$Count, total.seeds = 20, intervals = T1_95_30$Interval, partial = T)

head(T2_95_30)
GermPercent(germ.counts = T2_95_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_95_30$Count, intervals = T2_95_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_95_30$Count, intervals = T2_95_30$Interval, partial = T)
t50(germ.counts = T2_95_30$Count, intervals = T2_95_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_95_30$Count, intervals = T2_95_30$Interval, partial = T)
WeightGermPercent(germ.counts = T2_95_30$Count, total.seeds = 20, intervals = T2_95_30$Interval, partial = T)


head(WPS73_95_30)
GermPercent(germ.counts = WPS73_95_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_95_30$Count, intervals = WPS73_95_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_95_30$Count, intervals = WPS73_95_30$Interval, partial = T)
t50(germ.counts = WPS73_95_30$Count, intervals = WPS73_95_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_95_30$Count, intervals = WPS73_95_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_95_30$Count, total.seeds = 20, intervals = WPS73_95_30$Interval, partial = T)



head(T1_95_1)
GermPercent(germ.counts = T1_95_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_95_1$Count, intervals = T1_95_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_95_1$Count, intervals = T1_95_1$Interval, partial = T)
t50(germ.counts = T1_95_1$Count, intervals = T1_95_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_95_1$Count, intervals = T1_95_1$Interval, partial = T)
WeightGermPercent(germ.counts = T1_95_1$Count, total.seeds = 20, intervals = T1_95_1$Interval, partial = T)

head(T2_95_1)
GermPercent(germ.counts = T2_95_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_95_1$Count, intervals = T2_95_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_95_1$Count, intervals = T2_95_1$Interval, partial = T)
t50(germ.counts = T2_95_1$Count, intervals = T2_95_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_95_1$Count, intervals = T2_95_1$Interval, partial = T)
WeightGermPercent(germ.counts = T2_95_1$Count, total.seeds = 20, intervals = T2_95_1$Interval, partial = T)


head(WPS73_95_1)
GermPercent(germ.counts = WPS73_95_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_95_1$Count, intervals = WPS73_95_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_95_1$Count, intervals = WPS73_95_1$Interval, partial = T)
t50(germ.counts = WPS73_95_1$Count, intervals = WPS73_95_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_95_1$Count, intervals = WPS73_95_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_95_1$Count, total.seeds = 20, intervals = WPS73_95_1$Interval, partial = T)



head(T1_95_2)
GermPercent(germ.counts = T1_95_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_95_2$Count, intervals = T1_95_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_95_2$Count, intervals = T1_95_2$Interval, partial = T)
t50(germ.counts = T1_95_2$Count, intervals = T1_95_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_95_2$Count, intervals = T1_95_2$Interval, partial = T)
WeightGermPercent(germ.counts = T1_95_2$Count, total.seeds = 20, intervals = T1_95_2$Interval, partial = T)

head(T2_95_2)
GermPercent(germ.counts = T2_95_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_95_2$Count, intervals = T2_95_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_95_2$Count, intervals = T2_95_2$Interval, partial = T)
t50(germ.counts = T2_95_2$Count, intervals = T2_95_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_95_2$Count, intervals = T2_95_2$Interval, partial = T)
WeightGermPercent(germ.counts = T2_95_2$Count, total.seeds = 20, intervals = T2_95_2$Interval, partial = T)


head(WPS73_95_2)
GermPercent(germ.counts = WPS73_95_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_95_2$Count, intervals = WPS73_95_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_95_2$Count, intervals = WPS73_95_2$Interval, partial = T)
t50(germ.counts = WPS73_95_2$Count, intervals = WPS73_95_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_95_2$Count, intervals = WPS73_95_2$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_95_2$Count, total.seeds = 20, intervals = WPS73_95_2$Interval, partial = T)



head(T1_95_5)
GermPercent(germ.counts = T1_95_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_95_5$Count, intervals = T1_95_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_95_5$Count, intervals = T1_95_5$Interval, partial = T)
t50(germ.counts = T1_95_5$Count, intervals = T1_95_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_95_5$Count, intervals = T1_95_5$Interval, partial = T)
WeightGermPercent(germ.counts = T1_95_5$Count, total.seeds = 20, intervals = T1_95_5$Interval, partial = T)

head(T2_95_5)
GermPercent(germ.counts = T2_95_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_95_5$Count, intervals = T2_95_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_95_5$Count, intervals = T2_95_5$Interval, partial = T)
t50(germ.counts = T2_95_5$Count, intervals = T2_95_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_95_5$Count, intervals = T2_95_5$Interval, partial = T)
WeightGermPercent(germ.counts = T2_95_5$Count, total.seeds = 20, intervals = T2_95_5$Interval, partial = T)


head(WPS73_95_5)
GermPercent(germ.counts = WPS73_95_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_95_5$Count, intervals = WPS73_95_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_95_5$Count, intervals = WPS73_95_5$Interval, partial = T)
t50(germ.counts = WPS73_95_5$Count, intervals = WPS73_95_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_95_5$Count, intervals = WPS73_95_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_95_5$Count, total.seeds = 20, intervals = WPS73_95_5$Interval, partial = T)



head(T1_95_10)
GermPercent(germ.counts = T1_95_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_95_10$Count, intervals = T1_95_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_95_10$Count, intervals = T1_95_10$Interval, partial = T)
t50(germ.counts = T1_95_10$Count, intervals = T1_95_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_95_10$Count, intervals = T1_95_10$Interval, partial = T)
WeightGermPercent(germ.counts = T1_95_10$Count, total.seeds = 20, intervals = T1_95_10$Interval, partial = T)

head(T2_95_10)
GermPercent(germ.counts = T2_95_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_95_10$Count, intervals = T2_95_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_95_10$Count, intervals = T2_95_10$Interval, partial = T)
t50(germ.counts = T2_95_10$Count, intervals = T2_95_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_95_10$Count, intervals = T2_95_10$Interval, partial = T)
WeightGermPercent(germ.counts = T2_95_10$Count, total.seeds = 20, intervals = T2_95_10$Interval, partial = T)


head(WPS73_95_10)
GermPercent(germ.counts = WPS73_95_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_95_10$Count, intervals = WPS73_95_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_95_10$Count, intervals = WPS73_95_10$Interval, partial = T)
t50(germ.counts = WPS73_95_10$Count, intervals = WPS73_95_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_95_10$Count, intervals = WPS73_95_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_95_10$Count, total.seeds = 20, intervals = WPS73_95_10$Interval, partial = T)








head(T1_110_30)
GermPercent(germ.counts = T1_110_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_110_30$Count, intervals = T1_110_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_110_30$Count, intervals = T1_110_30$Interval, partial = T)
t50(germ.counts = T1_110_30$Count, intervals = T1_110_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_110_30$Count, intervals = T1_110_30$Interval, partial = T)
WeightGermPercent(germ.counts = T1_110_30$Count, total.seeds = 20, intervals = T1_110_30$Interval, partial = T)

head(T2_110_30)
GermPercent(germ.counts = T2_110_30$Count, total.seeds = 17)
PeakGermPercent(germ.counts = T2_110_30$Count, intervals = T2_110_30$Interval, total.seeds = 17)
PeakGermTime(germ.counts = T2_110_30$Count, intervals = T2_110_30$Interval, partial = T)
t50(germ.counts = T2_110_30$Count, intervals = T2_110_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_110_30$Count, intervals = T2_110_30$Interval, partial = T)
WeightGermPercent(germ.counts = T2_110_30$Count, total.seeds = 17, intervals = T2_110_30$Interval, partial = T)


head(WPS73_110_30)
GermPercent(germ.counts = WPS73_110_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_110_30$Count, intervals = WPS73_110_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_110_30$Count, intervals = WPS73_110_30$Interval, partial = T)
t50(germ.counts = WPS73_110_30$Count, intervals = WPS73_110_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_110_30$Count, intervals = WPS73_110_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_110_30$Count, total.seeds = 20, intervals = WPS73_110_30$Interval, partial = T)



head(T1_110_1)
GermPercent(germ.counts = T1_110_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_110_1$Count, intervals = T1_110_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_110_1$Count, intervals = T1_110_1$Interval, partial = T)
t50(germ.counts = T1_110_1$Count, intervals = T1_110_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_110_1$Count, intervals = T1_110_1$Interval, partial = T)
WeightGermPercent(germ.counts = T1_110_1$Count, total.seeds = 20, intervals = T1_110_1$Interval, partial = T)

head(T2_110_1)
GermPercent(germ.counts = T2_110_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_110_1$Count, intervals = T2_110_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_110_1$Count, intervals = T2_110_1$Interval, partial = T)
t50(germ.counts = T2_110_1$Count, intervals = T2_110_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_110_1$Count, intervals = T2_110_1$Interval, partial = T)
WeightGermPercent(germ.counts = T2_110_1$Count, total.seeds = 20, intervals = T2_110_1$Interval, partial = T)


head(WPS73_110_1)
GermPercent(germ.counts = WPS73_110_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_110_1$Count, intervals = WPS73_110_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_110_1$Count, intervals = WPS73_110_1$Interval, partial = T)
t50(germ.counts = WPS73_110_1$Count, intervals = WPS73_110_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_110_1$Count, intervals = WPS73_110_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_110_1$Count, total.seeds = 20, intervals = WPS73_110_1$Interval, partial = T)









head(WPS89_con)
GermPercent(germ.counts = WPS89_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_con$Count, intervals = WPS89_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_con$Count, intervals = WPS89_con$Interval, partial = T)
t50(germ.counts = WPS89_con$Count, intervals = WPS89_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_con$Count, intervals = WPS89_con$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_con$Count, total.seeds = 20, intervals = WPS89_con$Interval, partial = T)

head(WPS107_con)
GermPercent(germ.counts = WPS107_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_con$Count, intervals = WPS107_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_con$Count, intervals = WPS107_con$Interval, partial = T)
t50(germ.counts = WPS107_con$Count, intervals = WPS107_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_con$Count, intervals = WPS107_con$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_con$Count, total.seeds = 20, intervals = WPS107_con$Interval, partial = T)


head(WPS154_con)
GermPercent(germ.counts = WPS154_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_con$Count, intervals = WPS154_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_con$Count, intervals = WPS154_con$Interval, partial = T)
t50(germ.counts = WPS154_con$Count, intervals = WPS154_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_con$Count, intervals = WPS154_con$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_con$Count, total.seeds = 20, intervals = WPS154_con$Interval, partial = T)




head(WPS89_80_30)
GermPercent(germ.counts = WPS89_80_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_80_30$Count, intervals = WPS89_80_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_80_30$Count, intervals = WPS89_80_30$Interval, partial = T)
t50(germ.counts = WPS89_80_30$Count, intervals = WPS89_80_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_80_30$Count, intervals = WPS89_80_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_80_30$Count, total.seeds = 20, intervals = WPS89_80_30$Interval, partial = T)

head(WPS107_80_30)
GermPercent(germ.counts = WPS107_80_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_80_30$Count, intervals = WPS107_80_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_80_30$Count, intervals = WPS107_80_30$Interval, partial = T)
t50(germ.counts = WPS107_80_30$Count, intervals = WPS107_80_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_80_30$Count, intervals = WPS107_80_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_80_30$Count, total.seeds = 20, intervals = WPS107_80_30$Interval, partial = T)


head(WPS154_80_30)
GermPercent(germ.counts = WPS154_80_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_80_30$Count, intervals = WPS154_80_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_80_30$Count, intervals = WPS154_80_30$Interval, partial = T)
t50(germ.counts = WPS154_80_30$Count, intervals = WPS154_80_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_80_30$Count, intervals = WPS154_80_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_80_30$Count, total.seeds = 20, intervals = WPS154_80_30$Interval, partial = T)



head(WPS89_80_1)
GermPercent(germ.counts = WPS89_80_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_80_1$Count, intervals = WPS89_80_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_80_1$Count, intervals = WPS89_80_1$Interval, partial = T)
t50(germ.counts = WPS89_80_1$Count, intervals = WPS89_80_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_80_1$Count, intervals = WPS89_80_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_80_1$Count, total.seeds = 20, intervals = WPS89_80_1$Interval, partial = T)

head(WPS107_80_1)
GermPercent(germ.counts = WPS107_80_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_80_1$Count, intervals = WPS107_80_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_80_1$Count, intervals = WPS107_80_1$Interval, partial = T)
t50(germ.counts = WPS107_80_1$Count, intervals = WPS107_80_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_80_1$Count, intervals = WPS107_80_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_80_1$Count, total.seeds = 20, intervals = WPS107_80_1$Interval, partial = T)


head(WPS154_80_1)
GermPercent(germ.counts = WPS154_80_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_80_1$Count, intervals = WPS154_80_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_80_1$Count, intervals = WPS154_80_1$Interval, partial = T)
t50(germ.counts = WPS154_80_1$Count, intervals = WPS154_80_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_80_1$Count, intervals = WPS154_80_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_80_1$Count, total.seeds = 20, intervals = WPS154_80_1$Interval, partial = T)



head(WPS89_80_2)
GermPercent(germ.counts = WPS89_80_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_80_2$Count, intervals = WPS89_80_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_80_2$Count, intervals = WPS89_80_2$Interval, partial = T)
t50(germ.counts = WPS89_80_2$Count, intervals = WPS89_80_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_80_2$Count, intervals = WPS89_80_2$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_80_2$Count, total.seeds = 20, intervals = WPS89_80_2$Interval, partial = T)

head(WPS107_80_2)
GermPercent(germ.counts = WPS107_80_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_80_2$Count, intervals = WPS107_80_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_80_2$Count, intervals = WPS107_80_2$Interval, partial = T)
t50(germ.counts = WPS107_80_2$Count, intervals = WPS107_80_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_80_2$Count, intervals = WPS107_80_2$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_80_2$Count, total.seeds = 20, intervals = WPS107_80_2$Interval, partial = T)


head(WPS154_80_2)
GermPercent(germ.counts = WPS154_80_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_80_2$Count, intervals = WPS154_80_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_80_2$Count, intervals = WPS154_80_2$Interval, partial = T)
t50(germ.counts = WPS154_80_2$Count, intervals = WPS154_80_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_80_2$Count, intervals = WPS154_80_2$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_80_2$Count, total.seeds = 20, intervals = WPS154_80_2$Interval, partial = T)



head(WPS89_80_5)
GermPercent(germ.counts = WPS89_80_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_80_5$Count, intervals = WPS89_80_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_80_5$Count, intervals = WPS89_80_5$Interval, partial = T)
t50(germ.counts = WPS89_80_5$Count, intervals = WPS89_80_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_80_5$Count, intervals = WPS89_80_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_80_5$Count, total.seeds = 20, intervals = WPS89_80_5$Interval, partial = T)

head(WPS107_80_5)
GermPercent(germ.counts = WPS107_80_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_80_5$Count, intervals = WPS107_80_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_80_5$Count, intervals = WPS107_80_5$Interval, partial = T)
t50(germ.counts = WPS107_80_5$Count, intervals = WPS107_80_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_80_5$Count, intervals = WPS107_80_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_80_5$Count, total.seeds = 20, intervals = WPS107_80_5$Interval, partial = T)


head(WPS154_80_5)
GermPercent(germ.counts = WPS154_80_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_80_5$Count, intervals = WPS154_80_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_80_5$Count, intervals = WPS154_80_5$Interval, partial = T)
t50(germ.counts = WPS154_80_5$Count, intervals = WPS154_80_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_80_5$Count, intervals = WPS154_80_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_80_5$Count, total.seeds = 20, intervals = WPS154_80_5$Interval, partial = T)



head(WPS89_80_10)
GermPercent(germ.counts = WPS89_80_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_80_10$Count, intervals = WPS89_80_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_80_10$Count, intervals = WPS89_80_10$Interval, partial = T)
t50(germ.counts = WPS89_80_10$Count, intervals = WPS89_80_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_80_10$Count, intervals = WPS89_80_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_80_10$Count, total.seeds = 20, intervals = WPS89_80_10$Interval, partial = T)

head(WPS107_80_10)
GermPercent(germ.counts = WPS107_80_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_80_10$Count, intervals = WPS107_80_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_80_10$Count, intervals = WPS107_80_10$Interval, partial = T)
t50(germ.counts = WPS107_80_10$Count, intervals = WPS107_80_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_80_10$Count, intervals = WPS107_80_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_80_10$Count, total.seeds = 20, intervals = WPS107_80_10$Interval, partial = T)


head(WPS154_80_10)
GermPercent(germ.counts = WPS154_80_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_80_10$Count, intervals = WPS154_80_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_80_10$Count, intervals = WPS154_80_10$Interval, partial = T)
t50(germ.counts = WPS154_80_10$Count, intervals = WPS154_80_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_80_10$Count, intervals = WPS154_80_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_80_10$Count, total.seeds = 20, intervals = WPS154_80_10$Interval, partial = T)





head(WPS89_95_30)
GermPercent(germ.counts = WPS89_95_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_95_30$Count, intervals = WPS89_95_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_95_30$Count, intervals = WPS89_95_30$Interval, partial = T)
t50(germ.counts = WPS89_95_30$Count, intervals = WPS89_95_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_95_30$Count, intervals = WPS89_95_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_95_30$Count, total.seeds = 20, intervals = WPS89_95_30$Interval, partial = T)

head(WPS107_95_30)
GermPercent(germ.counts = WPS107_95_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_95_30$Count, intervals = WPS107_95_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_95_30$Count, intervals = WPS107_95_30$Interval, partial = T)
t50(germ.counts = WPS107_95_30$Count, intervals = WPS107_95_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_95_30$Count, intervals = WPS107_95_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_95_30$Count, total.seeds = 20, intervals = WPS107_95_30$Interval, partial = T)


head(WPS154_95_30)
GermPercent(germ.counts = WPS154_95_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_95_30$Count, intervals = WPS154_95_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_95_30$Count, intervals = WPS154_95_30$Interval, partial = T)
t50(germ.counts = WPS154_95_30$Count, intervals = WPS154_95_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_95_30$Count, intervals = WPS154_95_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_95_30$Count, total.seeds = 20, intervals = WPS154_95_30$Interval, partial = T)



head(WPS89_95_1)
GermPercent(germ.counts = WPS89_95_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_95_1$Count, intervals = WPS89_95_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_95_1$Count, intervals = WPS89_95_1$Interval, partial = T)
t50(germ.counts = WPS89_95_1$Count, intervals = WPS89_95_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_95_1$Count, intervals = WPS89_95_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_95_1$Count, total.seeds = 20, intervals = WPS89_95_1$Interval, partial = T)

head(WPS107_95_1)
GermPercent(germ.counts = WPS107_95_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_95_1$Count, intervals = WPS107_95_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_95_1$Count, intervals = WPS107_95_1$Interval, partial = T)
t50(germ.counts = WPS107_95_1$Count, intervals = WPS107_95_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_95_1$Count, intervals = WPS107_95_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_95_1$Count, total.seeds = 20, intervals = WPS107_95_1$Interval, partial = T)


head(WPS154_95_1)
GermPercent(germ.counts = WPS154_95_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_95_1$Count, intervals = WPS154_95_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_95_1$Count, intervals = WPS154_95_1$Interval, partial = T)
t50(germ.counts = WPS154_95_1$Count, intervals = WPS154_95_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_95_1$Count, intervals = WPS154_95_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_95_1$Count, total.seeds = 20, intervals = WPS154_95_1$Interval, partial = T)



head(WPS89_95_2)
GermPercent(germ.counts = WPS89_95_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_95_2$Count, intervals = WPS89_95_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_95_2$Count, intervals = WPS89_95_2$Interval, partial = T)
t50(germ.counts = WPS89_95_2$Count, intervals = WPS89_95_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_95_2$Count, intervals = WPS89_95_2$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_95_2$Count, total.seeds = 20, intervals = WPS89_95_2$Interval, partial = T)

head(WPS107_95_2)
GermPercent(germ.counts = WPS107_95_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_95_2$Count, intervals = WPS107_95_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_95_2$Count, intervals = WPS107_95_2$Interval, partial = T)
t50(germ.counts = WPS107_95_2$Count, intervals = WPS107_95_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_95_2$Count, intervals = WPS107_95_2$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_95_2$Count, total.seeds = 20, intervals = WPS107_95_2$Interval, partial = T)


head(WPS154_95_2)
GermPercent(germ.counts = WPS154_95_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_95_2$Count, intervals = WPS154_95_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_95_2$Count, intervals = WPS154_95_2$Interval, partial = T)
t50(germ.counts = WPS154_95_2$Count, intervals = WPS154_95_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_95_2$Count, intervals = WPS154_95_2$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_95_2$Count, total.seeds = 20, intervals = WPS154_95_2$Interval, partial = T)



head(WPS89_95_5)
GermPercent(germ.counts = WPS89_95_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_95_5$Count, intervals = WPS89_95_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_95_5$Count, intervals = WPS89_95_5$Interval, partial = T)
t50(germ.counts = WPS89_95_5$Count, intervals = WPS89_95_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_95_5$Count, intervals = WPS89_95_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_95_5$Count, total.seeds = 20, intervals = WPS89_95_5$Interval, partial = T)

head(WPS107_95_5)
GermPercent(germ.counts = WPS107_95_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_95_5$Count, intervals = WPS107_95_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_95_5$Count, intervals = WPS107_95_5$Interval, partial = T)
t50(germ.counts = WPS107_95_5$Count, intervals = WPS107_95_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_95_5$Count, intervals = WPS107_95_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_95_5$Count, total.seeds = 20, intervals = WPS107_95_5$Interval, partial = T)


head(WPS154_95_5)
GermPercent(germ.counts = WPS154_95_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_95_5$Count, intervals = WPS154_95_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_95_5$Count, intervals = WPS154_95_5$Interval, partial = T)
t50(germ.counts = WPS154_95_5$Count, intervals = WPS154_95_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_95_5$Count, intervals = WPS154_95_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_95_5$Count, total.seeds = 20, intervals = WPS154_95_5$Interval, partial = T)



head(WPS89_95_10)
GermPercent(germ.counts = WPS89_95_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_95_10$Count, intervals = WPS89_95_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_95_10$Count, intervals = WPS89_95_10$Interval, partial = T)
t50(germ.counts = WPS89_95_10$Count, intervals = WPS89_95_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_95_10$Count, intervals = WPS89_95_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_95_10$Count, total.seeds = 20, intervals = WPS89_95_10$Interval, partial = T)

head(WPS107_95_10)
GermPercent(germ.counts = WPS107_95_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_95_10$Count, intervals = WPS107_95_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_95_10$Count, intervals = WPS107_95_10$Interval, partial = T)
t50(germ.counts = WPS107_95_10$Count, intervals = WPS107_95_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_95_10$Count, intervals = WPS107_95_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_95_10$Count, total.seeds = 20, intervals = WPS107_95_10$Interval, partial = T)


head(WPS154_95_10)
GermPercent(germ.counts = WPS154_95_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_95_10$Count, intervals = WPS154_95_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_95_10$Count, intervals = WPS154_95_10$Interval, partial = T)
t50(germ.counts = WPS154_95_10$Count, intervals = WPS154_95_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_95_10$Count, intervals = WPS154_95_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_95_10$Count, total.seeds = 20, intervals = WPS154_95_10$Interval, partial = T)








head(WPS89_110_30)
GermPercent(germ.counts = WPS89_110_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_110_30$Count, intervals = WPS89_110_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_110_30$Count, intervals = WPS89_110_30$Interval, partial = T)
t50(germ.counts = WPS89_110_30$Count, intervals = WPS89_110_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_110_30$Count, intervals = WPS89_110_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_110_30$Count, total.seeds = 20, intervals = WPS89_110_30$Interval, partial = T)

head(WPS107_110_30)
GermPercent(germ.counts = WPS107_110_30$Count, total.seeds = 19)
PeakGermPercent(germ.counts = WPS107_110_30$Count, intervals = WPS107_110_30$Interval, total.seeds = 19)
PeakGermTime(germ.counts = WPS107_110_30$Count, intervals = WPS107_110_30$Interval, partial = T)
t50(germ.counts = WPS107_110_30$Count, intervals = WPS107_110_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_110_30$Count, intervals = WPS107_110_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_110_30$Count, total.seeds = 19, intervals = WPS107_110_30$Interval, partial = T)


head(WPS154_110_30)
GermPercent(germ.counts = WPS154_110_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_110_30$Count, intervals = WPS154_110_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_110_30$Count, intervals = WPS154_110_30$Interval, partial = T)
t50(germ.counts = WPS154_110_30$Count, intervals = WPS154_110_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_110_30$Count, intervals = WPS154_110_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_110_30$Count, total.seeds = 20, intervals = WPS154_110_30$Interval, partial = T)



head(WPS89_110_1)
GermPercent(germ.counts = WPS89_110_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_110_1$Count, intervals = WPS89_110_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_110_1$Count, intervals = WPS89_110_1$Interval, partial = T)
t50(germ.counts = WPS89_110_1$Count, intervals = WPS89_110_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_110_1$Count, intervals = WPS89_110_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_110_1$Count, total.seeds = 20, intervals = WPS89_110_1$Interval, partial = T)

head(WPS107_110_1)
GermPercent(germ.counts = WPS107_110_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_110_1$Count, intervals = WPS107_110_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_110_1$Count, intervals = WPS107_110_1$Interval, partial = T)
t50(germ.counts = WPS107_110_1$Count, intervals = WPS107_110_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_110_1$Count, intervals = WPS107_110_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_110_1$Count, total.seeds = 20, intervals = WPS107_110_1$Interval, partial = T)


head(WPS154_110_1)
GermPercent(germ.counts = WPS154_110_1$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_110_1$Count, intervals = WPS154_110_1$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_110_1$Count, intervals = WPS154_110_1$Interval, partial = T)
t50(germ.counts = WPS154_110_1$Count, intervals = WPS154_110_1$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_110_1$Count, intervals = WPS154_110_1$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_110_1$Count, total.seeds = 20, intervals = WPS154_110_1$Interval, partial = T)



head(WPS89_110_2)
GermPercent(germ.counts = WPS89_110_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_110_2$Count, intervals = WPS89_110_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_110_2$Count, intervals = WPS89_110_2$Interval, partial = T)
t50(germ.counts = WPS89_110_2$Count, intervals = WPS89_110_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_110_2$Count, intervals = WPS89_110_2$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_110_2$Count, total.seeds = 20, intervals = WPS89_110_2$Interval, partial = T)

head(WPS107_110_2)
GermPercent(germ.counts = WPS107_110_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_110_2$Count, intervals = WPS107_110_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_110_2$Count, intervals = WPS107_110_2$Interval, partial = T)
t50(germ.counts = WPS107_110_2$Count, intervals = WPS107_110_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_110_2$Count, intervals = WPS107_110_2$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_110_2$Count, total.seeds = 20, intervals = WPS107_110_2$Interval, partial = T)


head(WPS154_110_2)
GermPercent(germ.counts = WPS154_110_2$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_110_2$Count, intervals = WPS154_110_2$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_110_2$Count, intervals = WPS154_110_2$Interval, partial = T)
t50(germ.counts = WPS154_110_2$Count, intervals = WPS154_110_2$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_110_2$Count, intervals = WPS154_110_2$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_110_2$Count, total.seeds = 20, intervals = WPS154_110_2$Interval, partial = T)



head(WPS89_110_5)
GermPercent(germ.counts = WPS89_110_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_110_5$Count, intervals = WPS89_110_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_110_5$Count, intervals = WPS89_110_5$Interval, partial = T)
t50(germ.counts = WPS89_110_5$Count, intervals = WPS89_110_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_110_5$Count, intervals = WPS89_110_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_110_5$Count, total.seeds = 20, intervals = WPS89_110_5$Interval, partial = T)

head(WPS107_110_5)
GermPercent(germ.counts = WPS107_110_5$Count, total.seeds = 19)
PeakGermPercent(germ.counts = WPS107_110_5$Count, intervals = WPS107_110_5$Interval, total.seeds = 19)
PeakGermTime(germ.counts = WPS107_110_5$Count, intervals = WPS107_110_5$Interval, partial = T)
t50(germ.counts = WPS107_110_5$Count, intervals = WPS107_110_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_110_5$Count, intervals = WPS107_110_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_110_5$Count, total.seeds = 19, intervals = WPS107_110_5$Interval, partial = T)


head(WPS154_110_5)
GermPercent(germ.counts = WPS154_110_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_110_5$Count, intervals = WPS154_110_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_110_5$Count, intervals = WPS154_110_5$Interval, partial = T)
t50(germ.counts = WPS154_110_5$Count, intervals = WPS154_110_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_110_5$Count, intervals = WPS154_110_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_110_5$Count, total.seeds = 20, intervals = WPS154_110_5$Interval, partial = T)



head(WPS89_110_10)
GermPercent(germ.counts = WPS89_110_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_110_10$Count, intervals = WPS89_110_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_110_10$Count, intervals = WPS89_110_10$Interval, partial = T)
t50(germ.counts = WPS89_110_10$Count, intervals = WPS89_110_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_110_10$Count, intervals = WPS89_110_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_110_10$Count, total.seeds = 20, intervals = WPS89_110_10$Interval, partial = T)

head(WPS107_110_10)
GermPercent(germ.counts = WPS107_110_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_110_10$Count, intervals = WPS107_110_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_110_10$Count, intervals = WPS107_110_10$Interval, partial = T)
t50(germ.counts = WPS107_110_10$Count, intervals = WPS107_110_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_110_10$Count, intervals = WPS107_110_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_110_10$Count, total.seeds = 20, intervals = WPS107_110_10$Interval, partial = T)


head(WPS154_110_10)
GermPercent(germ.counts = WPS154_110_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_110_10$Count, intervals = WPS154_110_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_110_10$Count, intervals = WPS154_110_10$Interval, partial = T)
t50(germ.counts = WPS154_110_10$Count, intervals = WPS154_110_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_110_10$Count, intervals = WPS154_110_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_110_10$Count, total.seeds = 20, intervals = WPS154_110_10$Interval, partial = T)


# 2. Smoke treatments -----
smoke <- read.table('./00_Data/Optimisation_tests/smoke_tests.txt', header = T)
head(smoke); dim(smoke)

smoke$Time_spread <- smoke$Time_to_finish - smoke$Time_to_germ
smoke.seed <- smoke[, which(colnames(smoke) == "X2"): which(colnames(smoke) == "X35")]
smoke_cum.sum <- data.frame(t(apply(smoke.seed, 1, FUN = function(x) cumsum(unlist(x)))))
smoke_cum.prop <- smoke_cum.sum/smoke$Total_seeds_adjusted
smoke_cum.prop2 <- cbind(smoke[, c(1:4, 22)], smoke_cum.prop)
head(smoke_cum.prop2); dim(smoke_cum.prop2)
str(smoke_cum.prop2)


# Split the data by species for plotting 
lit_smoke <- subset(smoke_cum.prop2, smoke_cum.prop2$Species == "littoralis")
rownames(lit_smoke) <- 1:nrow(lit_smoke)
head(lit_smoke); dim(lit_smoke)

tor_smoke <- subset(smoke_cum.prop2, smoke_cum.prop2$Species == "torulosa")
rownames(tor_smoke) <- 1:nrow(tor_smoke)
head(tor_smoke); dim(tor_smoke)

# Extract just the proportion data
prop.smoke_lit <- lit_smoke[, which(colnames(lit_smoke) == "X2"): which(colnames(lit_smoke) == "X35")]
prop.smoke_tor <- tor_smoke[, which(colnames(tor_smoke) == "X2"): which(colnames(tor_smoke) == "X35")]


# In this case we want to split the plot by smoke type, it is too hard to visualise the durations and different types of smoke all on the one plot because of the number of lines. 
# So we want to subset the data further by smoke type for each species 
lit_allocas_smoke <- subset(lit_smoke, lit_smoke$Smoke_type == 'Allocas_leaf' | lit_smoke$Smoke_type == "Control")
tor_allocas_smoke <- subset(tor_smoke, tor_smoke$Smoke_type == "Allocas_leaf" | lit_smoke$Smoke_type == "Control")
lit_pine_smoke <- subset(lit_smoke, lit_smoke$Smoke_type == 'Pine_sawdust'| tor_smoke$Smoke_type == 'Control')
tor_pine_smoke <- subset(tor_smoke, tor_smoke$Smoke_type == "Pine_sawdust" | tor_smoke$Smoke_type == 'Control')



# Set colour palette to be used for plotting
# We want the shortest durations to be blue and light, with increasing duration meaning darker colours to red.
library(RColorBrewer)
brewer.pal(4, 'RdYlBu')
pal2 <- c('black',"#ABD9E9","#2C7BB6", "#FDAE61","#D7191C")
smoke.col <- data.frame(col = pal2, dur = c('0','5','10','20','30'))


# Create the plots for smoke
dev.new(height = 12, width = 14, dpi = 80, pointsize = 16, noRStudioGD = T)
par(mar = c(6,4,2,1), mgp = c(3,1,0), mfrow = c(2,2), oma = c(0,0,3,8))

plot(1:ncol(prop.smoke_lit), 1:ncol(prop.smoke_lit), ylim = c(0,1), type = 'n', las = 1, ylab = "Proportion germination", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.smoke_lit), labels = c(2,4,9,11,14,16,18,21,23,25,28,30,32,35), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(main = "Allocasuarina leaf")
for(i in 1:nrow(lit_allocas_smoke)){
  dat.thisrun <- lit_allocas_smoke[i, 6:ncol(lit_allocas_smoke)]
  duration.thisrun <- lit_allocas_smoke$Duration[i]
  col.thisrun <- smoke.col$col[smoke.col$dur == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}

plot(1:ncol(prop.smoke_lit), 1:ncol(prop.smoke_lit), ylim = c(0,1), type = 'n', las = 1, ylab = "", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.smoke_lit), labels = c(2,4,9,11,14,16,18,21,23,25,28,30,32,35), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(main = "Pine sawdust")
for(i in 1:nrow(lit_pine_smoke)){
  dat.thisrun <- lit_pine_smoke[i, 6:ncol(lit_pine_smoke)]
  duration.thisrun <- lit_pine_smoke$Duration[i]
  col.thisrun <- smoke.col$col[smoke.col$dur == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}

mtext(expression('(a) '~italic(Allocasuarina~littoralis)), side = 3, outer = T, font = 2)

mtext(expression('(b) '~italic(Allocasuarina~torulosa)), side = 3, line = -26, outer = T, font = 2)

plot(1:ncol(prop.smoke_tor), 1:ncol(prop.smoke_tor), ylim = c(0,1), type = 'n', las = 1, ylab = "Proportion germination", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.smoke_tor), labels = c(2,4,9,11,14,16,18,21,23,25,28,30,32,35), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(xlab = "Days", mgp = c(2.2, 1, 0))
title(main = "Allocasuarina leaf")
for(i in 1:nrow(tor_allocas_smoke)){
  dat.thisrun <- tor_allocas_smoke[i, 6:ncol(tor_allocas_smoke)]
  duration.thisrun <- tor_allocas_smoke$Duration[i]
  col.thisrun <- smoke.col$col[smoke.col$dur == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}


plot(1:ncol(prop.smoke_tor), 1:ncol(prop.smoke_tor), ylim = c(0,1), type = 'n', las = 1, ylab = "", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.smoke_tor), labels = c(2,4,9,11,14,16,18,21,23,25,28,30,32,35), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(xlab = "Days", mgp = c(2.2, 1, 0))
title(main = "Pine sawdust")
for(i in 1:nrow(tor_pine_smoke)){
  dat.thisrun <- tor_pine_smoke[i, 6:ncol(tor_pine_smoke)]
  duration.thisrun <- tor_pine_smoke$Duration[i]
  col.thisrun <- smoke.col$col[smoke.col$dur == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}

# Legend being produced incorrectly
par(xpd = NA)
legend(x = 15, y = 2.5, legend = c("Control", "5 min", '10 min', '20 min', '30 min'), col = smoke.col$col, lwd = 4, cex = 0.95, text.width = 0.2, title = 'Smoke exposure', bty = "n")
par(xpd = F)


# Pine sawdust seems to have a more consistent effect than allocasuarina leaf litter. 10 minutes looks best for smoke exposure. Otherwise if we do a second smoke I would choose 5 or 20 minutes but should look at some of the other comparison metrics to decide further


# Plot pine sawdust alone for each species
dev.new(height = 6, width = 12, dpi = 80, pointsize = 16, noRStudioGD = T)
par(mar = c(4,4,2,1), mgp = c(3,1,0), mfrow = c(1,2), oma = c(0,0,0,8))

plot(1:ncol(prop.smoke_lit), 1:ncol(prop.smoke_lit), ylim = c(0,1), type = 'n', las = 1, ylab = "Proportion germination", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.smoke_lit), labels = c(2,4,9,11,14,16,18,21,23,25,28,30,32,35), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(xlab = "Days", mgp = c(2.2, 0.7,0))
title(main = "(a) A. littoralis")
for(i in 1:nrow(lit_pine_smoke)){
  dat.thisrun <- lit_pine_smoke[i, 6:ncol(lit_pine_smoke)]
  duration.thisrun <- lit_pine_smoke$Duration[i]
  col.thisrun <- smoke.col$col[smoke.col$dur == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}

plot(1:ncol(prop.smoke_tor), 1:ncol(prop.smoke_tor), ylim = c(0,1), type = 'n', las = 1, ylab = "Proportion germination", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.smoke_tor), labels = c(2,4,9,11,14,16,18,21,23,25,28,30,32,35), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(xlab = "Days", mgp = c(2.2, 1, 0))
title(main = "(b) A. torulosa")
for(i in 1:nrow(tor_pine_smoke)){
  dat.thisrun <- tor_pine_smoke[i, 6:ncol(tor_pine_smoke)]
  duration.thisrun <- tor_pine_smoke$Duration[i]
  col.thisrun <- smoke.col$col[smoke.col$dur == duration.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}


par(xpd = NA)
legend(x = 15, y = 1, legend = c("Control", "5 min", '10 min', '20 min', '30 min'), col = smoke.col$col, lwd = 4, cex = 0.95, text.width = 0.2, title = 'Smoke exposure', bty = "n")
par(xpd = F)



# Thinking 10 minutes and/or 5 minutes smoke exposure would be better. 



# Take a look at some summary statistics
library(germinationmetrics)
smoke_long <- read.table('./00_Data/Optimisation_tests/smoke_long.txt', header = T)
head(smoke_long)

# To be able to work with the data we need to subset by seed lot, smoke type and duration
T1_con <- subset(smoke_long, smoke_long$Seed_lot == "T1" & smoke_long$Smoke_type == "Control")
T2_con <- subset(smoke_long, smoke_long$Seed_lot == "T2" & smoke_long$Smoke_type == "Control")
WPS73_con <- subset(smoke_long, smoke_long$Seed_lot == "WPS73" & smoke_long$Smoke_type == "Control")


T1_allo_5 <- subset(smoke_long, smoke_long$Seed_lot == "T1" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='5')
T1_allo_10 <- subset(smoke_long, smoke_long$Seed_lot == "T1" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='10')
T1_allo_20 <- subset(smoke_long, smoke_long$Seed_lot == "T1" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='20')
T1_allo_30 <- subset(smoke_long, smoke_long$Seed_lot == "T1" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='30')


T2_allo_5 <- subset(smoke_long, smoke_long$Seed_lot == "T2" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='5')
T2_allo_10 <- subset(smoke_long, smoke_long$Seed_lot == "T2" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='10')
T2_allo_20 <- subset(smoke_long, smoke_long$Seed_lot == "T2" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='20')
T2_allo_30 <- subset(smoke_long, smoke_long$Seed_lot == "T2" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='30')


WPS73_allo_5 <- subset(smoke_long, smoke_long$Seed_lot == "WPS73" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='5')
WPS73_allo_10 <- subset(smoke_long, smoke_long$Seed_lot == "WPS73" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='10')
WPS73_allo_20 <- subset(smoke_long, smoke_long$Seed_lot == "WPS73" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='20')
WPS73_allo_30 <- subset(smoke_long, smoke_long$Seed_lot == "WPS73" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='30')



T1_pine_5 <- subset(smoke_long, smoke_long$Seed_lot == "T1" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='5')
T1_pine_10 <- subset(smoke_long, smoke_long$Seed_lot == "T1" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='10')
T1_pine_20 <- subset(smoke_long, smoke_long$Seed_lot == "T1" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='20')
T1_pine_30 <- subset(smoke_long, smoke_long$Seed_lot == "T1" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='30')

T2_pine_5 <- subset(smoke_long, smoke_long$Seed_lot == "T2" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='5')
T2_pine_10 <- subset(smoke_long, smoke_long$Seed_lot == "T2" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='10')
T2_pine_20 <- subset(smoke_long, smoke_long$Seed_lot == "T2" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='20')
T2_pine_30 <- subset(smoke_long, smoke_long$Seed_lot == "T2" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='30')

WPS73_pine_5 <- subset(smoke_long, smoke_long$Seed_lot == "WPS73" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='5')
WPS73_pine_10 <- subset(smoke_long, smoke_long$Seed_lot == "WPS73" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='10')
WPS73_pine_20 <- subset(smoke_long, smoke_long$Seed_lot == "WPS73" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='20')
WPS73_pine_30 <- subset(smoke_long, smoke_long$Seed_lot == "WPS73" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='30')



WPS89_con <- subset(smoke_long, smoke_long$Seed_lot == "WPS89" & smoke_long$Smoke_type == "Control")
WPS107_con <- subset(smoke_long, smoke_long$Seed_lot == "WPS107" & smoke_long$Smoke_type == "Control")
WPS154_con <- subset(smoke_long, smoke_long$Seed_lot == "WPS154" & smoke_long$Smoke_type == "Control")


WPS89_allo_5 <- subset(smoke_long, smoke_long$Seed_lot == "WPS89" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='5')
WPS89_allo_10 <- subset(smoke_long, smoke_long$Seed_lot == "WPS89" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='10')
WPS89_allo_20 <- subset(smoke_long, smoke_long$Seed_lot == "WPS89" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='20')
WPS89_allo_30 <- subset(smoke_long, smoke_long$Seed_lot == "WPS89" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='30')


WPS107_allo_5 <- subset(smoke_long, smoke_long$Seed_lot == "WPS107" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='5')
WPS107_allo_10 <- subset(smoke_long, smoke_long$Seed_lot == "WPS107" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='10')
WPS107_allo_20 <- subset(smoke_long, smoke_long$Seed_lot == "WPS107" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='20')
WPS107_allo_30 <- subset(smoke_long, smoke_long$Seed_lot == "WPS107" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='30')


WPS154_allo_5 <- subset(smoke_long, smoke_long$Seed_lot == "WPS154" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='5')
WPS154_allo_10 <- subset(smoke_long, smoke_long$Seed_lot == "WPS154" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='10')
WPS154_allo_20 <- subset(smoke_long, smoke_long$Seed_lot == "WPS154" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='20')
WPS154_allo_30 <- subset(smoke_long, smoke_long$Seed_lot == "WPS154" & smoke_long$Smoke_type == "Allocas" & smoke_long$Duration =='30')



WPS89_pine_5 <- subset(smoke_long, smoke_long$Seed_lot == "WPS89" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='5')
WPS89_pine_10 <- subset(smoke_long, smoke_long$Seed_lot == "WPS89" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='10')
WPS89_pine_20 <- subset(smoke_long, smoke_long$Seed_lot == "WPS89" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='20')
WPS89_pine_30 <- subset(smoke_long, smoke_long$Seed_lot == "WPS89" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='30')

WPS107_pine_5 <- subset(smoke_long, smoke_long$Seed_lot == "WPS107" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='5')
WPS107_pine_10 <- subset(smoke_long, smoke_long$Seed_lot == "WPS107" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='10')
WPS107_pine_20 <- subset(smoke_long, smoke_long$Seed_lot == "WPS107" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='20')
WPS107_pine_30 <- subset(smoke_long, smoke_long$Seed_lot == "WPS107" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='30')

WPS154_pine_5 <- subset(smoke_long, smoke_long$Seed_lot == "WPS154" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='5')
WPS154_pine_10 <- subset(smoke_long, smoke_long$Seed_lot == "WPS154" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='10')
WPS154_pine_20 <- subset(smoke_long, smoke_long$Seed_lot == "WPS154" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='20')
WPS154_pine_30 <- subset(smoke_long, smoke_long$Seed_lot == "WPS154" & smoke_long$Smoke_type == "Pine" & smoke_long$Duration =='30')


# We are interested in the following metrics: 
# Germ Percent, PeakGermPercent, PeakGermTime, t50, MeanGermTime, VarGermTime, MeanGermRate, VarGermRate, WeightGermPercent, MeanGermNumber

# Peak germination percent is the mean daily germination of the most vigorous component of the seed lot

# Controls

head(T1_con)
GermPercent(germ.counts = T1_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_con$Count, intervals = T1_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_con$Count, intervals = T1_con$Interval, partial = T)
t50(germ.counts = T1_con$Count, intervals = T1_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_con$Count, intervals = T1_con$Interval, partial = T)
WeightGermPercent(germ.counts = T1_con$Count, total.seeds = 20, intervals = T1_con$Interval, partial = T)

head(T2_con)
GermPercent(germ.counts = T2_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_con$Count, intervals = T2_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_con$Count, intervals = T2_con$Interval, partial = T)
t50(germ.counts = T2_con$Count, intervals = T2_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_con$Count, intervals = T2_con$Interval, partial = T)
WeightGermPercent(germ.counts = T2_con$Count, total.seeds = 20, intervals = T2_con$Interval, partial = T)


head(WPS73_con)
GermPercent(germ.counts = WPS73_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_con$Count, intervals = WPS73_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_con$Count, intervals = WPS73_con$Interval, partial = T)
t50(germ.counts = WPS73_con$Count, intervals = WPS73_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_con$Count, intervals = WPS73_con$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_con$Count, total.seeds = 20, intervals = WPS73_con$Interval, partial = T)



head(WPS89_con)
GermPercent(germ.counts = WPS89_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_con$Count, intervals = WPS89_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_con$Count, intervals = WPS89_con$Interval, partial = T)
t50(germ.counts = WPS89_con$Count, intervals = WPS89_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_con$Count, intervals = WPS89_con$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_con$Count, total.seeds = 20, intervals = WPS89_con$Interval, partial = T)

head(WPS107_con)
GermPercent(germ.counts = WPS107_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_con$Count, intervals = WPS107_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_con$Count, intervals = WPS107_con$Interval, partial = T)
t50(germ.counts = WPS107_con$Count, intervals = WPS107_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_con$Count, intervals = WPS107_con$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_con$Count, total.seeds = 20, intervals = WPS107_con$Interval, partial = T)


head(WPS154_con)
GermPercent(germ.counts = WPS154_con$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_con$Count, intervals = WPS154_con$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_con$Count, intervals = WPS154_con$Interval, partial = T)
t50(germ.counts = WPS154_con$Count, intervals = WPS154_con$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_con$Count, intervals = WPS154_con$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_con$Count, total.seeds = 20, intervals = WPS154_con$Interval, partial = T)


# Allocas smoke 5 min

head(T1_allo_5)
GermPercent(germ.counts = T1_allo_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_allo_5$Count, intervals = T1_allo_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_allo_5$Count, intervals = T1_allo_5$Interval, partial = T)
t50(germ.counts = T1_allo_5$Count, intervals = T1_allo_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_allo_5$Count, intervals = T1_allo_5$Interval, partial = T)
WeightGermPercent(germ.counts = T1_allo_5$Count, total.seeds = 20, intervals = T1_allo_5$Interval, partial = T)

head(T2_allo_5)
GermPercent(germ.counts = T2_allo_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_allo_5$Count, intervals = T2_allo_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_allo_5$Count, intervals = T2_allo_5$Interval, partial = T)
t50(germ.counts = T2_allo_5$Count, intervals = T2_allo_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_allo_5$Count, intervals = T2_allo_5$Interval, partial = T)
WeightGermPercent(germ.counts = T2_allo_5$Count, total.seeds = 20, intervals = T2_allo_5$Interval, partial = T)


head(WPS73_allo_5)
GermPercent(germ.counts = WPS73_allo_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_allo_5$Count, intervals = WPS73_allo_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_allo_5$Count, intervals = WPS73_allo_5$Interval, partial = T)
t50(germ.counts = WPS73_allo_5$Count, intervals = WPS73_allo_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_allo_5$Count, intervals = WPS73_allo_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_allo_5$Count, total.seeds = 20, intervals = WPS73_allo_5$Interval, partial = T)



head(WPS89_allo_5)
GermPercent(germ.counts = WPS89_allo_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_allo_5$Count, intervals = WPS89_allo_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_allo_5$Count, intervals = WPS89_allo_5$Interval, partial = T)
t50(germ.counts = WPS89_allo_5$Count, intervals = WPS89_allo_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_allo_5$Count, intervals = WPS89_allo_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_allo_5$Count, total.seeds = 20, intervals = WPS89_allo_5$Interval, partial = T)

head(WPS107_allo_5)
GermPercent(germ.counts = WPS107_allo_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_allo_5$Count, intervals = WPS107_allo_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_allo_5$Count, intervals = WPS107_allo_5$Interval, partial = T)
t50(germ.counts = WPS107_allo_5$Count, intervals = WPS107_allo_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_allo_5$Count, intervals = WPS107_allo_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_allo_5$Count, total.seeds = 20, intervals = WPS107_allo_5$Interval, partial = T)


head(WPS154_allo_5)
GermPercent(germ.counts = WPS154_allo_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_allo_5$Count, intervals = WPS154_allo_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_allo_5$Count, intervals = WPS154_allo_5$Interval, partial = T)
t50(germ.counts = WPS154_allo_5$Count, intervals = WPS154_allo_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_allo_5$Count, intervals = WPS154_allo_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_allo_5$Count, total.seeds = 20, intervals = WPS154_allo_5$Interval, partial = T)


# allo smoke 10 min


head(T1_allo_10)
GermPercent(germ.counts = T1_allo_10$Count, total.seeds = 18)
PeakGermPercent(germ.counts = T1_allo_10$Count, intervals = T1_allo_10$Interval, total.seeds = 18)
PeakGermTime(germ.counts = T1_allo_10$Count, intervals = T1_allo_10$Interval, partial = T)
t50(germ.counts = T1_allo_10$Count, intervals = T1_allo_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_allo_10$Count, intervals = T1_allo_10$Interval, partial = T)
WeightGermPercent(germ.counts = T1_allo_10$Count, total.seeds = 18, intervals = T1_allo_10$Interval, partial = T)

head(T2_allo_10)
GermPercent(germ.counts = T2_allo_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_allo_10$Count, intervals = T2_allo_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_allo_10$Count, intervals = T2_allo_10$Interval, partial = T)
t50(germ.counts = T2_allo_10$Count, intervals = T2_allo_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_allo_10$Count, intervals = T2_allo_10$Interval, partial = T)
WeightGermPercent(germ.counts = T2_allo_10$Count, total.seeds = 20, intervals = T2_allo_10$Interval, partial = T)


head(WPS73_allo_10)
GermPercent(germ.counts = WPS73_allo_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_allo_10$Count, intervals = WPS73_allo_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_allo_10$Count, intervals = WPS73_allo_10$Interval, partial = T)
t50(germ.counts = WPS73_allo_10$Count, intervals = WPS73_allo_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_allo_10$Count, intervals = WPS73_allo_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_allo_10$Count, total.seeds = 20, intervals = WPS73_allo_10$Interval, partial = T)



head(WPS89_allo_10)
GermPercent(germ.counts = WPS89_allo_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_allo_10$Count, intervals = WPS89_allo_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_allo_10$Count, intervals = WPS89_allo_10$Interval, partial = T)
t50(germ.counts = WPS89_allo_10$Count, intervals = WPS89_allo_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_allo_10$Count, intervals = WPS89_allo_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_allo_10$Count, total.seeds = 20, intervals = WPS89_allo_10$Interval, partial = T)

head(WPS107_allo_10)
GermPercent(germ.counts = WPS107_allo_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_allo_10$Count, intervals = WPS107_allo_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_allo_10$Count, intervals = WPS107_allo_10$Interval, partial = T)
t50(germ.counts = WPS107_allo_10$Count, intervals = WPS107_allo_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_allo_10$Count, intervals = WPS107_allo_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_allo_10$Count, total.seeds = 20, intervals = WPS107_allo_10$Interval, partial = T)


head(WPS154_allo_10)
GermPercent(germ.counts = WPS154_allo_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_allo_10$Count, intervals = WPS154_allo_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_allo_10$Count, intervals = WPS154_allo_10$Interval, partial = T)
t50(germ.counts = WPS154_allo_10$Count, intervals = WPS154_allo_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_allo_10$Count, intervals = WPS154_allo_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_allo_10$Count, total.seeds = 20, intervals = WPS154_allo_10$Interval, partial = T)



# allo smoke 20 min

head(T1_allo_20)
GermPercent(germ.counts = T1_allo_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_allo_20$Count, intervals = T1_allo_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_allo_20$Count, intervals = T1_allo_20$Interval, partial = T)
t50(germ.counts = T1_allo_20$Count, intervals = T1_allo_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_allo_20$Count, intervals = T1_allo_20$Interval, partial = T)
WeightGermPercent(germ.counts = T1_allo_20$Count, total.seeds = 20, intervals = T1_allo_20$Interval, partial = T)

head(T2_allo_20)
GermPercent(germ.counts = T2_allo_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_allo_20$Count, intervals = T2_allo_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_allo_20$Count, intervals = T2_allo_20$Interval, partial = T)
t50(germ.counts = T2_allo_20$Count, intervals = T2_allo_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_allo_20$Count, intervals = T2_allo_20$Interval, partial = T)
WeightGermPercent(germ.counts = T2_allo_20$Count, total.seeds = 20, intervals = T2_allo_20$Interval, partial = T)


head(WPS73_allo_20)
GermPercent(germ.counts = WPS73_allo_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_allo_20$Count, intervals = WPS73_allo_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_allo_20$Count, intervals = WPS73_allo_20$Interval, partial = T)
t50(germ.counts = WPS73_allo_20$Count, intervals = WPS73_allo_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_allo_20$Count, intervals = WPS73_allo_20$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_allo_20$Count, total.seeds = 20, intervals = WPS73_allo_20$Interval, partial = T)


head(WPS89_allo_20)
GermPercent(germ.counts = WPS89_allo_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_allo_20$Count, intervals = WPS89_allo_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_allo_20$Count, intervals = WPS89_allo_20$Interval, partial = T)
t50(germ.counts = WPS89_allo_20$Count, intervals = WPS89_allo_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_allo_20$Count, intervals = WPS89_allo_20$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_allo_20$Count, total.seeds = 20, intervals = WPS89_allo_20$Interval, partial = T)

head(WPS107_allo_20)
GermPercent(germ.counts = WPS107_allo_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_allo_20$Count, intervals = WPS107_allo_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_allo_20$Count, intervals = WPS107_allo_20$Interval, partial = T)
t50(germ.counts = WPS107_allo_20$Count, intervals = WPS107_allo_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_allo_20$Count, intervals = WPS107_allo_20$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_allo_20$Count, total.seeds = 20, intervals = WPS107_allo_20$Interval, partial = T)


head(WPS154_allo_20)
GermPercent(germ.counts = WPS154_allo_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_allo_20$Count, intervals = WPS154_allo_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_allo_20$Count, intervals = WPS154_allo_20$Interval, partial = T)
t50(germ.counts = WPS154_allo_20$Count, intervals = WPS154_allo_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_allo_20$Count, intervals = WPS154_allo_20$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_allo_20$Count, total.seeds = 20, intervals = WPS154_allo_20$Interval, partial = T)

# allo smoke 30 min


head(T1_allo_30)
GermPercent(germ.counts = T1_allo_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_allo_30$Count, intervals = T1_allo_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_allo_30$Count, intervals = T1_allo_30$Interval, partial = T)
t50(germ.counts = T1_allo_30$Count, intervals = T1_allo_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_allo_30$Count, intervals = T1_allo_30$Interval, partial = T)
WeightGermPercent(germ.counts = T1_allo_30$Count, total.seeds = 20, intervals = T1_allo_30$Interval, partial = T)

head(T2_allo_30)
GermPercent(germ.counts = T2_allo_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_allo_30$Count, intervals = T2_allo_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_allo_30$Count, intervals = T2_allo_30$Interval, partial = T)
t50(germ.counts = T2_allo_30$Count, intervals = T2_allo_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_allo_30$Count, intervals = T2_allo_30$Interval, partial = T)
WeightGermPercent(germ.counts = T2_allo_30$Count, total.seeds = 20, intervals = T2_allo_30$Interval, partial = T)


head(WPS73_allo_30)
GermPercent(germ.counts = WPS73_allo_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_allo_30$Count, intervals = WPS73_allo_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_allo_30$Count, intervals = WPS73_allo_30$Interval, partial = T)
t50(germ.counts = WPS73_allo_30$Count, intervals = WPS73_allo_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_allo_30$Count, intervals = WPS73_allo_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_allo_30$Count, total.seeds = 20, intervals = WPS73_allo_30$Interval, partial = T)



head(WPS89_allo_30)
GermPercent(germ.counts = WPS89_allo_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_allo_30$Count, intervals = WPS89_allo_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_allo_30$Count, intervals = WPS89_allo_30$Interval, partial = T)
t50(germ.counts = WPS89_allo_30$Count, intervals = WPS89_allo_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_allo_30$Count, intervals = WPS89_allo_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_allo_30$Count, total.seeds = 20, intervals = WPS89_allo_30$Interval, partial = T)

head(WPS107_allo_30)
GermPercent(germ.counts = WPS107_allo_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_allo_30$Count, intervals = WPS107_allo_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_allo_30$Count, intervals = WPS107_allo_30$Interval, partial = T)
t50(germ.counts = WPS107_allo_30$Count, intervals = WPS107_allo_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_allo_30$Count, intervals = WPS107_allo_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_allo_30$Count, total.seeds = 20, intervals = WPS107_allo_30$Interval, partial = T)


head(WPS154_allo_30)
GermPercent(germ.counts = WPS154_allo_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_allo_30$Count, intervals = WPS154_allo_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_allo_30$Count, intervals = WPS154_allo_30$Interval, partial = T)
t50(germ.counts = WPS154_allo_30$Count, intervals = WPS154_allo_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_allo_30$Count, intervals = WPS154_allo_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_allo_30$Count, total.seeds = 20, intervals = WPS154_allo_30$Interval, partial = T)


# pine smoke 5 min

head(T1_pine_5)
GermPercent(germ.counts = T1_pine_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_pine_5$Count, intervals = T1_pine_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_pine_5$Count, intervals = T1_pine_5$Interval, partial = T)
t50(germ.counts = T1_pine_5$Count, intervals = T1_pine_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_pine_5$Count, intervals = T1_pine_5$Interval, partial = T)
WeightGermPercent(germ.counts = T1_pine_5$Count, total.seeds = 20, intervals = T1_pine_5$Interval, partial = T)

head(T2_pine_5)
GermPercent(germ.counts = T2_pine_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_pine_5$Count, intervals = T2_pine_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_pine_5$Count, intervals = T2_pine_5$Interval, partial = T)
t50(germ.counts = T2_pine_5$Count, intervals = T2_pine_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_pine_5$Count, intervals = T2_pine_5$Interval, partial = T)
WeightGermPercent(germ.counts = T2_pine_5$Count, total.seeds = 20, intervals = T2_pine_5$Interval, partial = T)


head(WPS73_pine_5)
GermPercent(germ.counts = WPS73_pine_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_pine_5$Count, intervals = WPS73_pine_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_pine_5$Count, intervals = WPS73_pine_5$Interval, partial = T)
t50(germ.counts = WPS73_pine_5$Count, intervals = WPS73_pine_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_pine_5$Count, intervals = WPS73_pine_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_pine_5$Count, total.seeds = 20, intervals = WPS73_pine_5$Interval, partial = T)



head(WPS89_pine_5)
GermPercent(germ.counts = WPS89_pine_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_pine_5$Count, intervals = WPS89_pine_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_pine_5$Count, intervals = WPS89_pine_5$Interval, partial = T)
t50(germ.counts = WPS89_pine_5$Count, intervals = WPS89_pine_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_pine_5$Count, intervals = WPS89_pine_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_pine_5$Count, total.seeds = 20, intervals = WPS89_pine_5$Interval, partial = T)

head(WPS107_pine_5)
GermPercent(germ.counts = WPS107_pine_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_pine_5$Count, intervals = WPS107_pine_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_pine_5$Count, intervals = WPS107_pine_5$Interval, partial = T)
t50(germ.counts = WPS107_pine_5$Count, intervals = WPS107_pine_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_pine_5$Count, intervals = WPS107_pine_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_pine_5$Count, total.seeds = 20, intervals = WPS107_pine_5$Interval, partial = T)


head(WPS154_pine_5)
GermPercent(germ.counts = WPS154_pine_5$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_pine_5$Count, intervals = WPS154_pine_5$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_pine_5$Count, intervals = WPS154_pine_5$Interval, partial = T)
t50(germ.counts = WPS154_pine_5$Count, intervals = WPS154_pine_5$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_pine_5$Count, intervals = WPS154_pine_5$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_pine_5$Count, total.seeds = 20, intervals = WPS154_pine_5$Interval, partial = T)


# pine smoke 10 min


head(T1_pine_10)
GermPercent(germ.counts = T1_pine_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_pine_10$Count, intervals = T1_pine_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_pine_10$Count, intervals = T1_pine_10$Interval, partial = T)
t50(germ.counts = T1_pine_10$Count, intervals = T1_pine_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_pine_10$Count, intervals = T1_pine_10$Interval, partial = T)
WeightGermPercent(germ.counts = T1_pine_10$Count, total.seeds = 20, intervals = T1_pine_10$Interval, partial = T)

head(T2_pine_10)
GermPercent(germ.counts = T2_pine_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_pine_10$Count, intervals = T2_pine_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_pine_10$Count, intervals = T2_pine_10$Interval, partial = T)
t50(germ.counts = T2_pine_10$Count, intervals = T2_pine_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_pine_10$Count, intervals = T2_pine_10$Interval, partial = T)
WeightGermPercent(germ.counts = T2_pine_10$Count, total.seeds = 20, intervals = T2_pine_10$Interval, partial = T)


head(WPS73_pine_10)
GermPercent(germ.counts = WPS73_pine_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_pine_10$Count, intervals = WPS73_pine_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_pine_10$Count, intervals = WPS73_pine_10$Interval, partial = T)
t50(germ.counts = WPS73_pine_10$Count, intervals = WPS73_pine_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_pine_10$Count, intervals = WPS73_pine_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_pine_10$Count, total.seeds = 20, intervals = WPS73_pine_10$Interval, partial = T)



head(WPS89_pine_10)
GermPercent(germ.counts = WPS89_pine_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_pine_10$Count, intervals = WPS89_pine_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_pine_10$Count, intervals = WPS89_pine_10$Interval, partial = T)
t50(germ.counts = WPS89_pine_10$Count, intervals = WPS89_pine_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_pine_10$Count, intervals = WPS89_pine_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_pine_10$Count, total.seeds = 20, intervals = WPS89_pine_10$Interval, partial = T)

head(WPS107_pine_10)
GermPercent(germ.counts = WPS107_pine_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_pine_10$Count, intervals = WPS107_pine_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_pine_10$Count, intervals = WPS107_pine_10$Interval, partial = T)
t50(germ.counts = WPS107_pine_10$Count, intervals = WPS107_pine_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_pine_10$Count, intervals = WPS107_pine_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_pine_10$Count, total.seeds = 20, intervals = WPS107_pine_10$Interval, partial = T)


head(WPS154_pine_10)
GermPercent(germ.counts = WPS154_pine_10$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_pine_10$Count, intervals = WPS154_pine_10$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_pine_10$Count, intervals = WPS154_pine_10$Interval, partial = T)
t50(germ.counts = WPS154_pine_10$Count, intervals = WPS154_pine_10$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_pine_10$Count, intervals = WPS154_pine_10$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_pine_10$Count, total.seeds = 20, intervals = WPS154_pine_10$Interval, partial = T)



# pine smoke 20 min

head(T1_pine_20)
GermPercent(germ.counts = T1_pine_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_pine_20$Count, intervals = T1_pine_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_pine_20$Count, intervals = T1_pine_20$Interval, partial = T)
t50(germ.counts = T1_pine_20$Count, intervals = T1_pine_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_pine_20$Count, intervals = T1_pine_20$Interval, partial = T)
WeightGermPercent(germ.counts = T1_pine_20$Count, total.seeds = 20, intervals = T1_pine_20$Interval, partial = T)

head(T2_pine_20)
GermPercent(germ.counts = T2_pine_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_pine_20$Count, intervals = T2_pine_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_pine_20$Count, intervals = T2_pine_20$Interval, partial = T)
t50(germ.counts = T2_pine_20$Count, intervals = T2_pine_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_pine_20$Count, intervals = T2_pine_20$Interval, partial = T)
WeightGermPercent(germ.counts = T2_pine_20$Count, total.seeds = 20, intervals = T2_pine_20$Interval, partial = T)


head(WPS73_pine_20)
GermPercent(germ.counts = WPS73_pine_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_pine_20$Count, intervals = WPS73_pine_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_pine_20$Count, intervals = WPS73_pine_20$Interval, partial = T)
t50(germ.counts = WPS73_pine_20$Count, intervals = WPS73_pine_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_pine_20$Count, intervals = WPS73_pine_20$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_pine_20$Count, total.seeds = 20, intervals = WPS73_pine_20$Interval, partial = T)


head(WPS89_pine_20)
GermPercent(germ.counts = WPS89_pine_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_pine_20$Count, intervals = WPS89_pine_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_pine_20$Count, intervals = WPS89_pine_20$Interval, partial = T)
t50(germ.counts = WPS89_pine_20$Count, intervals = WPS89_pine_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_pine_20$Count, intervals = WPS89_pine_20$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_pine_20$Count, total.seeds = 20, intervals = WPS89_pine_20$Interval, partial = T)

head(WPS107_pine_20)
GermPercent(germ.counts = WPS107_pine_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_pine_20$Count, intervals = WPS107_pine_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_pine_20$Count, intervals = WPS107_pine_20$Interval, partial = T)
t50(germ.counts = WPS107_pine_20$Count, intervals = WPS107_pine_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_pine_20$Count, intervals = WPS107_pine_20$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_pine_20$Count, total.seeds = 20, intervals = WPS107_pine_20$Interval, partial = T)


head(WPS154_pine_20)
GermPercent(germ.counts = WPS154_pine_20$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_pine_20$Count, intervals = WPS154_pine_20$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_pine_20$Count, intervals = WPS154_pine_20$Interval, partial = T)
t50(germ.counts = WPS154_pine_20$Count, intervals = WPS154_pine_20$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_pine_20$Count, intervals = WPS154_pine_20$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_pine_20$Count, total.seeds = 20, intervals = WPS154_pine_20$Interval, partial = T)

# pine smoke 30 min


head(T1_pine_30)
GermPercent(germ.counts = T1_pine_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T1_pine_30$Count, intervals = T1_pine_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T1_pine_30$Count, intervals = T1_pine_30$Interval, partial = T)
t50(germ.counts = T1_pine_30$Count, intervals = T1_pine_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T1_pine_30$Count, intervals = T1_pine_30$Interval, partial = T)
WeightGermPercent(germ.counts = T1_pine_30$Count, total.seeds = 20, intervals = T1_pine_30$Interval, partial = T)

head(T2_pine_30)
GermPercent(germ.counts = T2_pine_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = T2_pine_30$Count, intervals = T2_pine_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = T2_pine_30$Count, intervals = T2_pine_30$Interval, partial = T)
t50(germ.counts = T2_pine_30$Count, intervals = T2_pine_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = T2_pine_30$Count, intervals = T2_pine_30$Interval, partial = T)
WeightGermPercent(germ.counts = T2_pine_30$Count, total.seeds = 20, intervals = T2_pine_30$Interval, partial = T)


head(WPS73_pine_30)
GermPercent(germ.counts = WPS73_pine_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS73_pine_30$Count, intervals = WPS73_pine_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS73_pine_30$Count, intervals = WPS73_pine_30$Interval, partial = T)
t50(germ.counts = WPS73_pine_30$Count, intervals = WPS73_pine_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS73_pine_30$Count, intervals = WPS73_pine_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS73_pine_30$Count, total.seeds = 20, intervals = WPS73_pine_30$Interval, partial = T)



head(WPS89_pine_30)
GermPercent(germ.counts = WPS89_pine_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS89_pine_30$Count, intervals = WPS89_pine_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS89_pine_30$Count, intervals = WPS89_pine_30$Interval, partial = T)
t50(germ.counts = WPS89_pine_30$Count, intervals = WPS89_pine_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS89_pine_30$Count, intervals = WPS89_pine_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS89_pine_30$Count, total.seeds = 20, intervals = WPS89_pine_30$Interval, partial = T)

head(WPS107_pine_30)
GermPercent(germ.counts = WPS107_pine_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS107_pine_30$Count, intervals = WPS107_pine_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS107_pine_30$Count, intervals = WPS107_pine_30$Interval, partial = T)
t50(germ.counts = WPS107_pine_30$Count, intervals = WPS107_pine_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS107_pine_30$Count, intervals = WPS107_pine_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS107_pine_30$Count, total.seeds = 20, intervals = WPS107_pine_30$Interval, partial = T)


head(WPS154_pine_30)
GermPercent(germ.counts = WPS154_pine_30$Count, total.seeds = 20)
PeakGermPercent(germ.counts = WPS154_pine_30$Count, intervals = WPS154_pine_30$Interval, total.seeds = 20)
PeakGermTime(germ.counts = WPS154_pine_30$Count, intervals = WPS154_pine_30$Interval, partial = T)
t50(germ.counts = WPS154_pine_30$Count, intervals = WPS154_pine_30$Interval, partial = T, method = c( 'coolbear', 'farooq'))
MeanGermTime(germ.counts = WPS154_pine_30$Count, intervals = WPS154_pine_30$Interval, partial = T)
WeightGermPercent(germ.counts = WPS154_pine_30$Count, total.seeds = 20, intervals = WPS154_pine_30$Interval, partial = T)
