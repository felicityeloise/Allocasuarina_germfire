# Written by Felicity Charles
# 4/4/2024
# Caveat emptor

# Thermogradient bar simple analysis of results
# https://academic.oup.com/aob/article/129/7/787/6536948 - package for analysing germination results

dat <- read.table("./00_Data/Thermogradient_test/Cumulative_thermogradient_results.txt", header = T)
head(dat); dim(dat)

dat$Time_spread <- dat$Time_to_finish - dat$Time_to_germ


dat.seed <- dat[,which(colnames(dat) == "X2"): which(colnames(dat) == "X26")]

cum.sum <- data.frame(t(apply(dat.seed, 1, FUN = function(x) cumsum(unlist(x)))))
cum.prop <- cum.sum/dat$Tot_seed
cum.prop2<- cbind(dat[,c(1:3, 18)], cum.prop)
head(cum.prop2); dim(cum.prop2)
str(cum.prop2)

# Split the data by species for plotting

cum.prop2_lit <- subset(cum.prop2, cum.prop2$Species == "littoralis")
rownames(cum.prop2_lit) <- 1:nrow(cum.prop2_lit)
head(cum.prop2_lit); dim(cum.prop2_lit)



cum.prop2_tor <- subset(cum.prop2, cum.prop2$Species == "torulosa")
rownames(cum.prop2_tor) <- 1:nrow(cum.prop2_tor)
head(cum.prop2_tor); dim(cum.prop2_tor)




library(RColorBrewer)

pal <- brewer.pal(10, "RdYlBu")
pal <- c("#A50026","#D73027","#F46D43","#FDAE61", "#FEE090", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4","#313695")
chamber.col<- data.frame(col = rev(pal), chamber = unique(cum.prop2$Chamber))
prop.dat_lit <- cum.prop2_lit[,which(colnames(cum.prop2_lit) == "X2"): which(colnames(cum.prop2_lit) == "X26")]
prop.dat_tor <- cum.prop2_tor[,which(colnames(cum.prop2_tor) == "X2"): which(colnames(cum.prop2_tor) == "X26")]



head(cum.prop2_lit);dim(cum.prop2_lit)
head(cum.prop2_tor);dim(cum.prop2_tor)



dev.new(height = 6, width = 12, dpi = 80, pointsize = 16, noRStudioGD = T)
par(mar = c(4,4,2,1), mgp = c(3,1,0), mfrow = c(1,2), oma = c(0,0,0, 8))
plot(1:ncol(prop.dat_lit), 1:ncol(prop.dat_lit), ylim = c(0,1), type = "n", las = 1, ylab = "Proportion germinated", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.dat_lit), labels = c(2,5,7,9,12,14,20,21,23,26,28), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(xlab = "Days", mgp = c(2.2, 1, 0))
title(main = "(a) Allocasuarina littoralis")

for(i in 1:nrow(cum.prop2_lit)){
  dat.thisrun <- cum.prop2_lit[i,5:ncol(cum.prop2_lit)]
  chamber.thisrun <- cum.prop2_lit$Chamber[i]
  col.thisrun <- chamber.col.lit$col[chamber.col.lit$chamber == chamber.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}


plot(1:ncol(prop.dat_tor), 1:ncol(prop.dat_tor), ylim = c(0,1), type = "n", las = 1, ylab = "Proportion germinated", xaxt = "n", xlab = "")
axis(side = 1, at = 1:ncol(prop.dat_tor), labels = c(2,5,7,9,12,14,20,21,23,26,28), cex.axis = 0.8, mgp = c(2.2, 0.7, 0))
title(xlab = "Days", mgp = c(2.2, 1, 0))
title(main = "(b) Allocasuarina torulosa")
for(i in 1:nrow(cum.prop2_tor)){
  dat.thisrun <- cum.prop2_tor[i,5:ncol(cum.prop2_tor)]
  chamber.thisrun <- cum.prop2_tor$Chamber[i]
  col.thisrun <- chamber.col$col[chamber.col$chamber == chamber.thisrun]
  lines(1:length(dat.thisrun), dat.thisrun, col = col.thisrun, lwd = 3)
}


par(xpd = NA)
legend(x = 10.5, y= 1, legend = c("6°C", "10°C", "14°C", "17°C", "20°C", "23°C", "26°C", "29°C", "32°C", "36°C"), col = chamber.col$col, lty = 1, lwd = 4, cex = 0.95, text.width = 0.2, title = "Chamber temperature", bty = "n")
par(xpd = F)





# Combine the data from cumulative proportions with the data on time spread of germination

dat2 <- cbind(dat[,c(1:4, 18:21)], cum.prop)
# Check this looks right
head(dat2)
head(cum.prop2_tor)

write.table(dat2, "./00_Data/Thermogradient_test/Cumulative_thermogradient_results.txt", col.names = T, row.names = F, sep = ",")





