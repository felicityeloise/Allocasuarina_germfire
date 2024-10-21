# Written by Felicity Charles
# 10/09/2024


viable <- read.table('./00_Data/Full_experiment/Viability.txt', header = T, stringsAsFactors = T)
head(viable)

library(dplyr)
 

library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)
RColorBrewer::brewer.pal(11, "RdYlBu")

pal1 <- c("#D73027", "#FDAE61", "#4575B4")



cor.test(viable$X.ray_proportion_viable[viable$Species == "littoralis"], viable$TTC_proportion_viable[viable$Species == "littoralis"])
cor.test(viable$X.ray_proportion_viable[viable$Species == "torulosa" & viable$Group == "lowfi"], viable$TTC_proportion_viable[viable$Species == "torulosa" & viable$Group == "lowfi"])
cor.test(viable$X.ray_proportion_viable[viable$Species == "torulosa" & viable$Group == "hifi"], viable$TTC_proportion_viable[viable$Species == "torulosa" & viable$Group == "hifi"])



# Viability plots



dev.new(height = 7, width = 25, dpi = 80, pointsize = 18, noRStudioGD = T)
par(mar = c(5,5,2,3), mgp = c(3,1,0), mfrow = c(1, 3), oma = c(0,0,0,9), cex.axis = 1.6, cex.lab = 1.8, cex.main = 2)
 
plot(viable$X.ray_proportion_viable[viable$Species == "littoralis"], viable$TTC_proportion_viable[viable$Species == "littoralis"], pch = NA, xlab = expression(bold("X-ray")), ylab = expression(bold("Tetrazolium")), xlim = c(0,1), ylim = c(0,1), las = 1)
points(viable$X.ray_proportion_viable[viable$Species == "littoralis" & viable$Replicate == 1], viable$TTC_proportion_viable[viable$Species == "littoralis" & viable$Replicate == 1], col = "#D73027", pch = 20, cex = 2)
points(viable$X.ray_proportion_viable[viable$Species == "littoralis" & viable$Replicate == 2], viable$TTC_proportion_viable[viable$Species == "littoralis" & viable$Replicate == 2], col = "#FDAE61", pch = 20, cex = 2)
points(viable$X.ray_proportion_viable[viable$Species == "littoralis" & viable$Replicate == 3], viable$TTC_proportion_viable[viable$Species == "littoralis" & viable$Replicate == 3], col = "#4575B4", pch = 20, cex = 2)
text(0.05, 1, labels = "Pearson's r = 0.71", adj = 0, cex = 1.5)
text(0.16, 0.9, labels = "p < 0.001", cex = 1.5)
title('(a)', adj = 0, cex = 3)

plot(viable$X.ray_proportion_viable[viable$Species == "torulosa" & viable$Group == "lowfi"], viable$TTC_proportion_viable[viable$Species == "torulosa" & viable$Group == "lowfi"], pch = NA, xlab = expression(bold("X-ray")), ylab = "", xlim = c(0,1), ylim = c(0,1), las = 1)
points(viable$X.ray_proportion_viable[viable$Species == "torulosa" & viable$Group == "lowfi" & viable$Replicate == 1], viable$TTC_proportion_viable[viable$Species == "torulosa" & viable$Group == "lowfi" &viable$Replicate == 1], col = "#D73027", pch = 20, cex = 2)
points(viable$X.ray_proportion_viable[viable$Species == "torulosa" & viable$Group == "lowfi" & viable$Replicate == 2], viable$TTC_proportion_viable[viable$Species == "torulosa" & viable$Group == "lowfi" & viable$Replicate == 2], col = "#FDAE61", pch = 20, cex = 2)
points(viable$X.ray_proportion_viable[viable$Species == "torulosa" & viable$Group == "lowfi" & viable$Replicate == 3], viable$TTC_proportion_viable[viable$Species == "torulosa" & viable$Group == "lowfi" & viable$Replicate == 3], col = "#4575B4", pch = 20, cex = 2)
text(0.05, 1, labels = "Pearson's r = 0.81", adj = 0, cex = 1.5)
text(0.16, 0.9, labels = "p < 0.001", cex = 1.5)
title('(b)', adj = 0, cex = 3)



plot(viable$X.ray_proportion_viable[viable$Species == "torulosa" & viable$Group == "hifi"], viable$TTC_proportion_viable[viable$Species == "torulosa" & viable$Group == "hifi"], pch = NA, xlab = expression(bold("X-ray")), ylab = "", xlim = c(0,1), ylim = c(0,1), las = 1)
points(viable$X.ray_proportion_viable[viable$Species == "torulosa" & viable$Group == "hifi" & viable$Replicate == 1], viable$TTC_proportion_viable[viable$Species == "torulosa" & viable$Group == "hifi" &viable$Replicate == 1], col = "#D73027", pch = 20, cex = 2)
points(viable$X.ray_proportion_viable[viable$Species == "torulosa" & viable$Group == "hifi" & viable$Replicate == 2], viable$TTC_proportion_viable[viable$Species == "torulosa" & viable$Group == "hifi" & viable$Replicate == 2], col = "#FDAE61", pch = 20, cex = 2)
points(viable$X.ray_proportion_viable[viable$Species == "torulosa" & viable$Group == "hifi" & viable$Replicate == 3], viable$TTC_proportion_viable[viable$Species == "torulosa" & viable$Group == "hifi" & viable$Replicate == 3], col = "#4575B4", pch = 20, cex = 2)
text(0.05, 1, labels = "Pearson's r = 0.61", adj = 0, cex = 1.5)
text(0.16, 0.9, labels = "p < 0.001", cex = 1.5)
title('(c)', adj = 0, cex = 3)

par(xpd = NA)
legend(x = 1.1, y = 0.70, legend = c("Replicate 1", "Replicate 2", "Replicate 3"), col = c("#D73027", "#FDAE61", "#4575B4"), pch = 20, title = expression(bold("Replicate")), cex = 1.5)








# Lets change the structure of the data
levels(viable$Species) <- c("Allocasuarina littoralis", "Allocasuarina torulosa")
levels(viable$Group) <- c("high fire", "low fire")
viable$Treatment <- factor(viable$Treatment, levels = c("C", "H80", "H95", "S", "H80+S", "H95+S"))

str(viable)
head(viable)

library(dplyr)
options(scipen = 999)
# Extract the correlations for each species group
Viability_averages <- viable %>% 
  group_by(Species, Group) %>% 
  summarise(mean_TTC_via = mean(TTC_proportion_viable), X_ray_via = mean(X.ray_proportion_viable)) %>% 
  as.data.frame()
Viability_averages

Viability_correlations <- viable %>% 
  group_by(Species, Group) %>% 
  summarise(mean_TTC_viability = mean(TTC_proportion_viable), 
            X_ray_viability = mean(X.ray_proportion_viable),
            correlation = stats::cor.test(TTC_proportion_viable, X.ray_proportion_viable)$estimate,
            pvalue = stats::cor.test(TTC_proportion_viable, X.ray_proportion_viable)$p.value) %>% 
  as.data.frame()
Viability_correlations



# Now lets look into whether viability is changed by the treatment
Treatment_viability_averages_ind <- viable %>% 
  group_by(Species, Group, Individual, Treatment) %>% 
  summarise(X_ray_via = mean(X.ray_proportion_viable),
            mean_TTC_via = mean(TTC_proportion_viable),
            min_TTC = min(TTC_proportion_viable), 
            max_TTC = max(TTC_proportion_viable)) %>% 
  as.data.frame()
Treatment_viability_averages_ind



# This returns the correlation of x ray viability compared to TTC viability for each treatment per species grouping
            
Treatment_viability_cor <- Treatment_viability_averages_ind %>% 
  group_by(Species, Group, Treatment) %>% 
  summarise(min_TTC = min(min_TTC),
            max_TTC = max(max_TTC),
            mean_TTC = mean(mean_TTC_via),
            x_ray = mean(X_ray_via),
            correlation = stats::cor.test(mean_TTC_via, X_ray_via)$estimate,
            pvalue = stats::cor.test(mean_TTC_via, X_ray_via)$p.value)%>% 
  as.data.frame()
Treatment_viability_cor

Treatment_viability <- Treatment_viability_cor



Treatment_viability$mean_TTC <- round(Treatment_viability$mean_TTC, digits = 3)
Treatment_viability$x_ray <- round(Treatment_viability$x_ray, digits = 3)
Treatment_viability$correlation <- round(Treatment_viability$correlation, digits = 3)
Treatment_viability$pvalue <- round(Treatment_viability$pvalue, digits = 3)



options(scipen = 0)

str(Treatment_viability_cor)
levels(Treatment_viability_cor$Treatment)
# Lets plot these results, we want to facet by species grouping and colour the points by treatment type.
# First split the data

library(RColorBrewer)
pal <- c("#FEE090", "#F46D43", "#A50026", "#ABD9E9","#74ADD1","#313695")
  
brewer.pal(11, "RdYlBu")



cor.test(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis"])

cor.test(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire"])


cor.test(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire"])



 
dev.new(height = 4, width = 14, dpi = 80, pointsize = 16, noRStudioGD = T)
par(mar = c(5,5,2,3), mgp = c(3,1,0), mfrow = c(1,3), oma = c(0,0,0, 10), cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)



plot(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis"], pch = NA, xlim = c(0,1), ylim = c(0,1), xlab = expression(bold("X-ray")), ylab = expression(bold("Tetrazolium")), las = 1)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "C"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "C"], col = "#FEE090", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "H80"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "C"], col = "#F46D43", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "H95"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "C"], col = "#A50026", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "S"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "C"], col = "#ABD9E9", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "H80+S"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "C"], col = "#74ADD1", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "H95+S"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == "Allocasuarina littoralis" & Treatment_viability_averages_ind$Treatment == "C"], col = "#313695", pch = 20, cex = 2)
text(0.4, 0.4, labels = "Pearson's r = 0.85", adj = 0, cex = 1.2)
text(0.4, 0.3, labels = "p < 0.001", adj = 0, cex = 1.2)
title('(a)', adj = 0, cex = 2)





plot(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire"], pch = NA, xlim = c(0,1), ylim = c(0,1), xlab = expression(bold("X-ray")), ylab = "", las = 1)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "C"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "C"], col = "#FEE090", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "H80"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "H80"], col = "#F46D43", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "H95"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "H95"], col = "#A50026", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "S"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "S"], col = "#ABD9E9", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "H80+S"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "H80+S"], col = "#74ADD1", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "H95+S"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "low fire" & Treatment_viability_averages_ind$Treatment == "H95+S"], col = "#313695", pch = 20, cex = 2)
text(0.4, 0.2, labels = "Pearson's r = 0.88", adj = 0, cex = 1.2)
text(0.4, 0.1, labels = "p < 0.001", adj = 0, cex = 1.2)
title('(b)', adj = 0, cex = 2)






plot(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire"], pch = NA, xlim = c(0,1), ylim = c(0,1), xlab = expression(bold("X-ray")), ylab = "", las = 1)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "C"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "C"], col = "#FEE090", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "H80"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "H80"], col = "#F46D43", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "H95"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "H95"], col = "#A50026", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "S"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "S"], col = "#ABD9E9", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "H80+S"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "H80+S"], col = "#74ADD1", pch = 20, cex = 2)
points(Treatment_viability_averages_ind$mean_TTC_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "H95+S"], Treatment_viability_averages_ind$X_ray_via[Treatment_viability_averages_ind$Species == 'Allocasuarina torulosa' & Treatment_viability_averages_ind$Group == "high fire" & Treatment_viability_averages_ind$Treatment == "H95+S"], col = "#313695", pch = 20, cex = 2)
text(0.4, 0.2, labels = "Pearson's r = 0.64", adj = 0, cex = 1.2)
text(0.4, 0.1, labels = "p < 0.001", adj = 0, cex = 1.2)
title('(c)', adj = 0, cex = 2)

par(xpd = NA)
legend(x = 1.2, y = 1, legend = c("Control", "80°C", "95°C", "Smoke", "80°C+Smoke", "95°C+Smoke"), col = pal, title = expression(bold("Treatment")), pch = 19, cex = 1.2)


