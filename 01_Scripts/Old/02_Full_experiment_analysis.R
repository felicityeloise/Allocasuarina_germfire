# Written by Felicity Charles
# Caveat emptor
# Date: 30th October 2024

# Full experiment analysis


#1. Load packages
# Having a problem with this germinationmetrics::germination.indices(). A warning message was updated but does not seem to be showing the update message when I run the function. Cannot install from github.

library(dplyr)
library(lme4)
library(terra)
library(sf)
library(AICcmodavg)


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



# 4.  Germination metrics ----
head(dat_cum.prop)

# We cannot get the germinationmetrics::germination.indices functionality to work as it either fails because of some incorrect error thrown when t50 is included, or it decides that the resulting vector would be too long, despite working for the a control only group. We will instead create our own functions to run the calculations as even doing a for loop with the germinationmetrics standalone functions was not working correctly in some instances such as t50. 

# Specify the number of intervals
int <- 1:61

# Caculate total proportion germination
head(dat_cum.prop); dim(dat_cum.prop)
dat_cum.prop$Proportion_germ <- dat_cum.prop[, 73]
unique(dat_cum.prop$Proportion_germ) # Replace NAs with 0
dat_cum.prop$Proportion_germ[is.na(dat_cum.prop$Proportion_germ)] <- 0
unique(dat_cum.prop$Proportion_germ) # Check
head(dat_cum.prop); dim(dat_cum.prop)

# Calculate 50% germination 
dat_cum.prop$Perc50 <- dat_cum.prop$Proportion_germ/2



# Run calculations for each germination metric for each row of data
for(i in 1:nrow(dat_cum.prop)){
  dat_cum.prop$t50[i] <- paste(colnames(dat_cum.prop[i, 13:73])[which(dat_cum.prop[i, 13:73] >= dat_cum.prop$Perc50)])[1]
  dat_cum.prop$t50 <- as.numeric(sub("Day", "", dat_cum.prop$t50)) 
  dat_cum.prop$firstgerm[i] <- int[min(which(dat_cum.prop[i, 13:73] != 0))]
  dat_cum.prop$lastgerm[i] <- int[min(which(dat_cum.prop[i, 13:73] == dat_cum.prop$Total_germination))]
  dat_cum.prop$timespread[i] <-  int[min(which(dat_cum.prop[i, 13:73] == dat_cum.prop$Total_germination))] - int[min(which(dat_cum.prop[i, 13:73] != 0))]
  dat_cum.prop$germspeed[i] <- sum((dat_cum.prop[i, 13:73]/dat_cum.prop$Total_seeds)/int)
  dat_cum.prop[is.na(dat_cum.prop)] <- 0
} 

head(dat_cum.prop); tail(dat_cum.prop); dim(dat_cum.prop)

str(dat_cum.prop)
hist(dat_cum.prop$Proportion_germ[dat_cum.prop$Species == "littoralis"])
hist(dat_cum.prop$Proportion_germ[dat_cum.prop$Species == "torulosa"])


summary(dat_cum.prop$Proportion_germ)
View(dat_cum.prop[, c(1, 72:ncol(dat_cum.prop))])


# For each analysis we want to include a random effect for individual and replicate (set) accounting for any differences between individuals and germination between replicates. 



# 5. Add extra information columns for analyses ----
# Add information on seed weight, and fire frequency
seed_charact <- read.csv('./00_Data/Seed characteristics and sample information.csv', stringsAsFactors = T, header = T)
head(seed_charact)
seed_charact$seed_weight <- seed_charact$Seed_weight_subset_of_10_seeds/10 # Calculate weight per seed
seed_charact_s <- seed_charact[, c(5,23)]
head(seed_charact_s)



# 5.1 Load fire data ----
# We need to make decisions about what fire data to use, whether we need to keep QPWS data where it is available and whether to use the GLM or GAM predictions. We did conclude that the GAM predictions were better for investigations of fire frequency but lets just confirm this.
pts <- vect('D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/QPWS_random.gpkg')
GLM_pred <- rast('D:/PhD/R_analysis/Fire_freq/04_Results/Prediction_rasters/GLM_pred.tif')
GAM_pred <- rast('D:/PhD/R_analysis/Fire_freq/04_Results/Prediction_rasters/GAM_pred.tif')
QPWS <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_raster.tif')



# 5.2 Extract fire frequencies to points for seed lots and transects ----
transects <- read.csv('C:/Users/s4590925/OneDrive - The University of Queensland/Desktop/GitHub/Fire_recruit/00_Data/Transect_location_data.csv', stringsAsFactors = T, header = T)
transect_sf <- st_as_sf(transects, coords = c('Longitude', 'Latitude'), crs = 'EPGS:4326')
head(transect_sf)
transect_sf <- transect_sf[-6,] # Remove the row for leuhmannii
transect_sf
transect_sf <- sf::st_set_crs(transect_sf, 'EPSG:4326')
transect <- vect(transect_sf) %>% 
  project('EPSG:3577')


dat_sf <- st_as_sf(seed_charact, coords = c('Longitude', 'Latitude'), crs = 'EPSG:4326')
head(dat_sf)
dat_sf <- st_make_valid(dat_sf)
germ <- vect(dat_sf) %>% 
  project('EPSG:3577') 
germ
plet(germ)


GLM_g <- extract(GLM_pred, germ)
colnames(GLM_g) <- c("ID", "Fire_freq")
GLM_g$Data <- "GLM"
GLM_g$transect <- germ$Transect
GLM_g$Species <- germ$Species

GLM_t <- extract(GLM_pred, transect)
colnames(GLM_t) <- c("ID", "Fire_freq")
GLM_t$Data <- "GLM"
GLM_t$transect <- transect$Transect
GLM_t$Species <- transect$Species


GAM_g <- extract(GAM_pred, germ)
colnames(GAM_g) <- c("ID", "Fire_freq")
GAM_g$Data <- "GAM"
GAM_g$transect <- germ$Transect
GAM_g$Species <- germ$Species

GAM_t <- extract(GAM_pred, transect)
colnames(GAM_t) <- c("ID", "Fire_freq")
GAM_t$Data <- "GAM"
GAM_t$transect <- transect$Transect
GAM_t$Species <- transect$Species

QPWS_g <- extract(QPWS, germ)
colnames(QPWS_g) <- c("ID", "Fire_freq")
QPWS_g[is.na(QPWS_g)] <- 0
QPWS_g$Data <- "QPWS"
QPWS_g$transect <- germ$Transect
QPWS_g$Species <- germ$Species

QPWS_t <- extract(QPWS, transect)
colnames(QPWS_t) <- c("ID", "Fire_freq")
QPWS_t[is.na(QPWS_t)] <- 0
QPWS_t$Data <- "QPWS"
QPWS_t$transect <- transect$Transect
QPWS_t$Species <- transect$Species



# Combine into one dataframe
transect_fire <- rbind(GLM_t, GAM_t, QPWS_t)
tor_tran_fire <- transect_fire[transect_fire$Species == "torulosa",]
lit_tran_fire <- transect_fire[transect_fire$Species == "littoralis",]

germ_fire <- rbind(GLM_g, GAM_g, QPWS_g)
tor_germ_fire <- germ_fire[germ_fire$Species == "torulosa",]
lit_germ_fire <- germ_fire[germ_fire$Species == "littoralis",]


# 5.3 Produce and plot histograms for each species and dataset ----
QP_tt <- hist(tor_tran_fire$Fire_freq[tor_tran_fire$Data == "QPWS"], breaks = seq(-1, 7, 1))
GL_tt <- hist(tor_tran_fire$Fire_freq[tor_tran_fire$Data == "GLM"], breaks = seq(-1, 7, 1))
GA_tt <- hist(tor_tran_fire$Fire_freq[tor_tran_fire$Data == "GAM"], breaks = seq(-1, 7, 1))

QP_tg <- hist(tor_germ_fire$Fire_freq[tor_germ_fire$Data == "QPWS"], breaks = seq(-1, 7, 1))
GL_tg <- hist(tor_germ_fire$Fire_freq[tor_germ_fire$Data == "GLM"], breaks = seq(-1, 7, 1))
GA_tg <- hist(tor_germ_fire$Fire_freq[tor_germ_fire$Data == "GAM"], breaks = seq(-1, 7, 1))



QP_lt <- hist(lit_tran_fire$Fire_freq[lit_tran_fire$Data == "QPWS"], breaks = seq(-1,7,1))
GL_lt <- hist(lit_tran_fire$Fire_freq[lit_tran_fire$Data == "GLM"], breaks = seq(-1,7,1))
GA_lt <- hist(lit_tran_fire$Fire_freq[lit_tran_fire$Data == "GAM"], breaks = seq(-1,7,1))

QP_lg <- hist(lit_germ_fire$Fire_freq[lit_germ_fire$Data == "QPWS"], breaks = seq(-1,7,1))
GL_lg <- hist(lit_germ_fire$Fire_freq[lit_germ_fire$Data == "GLM"], breaks = seq(-1,7,1))
GA_lg <- hist(lit_germ_fire$Fire_freq[lit_germ_fire$Data == "GAM"], breaks = seq(-1,7,1))




# Plot histograms together for comparisons
par(mfrow = c(1,2))

# Torulosa
# Transects
plot(QP_tt, col = 'gray', ylim = c(0,30), main = "QPWS + GLM Transects", xlab = "Fire frequency", xaxt = "n", las = 1)
plot(GL_tt, col = rgb(73/255, 32/255, 80/255, 0.6), add = T)
axis(side = 1, at = seq(-1,7,1))

plot(QP_tt, col = 'gray', ylim = c(0,30),  main = "QPWS + GAM Transects", xlab = "Fire frequency", xaxt = "n", las = 1)
plot(GA_tt, col = rgb(170/255, 169/255, 112/255, 0.6), add = T)
axis(side = 1, at = seq(-1,7,1))


# Seed lots 
plot(QP_tg, col = 'gray', ylim = c(0,40), main = "QPWS + GLM Seed lots", xlab = "Fire frequency", xaxt = "n", las = 1)
plot(GL_tg, col = rgb(73/255, 32/255, 80/255, 0.6), add = T)
axis(side = 1, at = seq(-1,7,1))

plot(QP_tg, col = 'gray', ylim = c(0,40),  main = "QPWS + GAM Seed lots", xlab = "Fire frequency", xaxt = "n", las = 1)
plot(GA_tg, col = rgb(170/255, 169/255, 112/255, 0.6), add = T)
axis(side = 1, at = seq(-1,7,1))



# Littoralis
# Transects
plot(QP_lt, col = 'gray', ylim = c(0,30), main = "QPWS + GLM Transects", xlab = "Fire frequency", xaxt = "n", las = 1)
plot(GL_lt, col = rgb(73/255, 32/255, 80/255, 0.6), add = T)
axis(side = 1, at = seq(-1,7,1))

plot(QP_lt, col = 'gray', ylim = c(0,30),  main = "QPWS + GAM Transects", xlab = "Fire frequency", xaxt = "n", las = 1)
plot(GA_lt, col = rgb(170/255, 169/255, 112/255, 0.6), add = T)
axis(side = 1, at = seq(-1,7,1))


# Seed lots 
plot(QP_lg, col = 'gray', ylim = c(0,40), main = "QPWS + GLM Seed lots", xlab = "Fire frequency", xaxt = "n", las = 1)
plot(GL_lg, col = rgb(73/255, 32/255, 80/255, 0.6), add = T)
axis(side = 1, at = seq(-1,7,1))

plot(QP_lg, col = 'gray', ylim = c(0,40),  main = "QPWS + GAM Seed lots", xlab = "Fire frequency", xaxt = "n", las = 1)
plot(GA_lg, col = rgb(170/255, 169/255, 112/255, 0.6), add = T)
axis(side = 1, at = seq(-1,7,1))



# So we want to use QPWS data wherever it is available, where we have NA values for QPWS we want to use the GAM as it is more likely to capture higher fire activity and reduces the overall shift of high fire frequencies to lower fire frequencies. 




# 5.4 Extract Time since fire information ----
# NOTE: We can only extract date of last fire for QPWS and unmodelled satellite data where QPWS unavailable. So we need to start with QPWS data and then supplement this with satellite data. 

# To extract time since fire, we will convert the QPWS fire history data to raster, with raster cell values assigned as the maximum value for year burnt. 
#QPWS_hist <- st_read('D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/QPWS_fire_hist_1987.gpkg') %>% 
  st_transform('EPSG:3577')
#QPWS_hist <- st_crop(QPWS_hist, GAM_pred)
#rtemp <- raster::raster(xmn = 1902033, xmx = 2111776, ymn = -3257627, ymx = -2954985, res = 5, crs = 'EPSG:3577')
#QPWS_SEQ_freq_rast <- fasterize::fasterize(QPWS_hist, rtemp, field = 'OUTYEAR', fun = 'max')
#raster::writeRaster(QPWS_SEQ_freq_rast, 'D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/SEQ/QPWS_TSF.tif')


# Extract QPWS fire history for transects and seed lots.
QPWS_hist_rast <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/SEQ/QPWS_TSF.tif')
transect_TSF <- extract(QPWS_hist_rast, transect, xy = T)
germ_TSF <- extract(QPWS_hist_rast, germ, xy = T)

# Find which points have NAs to help with replacement
t_TSF_na <-  transect_TSF[is.na(transect_TSF$QPWS_TSF),]
g_TSF_na <- germ_TSF[is.na(germ_TSF$QPWS_TSF), ]

sent <- rast('D:/ADATA/Fire_data/QLD/Landsat_Annual_Fire_Scars/TERN_Sentinel2/cvmsre_qld_2023_afma2(1).tif')
t_TSF2 <- extract(sent, transect)
t_TSF2

g_TSF2 <- extract(sent, germ)
g_TSF2
#96,97,100,102,103, = 2011, 78 = 2006, 37,38,39,68,69,70,71,72,73 = 2004, 74 = 2003, 75,76,77 = 1999


# Replace NA values with year burnt from Sentinel data, anything that was not burnt between 1987-2023 will be assigned the value of 1986
transect_TSF[22, 2] <- 2004
transect_TSF[30, 2] <- 2020
transect_TSF[31, 2] <- 2004
transect_TSF[32, 2] <- 2003
transect_TSF[is.na(transect_TSF$QPWS_TSF), 2] <- 1986
transect_TSF

germ_TSF[c(96,97,100,102,103), 2] <- 2011
germ_TSF[78, 2] <- 2006
germ_TSF[c(37,38,39,68,70,71,72,73), 2] <- 2004
germ_TSF[74, 2] <- 2003
germ_TSF[c(75,76,77), 2] <- 1999
germ_TSF[is.na(germ_TSF$QPWS_TSF), 2] <- 1986
germ_TSF

# Now calculate time since fire
transect_TSF$Time_SF <- 2023-transect_TSF$QPWS_TSF
transect_TSF
transect_TSF$Transect <- transect$Transect
transect_TSF$Species <- transect$Species


germ_TSF 





hist(transect_TSF$Time_SF)
transect_TSF_l <- transect_TSF[transect_TSF$Species == "littoralis",]
transect_TSF_t <- transect_TSF[transect_TSF$Species == "torulosa", ]
hist(transect_TSF_l$Time_SF)
hist(transect_TSF_t$Time_SF)

# Then we want to look at boxplots of germination rates and seed size
lit <- seed_charact[seed_charact$Species == "littoralis", ]
tor <- seed_charact[seed_charact$Species == "torulosa", ]

# Seed weight
hist(lit$seed_weight)
hist(tor$seed_weight)


# Germination rate - time to 50% germination

germ_lit <- dat_cum.prop[dat_cum.prop$Species == "littoralis", ]
germ_tor <- dat_cum.prop[dat_cum.prop$Species == "torulosa", ]

hist(germ_lit$t50)
hist(germ_tor$t50)

hist(germ_lit$germspeed)
hist(germ_tor$germspeed)

hist(germ_lit$Proportion_germ)
hist(germ_tor$Proportion_germ)

















# Need to decide whether to keep information from QPWS as well or if this doesn't matter as much for the analyses as we are analysing all the individuals for a species together
pred_fire <- extract(pred_ff, germ)
seed_charact_s$Fire_freq <- pred_fire$lyr1
head(seed_charact_s)
unique(seed_charact_s$Fire_freq)
unique(is.na(seed_charact_s))
seed_charact_s$Fire_freq <- round(seed_charact_s$Fire_freq)
head(seed_charact_s)











# Check how the correlation looks
pred_pts <- extract(pred_ff, pts)
QPWS <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped_reproj.tif')

QPWS_pts <- extract(QPWS, pts)
GLM_pts <- extract(GLM_pred, pts)
GAM_pts <- extract(GAM_pred, pts)

cor.test(pred_pts$lyr1, QPWS_pts$QPWS_SEQ_freq_raster) # This correlation is lower than that of the GLM predictions but is higher than the GAM predictions and we already know the GAM produces values more similar to QPWS at higher fire frequencies.
cor.test(GLM_pts$lyr1, QPWS_pts$QPWS_SEQ_freq_raster)
cor.test(GAM_pts$lyr1, QPWS_pts$QPWS_SEQ_freq_raster)

# Upon further investigations, due to categorising torulosa individuals as high fire or low fire, this was based on fire frequency from QPWS. To preserve these categories we need to use QPWS data as well where it is available

QPWS_ff <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped_reproj.tif')
fire_freq <- mosaic(pred_ff, QPWS_ff, fun = 'max')
plot(fire_freq)


dat_sf <- st_as_sf(seed_charact, coords = c('Longitude', 'Latitude'), crs = 'EPSG:4326')
head(dat_sf)
dat_sf <- st_make_valid(dat_sf)
tmap::qtm(dat_sf)

# Transform to terra and reproject
germ <- vect(dat_sf) %>% 
  project('EPSG:3577') 
germ
plet(germ)

# Need to decide whether to keep information from QPWS as well or if this doesn't matter as much for the analyses as we are analysing all the individuals for a species together
pred_fire <- extract(pred_ff, germ)
seed_charact_s$Fire_freq <- pred_fire$lyr1
head(seed_charact_s)
unique(seed_charact_s$Fire_freq)
unique(is.na(seed_charact_s))
seed_charact_s$Fire_freq <- round(seed_charact_s$Fire_freq)
head(seed_charact_s)


dat_cum.prop <- left_join(dat_cum.prop, seed_charact_s, by = "Individual")
head(dat_cum.prop)
unique(dat_cum.prop$Fire_freq)


head(dat_cum.prop[1, ])

View(dat_cum.prop[dat_cum.prop$Species == "torulosa" & dat_cum.prop$Fire_freq >= 4, c(1,2, 5,82)])



# Notes on analyses:
# As zeros are "non-positive" values, GLMER will not accept a Gamma distribution for this data, so we must instead use a binomial family instead. 

# We use Akaike's information criterion corrected for small sample sizes to rank models with the best model chosen as the model with the lowest AICc, which improves model fit over the null model by a change in AICc > 2. 

# For all models as we are investigating two different species, with different modes of reproduction so we include a model term for specie in all models. 



save.image('./02_Workspaces/Full_experiment_analysis.RData')
