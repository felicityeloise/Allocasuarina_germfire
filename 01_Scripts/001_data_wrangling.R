# Written by Felicity Charles
# Caveat emptor
# Date: 30th October 2024

# Data wrangling


#1. Load packages ----
# Having a problem with this germinationmetrics::germination.indices(). A warning message was updated but does not seem to be showing the update message when I run the function. Cannot install from github.

library(dplyr)
library(terra)
library(sf)
library(gdalUtilities)


# 2. Read in data and combine by species grouping ----
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

# 3. Calculate the cumulative proportions of germination ----
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
dat_cum.prop$Proportion_germ <- dat_cum.prop[, ncol(dat_cum.prop)]
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
  dat_cum.prop[is.na(dat_cum.prop)] <- 0
} 

head(dat_cum.prop); tail(dat_cum.prop); dim(dat_cum.prop)

str(dat_cum.prop)
hist(dat_cum.prop$Proportion_germ[dat_cum.prop$Species == "littoralis"])
hist(dat_cum.prop$Proportion_germ[dat_cum.prop$Species == "torulosa"])


summary(dat_cum.prop$Proportion_germ)
#View(dat_cum.prop[, c(1, 72:ncol(dat_cum.prop))])


# For each analysis we want to include a random effect for individual and replicate (set) accounting for any differences between individuals and germination between replicates. 



# 5. Add extra information columns for analyses ----
# Add information on seed weight, and fire frequency
seed_charact <- read.csv('./00_Data/Seeds_data/Seed_characteristics_sample_information.csv', stringsAsFactors = T, header = T)
head(seed_charact)
seed_charact$seed_weight <- (seed_charact$Seed_weight_subset_of_10_seeds/10) # Calculate weight per seed
seed_charact_s <- seed_charact[, c(5,9,23,24)]
head(seed_charact_s)
seed_charact_s$seed_wt_mg <- seed_charact_s$seed_weight*1000


# 5.1 Load fire data ----
# We need to make decisions about what fire data to use, whether we need to keep QPWS data where it is available and whether to use the GLM or GAM predictions. We did conclude that the GAM predictions were better for investigations of fire frequency but lets just confirm this.
pts <- vect('D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/QPWS_random.gpkg')
GLM_pred <- rast('D:/PhD/R_analysis/Fire_freq/04_Results/Prediction_rasters/GLM_pred.tif')
GAM_pred <- rast('D:/PhD/R_analysis/Fire_freq/04_Results/Prediction_rasters/GAM_pred.tif')
QPWS <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_raster.tif')



# 5.2 Extract fire frequencies to points for seed lots and transects  ----
transects <- read.csv('C:/Users/s4590925/OneDrive - The University of Queensland/Desktop/GitHub/Allocasuarina_germfire/00_Data/Transect_data/Transect_level.csv', stringsAsFactors = T, header = T)
head(transects)
# Note for transects - we measured Gatton transect GB1 first, with fewer measurements taken from each individual and going out to 100m. As all other transects were reduced, we have limited data for this transect to 50m.

transect_sf <- st_as_sf(transects, coords = c('Longitude', 'Latitude'), crs = 'EPGS:4326')
head(transect_sf)
transect_sf
transect_sf <- sf::st_set_crs(transect_sf, 'EPSG:4326')
transect_v <- vect(transect_sf) %>% 
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

GLM_t <- extract(GLM_pred, transect_v)
colnames(GLM_t) <- c("ID", "Fire_freq")
GLM_t$Data <- "GLM"
GLM_t$transect <- transect_v$Transect
GLM_t$Species <- transect_v$Species


GAM_g <- extract(GAM_pred, germ)
colnames(GAM_g) <- c("ID", "Fire_freq")
GAM_g$Data <- "GAM"
GAM_g$transect <- germ$Transect
GAM_g$Species <- germ$Species

GAM_t <- extract(GAM_pred, transect_v)
colnames(GAM_t) <- c("ID", "Fire_freq")
GAM_t$Data <- "GAM"
GAM_t$transect <- transect_v$Transect
GAM_t$Species <- transect_v$Species

QPWS_g <- extract(QPWS, germ)
colnames(QPWS_g) <- c("ID", "Fire_freq")
QPWS_g[is.na(QPWS_g)] <- 0
QPWS_g$Data <- "QPWS"
QPWS_g$transect <- germ$Transect
QPWS_g$Species <- germ$Species

QPWS_t <- extract(QPWS, transect_v)
colnames(QPWS_t) <- c("ID", "Fire_freq")
QPWS_t[is.na(QPWS_t)] <- 0
QPWS_t$Data <- "QPWS"
QPWS_t$transect <- transect_v$Transect
QPWS_t$Species <- transect_v$Species



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
  #st_transform('EPSG:3577')
#QPWS_hist <- st_crop(QPWS_hist, GAM_pred)
#rtemp <- raster::raster(xmn = 1902033, xmx = 2111776, ymn = -3257627, ymx = -2954985, res = 5, crs = 'EPSG:3577')
#QPWS_SEQ_freq_rast <- fasterize::fasterize(QPWS_hist, rtemp, field = 'OUTYEAR', fun = 'max')
#raster::writeRaster(QPWS_SEQ_freq_rast, 'D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/SEQ/QPWS_TSF.tif')


# 5.5 Extract QPWS fire history for transects and seed lots ----
QPWS_hist_rast <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/SEQ/QPWS_TSF.tif')
transect_TSF <- extract(QPWS_hist_rast, transect_v, xy = T)
germ_TSF <- extract(QPWS_hist_rast, germ, xy = T)

# Find which points have NAs to help with replacement
t_TSF_na <-  transect_TSF[is.na(transect_TSF$QPWS_TSF),]
g_TSF_na <- germ_TSF[is.na(germ_TSF$QPWS_TSF), ]

sent <- rast('D:/ADATA/Fire_data/QLD/Landsat_Annual_Fire_Scars/TERN_Sentinel2/cvmsre_qld_2023_afma2(1).tif')
t_TSF2 <- extract(sent, transect_v)
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

# 5.5 Calculate time since fire ----
transects_fire <- as.data.frame(transects)
transects_fire$Fire_freq <- QPWS_t[, 2]
transects_fire$Fire_freq <- ifelse(is.na(transects_fire$Fire_freq), GAM_t$Fire_freq, transects_fire$Fire_freq)
transects_fire$Last_fire <- transect_TSF$QPWS_TSF
transects_fire$TSF <- transects$Year_measured - transects_fire$Last_fire
head(transects_fire)
unique(is.na(transects_fire))


seed_charact_s$Fire_freq <- QPWS_g$Fire_freq
seed_charact_s$Fire_freq <- ifelse(is.na(seed_charact_s$Fire_freq), GAM_g$Fire_freq, seed_charact_s$Fire_freq)
seed_charact_s$Last_fire <- germ_TSF$QPWS_TSF
seed_charact_s$TSF <- seed_charact_s$Collection_year - seed_charact_s$Last_fire
head(seed_charact_s)
unique(is.na(seed_charact_s))




# 5.6 Produce histograms to investigate time since fire ----
tor_tran_TSF <- hist(transects_fire$TSF[transects_fire$Species == "torulosa"], breaks = seq(-1,37,1), main = "torulosa", las = 1, xlab = "Fire_frequency")
lit_tran_TSF <- hist(transects_fire$TSF[transects_fire$Species == "littoralis"], breaks = seq(-1,37,1), main= "littoralis", las = 1, xlab = "Fire frequency")

tor_germ_TSF <- hist(dat_cum.prop$TSF[dat_cum.prop$Species == "torulosa"], breaks = seq(-1,37,1), main = "torulosa", las = 1, xlab = "Fire frequency")
lit_germ_TSF <- hist(dat_cum.prop$TSF[dat_cum.prop$Species == "littoralis"], breaks = seq(-1,37,1), main = "littoralis", las = 1, xlab = "Fire frequency")




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

hist(germ_lit$Proportion_germ)
hist(germ_tor$Proportion_germ)



# 6. Calculate proportion of seedlings to adults, saplings to adults, and seedlings+saplings to adults for transect data -----
transects_fire$Proportion_seedlings <- transects_fire$Number_seedlings/(transects_fire$Number_seedlings + transects_fire$Number_mature)
transects_fire$Proportion_saplings <- transects_fire$Number_saplings/(transects_fire$Number_saplings + transects_fire$Number_mature)
transects_fire$Proportion_recruits <- (transects_fire$Number_seedlings + transects_fire$Number_saplings)/(transects_fire$Number_seedlings + transects_fire$Number_saplings + transects_fire$Number_mature)

head(transects_fire)

# 7. Read in tree level transect data ----
treedat <- read.csv('./00_Data/Transect_data/Transect_tree_level.csv', header = T, stringsAsFactors = T)
head(treedat);dim(treedat)
str(treedat)
# Note for tree level transect data - we measured Gatton transect GB1 first, with fewer measurements taken from each individual and going out to 100m. As all other transects were reduced, we have limited data for this transect to 50m. However, we still lack a lot of data for the Gatton transect which will limit its inclusion in some analyses.

# 8. Convert tree level transect data to spatial ----
treedat_sf <- st_as_sf(treedat, coords = c('Longitude', 'Latitude'), crs = 'EPSG:4326')
treedat_v <- vect(treedat_sf) %>% 
  project('EPSG:3577')
plet(treedat_v)


# 9. Extract fire information -----
# Extract fire frequency, last year burnt and then time since last fire. As we move along the transect fire frequency may change, if it does then we will need to extract new TSF information

tree_ff <- extract(QPWS, treedat_v)
tree_ff$transect <- treedat_v$Transect
table(tree_ff$transect, tree_ff$QPWS_SEQ_freq_raster) # A couple of locations have multiple fire frequencies. 

tree_sent <- extract(GAM_pred, treedat_v)

treedat_v$Fire_freq <- tree_ff$QPWS_SEQ_freq_raster
treedat_v$Fire_freq <- ifelse(is.na(treedat_v$Fire_freq), round(tree_sent$lyr1), treedat_v$Fire_freq)
unique(treedat_v$Fire_freq)

tree_burn <- extract(QPWS_hist_rast, treedat_v)
is.na(tree_burn)

sent <- rast('D:/ADATA/Fire_data/QLD/Landsat_Annual_Fire_Scars/QLD_spatial_Landsat/1987/IMG_QLD_LANDSAT_FIRESCARS_1987.tif')
tree_burn2 <- extract(sent, treedat_v) 
  colnames(tree_burn2) <- c("ID", "sent")

tree_burn2[tree_burn2$sent !=0, ]

#tree_burn$QPWS_TSF <- ifelse(is.na(tree_burn$QPWS_TSF), tree_burn2$sent, tree_burn$QPWS_TSF)
unique(is.na(tree_burn))
unique(tree_burn$QPWS_TSF)
tree_burn[c(889:993, 909:914), 2] <- 2011
tree_burn[c(433:490, 761:765, 770:773), 2] <- 2004
tree_burn[716:726, 2] <- 2003
tree_burn[727:760, 2] <- 1999


treedat_v$Last_burn <- tree_burn$QPWS_TSF
unique(is.na(treedat_v$Last_burn)) # Still have some NAs
treedat_v$Last_burn[is.na(treedat_v$Last_burn)] <- 1986
unique(is.na(treedat_v$Last_burn)) # Still have some NAs
treedat_v$TSF <- treedat_v$Year_measured - treedat_v$Last_burn

treedat$Fire_freq <- treedat_v$Fire_freq
treedat$Last_burn <- treedat_v$Last_burn
treedat$TSF <- treedat$Year_measured - treedat$Last_burn
head(treedat);dim(treedat)
  




# 10. Extract environmental information -----
e <- ext(1936841, 2086019, -3241253, -3105183)


# 10.1 Foliage projective cover ----
# We want information relating to foliage projective cover and precipitation.
FPC14 <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/FPC/DP_QLD_FPC2014.tif')
unique(FPC14$DP_QLD_FPC2014)
# Need to do some more adjustments to this data - metadata states that data ranges between 100-200 which is equivalent to 0-100% FPC. values erroneously predicted above 100% or below 0% have been classed as above 200 and below 100 respectively. Zero values indicate NULL data. The data actually seems to be ranging between 88-213. Post 2014, values range between 0-100 which would denote the % cover without any further changes being required. Let's take a look at the data in ArcGIS as well to make sure this is true for the 2014 dataset.


# Create matrices for reclassification
A = matrix(
  c(88:99, 201:213),
  nrow = 25,
  ncol = 2)
A[,2] <- 0

B = matrix(
  c(100:200),
  nrow = 101,
  ncol = 1
)
B <- cbind(B, 0:100)


reclas <- rbind(A, B)

# Now reclassify FPC14
FPC14r <- classify(FPC14, rcl = reclas)
FPC14r # Check how this looks
plot(FPC14r)

FPC14_seq <- crop(FPC14r, e)
FPC14seq <- project(FPC14_seq, 'EPSG:3577')
FPC14seq
plot(FPC14seq)

writeRaster(FPC14seq, './00_Data/Spatial_data/Environmental_data/FPC/FPC14seq.tif')


# 1.3.1 Add in the data from more recent years post 2014 ----
# Firstly, look at the new data to see what needs to be changed
FPC18 <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/FPC/DP_QLD_S2_WOODY_FPC_2018.tif')

# Need to aggregate the data to a coarser resolution, from 10m to 30m, and then we also need to crop the data to SEQ
FPC18_seq <- crop(FPC18, e)

FPC18seq <- terra::aggregate(FPC18_seq, fact = 3)
FPC18seq # Check how this looks
plot(FPC18seq)

writeRaster(FPC18seq, './00_Data/Spatial_data/Environmental_data/FPC/FPC18_SEQ.tif')


FPC19 <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/FPC/DP_QLD_S2_WOODY_FPC_2019.tif')
FPC19_seq <- crop(FPC19, e)
FPC19seq <- terra::aggregate(FPC19_seq, fact = 3)
FPC19seq
plot(FPC19seq)
writeRaster(FPC19seq, './00_Data/Spatial_data/Environmental_data/FPC/FPC19_SEQ.tif')


FPC20 <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/FPC/DP_QLD_S2_WOODY_FPC_2020.tif')
FPC20_seq <- crop(FPC20, e)
FPC20seq <- terra::aggregate(FPC20_seq, fact = 3)
FPC20seq
plot(FPC20seq)
writeRaster(FPC20seq, './00_Data/Spatial_data/Environmental_data/FPC/FPC20_SEQ.tif')


FPC21 <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/FPC/DP_QLD_S2_FPC_2021.tif') # The coordinate reference system has not been read in the same manner as the others so we will need to fix this

FPC21_seq <- crop(FPC21, e)
FPC21seq <- project(FPC21_seq, 'EPSG:3577')
FPC21seq <- terra::aggregate(FPC21_seq, fact = 3)
FPC21seq
plot(FPC21seq)
writeRaster(FPC21seq, './00_Data/Spatial_data/Environmental_data/FPC/FPC21_SEQ.tif')


# 1.3.2 Combine the FPC data into one raster ----
FPC14seq <- resample(FPC14seq, FPC18seq) # Need the extents to match
FPC <- terra::mean(FPC14seq, FPC18seq, FPC19seq, FPC20seq, FPC21seq)
FPC
plot(FPC)


writeRaster(FPC, './00_Data/Spatial_data/Environmental_data/FPC/FPC_all.tif')
FPC <- rast('./00_Data/Spatial_data/Environmental_data/FPC/FPC_all.tif')

tree_FPC <- extract(FPC, treedat_v)
tran_FPC <- extract(FPC, transect_v)
seed_FPC <- extract(FPC, germ)

treedat_v$FPC <- tree_FPC$DP_QLD_FPC2014
transect_v$FPC <- tran_FPC$DP_QLD_FPC2014

treedat$FPC <- tree_FPC$DP_QLD_FPC2014
transects_fire$FPC <- tran_FPC$DP_QLD_FPC2014
seed_charact_s$FPC <- seed_FPC$DP_QLD_FPC2014

# 10.2 Precipitation seasonality ----
# Precipitation seasonality captures the extremes which may be more what we want than average precipitation. 

gdalwarp(srcfile = 'D:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/BioClim/wc2.1_30s_bio_15.tif',
         dstfile = './00_Data/Spatial_data/Environmental_data/BioClim/precipseason.tif',
         t_srs = 'EPSG:3577')

# Need to crop and change resolution
precipseason <- rast('./00_Data/Spatial_data/Environmental_data/BioClim/precipseason.tif')
precipseason
precip <- crop(precipseason, e)
precipr <- disagg(precip, fact = 30)

writeRaster(precipr, '00_Data/Spatial_data/Environmental_data/BioClim/precipseason_SEQ.tif')
precipr <- rast('./00_Data/Spatial_data/Environmental_data/BioClim/precipseason_SEQ.tif')

tree_precip <- extract(precipr, treedat_v)
tran_precip <- extract(precipr, transect_v)
seed_precip <- extract(precipr, germ)
treedat_v$Precip <- tree_precip$precipseason
transect_v$Precip <- tran_precip$precipseason

treedat$Precip <- tree_precip$precipseason
transects_fire$Precip <- tran_precip$precipseason
seed_charact_s$Precip <- seed_precip$precipseason



# 10.3 Temperature seasonality ---- 
# As above want to capture the affect of extremes better.

gdalwarp(srcfile = 'D:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/BioClim/wc2.1_30s_bio_4.tif',
         dstfile = './00_Data/Spatial_data/Environmental_data/BioClim/tempseason.tif',
         t_srs = 'EPSG:3577')

tempseason <- rast('./00_Data/Spatial_data/Environmental_data/BioClim/tempseason.tif')
temp <- crop(tempseason, e)
tempr <- disagg(temp, fact = 30)

writeRaster(tempr, './00_Data/Spatial_data/Environmental_data/BioClim/tempseason_SEQ.tif')
tempr <- rast('./00_Data/Spatial_data/Environmental_data/BioClim/tempseason_SEQ.tif')

tree_temp <- extract(tempr, treedat_v)
tran_temp <- extract(tempr, transect_v)
seed_temp <- extract(tempr, germ)

treedat$Temp <- tree_temp$tempseason
transects_fire$Temp <- tran_temp$tempseason
seed_charact_s$Temp <- seed_temp$tempseason


head(treedat_v)
head(treedat)
head(transects_fire)
head(transect_v)
head(seed_charact_s)





save.image('./02_Workspaces/001_data_wrangling.RData')

write.csv(treedat, './00_Data/Transect_data/Tree_level_enviro.csv')
write.csv(transects_fire, './00_Data/Transect_data/Transect_level_enviro.csv')

dat_cum.prop <- left_join(dat_cum.prop, seed_charact_s, by = "Individual")
head(dat_cum.prop)
unique(dat_cum.prop$Fire_freq)

write.csv(dat_cum.prop, './00_Data/Full_experiment/Full_experiment_cumulative_germ.csv')
