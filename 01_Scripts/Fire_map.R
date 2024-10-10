@ -0,0 +1,145 @@
  
  
  # Mape of fire frequency and sample locations for germination experiment
  # Written by Felicity Charles
  # Caveat emptor
  # Date: 10th October 2024
  
  # 1. Load packages -----
library(sf)
library(terra)
library(RColorBrewer)
library(dplyr)
library(leaflet)
library(ggplot2)
library(tidyterra)
library(ggspatial)

sf::sf_use_s2(FALSE)

# 2. Read and format sample data ----
sdat <- read.csv('./00_Data/Sample_information.csv', header  = T, stringsAsFactors = T)
head(sdat)
sdat <- na.omit(sdat) # Remove any samples with NA as they were not included in the experiment
dim(sdat) # 114 rows
sdat <- sdat[-which(sdat$Number_of_seeds_estimate_from_weight <60),]
dim(sdat)
sdat <- sdat[, c(1:5,9,10)]
head(sdat); dim(sdat)

# Convert this to a spatial features dataframe so we can plot this on a map
sdat_sf <- st_as_sf(sdat, coords = c('Longitude', 'Latitude'), crs = 'EPSG:4326')
head(sdat_sf)

# Convert to a SpatVetor and change the projection
sdat_r <- vect(sdat_sf) %>% 
  project('EPSG:3577') 

# Check how this looks
sdat_r
plet(sdat_r)

# 3. Read and format environmental spatial data ----
# Before reading the data in, create an extent object that captures the study area
e <- ext(1935620, 2021296, -3241847, -3126599)
protected_areas <- vect('E:/ADATA/QLD/QLD_Protected_areas/QSC_Extracted_Data_20220221_143858585000-64332/Protected_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(e)



# 4. Calculate fire frequency for QPWS estates ----
# These steps have already been run in a different project. Code for producing this file is copied below but we will just read in the final output.
#QPWS_fire_hist <- vect('E:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/QPWS_fire_history/Fire_history___QPWS.shp') %>% 
#project('EPSG:3577')
# Subset the fire history data 
# Keep only the data from 1987 - 2023
#QPWS_fire_hist_1987 <- subset(QPWS_fire_hist, QPWS_fire_hist$OUTYEAR >= 1987 & QPWS_fire_hist$OUTYEAR <2024)
#unique(QPWS_fire_hist_1987$OUTYEAR) # Check that this has worked

# Convert to raster and calculate fire frequency
# Fasterize only works with data in sf format so convert the data frame
#QPWS_fire <- st_as_sf(QPWS_fire_hist_1987)
#rtemp <- raster::raster(xmn = 1902033, xmx = 2111776, ymn = -3257627, ymx = -2954985, res = 5, crs = 'EPSG:3577')
#QPWS_SEQ_freq_rast <- fasterize::fasterize(QPWS_fire, rtemp, field = 'OUTYEAR', fun = 'count')


# Before reading the data in, create an extent object that captures the study area
QPWS_fire <- rast('E:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_raster.tif') %>% 
  crop(e)
plet(QPWS_fire)




# 5. Categorise the individuals by their fire frequency
sdat$ID <- 1:nrow(sdat)

# Extract the fire frequency for each sample
samp_fire <- extract(QPWS_fire, sdat_r)
sdat2 <- merge(sdat, samp_fire, by = "ID")
sdat2$fire_cat <- ifelse(sdat2$QPWS_SEQ_freq_raster >=4, "high", "low")
unique(sdat2$fire_cat)
sdat2$fire_cat <- ifelse(is.na(sdat2$fire_cat), "low", sdat2$fire_cat)
unique(sdat2$fire_cat)
head(sdat2);dim(sdat2)
sdat2$fire_cat <- paste0(sdat2$fire_cat, sdat2$Species)
sdat2$fire_cat <- ifelse(sdat2$fire_cat == "highlittoralis" | sdat2$fire_cat == "lowlittoralis", "littoralis", sdat2$fire_cat)
unique(sdat2$fire_cat)
head(sdat2); tail(sdat2); dim(sdat2)
str(sdat2)
sdat2$fire_cat <- plyr::revalue(sdat2$fire_cat, c("littoralis" = "Allocasuarina littoralis", "lowtorulosa" = "Low Allocasuarina torulosa", "hightorulosa" = "High Allocasuarina torulosa"))
head(sdat2); tail(sdat2); dim(sdat2)

# Convert this back to spatial data
sdat_r2 <- st_as_sf(sdat2, coords = c('Longitude', 'Latitude'), crs = 'EPSG:4326')
head(sdat_r2)

# Convert to a SpatVetor and change the projection
sdat_r2 <- vect(sdat_r2) %>% 
  project('EPSG:3577') 

# 6. Produce map ----

p1 <- ggplot() + 
  geom_spatraster(data = QPWS_fire) +
  theme_minimal()+
  scale_fill_continuous(na.value = "transparent", breaks = c(1,11),low = "#FFF5F0", high = "darkred")+  #Note we can change na.value to be "white", this would allow us to put any 0s as NAs and easily show these areas as white
  geom_spatvector(data = protected_areas, fill = NA, colour = "#4D60A9") +
  geom_spatvector(data = sdat_r2, aes(col = fire_cat, shape = fire_cat), size = 2)+
  scale_shape_manual(values = c(19, 17, 17)) +
  scale_color_manual(values = c("mediumblue",  "steelblue4", "dodgerblue2"))+
  labs(shape = "Species fire category", color = "Species fire category", fill = "Fire frequency")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(1, "cm"), width = unit(0.75, "cm"), pad_y = unit(0.5, "cm"), style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl")+
  theme(legend.text = element_text(size = rel(0.9)))


# Create the inset map
# Read in the data required for the inset ----
Aus <- download.file("https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/STE_2021_AUST_SHP_GDA2020.zip", destfile = './00_Data/Spatial data/Australia.zip', mode = "wb", cacheOK = F)
unzip(zipfile = './00_Data/Spatial data/Australia.zip', exdir = './00_Data/Spatial data/Australia')
QLD <- vect('./00_Data/Spatial data/Australia/STE_2021_AUST_GDA2020.shp') %>%
  project("EPSG:3577")
QLD <- subset(QLD, QLD$STE_NAME21 == "Queensland")


inset <- ggplot()+
  geom_spatvector(data = QLD, fill = NA)+
  geom_spatraster(data = QPWS_fire) +
  scale_fill_continuous(na.value = "transparent", breaks = c(1,11),low = "#FFF5F0", high = "darkred")+
  theme_void() +
  theme(legend.position = "none") +
  geom_rect(aes(xmin = 1935620, xmax = 2021296, ymin = -3241847, ymax = -3126599), alpha = 0, colour = "black")+
  annotation_scale(location = "bl")



sample_plot <- ggdraw()+
  draw_plot(p1)+
  draw_plot(inset, x = 0.7, y = 0.7, width = 0.25, height = 0.25)


ggsave("./03_Results/Sample_map.pdf", sample_plot, width = 10, height = 10)
