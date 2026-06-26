
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
library(cowplot)
library(ggnewscale)
library(gridGraphics)
sf::sf_use_s2(FALSE)

# 2. Read and format sample data ----
sdat <- read.csv('./00_Data/Transect_data/Transect_location_data.csv', header  = T, stringsAsFactors = T)
head(sdat)
sdat <- sdat[-6,]
sdat <- droplevels(sdat)
str(sdat)



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
# Create new extent for the study area
e <- ext(1936841, 2086019, -3241253, -3105183)



protected_areas <- vect('/Volumes/Extreme SSD/ADATA/QLD/QLD_Protected_areas/QSC_Extracted_Data_20220221_143858585000-64332/Protected_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(e)

# Download nature refuge files from internet and unzip
Bulimbah <- download.file("https://wetlandinfo.des.qld.gov.au/resources/wetland-summary/area/nature-refuge/kml/nature-refuge-bulimbah-nature-refuge.kmz", destfile = './00_Data/Spatial_data/Nature_refuges/Bulimbah_nature_refuge.kmz', mode = "wb", cacheOK = F)
unzip(zipfile = './00_Data/Spatial_data/Nature_refuges/Bulimbah_nature_refuge.kmz', exdir = './00_Data/Spatial_data/Nature_refuges')

Gillies <- download.file("https://wetlandinfo.des.qld.gov.au/resources/wetland-summary/area/nature-refuge/kml/nature-refuge-gillies-ridge-nature-refuge.kmz", destfile = './00_Data/Spatial_data/Nature_refuges/Gillies_nature_refuge.kmz', mode = "wb", cacheOK = F)
unzip(zipfile = './00_Data/Spatial_data/Nature_refuges/Gillies_nature_refuge.kmz', exdir = './00_Data/Spatial_data/Nature_refuges')

Bartopia <- download.file("https://wetlandinfo.des.qld.gov.au/resources/wetland-summary/area/nature-refuge/kml/nature-refuge-bartopia-nature-refuge.kmz", destfile = './00_Data/Spatial_data/Nature_refuges/Bartopia.kmz', mode = "wb", cacheOK = F)
unzip(zipfile = './00_Data/Spatial_data/Nature_refuges/Bartopia.kmz', exdir = './00_Data/Spatial_data/Nature_refuges/')

unzip(zipfile = './00_Data/Spatial_data/Nature_refuges/Entire_OHV_perimeter.kmz', exdir='./00_Data/Spatial_data/Nature_refuges/') # Note: this needs to be renamed

# Read in nature refuge files
HV <- vect('./00_Data/Spatial_data/Nature_refuges/doc.kml') %>% 
  project('EPSG:3577')
Gillies <- vect('./00_Data/Spatial_data/Nature_refuges/nature-refuge-gillies-ridge-nature-refuge.kml') %>% 
  project('EPSG:3577')
Bulimbah <- vect('./00_Data/Spatial_data/Nature_refuges/nature-refuge-bulimbah-nature-refuge.kml') %>% 
  project('EPSG:3577')
Bartopia <- vect('./00_Data/Spatial_data/Nature_refuges/nature-refuge-bartopia-nature-refuge.kml') %>% 
  project('EPSG:3577')

# 4. Read in fire frequency data for QPWS estates  ----
QPWS_fire <- rast('/Volumes/Extreme SSD/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped_reproj.tif') %>% 
  crop(e)
plet(QPWS_fire)

# 5. Read in modelled satellite fire frequency for outside QPWS estates ----
# This data was produced as part of a different project, refer to https://github.com/felicityeloise/Fire_freq_pred_modelling.git
mod_satellite_fire <- rast('/Volumes/Extreme SSD/PhD/R_analysis/Fire_freq/04_Results/Prediction_rasters/GAM_pred.tif') %>% 
  crop(e)
plet(mod_satellite_fire)


mod_satellite_fire <- mask(mod_satellite_fire, QPWS_fire, inverse = T)
plet(mod_satellite_fire)


# 6. Read in the data required for the inset ----
Aus <- download.file("https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/STE_2021_AUST_SHP_GDA2020.zip", destfile = './00_Data/Spatial_data/Australia.zip', mode = "wb", cacheOK = F)
unzip(zipfile = './00_Data/Spatial_data/Australia.zip', exdir = './00_Data/Spatial_data/Australia')
Australia <- vect('./00_Data/Spatial_data/Australia/STE_2021_AUST_GDA2020.shp') %>%
  project("EPSG:3577")
QLD <- subset(Australia, Australia$STE_NAME21 == "Queensland")
coast <- crop(QLD, e)


# 7. Produce map ----

p1 <- 
ggplot() + 
  geom_spatvector(data = coast, col = "black", fill = 'transparent')+
  geom_spatraster(data = QPWS_fire) +
  scale_fill_continuous(na.value = "transparent", limits = c(0,7), breaks = c(0,5,7), low = "#FFF5F0", high = "darkred", name = expression(bold("Fire frequency")))+
  geom_spatraster(data = mod_satellite_fire) +
  ggnewscale::new_scale_fill() +
  theme_cowplot(font_size = 26) +
  labs(fill = expression(bold("Fire frequency"))) +
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.25, "cm"), pad_y = unit(0.2, "cm"),pad_x = unit(23, 'cm'), style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", pad_x = unit(14, "cm"), text_cex = 1.75)+
  theme(legend.text = element_text(size = rel(1.2)),
        legend.key.height = unit(1.5, "cm"),
        legend.title = element_text(size = rel(1.4)),
        legend.justification = "bottom")+
  geom_spatvector(data = protected_areas, fill = NA, color = "black") +
  geom_spatvector(data = Gillies, fill = NA, color = 'black') +
  geom_spatvector(data = HV, fill = NA, color = 'black')+
  geom_spatvector(data = Bartopia, fill = NA, color = 'black')+
  geom_spatvector(data = Bulimbah, fill = NA, color = 'black')+
  geom_spatvector(data = sdat_r, aes(fill = Species, shape = Species), size = 4) +
  scale_shape_manual(values = c(21, 23), labels = c(expression(italic("Allocasuarina littoralis")), expression(italic("Allocasuarina torulosa"))), name = expression(bold("Species"))) +
  scale_fill_manual(values = c("#7B68EE", "#1E90FF"), labels = c(expression(italic("Allocasuarina littoralis")), expression(italic("Allocasuarina torulosa"))), name = expression(bold("Species"))) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 23), size = 6, alpha = 1))) +
  theme(legend.text = element_text(size = rel(1.2)),
        legend.key.height = unit(0.8, "cm"),
        legend.title = element_text(size = rel(1)),
        legend.position = "right",
        legend.justification = "center")


# 8. Create the inset map
inset <-
  ggplot()+
  geom_spatvector(data = Australia, fill = NA)+
  geom_spatraster(data = QPWS_fire) +
  geom_spatraster(data = mod_satellite_fire) +
  scale_fill_continuous(na.value = "transparent", breaks = c(1,13),low = "darkred", high = "darkred")+
  theme_void() +
  theme(legend.position = "none") +
  geom_rect(aes(xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4]), alpha = 0, colour = "black", linewidth = 1.5)+
  annotation_scale(location = "bl", text_cex = 1.3, pad_y = unit(0.05, "cm"), pad_x = unit(5, 'cm'))


sample_plot <- ggdraw()+
  draw_plot(p1)+
  draw_plot(inset, x = 0.55, y = 0.7, width = 0.4, height = 0.3) 
sample_plot
ggsave("./03_Results/Sample_map.pdf", sample_plot, width = 19, height = 10)



# Get information on the known distributions of A. torulosa and A. littoralis known occurrences in relation to fire

lit <- read.csv('/Volumes/Extreme SSD/ADATA/Observation data/ALA_4:3:2022/Allocasuarina_littoralis_records-2022-03-04/Allocasuarina_littoralis_records-2022-03-04.csv')
head(lit)
lit <- lit[!is.na(lit$Decimal.latitude..WGS84.),]
lit <- lit[!is.na(lit$Decimal.longitude..WGS84.),]
unique(is.na(lit))
lit_sf <- st_as_sf(lit, coords = c('Decimal.longitude..WGS84.', 'Decimal.latitude..WGS84.'), crs = 'EPSG:4326')
lit_v <- vect(lit_sf) %>% 
  project('EPSG:3577')

tor <- read.csv('/Volumes/Extreme SSD/ADATA/Observation data/ALA_4:3:2022/Allocasuarina_torulosa_records-2022-03-04/Allocasuarina_torulosa_records-2022-03-04.csv')
tor <- tor[!is.na(tor$Decimal.latitude..WGS84.),]
tor <- tor[!is.na(tor$Decimal.longitude..WGS84.),]
tor_sf <- st_as_sf(tor, coords = c('Decimal.longitude..WGS84.', 'Decimal.latitude..WGS84.'), crs = 'EPSG:4326')
tor_v <- vect(tor_sf) %>% 
  project('EPSG:3577')

# Extract fire frequency information for each Allocasuarina
fire <- c(QPWS_fire, mod_satellite_fire)
fire
lit_fire <- extract(fire, lit_v)
tor_fire <- extract(fire, tor_v)

lit_v$ff <- lit_fire$QPWS_SEQ_freq_raster
lit_v$ff <- ifelse(is.na(lit_fire$QPWS_SEQ_freq_raster), lit_fire$lyr1, lit_v$ff)
unique(lit_v$ff)

tor_v$ff <- tor_fire$QPWS_SEQ_freq_raster
tor_v$ff <- ifelse(is.na(tor_fire$QPWS_SEQ_freq_raster), tor_fire$lyr1, tor_v$ff)

lit_v <- crop(lit_v, e)
plet(lit_v)
round(unique(lit_v$ff))

tor_v <- crop(tor_v, e)

lit_c <- col2rgb("#483D8B", alpha = T)
lit_c
tor_c <- col2rgb("steelblue4", alpha = T)
tor_c



dev.new()
hist(round(tor_v$ff), breaks = seq(0,7,1), ylim = c(0,140), las = 1, col = rgb(54/255, 100/255, 139/255, 0.7), border = 'white', xlim = c(0,7), xaxt = "n", yaxt = "n", main = "", xlab = "", ylab = "")
hist(round(lit_v$ff), breaks = seq(0,7,1), ylim = c(0, 140), las = 1, col = rgb(72/255, 61/255, 139/255, 0.7), border = 'white', xaxt = "n", yaxt = "n", add = T)
axis(side = 1, at = seq(0,6,1), line = -1.1, tick = F, hadj = -3, cex = 1)
axis(side = 2, at = seq(0,140,20), las = 1, line = -0.7, cex = 1)
axis(side = 1, at = seq(0, 7,1), labels = F, line = -0.3)
mtext(expression(bold("Fire frequency")), side = 1, line = 1.5, cex = 1.5)
mtext(expression(bold("Number of occurrences")), side = 2, cex = 1.5, line = 1.6)

p2 <- recordPlot()
dev.off()
p2g <- as_grob(p2)

dev.new(width = 20, height = 15)
sample_plot <- ggdraw()+
  draw_plot(p1)+
  draw_plot(inset, x = 0.55, y = 0.7, width = 0.4, height = 0.3) +
  draw_plot(p2g, x = 0.55, y = 0.6, width = 0.4, height = 0.3)
sample_plot
ggsave("./03_Results/Sample_map_hist.pdf", sample_plot, width = 19, height = 10)
ggsave('./03_Results/Sample_map.png', sample_plot, width = 19, height = 10)

