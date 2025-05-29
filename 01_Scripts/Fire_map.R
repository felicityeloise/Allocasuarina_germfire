
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



protected_areas <- vect('D:/ADATA/QLD/QLD_Protected_areas/QSC_Extracted_Data_20220221_143858585000-64332/Protected_areas.shp') %>% 
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
QPWS_fire <- rast('D:/PhD/R_analysis/Fire_freq/00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped_reproj.tif') %>% 
  crop(e)
plet(QPWS_fire)

# 5. Read in modelled satellite fire frequency for outside QPWS estates ----
# This data was produced as part of a different project, refer to https://github.com/felicityeloise/Fire_freq_pred_modelling.git
mod_satellite_fire <- rast('D:/PhD/R_analysis/Fire_freq/04_Results/Prediction_rasters/GAM_pred.tif') %>% 
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
  scale_fill_continuous(na.value = "transparent", limits = c(1,13), breaks = seq(1,13,1), low = "#FFF5F0", high = "darkred")+
  geom_spatraster(data = mod_satellite_fire) +
  theme_minimal() +
  theme_cowplot(font_size = 20) +
  labs(fill = expression(bold("Fire frequency")), col = expression(bold("Species")), shape = expression(bold("Species"))) +
  annotation_north_arrow(location = "bl", which_north = T, height = unit(1.5, "cm"), width = unit(1.25, "cm"), pad_y = unit(0.2, "cm"),pad_x = unit(23, 'cm'), style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", pad_x = unit(14, "cm"), text_cex = 1.25)+
  theme(legend.text = element_text(size = rel(0.9)),
        legend.key.height = unit(1.5, "cm"),
        legend.title = element_text(size = rel(1)),
        legend.justification = "bottom")+
  geom_spatvector(data = protected_areas, fill = NA, colour = "black") +
  geom_spatvector(data = Gillies, fill = NA, colour = 'black') +
  geom_spatvector(data = HV, fill = NA, colour = 'black')+
  geom_spatvector(data = Bartopia, fill = NA, colour = 'black')+
  geom_spatvector(data = Bulimbah, fill = NA, colour = 'black')+
  geom_spatvector(data = sdat_r, aes(col = Species, shape = Species), size = 4) +
  scale_shape_manual(values = c(16,18), labels = c(expression(italic("Allocasuarina littoralis")), expression(italic("Allocasuarina torulosa")))) +
  scale_colour_manual(values = c("slateblue", "steelblue4"), labels = c(expression(italic("Allocasuarina littoralis")), expression(italic("Allocasuarina torulosa"))))+
  guides(color = guide_legend(override.aes = list(size = 4)))


    

# 8. Create the inset map
inset <-
  ggplot()+
  geom_spatvector(data = Australia, fill = NA)+
  geom_spatraster(data = QPWS_fire) +
  geom_spatraster(data = mod_satellite_fire) +
  scale_fill_continuous(na.value = "transparent", breaks = c(1,13),low = "darkred", high = "darkred")+
  theme_void() +
  theme(legend.position = "none") +
  geom_rect(aes(xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4]), alpha = 0, colour = "black")+
  annotation_scale(location = "bl", text_cex = 1, pad_y = unit(2.6, "cm"), pad_x = unit(6, 'cm'))



sample_plot <- ggdraw()+
  draw_plot(p1)+
  draw_plot(inset, x = 0.55, y = 0.5, width = 0.6, height = 0.5)
sample_plot




#ggsave("./03_Results/Sample_map.pdf", sample_plot, width = 10, height = 10)
ggsave('./03_Results/Sample_map.png', sample_plot, width = 19, height = 10)
