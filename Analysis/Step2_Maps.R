# Maps------------

## Libraries -----
library(ggplot2)
library(sf)
library(geosphere)
library(viridis)
library(tidyverse)
library(lwgeom)
library(ggspatial)



## Shape files ----------------- 
# These files should be located on your own computer. Replace directories/files as needed
### LML-----
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/LML_shape/")
LML_shape <- st_read("World_Lakes.shp") 

target_crs <- st_crs(LML_shape)

### FBL ----------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/FBL_shape/")
FBL_shape <- st_read("World_Lakes.shp") %>%
  st_transform(target_crs)
### ETL -----------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/kml/")
east_shape = st_read("east.kml")%>%
  st_transform(target_crs) %>%
  st_zm(drop = TRUE, what = "ZM")
### GNL ----------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files//ALC/GNL/")
GNL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### PRL ------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/PRL/")
PRL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### POL ------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/POL/")
POL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### CAL ----------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/CAL/")
CAL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### FBL----------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/FBL/")
FBL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### SBL ---------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/SBL/")
SBL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### TBL ---------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/TBL/")
TBL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### FOB --------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/FOB/")
FOB_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### CSL ------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/CSL/")
CSL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### USP ----
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/USP/")
USP_shape = st_read("Upper Sylvan.kml")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### LSP ----
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/LSP/")
LSP_shape = st_read("Lower Sylvan.kml")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### COM ----- 
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/COM/")
COM_shape = st_read("Combs.kml")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### TRP----- 
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/TRP/")
TRP_shape = st_read("taylor.kml")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### MNP ----- 
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/MNP/")
MNP_shape = st_read("Mountain.kml")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### PEP ----- 
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/PEP/")
PEP_shape = st_read("Pinchnose.kml")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### WLL ----- 
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/WLL/")
WLL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### SDL ----- 
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/SDL/")
SDL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### RKP ----- 
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/RKP/")
RKP_shape = st_read("Rock Pond.kml")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### HAL ----- 
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/HAL/")
HAL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### GEL ----- 
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/GEL/")
GEL_shape = st_read("GEL.kml")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### ORL ----- 
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/ORL/")
ORL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")
### JSL ----- 
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP/shape_files/ALC/JEL/")
JSL_shape = st_read("World_Lakes.shp")%>%
  st_transform(target_crs) %>%
  st_zm( drop = TRUE, what = "ZM")


## ALC Map
ACL_Map = ggplot() +
  theme_minimal(base_size = 12) +
  # Add North arrow
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
                         style = north_arrow_fancy_orienteering) +
  
  geom_sf(data = east_shape, fill = "#154c79") +
  geom_sf(data = PRL_shape, fill ="#154c79" ) +
  geom_sf(data = GNL_shape, fill = "#154c79" ) + 
  geom_sf(data = LML_shape, fill = "#154c79") +
  geom_sf(data = POL_shape, fill = "#154c79" ) + 
  geom_sf(data = CAL_shape, fill = "#154c79" ) + 
  geom_sf(data = FBL_shape, fill = "#154c79") + 
  geom_sf(data = SBL_shape, fill = NA ) + 
  geom_sf(data = TBL_shape, fill = NA ) + 
  geom_sf(data = FOB_shape, fill = "#154c79" ) + 
  geom_sf(data = CSL_shape, fill = NA ) + 
  geom_sf(data = MNP_shape, fill = "#154c79" ) +
  geom_sf(data = TRP_shape, fill = NA ) +
  geom_sf(data = PEP_shape, fill = "#154c79" ) +
  geom_sf(data = COM_shape, fill = NA ) +
  geom_sf(data = USP_shape, fill = NA ) + 
  geom_sf(data = LSP_shape, fill = "#154c79" ) +
  geom_sf(data = WLL_shape, fill = NA ) + 
  geom_sf(data = SDL_shape, fill = NA ) + 
  geom_sf(data = RKP_shape, fill = "#154c79" ) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/thermal-guild-otolith")
ggsave("Graphics/ALC_map.pdf", width = 5, height = 4)

## Upload TDo5 Data

afrp_lakes = read.csv("Data/TempDO Data/lake_metadata.csv") %>% 
  group_by(afrp_abbrev) %>%
  slice(1)

all_lakes = bind_rows( 
  #HAL_shape   %>% mutate(lake = "HAL"),
  east_shape %>% mutate(lake = "ETL"),
  PRL_shape   %>% mutate(lake = "PRL"),
  GNL_shape   %>% mutate(lake = "GNL"),
  LML_shape   %>% mutate(lake = "LML"),
  POL_shape   %>% mutate(lake = "POL"),
  CAL_shape   %>% mutate(lake = "CAL"),
  FBL_shape   %>% mutate(lake = "FBL"),
  SBL_shape   %>% mutate(lake = "SBL"),
  TBL_shape   %>% mutate(lake = "TBL"),
  FOB_shape   %>% mutate(lake = "FOB"),
  CSL_shape   %>% mutate(lake = "CSL"),
  MNP_shape   %>% mutate(lake = "MNP"),
  TRP_shape   %>% mutate(lake = "TRP"),
  PEP_shape   %>% mutate(lake = "PEP"),
  COM_shape   %>% mutate(lake = "COM"),
  USP_shape   %>% mutate(lake = "USP"),
  LSP_shape   %>% mutate(lake = "LSP"),
  WLL_shape   %>% mutate(lake = "WLL"),
  SDL_shape   %>% mutate(lake = "SDL"),
  RKP_shape   %>% mutate(lake = "RKP")
) %>% 
  left_join(afrp_lakes, by = c("lake" = "afrp_abbrev"))

tdo5_alc = ggplot(data = all_lakes, aes(fill = tdo5_avg)) + 
  geom_sf()+
  theme_minimal(base_size = 12) +
  # Add North arrow
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_viridis_c() + 
  labs(fill = "TDO5")

ggsave(tdo5_alc, file = "Graphics/tdo5_alc.png", width = 5, height = 5, dpi = 500)

ggplot() +
  theme_minimal(base_size = 12) +
  # Add North arrow
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
                         style = north_arrow_fancy_orienteering) +
  
  geom_sf(data = east_shape, fill = "#154c79") +
  geom_sf(data = PRL_shape, fill ="#154c79" ) +
  geom_sf(data = GNL_shape, fill = "#154c79" ) + 
  geom_sf(data = LML_shape, fill = "#154c79") +
  geom_sf(data = POL_shape, fill = "#154c79" ) + 
  geom_sf(data = CAL_shape, fill = "#154c79" ) + 
  geom_sf(data = FBL_shape, fill = "#154c79") + 
  geom_sf(data = SBL_shape, fill = NA ) + 
  geom_sf(data = TBL_shape, fill = NA ) + 
  geom_sf(data = FOB_shape, fill = "#154c79" ) + 
  geom_sf(data = CSL_shape, fill = NA ) + 
  geom_sf(data = MNP_shape, fill = "#154c79" ) +
  geom_sf(data = TRP_shape, fill = NA ) +
  geom_sf(data = PEP_shape, fill = "#154c79" ) +
  geom_sf(data = COM_shape, fill = NA ) +
  geom_sf(data = USP_shape, fill = NA ) + 
  geom_sf(data = LSP_shape, fill = "#154c79" ) +
  geom_sf(data = WLL_shape, fill = NA ) + 
  geom_sf(data = SDL_shape, fill = NA ) + 
  geom_sf(data = RKP_shape, fill = "#154c79" ) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))

ggplot()+
  theme_minimal(base_size = 10) +
  # Add North arrow
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
                         style = north_arrow_fancy_orienteering) + 
  geom_sf(data = LML_shape, fill = "#154c79") -> LML.map
ggsave(LML.map, file = "Graphics/LML.map.png", width = 6, height = 5, dpi = 500)
