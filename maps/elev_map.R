# spatial packages
library(raster) #for processing some spatial data
library(rnaturalearth) #for downloading shapefiles
library(sf) #for processing shapefiles
#library("maps")
library(maps)
library(SpatialEpi)
library(concaveman) # create range edge mapping
library(geosphere) # outlier distance calculation
library(raster) # importing SDM raster

# data management packages
library(dplyr) #to keep things tidy
library(plyr)
library(devtools)
library(readxl)

# plotting packages
library(ggplot2) # for tidy plotting
library(ggpubr) # for easy selection of symbols
library(ggnewscale) # for color schemes in more complex figures 

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/")

### SCRIPT COMPONENTS
# 1. SET UP - create functions, import data and shape files
# 2. RANGE BOUNDARY = creates an estimated range edge using a concave polygon
# 3. PLOTTING - common gardens, populations (LA and CROSSES)

### 1. SET UP
# custom function for range boundary outlier removal
cluster.func <- function(x) {
  sorted_x <- sort(x)
  return(sorted_x[cluster.size])
}
cluster.size = 5 # number of neighbors within given distance

# shapefiles
world <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))

# geospatial data
pops <- read.csv("maps/data/HET_pops.csv") %>% mutate(lat.2 = if_else(pop.ID=="EGG", pop_lat-0.30, pop_lat)) # heterosis crosses
locations <- as.data.frame(read_excel("maps/data/LA_pop_garden.xlsx")) # CG and LA pops
coordinates.full <- read.csv("maps/data/iNaturalist_coords.csv", sep=',') # inaturalist pops (for range boundary)
# alt <- raster("~/Downloads/wc2.1_30s_elev.tif") # elevation -- cutting because we're using SDM predictions now
sdm.pred <- raster("SDM/data/models/rf/rf_100avg.tif")
sdm.pred.local.exp <- read.csv("SDM/data/connectivity_rf_local_lat_extended_fix1.csv")
sdm.pred.local.app <- read.csv("SDM/data/connectivity_rf_local_app_extended_fix1.csv")

# slightly move populations to reduce overlap
locations <- locations %>%
  mutate_at(vars(pop, gar.pop, env.var, env.grp, lineage), as.factor) %>% 
  mutate(lat.2 = if_else(pop=="VA5", lat-0.30, lat)) %>%
  mutate(lat.2 = if_else(garden.site == "M.L.B.S.", lat.2+0.1, lat.2)) %>%
  mutate(lat.2 = if_else(garden.site == "KENTLAND FARM", lat.2-0.1, lat.2)) %>% 
  filter(garden.site != "BLANDY")

# modify SDM tif for later plotting
df.sdm <- as.data.frame(sdm.pred, xy = TRUE)



### 2. SET UP RANGE BOUNDARY
#Import set of coordinates, with columns "ID", "latitude" and "longitude"
# dist.mat <- distm(coordinates.full %>% dplyr::select("longitude", "latitude"), fun=distHaversine)
# diag(dist.mat) <- NA # remove distances to self (dist=0)
# nearest.neighbor <- apply(dist.mat, 1, cluster.func) # grabs nth nearest neighbor from matrix
# coordinates.full$distance <- nearest.neighbor
# # remove points further than threshold
# # threshold is decided as the 5th percentile distance of the nth neighbor
# coordinates.full <- coordinates.full %>% 
#   filter(distance <= (quantile(nearest.neighbor, 0.05) %>% as.numeric)*1e3) # distance is multiplied by 1e3 to convert to meters
# remove points in a botanical garden in NY

coordinates.full <- coordinates.full %>% 
  filter(id != 49258040) %>%
  filter(id != 96591108) %>%
  filter(id != 9655474) %>%
  filter(id != 53173885) %>%
  filter(id != 51020544) %>% 
  as.data.frame()
# create hull
sf.df <- st_as_sf(coordinates.full, coords=c("longitude", "latitude"), crs=4326)
concave.hull <- concaveman(sf.df) # sets up concave polygon
buffer.hull <- st_buffer(concave.hull, dist = 5*1e3) # 50km buffer to concave polygon

# clip elevation data for easier handling -- cutting out elevation within map
max.lat <- ceiling(max(coordinates.full$latitude))
min.lat <- floor(min(coordinates.full$latitude))
max.lon <- ceiling(max(coordinates.full$longitude))
min.lon <- floor(min(coordinates.full$longitude))

sdm.pred.local.app$region <- "app"
sdm.pred.local.exp$region <- "exp"

# joining local SDM df
sdm_combo <- rbind(sdm.pred.local.app, sdm.pred.local.exp) %>% 
  dplyr::select(-X) %>% 
  mutate(x = round(x, 5), y = round(y, 4)) %>% 
  group_by(x, y) %>% 
  reframe(pred = max(pred.cell))




### 3. PLOTTING
df.sdm <- df.sdm %>% mutate(bins = cut(layer, breaks = c(-Inf,0.2,0.4,0.6,0.8,Inf))) %>% as.data.frame()
df.sdm <- df.sdm %>% mutate(bin.round = round(layer/0.2)*0.2) %>% as.data.frame()

# SDM WITH GARDENS WITH SDM PREDICTIONS 1 GRAY PER MODEL
jpeg("maps/plots/grayscale_sdm/CG_LOCATIONS_sdm_local_rangepos_color_ext_fix1.jpeg", width=5, height=5, res=500, units='in')
ggplot(data=world)+
  geom_tile(data=sdm.pred.local.exp, aes(x=x, y=y), fill="gray75", alpha=1)+ # for continuous
  # scale_fill_gradient(low = "white", high = "orchid4", na.value = NA)+ # for continuous COLOR
  # ggnewscale::new_scale_fill()+
  geom_tile(data=sdm.pred.local.app, aes(x=x, y=y), fill="gray50", alpha=1)+ # for continuous
  # scale_fill_gradient(low = "white", high = "darkgreen", na.value = NA)+ # for continuous COLOR
  # ggnewscale::new_scale_fill()+
  geom_sf(data=states, fill=NA, color="black")+ #adds state boundaries
  geom_sf(data=buffer.hull, fill=NA, color="black", linetype=3, linewidth=0.85)+ # buffer hull (used in masking)
  geom_sf(data=lakes, fill="white")+ #adds the great lakes
  coord_sf(xlim=c(min(locations$long)-1,max(locations$long)-1), ylim=c(34, 46), expand=FALSE)+
  xlab("Latitude")+
  ylab("Longitude")+
  ggnewscale::new_scale_fill()+
  geom_point(data=subset(locations, pop == "non"), aes(x=long, y=lat.2, shape=env.var, fill=env.grp), size=5)+ #adds lab collections
  # geom_point(data=subset(locations, pop == "non"), aes(x=long, y=lat.2, shape=env.var), color="gray30", size=4.3)+ #adds lab collections
  scale_shape_manual(values=c(lat=21, elev=24))+
  scale_fill_manual(values=c(low="firebrick1", mid="orchid2", high="skyblue3"))+
  theme_bw()+
  theme(text = element_text(size=15), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

# SDM WITH POPS WITH SDM PREDICTIONS 1 GRAY PER MODEL
jpeg("maps/plots/grayscale_sdm/POP_SDM_local_rangepos_color_ext_fix1.jpeg", width=5, height=5, res=500, units='in')
ggplot(data=world)+
  geom_tile(data=sdm.pred.local.exp, aes(x=x, y=y), fill="gray75", alpha=1)+ # for continuous
  # scale_fill_gradient(low = "white", high = "orchid4", na.value = NA)+ # for continuous COLOR
  # ggnewscale::new_scale_fill()+
  geom_tile(data=sdm.pred.local.app, aes(x=x, y=y, fill="gray50"), alpha=1)+ # for continuous
  # scale_fill_gradient(low = "white", high = "darkgreen", na.value = NA)+ # for continuous COLOR
  # ggnewscale::new_scale_fill()+
  geom_sf(data=states, fill=NA, color="black")+ #adds state boundaries
  geom_sf(data=buffer.hull, fill=NA, color="black", linetype=3, linewidth=0.85)+ # buffer hull (used in masking)
  geom_sf(data=lakes, fill="white")+ #adds the great lakes
  coord_sf(xlim=c(min(locations$long)-1,max(locations$long)+1), ylim=c(34, 46), expand=FALSE)+
  xlab("Latitude")+
  ylab("Longitude")+
  # geom_polygon(data=pops, aes(x=pop_long, y=lat.2, group=cross_group), color="gray70", alpha = 0)+ # FOR HET CROSS MAP
  # geom_point(data=subset(locations, pop != "non" & cg_23 == 1), aes(x=long, y=lat.2, shape=env.var), color="black", size=4)+ #adds lab collections
  geom_point(data=subset(locations, pop != "non" & cg_23 == 1), aes(x=long, y=lat.2, fill=env.grp, shape=env.var), size=5)+ #adds lab collections
  scale_shape_manual(values=c(lat=21, elev=24))+
  # scale_colour_manual(values=c(W="purple2", E="purple2", A="green4"))+
  scale_fill_manual(values=c(low="firebrick1", mid="orchid2", high="skyblue3"))+
  theme_bw()+
  theme(text = element_text(size=12), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_rect(fill="white"), axis.line = element_line(colour = "black"))
dev.off()
