
# GIS-related packages
library(SpatialEpi) # for km lat/long grid
library(sp) # format conversion, general GIS
library(sf) # format conversion, general GIS
library(geodata) # Worldclim 2.1 bioclim download
library(rbioclim) # LGM data
library(geosphere) # outlier distance calculation
library(concaveman) # concave polygon for pseudo-absence generation now that rangemap is defunct
library(raster) # masking, general GIS
# library(rgeos) # for clipping/cropping
# library(maptools) # world map for simplest clipping/cropping
library(spatialEco) # pseudo-absence generation
library(terra) # bioclim rastering

# Data manipulation packages
library(dplyr) # filter and piping 
library(plyr) # round_any
library(readxl) # data import

# plotting
library(maps)
library(rnaturalearth)
library(ggplot2)

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/")
set.seed(1066) #set seed for ML random forest 

# Declare variables
grid_size = 10 # km resolution for grid cells
cluster.size = 5 # number of neighbors within given distance. 
buffer.dist.N = 5 # buffer distance for restricting presence points by phylogenetic lineage
resolution = 2.5 # worldclim resolution. finer resolution takes longer (0.5, 2.5, etc.)
pa.rep = 1 # multiplier for creating pseudo-absence data
sanity.plot = T # whether to plot sanity checks
size.t = 16 # size of text in plots
lat.cutoff = 33.5 # latitude cutoff for determining when buffer is used and when latitude is used
lat.cutoff.TF = FALSE # if you want to use the latitude cut off above. my recommendation is to set to F, as you only gain 2 remote points
region = "expansion" # "rear", "mid.lat", "leading", "expansion" (mid.lat and leading), "Appalachian"
bio.vers = "2.1" # which version of bioclim to use. 2.1 is higher resolution but does not work for LGM projection


### DATA PREPARATION

coordinates.pts <- read.csv("maps/data/iNaturalist_coords.csv", sep=',') # data from iNaturalist to May 2022

# expansion region (mid.lat and leading)
# mid-latitude (KY refugium)
coordinates.full <- coordinates.pts

# remove points from Appalachian and Smokies lineages by elevational cut-off
elev <- elevation_30s("USA", path="~/Desktop/Documents/Research/Q3/SDM/geodata/")
coordinates.full$elevation <- extract(elev, coordinates.full[,c("longitude", "latitude")])$USA_elv_msk
coordinates.full <- coordinates.full %>% filter(elevation <= quantile(coordinates.full$elevation, 0.85, na.rm=T))

# get pops close to sequenced populations of known structure
pop.lim <- read.csv("SDM/data/cam_pa2m30s_100km_expansion_1x_extended_fix1.csv")
# pop.lim <- pop.lim %>% filter(presence == 1)

# generate a buffered concave polygon around the desired samples
gen.sf <- st_as_sf(pop.lim, coords=c("longitude", "latitude"), crs=4326)
gen.hull <- concaveman(gen.sf, concavity=3)
gen.buffer <- st_buffer(gen.hull, dist = buffer.dist.N*1e3)

# apply buffer around sampling to grab points of interest
coordinates.sf <- st_as_sf(coordinates.full[,c("longitude", "latitude")], coords=c("longitude", "latitude"), crs=4326)
coordinates.buffer <- st_within(coordinates.sf, gen.buffer, sparse=F)
coordinates.full <- coordinates.sf[coordinates.buffer, ]
coordinates.full <- st_coordinates(coordinates.full) %>% as.data.frame()
names(coordinates.full) <- c("longitude", "latitude")



### appalachian
coordinates.full2 <- coordinates.pts
pop.lim2 <- read.csv("SDM/data/cam_pa2m30s_100km_Appalachian_1x_extended_fix1.csv")
# pop.lim2 <- pop.lim2 %>% filter(presence == 1)

# generate a buffered concave polygon around the desired samples
gen.sf2 <- st_as_sf(pop.lim2, coords=c("longitude", "latitude"), crs=4326)
gen.hull2 <- concaveman(gen.sf2, concavity = 2.25)
gen.buffer2 <- st_buffer(gen.hull2, dist = buffer.dist.N*1e3)

# apply buffer around sampling to grab points of interest
coordinates.sf2 <- st_as_sf(coordinates.full2[,c("longitude", "latitude")], coords=c("longitude", "latitude"), crs=4326)
coordinates.buffer2 <- st_within(coordinates.sf2, gen.buffer2, sparse=F)
coordinates.full2 <- coordinates.sf2[coordinates.buffer2, ]
coordinates.full2 <- st_coordinates(coordinates.full2) %>% as.data.frame()
names(coordinates.full2) <- c("longitude", "latitude")


# i. mask hulls. takes ~1 minute (plotting time not included)
# in case sanity check in 2.iv is skipped, load relevant data
world.sf <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))

# get data for north american countries of interest
world.sp <- ne_countries(scale = "large", returnclass = "sf")
world.crop <- crop(as_Spatial(world.sp), extent(min(coordinates.full$longitude)-10, max(coordinates.full$longitude)+10,
                                    min(coordinates.full$latitude)-10, max(coordinates.full$latitude)+10))

r <- raster(ncol=2e3, nrow=2e3) # create raster object to fill into. bigger gives better resolution but takes longer

# mask out oceans
r.mask <- rasterize(world.crop, r, mask=F) 

# mask out lakes and buffer zone
lakes.sp <- as(lakes, "Spatial")
r.mask <- raster::mask(r.mask, lakes.sp, inverse=T)

buffer.sp <- as(gen.buffer, "Spatial")
r.mask <- raster::mask(r.mask, buffer.sp, inverse=F)

raster.pts <- rasterToPoints(r.mask)
r.sf <- st_as_sf(data.frame(raster.pts), coords = c("x", "y"), crs = st_crs(r.mask))



# get data for north american countries of interest
world.crop2 <- crop(as_Spatial(world.sp), extent(min(coordinates.full2$longitude)-10, max(coordinates.full2$longitude)+10,
                                                min(coordinates.full2$latitude)-10, max(coordinates.full2$latitude)+10))

r2 <- raster(ncol=2e3, nrow=2e3) # create raster object to fill into. bigger gives better resolution but takes longer

# mask out oceans
r.mask2 <- rasterize(world.crop2, r2, mask=F) 

# mask out lakes and buffer zone
lakes.sp <- as(lakes, "Spatial")
r.mask2 <- raster::mask(r.mask2, lakes.sp, inverse=T)

buffer.sp2 <- as(gen.buffer2, "Spatial")
r.mask2 <- raster::mask(r.mask2, buffer.sp2, inverse=F)

raster.pts2 <- rasterToPoints(r.mask2)
r.sf2 <- st_as_sf(data.frame(raster.pts2), coords = c("x", "y"), crs = st_crs(r.mask2))



expand <- read.csv("SDM/data/cam_pa2m30s_100km_expansion_1x_extended_fix1.csv")
app <- read.csv("SDM/data/cam_pa2m30s_100km_Appalachian_1x_extended_fix1.csv")

# save plot -- CHANGE RESOLUTION
jpeg(sprintf("SDM/plots/data_prep2m30s_100km_combined_extApp_fix1.jpeg", buffer.dist.N, region, pa.rep), res=5e2, units="in", height=5, width=6)
ggplot(data=world.sf)+
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  geom_sf(data=lakes, fill="white")+ #adds the great lakes
  geom_sf(data=r.sf, fill=NA, color="gray", alpha=0.1)+ # full mask
  geom_sf(data=gen.buffer, fill=NA, color="purple4", linewidth=0.75)+ # buffer hull (used in masking)
  geom_sf(data=r.sf2, fill=NA, color="gray", alpha=0.1)+ # full mask
  geom_sf(data=gen.buffer2, fill=NA, color="forestgreen", linewidth=0.75)+ # buffer hull (used in masking)
  # geom_point(data=pa.grid, aes(x=longitude,y=latitude), color="coral1", size=0.25)+
  ggnewscale::new_scale_color()+
  geom_point(data=expand, aes(x=longitude, y=latitude, color=as.factor(presence)), size=0.35)+
  scale_color_manual(values=c("red", "green4"))+
  geom_point(data=app, aes(x=longitude, y=latitude, color=as.factor(presence)), size=0.75, shape=17)+
  coord_sf(xlim=c(-72, -100), ylim=c(25, 47), expand=FALSE)+
  xlab("Latitude")+
  ylab("Longitude")+
  theme_bw()+
  theme(text = element_text(size = 12), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_rect(fill="white"), axis.line = element_line(colour = "black"))
dev.off()

jpeg(sprintf("SDM/plots/data_prep2m30s_100km_combined_APP_fix1.jpeg", buffer.dist.N, region, pa.rep), res=5e2, units="in", height=5, width=6)
ggplot(data=world.sf)+
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  geom_sf(data=lakes, fill="white")+ #adds the great lakes
  # geom_sf(data=r.sf, fill=NA, color="gray", alpha=0.1)+ # full mask
  # geom_sf(data=gen.buffer, fill=NA, color="purple4", linewidth=0.75)+ # buffer hull (used in masking)
  geom_sf(data=r.sf2, fill=NA, color="gray", alpha=0.1)+ # full mask
  geom_sf(data=gen.buffer2, fill=NA, color="forestgreen", linewidth=0.75)+ # buffer hull (used in masking)
  # geom_point(data=pa.grid, aes(x=longitude,y=latitude), color="coral1", size=0.25)+
  ggnewscale::new_scale_color()+
  # geom_point(data=expand, aes(x=longitude, y=latitude, color=as.factor(presence)), size=0.35)+
  scale_color_manual(values=c("red", "green4"))+
  geom_point(data=app, aes(x=longitude, y=latitude, color=as.factor(presence)), size=0.75, shape=17)+
  coord_sf(xlim=c(-72, -100), ylim=c(25, 47), expand=FALSE)+
  xlab("Latitude")+
  ylab("Longitude")+
  theme_bw()+
  theme(text = element_text(size = 12), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_rect(fill="white"), axis.line = element_line(colour = "black"))
dev.off()

jpeg(sprintf("SDM/plots/data_prep2m30s_100km_combined_EXP_fix1.jpeg", buffer.dist.N, region, pa.rep), res=5e2, units="in", height=5, width=6)
ggplot(data=world.sf)+
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  geom_sf(data=lakes, fill="white")+ #adds the great lakes
  geom_sf(data=r.sf, fill=NA, color="gray", alpha=0.1)+ # full mask
  geom_sf(data=gen.buffer, fill=NA, color="purple4", linewidth=0.75)+ # buffer hull (used in masking)
  # geom_sf(data=r.sf2, fill=NA, color="gray", alpha=0.1)+ # full mask
  # geom_sf(data=gen.buffer2, fill=NA, color="forestgreen", linewidth=0.75)+ # buffer hull (used in masking)
  # geom_point(data=pa.grid, aes(x=longitude,y=latitude), color="coral1", size=0.25)+
  ggnewscale::new_scale_color()+
  geom_point(data=expand, aes(x=longitude, y=latitude, color=as.factor(presence)), size=0.35)+
  scale_color_manual(values=c("red", "green4"))+
  # geom_point(data=app, aes(x=longitude, y=latitude, color=as.factor(presence)), size=0.75, shape=17)+
  coord_sf(xlim=c(-72, -100), ylim=c(25, 47), expand=FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  theme(text = element_text(size = 12), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_rect(fill="white"), axis.line = element_line(colour = "black"))
dev.off()

