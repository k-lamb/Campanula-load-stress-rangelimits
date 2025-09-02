
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
set.seed(42) 

# Declare variables
grid_size = 25 # km resolution for grid cells
buffer.dist = 200 # buffer distance for restricting pseudo-absences
buffer.dist.N = 100 # buffer distance for restricting presence points by phylogenetic lineage
resolution = 2.5 # worldclim resolution. finer resolution takes longer (0.5, 2.5, etc.)
pa.rep = 1 # multiplier for creating pseudo-absence data
sanity.plot = T # whether to plot sanity checks
neighbor.dist = 200 # distance from nearest neighbor in km
size.t = 16 # size of text in plots
lat.cutoff = 33.5 # latitude cutoff for determining when buffer is used and when latitude is used
lat.cutoff.TF = FALSE # if you want to use the latitude cut off above. my recommendation is to set to F, as you only gain 2 remote points
region = "expansion" # "rear", "mid.lat", "leading", "expansion" (mid.lat and leading), "Appalachian"
bio.vers = "2.1" # which version of bioclim to use. 2.1 is higher resolution but does not work for LGM projection
method.pa.prune = "cut.out.inside" # "cut.around.presence" | "cut.out.inside" <-- whether to allow pseudo-absences to generate within the range or whether to use a hull to exclude them

### DATA PREPARATION

## Actions performed in this code:
# 1 - Import iNaturalist data
# 2 - Spatial down-sampling of presence data
# 3 - Generate pseudo-absence data bound by some reasonable proximity to known points
# 4 - Download/import WorldClim 2.1 data
# 5 - Save data


## 1. Import set of coordinates, with columns "id", "latitude" and "longitude"
# i - import data
# ii - limit presence data based on the region of interest

# i. import data
# iNaturalist data
coordinates.pts <- read.csv("maps/data/iNaturalist_coords.csv", sep=',') # data from iNaturalist to May 2022

# expansion region (mid.lat and leading)
# mid-latitude (KY refugium)
if (region == "expansion"){
  
  coordinates.full <- coordinates.pts
  pa.grid_size = 50 # km resolution for grid cells for p-pa distance
  
  # remove points from Appalachian and Smokies lineages by elevational cut-off
  elev <- elevation_30s("USA", path="~/Desktop/Documents/Research/Q3/SDM/geodata/")
  coordinates.full$elevation <- extract(elev, coordinates.full[,c("longitude", "latitude")])$USA_elv_msk
  coordinates.full <- coordinates.full %>% filter(elevation <= quantile(coordinates.full$elevation, 0.85, na.rm=T))
  
  # get pops close to sequenced populations of known structure
  pop.lim <- read_excel("./SDM/data//Samples detail.xlsx") %>% 
    as.data.frame() %>% 
    mutate_at(vars(POP...1, Clade), as.factor) %>% 
    mutate_at(vars(long, lat), as.numeric) %>% 
    dplyr::select(POP...1, Clade, long, lat)
  names(pop.lim) <- c("pop", "clade", "longitude", "latitude")
  
  # pop.lim <- pop.lim %>% filter(pop == "KY5" | pop == "KY1" | pop == "TN34" | pop == "OH119" | pop == "IN5" | pop == "IN7" | 
  #                                 pop == "MI2" | pop == "IN2" | pop == "TN3" | pop == "IA17" | pop == "IL10" | pop == "WI4" | 
  #                                 pop == "MO2" | pop == "IA12" | pop == "MN8" | pop == "MN117" | pop == "KS60" | pop == "OK1" | 
  #                                 pop == "OK61" | pop == "AR2")
  
  pop.lim <- pop.lim %>% filter(latitude > 35 & longitude < -83 & clade == "Western")
  
  # generate a buffered concave polygon around the desired samples
  gen.sf <- st_as_sf(pop.lim, coords=c("longitude", "latitude"), crs=4326)
  gen.hull <- concaveman(gen.sf)
  gen.buffer <- st_buffer(gen.hull, dist = buffer.dist.N*1e3)
  
  # apply buffer around sampling to grab points of interest
  coordinates.sf <- st_as_sf(coordinates.full[,c("longitude", "latitude")], coords=c("longitude", "latitude"), crs=4326)
  coordinates.buffer <- st_within(coordinates.sf, gen.buffer, sparse=F)
  coordinates.full <- coordinates.sf[coordinates.buffer, ]
  coordinates.full <- st_coordinates(coordinates.full) %>% as.data.frame()
  names(coordinates.full) <- c("longitude", "latitude")
  crs <- 32715
  
}

# expansion region (mid.lat and leading)
# mid-latitude (KY refugium)
if (region == "Appalachian"){
  
  coordinates.full <- coordinates.pts
  pa.grid_size = 50 # km resolution for grid cells for p-pa distance
  
  pop.lim <- read_excel("./SDM/data/C americana populations.xlsx") %>% as.data.frame() 
  pop.lim <- pop.lim[,c(5,9,7,6,10)] # pop, lineage, longitude, latitude
  pop.lim <- pop.lim[-1,]
  names(pop.lim) <- c("pop", "clade", "longitude", "latitude", "state")
  pop.lim <- pop.lim %>% filter(clade == "Appalachian" | state == "North Carolina" | (state == "Georgia" & latitude > 34))
  
  alt <- terra::rast("~/Downloads/wc2.1_30s_elev.tif")
  alt.coord <- terra::extract(alt, st_as_sf(pop.lim %>% na.omit(), coords=c("longitude", "latitude"), crs=4326))
  
  # generate a buffered concave polygon around the desired samples
  gen.sf <- st_as_sf(pop.lim, coords=c("longitude", "latitude"), crs=4326)
  gen.hull <- concaveman(gen.sf)
  gen.buffer <- st_buffer(gen.hull, dist = buffer.dist.N*1e3)
  
  # apply buffer around sampling to grab points of interest
  coordinates.sf <- st_as_sf(coordinates.full[,c("longitude", "latitude")], coords=c("longitude", "latitude"), crs=4326)
  coordinates.buffer <- st_within(coordinates.sf, gen.buffer, sparse=F)
  coordinates.full <- coordinates.sf[coordinates.buffer, ]
  coordinates.full <- st_coordinates(coordinates.full) %>% as.data.frame()
  names(coordinates.full) <- c("longitude", "latitude")
  crs <- 32717
  
}


## 2. Spatial down-sample of presence data
# i - create spatial grids
# ii - remove outlier points (bad observations?)
# iii - remove duplicate points (down-sample to 1 point per grid cell)
# iv - sanity check plot

# i. create spatial grids (10km*10km cell size)
grid <- round(latlong2grid(coordinates.full[,c("longitude","latitude")])) # gets km grid and rounds to nearest whole km
grid$x <- round_any(grid$x, grid_size)
grid$y <- round_any(grid$y, grid_size)
coordinates.full <- cbind(coordinates.full, grid)

# ii. remove outlier points ≥200km away from the nearest neighboring point (<1min to run)
dist.mat <- distm(coordinates.full %>% dplyr::select("longitude", "latitude"), fun=distHaversine)
diag(dist.mat) <- NA # remove distances to self (dist=0)

# coordinates <- coordinates.full
grid_m <- grid_size * 1000
coordinates <- data.frame(longitude=-50, latitude=0, x=50, y=0) # fake coordinate so there is at least one FALSE
coordinates <- st_as_sf(coordinates, coords = c("longitude", "latitude"), crs = 4326)
coordinates <- st_transform(coordinates, crs) # crs for illinois
min.distances <- c()

save <- coordinates.full
temp <- st_as_sf(save, coords = c("longitude", "latitude"), crs = 4326)
coordinates.full <- st_transform(temp, crs)

for (i in 2:nrow(coordinates.full)) {
  temp <- rbind(coordinates, coordinates.full[i,])
  distances <- st_distance(temp[-nrow(temp),], temp[nrow(temp),]) %>% as.numeric()
  
  if(all(distances >= grid_m)) {
    coordinates <- rbind(coordinates, coordinates.full[i,])
    min.distances <- c(min.distances, min(distances))
  }
  
  if (i %% 500 == 0) {
    print(paste("working on point:", i))
  }
}

coordinates <- coordinates[-1,] # removes fake coordinate point
coordinates <- coordinates %>% st_transform(., 4326) %>% st_coordinates() %>% as.data.frame()
names(coordinates) <- c("longitude", "latitude")
dim(coordinates)

coordinates.full <- coordinates.full %>% st_transform(., 4326) %>% st_coordinates() %>% as.data.frame()
names(coordinates.full) <- c("longitude", "latitude")

# iv. sanity check plotting
# set up map attributes
world.sf <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(maps::map("lakes", plot = FALSE, fill = TRUE))

if (sanity.plot == T) {
  ggplot(data=world.sf)+
    geom_sf()+
    geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
    geom_sf(data=lakes, fill="white")+ #adds the great lakes
    coord_sf(xlim=c(min(coordinates$longitude)-3, max(coordinates$longitude)+3), 
             ylim=c(min(coordinates$latitude)-3, max(coordinates$latitude)+3), expand=FALSE)+
    geom_point(data=coordinates, aes(x=longitude, y=latitude), size=0.5)+
    # geom_text(data=coordinates, aes(x=longitude, y=latitude, label=id))+
    xlab("Latitude")+
    ylab("Longitude")+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}


## 3. Generate pseudo-absence data
# i - mask the species range to avoid generating points outside buffer, in ocean, or in lakes
# ii - generate pseudo-absence data within mask and away from ANY known presence observation
# iii - merge pseudo-absence and presence data

# i. mask hulls. takes ~1 minute (plotting time not included)

# get data for north american countries of interest
world.sp <- ne_countries(scale = "large", returnclass = "sf")
world.crop <- crop(as_Spatial(world.sp), extent(min(coordinates.full$longitude)-10, max(coordinates.full$longitude)+10,
                                                min(coordinates.full$latitude)-10, max(coordinates.full$latitude)+10))

r <- raster(ncol=2e3, nrow=2e3) # create raster object to fill into. bigger gives better resolution but takes longer
r.mask <- rasterize(world.crop, r, mask=F) # mask out oceans
lakes.sp <- as(lakes, "Spatial")
r.mask <- raster::mask(r.mask, lakes.sp, inverse=T) # mask out lakes and buffer zone



## ii. generate pseudo-absence data

coord.sp <- SpatialPoints(coordinates.full[,c("longitude", "latitude")])
gen.sf <- st_as_sf(coordinates, coords=c("longitude", "latitude"), crs=4326)
gen.hull <- concaveman(gen.sf)
gen.buffer <- st_buffer(gen.hull, dist = buffer.dist*1e3)

buffer.sp <- as(gen.buffer, "Spatial")
r.mask <- raster::mask(r.mask, buffer.sp, inverse=F)



pa.grid <- expand.grid(long=seq(-100,-70, by=0.2), lat=seq(30,50, by=0.2))
pa.grid$long <- pa.grid$long + runif(nrow(pa.grid), -0.1, 0.1)
pa.grid$lat  <- pa.grid$lat + runif(nrow(pa.grid), -0.1, 0.1)

# remove pseudo-absences that are too far away and too close
coord.sf <- st_as_sf(coordinates, coords=c("longitude", "latitude"), crs=4326)
coord.hull <- concaveman(coord.sf)
coord.buffer <- st_buffer(coord.hull, dist = pa.grid_size*1e3)
coord.buffer.2 <- st_buffer(coord.hull, dist = buffer.dist*1e3)
coord.buffer.sp <- as(coord.buffer, "Spatial")
coord.buffer.2.sp <- as(coord.buffer.2, "Spatial")

pa.sf <- st_as_sf(pa.grid, coords=c("long", "lat"), crs=4326)


if (method.pa.prune == "cut.out.inside") {
  r.mask <- raster::mask(r.mask, coord.buffer.sp, inverse=T) # cuts out inside
}

r.mask <- raster::mask(r.mask, coord.buffer.2.sp, inverse=F)
values <- raster::extract(r.mask, st_coordinates(pa.sf))
pa.grid <- pa.grid[!is.na(values), ] %>% as.data.frame()
names(pa.grid) <- c("longitude", "latitude")

ggplot(data=pa.grid, aes(x=longitude, y=latitude))+
  geom_point()

# remove close points in absences
grid_m <- grid_size * 1e3
# if (region == "expansion") { 
#   grid_m <- grid_m * 2 # multiplying for expansion only since habitat is more uniform along latitude than in the mountains
# }

pa.start <- data.frame(longitude=-50, latitude=0) # fake coordinate so there is at least one FALSE
pa.start <- st_as_sf(pa.start, coords = c("longitude", "latitude"), crs = 4326)
pa.start <- st_transform(pa.start, crs) # crs for illinois
min.distances <- c()

temp <- st_as_sf(pa.grid, coords = c("longitude", "latitude"), crs = 4326)
pa.grid <- st_transform(temp, crs)

for (i in 2:nrow(pa.grid)) {
  temp <- rbind(pa.start, pa.grid[i,])
  # distances <- distm(temp[-nrow(temp),c("longitude", "latitude")] %>% as.matrix(), temp[nrow(temp),c("longitude", "latitude")], fun=distVincentyEllipsoid)
  distances <- st_distance(temp[-nrow(temp),], temp[nrow(temp),]) %>% as.numeric()
  
  if(all(distances >= grid_m)) {
    pa.start <- rbind(pa.start, pa.grid[i,])
    min.distances <- c(min.distances, min(distances))
  }
  
  if (i %% 1e3 == 0) {
    print(paste("working on point:", i))
  }
}

pa.start <- pa.start[-1,] # removes fake coordinate point
pa.start <- pa.start %>% st_transform(., 4326) %>% st_coordinates() %>% as.data.frame()
names(pa.start) <- c("longitude", "latitude")
dim(pa.start)

# double checking it worked
ggplot(data=pa.start, aes(x=longitude, y=latitude))+
  geom_point()

pa.grid <- pa.start
coord.pa <- coordinates

if(method.pa.prune == "cut.around.presence") {
  coord.pa <- st_as_sf(coord.pa, coords = c("longitude", "latitude"), crs = 4326)
  pa.grid <- st_as_sf(pa.grid, coords = c("longitude", "latitude"), crs = 4326)
  for (i in 1:nrow(pa.grid)) {
    temp <- rbind(coord.pa, pa.grid[i,])
    # distances <- distm(temp[-nrow(temp),c("longitude", "latitude")] %>% as.matrix(), temp[nrow(temp),c("longitude", "latitude")], fun=distVincentyEllipsoid)
    distances <- st_distance(temp[-c((nrow(coordinates)+1):nrow(temp)),], temp[nrow(temp),]) %>% as.numeric()
    
    if(all(distances >= pa.grid_size*1e3)) {
      coord.pa <- rbind(coord.pa, pa.grid[i,])
      min.distances <- c(min.distances, min(distances))
    }
    
    if (i %% 500 == 0) {
      print(paste("working on point:", i))
    }
  }
  
  coord.pa <- coord.pa[-1,] # removes fake coordinate point
  coord.pa <- coord.pa %>% st_coordinates() %>% as.data.frame()
  names(coord.pa) <- c("longitude", "latitude")
  pa.grid <- coord.pa
}


coordinates$presence <- 1
pa.grid$presence <- 0
c.americana <- rbind(coordinates %>% dplyr::select(longitude, latitude, presence), pa.grid)
c.americana <- c.americana %>% group_by(presence) %>% slice_sample(n=round((nrow(coordinates)*1))) %>% ungroup()
names(c.americana)[1:2] <- c("longitude", "latitude")

# sanity check plot
raster.pts <- rasterToPoints(r.mask)
r.sf <- st_as_sf(data.frame(raster.pts), coords = c("x", "y"), crs = st_crs(r.mask))

if (sanity.plot == T) {
  ggplot(data=world.sf)+
    geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
    geom_sf(data=lakes, fill="white")+ #adds the great lakes
    geom_sf(data=r.sf, fill=NA, color="gray", alpha=0.3)+ # full mask
    geom_sf(data=gen.buffer, fill=NA, color="blue3")+ # buffer hull (used in masking)
    geom_point(data=c.americana, aes(x=longitude, y=latitude, color=as.factor(presence)), size=0.5)+
    scale_color_manual(values=c("coral3", "green4"))+
    coord_sf(xlim=c(-72, -100), ylim=c(25, 47), expand=FALSE)+
    xlab("Latitude")+
    ylab("Longitude")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
}

## 4. download and prep bioclim data from WorldClim 2.1
# i - download and extract data (contemporary and LGM)
# ii - combine data 

### altitude for appalachian only
if (region == "Appalachian") {
  alt <- terra::rast("~/Downloads/wc2.1_30s_elev.tif")
  alt.coord <- read.csv("SDM/data/cam_pa2m30s_100km_Appalachian_1x_extended_fix1.csv")
  alt.coord <- terra::extract(alt, st_as_sf(alt.coord %>% filter(presence==1), coords=c("longitude", "latitude"), crs=4326))
  range(alt.coord$wc2.1_30s_elev)
  
}

# CONTEMPORARY data
# convert points for extracting bioclim data
cam.sf <- st_as_sf(c.americana, coords=c("longitude", "latitude"), crs=4326)

# download spatial data (1km^2 resolution)
bio.files <- list.files(path="~/Desktop/Documents/Research/Q3/SDM/wc2-5/", pattern = "^bio.*\\.bil$", full.names = TRUE)
bio.world <- terra::rast(bio.files)
cam.bio <- terra::extract(bio.world, cam.sf) 

cam.bio <- cam.bio[,-1]
names(cam.bio) <- paste0("bio", 1:19) # rename variable
c.americana <- cbind(c.americana, cam.bio) # combine data


# leading edge
if (resolution == 0.5) {
  write.csv(c.americana, sprintf("SDM/data/cam_pa30s_%skm_%s_%sx.csv", buffer.dist.N, region, pa.rep), row.names = F)
}

if (resolution == 2.5) {
  write.csv(c.americana, sprintf("SDM/data/cam_pa2m30s_%skm_%s_%sx_extended_fix1.csv", buffer.dist.N, region, pa.rep), row.names = F)
}

# save plot -- CHANGE RESOLUTION
jpeg(sprintf("SDM/plots/data_prep2m30s_%skm_%s_%sx_extended_fix1.jpeg", buffer.dist.N, region, pa.rep), res=5e2, units="in", height=5, width=6)
ggplot(data=world.sf)+
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  geom_sf(data=lakes, fill="white")+ #adds the great lakes
  geom_sf(data=r.sf, fill=NA, color="gray", alpha=0.3)+ # full mask
  geom_sf(data=gen.buffer, fill=NA, color="blue3")+ # buffer hull (used in masking)
  geom_point(data=c.americana, aes(x=longitude, y=latitude, color=as.factor(presence)), size=0.5)+
  scale_color_manual(values=c("coral3", "green4"))+
  coord_sf(xlim=c(-72, -100), ylim=c(25, 47), expand=FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = "none")
dev.off()

