
# libraries

# modeling
library(dismo) # max ent
library(rJava) # required by dismo
library(dplyr) # pipes
library(splitTools) # partitioning data
library(stringr)

# map projection
library(geodata) # worldclim spatrastering
library(sp) # cropping/clipping
library(sf) # cropping/clipping
library(raster) # extent clipping
library(terra) # dealing with SpatRasters
library(rnaturalearth) # setting up coordinate grid for prediction
library(maps) # for lakes and states

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/SDM/")

# define variables
set.seed(1066)
pca.perc = 0.95 # percent cut off for number of principal components to retain for modeling
region = "expansion" # expansion | Appalachian
buffer = 10 # lat/long buffer on world crop
resolution = 2.5 # worldclim data resolution. 2.5 recommended
bio.vers = "2.1"
buffer.dist = 100 # from dat prep

## Actions performed in this code:
# 1 - read in data and transform with PCA
# 2 - compute random forest model using contemporary climate data (1970-2000)
# 3 - predict N American presence/absence and save prediction raster

# getting full range data so PCA is calibrated to full species range
c.americana <- read.csv("~/Desktop/Documents/Research/Paper Code/adaptation_load/maps/data/iNaturalist_coords.csv")
bio_files <- list.files("~/Desktop/Documents/Research/Q3/ME_locus_detection/data/bioclimate/wc2.1_2.5m_bio/", pattern = "\\.tif$", full.names = TRUE) # List all the .tiff files in the directory
bioclim <- rast(bio_files) # Load all .tiff files as a SpatRaster object
c.am.bio <- terra::extract(bioclim, st_as_sf(c.americana, coords = c("longitude", "latitude"), crs = 4326))
c.am.bio.names <- names(c.am.bio)
names(c.am.bio) <- str_replace(c.am.bio.names, "wc2.1_2.5m_bio_", "bio")

c.am.full <- cbind(c.americana, c.am.bio[,paste0("bio", 1:19)])
cam.pca <- prcomp(c.am.full[,paste0("bio", 1:19)], center=T, scale.=T)

# extract the number of PC which explain at least N% of the data
cumsum <- summary(cam.pca)$importance["Cumulative Proportion",]
pc.num <- which(cumsum >= pca.perc)[1] %>% as.numeric()

## 1. read in data and transform with PCA (data generated in dat_prep.R)
if (region == "expansion") {
  c.americana <- read.csv(sprintf("data/cam_pa2m30s_%skm_expansion_1x_extended_fix1.csv", buffer.dist)) %>% na.omit() # 2 rows of NA's (coordinates are in ocean by NYC which escaped mask filter (resolution issues?)  
}

if (region == "Appalachian") {
  c.americana <- read.csv(sprintf("data/cam_pa2m30s_%skm_Appalachian_1x_extended_fix1.csv", buffer.dist)) %>% na.omit() # 2 rows of NA's (coordinates are in ocean by NYC which escaped mask filter (resolution issues?)  
}

cam.ex.pca <- predict(cam.pca, c.americana[,paste0("bio", 1:19)])


# extract and combine data
cam.pcN <- cam.ex.pca
c.americana <- cbind(c.americana[,c("presence", "latitude", "longitude")], cam.pcN)

## 2. compute random forest model using data obtained in dat_prep.R
# i - split data in testing and training to evaluate model accuracy/type I & II errors
# ii - run model on total aggregate data set (100 model)
# iii - set up worldclim data for broader prediction
# iv - predict full raster from each model and average across cells
# v - plot predictions and save contemporary prediction as raster

# i. split data into training (70%) and testing (30%) data sets
inds <- partition(c.americana$presence, p=c(train = 0.7, test = 0.3))
c.am.training <- c.americana[inds$train, ]
c.am.testing <- c.americana[inds$test, ]

# baseline model with 70% training set
set.seed(216) # set random set so model is replicable 

# run model on the number of PCA identified earlier as explaining 95% of climate variance (pc.num)
m_me <- maxent(c.am.training[,paste0("PC", 1:pc.num)], c.am.training[,"presence"] %>% as.factor())

# test error rates
prediction.table <- predict(m_me, c.am.testing[,paste0("PC", 1:pc.num)]) %>% round() #create a prediction table
con_matrix <- table(observed=c.am.testing[,"presence"],predicted=prediction.table) %>% as.data.frame()
false.neg <- con_matrix %>% filter(observed == 1 & predicted == 0) %>% dplyr::select(Freq) %>% as.numeric() / con_matrix %>% filter(predicted == 0) %>% dplyr::select(Freq) %>% sum()
false.pos <- con_matrix %>% filter(observed == 0 & predicted == 1) %>% dplyr::select(Freq) %>% as.numeric / con_matrix %>% filter(predicted == 1) %>% dplyr::select(Freq) %>% sum()

false.neg # extended_fix: expansion=0.113; Appalachian=0.207
false.pos # extended_fix: expansion=0.075; Appalachian=0.327

# ii. run full model 
m_me <- maxent(c.americana[,paste0("PC", 1:pc.num)], c.americana[,"presence"] %>% as.factor())

# geom_histogram()# iii. set up worldclim data for broader prediction

# download world data to crop extent of full bioclim data
world.sp <- ne_countries(scale = "large", returnclass = "sf")
world.crop <- crop(as_Spatial(world.sp), extent(min(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))-buffer-4, 
                                    max(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))+buffer+10, 
                                    min(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))-buffer, 
                                    max(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))+buffer+4))
crop.sf <- st_as_sf(world.crop)
coord.sf <- st_coordinates(crop.sf) %>% as.data.frame() %>% dplyr::select(X,Y)
names(coord.sf) <- c("lon", "lat")

if (bio.vers == "2.1") {
  bio.files <- list.files(path="~/Desktop/Documents/Research/Q3/SDM/wc2-5/", pattern = "^bio.*\\.bil$", full.names = TRUE)
  bio.world <- terra::rast(bio.files)
}

bio.crop <- terra::extract(bio.world, vect(world.crop), cells=T, xy=T)


bio.crop <- bio.crop[,-1] # gets rid of useless ID column
names(bio.crop) <- c(paste0("bio", 1:19), "cell", "longitude", "latitude") # make names match cam.pca

bio.bio <- bio.crop[,paste0("bio", 1:19)]
bio.xy <- bio.crop[,c("longitude", "latitude")]
bio.pca <- predict(cam.pca, bio.bio) # predict values of world data in existing PCA so as to not change PCA values. takes a while (~1-2 min.)
bio.pcN <- bio.pca[,1:pc.num]  # limit to same number of PC as before

# determine which cells have missing data and remove
bio.NAomit <- cbind(bio.xy, bio.pcN)
bio.NAomit <- bio.NAomit %>% na.omit()
coord.sf <- bio.NAomit %>% dplyr::select(longitude, latitude) %>% as.data.frame()
bio.pcN <- bio.NAomit %>% dplyr::select(paste0("PC", 1:pc.num)) %>% as.data.frame()

# iv. generate predictions from all models for range extent +/- 3 lat/long
pred.cell <- predict(m_me, bio.pcN) # predict(m_rf, bio.pcN, "prob")

# combine predictions with spatial data
pred.grid <- cbind(pred.cell, coord.sf)

# v. plot check
pred.grid <- pred.grid %>% mutate(bins = cut(pred.cell, breaks = c(-Inf,0.2,0.4,0.6,0.8,Inf))) %>% as.data.frame()
# pred.grid <- pred.grid %>% mutate(bins = cut(pred.cell, breaks = c(-Inf,seq(0,1,0.1),Inf))) %>% as.data.frame()

world.sf <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))

p <- ggplot()+
  geom_sf(data=world.sf, fill=NA)+
  geom_sf(data=states, fill=NA)+ #adds state boundaries
  geom_tile(data=pred.grid, aes(x=longitude, y=latitude, color=bins), size=1)+ #for crop
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  scale_color_manual(values = c("steelblue4",  "yellow3", "yellow", "chartreuse3", "green4"))+
  # scale_color_manual(values = c("dodgerblue4", "steelblue4", "steelblue", "gold4", "yellow3", "yellow",
  #                               "darkolivegreen2","darkolivegreen3", "chartreuse3", "green4", "darkgreen"))+
  geom_sf(data=lakes, fill="white")+ #adds the great lakes
  coord_sf(xlim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))-3, 
                  max(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))+3), 
           ylim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))-3, 
                  max(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))+3), expand=FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  # geom_point(data=subset(c.americana, presence==1), aes(x=longitude, y=latitude), colour="black", size=0.25)+ #adds lab collections
  theme_bw()+
  theme(text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_rect(fill="white"), axis.line = element_line(colour = "black"))

if (region == "expansion") {
  p <- p + coord_sf(xlim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))-3, 
                      max(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))+15), 
               ylim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))-3, 
                      max(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))+3), expand=FALSE)
}

jpeg(sprintf("plots/maxent/local/m_me_%s_%skm_extended_fix1.jpeg", region, buffer.dist), width=10, height=10, units="in", res=500)
print(p)
dev.off()

# save files for export
write.csv(pred.grid, sprintf("data/models/maxent/local/m_me_%skm_%s_extended_fix1.csv", buffer.dist, region), row.names = F)

# save raster
pred.raster <- rasterFromXYZ(pred.grid[,c("longitude", "latitude", "pred.cell")],
                             crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(pred.raster) # check that it worked
writeRaster(pred.raster, sprintf('data/models/maxent/local/m_me_%skm_%s_extended_fix1.tif', buffer.dist, region), overwrite=T)
