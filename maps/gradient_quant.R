library(raster) #for mapping and creating raster layers and bricks
# library(rgdal)
library(sp) #spatial libraries for mapping
library(spatialEco) #for pseudo-absence data generation
library(ggplot2) #graphing maps
library(randomForest) #RF model generation
library(RColorBrewer) #for transparent color palette generation
library(mecofun) #for predicting p/a data from RF model
# library(maptools)
library(rnaturalearth) #for downloading shapefiles
library(sf) #for processing shapefiles
library(dplyr) #to keep things tidy
library(ggspatial) #for scale bars and arrows
library(ggplot2) #for tidy plotting
library(SpatialEpi) #for converting latlong to km grid
library(devtools) #for image exporting 
library(FRK) #for converting dfs to spatial polygons 
library(maps) #for maps in graphs
library(plyr) #for round_any / km rounding from latlong coords
library(ggbiplot) #for PCA graphing
library(factoextra) #for PCA variance percentages
library(corrplot) #for correlation plot
library(data.table) #for duplicate removal
library(geosphere)
library(terra)
library(sf)
library(sp)
library(ggplot2)
library(geosphere)
library(reshape2)

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/maps/")

coordinates <- read.csv("data/iNaturalist_coords.csv", sep=',')
ID <- coordinates[,1]
sp <- SpatialPoints(coordinates[,c(3,2)])

#download world clim datasets
elev <- rast("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Documents/Research/Q3/ME_locus_detection/data/bioclimate/wc2.1_30s_elev.tif")
ex <- terra::extract(elev, st_as_sf(sp))
df_presence <- cbind.data.frame(ID,coordinates(sp),ex) # cbind ID, spatial points, and climate data together

# get environmental data associated with each point from an unbiased PCA
bio_files <- list.files("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Documents/Research/Q3/ME_locus_detection/data/bioclimate/wc2.1_2.5m_bio", pattern = "\\.tif$", full.names = TRUE) # List all the .tiff files in the directory
bio <- rast(bio_files) # load all .tiff files as a SpatRaster object
bio_dat <- terra::extract(bio, st_as_sf(sp)) %>% as.data.frame()

# global PCA 
c.am.pca <- prcomp(bio_dat, center=T, scale.=T)[["x"]]
summary(prcomp(bio_dat, center=T, scale.=T)) # take first 5 PCA 96.4%
all_dat <- cbind(coordinates, c.am.pca[,1:5])
names(all_dat)[1] <- "ID"

# # population data

# # for drift chapter populations
pops <- read.csv("data/pop_levels.csv", stringsAsFactors = T) 
pops <- pops %>% dplyr::filter(lineage != "E")

pops <- pops %>%
  dplyr::filter(AP == 0 & lineage != "E" & alt.pop != "RCA" & pop != "ALBG" & pop != "WV7") %>%
  dplyr::mutate(lat.2 = dplyr::case_when(
    pop == "VA5" ~ lat-0.25,
    pop != "VA5" ~ lat
  ))

pops$env.grp2 <- factor(pops$env.grp, levels=c("low", "mid", "high"))

# function for pairwies geographic distances
calculate_pairwise_distances <- function(lat_lon_data) {
  n <- nrow(lat_lon_data)
  distances <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      distances[i, j] <- distGeo(lat_lon_data[i,], lat_lon_data[j,])
      distances[j, i] <- distances[i, j] # distances matrix is symmetric
    }
  }
  rownames(distances) <- colnames(distances) <- 1:n
  return(distances)
}

###################################################################################################################
#
# ELEVATIONS ONLY
#
###################################################################################################################


# restricting to mountainous populations
elev_pops <- pops %>% dplyr::filter(elev_grp == 1) %>% dplyr::filter(pop != "PA104")
elev_pops <- rbind(elev_pops, c("VA3", NA, 0, 38.99, -77.99, 1, 0, "A", "low", "e4", NA, 38.99, "low")) %>% as.data.frame()
elev_pops <- elev_pops %>% 
  dplyr::mutate_at(dplyr::vars(lat, long, elev_m, lat.2), as.numeric) %>% 
  dplyr::mutate_at(dplyr::vars(elev_grp, lat_grp), as.factor)

# elevation PCA values from global
# sp <- SpatialPoints(elev_pops[,c("long", "lat")]) # for sequenced pops
sp <- SpatialPoints(elev_pops[,c("long", "lat")])
ex <- terra::extract(bio, st_as_sf(sp)) %>% as.data.frame()

elev_pca <- predict(prcomp(bio_dat, center=T, scale.=T), ex) %>% as.data.frame()
elev_pops <- cbind(elev_pops, elev_pca[,1:5])


# find distances in env/PCA space and geographic distance and reduce to keep comparisons within 50km of each other
# geographic distance
# tt <- calculate_pairwise_distances(elev_pops[,c("long", "lat")]) # for seq pops
tt <- calculate_pairwise_distances(elev_pops[,c("long", "lat")]) # for drift pops
rownames(tt) <- elev_pops$pop # pop for seq
colnames(tt) <- elev_pops$pop # pop for seq
tt <- tt/1000 # convert to km from m
elev_gdist <- melt(tt)
names(elev_gdist) <- c("pop.1", "pop.2", "dist_km")

# environmental distance 
env_distances <- dist(elev_pops[,c("PC1", "PC2")]) %>% as.matrix()
rownames(env_distances) <- elev_pops$pop # pop for seq
colnames(env_distances) <- elev_pops$pop # pop for seq
elev_edist <- melt(env_distances)
names(elev_edist) <- c("pop.1", "pop.2", "env_dist_pca")
elev_dists <- merge(elev_gdist, elev_edist, by=c("pop.1", "pop.2"))

# re-add metadata
elev_pops$grad.env <- paste0(elev_pops$lineage, "_", elev_pops$elev_grp, "_", elev_pops$env.grp)
elev_pops$pop.1 <- elev_pops$pop
elev_pops$pop.2 <- elev_pops$pop

test2 <- merge(elev_dists, elev_pops[,c("lineage", "elev_grp", "grad.env", "env.grp", "transect", "pop.1")], by="pop.1")
test2 <- merge(test2, elev_pops[,c("lineage", "elev_grp", "grad.env", "env.grp", "transect", "pop.2")], by="pop.2")
head(test2)

# filter to only pops of interest along transects
# for drift
test2$position.comp <- paste0(test2$env.grp.x, "_", test2$env.grp.y)

elev_dists2 <- test2 %>% 
  dplyr::filter(pop.1 != pop.2) %>% 
  dplyr::filter(elev_grp.x == elev_grp.y & lineage.x == lineage.y & transect.x == transect.y) %>% 
  dplyr::filter(elev_grp.x == 1) %>% 
  dplyr::filter(position.comp == "high_low")

# median gradient
e_grad <- elev_dists2 %>% 
  # select(pop.1, pop.2, dist_km, env_dist_pca) %>%
  dplyr::filter(transect.x == transect.y) %>%
  dplyr::mutate(gradient = env_dist_pca / dist_km)




###################################################################################################################
#
# LATITUDES ONLY
#
###################################################################################################################

# restricting to mountainous populations
lat_pops <- pops %>% dplyr::filter(elev_grp == 0)

# map it out
world <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))

ggplot(data=world)+
  geom_sf()+
  geom_sf(data=states, fill=NA)+ #adds state boundaries
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  geom_sf(data=lakes, fill="gray40")+ #adds the great lakes
  coord_sf(xlim=c(-98,-83), ylim=c(36,47), expand=FALSE)+ #sets boundaries to extent of model
  xlab("Latitude")+
  ylab("Longitude")+
  # geom_point(data=lat_pops, aes(x=long, y=lat), colour="black", size=1.5)+ # sequenced pops
  geom_point(data=lat_pops, aes(x=long, y=lat, group=transect), size=4)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


# latation PCA values from global
# sp <- SpatialPoints(lat_pops[,c("long", "lat")]) # for sequenced pops
sp <- SpatialPoints(lat_pops[,c("long", "lat")])
ex <- terra::extract(bio, st_as_sf(sp)) %>% as.data.frame()

lat_pca <- predict(prcomp(bio_dat, center=T, scale.=T), ex) %>% as.data.frame()
lat_pops <- cbind(lat_pops, lat_pca[,1:5])


# find distances in env/PCA space and geographic distance and reduce to keep comparisons within 50km of each other
# geographic distance
# tt <- calculate_pairwise_distances(lat_pops[,c("long", "lat")]) # for seq pops
tt <- calculate_pairwise_distances(lat_pops[,c("long", "lat")]) # for drift pops
rownames(tt) <- lat_pops$pop # pop for seq
colnames(tt) <- lat_pops$pop # pop for seq
tt <- tt/1000 # convert to km from m
lat_gdist <- melt(tt)
names(lat_gdist) <- c("pop.1", "pop.2", "dist_km")

# environmental distance 
env_distances <- dist(lat_pops[,c("PC1", "PC2")]) %>% as.matrix()
rownames(env_distances) <- lat_pops$pop # pop for seq
colnames(env_distances) <- lat_pops$pop # pop for seq
lat_edist <- melt(env_distances)
names(lat_edist) <- c("pop.1", "pop.2", "env_dist_pca")
lat_dists <- merge(lat_gdist, lat_edist, by=c("pop.1", "pop.2"))

# re-add metadata
lat_pops$grad.env <- paste0(lat_pops$lineage, "_", lat_pops$elev_grp, "_", lat_pops$env.grp)
lat_pops$pop.1 <- lat_pops$pop
lat_pops$pop.2 <- lat_pops$pop

test <- merge(lat_dists, lat_pops[,c("lineage", "elev_grp", "grad.env", "env.grp", "transect", "pop.1")], by="pop.1")
test <- merge(test, lat_pops[,c("lineage", "elev_grp", "grad.env", "env.grp", "transect", "pop.2")], by="pop.2")
head(test)

# filter to only pops of interest along transects
# for drift
test$position.comp <- paste0(test$env.grp.x, "_", test$env.grp.y)
head(test)

lat_dists2 <- test %>% 
  dplyr::filter(pop.1 != pop.2) %>% 
  dplyr::filter(elev_grp.x == elev_grp.y & lineage.x == lineage.y & transect.x == transect.y) %>% 
  dplyr::filter(elev_grp.x == 0) %>% 
  dplyr::filter(position.comp == "high_mid" | position.comp == "mid_low")

# median gradient
l_grad <- lat_dists2 %>% 
  # select(pop.1, pop.2, dist_km, env_dist_pca) %>% 
  dplyr::mutate(gradient = env_dist_pca / dist_km)

median(l_grad$gradient) # mean=0.0083, median=0.0083

l_grad <- l_grad %>% 
  dplyr::mutate_at(vars(elev_grp.x, elev_grp.y), as.factor)
grad <- rbind(e_grad, l_grad)

grad <- grad %>% 
  # dplyr::group_by(elev_grp.x, lineage.x) %>% 
  dplyr::group_by(elev_grp.x) %>% 
  dplyr::reframe(mean_grad = mean(gradient),
                 sd_grad = sd(gradient),
                 se_grad = sd(gradient)/dplyr::n())

grad$se_low <- grad$mean_grad-grad$se_grad
grad$se_high <- grad$mean_grad+grad$se_grad

grad <- grad %>% 
  dplyr::mutate(gradient = dplyr::case_when(
    elev_grp.x == 1 ~ "Elevation",
    elev_grp.x == 0 ~ "Latitude"
  )) %>% 
  dplyr::mutate(gradient.n = dplyr::case_when(
    elev_grp.x == 0 ~ "Shallow",
    elev_grp.x == 1 ~ "Steep"
  ))

jpeg("plots/GRADIENT_combined.jpeg", res=500, units="in", height=3.9, width=5)
ggplot()+
  geom_errorbar(data=grad, aes(x=gradient.n, ymin=se_low, ymax=se_high), width = 0.2)+
  geom_point(data=grad, aes(x=gradient.n, y=mean_grad, shape=gradient.n), color="black", size=4)+
  geom_point(data=grad, aes(x=gradient.n, y=mean_grad, shape=gradient.n), size=4)+
  scale_shape_manual(values=c(21, 24))+
  # scale_color_manual(values=c("green4", "purple2"))+
  ylim(c(0,NA))+
  theme_bw()+
  ylab("β Environment")+
  # ggtitle("Median steepness of environmental gradient by lineage and region")+
  xlab("Gradient")+
  theme(text = element_text(size = 11), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio=3)
dev.off()
