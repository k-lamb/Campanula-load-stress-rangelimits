
# data prep
library(dplyr)
# library(tidyverse)

# geospatial
library(geodata) # Worldclim 2.1 bioclim download
library(terra)
library(raster)
library(sp)
library(sf)

# plotting
library(ggplot2)
library(ggfortify)

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/")

# data import
coord.site <- read.csv("local_adaptation/data/CG_locations.csv")
coord.pop <- read.csv("local_adaptation/data/pops_mapping.csv")
coord.site <- subset(coord.site, range=="A")
coord.pop <- subset(coord.pop,  range ==  "Appalachia")

cg.no_herb <- read.csv("local_adaptation/data/elevation/elev_rel_abs_NO_HERB.csv", stringsAsFactors = T)
het <- read.csv("load/data/heterosis/heterosis_all.csv") %>% dplyr::select(-c("lat", "long", "transect", "env.var", "env.group"))

sdm.pred <- raster("SDM/data/models/rf/rf_100avg.tif") # deprecated -- preserving so column numbers don't change
me.sdm.pred <- raster("SDM/data/models/maxent/m_me_2m30s.tif") # deprecated -- preserving so column numbers don't change

rf.sdm.app.pred <- raster("SDM/data/models/rf/local/m_rf_100km_Appalachian_extended_fix1.tif")
me.sdm.app.pred <- raster("SDM/data/models/maxent/local/m_me_100km_Appalachian_extended_fix1.tif")


pc.dat <- read.csv("maps/data/iNaturalist_coords.csv")
# bio.world <- worldclim_global(var="bio", res=2.5, path="~/Desktop/Documents/Research/Q3/SDM/geodata/", version="2.1")
bio.files <- list.files(path="~/Desktop/Documents/Research/Q3/SDM/wc2-5/", pattern = "^bio.*\\.bil$", full.names = TRUE)
bio.world <- terra::rast(bio.files)

### DATA PREP FOR CLIMATE PCA
ID.site <- coord.site[,"common.gardens"]
ID.pop <- coord.pop[,"pop.ID"]

sp.site <- SpatialPoints(coord.site[,c(5,4)])
sp.pop <- SpatialPoints(coord.pop[,c(13,12)])
sp <- SpatialPoints(pc.dat[,c(3,2)])

# extract climate data
ex <- terra::extract(bio.world, st_as_sf(sp)) 
ex.site <- terra::extract(bio.world, st_as_sf(sp.site))
ex.pop <- terra::extract(bio.world, st_as_sf(sp.pop))

df <- cbind.data.frame(coordinates(sp), ex) # cbind ID, spatial points, and climate data together
df.site <- cbind.data.frame(ID.site, coordinates(sp.site), ex.site)
df.pop <- cbind.data.frame(ID.pop,coordinates(sp.pop), ex.pop)

# PCA from all iNaturalist first
pca <- prcomp(df[,paste0("bio", 1:19)],scale=T, center=T)
summary(pca)


# project population and site data to PCA
pr.pop <- predict(pca, df.pop[,paste0("bio", 1:19)])
pc1.pop <- as.data.frame(pr.pop)$PC1
pc2.pop <- as.data.frame(pr.pop)$PC2

pr.site <- predict(pca, df.site[,paste0("bio", 1:19)])
pc1.site <- as.data.frame(pr.site)$PC1
pc2.site <- as.data.frame(pr.site)$PC2

# add data by groups
cg.no_herb <- cg.no_herb %>% mutate(pc1.sites = case_when(
  full.garden == "MLBS" ~ pc1.site[2],
  full.garden == "kentland_farm" ~ pc1.site[3],
  full.garden == "highlands" ~ pc1.site[4],
  full.garden == "clemson" ~ pc1.site[5]
))

cg.no_herb <- cg.no_herb %>% mutate(pc2.sites = case_when(
  full.garden == "MLBS" ~ pc2.site[2],
  full.garden == "kentland_farm" ~ pc2.site[3],
  full.garden == "highlands" ~ pc2.site[4],
  full.garden == "clemson" ~ pc2.site[5]
))

cg.no_herb <- cg.no_herb %>% mutate(pc1.pop = case_when(
  pop == "SK" ~ pc1.pop[1],
  pop == "SM" ~ pc1.pop[2],
  pop == "VA73" ~ pc1.pop[3],
  pop == "Egg" ~ pc1.pop[4],
  pop == "GA2" ~ pc1.pop[5],
  pop == "GA1" ~ pc1.pop[6]
))

cg.no_herb <- cg.no_herb %>% mutate(pc2.pop = case_when(
  pop == "SK" ~ pc2.pop[1],
  pop == "SM" ~ pc2.pop[2],
  pop == "VA73" ~ pc2.pop[3],
  pop == "Egg" ~ pc2.pop[4],
  pop == "GA2" ~ pc2.pop[5],
  pop == "GA1" ~ pc2.pop[6]
))

# calculate difference between garden and origin
cg.no_herb$delta_pc1 <- cg.no_herb$pc1.sites - cg.no_herb$pc1.pop
cg.no_herb$delta_pc2 <- cg.no_herb$pc2.sites - cg.no_herb$pc2.pop


temp <- cbind(pc.dat[,c(3,2)], pca[["x"]])

# autoplot(pca, loadings=T, alpha=0, loadings.label=F)+
ggplot(temp, aes(x=PC1, y=PC2))+
  theme_bw()+
  geom_point(size=0.1, color="gray70")+
  # geom_label(data=cg.no_herb, aes(x=pc1.sites, y=pc2.sites, label=full.garden), color="darkorchid", size=4, alpha=1)+
  # geom_point(data=cg.no_herb, aes(x=pc1.pop, y=pc2.pop), color="black", shape=17, size=3.5, alpha=1)+
  geom_point(data=cg.no_herb, aes(x=pc1.pop, y=pc2.pop, fill=pop.group), shape=24, size=5, alpha=1)+
  # geom_point(data=cg.no_herb, aes(x=pc1.sites, y=pc2.sites), color="black", shape=15, size=4.5, alpha=1)+
  geom_point(data=cg.no_herb, aes(x=pc1.sites, y=pc2.sites, fill=garden.group), shape=22, size=5, alpha=1)+
  # geom_label(data=cg.no_herb, aes(x=pc1.pop+1.25, y=pc2.pop, label=pop, color=pop.group), size=4, alpha=1)+
  scale_fill_manual(values=c("low" = "red", "mid"="purple", "high" = "blue"))+
  # stat_bin_hex(bins=20)+
  theme(text = element_text(size=20), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio=1)
# ggsave("LA_HET/plots/ELEV_PC.jpeg")


### DATA PREP FOR SDM 

# populations
coord.pop <- coord.pop[,c("pop.ID", "pop_lat", "pop_long")]
spp <- SpatialPoints(coord.pop[,c(3,2)])

sdm.pop <- raster::extract(sdm.pred, spp)
me.sdm.pop <- raster::extract(me.sdm.pred, spp)
rf.sdm.app.pop <- raster::extract(rf.sdm.app.pred, spp)
me.sdm.app.pop <- raster::extract(me.sdm.app.pred, spp)

coord.pop$pop.ll <- sdm.pop
coord.pop$me.pop.ll <- me.sdm.pop
coord.pop$pop.ll.app <- rf.sdm.app.pop
coord.pop$me.pop.ll.app <- me.sdm.app.pop
names(coord.pop) <- c("pop", "pop.lat", "pop.long", "rf.pop.ll", "me.pop.ll", "rf.pop.app.ll", "me.pop.app.ll")

# garden sites
coord.sites <- coord.site[,5:4]
sps <- SpatialPoints(coord.sites)

sdm.site <- raster::extract(sdm.pred, sps)
me.sdm.site <- raster::extract(me.sdm.pred, sps)
rf.sdm.app.site <- raster::extract(rf.sdm.app.pred, sps)
me.sdm.app.site <- raster::extract(me.sdm.app.pred, sps)

coord.site$rf.site.ll <- sdm.site
coord.site$me.site.ll <- me.sdm.site
coord.site$rf.site.app.ll <- rf.sdm.app.site
coord.site$me.site.app.ll <- me.sdm.app.site

coord.site <- coord.site[,c(1,4:ncol(coord.site))]
names(coord.site) <- c("full.garden", "site.lat", "site.long", "rf.site.ll", "me.site.ll", "rf.site.app.ll", "me.site.app.ll")

### MERGE EVERYTHING TOGETHER
het <- het %>% mutate(pop = if_else(pop=="EGG", "Egg", pop))
cg.het <- merge(cg.no_herb, het, by="pop") %>% dplyr::select(-c(index, g.index, X.x, X.y))
cg.het <- merge(cg.het, coord.pop, by="pop")
cg.het <- merge(cg.het, coord.site, by="full.garden")

write.csv(cg.het, "LA_HET/data/elev_LA_HET_NOHERB_rfme_ext_fix1.csv", row.names = F)
