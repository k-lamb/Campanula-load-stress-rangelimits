library(ggnewscale)
library(raster)
library(sp)
library(rnaturalearth) # setting up coordinate grid for prediction
library(ggpubr)
library(ggplot2)
setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/SDM/")

# # no longer true:
# requires pred.grid objects (pred.grid.app and pred.grid.exp) from either modeling script (maxent_model.R or rf_model.R)

model.type = "rf" # "maxent" or "rf"
region = "ALL" # expansion | Appalachian
buffer.dist = 100 # lat/long buffer on world crop
connectivity.pred.cutoff <- 0.5

if (model.type == "rf") {
  # pred.grid.app <- pred.grid.app.rf
  # pred.grid.exp <- pred.grid.exp.rf
  pred.grid.app <- read.csv("data/models/rf/local/m_rf_100km_Appalachian_extended_fix1.csv")
  pred.grid.exp <- read.csv("data/models/rf/local/m_rf_100km_expansion_extended_fix1.csv")
}

# trying to pull from existing data frames (.csv's) below instead of regenerating each time
if (model.type == "maxent") {
  # pred.grid.app <- pred.grid.app.me
  # pred.grid.exp <- pred.grid.exp.me
  pred.grid.app <- read.csv("data/models/maxent/local/m_me_100km_Appalachian_extended_fix1.csv")
  pred.grid.exp <- read.csv("data/models/maxent/local/m_me_100km_expansion_extended_fix1.csv")
}

# for APP:
xy <- c(-81, 37) # in NC... should work for both
app.r <- rasterFromXYZ(pred.grid.app %>% dplyr::select(longitude, latitude, pred.cell) %>% filter(pred.cell >= connectivity.pred.cutoff))
cell_coords <- cellFromXY(app.r, xy)
connected.app <- clump(app.r)
region_number <- connected.app[cell_coords]
mask <- connected.app == region_number
app.r[!mask] <- NA  # Set unconnected cells to NA or modify them as needed
# plot(app.r)

pred.grid.app2 <- as.data.frame(app.r, xy=TRUE) %>% filter(is.na(pred.cell) == F)
pred.grid.app2 <- pred.grid.app2 %>% mutate(bins = cut(pred.cell, breaks = c(-Inf,seq(0,1,0.1),Inf))) %>% as.data.frame()

# for WEST:
xy <- c(-87, 40) # in NC... should work for both
exp.r <- rasterFromXYZ(pred.grid.exp %>% dplyr::select(longitude, latitude, pred.cell) %>% filter(pred.cell >= connectivity.pred.cutoff))
cell_coords <- cellFromXY(exp.r, xy)
connected.exp <- clump(exp.r)
region_number <- connected.exp[cell_coords]
mask <- connected.exp == region_number
exp.r[!mask] <- NA  # Set unconnected cells to NA or modify them as needed
# plot(exp.r)

pred.grid.exp2 <- as.data.frame(exp.r, xy=TRUE) %>% filter(is.na(pred.cell) == F)
pred.grid.exp2 <- pred.grid.exp2 %>% mutate(bins = cut(pred.cell, breaks = c(-Inf,seq(0,1,0.1),Inf))) %>% as.data.frame()
world.sf <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))
t.size=20 # was 20


p.con <- 
  ggplot()+
    geom_sf(data=world.sf, fill=NA)+
    geom_sf(data=states, fill=NA)+ #adds state boundaries
    geom_tile(data=pred.grid.exp2, aes(x=x, y=y, fill=bins), size=0.2)+
    scale_fill_manual(values = c("azure2", "plum1", "violet","orchid", "orchid3", "orchid4", "darkorchid4", "purple4"))+
    ggnewscale::new_scale_fill()+
    geom_tile(data=pred.grid.app2, aes(x=x, y=y, fill=bins), size=1, alpha=1)+ #for crop
    scale_fill_manual(values = c("azure2", "darkolivegreen1", "darkolivegreen2","darkolivegreen3", "chartreuse3", "seagreen", "green4", "darkgreen"))+
    geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
    # scale_color_manual(values = c("steelblue4",  "yellow3", "yellow", "chartreuse3", "green4"))+
    geom_sf(data=lakes, fill="white")+ #adds the great lakes
    coord_sf(xlim=c(-100,-72), ylim=c(32,48), expand=FALSE)+
    xlab("Longitude")+
    ylab("Latitude")+
    # geom_point(data=subset(c.americana, presence==1), aes(x=longitude, y=latitude), colour="black", size=0.25)+ #adds lab collections
    theme_bw()+
    theme(text = element_text(size=t.size), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_rect(fill="white"), axis.line = element_line(colour = "black"), legend.position = "none")

p.app <- 
  ggplot()+
  geom_sf(data=world.sf, fill=NA)+
  geom_sf(data=states, fill=NA)+ #adds state boundaries
  # geom_tile(data=pred.grid.exp2 %>% filter(pred.cell >= 0.5), aes(x=x, y=y, fill=bins), size=0.2)+
  # scale_fill_manual(values = c("plum1", "violet","orchid", "orchid3", "orchid4", "darkorchid4", "purple4"))+
  ggnewscale::new_scale_fill()+
  geom_tile(data=pred.grid.app2, aes(x=x, y=y, fill=bins), size=1, alpha=1)+ #for crop
  scale_fill_manual(values = c("azure2","darkolivegreen1", "darkolivegreen2","darkolivegreen3", "chartreuse3", "seagreen", "green4", "darkgreen"))+
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  # scale_color_manual(values = c("steelblue4",  "yellow3", "yellow", "chartreuse3", "green4"))+
  geom_sf(data=lakes, fill="white")+ #adds the great lakes
  coord_sf(xlim=c(-100,-72), ylim=c(32,48), expand=FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  # geom_point(data=subset(c.americana, presence==1), aes(x=longitude, y=latitude), colour="black", size=0.25)+ #adds lab collections
  theme_bw()+
  theme(text = element_text(size=t.size), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_rect(fill="white"), axis.line = element_line(colour = "black"), legend.position = "none")

p.west <-
  ggplot()+
  geom_sf(data=world.sf, fill=NA)+
  geom_sf(data=states, fill=NA)+ #adds state boundaries
  # geom_tile(data=pred.grid.exp2 %>% filter(pred.cell >= connectivity.pred.cutoff), aes(x=x, y=y, fill=bins), size=0.2)+
  geom_tile(data=pred.grid.exp2, aes(x=x, y=y, fill=bins), size=0.2)+
  scale_fill_manual(values = c("azure2","plum1", "violet","orchid", "orchid3", "orchid4", "darkorchid4", "purple4"))+
  ggnewscale::new_scale_fill()+
  # geom_tile(data=pred.grid.app2, aes(x=x, y=y, fill=bins), size=1, alpha=1)+ #for crop
  # scale_fill_manual(values = c("darkolivegreen1", "darkolivegreen2","darkolivegreen3", "chartreuse3", "seagreen", "green4", "darkgreen"))+
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  # scale_color_manual(values = c("steelblue4",  "yellow3", "yellow", "chartreuse3", "green4"))+
  geom_sf(data=lakes, fill="white")+ #adds the great lakes
  coord_sf(xlim=c(-100,-72), ylim=c(32,48), expand=FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  # geom_point(data=subset(c.americana, presence==1), aes(x=longitude, y=latitude), colour="black", size=0.25)+ #adds lab collections
  theme_bw()+
  theme(text = element_text(size=t.size), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_rect(fill="white"), axis.line = element_line(colour = "black"), legend.position = "none")


p.agg <- ggarrange(p.west, p.app, ncol=1)
p.agg <- ggarrange(p.con, p.agg, ncol=2)

p.westapp <- ggarrange(p.west, p.app, ncol=2)

if (model.type == "maxent") {
  jpeg(sprintf("plots/maxent/local/connected_cells/m_me_ALL_50percpred_connectedcellsONLY_%s_%skm_AGG_extended_fix1.jpeg", region, buffer.dist), width=30, height=10, units="in", res=500)
  print(p.agg)
  dev.off()
  
  jpeg(sprintf("plots/maxent/local/connected_cells/m_me_ALL_50percpred_connectedcellsONLY_%s_%skm_noOverlapPanel_extended_fix1.jpeg", region, buffer.dist), width=30, height=15, units="in", res=500)
  print(p.westapp)
  dev.off()
  
  write.csv(pred.grid.app2, "data/connectivity_me_local_app_extended_fix1.csv")
  write.csv(pred.grid.exp2, "data/connectivity_me_local_lat_extended_fix1.csv")
  
}

if (model.type == "rf") {
  jpeg(sprintf("plots/rf/local/connected_cells/m_rf_ALL_50percpred_connectedcellsONLY_%s_%skm_AGG_extended_fix1.jpeg", region, buffer.dist), width=30, height=15, units="in", res=500)
  print(p.agg)
  dev.off()
  
  jpeg(sprintf("plots/rf/local/connected_cells/m_rf_ALL_50percpred_connectedcellsONLY_%s_%skm_noOverlapPanel_extended_fix1.jpeg", region, buffer.dist), width=30, height=15, units="in", res=500)
  print(p.westapp)
  dev.off()
  
  write.csv(pred.grid.app2, "data/connectivity_rf_local_app_extended_fix1.csv")
  write.csv(pred.grid.exp2, "data/connectivity_rf_local_lat_extended_fix1.csv")
}

