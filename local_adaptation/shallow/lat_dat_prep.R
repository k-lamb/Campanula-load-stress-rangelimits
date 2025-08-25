### effect size Leimu and Fischer 2008 local adaptation test ###

# libraries
library(readxl)
library(dplyr)
library(plyr)
library(ggplot2)
library(devtools)

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/local_adaptation/")

### PIPELINE
# 1. DATA IMPORT AND PREPARATION
#    - convert strings to factors and numerics
#    - calculate relativized metrics to control for site conditions
#    - filter data sets to control for herbivory issues
#    - calculate Hedge's D for HERB and NO HERB

### 1. DATA IMPORT AND PREPARATION
cg <- read_excel("data/latitude/western_cg_2023.xlsx") %>% as.data.frame()

cg.dat <- read.csv("data/CG_locations.csv")
names(cg.dat)[c(1,3)] <- c("full.garden","garden.position")

# convert to factors and fix NA calling issue
cg[cg == "NA"] <- NA # convert character NA to NA class
fac_list <- c("pop", "pop.group", "label.ID", "block", "plant", "new.ID", "garden", "full.garden", 
              "garden.group", "stage", "home_group", "pop.transect", "garden.transect")
num_list<- c("live", "herb", "fl.fr", "buds", "home", "stem.d")

cg <- cg %>% 
  mutate_at(fac_list, as.factor) %>%
  mutate_at(num_list, as.numeric) %>%
  mutate(theo.repro = fl.fr+buds) %>% # creates var with theoretical maximum reproduction in the given season
  group_by(full.garden) %>% 
  mutate(rel.stem = stem.d / mean(stem.d, na.rm=T)) %>% # creates relative variables within gardens to control for site conditions
  ungroup() %>% 
  filter(is.na(live) == F) %>% 
  as.data.frame()

# dealing with herbivory issues
cg.herb <- cg %>% 
  filter(herb == 0) %>% # remove herbivorized plants
  group_by(full.garden) %>%
  mutate(rel.fr = fl.fr / mean(fl.fr, na.rm=T), # now we can calculate other relative metrics
         rel.theo = theo.repro / mean(theo.repro, na.rm=T)) %>%
  ungroup() %>% 
  as.data.frame()

# # write data
write.csv(cg.herb, "data/latitude/lat_rel_abs_NO_HERB.csv")
write.csv(cg, "data/latitude/lat_rel_abs_HERB.csv")



# HEDGE'S D FOR NO HERBIVORY (LOCAL V FOREIGN)
# calculate hedges D individual (by garden) adaptation
mean_home <- cg.herb %>%
  group_by(full.garden) %>%
  subset(home == 1) %>%
  dplyr::summarise(mean_home = mean(fl.fr, na.rm=T))

mean_foreign <- cg.herb %>%
  group_by(full.garden) %>%
  subset(home == 0) %>%
  dplyr::summarise(mean_foreign = mean(fl.fr, na.rm=T))

sd_all <- cg.herb %>%
  group_by(full.garden) %>%
  dplyr::summarise(sd_all = sd(fl.fr, na.rm=T))

temp <- merge(as.data.frame(mean_home), as.data.frame(mean_foreign), by="full.garden")
cg.herb.1 <- merge(temp, as.data.frame(sd_all), by="full.garden")
cg.herb.1$hedges_d_NOHERB <- (cg.herb.1$mean_home - cg.herb.1$mean_foreign) / cg.herb.1$sd_all

# build out data attributes and save data
fin <- cg.herb.1 %>% dplyr::select(-c(mean_home, mean_foreign))
fin <- merge(fin, cg.dat, by="full.garden")


### HOME VS AWAY NO HERBIVORY
home.away.fl <- cg.herb %>%
  group_by(pop, home) %>%
  # filter(pop.transect==garden.transect) %>% 
  dplyr::summarise(mean=mean(fl.fr))
home.away.sd <- cg.herb %>%
  group_by(pop) %>%
  # filter(pop.transect==garden.transect) %>% 
  dplyr::summarise(sd=sd(fl.fr))

home.away <- merge(home.away.sd, home.away.fl, by="pop")
away <- subset(home.away, home==0 & pop != "SK" & pop != "SM" & pop != "SH") # populations with no homes assigned
home <- subset(home.away, home==1)

hedge.d <- (home$mean-away$mean)/home$sd
home$hedge_homeaway <- (home$mean - away$mean) / home$sd

lim <- unique(cg.herb[,c("pop", "full.garden", "home")])
lim <- subset(lim, home==1)
limx <- merge(lim, home, by=c("pop", "home"))
limx <- limx[,c("pop", "full.garden", "hedge_homeaway")]

Hedge.NOHERB <- merge(limx, fin, by="full.garden")

Hedge.NOHERB <- Hedge.NOHERB %>% 
  mutate(garden.position = case_when(
    garden.position == "low" ~ "core",
    garden.position == "mid" ~ "mid",
    garden.position == "high" ~ "edge"
  )) %>% 
  dplyr::select(-sd_all)

# write
write.csv(Hedge.NOHERB, "data/latitude/lat_hedgesD_NO_HERB.csv")



# HEDGE'S D FOR HERBIVORY (LOCAL V FOREIGN)
# calculate hedges D individual (by garden) adaptation
mean_home <- cg %>%
  group_by(full.garden) %>%
  subset(home == 1) %>%
  dplyr::summarise(mean_home = mean(stem.d, na.rm=T))

mean_foreign <- cg %>%
  group_by(full.garden) %>%
  subset(home == 0) %>%
  dplyr::summarise(mean_foreign = mean(stem.d, na.rm=T))

sd_all <- cg %>%
  group_by(full.garden) %>%
  dplyr::summarise(sd_all = sd(stem.d, na.rm=T))

temp <- merge(as.data.frame(mean_home), as.data.frame(mean_foreign), by="full.garden")
cg.1 <- merge(temp, as.data.frame(sd_all), by="full.garden")
cg.1$hedges_d_HERB <- (cg.1$mean_home - cg.1$mean_foreign) / cg.1$sd_all

fin <- cg.1 %>% dplyr::select(-c(mean_home, mean_foreign))
fin <- merge(fin, cg.dat, by="full.garden")


### HOME VS AWAY HERBIVORY
home.away.fl <- cg %>%
  group_by(pop, home) %>%
  # filter(pop.transect==garden.transect) %>% 
  dplyr::summarise(mean=mean(stem.d, na.rm = T))
home.away.sd <- cg %>%
  group_by(pop) %>%
  # filter(pop.transect==garden.transect) %>% 
  dplyr::summarise(sd=sd(stem.d, na.rm = T))

home.away <- merge(home.away.sd, home.away.fl, by="pop")
away <- subset(home.away, home==0 & pop != "SK" & pop != "SM" & pop != "SH") # populations with no homes assigned
home <- subset(home.away, home==1)

hedge.d <- (home$mean-away$mean)/home$sd
home$hedge_homeaway <- (home$mean - away$mean) / home$sd

lim <- unique(cg[,c("pop", "full.garden", "home")])
lim <- subset(lim, home==1)
limx <- merge(lim, home, by=c("pop", "home"))
limx <- limx[,c("pop", "full.garden", "hedge_homeaway")]

Hedge.HERB <- merge(limx, fin, by="full.garden")

Hedge.HERB <- Hedge.HERB %>% 
  mutate(garden.position = case_when(
    garden.position == "low" ~ "core",
    garden.position == "mid" ~ "mid",
    garden.position == "high" ~ "edge"
  )) %>% 
  dplyr::select(-sd_all)

# write
write.csv(Hedge.HERB, "data/latitude/lat_hedgesD_HERB.csv")

