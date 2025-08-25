
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggbreak)

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/LA_HET/")

elev_dat <- read.csv("data/elev_LA_HET_NOHERB_rfme_ext_fix1.csv")
lat_dat <- read.csv("data/lat_LA_HET_NOHERB_rfme_ext_fix1.csv")

local.global = "local" # global or local -- global now deprecated

### latitude ####
if (local.global == "local") {
  dat.pop <- lat_dat %>% 
    group_by(pop.group) %>% 
    reframe(mean.pop=mean(rf.pop.lat.ll), se.pop=1.96*(sd(rf.pop.lat.ll)/sqrt(n()))) %>% 
    mutate(lower.pop=mean.pop-se.pop, upper.pop=mean.pop+se.pop) %>% 
    mutate(pop.group = case_when(pop.group == "low" ~ "Core",
                                 pop.group == "high" ~ "Edge",
                                 pop.group == "mid" ~ "Mid"))
  dat.lat.cg <- lat_dat %>% 
    group_by(garden.group) %>% 
    reframe(mean=mean(rf.site.lat.ll), se=1.96*(sd(rf.site.lat.ll)/sqrt(n()))) %>% 
    mutate(lower=mean-se, upper=mean+se) %>% 
    mutate(garden.group = case_when(garden.group == "low" ~ "Core",
                                    garden.group == "high" ~ "Edge",
                                    garden.group == "mid" ~ "Mid"))
}

if (local.global == "global") {
  dat.pop <- lat_dat %>% 
    group_by(pop.group) %>% 
    reframe(mean.pop=mean(rf.pop.ll), se.pop=1.96*(sd(rf.pop.ll)/sqrt(n()))) %>% 
    mutate(lower.pop=mean.pop-se.pop, upper.pop=mean.pop+se.pop) %>% 
    mutate(pop.group = case_when(pop.group == "low" ~ "Core",
                                 pop.group == "high" ~ "Edge",
                                 pop.group == "mid" ~ "Mid"))
  dat.lat.cg <- lat_dat %>% 
    group_by(garden.group) %>% 
    reframe(mean=mean(rf.site.ll), se=1.96*(sd(rf.site.ll)/sqrt(n()))) %>% 
    mutate(lower=mean-se, upper=mean+se) %>% 
    mutate(garden.group = case_when(garden.group == "low" ~ "Core",
                                    garden.group == "high" ~ "Edge",
                                    garden.group == "mid" ~ "Mid"))
}


dat.pop$pop.group <- factor(dat.pop$pop.group, levels = c("Core", "Mid", "Edge"))
dat.lat.cg$garden.group <- factor(dat.lat.cg$garden.group, levels = c("Core", "Mid", "Edge"))

lat_gg <- ggplot()+
  geom_point(data=dat.pop, aes(x=pop.group, y=mean.pop), shape=15, size=2)+
  geom_line(data=dat.pop, aes(x=pop.group, y=mean.pop, group=NA))+
  geom_errorbar(data=dat.pop, aes(x=pop.group, y=mean.pop, ymin=lower.pop, ymax=upper.pop), width=0.1)+
  geom_errorbar(data=dat.lat.cg, aes(x=garden.group, y=mean, ymin=lower, ymax=upper), width=0.1, color="blue")+
  geom_point(data=dat.lat.cg, aes(x=garden.group, y=mean), shape=15, size=2, color="blue")+
  geom_line(data=dat.lat.cg, aes(x=garden.group, y=mean, group=NA), color="blue", linetype=2)+
  theme_bw()+
  ylim(c(0.5, 1))+
  ylab("suitability")+
  xlab("population range position")+
  theme(text = element_text(size=12), panel.grid.major = element_blank(), aspect.ratio=2,
        panel.grid.minor = element_blank(), plot.background = element_rect(fill="white"), 
        axis.line = element_line(colour = "black"))


### elevation ###
if (local.global == "local") {
  dat <- elev_dat %>% 
    group_by(pop.group) %>% 
    reframe(mean=mean(rf.pop.app.ll), se=1.96*(sd(rf.pop.app.ll)/sqrt(n()))) %>% 
    mutate(lower=mean-se, upper=mean+se) %>% 
    mutate(pop.group = case_when(pop.group == "low" ~ "Core",
                                 pop.group == "high" ~ "Edge"))
  dat.cg <- elev_dat %>% 
    group_by(garden.group) %>% 
    reframe(mean=mean(rf.site.app.ll), se=1.96*(sd(rf.site.app.ll)/sqrt(n()))) %>% 
    mutate(lower=mean-se, upper=mean+se) %>% 
    mutate(garden.group = case_when(garden.group == "low" ~ "Core",
                                    garden.group == "high" ~ "Edge"))
}

if (local.global == "global") {
  dat <- elev_dat %>% 
    group_by(pop.group) %>% 
    reframe(mean=mean(me.pop.ll), se=1.96*(sd(me.pop.ll)/sqrt(n()))) %>% 
    mutate(lower=mean-se, upper=mean+se) %>% 
    mutate(pop.group = case_when(pop.group == "low" ~ "Core",
                                 pop.group == "high" ~ "Edge"))
  dat.cg <- elev_dat %>% 
    group_by(garden.group) %>% 
    reframe(mean=mean(me.site.ll), se=1.96*(sd(me.site.ll)/sqrt(n()))) %>% 
    mutate(lower=mean-se, upper=mean+se) %>% 
    mutate(garden.group = case_when(garden.group == "low" ~ "Core",
                                    garden.group == "high" ~ "Edge"))
}

dat$pop.group <- factor(dat$pop.group, levels = c("Core", "Edge"))
dat.cg$garden.group <- factor(dat.cg$garden.group, levels = c("Core", "Edge"))


elev_gg <- ggplot()+
  geom_point(data=dat, aes(x=pop.group, y=mean), shape=15, size=2)+
  geom_line(data=dat, aes(x=pop.group, y=mean, group = 1))+
  geom_errorbar(data=dat, aes(x=pop.group, y=mean, ymin=lower, ymax=upper), width=0.1)+
  geom_point(data=dat.cg, aes(x=garden.group, y=mean), shape=15, size=2, color="blue")+
  geom_line(data=dat.cg, aes(x=garden.group, y=mean, group=NA), color="blue", linetype=2)+
  geom_errorbar(data=dat.cg, aes(x=garden.group, y=mean, ymin=lower, ymax=upper), width=0.1, color="blue")+
  theme_bw()+
  ylim(c(0.5, 1))+
  ylab("")+
  xlab("population range position")+
  theme(text = element_text(size=12), panel.grid.major = element_blank(), aspect.ratio=2,
        panel.grid.minor = element_blank(), plot.background = element_rect(fill="white"), 
        axis.line = element_line(colour = "black"))
elev_gg

### plotted together ###
ggarrange(lat_gg, elev_gg, ncol=2, align="hv")


s.size = 4.5
west.color = "purple2"
mtn.color = "green4"

temp <- 
ggplot()+
  # latitude
  # geom_point(data=dat.pop, aes(x=pop.group, y=mean.pop), shape=16, size=s.size, color="black")+
  geom_point(data=dat.pop, aes(x=pop.group, y=mean.pop), shape=21, size=s.size-0.5, color="black")+
  geom_line(data=dat.pop, aes(x=pop.group, y=mean.pop, group=NA), color="black")+
  geom_errorbar(data=dat.pop, aes(x=pop.group, y=mean.pop, ymin=lower.pop, ymax=upper.pop), width=0.1, color="black")+
  
  # geom_point(data=dat.lat.cg, aes(x=garden.group, y=mean), shape=16, size=s.size, color="black")+
  geom_point(data=dat.lat.cg, aes(x=garden.group, y=mean), shape=21, size=s.size-0.5, color="black")+
  geom_line(data=dat.lat.cg, aes(x=garden.group, y=mean, group=NA), color="black", linetype=2)+
  geom_errorbar(data=dat.lat.cg, aes(x=garden.group, y=mean, ymin=lower, ymax=upper), width=0.1, color="black")+
  # elevation
  # geom_point(data=dat, aes(x=pop.group, y=mean+0.01), shape=17, size=s.size, color="black")+
  geom_point(data=dat, aes(x=pop.group, y=mean+0.01), shape=24, size=s.size-0.5, color="black",
             position = position_dodge2(width = 0.1))+
  geom_line(data=dat, aes(x=pop.group, y=mean+0.01, group=NA), color="black", position = position_dodge2(width = 0.1))+
  geom_errorbar(data=dat, aes(x=pop.group, y=mean+0.01, ymin=lower+0.01, ymax=upper+0.01, group=NA), 
                width=0.1, color="black", position = position_dodge2(width = 0.1))+
  
  # geom_point(data=dat.cg, aes(x=garden.group, y=mean), shape=17, size=s.size, color="black")+
  geom_point(data=dat.cg, aes(x=garden.group, y=mean), shape=24, size=s.size-0.5, color="black")+
  geom_line(data=dat.cg, aes(x=garden.group, y=mean, group=NA), color="black", linetype=2)+
  geom_errorbar(data=dat.cg, aes(x=garden.group, y=mean, ymin=lower, ymax=upper), width=0.1, color="black")+
  
  theme_bw()+
  # ylim(c(0.15, 0.75))+
  scale_y_continuous(breaks=seq(0.0,1,by=0.1), limits=c(0.7,1.01))+
  ylab("Suitability")+
  xlab("Population Range Position")+
  theme(text = element_text(size=11), panel.grid.major = element_blank(), aspect.ratio=2,
        panel.grid.minor = element_blank(), plot.background = element_rect(fill="white"), 
        axis.line = element_line(colour = "black"), legend.position = "right")
temp

if (local.global == "local") {
  loc <- temp
}
if (local.global == "global") {
  glob <- temp
}


jpeg("plots/pop_suitability_appextended_rf_fix1.jpeg", res=500, units="in", height=3.9, width=5)
temp
dev.off()
