
# data prep
library(dplyr)

# modeling
library(glmmTMB)
library(DHARMa)
library(car)

# plotting
library(visreg)
library(sjPlot)
library(ggplot2)

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/LA_HET/")

### Read in data formatted by lat_het_la_dat_prep.R
dat <- read.csv("data/lat_LA_HET_NOHERB_rfme_ext_fix1.csv", stringsAsFactors = T)
dat <- dat %>% 
  mutate_at(vars(block), as.factor) %>% 
  as.data.frame()
dat$delta.me.ll <- dat$me.site.ll-dat$me.pop.ll
dat$delta.rf.ll <- dat$rf.site.ll-dat$rf.pop.ll

dat$pc.dist <- sqrt((dat$delta_pc1)^2 + (dat$delta_pc2)^2)

### Read in data formatted by lat_het_la_dat_prep.R
dat2 <- read.csv("data/elev_LA_HET_NOHERB_rfme_ext_fix1.csv", stringsAsFactors = T)

dat2 <- dat2 %>% 
  mutate_at(vars(block), as.factor) %>% 
  filter(life_het >= -0.05) %>% # 5% lower fitness of hybrids than mid-parent
  as.data.frame()

dat2$delta.me.ll <- dat2$me.site.ll-dat2$me.pop.ll
dat2$delta.rf.ll <- dat2$rf.site.ll-dat2$rf.pop.ll

dat2$pc.dist <- sqrt((dat2$delta_pc1)^2 + (dat2$delta_pc2)^2)

### combine
dat <- dat %>% dplyr::select(rel.fr, life_het, me.site.lat.ll, me.site.ll, rf.site.lat.ll, rf.site.ll, full.garden, block, pop, fl.fr, home)
dat2 <- dat2 %>% dplyr::select(rel.fr, life_het, me.site.app.ll, me.site.ll, rf.site.app.ll, rf.site.ll, full.garden, block, pop, fl.fr, home)
names(dat)[c(3,5)] <- c("me.site.local.ll","rf.site.local.ll")
names(dat2)[c(3,5)] <- c("me.site.local.ll", "rf.site.local.ll")
dat$grad <- "shallow"
dat2$grad <- "steep"

dat.tot <- rbind(dat,dat2)
dat.tot$grad <- as.factor(dat.tot$grad)

dat.tot <- dat.tot %>% 
  group_by(grad) %>% 
  mutate(local.me.ll.scale = scale(me.site.local.ll),
         global.me.ll.scale = scale(me.site.ll),
         local.rf.ll.scale = scale(rf.site.local.ll),
         global.rf.ll.scale = scale(rf.site.ll),
         fl.fr.sqrt = sqrt(fl.fr))


### local sdm models
### MAXENT
{
  het.site.loc.lm.freq.scale <- glmmTMB(rel.fr ~ life_het * poly(local.me.ll.scale,1) * grad + (1|full.garden/block), 
                                  family = tweedie(link = "log"), data = dat.tot)
  Anova(het.site.loc.lm.freq.scale)

  het.site.loc.lm.freq.scale.onemod <- glmmTMB(rel.fr ~ life_het * poly(local.me.ll.scale,1) + (1|full.garden/block), 
                                        family = tweedie(link = "log"), data = dat.tot)
  Anova(het.site.loc.lm.freq.scale.onemod)

  het.site.loc.lm.freq.scale.onemod.grad <- glmmTMB(rel.fr ~ life_het * poly(local.me.ll.scale,1) + grad + (1|full.garden/block),
                                               family = tweedie(link = "log"), data = dat.tot)
  Anova(het.site.loc.lm.freq.scale.onemod.grad)
}

### RANDOM FOREST
{
  het.site.loc.lm.freq.scale.gradint <- glmmTMB(rel.fr ~ life_het * local.rf.ll.scale * grad + (1|full.garden/block),
                                        family = tweedie(link = "log"), data = dat.tot)
  Anova(het.site.loc.lm.freq.scale.gradint)
  
  het.site.loc.lm.freq.scale.nograd <- glmmTMB(rel.fr ~ life_het * local.rf.ll.scale + (1|full.garden/block), 
                                               family = tweedie(link = "log"), data = dat.tot)
  Anova(het.site.loc.lm.freq.scale.nograd)
  
  het.site.loc.lm.freq.scale.fixgrad <- glmmTMB(rel.fr ~ life_het * local.rf.ll.scale + grad + (1|full.garden/block),
                                                    family = tweedie(link = "log"), data = dat.tot)
  Anova(het.site.loc.lm.freq.scale.fixgrad)
}

anova(het.site.loc.lm.freq.scale.gradint, het.site.loc.lm.freq.scale.nograd, test = "Chisq")
BIC(het.site.loc.lm.freq.scale.gradint) - BIC(het.site.loc.lm.freq.scale.nograd)

anova(het.site.loc.lm.freq.scale.fixgrad, het.site.loc.lm.freq.scale.nograd, test = "Chisq")
BIC(het.site.loc.lm.freq.scale.fixgrad) - BIC(het.site.loc.lm.freq.scale.nograd)

summary(het.site.loc.lm.freq.scale.nograd)
performance::performance(het.site.loc.lm.freq.scale.nograd)
simulationOutput <- DHARMa::simulateResiduals(het.site.loc.lm.freq.scale.nograd, plot = T)

### plot
visreg2d(het.site.loc.lm.freq.scale.nograd, "local.rf.ll.scale", "life_het", plot.type="image", type=c("conditional"), 
         scale="response", xlab="Scaled Habitat Suitability", ylab="Genetic Load", plot=T)

# new_data.scale <- expand.grid(
#   life_het = seq(min(dat.tot$life_het)-0.025, max(dat.tot$life_het)+0.025, by=0.01),
#   local.rf.ll.scale = seq(min(dat.tot$local.rf.ll.scale)-0.075, max(dat.tot$local.rf.ll.scale)+0.1, by=0.05),
#   full.garden = unique(dat.tot$full.garden),
#   pop = unique(dat.tot$pop),
#   block = unique(dat.tot$block))
# new_data.scale$predicted <- predict(het.site.loc.lm.freq.scale.nograd, newdata = new_data.scale, type="response")

new_data.scale <- expand.grid(
  life_het = seq(min(dat.tot$life_het)-0.025, max(dat.tot$life_het)+0.025, by=0.01),
  local.rf.ll.scale = seq(min(dat.tot$local.rf.ll.scale)-0.075, max(dat.tot$local.rf.ll.scale)+0.1, by=0.05))
new_data.scale$predicted <- predict(het.site.loc.lm.freq.scale.nograd, newdata = new_data.scale, 
                                    type="response", re.form = NA, allow.new.levels = TRUE)



# # no longer major ugo
# new_data.scale.2 <- new_data.scale %>% 
#   mutate(life_het = round(life_het, 2), local.rf.ll.scale = round(local.rf.ll.scale/0.1)*0.1) %>% 
#   group_by(pop, full.garden) %>% 
#   reframe(predicted=mean(predicted, na.rm=T))
# mid <- quantile(new_data.scale$predicted)[3] %>% as.numeric()

new_data.scale.2 <- new_data.scale %>% 
  mutate(life_het = round(life_het, 2), local.rf.ll.scale = round(local.rf.ll.scale/0.1)*0.1) %>% 
  group_by(life_het, local.rf.ll.scale) %>% 
  reframe(predicted=mean(predicted, na.rm=T))
mid <- quantile(new_data.scale$predicted)[3] %>% as.numeric()

n.contour <- 12
la.het.plot <- 
ggplot() +
  # geom_tile(new_data.scale.2, aes(x=local.rf.ll.scale, y=life_het, fill=log10(predicted))) +
  geom_contour_filled(data=new_data.scale.2 %>% ungroup(), 
                      aes(x=local.rf.ll.scale, y=life_het, z=log10(predicted)), bins=n.contour) + # log10(predicted)
  geom_contour(data=new_data.scale.2 %>% ungroup(), 
               aes(x=local.rf.ll.scale, y=life_het, z=log10(predicted)),
               bins=n.contour, size=0.25, color="white") +
  geom_contour(data=new_data.scale.2 %>% ungroup(), 
               aes(x=local.rf.ll.scale, y=life_het, z=log10(predicted)),
               bins=n.contour, size=0.1, color="black") +
  geom_point(data=dat.tot %>% group_by(pop, full.garden) %>% slice_head(n=1), 
             aes(x=local.rf.ll.scale, y=life_het), shape=21, size=5, alpha=0.5, color="gray45")+
  # facet_wrap(~grad) +
  # scale_color_gradient2(low="blue", mid="white", high="red", midpoint = -0.24)+
  # scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = log10(mid))+
  scale_fill_manual(values=(colorRampPalette(c("tomato1", "white", "steelblue1"))(n.contour)),
                    labels=sprintf("%.1f - %.1f", seq(0,max(new_data.scale.2$predicted), 
                                                      by=(max(new_data.scale.2$predicted)/n.contour))[-(n.contour+1)], 
                                   seq(0,max(new_data.scale.2$predicted), 
                                       by=(max(new_data.scale.2$predicted)/n.contour))[-1]))+
  ylim(c(min(dat.tot$life_het)-0.02, max(dat.tot$life_het)+0.02))+
  # xlim(c(min(dat.tot$local.rf.ll.scale)-0.1, max(dat.tot$local.rf.ll.scale)+0.1))+
  ylab("Genetic Load")+
  xlab("Scaled Garden Site Suitability")+
  labs(fill="Relative Fitness")+
  theme_bw()+
  theme(text = element_text(size=15), panel.grid.major = element_blank(), aspect.ratio=1,
        panel.grid.minor = element_blank(), plot.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black"), legend.position = "right")+
  guides(fill = guide_legend(reverse = TRUE))
la.het.plot

jpeg("plots/ensemble_scale_fix1_points_xlab.jpeg", units="in", height=8, width=8, res=300)
la.het.plot
dev.off()






### Site quality
lat.dat <- read.csv("~/Desktop/Documents/Research/Paper Code/adaptation_load/LA_HET/data/lat_LA_HET_NOHERB_rfme_ext_fix1.csv", stringsAsFactors = T)
app.dat <- read.csv("~/Desktop/Documents/Research/Paper Code/adaptation_load/LA_HET/data/elev_LA_HET_NOHERB_rfme_ext_fix1.csv", stringsAsFactors = T)

lat.dat <- lat.dat %>% 
  dplyr::select(garden.group, full.garden, pop, pop.group, block, new.ID, fl.fr, live, rel.fr, rf.site.lat.ll, live, home) %>% 
  dplyr::rename("rf.site.ll"="rf.site.lat.ll") %>% 
  mutate(rf.site.ll.scale=scale(rf.site.ll))
lat.dat$grad <- "shallow"
app.dat <- app.dat %>% 
  dplyr::select(garden.group, full.garden, pop, pop.group, block, new.ID, fl.fr, live, rel.fr, rf.site.app.ll, live, home) %>% 
  dplyr::rename("rf.site.ll"="rf.site.app.ll") %>% 
  mutate(rf.site.ll.scale=scale(rf.site.ll))
app.dat$grad <- "steep"

dat <- rbind(lat.dat, app.dat)
head(dat)

dat.flfr <- dat %>% 
  ungroup() %>% 
  group_by(pop) %>% 
  mutate(pop.flfr.mu = mean(fl.fr[home == 0], na.rm = TRUE)) %>%
  group_by(full.garden, pop) %>% 
  mutate(pop.flfr.norm = mean(fl.fr)/pop.flfr.mu) %>% 
  ungroup() %>%
  group_by(garden.group, full.garden, rf.site.ll, rf.site.ll.scale, grad) %>%
  reframe(garden.quality = mean(pop.flfr.norm))

names(dat.flfr)[4] <- "rf.site.ll.scale"
dat.flfr$rf.site.ll.scale <- as.numeric(dat.flfr$rf.site.ll.scale)
dat.flfr$grad <- as.factor(dat.flfr$grad)

dat.flfr$log10.garden.quality <- log10(dat.flfr$garden.quality)
# qual.lm <- lm(garden.quality~rf.site.ll.scale, data=dat.flfr)
qual.lm <- lm(log10.garden.quality~rf.site.ll+grad, data=dat.flfr)
Anova(qual.lm, type=3)

eff <- ggeffects::ggpredict(qual.lm, terms=c("rf.site.ll", "grad")) %>% as.data.frame()

jpeg("plots/site_quality_suitability.jpeg", width=6, height=5, units="in", res=5e2)
ggplot(data=eff %>% 
         mutate(x = if_else(group=="shallow" & x <= ((dat.flfr %>% filter(grad=="shallow"))$rf.site.ll %>% min())-0.05, NA, x)) %>% 
         mutate(x = if_else(group=="steep" & x >= ((dat.flfr %>% filter(grad=="steep"))$rf.site.ll %>% max())+0.01, NA, x)), 
       aes(x=x, y=10^predicted, color=group))+
  geom_point(data=dat.flfr, aes(x=rf.site.ll, y=10^log10.garden.quality, color=grad, shape=grad), size=2)+
  geom_line()+
  geom_ribbon(aes(ymin = 10^conf.low, ymax = 10^conf.high),
              alpha=0.05, linetype=2, linewidth=0.1)+
  theme_bw()+
  xlab("Predicted Garden Suitability")+
  scale_color_manual(values=c("black", "black"))+
  scale_shape_manual(values=c(21,24))+
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits=c(0.0001,100))+
  ylab("Garden Site Quality")+
  theme(text = element_text(size=15), panel.grid.major = element_blank(), aspect.ratio=1,
        panel.grid.minor = element_blank(), plot.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black"), legend.position = "right")+
  guides(fill = guide_legend(reversxe = TRUE))+
  labs(color="gradient", shape="gradient")
dev.off()



### Site quality (a la Bontrager et al., 2021) and genetic load
dat.qual <- dat.tot %>% 
  ungroup() %>% 
  group_by(pop) %>% 
  mutate(pop.flfr.mu = mean(fl.fr[home == 0], na.rm = TRUE)) %>%
  group_by(full.garden, pop) %>% 
  mutate(pop.flfr.norm = mean(fl.fr)/pop.flfr.mu) %>% 
  ungroup() %>%
  group_by(full.garden, rf.site.local.ll, local.rf.ll.scale, grad) %>%
  mutate(garden.quality = mean(pop.flfr.norm)) %>% 
  mutate(log10.garden.quality = log10(garden.quality)) %>% 
  group_by(grad) %>% 
  mutate(log10.garden.quality.scale = scale(log10.garden.quality)) %>% 
  ungroup()

qual.mod.nograd <- glmmTMB(rel.fr ~ life_het * log10.garden.quality.scale + (1|full.garden/block), family = tweedie(link = "log"), 
                           data = dat.qual)
Anova(qual.mod.nograd)

new_data.scale <- expand.grid(
  life_het = seq(min(dat.qual$life_het)-0.025, max(dat.qual$life_het)+0.025, by=0.01),
  log10.garden.quality.scale = seq(min(dat.qual$log10.garden.quality.scale)-0.075, max(dat.qual$log10.garden.quality.scale)+0.1, by=0.05))
new_data.scale$predicted <- predict(qual.mod.nograd, newdata = new_data.scale, 
                                    type="response", re.form = NA, allow.new.levels = TRUE)

new_data.scale.2 <- new_data.scale %>% 
  mutate(life_het = round(life_het, 2), log10.garden.quality.scale = round(log10.garden.quality.scale/0.1)*0.1) %>% 
  group_by(life_het, log10.garden.quality.scale) %>% 
  reframe(predicted=mean(predicted, na.rm=T))
mid <- quantile(new_data.scale$predicted)[3] %>% as.numeric()

n.contour <- 12
la.het.qual.plot <- 
  ggplot() +
  # geom_tile(new_data.scale.2, aes(x=local.rf.ll.scale, y=life_het, fill=log10(predicted))) +
  geom_contour_filled(data=new_data.scale.2 %>% ungroup(), 
                      aes(x=log10.garden.quality.scale, y=life_het, z=log10(predicted)), bins=n.contour) + # log10(predicted)
  geom_contour(data=new_data.scale.2 %>% ungroup(), 
               aes(x=log10.garden.quality.scale, y=life_het, z=log10(predicted)),
               bins=n.contour, size=0.25, color="white") +
  geom_contour(data=new_data.scale.2 %>% ungroup(), 
               aes(x=log10.garden.quality.scale, y=life_het, z=log10(predicted)),
               bins=n.contour, size=0.1, color="black") +
  geom_point(data=dat.qual %>% group_by(pop, full.garden) %>% slice_head(n=1), 
             aes(x=log10.garden.quality.scale, y=life_het), shape=21, size=5, alpha=0.5, color="gray45")+
  # facet_wrap(~grad) +
  # scale_color_gradient2(low="blue", mid="white", high="red", midpoint = -0.24)+
  # scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = log10(mid))+
  scale_fill_manual(values=(colorRampPalette(c("tomato1", "white", "steelblue1"))(n.contour)),
                    labels=sprintf("%.1f - %.1f", seq(0,max(new_data.scale.2$predicted), 
                                                      by=(max(new_data.scale.2$predicted)/n.contour))[-(n.contour+1)], 
                                   seq(0,max(new_data.scale.2$predicted), 
                                       by=(max(new_data.scale.2$predicted)/n.contour))[-1]))+
  ylim(c(min(dat.tot$life_het)-0.02, max(dat.tot$life_het)+0.02))+
  # xlim(c(min(dat.tot$local.rf.ll.scale)-0.1, max(dat.tot$local.rf.ll.scale)+0.1))+
  ylab("Genetic Load")+
  xlab("Scaled Garden Quality")+
  labs(fill="Relative Fitness")+
  theme_bw()+
  theme(text = element_text(size=15), panel.grid.major = element_blank(), aspect.ratio=1,
        panel.grid.minor = element_blank(), plot.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black"), legend.position = "right")+
  guides(fill = guide_legend(reverse = TRUE))
la.het.qual.plot

