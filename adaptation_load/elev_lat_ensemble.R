
# data prep
library(dplyr)

# modeling
library(glmmTMB)
library(DHARMa)
library(car)

# plotting
library(visreg)
library(sjPlot)

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
dat <- dat %>% dplyr::select(rel.fr, life_het, me.site.lat.ll, me.site.ll, rf.site.lat.ll, rf.site.ll, full.garden, block, pop, fl.fr)
dat2 <- dat2 %>% dplyr::select(rel.fr, life_het, me.site.app.ll, me.site.ll, rf.site.app.ll, rf.site.ll, full.garden, block, pop, fl.fr)
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
performance(het.site.loc.lm.freq.scale.nograd)
simulationOutput <- DHARMa::simulateResiduals(het.site.loc.lm.freq.scale.nograd, plot = T)

### plot
visreg2d(het.site.loc.lm.freq.scale.nograd, "local.rf.ll.scale", "life_het", plot.type="image", type=c("conditional"), 
         scale="response", xlab="Scaled Habitat Suitability", ylab="Genetic Load", plot=T)

new_data.scale <- expand.grid(
  life_het = seq(min(dat.tot$life_het), max(dat.tot$life_het), by=0.01),
  local.rf.ll.scale = seq(min(dat.tot$local.rf.ll.scale), max(dat.tot$local.rf.ll.scale), by=0.05),
  # grad = unique(dat.tot$grad),
  full.garden = unique(dat.tot$full.garden),
  pop = unique(dat.tot$pop),
  block = unique(dat.tot$block))
new_data.scale$predicted <- predict(het.site.loc.lm.freq.scale.nograd, newdata = new_data.scale, type="response")

# no longer major ugo
new_data.scale.2 <- new_data.scale %>% 
  mutate(life_het = round(life_het, 2), local.rf.ll.scale = round(local.rf.ll.scale/0.2)*0.2) %>% 
  group_by(life_het, local.rf.ll.scale) %>% 
  reframe(predicted=mean(predicted, na.rm=T))
mid <- quantile(new_data.scale$predicted)[3] %>% as.numeric()

jpeg("plots/ensemble_scale_fix1.jpeg", units="in", height=8, width=8, res=300)
ggplot() +
  # geom_tile(new_data.scale.2, aes(x=local.rf.ll.scale, y=life_het, fill=log10(predicted))) +
  geom_contour_filled(data=new_data.scale.2 %>% ungroup(), aes(x=local.rf.ll.scale, y=life_het, 
                                                 z=log10(predicted)), bins=20) + # log10(predicted)
  geom_contour(data=new_data.scale.2 %>% ungroup(), aes(x=local.rf.ll.scale, y=life_het,
                                                        z=log10(predicted)),
                      bins=20, size=0.25, color="white") +
  geom_contour(data=new_data.scale.2 %>% ungroup(), aes(x=local.rf.ll.scale, y=life_het,
                                                        z=log10(predicted)),
               bins=20, size=0.05, color="black") +
  # facet_wrap(~grad) +
  # scale_color_gradient2(low="blue", mid="white", high="red", midpoint = -0.24)+
  # scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = log10(mid))+
  scale_fill_manual(values=(colorRampPalette(c("steelblue1", "white", "tomato1"))(20)),
                    labels=sprintf("%.3f - %.3f", seq(0,3.3, by=0.165)[-21], seq(0,3.3, by=0.165)[-1]))+
  ylab("Genetic Load")+
  xlab("Scaled Site Suitability")+
  labs(fill="Relative Fitness")+
  theme_bw()+
  theme(text = element_text(size=15), panel.grid.major = element_blank(), aspect.ratio=1,
        panel.grid.minor = element_blank(), plot.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black"), legend.position = "right")+
  guides(fill = guide_legend(reverse = TRUE))
dev.off()

