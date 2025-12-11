
# modeling
library(lme4)
library(glmmTMB)
library(car)
library(emmeans)
library(effects)
library(pscl)
library(DHARMa)
library(ggResidpanel)
library(blme)
library(glmmTMB)

# data manipulation
library(dplyr)
# library(tidyverse)

# plotting
library(ggplot2)
library(sjPlot)
library(ggeffects)
library(ggpubr)

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/local_adaptation/")

# read in data and adjust variables to factors (and their order)
cg <- read.csv("data/elevation/elev_rel_abs_HERB.csv", stringsAsFactors = T)[,-c(1:3)]
cg.no_herb <- read.csv("data/elevation/elev_rel_abs_NO_HERB.csv", stringsAsFactors = T)[,-1]

cg <- cg %>% mutate_at(vars(block), as.factor) %>% mutate(pop.transect = if_else(pop == "SM", "e4", pop.transect))
cg.no_herb <- cg.no_herb %>% mutate_at(vars(block), as.factor) %>% mutate(pop.transect = if_else(pop == "SM", "e4", pop.transect))

cg$garden.group <- factor(cg$garden.group, c("low", "high"))
cg$pop.group <- factor(cg$pop.group, c("low", "high"))

cg.no_herb$garden.group <- factor(cg.no_herb$garden.group, c("low", "high"))
cg.no_herb$pop.group <- factor(cg.no_herb$pop.group, c("low", "high"))

cg.no_herb <- cg.no_herb %>% 
  group_by(garden.transect) %>% mutate(fl.fr.scale.gr=scale(fl.fr)) %>% 
  ungroup()

### RELATIVE MODELS--NO HERBIVORY

# relative flowering & fruiting

# # base LMER... does not account for zero-bounding
# cg.lme.tot <- blmer(rel.fr ~ garden.group*pop.group + (1-garden.group|full.garden/block) + (1-pop.group|pop),
#                     control=lmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=2e5)),
#                     data=cg.no_herb)
cg.lme.tot <- glmmTMB(rel.fr ~ garden.group*pop.group + (1-garden.group|full.garden/block) + (1-pop.group|pop),
                      family = tweedie(link = "log"), data=cg.no_herb)


# check assumptions
# simulationOutput <- simulateResiduals(fittedModel = cg.lme.tot, plot = T)
plot(residuals(cg.lme.tot)~cg.no_herb$pop.group)
plot(residuals(cg.lme.tot)~cg.no_herb$garden.group)

Anova(cg.lme.tot, type=3)

# plotting standard errors and estimated means
# cg.se <- effects::effect("garden.group*pop.group", cg.lme.tot) %>% as.data.frame()
cg.se <- emmeans(cg.lme.tot, ~ garden.group * pop.group, type="response") %>% as.data.frame()

cg.se <- cg.se %>% 
  mutate(garden.group = case_when(garden.group == "low" ~ "Core",
                                  garden.group == "high" ~ "Edge")) %>% 
  mutate(pop.group = case_when(pop.group == "low" ~ "Core",
                               pop.group == "high" ~ "Edge"))
cg.se$pop.group <- factor(cg.se$pop.group, levels = c("Core","Edge"))
cg.se$garden.group <- factor(cg.se$garden.group, levels = c("Core", "Edge"))

# b <-
s.size=6
dodge=0.5
ggplot()+
  geom_line(data=cg.se, aes(x=garden.group, y=response, group=pop.group, color=pop.group),
            position = position_dodge(width=dodge), linetype=3)+
  geom_errorbar(data=cg.se, aes(x=garden.group, ymin=(response-SE), ymax=(response+SE), color=pop.group), 
                width=0.25, position = position_dodge(width=dodge))+
  
  geom_point(data=cg.se, aes(x=garden.group, y=response, fill=pop.group), 
             color="black", size=s.size, shape=24, position = position_dodge(width=dodge))+
  
  theme_bw()+
  ylim(c(0,2.5))+
  scale_color_manual(values=c("Core" = "firebrick1", "Mid"="orchid2", "Edge" = "skyblue3"))+
  scale_fill_manual(values=c("Core" = "firebrick1", "Mid"="orchid2", "Edge" = "skyblue3"))+
  ylab("Relative fitness")+
  xlab("Garden range position")+
  geom_hline(yintercept=1, linetype=2)+
  theme(text = element_text(size = 28), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio = 2)
ggsave("plots/elevation/rel_fr_LA_color_fix.jpg")

# ggpubr::ggarrange(a,b, ncol=2)
# ggsave("plots/elevation/elev_fl_fr_LA.jpg")


