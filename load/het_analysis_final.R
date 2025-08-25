
# modeling
library(lme4)
library(car)
library(blme)

# data manipulation
library(reshape)
library(dplyr)
library(tidyr)

# plotting
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(sjPlot)
library(devtools)

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/load/")


het.df <- read.csv("data/heterosis/heterosis_all.csv", stringsAsFactors = T)[,-1]
het.df$env.group <- factor(het.df$env.group, levels=c("low", "mid", "high")) # fixes order issue in graphing later
het.df <- het.df %>% mutate_at(vars(transect), as.factor)

# filtering out populations with negative estimates of heterosis
het.boot <- read.csv("data/heterosis/heterosis_lifetime_bootstrap_samplewithreplacement.csv", stringsAsFactors = T)[,-1]
het.boot <- het.boot %>% 
  group_by(pop, lineage, env.var, env.group) %>% 
  reframe(mean = mean(pop_hets), median=mean(pop_hets),
          se = sd(pop_hets)/sqrt(n())) %>% 
  mutate(ci95 = mean+(3*se))
pop.filter <- (het.boot %>% filter(ci95 < 0))$pop %>% droplevels() %>% unique()

`%ni%` <- Negate(`%in%`)
het.df <- het.df %>% filter(pop %ni% pop.filter)

het.lat.lmer <- lm(life_het ~ env.group, data=het.df %>% filter(env.var=="lat"))
car::Anova(het.lat.lmer, type=3)

het.elev.lmer <- lm(life_het ~ env.group, data=het.df %>% filter(env.var=="elev"))
car::Anova(het.elev.lmer, type=3)

# all together
het.all.lm2 <- lm(life_het ~ env.var * env.group, data=het.df %>% filter(env.group != "mid"))
het.all.lm1 <- lm(life_het ~ env.var + env.group, data=het.df %>% filter(env.group != "mid")) # best model
het.all.lm <- lm(life_het ~ env.group, data=het.df %>% filter(env.group != "mid"))

car::Anova(het.all.lm, type=3)
car::Anova(het.all.lm1, type=3)
car::Anova(het.all.lm2, type=3)

BIC(het.all.lm2) - BIC(het.all.lm)
BIC(het.all.lm2) - BIC(het.all.lm1)
BIC(het.all.lm1) - BIC(het.all.lm)

sjPlot::plot_model(het.all.lm1, "pred", terms = c("env.group"))+theme_bw()

## LMM version
# het.all.lmer <- lmer(life_het ~ env.var + (1-env.var|env.group), data=het.df %>% filter(env.group != "mid"))
# het.all.lmer <- lmer(life_het ~ env.group + (1|env.var), data=het.df %>% filter(env.group != "mid"))
# car::Anova(het.all.lmer, type=3)
# sjPlot::plot_model(het.all.lmer, type="pred", pred.type="fe")+theme_bw()




### plotting

temp <- het.df
temp$env.var.lin <- paste0(temp$env.var, "_", temp$lineage)

temp %>%
  group_by(env.var) %>%
  dplyr::summarise(life_het.mu = mean(life_het, na.rm=T),
                   life_het.se.low = life_het.mu - (sd(life_het, na.rm=T)/sqrt(n())),
                   life_het.se.high = life_het.mu + (sd(life_het, na.rm=T)/sqrt(n())))


temp <- temp %>%
  group_by(env.var, env.group) %>%
  dplyr::summarise(life_het.mu = mean(life_het, na.rm=T),
                   life_het.se.low = life_het.mu - (sd(life_het, na.rm=T)/sqrt(n())),
                   life_het.se.high = life_het.mu + (sd(life_het, na.rm=T)/sqrt(n())))

temp <- temp %>% 
  mutate(range.pos = case_when(
    env.group == "low" ~ "Core",
    env.group == "mid" ~ "Mid",
    env.group == "high" ~ "Edge"
  ))


s.size=6
dodge=0.5

temp2 <- rbind(temp, data.frame(env.var="elev", env.group="mid", life_het.mu=NA, life_het.se.low=NA, life_het.se.high=NA, range.pos="mid"))
temp3 <- temp %>% filter(env.var=="elev")
temp$range.pos <- factor(temp$range.pos, c("Core", "Mid", "Edge"))
temp2$range.pos <- factor(temp2$range.pos, c("Core", "Mid", "Edge"))
temp3$range.pos <- factor(temp3$range.pos, c("Core", "Mid", "Edge"))

temp2 <- temp2 %>% mutate(env.var = if_else(env.var=="elev", "Steep", "Shallow"))

ggplot(temp2 %>% na.omit(), aes(x=range.pos, y=life_het.mu, group=env.var))+
  geom_line(data=temp2 %>% na.omit(), aes(x=range.pos, y=life_het.mu, group=env.var), 
            position = position_dodge(width=dodge), linetype=3)+
  # geom_line(data=temp3, aes(x=range.pos, y=life_het.mu, group=env.var),
  #           position = position_dodge(width=dodge), linetype=3)+
  geom_errorbar(data=temp2 %>% na.omit(), aes(x=range.pos, ymin=life_het.se.low, ymax=life_het.se.high, y=life_het.mu, color=range.pos), 
                width=0.25, position = position_dodge(width=dodge))+
  geom_point(aes(shape=env.var, fill=range.pos), 
             size=s.size, color="black", position = position_dodge(width=dodge))+
  ylim(c(0, 0.5))+
  ylab("Genetic load")+
  xlab("Population range position")+
  # geom_hline(yintercept=0,linetype=2)+
  scale_color_manual(values=c("Core" = "firebrick1", "Mid"="orchid2", "Edge" = "skyblue3"))+
  scale_fill_manual(values=c("Core" = "firebrick1", "Mid"="orchid2", "Edge" = "skyblue3"))+
  scale_shape_manual(values = c("Steep" = 24, "Shallow" = 21))+
  theme_bw()+
  theme(text = element_text(size = 28), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_rect(fill="white"), axis.line = element_line(colour = "black"),
        aspect.ratio = 1.5)+
  labs(color="Population range position", fill="Population range position", shape="Environmental gradient")
ggsave("plots/het_bootstrap_color_legend.jpg")

