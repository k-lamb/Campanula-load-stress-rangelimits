
library(dplyr)
library(car)
library(betareg)

# data frame that matches crosses performed per population
# dat <- read.csv("~/Desktop/Documents/Research/Q2/data/status1_crosses.csv")
# 
# dat <- dat %>%
#   group_by(Family) %>% 
#   dplyr::summarise(count=sum(status))
# 
# write.csv(dat, "~/Desktop/Documents/Research/Q2/data/cross_counts_by_pop.csv")

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/load/data/")
flw.dat <- read.csv("flw_reshape_final.csv", stringsAsFactors = T) # this data has been amended to this format from earlier melt & adds column for cross type, elerv group, and env type (lat or elev)
germ.dat <- read.csv("germ_reshape_final.csv", stringsAsFactors = T) # this data has been amended to this format from earlier melt & adds column for cross type, elerv group, and env type (lat or elev)

flw.test <- flw.dat %>% 
  group_by(FULL_ID, Cohort, Rep, Cross.Type, Lineage, Env.Var, Env.Group) %>% 
  dplyr::reframe(survival = survival,
                   bolt=bolt,
                   flw=sum(flw.count, na.rm=T)) %>%
  group_by(FULL_ID, Cohort, Rep, Cross.Type, Lineage, Env.Var, Env.Group) %>% 
  slice(1)

germ.test <- germ.dat %>% 
  group_by(FULL_ID, Cohort, Rep, Cross.Type, Lineage, Env.Var, Env.Group) %>% 
  dplyr::reframe(germ=max(Germ.Count/Seed.Total, na.rm=T)) %>%
  group_by(FULL_ID, Cohort, Rep, Cross.Type, Lineage, Env.Var, Env.Group) %>% 
  slice(1)

germ.test <- germ.test %>% 
  mutate(germ.prop = if_else(germ == 1, 0.9999, if_else(germ == 0, 0.0001, germ)))

# components separately
m1 <- lm(flw ~ Cohort*Cross.Type*Env.Var, data=flw.test)
car::Anova(m1, type=3)

m2 <- betareg::betareg(germ.prop ~ Cohort*Cross.Type*Env.Var, data=germ.test, link = "logit")
car::Anova(m2, type=3)




### lifetime all together model

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/load/")


#######################################################################################################
#
### adding germination to the heterosis data for lifetime heterosis calculations ###
#
#######################################################################################################

# read in data
germ.dat <- read.csv("data/germ_reshape_final.csv") #this data has been amended to this format from earlier melt & adds column for cross type, elerv group, and env type (lat or elev)
flw.dat <- read.csv("data/flw_reshape_final.csv") #this data has been amended to this format from earlier melt & adds column for cross type, elerv group, and env type (lat or elev)

### data prep germination
Factors<-c("Cross","MatID","PatID","Cohort","Cross.Type","Env.Var","Env.Group", "Rep", "Lineage")
germ.dat[Factors]<-lapply(germ.dat[Factors],factor) # converting vars to factors
germ.dat$Env.Group.fac <- factor(germ.dat$Env.Group, levels=c("low", "mid", "high")) # fixes order issue in graphing later
germ.dat$Germ.Prop <- germ.dat$Germ.Count/germ.dat$Seed.Total
germ.dat.sub <- subset(germ.dat, DSP == 25 | DSP == 26)
germ.dat.sub <- na.omit(germ.dat.sub)
germ.dat.sub <- germ.dat.sub[,c("Cross", "MatID", "PatID", "Rep", "Cohort", "FULL_ID", "Germ.Prop")]

### data prep flowering
Factors<-c("Cross","MatID","PatID","Cohort","Cross.Type","Env.Var","Env.Group", "Rep", "Lineage", "FULL_ID")
flw.dat[Factors]<-lapply(flw.dat[Factors],factor) # converting vars to factors
flw.dat$Env.Group.fac <- factor(flw.dat$Env.Group, levels=c("low", "mid", "high")) # fixes order issue in graphing later
flw.dat.sub <- aggregate(flw.dat, flw.count ~ Cross + MatID + PatID + Rep + Cohort + bolt + survival + Lineage + Cross.Type + Env.Var + Env.Group + Env.Group.fac + FULL_ID, FUN = "sum")

### merge flws and germ
dat <- merge(germ.dat.sub, flw.dat.sub, by=c("Cross", "MatID", "PatID", "Cohort", "Rep", "FULL_ID"), all=T)
dat <- dat %>% drop_na(Germ.Prop)

# add 0's to any individuals which didn't survive germination (NA's right now)
for (i in 1:dim(dat)[1]) {
  if ((dat$Germ.Prop[i] == 0) == TRUE) {
    dat$bolt[i] <- 0
    dat$survival[i] <- 0
    dat$flw.count[i] <- 0
  }
}

### rerun data prep from earlier to assign missing lineage, etc.
# Create Lineage
dat$Lineage <- NA
dat <- dat %>%
  mutate(Lineage = case_when(
    Cross == "ALBG" | Cross == "FC" | Cross == "FT" | Cross == "GA1" | Cross == "GA2" | 
      Cross == "HB" | Cross == "HP" | Cross == "IJ" | Cross == "MDP" | Cross == "NC130" | 
      Cross == "NC16" | Cross == "QH" | Cross == "SH" | Cross == "TN34" | Cross == "V" |
      Cross == "ALBGxGA1" | Cross == "GA1xTN34" | Cross == "GA2xNC130" | Cross == "NC130xNC16" | 
      Cross == "NC16xGA2" | Cross == "TN34xALBG" | Cross == "FCxQH" | Cross == "FTxMDP" | 
      Cross == "HBxFC" | Cross == "HPxIJ" | Cross == "IJxSH" | Cross == "MDPxV" | 
      Cross == "QHxHB" | Cross == "SHxHP" | Cross == "VxFT"
    ~ "Western",
    Cross == "EGG" | Cross == "NC91" | Cross == "PA104" | Cross == "SK" | Cross == "TN92" | 
      Cross == "VA73" |
      Cross == "EGGxTN92" | Cross == "NC91xSK" | Cross == "PA104xEGG" | Cross == "SKxVA73"|
      Cross == "TN92xPA104" | Cross == "VA73xNC91"
    ~ "Appalachian",
    Cross == "SM" | Cross == "VA112" |
      Cross == "SMxVA112" 
    ~ "Eastern"))

# Create Cross.Type
dat$Cross.Type <- NA
dat <- dat %>%
  mutate(Cross.Type = case_when(
    Cross == "ALBG" | Cross == "EGG" | Cross == "FC" | Cross == "FT" | Cross == "GA1" | 
      Cross == "GA2" | Cross == "HB" | Cross == "HP" | Cross == "IJ" | Cross == "MDP" | 
      Cross == "NC130" | Cross == "NC16" | Cross == "NC91" | Cross == "PA104" | 
      Cross == "QH" | Cross == "SH" | Cross == "SK" | Cross == "SM" | Cross == "TN34" | 
      Cross == "TN92" | Cross == "V" | Cross == "VA112" | Cross == "VA73" ~ "WI")) # adds WI column

dat[["Cross.Type"]][is.na(dat[["Cross.Type"]])] <- "BW" # add BW to column

# Create Env.Var
dat$Env.Var <- NA
dat <- dat %>%
  mutate(Env.Var = case_when(
    Cross == "ALBG" | Cross == "EGG" | Cross == "GA1" | Cross == "GA2" | Cross == "NC130" | 
      Cross == "NC16" | Cross == "NC91" | Cross == "PA104" | Cross == "SK" | Cross == "SM" | 
      Cross == "TN34" | Cross == "TN92" | Cross == "VA112" | Cross == "VA73" |
      Cross == "ALBGxGA1" | Cross == "EGGxTN92" | Cross == "GA1xTN34" | Cross == "GA2xNC130" |
      Cross == "NC130xNC16" | Cross == "NC16xGA2" | Cross == "NC91xSK" | Cross == "PA104xEGG" | 
      Cross == "SKxVA73" | Cross == "TN34xALBG" | Cross == "TN92xPA104" | Cross == "VA73xNC91" |
      Cross == "SMxVA112"
    ~ "elev",
    Cross == "FC" | Cross == "FT" | Cross == "GA1" | Cross == "HB" | Cross == "HP" | 
      Cross == "IJ" | Cross == "MDP" | Cross == "QH" | Cross == "SH" | Cross == "V" |
      Cross == "FCxQH" | Cross == "FTxMDP" | Cross == "HBxFC" | Cross == "HPxIJ" | 
      Cross == "IJxSH" | Cross == "MDPxV" | Cross == "QHxHB" | Cross == "SHxHP" | 
      Cross == "VxFT"
    ~ "lat")) # adds environmental delineater column

# Create Env.Group
dat$Env.Group <- NA
dat <- dat %>%
  mutate(Env.Group = case_when(
    Cross == "ALBG" | Cross == "EGG" | Cross == "GA1" | Cross == "HP" | Cross == "IJ" | 
      Cross == "PA104" | Cross == "SH" | Cross == "SM" | Cross == "TN34" | Cross == "TN92" | 
      Cross == "VA112" |
      Cross == "ALBGxGA1" | Cross == "EGGxTN92" | Cross == "GA1xTN34" | Cross == "PA104xEGG" | 
      Cross == "TN34xALBG" | Cross == "TN92xPA104" | Cross == "IJxSH" | Cross == "SHxHP" |
      Cross == "HPxIJ" | Cross == "SMxVA112"
    ~ "low",
    Cross == "FC" | Cross == "HB" | Cross == "QH" |
      Cross == "FCxQH" | Cross == "QHxHB" | Cross == "HBxFC"
    ~ "mid",
    Cross == "FT" | Cross == "GA2" | Cross == "MDP" | Cross == "NC130" | Cross == "NC16" | 
      Cross == "NC91" | Cross == "SK" | Cross == "V" | Cross == "VA73" |
      Cross == "MDPxV" | Cross == "VxFT" | Cross == "FTxMDP" | Cross == "GA2xNC130" |
      Cross == "NC130xNC16" | Cross == "NC16xGA2" | Cross == "NC91xSK" |
      Cross == "SKxVA73" | Cross == "VA73xNC91"
    ~ "high"
  )) # adds 3 tiers of environmental group

dat <- na.omit(dat)

dat$lifetime <- dat$Germ.Prop*dat$bolt*dat$survival*dat$flw.count # generate lifetime fitness metric
m3 <- lm(lifetime ~ Cohort*Cross.Type*Env.Var, data=dat)
car::Anova(m3, type=3)

dat$combo <- paste0(dat$MatID,"_",dat$PatID)


lme <- lmer(lifetime ~ Cohort*Cross.Type*Env.Var + (1-Cross.Type|Cross/combo), data=dat)
Anova(lme, type=3)
summary(lme)

glmm <- glmmTMB(lifetime ~ Cohort*Cross.Type*Env.Var + (1-Cross.Type|Cross/combo), data=dat,
                family = tweedie(link = "log"))
Anova(glmm, type=3)
summary(glmm)


lme <- lmer(lifetime ~ Cohort + (1-Cross|combo), data=dat)
Anova(lme, type=3)
summary(lme)

plot_model(lme, "pred")+theme_bw()
