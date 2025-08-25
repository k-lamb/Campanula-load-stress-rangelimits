library(reshape)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(lme4)
library(devtools)
library(car)
library(lmerTest)
library(blme)
library(dplyr)
library(ggeffects)
library(sjPlot)
library(optimx)
library(dplyr)
library(tidyr)

setwd("~/Desktop/Documents/Research/Paper Code/adaptation_load/load/")


#######################################################################################################
#
### adding germination to the heterosis data for lifetime heterosis calculations ###
#
#######################################################################################################

# read in data
germ.dat <- read.csv("data/germ_reshape_final.csv") # this data has been amended to this format from earlier melt & adds column for cross type, elerv group, and env type (lat or elev)
flw.dat <- read.csv("data/flw_reshape_final.csv") # this data has been amended to this format from earlier melt & adds column for cross type, elerv group, and env type (lat or elev)

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

### generate a table of everything ###
temp <- dat
t1 <- table(dat$Cross)
temp$uniqueID <- paste(temp$MatID, temp$PatID, sep="_")
temp <- temp %>% group_by(uniqueID) %>% slice_head(n=1)
t2 <- table(temp$Cross)
temp <- merge(t1 %>% as.data.frame(), t2 %>% as.data.frame(), by="Var1")



### lifetime heterosis ###
### APPALACHIAN ELEVATION ###

#ALBG
b1 <- mean(dat[dat$Cross == "ALBGxGA1",]$lifetime)
b2 <- mean(dat[dat$Cross == "TN34xALBG",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "ALBG",]$lifetime), mean(dat[dat$Cross == "GA1",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "ALBG",]$lifetime), mean(dat[dat$Cross == "TN34",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
ALBG_het <- mean(c(h1,h2))

#GA1
b1 <- mean(dat[dat$Cross == "ALBGxGA1",]$lifetime)
b2 <- mean(dat[dat$Cross == "GA1xTN34",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "ALBG",]$lifetime), mean(dat[dat$Cross == "GA1",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "GA1",]$lifetime), mean(dat[dat$Cross == "TN34",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
GA1_het <- mean(c(h1,h2))

#TN34
b1 <- mean(dat[dat$Cross == "TN34xALBG",]$lifetime)
b2 <- mean(dat[dat$Cross == "GA1xTN34",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "ALBG",]$lifetime), mean(dat[dat$Cross == "TN34",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "GA1",]$lifetime), mean(dat[dat$Cross == "TN34",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
TN34_het <- mean(c(h1,h2))

#GA2
b1 <- mean(dat[dat$Cross == "GA2xNC130",]$lifetime)
b2 <- mean(dat[dat$Cross == "NC16xGA2",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "GA2",]$lifetime), mean(dat[dat$Cross == "NC130",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "GA2",]$lifetime), mean(dat[dat$Cross == "NC16",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
GA2_het <- mean(c(h1,h2))

#NC130
b1 <- mean(dat[dat$Cross == "GA2xNC130",]$lifetime)
b2 <- mean(dat[dat$Cross == "NC130xNC16",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "GA2",]$lifetime), mean(dat[dat$Cross == "NC130",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "NC130",]$lifetime), mean(dat[dat$Cross == "NC16",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
NC130_het <- mean(c(h1,h2))

#NC16
b1 <- mean(dat[dat$Cross == "NC16xGA2",]$lifetime)
b2 <- mean(dat[dat$Cross == "NC130xNC16",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "GA2",]$lifetime), mean(dat[dat$Cross == "NC16",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "NC130",]$lifetime), mean(dat[dat$Cross == "NC16",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
NC16_het <- mean(c(h1,h2))

#EGG
b1 <- mean(dat[dat$Cross == "EGGxTN92",]$lifetime)
b2 <- mean(dat[dat$Cross == "PA104xEGG",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "EGG",]$lifetime), mean(dat[dat$Cross == "TN92",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "EGG",]$lifetime), mean(dat[dat$Cross == "PA104",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
EGG_het <- mean(c(h1,h2))

#TN92
b1 <- mean(dat[dat$Cross == "EGGxTN92",]$lifetime)
b2 <- mean(dat[dat$Cross == "TN92xPA104",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "EGG",]$lifetime), mean(dat[dat$Cross == "TN92",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "TN92",]$lifetime), mean(dat[dat$Cross == "PA104",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
TN92_het <- mean(c(h1,h2))

#PA104
b1 <- mean(dat[dat$Cross == "PA104xEGG",]$lifetime)
b2 <- mean(dat[dat$Cross == "TN92xPA104",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "EGG",]$lifetime), mean(dat[dat$Cross == "PA104",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "TN92",]$lifetime), mean(dat[dat$Cross == "PA104",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
PA104_het <- mean(c(h1,h2))

#NC91
b1 <- mean(dat[dat$Cross == "NC91xSK",]$lifetime)
b2 <- mean(dat[dat$Cross == "VA73xNC91",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "NC91",]$lifetime), mean(dat[dat$Cross == "SK",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "NC91",]$lifetime), mean(dat[dat$Cross == "VA73",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
NC91_het <- mean(c(h1,h2))

#SK
b1 <- mean(dat[dat$Cross == "NC91xSK",]$lifetime)
b2 <- mean(dat[dat$Cross == "SKxVA73",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "NC91",]$lifetime), mean(dat[dat$Cross == "SK",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "SK",]$lifetime), mean(dat[dat$Cross == "VA73",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
SK_het <- mean(c(h1,h2))

#VA73
b1 <- mean(dat[dat$Cross == "VA73xNC91",]$lifetime)
b2 <- mean(dat[dat$Cross == "SKxVA73",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "NC91",]$lifetime), mean(dat[dat$Cross == "VA73",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "SK",]$lifetime), mean(dat[dat$Cross == "VA73",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
VA73_het <- mean(c(h1,h2))


### WESTERN LATITUDE ###

#IJ
b1 <- mean(dat[dat$Cross == "IJxSH",]$lifetime)
b2 <- mean(dat[dat$Cross == "HPxIJ",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "IJ",]$lifetime), mean(dat[dat$Cross == "SH",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "IJ",]$lifetime), mean(dat[dat$Cross == "HP",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
IJ_het <- mean(c(h1,h2))

#SH
b1 <- mean(dat[dat$Cross == "IJxSH",]$lifetime)
b2 <- mean(dat[dat$Cross == "SHxHP",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "IJ",]$lifetime), mean(dat[dat$Cross == "SH",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "SH",]$lifetime), mean(dat[dat$Cross == "HP",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
SH_het <- mean(c(h1,h2))

#HP
b1 <- mean(dat[dat$Cross == "HPxIJ",]$lifetime)
b2 <- mean(dat[dat$Cross == "SHxHP",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "IJ",]$lifetime), mean(dat[dat$Cross == "HP",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "SH",]$lifetime), mean(dat[dat$Cross == "HP",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
HP_het <- mean(c(h1,h2))

#FC
b1 <- mean(dat[dat$Cross == "FCxQH",]$lifetime)
b2 <- mean(dat[dat$Cross == "HBxFC",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "FC",]$lifetime), mean(dat[dat$Cross == "QH",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "FC",]$lifetime), mean(dat[dat$Cross == "HB",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
FC_het <- mean(c(h1,h2))

#QH
b1 <- mean(dat[dat$Cross == "FCxQH",]$lifetime)
b2 <- mean(dat[dat$Cross == "QHxHB",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "FC",]$lifetime), mean(dat[dat$Cross == "QH",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "QH",]$lifetime), mean(dat[dat$Cross == "HB",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
QH_het <- mean(c(h1,h2))

#HB
b1 <- mean(dat[dat$Cross == "HBxFC",]$lifetime)
b2 <- mean(dat[dat$Cross == "QHxHB",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "FC",]$lifetime), mean(dat[dat$Cross == "HB",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "QH",]$lifetime), mean(dat[dat$Cross == "HB",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
HB_het <- mean(c(h1,h2))

#FT
b1 <- mean(dat[dat$Cross == "FTxMDP",]$lifetime)
b2 <- mean(dat[dat$Cross == "VxFT",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "FT",]$lifetime), mean(dat[dat$Cross == "MDP",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "FT",]$lifetime), mean(dat[dat$Cross == "V",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
FT_het <- mean(c(h1,h2))

#MDP
b1 <- mean(dat[dat$Cross == "FTxMDP",]$lifetime)
b2 <- mean(dat[dat$Cross == "MDPxV",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "FT",]$lifetime), mean(dat[dat$Cross == "MDP",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "MDP",]$lifetime), mean(dat[dat$Cross == "V",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
MDP_het <- mean(c(h1,h2))

#V
b1 <- mean(dat[dat$Cross == "VxFT",]$lifetime)
b2 <- mean(dat[dat$Cross == "MDPxV",]$lifetime)
m1 <- mean(c(mean(dat[dat$Cross == "FT",]$lifetime), mean(dat[dat$Cross == "V",]$lifetime)))
m2 <- mean(c(mean(dat[dat$Cross == "MDP",]$lifetime), mean(dat[dat$Cross == "V",]$lifetime)))
h1 <- (b1 - m1)/m1
h2 <- (b2 - m2)/m2
V_het <- mean(c(h1,h2))

### DATA FRAME OF HETEROSIS ###

#set up df
pop_hets <- c(ALBG_het, EGG_het, FC_het, FT_het, GA1_het, GA2_het, HB_het, HP_het, IJ_het, MDP_het, NC130_het, NC16_het, NC91_het,
              PA104_het, QH_het, SH_het, SK_het, TN34_het, TN92_het, V_het, VA73_het) #21
pop_list <- c("ALBG", "EGG", "FC", "FT", "GA1", "GA2", "HB", "HP", "IJ", "MDP", "NC130", "NC16", "NC91", "PA104", "QH", "SH", "SK",
              "TN34", "TN92", "V", "VA73")
het.df <- as.data.frame(pop_hets, pop_list)
het.df <- cbind(pop = rownames(het.df), het.df)
rownames(het.df) <- 1:nrow(het.df)

#add specific columns we need for analysis: lineage, env var and env group

#Create lineage
het.df$lineage <- NA
het.df <- het.df %>%
  mutate(lineage = case_when(
    pop == "ALBG" | pop == "FC" | pop == "FT" | pop == "GA1" | pop == "GA2" | 
      pop == "HB" | pop == "HP" | pop == "IJ" | pop == "MDP" | pop == "NC130" | 
      pop == "NC16" | pop == "QH" | pop == "SH" | pop == "TN34" | pop == "V" 
    ~ "Western",
    pop == "EGG" | pop == "NC91" | pop == "PA104" | pop == "SK" | pop == "TN92" | 
      pop == "VA73"
    ~ "Appalachian",
    pop == "SM" | pop == "VA112"
    ~ "Eastern"))

#env var
het.df$env.var <- NA
het.df <- het.df %>%
  mutate(env.var = case_when(
    pop == "ALBG" | pop == "EGG" | pop == "GA1" | pop == "GA2" | pop == "NC130" | 
      pop == "NC16" | pop == "NC91" | pop == "PA104" | pop == "SK" | pop == "SM" | 
      pop == "TN34" | pop == "TN92" | pop == "VA112" | pop == "VA73"
    ~ "elev",
    pop == "FC" | pop == "FT" | pop == "GA1" | pop == "HB" | pop == "HP" | 
      pop == "IJ" | pop == "MDP" | pop == "QH" | pop == "SH" | pop == "V"
    ~ "lat")) #adds environmental delineater column

#env group
#Create env.group
het.df$env.group <- NA
het.df <- het.df %>%
  mutate(env.group = case_when(
    pop == "ALBG" | pop == "EGG" | pop == "GA1" | pop == "HP" | pop == "IJ" | 
      pop == "PA104" | pop == "SH" | pop == "SM" | pop == "TN34" | pop == "TN92" | 
      pop == "VA112"
    ~ "low",
    pop == "FC" | pop == "HB" | pop == "QH" 
    ~ "mid",
    pop == "FT" | pop == "GA2" | pop == "MDP" | pop == "NC130" | pop == "NC16" | 
      pop == "NC91" | pop == "SK" | pop == "V" | pop == "VA73" 
    ~ "high"
  )) #adds 3 tiers of environmental group


write.csv(het.df, "data/heterosis/heterosis_lifetime.csv")


### combine all heterosis data into single file
het.df1 <- read.csv("data/heterosis/heterosis_lifetime.csv")
het.df2 <- read.csv("data/heterosis/heterosis_flw.csv")
het.df3 <- read.csv("data/heterosis/heterosis_germ.csv")

het.df <- cbind(het.df1[,c(2,4,5,6,3)], het.df2[,3], het.df3[,3])
colnames(het.df)[5:7] <- c("life_het", "flw_het", "germ_het")

coords <- read.csv("data/pop_coords.csv")
het.df <- merge(het.df, coords, by="pop")

# hard coding transect info
het.df <- het.df %>% 
  mutate(transect = case_when(
    pop == "EGG" | pop == "VA73" ~ 1,
    pop == "PA104" | pop == "SK" ~ 2,
    pop == "TN92" | pop == "NC91" ~ 3,
    pop == "GA1" | pop == "GA2" ~ 4,
    pop == "ALBG" | pop == "NC16" ~ 5,
    pop == "TN34" | pop == "NC130" ~ 6,
    pop == "HP" | pop == "FC" | pop == "FT" ~ 7,
    pop == "SH" | pop == "HB" | pop == "V" ~ 8,
    pop == "IJ" | pop == "QH" | pop == "MDP" ~ 9
  ))

write.csv(het.df, "data/heterosis/heterosis_all.csv")

