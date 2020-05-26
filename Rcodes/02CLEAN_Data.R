## Clean the data, calculate diversity, abundances, predation rates
## store the outputs (final datasets ready for analyses)


## I. Main Analysis ----------------------------------------------------

## Functions
library(dplyr)
library(vegan)
library(tidyr)

## 1. Natural enemies-------------------------------------------------------

## Load data-------------------------------------
NE <- read.csv("Output/NatEnemies_raw.csv")

## Rename variables------------------------------
NE$Treatment <- factor(ifelse(NE$bande=="BE", "Low Div", "High Div"))
NE$Guild <- factor(ifelse(NE$piege == "B", "Vine",
                          ifelse(NE$piege == "F", "Vegetation",
                                 "Soil")))
NE$Distance <- ifelse(NE$mod == "00m", 0,
                      ifelse(NE$mod  == "05m", 15,
                             30))
NE$Site <- factor(NE$couple)

## Add landscape--------------------------------
Ldscp <- read.csv("Output/LandscapeVars.csv") # proportion of SNH in landscape
NE_lddf <- left_join(NE, Ldscp, by = "couple")

# landscape variable is proportion of SNH in 1000 m radius
NE$Ldscp <- NE_lddf$HSN1000

## CHECK Data for missing values, NAs and zero abundances---------------------------------
NE %>% group_by(Site,Treatment,Guild,Distance, session) %>% summarize(n()) # no. obs
NE %>% group_by(Site,Treatment,Guild,Distance, session) %>% summarize(sum(eff)) # sum abdc has NAs

## Nas for eff: site 9 and 10 only
filter(NE, is.na(eff)) %>% group_by(Site, Treatment, Guild, Distance, session) %>% summarize(n())

# Based on planning of sampling, should have NA for Vegetation in high div and low div treatments in site 9, 
# and in site 10, should be NA for all Vine and Vegetation samplings 
filter(NE, Site=="9") %>% group_by(Treatment, Guild, Distance, session) %>% summarize(mean(eff), min(eff), max(eff)) %>% View()

with(filter(NE, Site=="9"), table(Treatment, Guild, Distance, session))
# Missing rows: Vegetation at 15 and 30m for all sessions and all treatments
# Check in initial table (raw data) if values are present: no

## NAs were handled inconsistently in NE dataset: some are NAs and some are missing
## Remove all NAs from dataset
NE <- na.omit(NE)

## Missing values are not random: table shows that for Vegetation sampling 
# none of the sites have data for Distance 15 and 30 m
# vegetation was absent due to drought (except in strips) and sampling was impossible
NE %>% group_by(Guild, Distance) %>% summarize(n())

NE %>% group_by(Site, Treatment, Guild, Distance, session) %>% summarize(n())

## 1. Calculate abundance-------------------

## Total abundance per guild, distance and session----------------------------------------
Abundance <- NE %>% group_by(Site, Treatment, Guild, Distance, session) %>% 
  summarize(Total = sum(eff)) %>% ungroup()

## Missing values--------------------------
# we have 366 observations
nrow(Abundance)
# complete combination would be 648 observations
nrow(expand_grid(Site = c(1:12), 
                 Treatment = c("high", "low"), 
                 Guild=c("Soil", "Vegetation", "Vine"), 
                 Distance = c(0,15,30), 
                 session=c(1:3)))

## Save dataset
write.csv(Abundance, "Output/AbundanceClean.csv")


##### STOPPED HERE : TO DO : diversity calculation-------------------------------
## 2. Calculate rarefied richness-----------
Diversity <- NE %>% group_by(Site, Treatment, Guild, Distance, session) %>% 
  summarize()

## 3. Add explanatory variables and remove unnecessary columns----------------------



## Save cleaned data----------------------
write.csv(NE, "Output/NatEnemiesDiv_clean.csv")

# 2. Predation of Lobesia eggs-------------------------------------------------------

## Load data--------------------------------------------
Pred <- read.csv("Output/PredationPest_raw.csv")

## 1. Calculate predation rates---------------------------------

## 2. Add explanatory variables---------------------------------



# Save cleaned data---------------------------------
write.csv("Output/PredationPest_clean.csv")


## Supplementary Analysis - Year effect