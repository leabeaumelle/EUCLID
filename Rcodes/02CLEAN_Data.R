## Clean the data, calculate diversity, abundances, predation rates
## store the outputs (final datasets ready for analyses)


## I. Main Analysis ----------------------------------------------------

## Functions
library(dplyr)
library(vegan)

## 1. Natural enemies-------------------------------------------------------

## Load data-------------------------------------
NE <- read.csv("Output/NatEnemies_raw.csv")


## 1. Calculate abundance-------------------

## 2. Calculate rarefied richness-----------


## 3. Add explanatory variables and remove unnecessary columns----------------------
Treatment <- factor(levels=c("low plant div", "high plant div"))

Ldscp <- read.csv("Output/LanscapeVars.csv") # proportion of SNH in landscape

Community <- factor(levels = c("Vine", "Soil", "Vegetation")) # guild of NE 

Dist <- c(0,15,30) # distance in meters to the center of the flower strip

Site <- levels(as.factor(NE$couple))

Session <- levels(as.factor(NE$session)) # session nested within site?


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