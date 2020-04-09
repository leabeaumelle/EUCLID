## Create figures for the manuscript

## The script creates and saves three main figures

## Functions----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)

## Data----------------------------------------------------------------------------
NEData <- read.csv("Output/NatEnemies_clean.csv")
Preddata <- read.csv("Output/PredationPest_clean.csv")



## Figure 1 - overall effect of flower strips on natural enemies and predation
# three panels with overall flower strip effect across landscapes and communitiy types


## Figure 2 - response of different guilds at different distances of flower strips
# three to nine panels showing how response depends on community type and distance to strip


## Figure 3 - landscape effect
# figure showing for one selected outcome how landscape modulates flower strip effects
