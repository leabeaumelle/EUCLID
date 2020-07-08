## Models for predation rates

# The script tests the effects of flower strips, landscape on Lobesia eggs predation
# Model selection, model assumption checks, and results stored into output folder.

## Main Analysis - Landscape*FlowerStrips + NatEnemies

# Y ~ Landscp*Treatment + Treatment*Distc + NatEnemiesDensity + (1|Site/Session)

# Landscp = cont. var: gradient of proportion of SNH
# Treatment = factor: flower strip vs. grassy strip
# Distc = continuous variable taking values 0, 15 and 30 m


# site = pairs of plots with same Landscp
# Session = multiple sessions per plot

## Functions---------------------------------------------------------------------------
library(dplyr)
library(lme4)
library(MuMIn)
library(DHARMa)
library(MASS)
library(sjPlot)
library(lattice)
library(optimx)
library(car)

## Load data---------------------------------------------------------------------------
Pred <- read.csv("Output/PredationPest_clean.csv")
Pred$Site <- as.factor(as.character(Pred$Site))

Ldscp <- read.csv("Output/Landscapevars.csv")
Ldscp$Site <- as.factor(as.character(Ldscp$couple))
Ldscp$Ldscp <- Ldscp$HSN1000

Pred <- left_join(Pred, Ldscp, by = "Site")
Pred$Site <- as.factor(as.character(Pred$Site))

## Data exploration---------------------------------------------------------------------------

summary(Pred)

# NAs predation rates: sites 9, 6 and 11, in high div and low div treatments, at distance 0m and in sessions 1 and 2
Pred[is.na(Pred$PredRate),]

# bimodal distribution: lots of zeros and ones
hist(Pred$PredRate)

# relationships
plot(Pred$PredRate ~ Pred$Ldscp)
plot(Pred$PredRate ~ Pred$Treatment)
plot(Pred$PredRate ~ Pred$Distance)

# random effects
plot(Pred$PredRate ~ Pred$Site)
plot(Pred$PredRate ~ Pred$Session)


## Run model---------------------------------------------------------------------------


## Check assumptions ---------------------------------------------------------------------


## Store model results as output ---------------------------------------------------------

