## Models for natural enemies abundance

# The scripts tests the effects of flower strips, landscape on Nat enemies densities
# Model selection, model assumption checks, and results stored into output folder.

## Main Analysis - Landscape*FlowerStrips*Community

# Y ~ Landscp*Treatment*Community + Treatment*Distc*Community + (1|Site/Session)

# Landscp = cont. var: gradient of proportion of SNH
# Treatment = factor: flower strip vs. grassy strip
# Community = factor: vine, soil, vegetation
# Distc = factor: 0m, 5m, 20 m or continuous

# site = pairs of plots with same Landscp
# Session = multiple sessions per plot

## Functions---------------------------------------------------------------------------
library(dplyr)
library(lme4)
library(MuMIn)
library(DHARMa)
library(MASS)

# from previous codes
# library(lmerTest)
# library("blmeco")
# library("sjstats")

## Load data---------------------------------------------------------------------------
Abundance <- read.csv("Output/AbundanceClean.csv")
Abundance$Site <- as.factor(Abundance$Site)

Ldscp <- read.csv("Output/Landscapevars.csv")
Ldscp$Site <- as.factor(as.character(Ldscp$couple))

Abundance <- left_join(Abundance, Ldscp, by = "Site")
Abundance$Site <- as.factor(Abundance$Site)

## Descriptive stats and plot--------------------
summary(Abundance$Total)

hist(Abundance$Total)

par(mfrow=c(2,2))
plot(Abundance$Total ~ Abundance$Site)
plot(Abundance$Total ~ Abundance$Treatment)
plot(Abundance$Total ~ Abundance$Distance)
plot(Abundance$Total ~ Abundance$Guild)
par(mfrow=c(1,1))

# Missing values


## Modelling-----------------------------------
mod1 <- glmer.nb(Total ~ Ldscp*Treatment*Guild + Treatment*Distance*Guild + 
                    (1|Site/session) ,
                  data=Abundance) #, control=glmerControl(optimizer="bobyqa"))


## Run model---------------------------------------------------------------------------
# From Arthur
# NEData <- read.csv("Output/NatEnemies_clean.csv")

totEN <- glmer.nb(Y ~ Landscp*Treatment*Community + Treatment*Distc*Community + 
                    (1|Site/Session) ,
                  data=NEData , control=glmerControl(optimizer="bobyqa"))

standardize(totEN)
summary(totEN)
overdisp_fun(totEN)
plot(simulateResiduals(totEN))
set_totEN<-dredge(totEN,rank="AICc")
topmod_set_totEN<-get.models(set_totEN,subset=delta<2) 
topmod_set_totEN # Deux meilleurs modeles : bande + mod + (1|couple) & mod + (1|couple)
modavg_totEN<-model.avg(topmod_set_totEN) ; modavg_totEN
summary(modavg_totEN) ; confint(modavg_totEN)
r.squaredGLMM(topmod_set_totEN$'4')
# Significativement moins d'Ennemis naturels a 20m.

## Check assumptions ---------------------------------------------------------------------


## Store model results as output ---------------------------------------------------------


