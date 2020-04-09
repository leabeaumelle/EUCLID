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

# from previous codes
# library(lmerTest)
# library("blmeco")
# library("sjstats")

## Load data---------------------------------------------------------------------------
Preddata <- read.csv("Output/PredationPest_clean.csv")

## Run model---------------------------------------------------------------------------
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

