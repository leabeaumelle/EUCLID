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
library(sjPlot)
library(lattice)
library(optimx)

# from previous codes
# library(lmerTest)
# library("blmeco")
# library("sjstats")

## Load data---------------------------------------------------------------------------
Abundance <- read.csv("Output/AbundanceClean.csv")
Abundance$Site <- as.factor(Abundance$Site)

Ldscp <- read.csv("Output/Landscapevars.csv")
Ldscp$Site <- as.factor(as.character(Ldscp$couple))
Ldscp$Ldscp <- Ldscp$HSN1000

Abundance <- left_join(Abundance, Ldscp, by = "Site")
Abundance$Site <- as.factor(Abundance$Site)

## Data exploration--------------------
summary(Abundance$Total)

hist(Abundance$Total)

par(mfrow=c(2,2))
plot(Abundance$Total ~ Abundance$Site, varwidth = TRUE)
plot(Abundance$Total ~ Abundance$Treatment, varwidth = TRUE)
plot(Abundance$Total ~ factor(Abundance$Distance), varwidth = TRUE)
plot(Abundance$Total ~ Abundance$Guild, varwidth = TRUE)
par(mfrow=c(1,1))

# Missing values
table(Abundance$Guild, Abundance$Distance)
# vegetation samplings at distances 15 and 30 m of flower strips: for ALL sites and treatments and sessions

table(Abundance$Site, Abundance$session, Abundance$Guild)
# vine samplings in session 1 for all treatments at site 10, and for high div treatment at site 9: for ALL distances
# vegetation sampling at distance 0m at site 10 for both treatments, and at site 9 for high div treatment 

# Extreme values
dotchart(Abundance$Total)

# Pairplots
pairs(Abundance[,c("Total", "Treatment", "Ldscp", "Guild", "Distance", "Site", "session")])

# Landscape relationships across treatments
plot(Abundance$Total ~ Abundance$Ldscp)
xyplot(Total ~ Ldscp |factor(Treatment)*factor(Guild)*factor(Distance), data = Abundance)

# Random effect structure : measurements from the same session at the same site are not independent
boxplot(Abundance$Total ~ factor(Abundance$Site), varwidth = TRUE)
plot(Abundance$Total ~ factor(Abundance$session), varwidth = TRUE)
xyplot(Total ~ factor(Site):factor(session), data = Abundance)


## Modelling-----------------------------------
# full model
mod1 <- glmer.nb(Total ~ Ldscp*Treatment*Guild + Treatment*Distance*Guild + 
                    (1|Site/session) ,
                  data=Abundance, control=glmerControl(optimizer="bobyqa"))

# model complexity vs. sample size
# k = 19: 16 (fixed) + 2 (random) + 1 (error term); n = 366
logLik(mod1) # gives the df of the model (19)
nrow(Abundance)
# the ratio n/k should be between 3 and 10 (Harrison, 2018)
366/19

# Rescale and center continuous predictors: landscape and distance variables
numcols <- grep("Ldscp|Dist",names(Abundance))
Abs <- Abundance
Abs[,numcols] <- scale(Abs[,numcols])
mod1_sc <- update(mod1,data=Abs)

## Check model assumptions------------------

## Inspect residuals
## residuals vs. fitted
plot(mod1_sc)

# check residuals with Dharma
res <- simulateResiduals(mod1_sc, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(scale(Abundance$Ldscp), res$scaledResiduals, asFactor = FALSE, main = "Landscape")
plotResiduals(Abundance$Guild, res$scaledResiduals, main = "Guild")
plotResiduals(Abundance$Treatment, res$scaledResiduals, main = "Treatment")
plotResiduals(Abundance$Distance, res$scaledResiduals, main = "Distance")
par(mfrow = c(1,1))

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


