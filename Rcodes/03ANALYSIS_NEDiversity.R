## Models for natural enemies richness

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
Diversity <- read.csv("Output/DiversityClean.csv")
Diversity$Site <- as.factor(Diversity$Site)

Ldscp <- read.csv("Output/Landscapevars.csv")
Ldscp$Site <- as.factor(as.character(Ldscp$couple))
Ldscp$Ldscp <- Ldscp$HSN1000

Diversity <- left_join(Diversity, Ldscp, by = "Site")
Diversity$Site <- as.factor(Diversity$Site)
Diversity$session <- as.factor(as.character(Diversity$session))

## Data exploration--------------------
summary(Diversity$GenusR)
summary(Diversity$TaxaR)

par(mfrow=c(1,2))
hist(Diversity$GenusR, main = "Genus richness")
hist(Diversity$TaxaR, main = "Taxa richness")


par(mfrow=c(4,2))
plot(Diversity$GenusR ~ Diversity$Site, varwidth = TRUE, main = "Genus")
plot(Diversity$TaxaR ~ Diversity$Site, varwidth = TRUE, main = "Taxonomic")
plot(Diversity$GenusR ~ Diversity$Treatment, varwidth = TRUE)
plot(Diversity$TaxaR ~ Diversity$Treatment, varwidth = TRUE)
plot(Diversity$GenusR ~ factor(Diversity$Distance), varwidth = TRUE)
plot(Diversity$TaxaR ~ factor(Diversity$Distance), varwidth = TRUE)
plot(Diversity$GenusR ~ Diversity$Guild, varwidth = TRUE)
plot(Diversity$TaxaR ~ Diversity$Guild, varwidth = TRUE)
par(mfrow=c(1,1))

# Missing values
table(Diversity$Guild, Diversity$Distance)
# vegetation samplings at distances 15 and 30 m of flower strips: for ALL sites and treatments and sessions

table(Diversity$Site, Diversity$session, Diversity$Guild)
# vine samplings in session 1 for all treatments at site 10, and for high div treatment at site 9: for ALL distances
# vegetation sampling at distance 0m at site 10 for both treatments, and at site 9 for high div treatment 

# Extreme values
dotchart(Diversity$GenusR)
dotchart(Diversity$TaxaR)

# Pairplots
pairs(Diversity[,c("TaxaR","GenusR", "Treatment", "Ldscp", "Guild", "Distance", "Site", "session")])

# Landscape relationships across treatments
par(mfrow=c(1,2))
plot(Diversity$GenusR ~ Diversity$Ldscp)
plot(Diversity$TaxaR ~ Diversity$Ldscp)

xyplot(GenusR ~ Ldscp |factor(Treatment)*factor(Guild)*factor(Distance), data = Diversity)

# Random effect structure : measurements from the same session at the same site are not independent
boxplot(Diversity$TaxaR ~ factor(Diversity$Site), varwidth = TRUE, main = "Taxonomic richness")
boxplot(Diversity$GenusR ~ factor(Diversity$Site), varwidth = TRUE, main = "Genus richness")

plot(Diversity$TaxaR ~ factor(Diversity$session), varwidth = TRUE)
plot(Diversity$GenusR ~ factor(Diversity$session), varwidth = TRUE)

xyplot(TaxaR ~ factor(Site):factor(session), data = Diversity)
xyplot(GenusR ~ factor(Site):factor(session), data = Diversity)


## Modelling-----------------------------------
# full model
mod1 <- lmer(GenusR ~ Ldscp*Treatment*Guild + Treatment*Distance*Guild + 
                   (1|Site/session) ,
                 data=Diversity, control=lmerControl(optimizer="bobyqa"))

# model complexity vs. sample size
# k = 19: 16 (fixed) + 2 (random) + 1 (error term); n = 366
logLik(mod1) # gives the df of the model (19)
nrow(Diversity[!is.na(Diversity$GenusR),])
# the ratio n/k should be between 3 and 10 (Harrison, 2018)
364/19

# Rescale and center continuous predictors: landscape and distance variables
numcols <- grep("Ldscp|Dist",names(Diversity))
Div <- Diversity
Div[,numcols] <- scale(Div[,numcols])
mod1_sc <- update(mod1,data=Div)

## Check model assumptions------------------

## Inspect residuals
## residuals vs. fitted
plot(mod1_sc)

# check residuals with Dharma
res <- simulateResiduals(mod1_sc, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors : no strong sign of var heterogeneity
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
plotResiduals(scale(Diversity$Ldscp[!is.na(Diversity$GenusR)]), res$scaledResiduals, asFactor = FALSE, main = "Landscape")
plotResiduals(Diversity$Guild[!is.na(Diversity$GenusR)], res$scaledResiduals, main = "Guild")
plotResiduals(Diversity$Treatment[!is.na(Diversity$GenusR)], res$scaledResiduals, main = "Treatment")
plotResiduals(Diversity$Distance[!is.na(Diversity$GenusR)], res$scaledResiduals, main = "Distance")
par(op)

# residuals vs. random factors : plots suggest variance heterogeneity for the different sites
op2 <- par(mfrow=c(1,2), mar = c(4,4,2,2))
plotResiduals(Diversity$Site[!is.na(Diversity$GenusR)], res$scaledResiduals)
plotResiduals(Diversity$session[!is.na(Diversity$GenusR)], res$scaledResiduals)
par(op2)

# plot suggesting non-normality (blue line not so close to the red line)
plot_model(mod1_sc, type = "slope")


## Full model 2: add non-linear pattern with landscape
mod2 <- lmer(GenusR ~ poly(Ldscp, 2)*Treatment*Guild + Treatment*Distance*Guild + 
               (1|Site/session) ,
             data=Diversity, control=lmerControl(optimizer="bobyqa"))

# model complexity vs. sample size
# k = 25; n = 364
logLik(mod2) # gives the df of the model (25)
nrow(Diversity[!is.na(Diversity$GenusR),])
# the ratio n/k should be between 3 and 10 (Harrison, 2018)
364/25

# Rescale and center continuous predictors: landscape and distance variables
mod2_sc <- update(mod2,data=Div)

## Check model assumptions------------------

## Inspect residuals
## residuals vs. fitted
plot(mod2_sc)

# check residuals with Dharma
res2 <- simulateResiduals(mod2_sc, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(scale(Diversity$Ldscp[!is.na(Diversity$GenusR)]), res2$scaledResiduals, asFactor = FALSE, main = "Landscape")
plotResiduals(Diversity$Guild[!is.na(Diversity$GenusR)], res2$scaledResiduals, main = "Guild")
plotResiduals(Diversity$Treatment[!is.na(Diversity$GenusR)], res2$scaledResiduals, main = "Treatment")
plotResiduals(Diversity$Distance[!is.na(Diversity$GenusR)], res2$scaledResiduals, main = "Distance")
par(mfrow = c(1,1))

# residuals vs. random factors : plots suggest variance heterogeneity for the different sites
op3 <- par(mfrow=c(1,3), mar = c(4,4,2,2))
plotResiduals(Diversity$Site[!is.na(Diversity$GenusR)], res2$scaledResiduals)
plotResiduals(Diversity$session[!is.na(Diversity$GenusR)], res2$scaledResiduals)
plotResiduals(Diversity$session[!is.na(Diversity$GenusR)]:Diversity$Site[!is.na(Diversity$GenusR)], res2$scaledResiduals)
par(op3)

## Full model 3 : variance heterogeneity

