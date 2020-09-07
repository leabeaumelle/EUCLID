## Models for natural enemies richness

# The scripts tests the effects of flower strips, landscape on Nat enemies genus richness
# Model selection, model assumption checks, and results stored into output folder.

## Main Analysis - Landscape*PlantDiv*Guild

# Y ~ Landscp*Treatment*Guild + Treatment*Distc*Guild + (1|Site:Session)

# Landscp = cont. var: gradient of proportion of SNH
# Treatment = factor: flower strip (high div) vs. grassy strip (low div)
# Guild = factor: vine, soil, vegetation
# Distc = continuous: 0m, 15m, 30 m

# site = pairs of plots with same Landscp
# Session = multiple sessions per plot


## Functions ------------------------------
library(dplyr)
library(lme4)
library(MuMIn)
library(DHARMa)
library(MASS)
library(sjPlot)
library(lattice)
library(optimx)
library(car)

## Load data------------------------------

Diversity <- read.csv("Output/DiversityClean.csv")
Diversity$Site <- as.factor(Diversity$Site)

Ldscp <- read.csv("Output/Landscapevars.csv")
Ldscp$Site <- as.factor(as.character(Ldscp$couple))
Ldscp$Ldscp <- Ldscp$HSN1000

Diversity <- left_join(Diversity, Ldscp, by = "Site")
Diversity$Site <- as.factor(Diversity$Site)
Diversity$session <- as.factor(as.character(Diversity$session))

# Rescale and center continuous predictors: landscape and distance variables
numcols <- grep("Ldscp|Dist",names(Diversity))
Div <- Diversity
Div[,numcols] <- scale(Div[,numcols])

# Remove vegetation guild data
Div <- droplevels(Div[Div$Guild != "Vegetation",])


## 1. Data exploration-----------------------

summary(Diversity$GenusR)
summary(Diversity$TaxaR)

par(mfrow=c(1,2))
hist(Diversity$GenusR, main = "Genus richness")
plot(Diversity$GenusR, Diversity$TaxaR); abline(0,1)

par(mfrow=c(2,2), mar = c(4,4,2,2))
plot(Diversity$TaxaR ~ Diversity$Site, varwidth = TRUE, main = "Taxonomic")
plot(Diversity$TaxaR ~ Diversity$Treatment, varwidth = TRUE)
plot(Diversity$TaxaR ~ factor(Diversity$Distance), varwidth = TRUE)
plot(Diversity$TaxaR ~ Diversity$Guild, varwidth = TRUE)
par(mfrow=c(1,1))

# Missing values
table(Diversity$Guild, Diversity$Distance)
# vegetation samplings at distances 15 and 30 m of flower strips: for ALL sites and treatments and sessions

table(Diversity$Site, Diversity$session, Diversity$Guild)
# vine samplings in session 1 for all treatments at site 10, and for high div treatment at site 9: for ALL distances
# vegetation sampling at distance 0m at site 10 for both treatments, and at site 9 for high div treatment 

# Extreme values
dotchart(Diversity$TaxaR)

# Pairplots
pairs(Diversity[,c("TaxaR","GenusR", "Treatment", "Ldscp", "Guild", "Distance", "Site", "session")])

# Landscape relationships across treatments
plot(Diversity$TaxaR ~ Diversity$Ldscp)

xyplot(GenusR ~ Ldscp |factor(Treatment)*factor(Guild)*factor(Distance), data = Diversity)

# Random effect structure : measurements from the same session at the same site are not independent
par(mfrow = c(1,2))
boxplot(Diversity$TaxaR ~ factor(Diversity$Site), varwidth = TRUE, main = "Taxonomic richness")
plot(Diversity$TaxaR ~ factor(Diversity$session), varwidth = TRUE)

xyplot(TaxaR ~ factor(Site):factor(session), data = Diversity)

# Taxonomic richness, without vegetation guild
par(mfrow=c(2,2), mar = c(4,4,2,2))
plot(Div$TaxaR ~ Div$Site, varwidth = TRUE, main = "Taxonomic")
plot(Div$TaxaR ~ Div$Treatment, varwidth = TRUE)
plot(Div$TaxaR ~ factor(Div$Distance), varwidth = TRUE)
plot(Div$TaxaR ~ Div$Guild, varwidth = TRUE)
par(mfrow=c(1,1))


plot(Div$TaxaR~Div$Ldscp)
plot(Div$TaxaR~Div$session)
plot(Div$TaxaR~Div$Site)

## 2. Full models ---------------------------

# taxonomic richness is assuming that within each sample,
# individuals not identified to the genus belong to a genus not represented by individuals identified
# i.e. is a number of taxonomic units mixing different taxonomic resolution

## Full model1: taxonomic richness ----------
mod1_taxar <- lmer(TaxaR ~ Ldscp*Treatment*Guild + Distance*Treatment*Guild + 
                           (1|Site:session),
                         REML = TRUE,
                         data = Div,
                         control= lmerControl(optimizer="bobyqa"))

## residuals vs. fitted
plot(mod1_taxar)

# check residuals with Dharma
res <- simulateResiduals(mod1_taxar, plot = T)

# residuals vs. predictors
par(mfrow = c(2,2), mar = c(4, 4, 2, 2))
plotResiduals(res, Div$Ldscp, asFactor = FALSE, main = "Landscape")
plotResiduals(res, Div$Guild, main = "Guild")
plotResiduals(res, Div$Treatment, main = "Treatment")
plotResiduals(res, Div$Distance,  main = "Distance")
par(mfrow = c(1,1))

# residuals are uniform, but violate variance heterogeneity, and strong trend with landscape suggest non linear effect needed

## Full model2: Add non linear landscape effect------
mod1_taxar2 <- lmer(TaxaR ~ poly(Ldscp, 2)*Treatment*Guild + Distance*Treatment*Guild + 
                            (1|Site:session),
                          REML = TRUE,
                          data = Div,
                          control= lmerControl(optimizer="bobyqa"))

## residuals vs. fitted
plot(mod1_taxar2)

# check residuals with Dharma
res2 <- simulateResiduals(mod1_taxar2, plot = T)

testQuantiles(res2)

# deviation detected, but is not significantly deviating from uniform

# residuals vs. predictors
par(mfrow = c(2,2), mar = c(4, 4, 2, 2))
plotResiduals(res2, Div$Ldscp, asFactor = FALSE, main = "Landscape")
plotResiduals(res2, Div$Guild, main = "Guild")
plotResiduals(res2, Div$Treatment, main = "Treatment")
plotResiduals(res2, Div$Distance,  main = "Distance")
par(mfrow = c(1,1))

# variance heterogeneity is still present,

## Full model3:  Log transform variable -----
mod1_taxar3 <- lmer(log10(TaxaR+1) ~ Ldscp*Treatment*Guild + Distance*Treatment*Guild + 
                            (1|Site:session),
                          REML = TRUE,
                          data = Div,
                          control= lmerControl(optimizer="bobyqa"))

## residuals vs. fitted
plot(mod1_taxar3)

# check residuals with Dharma
res3 <- simulateResiduals(mod1_taxar3, plot = T)

# no deviation but several outliers

# residuals vs. predictors
par(mfrow = c(2,2), mar = c(4, 4, 2, 2))
plotResiduals(res3, Div$Ldscp, asFactor = FALSE, main = "Landscape")
plotResiduals(res3, Div$Guild, main = "Guild")
plotResiduals(res3, Div$Treatment, main = "Treatment")
plotResiduals(res3, Div$Distance,  main = "Distance")
par(mfrow = c(1,1))

# residuals look better, but landscape effect detected

## Full model4:  Log transform variable And non linear landscape effect---------------
modfin <- lmer(log10(TaxaR+1) ~ poly(Ldscp, 2)*Treatment*Guild + Distance*Treatment*Guild + 
                      (1|Site:session),
                    REML = TRUE,
                    data = Div,
                    control= lmerControl(optimizer="bobyqa"))

## residuals vs. fitted
plot(modfin)

# check residuals with Dharma
res4 <- simulateResiduals(modfin, plot = T)

# deviation from uniform detected

testQuantiles(res4)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res4, Div$Ldscp, asFactor = FALSE, main = "Landscape")
plotResiduals(res4, Div$Guild, main = "Guild")
plotResiduals(res4, Div$Treatment, main = "Treatment")
plotResiduals(res4, Div$Distance,  main = "Distance")
par(mfrow = c(1,1))

## Full model: new distribution---------

# trying out new distribution to fit the data: negative binomial chosen because
# 1. data are counts (number of species), integer values, positive
# 2. the variance is much higher than the mean taxa richness (overdispersion)
mean(Div$TaxaR); var(Div$TaxaR)

# Full model 1b: Negative Binomial --------------------------------------
mod1b_taxar <- glmer.nb(TaxaR ~ Ldscp*Treatment*Guild + Distance*Treatment*Guild + 
                     (1|Site:session),
                   data = Div,
                   control= glmerControl(optimizer="bobyqa"))

## residuals vs. fitted
plot(mod1b_taxar)

# check residuals with Dharma
resb <- simulateResiduals(mod1b_taxar, plot = T)

# residuals vs. predictors
par(mfrow = c(2,2), mar = c(4, 4, 2, 2))
plotResiduals(resb, Div$Ldscp, asFactor = FALSE, main = "Landscape")
plotResiduals(resb, Div$Guild, main = "Guild")
plotResiduals(resb, Div$Treatment, main = "Treatment")
plotResiduals(resb, Div$Distance,  main = "Distance")
par(mfrow = c(1,1))

# residuals are uniform, and have homogeneous variance, but there are patterns with landscape

# Full model 2b: NB with non linear landscape effect --------------------------------------
mod2b_taxar <- glmer.nb(TaxaR ~ poly(Ldscp,2)*Treatment*Guild + Distance*Treatment*Guild + 
                          (1|Site:session),
                        data = Div,
                        control= glmerControl(optimizer="bobyqa"))

## residuals vs. fitted
plot(mod2b_taxar)

# check residuals with Dharma
res2b <- simulateResiduals(mod2b_taxar, plot = T)

# residuals vs. predictors
par(mfrow = c(2,2), mar = c(4, 4, 2, 2))
plotResiduals(res2b, Div$Ldscp, asFactor = FALSE, main = "Landscape")
plotResiduals(res2b, Div$Guild, main = "Guild")
plotResiduals(res2b, Div$Treatment, main = "Treatment")
plotResiduals(res2b, Div$Distance,  main = "Distance")
par(mfrow = c(1,1))

# residuals are uniform, and have homogeneous variance, no pattern with landscape and other covariates

## SAVE Final Full model ---------------
tab_model(mod2b_taxar)

saveRDS(mod2b_taxar, file = "Output/NEDiversity_FullModel.rds")


## 3. Model selection taxo richness data--------

# Estimate with ML instead of REML: not possible with glmer.nb
# mod1_taxarFull <- update(modfin2, REML = FALSE)

# first step
drop1(mod2b_taxar, test = "Chisq")

# The interaction Treatment:Guild:Distance is ns
modsel1b_taxar <- update(mod2b_taxar, .~. -Treatment:Guild:Distance)
drop1(modsel1b_taxar, test = "Chisq")

# interaction landscape:treatment:guild ns (dropping three way interactions first)
modsel2b_taxar <- update(modsel1b_taxar, .~. -poly(Ldscp, 2):Treatment:Guild)
drop1(modsel2b_taxar, test= "Chisq")

# Treatment:distance is the least significant interaction
modsel3b_taxar <- update(modsel2b_taxar, .~. -Treatment:Distance)
drop1(modsel3b_taxar, test = "Chisq")

# Guild:distance is the least significant interaction
modsel4b_taxar <- update(modsel3b_taxar, .~. -Guild:Distance)
drop1(modsel4b_taxar, test = "Chisq")

# Ldscp:treatment is the least significant interaction
modsel5b_taxar <- update(modsel4b_taxar, .~. -poly(Ldscp, 2):Treatment)
drop1(modsel5b_taxar, test = "Chisq")

# Ldscp:guild is the least significant interaction
modsel6b_taxar <- update(modsel5b_taxar, .~. -poly(Ldscp, 2):Guild )
drop1(modsel6b_taxar, test = "Chisq")

# Distance is the least signficant
modsel7b_taxar <- update(modsel6b_taxar, .~. -Distance)
drop1(modsel7b_taxar, test = "Chisq")

# Treatment:guild is the least significant
modsel8b_taxar <- update(modsel7b_taxar, .~. -Treatment:Guild)
drop1(modsel8b_taxar, test = "Chisq")

# Treatment is least significant
modsel9b_taxar <- update(modsel8b_taxar, .~. -Treatment)
drop1(modsel9b_taxar, test = "Chisq")

# Landscape is least significant
modsel10b_taxar <- update(modsel9b_taxar, .~. -poly(Ldscp, 2))
drop1(modsel10b_taxar, test = "Chisq")

# no further deletion

## Validation of optimal model--------

modopt_taxar <- modsel10b_taxar

# Check model assumptions for the optimal model
## Inspect residuals
plot(modopt_taxar)

# check residuals with Dharma
res <- simulateResiduals(modopt_taxar, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# distribution not great, but K-S test is not statistically significant

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res, scale(Div$Ldscp), asFactor = FALSE, main = "Landscape")
plotResiduals(res, Div$Guild, main = "Guild")
plotResiduals(res, Div$Treatment, main = "Treatment")
plotResiduals(res, Div$Distance,  main = "Distance")
par(mfrow = c(1,1))

# trend with landscape and with guild...

# Save model results-------
tab_model(modopt_taxar)

saveRDS(modopt_taxar, file = "Output/NEDiversity_OptimalModel.rds")


# Plot showing the Guild effect------
plot_model(modopt_taxar, type = "pred", terms = c("Guild"))

# Anova tables---------
Anova(modopt_taxar)
Anova(mod2b_taxar)


# Save results of LRTs-------



