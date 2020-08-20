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
library(dfoptim)
library(car)

## Load data---------------------------------------------------------------------------
Abundance <- read.csv("Output/AbundanceClean.csv")
Abundance$Site <- as.factor(Abundance$Site)

Ldscp <- read.csv("Output/Landscapevars.csv")
Ldscp$Site <- as.factor(as.character(Ldscp$couple))
Ldscp$Ldscp <- Ldscp$HSN1000

Abundance <- left_join(Abundance, Ldscp, by = "Site")
Abundance$Site <- as.factor(Abundance$Site)
Abundance$session <- as.factor(as.character(Abundance$session))

# Rescale and center continuous predictors: landscape and distance variables
numcols <- grep("Ldscp|Dist",names(Abundance))
Abs <- Abundance
Abs[,numcols] <- scale(Abs[,numcols])

# Remove vegetation guild data to fit models
Abs <- droplevels(Abs[Abs$Guild != "Vegetation",])


## Data exploration----------------------------------------------------------
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
# vegetation samplings at distances 15 and 30 m: for ALL sites and treatments and sessions

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


## Is a random effect on the Site necessary? (https://stackoverflow.com/questions/53034261/warning-lme4-model-failed-to-converge-with-maxgrad)
# each site has a single landscape value for all observations
all(rowSums(with(Abs, table(Site, Ldscp))>0)==1)

## Full models -------------------------------------------------------------------------
# using dataset Abs where continuous variable are scaled (unscaled vars results in convergence issues)

mod1_sc <- glmer.nb(Total ~ Ldscp*Treatment*Guild + Treatment*Distance*Guild + 
                    (1|Site:session) ,
                  data=Abs, control=glmerControl(optimizer="bobyqa"))

# model complexity vs. sample size
# k = 14: 12 (fixed) + 1 (random) + 1 (error term); n = 315
logLik(mod1_sc) # gives the df of the model (14)
nrow(Abs)
# the minimum ratio n/k should be between 3 and 10 (Harrison, 2018)
315/14


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
plotResiduals(res, Abs$Ldscp, asFactor = FALSE, main = "Landscape")
plotResiduals(res, Abs$Guild, main = "Guild")
plotResiduals(res, Abs$Treatment, main = "Treatment")
plotResiduals(res, Abs$Distance,  main = "Distance")
par(mfrow = c(1,1))

## Patterns with landscape: need to add a non-linear effect of the landscape


## Full model with Non linear trend with landscape----

# optimizer bobyqa has convergence issues: so I use Nelder_Mead
mod1_nl <- glmer.nb(Total ~ poly(Ldscp, 2)*Treatment*Guild + Treatment*Distance*Guild + 
                      (1|Site:session) ,
                    data = Abs, control = glmerControl(optimizer = "bobyqa"))


# minimum n/k should be between 3-10: ok
nrow(Abs)/attr(logLik(mod1_nl),"df")

# LRT test of quadratic term: highly significant
anova(mod1_nl, mod1_sc)

# # Check model assumptions----

## residuals vs. fitted
plot(mod1_nl)

# check residuals with Dharma
res <- simulateResiduals(mod1_nl, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res$scaledResiduals, Abs$Ldscp,  asFactor = FALSE, main = "Landscape")
plotResiduals(res$scaledResiduals, Abs$Guild,  main = "Guild")
plotResiduals(res$scaledResiduals, Abs$Treatment,  main = "Treatment")
plotResiduals(res$scaledResiduals, Abs$Distance,  main = "Distance")
par(mfrow = c(1,1))

## Non linear effect of the landscape ameliorate model residuals and is significant according
# to ANOVA: new full model


## Results of Final Full model-------------------
tab_model(mod1_nl)

# anova type II (wald chi-square tests : least recommended solution according to Bolker (https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html))
Anova(mod1_nl)

saveRDS(mod1_nl, file = "Output/NEAbundance_FullModel.rds")


## Model selection------------------------
# step 1
drop1(mod1_nl, test = "Chisq")

#both three way interactions are ns : dropping treatment:guild:distance first
modsel1 <- update(mod1_nl, .~. -Treatment:Guild:Distance)

# setp 2
drop1(modsel1, test = "Chisq")

# Three way interaction Ld:Treatment:Guild ns
modsel2 <- update(modsel1, .~. -poly(Ldscp, 2):Treatment:Guild)

# step 3
drop1(modsel2, test = "Chisq")

# least significant term is the Treatment:Distance interaction (lowest LRT and highest P)
modsel3 <- update(modsel2, .~. -Treatment:Distance)

# step 4
drop1(modsel3, test = "Chisq")

# least significant term is the Guild:Distance interaction
modsel4 <- update(modsel3, .~. -Guild:Distance)

# step5
drop1(modsel4, test = "Chisq")

# Treatment:Guild ns
modsel5 <- update(modsel4, .~. -Treatment:Guild)

# step 6
drop1(modsel5, test = "Chisq")

# Ldscp:Guild is ns but p = 0.06...
modsel6 <- update(modsel5, .~. -poly(Ldscp, 2):Guild)

# no further deletion
drop1(modsel6, test = "Chisq")


modOpt <- modsel6

## Checking otimal model assumptions
## Inspect residuals
## residuals vs. fitted
plot(modOpt)

# check residuals with Dharma
res <- simulateResiduals(modOpt, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res$scaledResiduals, (Abs$Ldscp),  asFactor = FALSE, main = "Landscape")
plotResiduals(res$scaledResiduals, Abs$Guild, main = "Guild")
plotResiduals(res$scaledResiduals, Abs$Treatment, main = "Treatment")
plotResiduals(res$scaledResiduals, Abs$Distance,  main = "Distance")
par(mfrow = c(1,1))


## Results
tab_model(modOpt)
Anova(modOpt)


# different resposne to landscape depending on guild
plot_model(modOpt, type = "pred", terms = c("Treatment", "Guild"))
plot_model(modOpt, type = "pred", terms = c("Ldscp [all]", "Treatment"))
plot_model(modOpt, type = "pred", terms = c("Ldscp [all]", "Guild"))
plot_model(modOpt, type = "pred", terms = c("Distance [all]", "Guild"))

saveRDS(mod1_nl, file = "Output/NEAbundance_OptModel.rds")

