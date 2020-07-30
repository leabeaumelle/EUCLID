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
# library(MuMIn)
library(DHARMa)
library(sjPlot)
# library(lattice)
# library(optimx)
library(car)

## Load data---------------------------------------------------------------------------
Pred <- read.csv("Output/PredationPest_clean.csv")
Pred$Site <- as.factor(as.character(Pred$Site))

Ldscp <- read.csv("Output/Landscapevars.csv")
Ldscp$Site <- as.factor(as.character(Ldscp$couple))
Ldscp$Ldscp <- Ldscp$HSN1000

Pred <- left_join(Pred, Ldscp, by = "Site")
Pred$Site <- as.factor(as.character(Pred$Site))
Pred$Session <- as.factor(as.character(Pred$Session))

# Rescale and center continuous predictors: landscape and distance variables
numcols <- grep("Ldscp|Dist",names(Pred))
Pred_sc <- Pred
Pred_sc[,numcols] <- scale(Pred[,numcols])

## Data exploration---------------------------------------------------------------------------

summary(Pred)

# NAs predation rates: sites 9, 6 and 11, in high div and low div treatments, at distance 0m and in sessions 1 and 2
Pred[is.na(Pred$PredRate),]

# dataset structure
summary(Pred$Site) # Site 6 has only 90 obs versus 114 in the others
table(Pred$Site, Pred$Session) # SIte 6 has twice less obsrvations than other sites in session1
table(Pred$Site, Pred$Session, Pred$Distance) # these missing observations are spread across distances,
table(Pred$Site, Pred$Session, Pred$Distance, Pred$Treatment) # these missing observations are spread across distances,

# bimodal distribution: lots of zeros and ones
hist(Pred$PredRate)

# relationships
plot(Pred$PredRate ~ Pred$Ldscp)
plot(Pred$PredRate ~ Pred$Treatment)
plot(Pred$PredRate ~ Pred$Distance)

# random effects
plot(Pred$PredRate ~ Pred$Site)
plot(Pred$PredRate ~ Pred$Session)

# data structure
Pred %>% group_by(Site, Treatment, Distance, Session) %>% summarize(n()) %>% View()

## Full model---------------------------------------------------------------------------

modFull1 <- glmer(PredRate ~ Ldscp*Treatment + Treatment*Distance + (1|Site/Session),
                  family = binomial, data = Pred_sc,
                  weights = rep(10, nrow(Pred_sc)))  # there are 10 eggs per sentinel card

# model complexity vs. sample size
logLik(modFull1)  # gives the df of the model (8)
nrow(Pred[!is.na(Pred$PredRate),])/attr(logLik(modFull1), "df")  # the ratio n/k largely exceeds 10


# check residuals with Dharma: not good
res1 <- simulateResiduals(modFull1, plot = T)

# residuals vs. predictors
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
plotResiduals(res1, Pred_sc$Ldscp[!is.na(Pred_sc$PredRate)], asFactor = FALSE, main = "Landscape")
plotResiduals(res1, Pred_sc$Treatment[!is.na(Pred_sc$PredRate)], main = "Treatment")
plotResiduals(res1, Pred_sc$Distance[!is.na(Pred_sc$PredRate)],  main = "Distance")
par(op)


# residuals: not great: deviation from uniform distrib and overdispersion


## full model 2: try poisson and model the number of predated eggs (PredRate*10)
modFull2 <- glmer(PredRate*10 ~ Ldscp*Treatment + Treatment*Distance + (1|Site/Session),
                  family = poisson, data = Pred_sc) 

# check residuals with Dharma: not great
res2 <- simulateResiduals(modFull2, plot = T)

# residuals vs. predictors
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
plotResiduals(res2, Pred_sc$Ldscp[!is.na(Pred_sc$PredRate)], asFactor = FALSE, main = "Landscape")
plotResiduals(res2, Pred_sc$Treatment[!is.na(Pred_sc$PredRate)], main = "Treatment")
plotResiduals(res2, Pred_sc$Distance[!is.na(Pred_sc$PredRate)],  main = "Distance")
par(op)

# residuals still not great


## full model 3: try negative binomial
modFull3 <- glmer.nb(PredRate*10 ~ Ldscp*Treatment + Treatment*Distance + (1|Site/Session),
                     data = Pred_sc) 

# check residuals with Dharma: better but outlier and quantile deviation detected
res3 <- simulateResiduals(modFull3, plot = T)

# residuals vs. predictors : no strong sign of var heterogeneity
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
plotResiduals(res3, Pred_sc$Ldscp[!is.na(Pred_sc$PredRate)], asFactor = FALSE, main = "Landscape")
plotResiduals(res3, Pred_sc$Treatment[!is.na(Pred_sc$PredRate)], main = "Treatment")
plotResiduals(res3, Pred_sc$Distance[!is.na(Pred_sc$PredRate)],  main = "Distance")
par(op)


## full model 4: try non linear landscape effects
modFull4 <- glmer.nb(PredRate*10 ~ poly(Ldscp, 2)*Treatment + Treatment*Distance + (1|Site/Session),
                     data = Pred_sc) 

# check residuals with Dharma: much better. One outlier
res4 <- simulateResiduals(modFull4, plot = T)

# residuals vs. predictors : no strong sign of var heterogeneity
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
plotResiduals(res4, Pred_sc$Ldscp[!is.na(Pred_sc$PredRate)], asFactor = FALSE, main = "Landscape")
plotResiduals(res4, Pred_sc$Treatment[!is.na(Pred_sc$PredRate)], main = "Treatment")
plotResiduals(res4, Pred_sc$Distance[!is.na(Pred_sc$PredRate)],  main = "Distance")
par(op)


## Model with negative binomial and polynomial 2 effect landscape is the best so far
#BUT model is singular: parameters are on the boundary of feasible parameter space
# variances of linear combination of effects are zero
isSingular(modFull4)
isSingular(modFull4, tol = 1e-7) # no longer singular at tolerance 1e-7

summary(modFull4)
# shows that variance of random effect of the Site is zero:
# could mean that landscape effect is taking all the variance, none is left for site random effect?



## full model 5: try non linear landscape effects
modFull5 <- glmer(PredRate ~ poly(Ldscp, 2)*Treatment + Treatment*Distance + (1|Site/Session),
                     family = binomial, data = Pred_sc,
                    weights = rep(10, nrow(Pred_sc)))  # there are 10 eggs per sentinel card


# check residuals with Dharma: not great.
res5 <- simulateResiduals(modFull5, plot = T)


## Adding a non linear effect to the binomial distrib is not ameliorating residuals


## Save Full model results-------
modfull <- modFull4

# tab_model(modfull)
Anova(modfull)

saveRDS(modfull, file = "Output/PredRate_FullModel.rds")

## MODEL SELECTION--------------

# step 1
drop1(modfull, test = "Chisq")

# Treatment:Distance is not significant
modsel1 <- update(modfull,  .~. -Treatment:Distance)

# step 2
drop1(modsel1, test = "Chisq")

# Distance ns
modsel2 <- update(modsel1, .~. -Distance)

# step 3
drop1(modsel2, test = "Chisq")

# no further deletion, but singular warning is still on
modOpt <- modsel2
isSingular(modOpt)

summary(modOpt) # issue is the random effect structure: variance is zero
VarCorr(modOpt)

# check residuals with Dharma: outlier test significant (2 values > 1)
resOpt <- simulateResiduals(modOpt, plot = T)

# residuals vs. predictors : no strong sign of var heterogeneity nor patterns with predictors
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
plotResiduals(resOpt, Pred_sc$Ldscp[!is.na(Pred_sc$PredRate)], asFactor = FALSE, main = "Landscape")
plotResiduals(resOpt, Pred_sc$Treatment[!is.na(Pred_sc$PredRate)], main = "Treatment")
plotResiduals(resOpt, Pred_sc$Distance[!is.na(Pred_sc$PredRate)],  main = "Distance")
par(op)

## The model is singular because the SIte random effect has VAR = 0

## Find the optimal random effect structure ----------

# Following Zuur et al. Book protocol (p.121)
## Step 1
modre1 <- glm.nb(PredRate * 10 ~ poly(Ldscp, 2) * Treatment + Treatment * Distance, 
                 data = Pred_sc)

modre2 <- glmer.nb(PredRate * 10 ~ poly(Ldscp, 2) * Treatment + Treatment * Distance + (1 | Site),
                   data = Pred_sc)
modre3 <- glmer.nb(PredRate * 10 ~ poly(Ldscp, 2) * Treatment + Treatment * Distance + (1 | Site/Session),
                   data = Pred_sc)

# Step 2: Best model is with Site/session as random effects
AIC(modre1, modre2, modre3)
anova(modre3, modre2, modre1)

## And starting with Session random effect

## Step 1
modre1b <- glm.nb(PredRate * 10 ~ poly(Ldscp, 2) * Treatment + Treatment * Distance, 
                 data = Pred_sc)

modre2b <- glmer.nb(PredRate * 10 ~ poly(Ldscp, 2) * Treatment + Treatment * Distance + (1 | Session),
                   data = Pred_sc)
# Step 2: Best model is with Site/session as random effects
AIC(modre1b, modre2b, modre3)
anova(modre3, modre2b, modre1b)

# the model with only session as a random effect structure is the best, and has no singularity issue.


## New Full model with new random effect structure----------------------------------------
Pred_sc$PredCount <- Pred_sc$PredRate*10 # creates variable number of eggs predated to solve issue with tab_model


modFullnew <- glmer.nb(PredCount ~ poly(Ldscp, 2)*Treatment + Treatment*Distance + (1|Session),
                  data = Pred_sc)

# model complexity vs. sample size
logLik(modFullnew)  # gives the df of the model (10)
nrow(Pred[!is.na(Pred$PredRate),])/attr(logLik(modFullnew), "df")  # the ratio n/k largely exceeds 10


# check residuals with Dharma: good but outlier test is significant
res1new <- simulateResiduals(modFullnew, plot = T)

# residuals vs. predictors: good but at distance 15m, quartiles don't fully match uniform distrib
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 2))
plotResiduals(res1new, Pred_sc$Ldscp[!is.na(Pred_sc$PredRate)], asFactor = FALSE, main = "Landscape")
plotResiduals(res1new, Pred_sc$Treatment[!is.na(Pred_sc$PredRate)], main = "Treatment")
plotResiduals(res1new, Pred_sc$Distance[!is.na(Pred_sc$PredRate)],  main = "Distance")
plotResiduals(res1new, Pred_sc$Site[!is.na(Pred_sc$PredRate)],  main = "Site")
plotResiduals(res1new, Pred_sc$Session[!is.na(Pred_sc$PredRate)],  main = "Session")
par(op)


# residuals: look good, but quartiles of residuals per Site, DIstance and Session levels slightly differ from uniformity

## Save Full model with new random effect structure----------------------------------------
modfull <- modFullnew

# tab_model(modfull)
Anova(modfull)

saveRDS(modfull, file = "Output/PredRate_FullModelnew.rds")

## MODEL SELECTION with new random effect structure--------------

# step 1
drop1(modfull, test = "Chisq")

# Treatment:Distance is not significant
modsel1 <- update(modfull,  .~. -Treatment:Distance)

# step 2
drop1(modsel1, test = "Chisq")

# Distance ns
modsel2 <- update(modsel1, .~. -Distance)

# step 3
drop1(modsel2, test = "Chisq")

# no further deletion

## Optimal model selected ----------------------------------------------------------------

modOpt <- glmer.nb(PredCount ~ poly(Ldscp, 2)*Treatment + (1|Session),
                                data = Pred_sc)
isSingular(modOpt)

summary(modOpt)

# check residuals with Dharma: outlier test significant (2 values > 1)
resOpt <- simulateResiduals(modOpt, plot = T)

# residuals vs. predictors : no strong sign of var heterogeneity nor patterns with predictors
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 2))
plotResiduals(resOpt, Pred_sc$Ldscp[!is.na(Pred_sc$PredRate)], asFactor = FALSE, main = "Landscape")
plotResiduals(resOpt, Pred_sc$Treatment[!is.na(Pred_sc$PredRate)], main = "Treatment")
plotResiduals(resOpt, Pred_sc$Distance[!is.na(Pred_sc$PredRate)],  main = "Distance")
plotResiduals(resOpt, Pred_sc$Site[!is.na(Pred_sc$PredRate)],  main = "Site")
plotResiduals(resOpt, Pred_sc$Session[!is.na(Pred_sc$PredRate)],  main = "Session")
par(op)


## Save results of optimal model---------------------------------------------------------
saveRDS(modOpt, file = "Output/PredRate_OptimalModel.rds")


## Results ------------------------------------------------------------------------------
modOpt <- readRDS(file = "Output/PredRate_OptimalModel.rds")

tab_model(modOpt)

# Plot showing Landscape:Guild effect on log10(richness)
plot_model(modOpt, type = "pred", terms = c("Ldscp [all]", "Treatment"))
plot_model(modOpt, type = "pred", terms = c("Ldscp [all]"))
plot_model(modOpt, type = "pred", terms = c("Treatment"))





