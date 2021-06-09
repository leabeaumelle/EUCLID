## Models for predation rates

# The script tests the effects of flower strips, landscape on Lobesia eggs predation
# Model selection, model assumption checks, and results stored into output folder.

## Main Analysis -

# Y ~ Landscp*Treatment + Treatment*Distc + NatEnemiesDensity + (1|Site/Session)

# Landscp = cont. var: gradient of proportion of SNH
# Treatment = factor: flower strip vs. grassy strip
# Distc = continuous variable taking values 0, 15 and 30 m


# site = pairs of plots with same Landscp
# Session = multiple sessions per plot

## Functions---------------------------------------------------------------------------
library(dplyr)
library(lme4)
library(DHARMa)
library(sjPlot)
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
Pred$Session <- as.factor(as.character(Pred$Session))

# Rescale and center continuous predictors: landscape and distance variables
numcols <- grep("Ldscp|Dist",names(Pred))
Pred_sc <- Pred
Pred_sc[,numcols] <- scale(Pred[,numcols])
Pred_sc$PredatedEggs <- Pred_sc$PredRate*10

# GPS location of the sites
Locations <- read.csv("Data/EUCLID - Parcelles.csv")
Locations$Site <- factor(rep(1:12, each = 2))
Locations$Treatment <- rep(c("Low Div", "High Div"), 12)
Locations <- Locations[, c("x", "y", "Site", "Treatment")]

Pred_sc <- left_join(Pred_sc, Locations, by = c("Site", "Treatment"))

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

modFull1 <- glmer(PredRate ~ Ldscp*Treatment + Treatment*Distance + (1|Site:Session),
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
modFull2 <- glmer(PredatedEggs ~ Ldscp*Treatment + Treatment*Distance + (1|Site:Session),
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
modFull3 <- glmer.nb(PredatedEggs ~ Ldscp*Treatment + Treatment*Distance + (1|Site:Session),
                     data = Pred_sc) 

# check residuals with Dharma: better but outlier and quantile deviation detected
res3 <- simulateResiduals(modFull3, plot = T)

# residuals vs. predictors : no strong sign of var heterogeneity but deviation from uniform
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
plotResiduals(res3, Pred_sc$Ldscp[!is.na(Pred_sc$PredRate)], asFactor = FALSE, main = "Landscape")
plotResiduals(res3, Pred_sc$Treatment[!is.na(Pred_sc$PredRate)], main = "Treatment")
plotResiduals(res3, Pred_sc$Distance[!is.na(Pred_sc$PredRate)],  main = "Distance")
par(op)

# patterns with landscape detected

## full model 4: try non linear landscape effects
modFull4 <- glmer.nb(PredatedEggs ~ poly(Ldscp, 2)*Treatment + Treatment*Distance + (1|Site:Session),
                     data = Pred_sc) 

# check residuals with Dharma: much better. One outlier
res4 <- simulateResiduals(modFull4, plot = T)

# residuals vs. predictors : no strong sign of var heterogeneity
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
plotResiduals(res4, Pred_sc$Ldscp[!is.na(Pred_sc$PredRate)], asFactor = FALSE, main = "Landscape")
plotResiduals(res4, Pred_sc$Treatment[!is.na(Pred_sc$PredRate)], main = "Treatment")
plotResiduals(res4, Pred_sc$Distance[!is.na(Pred_sc$PredRate)],  main = "Distance")
par(op)


# Model with negative binomial and polynomial 2 effect landscape is the best so far

# Spatial autocorrelation ? 
plotResiduals(res4$scaledResiduals, Pred_sc$x[!is.na(Pred_sc$PredRate)])
plotResiduals(res4$scaledResiduals, Pred_sc$y[!is.na(Pred_sc$PredRate)])


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

# check residuals with Dharma: outlier test significant (2 values > 1)
resOpt <- simulateResiduals(modOpt, plot = T)

# residuals vs. predictors : no strong sign of var heterogeneity nor patterns with predictors
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
plotResiduals(resOpt, Pred_sc$Ldscp[!is.na(Pred_sc$PredRate)], asFactor = FALSE, main = "Landscape")
plotResiduals(resOpt, Pred_sc$Treatment[!is.na(Pred_sc$PredRate)], main = "Treatment")
plotResiduals(resOpt, Pred_sc$Distance[!is.na(Pred_sc$PredRate)],  main = "Distance")
par(op)


## Save results of optimal model---------------------------------------------------------
saveRDS(modOpt, file = "Output/PredRate_OptimalModel.rds")


tab_model(modOpt)

plot_model(modOpt, type = "pred", terms = c("Ldscp [all]", "Treatment"))
plot_model(modOpt, type = "pred", terms = c("Ldscp [all]"))
plot_model(modOpt, type = "pred", terms = c("Treatment"))


## Revisions -------------------

# R2: effect sizes seem small: is it relevant to farmers? 

# Calculate effect sizes for each landscape based on model predictions
library(ggeffects)
library(dplyr)

## Data
Pred <- read.csv("Output/PredationPest_clean.csv")
Ldscp <- read.csv("Output/Landscapevars.csv")
Ldscp$Site <- as.factor(as.character(Ldscp$couple))
Ldscp$Ldscp <- Ldscp$HSN1000

Pred$Site <- as.factor(as.character(Pred$Site))
Pred <- left_join(Pred, Ldscp, by = "Site")
Pred$Site <- as.factor(as.character(Pred$Site))
Pred$Session <- as.factor(as.character(Pred$Session))


# Rescale and center continuous predictors: landscape and distance variables
numcols <- grep("Ldscp|Dist",names(Pred))
Pred_sc <- Pred
Pred_sc[,numcols] <- scale(Pred[,numcols])
Pred_sc$PredatedEggs <- Pred_sc$PredRate*10

# Load model results from scripts 03_
modPred <- readRDS(file = "Output/PredRate_OptimalModel.rds")
modFullPred <- readRDS(file = "Output/PredRate_FullModel.rds")


# get the predictions with ggeffects
me <- ggemmeans(modFullPred, c("Ldscp [all]", "Treatment"))

# take the raw data from ggeffects and compute means and sds
raw <- attr(me, "rawdata")
raw2 <- raw %>% group_by(x, group) %>% summarise(means = mean(response, na.rm = TRUE), 
                                                 sdev = sd(response, na.rm = TRUE))


# effect sizes (LRR): 
LRR <- log(raw2$means[raw2$group=="High Div"]/raw2$means[raw2$group=="Low Div"])

# in percent change: 
PercentChange <- 100 * (exp(LRR) - 1)
