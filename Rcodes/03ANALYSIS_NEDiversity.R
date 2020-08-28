## Models for natural enemies richness

# The scripts tests the effects of flower strips, landscape on Nat enemies genus richness
# Model selection, model assumption checks, and results stored into output folder.

## Main Analysis - Landscape*FlowerStrips*Community

# Y ~ Landscp*Treatment*Community + Treatment*Distc*Community + (1|Site/Session)

# Landscp = cont. var: gradient of proportion of SNH
# Treatment = factor: flower strip vs. grassy strip
# Community = factor: vine, soil, vegetation
# Distc = factor: 0m, 5m, 20 m or continuous

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
hist(Diversity$TaxaR, main = "Taxa richness")


par(mfrow=c(4,2), mar = c(4,4,2,2))
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

# residuals are uniform, but visually violate variance heterpgeneity, and strong trend with landscape suggest non linear effect needed

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

# deviation due to the distance effect? higher residuals at distance = 15 m.

## Full model5:  Distance as a factor -------------
Div$DistanceFactor <- ifelse(Div$Distance == min(Div$Distance), "0m", 
                             ifelse(Div$Distance==max(Div$Distance), "30m", "15m"))
  
modfin2 <- lmer(log10(TaxaR+1) ~ poly(Ldscp, 2)*Treatment*Guild + DistanceFactor*Treatment*Guild + 
                 (1|Site:session),
               REML = TRUE,
               data = Div,
               control= lmerControl(optimizer="bobyqa"))

## residuals vs. fitted
plot(modfin2)

# check residuals with Dharma
res5 <- simulateResiduals(modfin2, plot = T)

testQuantiles(res5)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res5, Div$Ldscp, asFactor = FALSE, main = "Landscape")
plotResiduals(res5, Div$Guild, main = "Guild")
plotResiduals(res5, Div$Treatment, main = "Treatment")
plotResiduals(res5, Div$Distance,  main = "Distance")
par(mfrow = c(1,1))

# no deviation from uniform distribution detected, no pattern with covariates, 

## SAVE Final Full model ---------------
tab_model(modfin2)

saveRDS(modfin2, file = "Output/NEDiversity_FullModel_TaxoR.rds")



## 3. Model selection taxo richness data--------

# Estimate with ML instead of REML
mod1_taxarFull <- update(modfin2, REML = FALSE)

# first step
drop1(mod1_taxarFull, test = "Chisq")

# The interaction Treatment:Guild:Distance is ns
modsel1_taxar <- update(mod1_taxarFull, .~. -poly(Ldscp, 2):Treatment:Guild)

# step 2
drop1(modsel1_taxar, test = "Chisq")

# interaction Treatment:Guild:DistanceFactor ns
modsel2_taxar <- update(modsel1_taxar, .~. -Treatment:Guild:DistanceFactor)

drop1(modsel2_taxar, test= "Chisq")

# Treatment:distance is the least significant interaction
modsel3_taxar <- update(modsel2_taxar, .~. -Treatment:DistanceFactor)
drop1(modsel3_taxar, test = "Chisq")

# poly(Ldscp, 2):Treatment  is the least significant interaction
modsel4_taxar <- update(modsel3_taxar, .~. -poly(Ldscp, 2):Treatment )
drop1(modsel4_taxar, test = "Chisq")

# guild:distance is the least significant interaction
modsel5_taxar <- update(modsel4_taxar, .~. -Guild:DistanceFactor)
drop1(modsel5_taxar, test = "Chisq")

# Ldscp:Guild is the least significant interaction
modsel6_taxar <- update(modsel5_taxar, .~. -poly(Ldscp, 2):Guild )
drop1(modsel6_taxar, test = "Chisq")

# treatment is the least signficant
modsel7_taxar <- update(modsel6_taxar, .~. -Treatment:Guild)
drop1(modsel7_taxar, test = "Chisq")

# Landscape is the least significant
modsel8_taxar <- update(modsel7_taxar, .~. -poly(Ldscp, 2))
drop1(modsel8_taxar, test = "Chisq")

# Treatment is the least significant
modsel9_taxar <- update(modsel8_taxar, .~. -Treatment)
drop1(modsel9_taxar, test = "Chisq")


# no further deletion

## Validation of optimal model--------

modopt_taxar <- update(modsel9_taxar, REML = TRUE)

# Check model assumptions for the optimal model
## Inspect residuals
plot(modopt_taxar)

# check residuals with Dharma
res <- simulateResiduals(modopt_taxar, plot = T)

# not great

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

# trend with landscape and with guild: deviation from uniform: revise the full model

# 4. New Full models-----------------------------------

# trying out new distribution to fit the data: negative binomial chosen because
# 1. data are counts (number of species), integer values, positive
# 2. the variance is much higher than the mean taxa richness (overdispersion)
mean(Div$TaxaR); var(Div$TaxaR)

# Full model 1b: NB --------------------------------------
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

# residuals are uniform, and have homogeneous variance, no pattern with landscape, but some patternsfound with guild and distance...


## SAVE Final Full model ---------------
tab_model(mod2b_taxar)

saveRDS(mod2b_taxar, file = "Output/NEDiversity_FullModelNB_TaxoR.rds")


## 3. Model selection taxo richness data--------

# Estimate with ML instead of REML: not possible with glmer.nb
# mod1_taxarFull <- update(modfin2, REML = FALSE)

# first step
drop1(mod2b_taxar, test = "Chisq")

# The interaction Treatment:Guild:Distance is ns
modsel1b_taxar <- update(mod2b_taxar, .~. -Treatment:Guild:Distance)

# step 2
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

Anova(modopt_taxar)


# save model results
tab_model(modopt_taxar)

saveRDS(modopt_taxar, file = "Output/NEDiversity_OptimalModel_TaxoR.rds")

# Plot showing Landscape:Guild effect on log10(richness)
plot_model(modopt_taxar, type = "pred", terms = c("Ldscp [all]"))
plot_model(modopt_taxar, type = "pred", terms = c("Guild"))
plot_model(modopt_taxar, type = "pred", terms = c("Guild [all]"))





## Repeat the analysis with genus richness ----------

## Full model 1 -----
mod1_sc <- lmer(GenusR ~ Ldscp*Treatment*Guild + Treatment*Distance*Guild + 
                   (1|Site:session) ,
                 data=Div, control=lmerControl(optimizer="bobyqa"))

# model complexity vs. sample size
# k = 19: 16 (fixed) + 2 (random) + 1 (error term); n = 364
logLik(mod1) # gives the df of the model (19)
nrow(Diversity[!is.na(Diversity$GenusR),])
# the ratio n/k should be between 3 and 10 (Harrison, 2018)

## residuals vs. fitted
plot(mod1_sc)

# check residuals with Dharma
res1 <- simulateResiduals(mod1_sc, plot = T)

# Formal goodness of fit tests
testResiduals(res1)

# residuals vs. predictors : no strong sign of var heterogeneity
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
plotResiduals(res1, scale(Diversity$Ldscp[!is.na(Diversity$GenusR)]), asFactor = FALSE, main = "Landscape")
plotResiduals(res1, Diversity$Guild[!is.na(Diversity$GenusR)], main = "Guild")
plotResiduals(res1, Diversity$Treatment[!is.na(Diversity$GenusR)], main = "Treatment")
plotResiduals(res1, Diversity$Distance[!is.na(Diversity$GenusR)],  main = "Distance")
par(op)

# other plots suggesting non-normality (blue line not so close to the red line)
plot_model(mod1_sc, type = "slope")


## Full model 2: add non-linear pattern with landscape-----
mod2_sc <- lmer(GenusR ~ poly(Ldscp, 2)*Treatment*Guild + Treatment*Distance*Guild + 
               (1|Site/session) ,
             data=Div, control=lmerControl(optimizer="bobyqa"))

# model complexity vs. sample size
# k = 25; n = 364
logLik(mod2) # gives the df of the model (25)
nrow(Diversity[!is.na(Diversity$GenusR),])
# the ratio n/k should be between 3 and 10 (Harrison, 2018)
364/25

## residuals vs. fitted
plot(mod2_sc)

# check residuals with Dharma
res2 <- simulateResiduals(mod2_sc, plot = T)

# Formal goodness of fit tests
testResiduals(res2)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res2, scale(Diversity$Ldscp[!is.na(Diversity$GenusR)]), asFactor = FALSE, main = "Landscape")
plotResiduals(res2, Diversity$Guild[!is.na(Diversity$GenusR)], main = "Guild")
plotResiduals(res2, Diversity$Treatment[!is.na(Diversity$GenusR)], main = "Treatment")
plotResiduals(res2, Diversity$Distance[!is.na(Diversity$GenusR)], main = "Distance")
par(mfrow = c(1,1))

# residuals vs. random factors : plots suggest variance heterogeneity for the different sites
op3 <- par(mfrow=c(1,2), mar = c(4,4,2,2))
plotResiduals(res2, Diversity$Site[!is.na(Diversity$GenusR)])
plotResiduals(res2, Diversity$session[!is.na(Diversity$GenusR)])
par(op3)


## Full model 3 : transform response variable-------
mod3_sc <- lmer(log10(GenusR+1) ~ Ldscp*Treatment*Guild + Treatment*Distance*Guild + 
               (1|Site/session) ,
             data=Div, control=lmerControl(optimizer="bobyqa"))

## residuals vs. fitted
plot(mod3_sc)

# check residuals with Dharma
res3 <- simulateResiduals(mod3_sc, plot = T)

# Formal goodness of fit tests
testResiduals(res3)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res3, scale(Diversity$Ldscp[!is.na(Diversity$GenusR)]), asFactor = FALSE, main = "Landscape")
plotResiduals(res3, Diversity$Guild[!is.na(Diversity$GenusR)], main = "Guild")
plotResiduals(res3, Diversity$Treatment[!is.na(Diversity$GenusR)], main = "Treatment")
plotResiduals(res3, Diversity$Distance[!is.na(Diversity$GenusR)], main = "Distance")
par(mfrow = c(1,1))

# residuals vs. random factors : plots suggest variance heterogeneity for the different sites
op3 <- par(mfrow=c(1,2), mar = c(4,4,2,2))
plotResiduals(res3, Diversity$Site[!is.na(Diversity$GenusR)])
plotResiduals(res3, Diversity$session[!is.na(Diversity$GenusR)])
par(op3)

## Final Full Model ------------------
modfull_sc <- lmer(log10(GenusR+1) ~ poly(Ldscp,2)*Treatment*Guild + Treatment*Distance*Guild + 
               (1|Site/session) ,
             data=Div, control=lmerControl(optimizer="bobyqa"))

## residuals vs. fitted
plot(modfull_sc)

# check residuals with Dharma
res <- simulateResiduals(modfull_sc, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res, scale(Diversity$Ldscp[!is.na(Diversity$GenusR)]), asFactor = FALSE, main = "Landscape")
plotResiduals(res, Diversity$Guild[!is.na(Diversity$GenusR)], main = "Guild")
plotResiduals(res, Diversity$Treatment[!is.na(Diversity$GenusR)], main = "Treatment")
plotResiduals(res, Diversity$Distance[!is.na(Diversity$GenusR)], main = "Distance")
par(mfrow = c(1,1))

# residuals vs. random factors : plots suggest variance heterogeneity for the different sites
op3 <- par(mfrow=c(1,3), mar = c(4,4,2,2))
plotResiduals(res, Diversity$Site[!is.na(Diversity$GenusR)])
plotResiduals(res, Diversity$session[!is.na(Diversity$GenusR)])
plotResiduals(res, Diversity$session[!is.na(Diversity$GenusR)]:Diversity$Site[!is.na(Diversity$GenusR)])
par(op3)

## Save Full model results-------
tab_model(modfull_sc)
Anova(modfull_sc)

saveRDS(modfull_sc, file = "Output/NEDiversity_FullModel.rds")




##----------------------------------------------
## 3. Model selection------------------------------
##----------------------------------------------

# recode model with ML instead of REML
modfull_ml <- lme4::lmer(log10(GenusR+1) ~ poly(Ldscp,2)*Treatment*Guild + Treatment*Distance*Guild + 
                     (1|Site/session) ,
                   REML = FALSE,
                   data=Div[!is.na(Div$GenusR),], 
                   control=lmerControl(optimizer="bobyqa"))

drop1(modfull_ml, test = "Chisq")

# The interaction Ldscp:Treatment:Guild has highest P
modsel1 <- update(modfull_ml, .~. -poly(Ldscp, 2):Treatment:Guild)

# step 2
drop1(modsel1, test = "Chisq")

# three way interaction treatment:guild:distance is ns still
modsel2 <- update(modsel1, .~. -Treatment:Guild:Distance)
drop1(modsel2, test = "Chisq")

# guid:distance is the least significant interaction
modsel3 <- update(modsel2, .~. -Guild:Distance)
drop1(modsel3, test = "Chisq")

# Ldscp:treatment is the least significant interaction
modsel4 <- update(modsel3, .~. -poly(Ldscp, 2):Treatment )
drop1(modsel4, test = "Chisq")

# Treatment:distance is ns
modsel5 <- update(modsel4, .~. -Treatment:Distance)
drop1(modsel5, test = "Chisq")

# treatment;guild is the least significant interaction (F value and p)
modsel6 <- update(modsel5, .~. -Treatment:Guild)
drop1(modsel6, test = "Chisq")

# treatment is least significant
modsel7 <- update(modsel6, .~. -Treatment)
drop1(modsel7, test = "Chisq")

# distance is least significant
modsel8 <- update(modsel7, .~. -Distance)
drop1(modsel8, test = "Chisq")

# no further model simplification

## Model validation ------------------------------------------------

# Check model assumptions for the optimal model
## Inspect residuals
## residuals vs. fitted
modfin_lme4 <- lme4::lmer(log10(GenusR + 1) ~ poly(Ldscp, 2) + Guild + poly(Ldscp, 2):Guild  + (1 | Site/session),
             REML = TRUE,
             data = Div, 
             control=lmerControl(optimizer="bobyqa"))
modfin_lmerT <- lmer(log10(GenusR + 1) ~ poly(Ldscp, 2) + Guild + poly(Ldscp, 2):Guild  + (1 | Site/session),
                          REML = TRUE,
                          data = Div, 
                          control=lmerControl(optimizer="bobyqa"))

plot(modfin_lme4)

# check residuals with Dharma
res <- simulateResiduals(modfin_lme4, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res, scale(Diversity$Ldscp[!is.na(Diversity$GenusR)]), asFactor = FALSE, main = "Landscape")
plotResiduals(res, Diversity$Guild[!is.na(Diversity$GenusR)], main = "Guild")
plotResiduals(res, Diversity$Treatment[!is.na(Diversity$GenusR)], main = "Treatment")
plotResiduals(res, Diversity$Distance[!is.na(Diversity$GenusR)],  main = "Distance")
par(mfrow = c(1,1))




##----------------------------------------
## 4. Results ------------------------------
##----------------------------------------

tab_model(modfin_lme4)

# Plot showing Landscape:Guild effect on log10(richness)
plot_model(modfin_lme4, type = "pred", terms = c("Ldscp [all]", "Guild"))

# save model results----
saveRDS(modfin_lme4, file = "Output/NEDiversity_OptimalModel.rds")




