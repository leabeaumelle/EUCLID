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
                    (1|Site/session) ,
                  data=Abs, control=glmerControl(optimizer="bobyqa"))

# model complexity vs. sample size
# k = 19: 16 (fixed) + 2 (random) + 1 (error term); n = 366
logLik(mod1_sc) # gives the df of the model (19)
nrow(Abundance)
# the minimum ratio n/k should be between 3 and 10 (Harrison, 2018)
366/19


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
plotResiduals(res, scale(Abundance$Ldscp), asFactor = FALSE, main = "Landscape")
plotResiduals(res, Abundance$Guild, main = "Guild")
plotResiduals(res, Abundance$Treatment, main = "Treatment")
plotResiduals(res, Abundance$Distance,  main = "Distance")
par(mfrow = c(1,1))

## Patterns with landscape: need to add a non-linear effect of the landscape


## Full model with Non linear trend with landscape----

# optimizer bobyqa has convergence issues: so I use Nelder_Mead
mod1_nl <- glmer.nb(Total ~ poly(Ldscp, 2)*Treatment*Guild + Treatment*Distance*Guild + 
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="Nelder_Mead",
                                                   optCtrl=list(maxfun=1e4)))
# testing new re structure
# mod1_nl <- glmer.nb(Total ~ poly(Ldscp, 2)*Treatment*Guild + Treatment*Distance*Guild + 
#                       (1|Site:session) ,
#                     data=Abs, control=glmerControl(optimizer="Nelder_Mead",
#                                                    optCtrl=list(maxfun=1e4)))


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
modsel3 <- update(modsel2, .~. -Treatment:Distance) # convergence warning...

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

# Ldscp:Treatment is ns
modsel6 <- update(modsel5, .~. -poly(Ldscp, 2):Treatment)

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
plotResiduals(res$scaledResiduals, scale(Abundance$Ldscp),  asFactor = FALSE, main = "Landscape")
plotResiduals(res$scaledResiduals, Abundance$Guild, main = "Guild")
plotResiduals(res$scaledResiduals, Abundance$Treatment, main = "Treatment")
plotResiduals(res$scaledResiduals, Abundance$Distance,  main = "Distance")
par(mfrow = c(1,1))


## Results
tab_model(modOpt)
Anova(modOpt)


# different resposne to landscape depending on guild
plot_model(modOpt, type = "pred", terms = c("Treatment", "Guild"))
plot_model(modOpt, type = "pred", terms = c("Ldscp [all]", "Treatment"))
plot_model(modOpt, type = "pred", terms = c("Ldscp [all]", "Guild"))
plot_model(modOpt, type = "pred", terms = c("Distance [all]", "Guild"))



## Repeat the analysis removing vegetation data ------------------------------------------
# these data confound the distance effect because
# no vegetation sampling was carried out at distances 15 and 30m

data_novg <- droplevels(Abs[which(Abs$Guild != "Vegetation"),])


## Full model without vegetation data-----------
mod1_novg <- glmer.nb(Total ~ Ldscp*Treatment*Guild + Distance*Treatment*Guild + 
                        (1|Site/session),
                      data = data_novg,
                      control=glmerControl(optimizer="bobyqa"))

# Check model assumptions for the optimal model
## Inspect residuals
## residuals vs. fitted
plot(mod1_novg)

# check residuals with Dharma
res <- simulateResiduals(mod1_novg, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res$scaledResiduals, data_novg$Ldscp, asFactor = FALSE, main = "Landscape")
plotResiduals(res$scaledResiduals, data_novg$Guild,  main = "Guild")
plotResiduals(res$scaledResiduals, data_novg$Treatment,  main = "Treatment")
plotResiduals(res$scaledResiduals, data_novg$Distance,  main = "Distance")
par(mfrow = c(1,1))

# signs of non-linear effect of landscape in the residuals
mod1_novg2 <- glmer.nb(Total ~ poly(Ldscp,2)*Treatment*Guild + Distance*Treatment*Guild + 
                        (1|Site/session),
                      data = data_novg,
                      control=glmerControl(optimizer="optimx", optCtrl = list(method="nlminb")))

# check n/k ratio
nrow(data_novg)/attr(logLik(mod1_novg2), "df")

# Check model assumptions for the optimal model
## Inspect residuals
## residuals vs. fitted
plot(mod1_novg2)

# check residuals with Dharma
res <- simulateResiduals(mod1_novg2, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res$scaledResiduals, data_novg$Ldscp, asFactor = FALSE, main = "Landscape")
plotResiduals(res$scaledResiduals, data_novg$Guild,  main = "Guild")
plotResiduals(res$scaledResiduals, data_novg$Treatment,  main = "Treatment")
plotResiduals(res$scaledResiduals, data_novg$Distance,  main = "Distance")
par(mfrow = c(1,1))

# Final full model
modFull_novg <- mod1_novg2

# save model results
tab_model(modFull_novg)
Anova(modFull_novg)

# compare with full model with vegetation data
Anova(mod1_nl)

# store the model
saveRDS(modFull_novg, file = "Output/NEAbundance_FullModel_novg.rds")


## Model selection without vegetation dat------------------------

#step 1
drop1(modFull_novg, test = "Chisq")
# dropping the Treatment:Guild:Distance first [same as previous]
modsel1_novg <- update(modFull_novg, .~. -Treatment:Guild:Distance)

# step 2
drop1(modsel1_novg, test = "Chisq")
#  Ld:Treatment:Guild still not significant: dropping it [same as previous]
modsel2_novg <- update(modsel1_novg, .~. -poly(Ldscp, 2):Treatment:Guild)


# step 3
drop1(modsel2_novg, test = "Chisq")
# Treatment:distance interaction [same as before]
modsel3_novg <- update(modsel2_novg, .~. -Treatment:Distance)

# step 4
drop1(modsel3_novg, test = "Chisq")
# Guild:distance interaction ns [same as before]
modsel4_novg <- update(modsel3_novg, .~. -Guild:Distance)

# step 5
drop1(modsel4_novg, test = "Chisq")
# Treatment:Guild interaction ns [same as before]
modsel5_novg <- update(modsel4_novg, .~. -Treatment:Guild)

# step 6
drop1(modsel5_novg, test = "Chisq")
# Ldscp:Guild interaction ns [different than before!]
modsel6_novg <- update(modsel5_novg, .~. -poly(Ldscp, 2):Guild)

# step 7 
drop1(modsel6_novg, test = "Chisq")

# no further deletion

## Optimal model without vegetation
modOpt_novg <- modsel6_novg


# Check model assumptions for the optimal model
## Inspect residuals
## residuals vs. fitted
plot(modOpt_novg)

# check residuals with Dharma
res <- simulateResiduals(modOpt_novg, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res$scaledResiduals, scale(data_novg$Ldscp), asFactor = FALSE, main = "Landscape")
plotResiduals(res$scaledResiduals, data_novg$Guild, main = "Guild")
plotResiduals(res$scaledResiduals, data_novg$Treatment,  main = "Treatment")
plotResiduals(res$scaledResiduals, data_novg$Distance,  main = "Distance")
par(mfrow = c(1,1))


# save model results
tab_model(modOpt_novg)
Anova(modOpt_novg)
saveRDS(modOpt_novg, file = "Output/NEAbundance_OptimalModel_novg.rds")

# get mean values predicted
plot_model(modOpt_novg, type = "pred", terms = "Ldscp [all]")

# same resposne to landscape depending on guild
plot_model(modOpt_novg, type = "pred", terms = c("Ldscp [all]", "Guild", "Treatment"))
plot_model(modOpt_novg, type = "pred", terms = c("Distance [all]", "Guild"))
plot_model(modOpt_novg, type = "pred", terms = c("Treatment", "Guild"))
