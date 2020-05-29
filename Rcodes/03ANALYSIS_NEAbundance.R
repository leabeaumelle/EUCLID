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
Abundance$session <- as.factor(as.character(Abundance$session))

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

## LRT test global model
mod.null <-  glmer.nb(Total ~ 1 + (1|Site/session),
                      data=Abs, control=glmerControl(optimizer="bobyqa"))

anova(mod1_sc, mod.null)

# Save results of the full model
tab_model(mod1_sc)

saveRDS(mod1_sc, file = "Output/NEAbundance_FullModel.rds")


## Model selection------------------------
drop1(mod1_sc, test = "Chisq")

# Three way interactions are not significant: dropping the Treatment:Guild:Distance first
mod2_sc <- glmer.nb(Total ~ Ldscp*Treatment*Guild + Distance + 
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="bobyqa"))
isSingular(mod2_sc) # warning: singular fit , but isSingular indicates none of the random effects covariance matrices are singular

drop1(mod2_sc, test = "Chisq")

# Three way interaction Ld:Treatment:Guild still not significant: dropping it
mod3_sc <- glmer.nb(Total ~ Ldscp*Treatment+ Ldscp*Guild + Treatment*Guild + Distance + 
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="bobyqa"))
isSingular(mod3_sc) # warning: singular fit , but isSingular indicates none of the random effects covariance matrices are singular

drop1(mod3_sc, test = "Chisq")

# least significant term is the landscape: treatment interaction
mod4_sc <- glmer.nb(Total ~ Ldscp*Guild + Treatment*Guild + Distance + 
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="bobyqa"))
isSingular(mod4_sc)

summary(mod4_sc)

drop1(mod4_sc, test = "Chisq")

# guild:treatment ns
mod5_sc <- glmer.nb(Total ~ Ldscp*Guild + Treatment + Distance + 
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="bobyqa"))
isSingular(mod5_sc)

summary(mod5_sc)

drop1(mod5_sc, test = "Chisq")

# ldscp:guild interaction ns
mod6_sc <- glmer.nb(Total ~ Ldscp + Guild + Treatment + Distance + 
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="bobyqa"))
isSingular(mod6_sc)

summary(mod6_sc)

drop1(mod6_sc, test = "Chisq")

# ldscp ns
mod7_sc <- glmer.nb(Total ~  Guild + Treatment + Distance + 
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="bobyqa"))
isSingular(mod7_sc)

summary(mod7_sc)

drop1(mod7_sc, test = "Chisq")

# distance ns
mod8_sc <- glmer.nb(Total ~  Guild + Treatment +  
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="bobyqa"))
isSingular(mod8_sc)

summary(mod8_sc)

drop1(mod8_sc, test = "Chisq")

# no further variables to be dropped

# Check model assumptions for the optimal model
## Inspect residuals
## residuals vs. fitted
plot(mod8_sc)

# check residuals with Dharma
res <- simulateResiduals(mod8_sc, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(scale(Abundance$Ldscp), res$scaledResiduals, asFactor = FALSE, main = "Landscape")
plotResiduals(Abundance$Guild, res$scaledResiduals, main = "Guild")
plotResiduals(Abundance$Treatment, res$scaledResiduals, main = "Treatment")
plotResiduals(Abundance$Distance, res$scaledResiduals, main = "Distance")
par(mfrow = c(1,1))


# save model results
tab_model(mod8_sc)
# tab_model(mod8_sc, p.val = "kr") #more precise p-values but require model to be fitted with REML

saveRDS(mod8_sc, file = "Output/NEAbundance_OptimalModel.rds")

# get mean values predicted
plot_model(mod8_sc, type = "pred")


## Repeat the analysis removing vegetation data ------------------------------------------
# these data confound the distance effect because
# no vegetation sampling was carried out at distances 15 and 30m

## Full model without vegetation data-----------
data_novg <- Abs[which(Abs$Guild != "Vegetation"),]

mod1_novg <- glmer.nb(Total ~ Ldscp*Treatment*Guild + Distance*Treatment*Guild + 
                        (1|Site/session),
                      data = data_novg,
                      control=glmerControl(optimizer="bobyqa"))
isSingular(mod1_novg)

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
plotResiduals(data_novg$Ldscp, res$scaledResiduals, asFactor = FALSE, main = "Landscape")
plotResiduals(data_novg$Guild, res$scaledResiduals, main = "Guild")
plotResiduals(data_novg$Treatment, res$scaledResiduals, main = "Treatment")
plotResiduals(data_novg$Distance, res$scaledResiduals, main = "Distance")
par(mfrow = c(1,1))

# signs of non-linear effect of landscape in the residuals

# save model results
tab_model(mod1_novg)
# tab_model(mod1_novg, p.val = "kr") #more precise p-values but require model to be fitted with REML

saveRDS(mod1_novg, file = "Output/NEAbundance_FullModel_novg.rds")


## Model selection without vegetation dat------------------------
drop1(mod1_novg, test = "Chisq")

# Three way interactions are not significant: dropping the Treatment:Guild:Distance first [same as previous]
mod2_novg <- glmer.nb(Total ~ Ldscp*Treatment*Guild + Distance + 
                      (1|Site/session) ,
                    data=data_novg, control=glmerControl(optimizer="bobyqa"))
isSingular(mod2_novg) # warning: singular fit , but isSingular indicates none of the random effects covariance matrices are singular

drop1(mod2_novg, test = "Chisq")

# Three way interaction Ld:Treatment:Guild still not significant: dropping it [same as previous]
mod3_novg <- glmer.nb(Total ~ Ldscp*Treatment+ Ldscp*Guild + Treatment*Guild + Distance + 
                      (1|Site/session) ,
                    data=data_novg, control=glmerControl(optimizer="bobyqa"))
isSingular(mod3_novg) # warning: singular fit , but isSingular indicates none of the random effects covariance matrices are singular

drop1(mod3_novg, test = "Chisq")

# least significant term is the treatment:guild interaction [different than previous: dropped landscape:treatment]
mod4_novg <- glmer.nb(Total ~ Ldscp*Treatment + Ldscp*Guild + Distance + 
                      (1|Site/session) ,
                    data=data_novg, control=glmerControl(optimizer="bobyqa"))
isSingular(mod4_novg)

summary(mod4_novg)

drop1(mod4_novg, test = "Chisq")

# least significant term is the landscape:treatment interaction [same as previous now]
mod5_novg <- glmer.nb(Total ~ Ldscp*Guild + Treatment + Distance + 
                      (1|Site/session) ,
                    data=data_novg, control=glmerControl(optimizer="bobyqa"))
isSingular(mod5_novg)

summary(mod5_novg)

drop1(mod5_novg, test = "Chisq")

# ldscp:guild interaction ns [same as previous]
mod6_novg <- glmer.nb(Total ~ Ldscp + Guild + Treatment + Distance + 
                      (1|Site/session) ,
                    data=data_novg, control=glmerControl(optimizer="bobyqa"))
isSingular(mod6_novg)

summary(mod6_novg)

drop1(mod6_novg, test = "Chisq")

# ldscp ns [same as previous]
mod7_novg <- glmer.nb(Total ~  Guild + Treatment + Distance + 
                      (1|Site/session) ,
                    data=data_novg, control=glmerControl(optimizer="bobyqa"))
isSingular(mod7_novg)

summary(mod7_novg)

drop1(mod7_novg, test = "Chisq")

# distance ns
mod8_novg <- glmer.nb(Total ~  Guild + Treatment +  
                      (1|Site/session) ,
                    data=data_novg, control=glmerControl(optimizer="bobyqa"))
isSingular(mod8_novg)

summary(mod8_novg)

drop1(mod8_novg, test = "Chisq")

# no further variables to be dropped

# Check model assumptions for the optimal model
## Inspect residuals
## residuals vs. fitted
plot(mod8_novg)

# check residuals with Dharma
res <- simulateResiduals(mod8_novg, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(scale(data_novg$Ldscp), res$scaledResiduals, asFactor = FALSE, main = "Landscape")
plotResiduals(data_novg$Guild, res$scaledResiduals, main = "Guild")
plotResiduals(data_novg$Treatment, res$scaledResiduals, main = "Treatment")
plotResiduals(data_novg$Distance, res$scaledResiduals, main = "Distance")
par(mfrow = c(1,1))


# save model results
tab_model(mod8_novg)
# tab_model(mod8_novg, p.val = "kr") #more precise p-values but require model to be fitted with REML

saveRDS(mod8_novg, file = "Output/NEAbundance_OptimalModel_novg.rds")

# get mean values predicted
plot_model(mod8_novg, type = "pred")



## Explore Non linear trend with landscape------------------------------------------------

# adding a polynomial term for landscape
mod1_nl <- glmer.nb(Total ~ poly(Ldscp, 2)*Treatment*Guild + Treatment*Distance*Guild + 
                   (1|Site/session) ,
                 data=Abs, control=glmerControl(optimizer="bobyqa"))

# LRT test of quadratic term: highly significant
anova(mod1_nl, mod1_sc)

# # Check model assumptions: not great
## Inspect residuals
## residuals vs. fitted
plot(mod1_nl)

# check residuals with Dharma
res <- simulateResiduals(mod1_nl, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(Abs$Ldscp, res$scaledResiduals, asFactor = FALSE, main = "Landscape")
plotResiduals(Abs$Guild, res$scaledResiduals, main = "Guild")
plotResiduals(Abs$Treatment, res$scaledResiduals, main = "Treatment")
plotResiduals(Abs$Distance, res$scaledResiduals, main = "Distance")
par(mfrow = c(1,1))

## Model selection

drop1(mod1_nl, test = "Chisq")

#both three way interactions are ns : dropping treatment:guild:distance[same as before]
mod2_nl <- glmer.nb(Total ~ poly(Ldscp,2)*Treatment*Guild + Distance + 
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="bobyqa"))
isSingular(mod2_nl) # warning: singular fit , but isSingular indicates none of the random effects covariance matrices are singular

drop1(mod2_nl, test = "Chisq")

# Three way interaction Ld:Treatment:Guild still not significant: dropping it [same as before]
mod3_nl <- glmer.nb(Total ~ poly(Ldscp,2)*Treatment+ poly(Ldscp,2)*Guild + Treatment*Guild + Distance + 
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="bobyqa"))
isSingular(mod3_nl) # warning: singular fit , but isSingular indicates none of the random effects covariance matrices are singular

drop1(mod3_nl, test = "Chisq")

# least significant term is the Treatment:Guild interaction [different than before]
mod4_nl <- glmer.nb(Total ~ poly(Ldscp,2)*Treatment+ poly(Ldscp,2)*Guild + Distance + 
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="bobyqa"))
isSingular(mod4_nl)

drop1(mod4_nl, test = "Chisq")

# landscape:treatment effect is not significant (p=0.08) [different than before, but same as previous step]
mod5_nl <- glmer.nb(Total ~ poly(Ldscp,2)*Guild + Treatment + Distance + 
                      (1|Site/session) ,
                    data=Abs, control=glmerControl(optimizer="bobyqa"))
isSingular(mod5_nl)

summary(mod5_nl)

drop1(mod5_nl, test = "Chisq")

# no further variables to be dropped (although we are checking on the boundary and p = 0.03 for distance)

## Checking otimal model assumptions
## Inspect residuals
## residuals vs. fitted
plot(mod5_nl)

# check residuals with Dharma
res <- simulateResiduals(mod5_nl, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(scale(Abundance$Ldscp), res$scaledResiduals, asFactor = FALSE, main = "Landscape")
plotResiduals(Abundance$Guild, res$scaledResiduals, main = "Guild")
plotResiduals(Abundance$Treatment, res$scaledResiduals, main = "Treatment")
plotResiduals(Abundance$Distance, res$scaledResiduals, main = "Distance")
par(mfrow = c(1,1))


## Results
tab_model(mod5_nl)
anova(mod5_nl)


# different resposne to landscape depending on guild
plot_model(mod5_nl, type = "pred", terms = c("Ldscp", "Guild"))
plot_model(mod5_nl, type = "pred", terms = c("Distance", "Guild"))
plot_model(mod5_nl, type = "pred", terms = c("Treatment", "Guild"))

# From Arthur-----------------------------------------------------------------------------
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


