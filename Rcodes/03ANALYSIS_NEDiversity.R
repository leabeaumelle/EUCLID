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


## Functions---------------------------------------------------------------------------
library(dplyr)
library(lme4)
library(MuMIn)
library(DHARMa)
library(MASS)
library(sjPlot)
library(lattice)
library(optimx)
library(lmerTest)

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

# Rescale and center continuous predictors: landscape and distance variables
numcols <- grep("Ldscp|Dist",names(Diversity))
Div <- Diversity
Div[,numcols] <- scale(Div[,numcols])

## Data exploration--------------------
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


## Modelling-----------------------------------
## Full model 1 -----
mod1 <- lmer(GenusR ~ Ldscp*Treatment*Guild + Treatment*Distance*Guild + 
                   (1|Site/session) ,
                 data=Diversity, control=lmerControl(optimizer="bobyqa"))

# model complexity vs. sample size
# k = 19: 16 (fixed) + 2 (random) + 1 (error term); n = 364
logLik(mod1) # gives the df of the model (19)
nrow(Diversity[!is.na(Diversity$GenusR),])
# the ratio n/k should be between 3 and 10 (Harrison, 2018)
logLik(mod1)/nrow(Diversity[!is.na(Diversity$GenusR),])

mod1_sc <- update(mod1,data=Div)

## Inspect residuals-----
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

# residuals vs. random factors : plots suggest variance heterogeneity for the different sites
op2 <- par(mfrow=c(1,2), mar = c(4,4,2,2))
plotResiduals(res1, Diversity$Site[!is.na(Diversity$GenusR)])
plotResiduals(res1, Diversity$session[!is.na(Diversity$GenusR)])
par(op2)

# plot suggesting non-normality (blue line not so close to the red line)
plot_model(mod1_sc, type = "slope")

## Full model 2: add non-linear pattern with landscape-----
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

## Inspect residuals------
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
mod3 <- lmer(log10(GenusR+1) ~ Ldscp*Treatment*Guild + Treatment*Distance*Guild + 
               (1|Site/session) ,
             data=Diversity, control=lmerControl(optimizer="bobyqa"))

# Rescale and center continuous predictors: landscape and distance variables
mod3_sc <- update(mod3,data=Div)

## Inspect residuals-------
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
modfull <- lmer(log10(GenusR+1) ~ poly(Ldscp,2)*Treatment*Guild + Treatment*Distance*Guild + 
               (1|Site/session) ,
             data=Diversity, control=lmerControl(optimizer="bobyqa"))

# Rescale and center continuous predictors: landscape and distance variables
modfull_sc <- update(modfull,data=Div)

# Inspect residuals ------
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

saveRDS(modfull_sc, file = "Output/NEDiversity_FullModel.rds")


## Model selection-------------------------------------------------------------
# recode model with ML instead of REML
modfull_ml <- lmer(log10(GenusR+1) ~ poly(Ldscp,2)*Treatment*Guild + Treatment*Distance*Guild + 
                     (1|Site/session) ,
                   REML = FALSE,
                   data=Div, 
                   control=lmerControl(optimizer="bobyqa"))

drop1(modfull_ml, test = "Chisq")

# drop1 renders a F test based on type III analysis of variance, with Sattherthwaite's method for DF
anova(modsel1)
# shows that both methods give the same results
anova(modsel1, ddf = "Kenward-Roger")

# The interaction Treatment:Guild:Distance is ns, delta AIC is 2, and LRT is very small
modsel1 <- update(modfull_ml, .~. -Treatment:Guild:Distance)

# step 2
drop1(modsel1)

# three way interaction landscape:treatment:guild is ns still
modsel2 <- update(modsel1, .~. -poly(Ldscp, 2):Treatment:Guild)

anova(modsel2)
drop1(modsel2)

# guid:distance is the least significant interaction (F value and p)
modsel3 <- update(modsel2, .~. -Guild:Distance)
drop1(modsel3)

# treatment:distance is the least significant interaction (F value and p)
modsel4 <- update(modsel3, .~. -Treatment:Distance)
drop1(modsel4)

# landscape:treatment is the least significant interaction (F value and p)
modsel5 <- update(modsel4, .~. -poly(Ldscp, 2):Treatment)
drop1(modsel5)

# treatment;guild is the least significant interaction (F value and p)
modsel6 <- update(modsel5, .~. -Treatment:Guild)
drop1(modsel6)

# treatment is least significant
modsel7 <- update(modsel6, .~. -Treatment)
drop1(modsel7)

# distance is least significant
modsel8 <- update(modsel7, .~. -Distance)
drop1(modsel8)

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


# save model results
tab_model(modfin_lme4)
# tab_model(mod8_sc, p.val = "kr") #more precise p-values but require model to be fitted with REML

saveRDS(modfin_lme4, file = "Output/NEDiversity_OptimalModel.rds")

# Plot showing Landscape:Guild effect on log10(richness)
plot_model(modfin_lme4, type = "pred", terms = c("Ldscp [all]", "Guild"))


## Repeat the analysis removing vegetation data ------------------------------------------
# these data confound the distance effect because
# no vegetation sampling was carried out at distances 15 and 30m

DivNoVg <- Div %>% filter(Guild != "Vegetation")
DivNoVg <- droplevels(DivNoVg)


## Full model without vegetation data-----------
mod1_novg <- lme4::lmer(GenusR ~ poly(Ldscp, 2)*Treatment*Guild + Distance*Treatment*Guild + 
                        (1|Site/session),
                        REML = TRUE,
                      data = DivNoVg,
                      control= lmerControl(optimizer="bobyqa"))


# Check model assumptions for the optimal model
## Inspect residuals
## residuals vs. fitted
plot(mod1_novg)

# check residuals with Dharma
res <- simulateResiduals(mod1_novg, plot = T)

# signs of heteroscedasticity

# Log transform variable
mod1_novg2 <- lme4::lmer(log10(GenusR+1) ~ poly(Ldscp, 2)*Treatment*Guild + Distance*Treatment*Guild + 
                          (1|Site/session),
                        REML = TRUE,
                        data = DivNoVg,
                        control= lmerControl(optimizer="bobyqa"))


# Check model assumptions for the optimal model
## Inspect residuals
## residuals vs. fitted
plot(mod1_novg2)

# check residuals with Dharma
res <- simulateResiduals(mod1_novg2, plot = T)

# no signs of heteroscedasticity

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res, scale(DivNoVg$Ldscp[!is.na(DivNoVg$GenusR)]), asFactor = FALSE, main = "Landscape")
plotResiduals(res, DivNoVg$Guild[!is.na(DivNoVg$GenusR)], main = "Guild")
plotResiduals(res, DivNoVg$Treatment[!is.na(DivNoVg$GenusR)], main = "Treatment")
plotResiduals(res, DivNoVg$Distance[!is.na(DivNoVg$GenusR)],  main = "Distance")
par(mfrow = c(1,1))

# adjusted quantile test significant for landscape: strange pattern in residuals vs. landscape

# Remove non linear effect
mod1_novg3 <- lme4::lmer(log10(GenusR+1) ~ Ldscp*Treatment*Guild + Distance*Treatment*Guild + 
                           (1|Site/session),
                         REML = TRUE,
                         data = DivNoVg,
                         control= lmerControl(optimizer="bobyqa"))


# Check model assumptions for the optimal model
## Inspect residuals
## residuals vs. fitted
plot(mod1_novg3)

# check residuals with Dharma
res <- simulateResiduals(mod1_novg3, plot = T)

# no signs of heteroscedasticity, but slight deviation from uniform distrib (lower quantiles..)

# Formal goodness of fit tests
testResiduals(res)

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res, scale(DivNoVg$Ldscp[!is.na(DivNoVg$GenusR)]), asFactor = FALSE, main = "Landscape")
plotResiduals(res, DivNoVg$Guild[!is.na(DivNoVg$GenusR)], main = "Guild")
plotResiduals(res, DivNoVg$Treatment[!is.na(DivNoVg$GenusR)], main = "Treatment")
plotResiduals(res, DivNoVg$Distance[!is.na(DivNoVg$GenusR)],  main = "Distance")
par(mfrow = c(1,1))

# patterns increase if landscape effect is linear. Keep the non-linear trend

## Save No vegetation Full model results-------
tab_model(mod1_novg2)

saveRDS(mod1_novg2, file = "Output/NEDiversity_FullModel.rds")

## Model selection no vegetation data--------

# Estimate with ML instead of REML
mod1_novgFull <- update(mod1_novg2, REML = FALSE)

# first step
drop1(mod1_novgFull, test = "Chisq")

# The interaction Treatment:Guild:Distance is ns : SAME AS BEFORE
modsel1_novg <- update(mod1_novgFull, .~. -Treatment:Guild:Distance)

# step 2
drop1(modsel1_novg, test = "Chisq")

# interaction Guild:Distance not significant: DIFFERENT THAN BEFORE (landscape:treatment:guild)
modsel2_novg <- update(modsel1_novg, .~. -Guild:Distance)

drop1(modsel2_novg, test= "Chisq")

# landscape:treatment:guild is the least significant interaction: NOW SAME AS BEFORE
modsel3_novg <- update(modsel2_novg, .~. -poly(Ldscp, 2):Treatment:Guild)
drop1(modsel3_novg, test = "Chisq")

# treatment:distance is the least significant interaction (based on AIC and LRT values): SAME AS BEFORE
modsel4_novg <- update(modsel3_novg, .~. -Treatment:Distance)
drop1(modsel4_novg, test = "Chisq")

# landscape:treatment is the least significant interaction (based on LRT, not AIC): SAME AS BEFORE
modsel5_novg <- update(modsel4_novg, .~. -poly(Ldscp, 2):Treatment)
drop1(modsel5_novg, test = "Chisq")

# treatment;guild is the least significant interaction: SAME AS BEFORE
modsel6_novg <- update(modsel5_novg, .~. -Treatment:Guild)
drop1(modsel6_novg, test = "Chisq")

# treatment is least significant: SAME AS BEFORE
modsel7_novg <- update(modsel6_novg, .~. -Treatment)
drop1(modsel7_novg, test = "Chisq")

# distance is least significant: SAME AS BEFORE
modsel8_novg <- update(modsel7_novg, .~. -Distance)
drop1(modsel8_novg, test = "Chisq")

# no further deletion

## Validation of optimal model--------

modopt_novg <- update(modsel8_novg, REML = TRUE)

## Model validation ------------------------------------------------

# Check model assumptions for the optimal model
## Inspect residuals
plot(modopt_novg)

# check residuals with Dharma
res <- simulateResiduals(modopt_novg, plot = T)

# Formal goodness of fit tests
testResiduals(res)

# distribution not great, but K-S test is not statistically significant

# residuals vs. predictors
par(mfrow = c(2,2))
plotResiduals(res, scale(DivNoVg$Ldscp[!is.na(DivNoVg$GenusR)]), asFactor = FALSE, main = "Landscape")
plotResiduals(res, DivNoVg$Guild[!is.na(DivNoVg$GenusR)], main = "Guild")
plotResiduals(res, DivNoVg$Treatment[!is.na(DivNoVg$GenusR)], main = "Treatment")
plotResiduals(res, DivNoVg$Distance[!is.na(DivNoVg$GenusR)],  main = "Distance")
par(mfrow = c(1,1))


# save model results
tab_model(modopt_novg)

saveRDS(modopt_novg, file = "Output/NEDiversity_OptimalModel_NoVegetation.rds")

# Plot showing Landscape:Guild effect on log10(richness)
plot_model(modopt_novg, type = "pred", terms = c("Ldscp [all]", "Guild"))

## Compare both models with and without vegetation-------------------
tab_model(modopt_novg, modfin_lme4)
