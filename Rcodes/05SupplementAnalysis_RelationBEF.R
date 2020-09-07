## Script to investigate relationships between predation and natural enemies


## Functions---------------------------------------------------------------------------
library(dplyr)
library(lme4)
library(DHARMa)
library(sjPlot)
library(car)

## Data----------------------------------------------------------------------------
Abundance <- read.csv("Output/AbundanceClean.csv")
Abundance <- read.csv("Output/AbundanceClean.csv")
Abundance$Site <- as.factor(Abundance$Site)

Ldscp <- read.csv("Output/Landscapevars.csv")
Ldscp$Site <- as.factor(as.character(Ldscp$couple))
Ldscp$Ldscp <- Ldscp$HSN1000

Abundance <- left_join(Abundance, Ldscp, by = "Site")
Abundance$Site <- as.factor(Abundance$Site)
Abundance$session <- as.factor(as.character(Abundance$session))

Diversity <- read.csv("Output/DiversityClean.csv")
Diversity <- left_join(Diversity, Ldscp, by = "Site")
Diversity$Site <- as.factor(Diversity$Site)
Diversity$session <- as.factor(as.character(Diversity$session))

Pred <- read.csv("Output/PredationPest_clean.csv")
Pred <- left_join(Pred, Ldscp, by = "Site")
Pred$Site <- as.factor(as.character(Pred$Site))
Pred$Session <- as.factor(as.character(Pred$Session))

## Process data: average predation and abundance across different sessions----------------

# Abundance
datAb <- Abundance %>% group_by(Site, Treatment, Guild, Distance) %>% summarize(AverageAb = mean(Total, na.rm = TRUE))
datAb <- as.data.frame(datAb)

# Diversity
datDiv <- Diversity %>% group_by(Site, Treatment, Guild, Distance) %>% summarize(AverageDiv = mean(TaxaR, na.rm = TRUE))
datDiv <- as.data.frame(datDiv)

# Predation
datPred <- Pred %>% group_by(Site, Treatment, Distance) %>% summarize(AveragePred = mean(PredRate, na.rm = TRUE))
datPred <- as.data.frame(datPred)

## Assemble data
datAbVine <- datAb[datAb$Guild =="Vine",]
datDivVine <- datDiv[datDiv$Guild =="Vine",]

datbef <- left_join(datAbVine, datPred, by = c("Site", "Treatment", "Distance"))
datbef <- left_join(datbef, datDivVine, by = c("Site", "Treatment", "Distance"))
nrow(datbef)

## Relationships -------------------

plot(datbef)
with(datbef, plot(AverageAb, AveragePred, col = Treatment))
with(datbef, plot(AverageDiv, AveragePred, col = factor(Distance)))

## Relationships stats---------------

hist(datbef$AveragePred)
hist(datbef$AverageAb)
hist(datbef$AverageDiv)

# correlation
modAb <- lmer(AveragePred ~ AverageAb + (1|Site/Treatment),
                                  data = datbef)

# residual diagnostic plots
plot(modAb)
res1 <- simulateResiduals(modAb, plot = T)

# residuals vs. predictors: good but at distance 15m
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 2))
plotResiduals(res1, datbef$AverageAb, main = "Abdc")
plotResiduals(res1, datbef$AverageDiv, main = "Div")
plotResiduals(res1, datbef$Site, main = "Site")
plotResiduals(res1, datbef$Distance, main = "Distance")
plotResiduals(res1, datbef$Treatment, main = "Treatment")
par(op)

# patterns with Site and treatment detected (maybe normal because there are random effects?)

modAb2 <- lmer(log10(AveragePred) ~ AverageAb + (1|Site/Treatment),
              data = datbef)

# residual diagnostic plots
plot(modAb2)
res2 <- simulateResiduals(modAb2, plot = T)

# residuals vs. predictors: good but at distance 15m
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 2))
plotResiduals(res2, datbef$AverageAb, main = "Abdc")
plotResiduals(res2, datbef$AverageDiv, main = "Div")
plotResiduals(res2, datbef$Site, main = "Site")
plotResiduals(res2, datbef$Distance, main = "Distance")
plotResiduals(res2, datbef$Treatment, main = "Treatment")
par(op)

# note that the residuals are higher for the high diversity than the low diversity treatment on average

Anova(modAb2)
summary(modAb2)

# interaction barely significant
drop1(update(modAb2, REML = FALSE), test = "Chisq")

## Plot -------------------------
tab_model(modAb2)
plot_model(modAb2, type = "pred")


## What if we average at the plot scale? ------------------

datbef2 <- datbef %>% group_by(Site, Treatment) %>% summarize(AverageAb2 = mean(AverageAb),
                                                              AverageDiv2 = mean(AverageDiv),
                                                              AveragePred2 = mean(AveragePred))

nrow(datbef2)

# exploratory plots
plot(datbef2)
plot(datbef2$AverageAb2, datbef2$AveragePred2, col = datbef2$Treatment)
plot(datbef2$AverageDiv2, datbef2$AveragePred2, col = datbef2$Treatment)

# the treatment effect is huge
boxplot(datbef2$AveragePred2~ datbef2$Treatment)
boxplot(datbef2$AverageAb2~ datbef2$Treatment)
boxplot(datbef2$AverageDiv2~ datbef2$Treatment)

# 
hist(datbef2$AverageAb2)
hist(datbef2$AveragePred2)

#
modAb3 <- lm(log10(AveragePred2) ~ AverageAb2*Treatment, data = datbef2)

plot(modAb3)
summary(modAb3)
Anova(modAb3)

drop1(modAb3, test = "F")
modsel2 <- update(modAb3, .~. -AverageAb2:Treatment )
drop1(modsel2, test = "F")



# sous groupe? composition des communautés versus predation/ par groupe taxo? 
# NMDS/PErmanova taux de predation en variable explicative