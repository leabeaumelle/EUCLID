## Create table to show model results



## Functions----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)
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


# Rescale and center continuous predictors: landscape and distance variables
numcols <- grep("Ldscp|Dist",names(Abundance))
Abs <- Abundance
Abs[,numcols] <- scale(Abs[,numcols])

numcols <- grep("Ldscp|Dist",names(Diversity))
Div <- Diversity
Div[,numcols] <- scale(Div[,numcols])

numcols <- grep("Ldscp|Dist",names(Pred))
Pred_sc <- Pred
Pred_sc[,numcols] <- scale(Pred[,numcols])

# Load model results from scripts 03_
modAb <- readRDS(file = "Output/NEAbundance_OptimalModel.rds")
modDiv <- readRDS(file = "Output/NEDiversity_OptimalModel.rds")
modPred <- readRDS(file = "Output/PredRate_OptimalModel.rds")

modFullAb <- readRDS(file = "Output/NEAbundance_FullModel.rds")
modFullDiv <- readRDS(file = "Output/NEDiversity_FullModel.rds")
modFullPred <- readRDS(file = "Output/PredRate_FullModel.rds")

modFullAb_novg <- readRDS(file = "Output/NEAbundance_FullModel_novg.rds")
modFullDiv_novg <- readRDS(file = "Output/NEDiversity_FullModel_novg.rds")

modFullDiv_taxo <- readRDS(file = "Output/NEDiversity_FullModel_TaxoR.rds")
modOptDiv_taxo <- readRDS(file = "Output/NEDiversity_OptimalModel_TaxoR.rds")


## Make Table ---------------------------------------------------------------------------
Anova(modFullAb)
Anova(modFullDiv)
Anova(modFullPred)
Anova(modFullAb_novg)
Anova(modFullDiv_novg)
Anova(modFullDiv_taxo)
Anova(modOptDiv_taxo)

ColumnsT1 <- c("Variable", "Predictor", "Chi-sq", "df", "p-value")


## Store final table in Table folder----------------------------------------------------
write.csv(Table1, "Tables/Table1.csv")