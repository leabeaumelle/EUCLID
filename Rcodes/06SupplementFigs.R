## Creates supplementary figures

## Functions----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(sjPlot)

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


## Plot the effect of distance
# colors: get them from viridis
mycols <- c(scales::viridis_pal(begin = 0.5)(2)[2], 
            scales::viridis_pal(begin = 0.5)(2)[1])
mycols <- c("#F2DA02",scales::viridis_pal(begin = 0.5)(2)[1])

# set sizes of text in plots
sizetext <- 8
sizelegend <- 7

# Make plot
FigS1 <- plot_model(modFullAb, type = "pred", terms = c("Distance"),
                    colors = mycols, dot.size = 1.5, line.size = 1, 
                    show.data = TRUE)+
    scale_x_continuous(breaks = c(-1.0319298, 0.1670744, 1.3660785),
                       labels = c(0, 15, 30)
                      )+
  ylab("Abundance (individuals)")+
  xlab("Distance to the center of the strip (m)")+
  ggtitle("")+
  theme_bw()+
  theme(legend.position = "none",
        legend.text = element_text(size = sizelegend),
        legend.title = element_text(size = sizetext),
        axis.text.y=element_text(face = "bold", size = sizelegend),
        axis.text.x=element_text(face = "bold", size = sizelegend),
        axis.title.y = element_text(size=sizetext, face = "bold"),
        axis.title.x = element_text(size=sizetext, face = "bold"))
FigS1

# save a png with high res
ppi <- 300# final: 600 # resolution
w <- 20 # width in cm

png("Figures/FigSupp_DistanceAbd",
    width=w,
    height=w/3,
    units = "cm",
    res=ppi)

FigS1
dev.off()
