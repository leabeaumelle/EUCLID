## Create figures for the manuscript

## The script creates and saves three main figures

## Functions----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(sjPlot)
library(ggeffects)

## Data----------------------------------------------------------------------------
Abundance <- read.csv("Output/AbundanceClean.csv")
Diversity <- read.csv("Output/DiversityClean.csv")
Pred <- read.csv("Output/PredationPest_clean.csv")


## Figure 1 - overall effect of flower strips on natural enemies and predation------------
# three panels with overall flower strip effect across landscapes and communitiy types

# Panel A - Flower strips increase the abundance of natural enemies

# remove vegetation guild abundances
Abundance <- Abundance[Abundance$Guild != "Vegetation",]

# reorder levels of treatment
Abundance$Treatment <- relevel(Abundance$Treatment, ref = "Low Div")

# set sizes of text in plots
sizetext <- 10
sizelegend <- 8

F1A <- ggplot(data = Abundance, aes(y=Total, x=Treatment, fill = Treatment))+
  geom_boxplot()+
  geom_jitter(width = .25, alpha = .4, size = 0.4)+
  scale_fill_viridis_d(begin = 0.5)+
  # scale_fill_manual(values = c(pal[2], pal[3]))+
  ylab("Natural enemies abundance")+
  theme_bw()+
  theme(legend.position = "none", 
        axis.text.y=element_text(face = "bold", size = sizelegend),
        axis.text.x=element_text(size = sizetext),
        axis.title.y = element_text(size=sizetext, face = "bold"),
        axis.title.x = element_blank())
F1A

# Sample sizes
Abundance %>% filter(!is.na(Total)) %>% group_by(Treatment) %>% summarize(n())

# Panel B - Flower strips don't change the diversity of natural enemies (genus richness)

# remove vegetation guild abundances
Diversity <- Diversity[Diversity$Guild != "Vegetation",]

# reorder levels of treatment
Diversity$Treatment <- relevel(Diversity$Treatment, ref = "Low Div")

# set sizes of text in plots
sizetext <- 10
sizelegend <- 8

F2A <- ggplot(data = Diversity[!is.na(Diversity$TaxaR),], aes(y=TaxaR, x=Treatment, fill = Treatment))+
  geom_boxplot()+
  geom_jitter(width = .25, alpha = .4, size = 0.4)+
  scale_fill_viridis_d(begin = 0.5)+
  # scale_fill_manual(values = c(pal[2], pal[3]))+
  ylab("Natural enemies taxa richness")+
  theme_bw()+
  theme(legend.position = "none", 
        axis.text.y=element_text(face = "bold", size = sizelegend),
        axis.text.x=element_text(size = sizetext),
        axis.title.y = element_text(size=sizetext, face = "bold"),
        axis.title.x = element_blank())
F2A

# Sample sizes
Diversity %>% filter(!is.na(TaxaR)) %>% group_by(Treatment) %>% summarize(n())


# Panel C - Flower strips effect on predation depends on the landscape

# reorder levels of treatment
Pred$Treatment <- relevel(Pred$Treatment, ref = "Low Div")

# set sizes of text in plots
sizetext <- 10
sizelegend <- 8

F2C <- ggplot(data = Pred[!is.na(Pred$PredRate),], aes(y=PredRate, x=Treatment, fill = Treatment))+
  geom_boxplot()+
  geom_jitter(width = .25, alpha = .4, size = 0.4)+
  scale_fill_viridis_d(begin = 0.5)+
  # scale_fill_manual(values = c(pal[2], pal[3]))+
  ylab("Predation Rates")+
  theme_bw()+
  theme(legend.position = "none", 
        axis.text.y=element_text(face = "bold", size = sizelegend),
        axis.text.x=element_text(size = sizetext),
        axis.title.y = element_text(size=sizetext, face = "bold"),
        axis.title.x = element_blank())
F2C

# Sample sizes
Pred %>% filter(!is.na(PredRate)) %>% group_by(Treatment) %>% summarize(n())

# save a png with high res
ppi <- 300# final: 600 # resolution
w <- 20 # width in cm

png("Figures/Fig1.png",
    width=w,
    height=w/3,
    units = "cm",
    res=ppi)

F1A+F2A+F2C

dev.off()


## Figure 2 - landscape effect----------------------------------------------------
# figure showing how landscape modulates flower strip effects

## Data----------------------------------------------------------------------------
Abundance <- read.csv("Output/AbundanceClean.csv")
Diversity <- read.csv("Output/DiversityClean.csv")
Pred <- read.csv("Output/PredationPest_clean.csv")
Ldscp <- read.csv("Output/Landscapevars.csv")
Ldscp$Site <- as.factor(as.character(Ldscp$couple))
Ldscp$Ldscp <- Ldscp$HSN1000


Abundance$Site <- as.factor(Abundance$Site)
Abundance <- left_join(Abundance, Ldscp, by = "Site")
Abundance$Site <- as.factor(Abundance$Site)
Abundance$session <- as.factor(as.character(Abundance$session))

Diversity$Site <- as.factor(Diversity$Site)
Diversity <- left_join(Diversity, Ldscp, by = "Site")
Diversity$Site <- as.factor(Diversity$Site)
Diversity$session <- as.factor(as.character(Diversity$session))

Pred$Site <- as.factor(as.character(Pred$Site))
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
Pred_sc$PredatedEggs <- Pred_sc$PredRate*10

# Remove vegetation guild
Abs <- droplevels(Abs[Abs$Guild != "Vegetation",])
Div <- droplevels(Div[Div$Guild != "Vegetation",])

# Load model results from scripts 03_
modAb <- readRDS(file = "Output/NEAbundance_OptimalModel.rds")
modDiv <- readRDS(file = "Output/NEDiversity_OptimalModel.rds")
modPred <- readRDS(file = "Output/PredRate_OptimalModel.rds")

modFullAb <- readRDS(file = "Output/NEAbundance_FullModel.rds")
modFullDiv <- readRDS(file = "Output/NEDiversity_FullModel.rds")
modFullPred <- readRDS(file = "Output/PredRate_FullModel.rds")

# Make plots: prototype
plot_model(modFullAb, type = "pred", terms = c("Ldscp [all]","Treatment"))+
plot_model(modFullDiv, type = "pred", terms = c("Ldscp [all]","Treatment"))+
plot_model(modFullPred, type = "pred", terms = c("Ldscp [all]","Treatment"))


# colors: get them from viridis
mycols <- c(scales::viridis_pal(begin = 0.5)(2)[2], 
            scales::viridis_pal(begin = 0.5)(2)[1])

# set sizes of text in plots
sizetext <- 10
sizelegend <- 9
#colours
mycols <- c("#F2DA02",scales::viridis_pal(begin = 0.5)(2)[1])


# PANEL ABUNDANCE

# get the predictions with ggeffects
me <- ggemmeans(modFullAb, c("Ldscp [all]", "Treatment"))
# take the raw data from ggeffects and compute means and sds
raw <- attr(me, "rawdata")
raw2 <- raw %>% group_by(x, group) %>% summarise(means = mean(response, na.rm = TRUE), 
                                                 sdev = sd(response, na.rm = TRUE))



Fig2A <- ggplot(me, aes(x = x, y = predicted, colour = group), col = mycols) +
  geom_line(size = 1)+
  geom_errorbar(data = raw2, mapping = aes(x = x, y = means, ymin=means-sdev, ymax = means+sdev, group = group), colour = "black", 
                width = 0.1, position = position_dodge(width = 0.1), size = 0.1)+
  geom_point(data = raw2, mapping = aes(x = x, y = means, fill = group, shape = group), 
             size = 2.2, colour = "black", position = position_dodge(width = 0.1))+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_shape_manual(values = c(21,24))+
  geom_ribbon(inherit.aes = FALSE, 
              mapping = aes(x = x, y = predicted, group = group, fill = group,
                            ymin = conf.low, ymax = conf.high), alpha = 0.15)+
  # using jitter for observations here
  # scale_fill_manual(values=mycols)+
  ylab("Abundance (individuals)")+
  scale_x_continuous(breaks = c((30-mean(Abs$HSN1000))/sd(Abs$HSN1000), 
                                (40-mean(Abs$HSN1000))/sd(Abs$HSN1000), 
                                (50-mean(Abs$HSN1000))/sd(Abs$HSN1000), 
                                ((60-mean(Abs$HSN1000))/sd(Abs$HSN1000))),
                     labels = c(30, 40, 50, 60)
  )+
  xlab("Landscape complexity (%)")+
  ggtitle("")+
  theme_bw()+
  theme(legend.position = "none",
        legend.text = element_text(size = sizelegend),
        legend.title = element_blank(),
        axis.text.y=element_text(face = "bold", size = sizelegend),
        axis.text.x=element_text(face = "bold", size = sizelegend),
        axis.title.y = element_text(size=sizetext, face = "bold"),
        axis.title.x = element_text(size=sizetext, face = "bold"))

# PANEL RICHNESS

# get the predictions with ggeffects
me <- ggemmeans(modFullDiv, c("Ldscp [all]", "Treatment"))
# take the raw data from ggeffects and compute means and sds
raw <- attr(me, "rawdata")
raw2 <- raw %>% group_by(x, group) %>% summarise(means = mean(response, na.rm = TRUE), 
                                                 sdev = sd(response, na.rm = TRUE))

Fig2B <- ggplot(me, aes(x = x, y = predicted, colour = group), col = mycols) +
  geom_line(size = 1, linetype = 2)+
  geom_errorbar(data = raw2, mapping = aes(x = x, y = means, ymin=means-sdev, ymax = means+sdev, group = group), colour = "black", 
                width = 0.1, position = position_dodge(width = 0.1), size = 0.1)+
  geom_point(data = raw2, mapping = aes(x = x, y = means, fill = group, shape = group), 
             size = 2.2, colour = "black", position = position_dodge(width = 0.1))+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_shape_manual(values = c(21,24))+
  geom_ribbon(inherit.aes = FALSE, 
              mapping = aes(x = x, y = predicted, group = group, fill = group,
                            ymin = conf.low, ymax = conf.high), alpha = 0.15)+
  ylim(c(0,16))+
  ylab("Taxonomic richness (taxa)")+
  scale_x_continuous(breaks = c((30-mean(Abs$HSN1000))/sd(Abs$HSN1000), 
                                (40-mean(Abs$HSN1000))/sd(Abs$HSN1000), 
                                (50-mean(Abs$HSN1000))/sd(Abs$HSN1000), 
                                ((60-mean(Abs$HSN1000))/sd(Abs$HSN1000))),
                     labels = c(30, 40, 50, 60)
  )+
  xlab("Landscape complexity (%)")+
  ggtitle("")+
  theme_bw()+
  theme(legend.position = "none",
        legend.text = element_text(size = sizelegend),
        legend.title = element_blank(),
        axis.text.y=element_text(face = "bold", size = sizelegend),
        axis.text.x=element_text(face = "bold", size = sizelegend),
        axis.title.y = element_text(size=sizetext, face = "bold"),
        axis.title.x = element_text(size=sizetext, face = "bold"))

# PANEL PREDATION

# get the predictions with ggeffects
me <- ggemmeans(modFullPred, c("Ldscp [all]", "Treatment"))
# take the raw data from ggeffects and compute means and sds
raw <- attr(me, "rawdata")
raw2 <- raw %>% group_by(x, group) %>% summarise(means = mean(response, na.rm = TRUE), 
                                                 sdev = sd(response, na.rm = TRUE))

Fig2C <- ggplot(me, aes(x = x, y = predicted, colour = group), col = mycols) +
  geom_line(size = 1)+
  geom_errorbar(data = raw2, mapping = aes(x = x, y = means, ymin=means-sdev, ymax = means+sdev, group = group), colour = "black", 
                width = 0.1, position = position_dodge(width = 0.1), size = 0.1)+
  geom_point(data = raw2, mapping = aes(x = x, y = means, fill = group, shape = group), 
             size = 2.2, colour = "black", position = position_dodge(width = 0.1))+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_shape_manual(values = c(21,24))+
  geom_ribbon(inherit.aes = FALSE, 
              mapping = aes(x = x, y = predicted, group = group, fill = group,
                            ymin = conf.low, ymax = conf.high), alpha = 0.15)+
  # using jitter for observations here
  # scale_fill_manual(values=mycols)+
  ylab("Predation (eggs predated)")+
  scale_x_continuous(breaks = c((30-mean(Abs$HSN1000))/sd(Abs$HSN1000), 
                                (40-mean(Abs$HSN1000))/sd(Abs$HSN1000), 
                                (50-mean(Abs$HSN1000))/sd(Abs$HSN1000), 
                                ((60-mean(Abs$HSN1000))/sd(Abs$HSN1000))),
                     labels = c(30, 40, 50, 60)
  )+
  xlab("Landscape complexity (%)")+
  ggtitle("")+
  theme_bw()+
  theme(legend.position = "right",
        legend.text = element_text(size = sizelegend),
        legend.title = element_blank(),
        axis.text.y=element_text(face = "bold", size = sizelegend),
        axis.text.x=element_text(face = "bold", size = sizelegend),
        axis.title.y = element_text(size=sizetext, face = "bold"),
        axis.title.x = element_text(size=sizetext, face = "bold"))


# save a png with high res
ppi <- 300# final: 600 # resolution
w <- 25 # width in cm

png("Figures/Fig2.png",
    width=w,
    height=w/2.5,
    units = "cm",
    res=ppi)

Fig2A+labs(tag = "A")+
  Fig2B+labs(tag = "B")+
  Fig2C+labs(tag = "C")
dev.off()

# Get sample sizes for the caption
Abs %>% group_by(Treatment, HSN1000) %>% summarize(n())
Div %>% group_by(Treatment, HSN1000) %>% summarize(n())
Pred_sc %>% group_by(Treatment, HSN1000) %>% summarize(n())



## Figure 2 - response of different guilds at different distances of flower strips
# three to nine panels showing how response depends on community type and distance to strip
