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

png("Figures/FigSupp_DistanceAbd.png",
    width=w/3,
    height=w/3,
    units = "cm",
    res=ppi)

FigS1
dev.off()





## Alternative figure 1 showing different guilds
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

F1A <- ggplot(data = Abundance, aes(y=Total, x=Guild, fill = Treatment))+
  geom_boxplot()+
  stat_summary(fun="mean", position = position_dodge(width = 0.75), geom = "point", colour = "red")+
  geom_point(position = position_jitterdodge(dodge.width = 0.75), alpha = .4, size = 0.5)+
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
Abundance %>% filter(!is.na(Total)) %>% group_by(Treatment, Guild) %>% summarize(n())

# Panel B - Flower strips don't change the diversity of natural enemies (genus richness)

# remove vegetation guild abundances
Diversity <- Diversity[Diversity$Guild != "Vegetation",]

# reorder levels of treatment
Diversity$Treatment <- relevel(Diversity$Treatment, ref = "Low Div")

# set sizes of text in plots
sizetext <- 10
sizelegend <- 8

F2A <- ggplot(data = Diversity[!is.na(Diversity$TaxaR),], aes(y=TaxaR, x=Guild, fill = Treatment))+
  geom_boxplot()+
  stat_summary(fun="mean", position = position_dodge(width = 0.75), colour = "red")+
  geom_point(position = position_jitterdodge(jitter.height = .05), alpha = .4, size = 0.5)+
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
Diversity %>% filter(!is.na(TaxaR)) %>% group_by(Treatment, Guild) %>% summarize(n())


# Panel C - Flower strips effect on predation depends on the landscape

# reorder levels of treatment
Pred$Treatment <- relevel(Pred$Treatment, ref = "Low Div")

# set sizes of text in plots
sizetext <- 10
sizelegend <- 8

F2C <- ggplot(data = Pred[!is.na(Pred$PredRate),], aes(y=PredRate, x=Treatment, fill = Treatment))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0.95, jitter.height = 0.01), alpha = .4, size = 0.5)+
  stat_summary(fun="mean", position = position_dodge(width = 0.75), colour = "red", size = 0.7, shape = 16)+
  scale_fill_viridis_d(begin = 0.5)+
  ylab("Predation Rates")+
  theme_bw()+
  theme(legend.position = "right", 
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

png("Figures/Fig1_guilds.png",
    width=w,
    height=w/3,
    units = "cm",
    res=ppi)

F1A+F2A+F2C

dev.off()


## Ladnscape gradient across sites
Plot_df <- Ldscp[, names(Ldscp) %in% c("couple", "HSN1000")] %>% 
  rename("Lscp" = HSN1000, Sites = "couple")
Plot_df$Site <- as.factor(Plot_df$Site)

Plot_df2 <- droplevels(semi_join(Plot_df, Abundance, by = c("Site")))

Plot_df2$Sites = with(Plot_df2, reorder(Sites, Lscp, median))
p<-Plot_df2 %>% 
  ggplot(aes(x=Sites, y=Lscp, fill = Sites))+
  geom_point(size = 3)+
  geom_step()+
  ylab(expression(atop("Landscape complexity", paste("% semi-natural habitats"))))+
  ylim(c(0,100))+
  
  theme_bw()+
  theme(legend.position = "none", 
        axis.text.x = element_blank())

# save a png with high res
ppi <- 250# final: 600 # resolution
w <- 10 # width in cm

png("Figures/LandscapeGradient.png",
    width=w,
    height=w*0.8,
    units = "cm",
    res=ppi)

p
dev.off()
