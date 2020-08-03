## Create figures for the manuscript

## The script creates and saves three main figures

## Functions----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)

## Data----------------------------------------------------------------------------
Abundance <- read.csv("Output/AbundanceClean.csv")
Diversity <- read.csv("Output/DiversityClean.csv")
Pred <- read.csv("Output/PredationPest_clean.csv")


## Figure 1 - overall effect of flower strips on natural enemies and predation------------
# three panels with overall flower strip effect across landscapes and communitiy types

# Panel A - Flower strips increase the abundance of natural enemies

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

# reorder levels of treatment
Diversity$Treatment <- relevel(Diversity$Treatment, ref = "Low Div")

# set sizes of text in plots
sizetext <- 10
sizelegend <- 8

F2A <- ggplot(data = Diversity[!is.na(Diversity$GenusR),], aes(y=GenusR, x=Treatment, fill = Treatment))+
  geom_boxplot()+
  geom_jitter(width = .25, alpha = .4, size = 0.4)+
  scale_fill_viridis_d(begin = 0.5)+
  # scale_fill_manual(values = c(pal[2], pal[3]))+
  ylab("Natural enemies genus richness")+
  theme_bw()+
  theme(legend.position = "none", 
        axis.text.y=element_text(face = "bold", size = sizelegend),
        axis.text.x=element_text(size = sizetext),
        axis.title.y = element_text(size=sizetext, face = "bold"),
        axis.title.x = element_blank())
F2A

# Sample sizes
Diversity %>% filter(!is.na(GenusR)) %>% group_by(Treatment) %>% summarize(n())


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

## Figure 2 - response of different guilds at different distances of flower strips
# three to nine panels showing how response depends on community type and distance to strip


## Figure 3 - landscape effect
# figure showing for one selected outcome how landscape modulates flower strip effects
