## Clean the data, calculate diversity, abundances, predation rates
## store the outputs (final datasets ready for analyses)


## I. Main Analysis ----------------------------------------------------

## Functions
library(dplyr)
library(vegan)
library(tidyr)

## 1. Natural enemies-------------------------------------------------------

## Load data-------------------------------------
NE <- read.csv("Output/NatEnemies_raw.csv")

## Rename variables------------------------------
NE$Treatment <- factor(ifelse(NE$bande=="BE", "Low Div", "High Div"))
NE$Guild <- factor(ifelse(NE$piege == "B", "Vine",
                          ifelse(NE$piege == "F", "Vegetation",
                                 "Soil")))
NE$Distance <- ifelse(NE$mod == "00m", 0,
                      ifelse(NE$mod  == "05m", 15,
                             30))
NE$Site <- factor(NE$couple)


## Add landscape--------------------------------
Ldscp <- read.csv("Output/LandscapeVars.csv") # proportion of SNH in landscape
NE_lddf <- left_join(NE, Ldscp, by = "couple")

# landscape variable is proportion of SNH in 1000 m radius
NE$Ldscp <- NE_lddf$HSN1000

## CHECK Data for missing values, NAs and zero abundances---------------------------------
NE %>% group_by(Site,Treatment,Guild,Distance, session) %>% summarize(n()) # no. obs
NE %>% group_by(Site,Treatment,Guild,Distance, session) %>% summarize(sum(eff)) # sum abdc has NAs

## Nas for eff: site 9 and 10 only
filter(NE, is.na(eff)) %>% group_by(Site, Treatment, Guild, Distance, session) %>% summarize(n())

# Based on planning of sampling, should have NA for Vegetation in high div and low div treatments in site 9, 
# and in site 10, should be NA for all Vine and Vegetation samplings 
filter(NE, Site=="9") %>% group_by(Treatment, Guild, Distance, session) %>% summarize(mean(eff), min(eff), max(eff)) %>% View()

with(filter(NE, Site=="9"), table(Treatment, Guild, Distance, session))
# Missing rows: Vegetation at 15 and 30m for all sessions and all treatments
# Check in initial table (raw data) if values are present: no

## NAs were handled inconsistently in NE dataset: some are NAs and some are missing
## Remove all NAs from dataset
NE <- na.omit(NE)

## Missing values are not random: table shows that for Vegetation sampling 
# none of the sites have data for Distance 15 and 30 m
# vegetation was absent due to drought (except in strips) and sampling was impossible
NE %>% group_by(Guild, Distance) %>% summarize(n())

NE %>% group_by(Site, Treatment, Guild, Distance, session) %>% summarize(n())

## 1. Calculate abundance-------------------

## Total abundance per guild, distance and session----------------------------------------
Abundance <- NE %>% group_by(Site, Treatment, Guild, Distance, session) %>% 
  summarize(Total = sum(eff)) %>% ungroup()

## Missing values--------------------------
# we have 366 observations
nrow(Abundance)
# complete combination would be 648 observations
nrow(expand_grid(Site = c(1:12), 
                 Treatment = c("high", "low"), 
                 Guild=c("Soil", "Vegetation", "Vine"), 
                 Distance = c(0,15,30), 
                 session=c(1:3)))

## Save dataset
write.csv(Abundance, "Output/AbundanceClean.csv")


## 2. Calculate richness-----------
NE$session <- as.factor(NE$session)
# variables identifying the sample ID (to be used in dplyr:group_by)
WhichSampleID <- c("Site", "Treatment", "Guild", "Distance", "session")

# 2.1. Quantify data availability at different taxonomic resolution---------
# Makes table with the total number of individuals and those identified at each taxo resolution
TaxoReso <- c(
  TotalInd <-  sum(NE$eff),
  IndSpecies <-  sum(NE$eff[NE$Rspec=="x"]),
  IndGenus <-  sum(NE$eff[NE$Rgen=="x"&NE$Rspec==""]),
  IndFamily <-  sum(NE$eff[NE$Rgen=="" & NE$family != ""]),
  IndOrder <-  TotalInd - (IndSpecies+IndGenus+IndFamily))
TaxoReso

# Proportion of individuals id to the species vs genus level
TaxoReso[2]/TaxoReso[1] # 42% at species
(TaxoReso[2]+TaxoReso[3])/TaxoReso[1] # 77% at genus

## Taxonomic resolution available across all samples
TaxoResoPerSample <- NE %>% group_by_at(WhichSampleID) %>% 
  summarize(TotalInd = sum(eff),
            IndSpecies = sum(eff[Rspec=="x"]),
            IndGenus = sum(eff[Rgen=="x"]),
            PropInfoLostSpecies = (TotalInd - IndSpecies)/TotalInd,
            PropInfoLostGenus = (TotalInd - IndGenus)/TotalInd)

# 49 samples have no information at the species level
nrow(TaxoResoPerSample[TaxoResoPerSample$IndSpecies==0 & 
                         TaxoResoPerSample$TotalInd != 0,])
# 8 samples have no information at the genus level either
nrow(TaxoResoPerSample[TaxoResoPerSample$IndGenus==0 &
                         TaxoResoPerSample$TotalInd != 0,])
# But among those 8 samples, most have 1 or 2 individuals only
# which means, we can easily infer how many taxa are present: 
# 1 individual = 1 species; 2 individuals = 2 species if they belong to different families
TaxoResoPerSample$TotalInd[TaxoResoPerSample$IndGenus==0 & 
                    TaxoResoPerSample$IndSpecies == 0 &
                    TaxoResoPerSample$TotalInd != 0]

# How representative of each sample would taxonomic richness calculation be?
# On average, percentage of un-id individuals per sample
summary(TaxoResoPerSample$PropInfoLostSpecies)# on average our species richness estimates would miss 60% of individuals
hist(TaxoResoPerSample$PropInfoLostSpecies) # for most samples, we would miss 80% of individuals

# versus for genus richness, around 20% of individuals per sample are lost
summary(TaxoResoPerSample$PropInfoLostGenus)
hist(TaxoResoPerSample$PropInfoLostGenus)# and most samples miss less than 30% of individuals

## Data table NE contains information about families that are not represented in each sample
# i.e. 958 rows at the family level have counts = 0 (variable NE$eff)
length(NE$eff[NE$eff=="0"])


# 2.2. Compare richness based on different taxonomic scale------

## Taxonomic richness (across taxonomic resolutions) ----------------
# sum number of individuals per taxa (species, genus, family, order confounded)
TaxaList <- NE %>% group_by(Site, Treatment, Guild, Distance, session, taxon) %>% 
  summarize(total = sum(eff, na.rm = TRUE)) %>% ungroup()

# make wide table to calculate richness
TaxaList_wide <- as.data.frame(spread(TaxaList, taxon, total, fill = 0))

# remove categorical variables
JustTaxa <- subset(TaxaList_wide, 
                      select = -c(Site, Treatment, Guild, Distance, session))

TaxaList_wide$TR <- specnumber(JustTaxa)
summary(TaxaList_wide$TR) # between 0 and 20 taxa per sample (including various taxo scales)


## Species richness -------------
##Remove data for individuals that were not id to the species level
SpeciesList <- NE %>% 
  filter(Rspec =="x") %>% 
  group_by(Site, Treatment, Guild, Distance, session, taxon) %>% 
  summarize(total = sum(eff, na.rm = TRUE)) %>% ungroup()

# make wide table to calculate richness
SpeciesList_wide <- as.data.frame(spread(SpeciesList, taxon, total, fill = 0))

# remove categorical vars
JustSpecies <- subset(SpeciesList_wide, 
                       select = -c(Site, Treatment, Guild, Distance, session))
# richness calculation
SpeciesList_wide$SR <- specnumber(JustSpecies)


## Genus richness ------------------
# Only those individuals that were id to the genus
GenusList <- NE %>% 
  filter(Rgen == "x") %>% 
  group_by(Site, Treatment, Guild, Distance, session, Rgen_taxon) %>% 
  summarize(total = sum(eff, na.rm = TRUE)) %>% ungroup()

# Make wide table to calculate genus richness
GenusList_wide <- as.data.frame(spread(GenusList, Rgen_taxon, total, fill = 0))

# remove categorical vars
JustGenus <- subset(GenusList_wide, 
                      select = -c(Site, Treatment, Guild, Distance, session))
# richness calculation
GenusList_wide$GR <- specnumber(JustGenus)

## Tableau pour RDA avec landscape
Ldscp <- read.csv("Output/Landscapevars.csv")
Ldscp$Site <- as.factor(as.character(Ldscp$couple))
Ldscp$Ldscp <- Ldscp$HSN1000

TaxaList_Adrien <- dplyr::left_join(TaxaList_wide, Ldscp, by = "Site")
TaxaList_Adrien <- TaxaList_Adrien[TaxaList_Adrien$Guild != "Vegetation",]
write.csv(TaxaList_Adrien, "Output/TaxaList_forAdrien.csv")

## Make a diversity table with different richness estimates-----

# Create diversity dataframe, and adds taxonomic richness
DivDf <- subset(TaxaList_wide, 
                select = c(Site, Treatment, Guild, Distance, session, TR))

# add species richness when possible
SpDf <- subset(SpeciesList_wide, 
               select = c(Site, Treatment, Guild, Distance, session, SR))

DivDf <- left_join(DivDf, SpDf, by = WhichSampleID)

## add genus richness when possible
GnDf <- subset(GenusList_wide, 
               select = c(Site, Treatment, Guild, Distance, session, GR))

DivDf <- left_join(DivDf, GnDf, by = WhichSampleID)

## add abundance data
Diversity <- left_join(DivDf, Abundance, by = WhichSampleID)
head(Diversity)
summary(Diversity) # many NAs

## When abundance was 0 and 1, richness should be 0 and 1 too
Div_clean <- Diversity
Div_clean$SR[Div_clean$Total <= 1] <-  Div_clean$GR[Div_clean$Total <= 1] <- 
  Div_clean$TR[Div_clean$Total <= 1]
  
summary(Div_clean) # less NAs

## What about remaining NAs for genus richness? can we make educated guess about how many genus?
Div_clean[is.na(Div_clean$GR),]
# look at the raw data for those NAs,
left_join(Div_clean[is.na(Div_clean$GR),], NE, by = WhichSampleID) %>% filter(eff != 0)

# At Site 1, 2 individuals, belonging to different families: GR = TR = 2
Div_clean$GR[is.na(Div_clean$GR) & Div_clean$Site==1] <- Div_clean$TR[is.na(Div_clean$GR) & Div_clean$Site==1]

# At Site 6, there were 3 individuals, belonging to 2 distinct families
# At Site 11, there were 4 individuals, belong to 2 distinct families as well
## for those, we could assume at least 2 species are present, but then we would be
# treating data differently (special case for the NAs versus the rest of the data)

## Compare richness estimates----
# Accumulation curves
par(mfrow=c(1,3))
plot(Div_clean$Total, Div_clean$TR, main = "taxonomic")
plot(Div_clean$Total, Div_clean$GR, main = "genus")
plot(Div_clean$Total, Div_clean$SR, main = "species")
par(mfrow=c(1,1))
# Genus richness vs. other richness
par(mfrow=c(1,2))
plot(Div_clean$TR, Div_clean$GR, main = "genus vs. taxonomic")
abline(0,1)
plot(Div_clean$SR, Div_clean$GR, main = "genus vs. species")
abline(0,1)
par(mfrow=c(1,1))


# 2.3. Rarefy richness based on total abundance------------
# rarefied richness: some samples have zero abundances= not possible to rarefy
(raremaxT <- min(rowSums(JustTaxa)))
(raremaxS <- min(rowSums(JustSpecies)))
(raremaxG <- min(rowSums(JustGenus)))


# 2.4. Store richness data table------------------
## Rename cols
Div_final <- subset(Div_clean, select = -c(SR, Total)) %>% 
  rename(TaxaR = TR, GenusR = GR)

## Save dataset
write.csv(Div_final, "Output/DiversityClean.csv")




# 2. Predation of Lobesia eggs-------------------------------------------------------

## Load data--------------------------------------------
Pred <- read.csv("Output/PredationPest_raw.csv")

## Rename variables------------------------------
Pred$Treatment <- factor(ifelse(Pred$bande=="BE", "Low Div", "High Div"))

Pred$Distance <- ifelse(Pred$mod == "00m", 0,
                      ifelse(Pred$mod  == "05m", 15,
                             30))
Pred$Site <- factor(Pred$couple)
Pred$Session <- factor(Pred$session)

## Check missing values------------------------

table(Pred$Site, Pred$Treatment) # on site 6, there are missing values
table(Pred$Site, Pred$Session) # on site 6, half the data for session 1 are missing
table(Pred$Site, Pred$Session) # session 3 as only 18 observations for all sites
table(Pred$Site, Pred$Distance) # All distances were sampled equally, even site 6


## Predation rates------------------------------

Pred$PredRate <- Pred$Tx_pred # proportion of eggs predated per card

# Subset dataframe-------------------------------

PredClean <- subset(Pred, select = c("Site", "Treatment", "Distance", "Session", "PredRate"))

# Save cleaned data---------------------------------
write.csv(PredClean, "Output/PredationPest_clean.csv")

