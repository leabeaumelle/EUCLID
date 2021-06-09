## Loads the data 

# The script loads the raw data from EUCLID project, creates metadata list for natural enemies data
# raw data processed are then stored into the Output folder for further analyses

# Script is organized according to the paper: with a main analysis (I) on natural enemies and predation rates
# and supplementary analyses (II) for testing the year effect, and predation of larvae of Lobesia to complement
# the eggs predation data.

## Functions and packages---------------------------------------------------------------
library(dplyr)


## I. Main analysis on 2018 dataset---------------------------------------------------------

## 1. Natural enemies--------------------------------------------------------------------

# TOTAL INVENTORY of Natural enemies (EN) 2018
EN_2018<-read.table("Data/EUCLID - inventaire 2018.txt",sep="\t",dec=",",header=T)

# Save raw data in the output folder
write.csv(EN_2018, "Output/NatEnemies_raw.csv")

# save metadata
MetaDataEN <- list(
  title="Metadata for natural enemies raw dataset",
  
  README="Arthropods were sampled three times (sessions) in June, August, September, using three sampling methods
at 9 vineyard sites were two treatments (flower strips vs. grassy strips) were applied to separate plots.
Natural enemies belonging to 4 orders (Araneae, Neuroptera, Dermaptera, Opiliones) were counted
and identified to the lowest possible taxonomic resolution. Abundances are sums of number of individuals
collected from several sampling units (pitfall traps) per sampling method per session (HOW MANY?). Each row is a taxonomic group,
when there were e.g. no Araneae, a row for Araneae is created with abundance = 0. When there were several species of
  Araneae, several rows for each species and their respective abundance is created. There are NAs when the sampling
could not be done (especially for sweaping (fauchage) and beating (battage) when no vegetation was present 
(due to drought and/or farmer had to remove vegetation. ",
  
  data.frame(columns=names(EN_2018), 
             values=c("name of the site", 
                      "treatment modality (grassy or flower strip)",
                      "distance (m) to the right border of the strip where sampling conducted",
                      "sampling method (pitfall, sweaping or beating)",
                      "",
                      "french date of sampling",
                      "","","","",
                      "number of individuals for the taxonomic group considered, total abundance across several sampling units that were pooled per session",
                      "sampling session (three sessions were conducted)",
                      "x marked rows are relevant to compute species richness: individuals are id to the species",
                      "x marked rows are relevant to compute genera richess: individuals are id to the genera",
                      "full name of the taxonomic group for species richness calculations (family_genera_species)",
                      "full name of the taxonomic group for genera richness calculations (family_genera)",
                      "", "", ""))
)
# save metadata with raw data
cat(capture.output(print(MetaDataEN), file="Data/MetaData_NatEnemies.txt"))


## Create Raw data tables for natural enemies-----------------------------------
# As factor
EN_2018$couple<-as.factor(EN_2018$couple)     ;     EN_2018$mod<-as.factor(EN_2018$mod)
EN_2018$session<-as.factor(EN_2018$session)   ;     EN_2018$Rspec<-as.factor(EN_2018$Rspec)
EN_2018$Rgen<-as.factor(EN_2018$Rgen)
# Add "genus sp" column : concatenation of "genus" and "species.
EN_2018[,"genus_sp"] <- apply(EN_2018[,c("genus","species")] , 1 , paste , collapse=" ")
EN_2018$genus_sp<-as.factor(EN_2018$genus_sp)
# Add a code2_session column
EN_2018 <- EN_2018 %>% mutate(code2_session = paste(code2,session))
# "taxon" column creation : order_family_genus_species
EN_2018 <- EN_2018 %>% mutate(taxon = paste(order,family,genus,species,sep="_"))
# Column "codeRs" : if a S_ is added before the name, it means that the individual is identified at a species level,
# a G_ at the genus level, a F_ at the family level and a 0_ at the order level.
EN_2018 <- EN_2018 %>% mutate(codeRs = case_when(
  str_detect(taxon , "___") == TRUE ~ str_c("O_" , taxon), # Order level
  str_detect(taxon , "idae__") == TRUE ~ str_c("F_" , taxon), # Family level
  species == "sp." & nchar(as.character(species)) <= 3 ~ str_c("G_" , taxon), # Genus level
  genus != "" & nchar(as.character(species)) > 3 ~ str_c("S_" , taxon) # Species level
))


# save raw data in output
write.csv(EN_2018, "Output/NatEnemies_raw.csv")

## 2. Predation rates---------------------------------------------------------
# EGG CARDS IMPORTATION (data 2017 & 2018)
CO<-read.table("Data/EUCLID - CO.txt",sep="\t",dec=".",header=T,na.strings="NA")
CO$couple<-as.factor(CO$couple)
CO$mod<-as.factor(CO$mod)
CO$session<-as.factor(CO$session)
CO$annee<-as.factor(CO$annee)

# code2-session column
CO <- CO %>% mutate(code2_session = paste(code2,session))

# Keep only 2018 data
CO2018 <- subset(CO, annee == "2018")

# save raw data in output
write.csv(CO2018, "Output/PredationPest_raw.csv")


## 3. Landscape variables----------------------------------------------------------------
# IMPORTATION DATA SNH (Semi-natural-habitats)
HSN<-read.table("Data/EUCLID - HSN.txt",sep="\t",dec=".",header=T,na.strings="NA")
HSN$couple<-as.factor(HSN$couple)

# save raw data in output
write.csv(HSN, "Output/LandscapeVars.csv")

