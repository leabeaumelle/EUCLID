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

