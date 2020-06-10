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


## 2. Calculate rarefied richness-----------

# 2.1. Make a species matrix-------------------
# create unique identifier for each sample
NE$SampleID <- factor(paste0(NE$Site, NE$Treatment, NE$Guild, NE$Distance, NE$session))

## PREVIOUS CODE FROM ARTHUR:
### SPECIFIC RICHNESS DATAFRAME CREATION 
WhichVars <- c("Site","Treatment","Distance","Guild","session",
                    "SampleID","order","family","genus","species",
                    "eff","taxon","codeRs")

# A final df where data for all samples will be put
NE_R_fin <- data.frame(matrix(ncol=length(WhichVars)))
colnames(NE_R_fin) <- WhichVars

# An intermediate df to store abundance data for each sample in the following loop
NE_R <- data.frame(matrix(ncol=length(WhichVars)))


for (j in levels(NE$SampleID)) {  # For each sample (SampleID)
  
  # Subset count data
  NE_R <- NE[NE$SampleID==j, WhichVars] ;  rownames(NE_R) <- NULL
  
  # SPECIES : creates a df to store count data for individuals identified at species level
  NE_RS <- data.frame(matrix(ncol=length(WhichVars))) ; colnames(NE_RS) <- colnames(NE_R)
  
  for (i in levels(NE_R$codeRs)){ # For each taxon (codeRs)
    
    # add in the dataframe all individuals identified at the species level.
    if(grepl("S_", i) == TRUE) ## Here I changed Arthur code from str_detect to grepl
      NE_RS <- rbind(NE_RS , NE_R[NE_R$codeRs ==i,])
    
  }
  NE_RS <- NE_RS[-1,] ; rownames(NE_RS) <- NULL
  
  
  # GENUS : creates a df to store count data for individuals indentified at genus level
  NE_RG <- data.frame(matrix(ncol=length(WhichVars))) ; colnames(NE_RG) <- colnames(NE_R)
  for (i in 1:nrow(NE_R)){
    
    # add in the dataframe all individuals identified at the genus level
    if(grepl("G_", NE_R[i, "codeRs"]) == TRUE & 
       
       # AND whose genus isn't present in the Species-level dataframe (EN_Rs).
       as.character(NE_R[i,"genus"]) %in% NE_RS$genus == FALSE)
      
      NE_RG <- rbind(NE_RG , NE_R[i,])
    
  }
  NE_RG <- NE_RG[-1,] ; rownames(NE_RG) <- NULL
  
  
  # FAMILY
  NE_RF <- data.frame(matrix(ncol=length(WhichVars))) ; colnames(NE_RF) <- colnames(NE_R)
  for (i in 1:nrow(NE_R)){
    
    # add in the dataframe all individuals identified at the family level...
    if(grepl("F_", NE_R[i, "codeRs"]) == TRUE & 
       
       # AND whose family isn't present in the Species-level dataframe (EN_Rs)...
       as.character(NE_R[i,"family"]) %in% NE_RS$family == FALSE &
       
       # AND whose family isn't present in the Genus-level dataframe (EN_Rg).
       as.character(NE_R[i,"family"]) %in% NE_RG$family == FALSE)
      
      NE_RF <- rbind(NE_RF , NE_R[i,])
    
  }
  NE_RF <- NE_RF[-1,] ; rownames(NE_RF) <- NULL
  
  # ORDER : this part of the code looks strange to me, there are not so many orders, why 
  # the loop is excluding orders that are already present in previous df? There are only 5 orders
  NE_RO <- data.frame(matrix(ncol=length(WhichVars))) ; colnames(NE_RO) <- colnames(NE_R)
  for (i in 1:nrow(NE_R)){
    
    # ajoute dans le tableau les individus identifies a l'ordre ...
    if(grepl("O_", NE_R[i, "codeRs"]) == TRUE & 
       
       # AND whose order isn't present in the Species-level dataframe (EN_Rs)...
       as.character(NE_R[i,"order"]) %in% NE_RS$order == FALSE &
       
       # AND whose order isn't present in the Genus-level dataframe (EN_Rg)...
       as.character(NE_R[i,"order"]) %in% NE_RG$order == FALSE &
       
       # AND whose order isn't present in the Family-level dataframe (EN_Rf)...
       as.character(NE_R[i,"order"]) %in% NE_RF$order == FALSE)
      
      NE_RO <- rbind(NE_RO , NE_R[i,])
    
  }
  NE_RO <- NE_RO[-1,] ; rownames(NE_RO) <- NULL
  
  
  # The 4 dataframes at each taxonomic levels are grouped in the EN_R_fin dataframe.
  NE_R_fin <- rbind(NE_R_fin, NE_RS , NE_RG , NE_RF , NE_RO)
  
  # EN_R_fin est donc le tableau qui permettra de calculer les richesses
  # I don't understand how : this table is just reordered according to taxonomic level??
}
rm(NE_RS) ; rm(NE_RG) ; rm(NE_RF) ; rm(NE_RO)
NE_R_fin <- NE_R_fin[-1,]
NE_R_fin$taxon <- as.factor(NE_R_fin$taxon)
NE_R_fin$SampleID <- as.factor(NE_R_fin$SampleID)



## EN RICHNESS (at a taxonomic level)
datalong <- NE_R_fin %>% group_by(SampleID, taxon) %>% 
  summarize(total = sum(eff, na.rm = TRUE)) %>% ungroup()

EUCLID_Rstot <- as.data.frame(spread(datalong, taxon, total, fill = 0))

# Now checking if NE_R_fin dataset is necessary
test1 <- NE %>% group_by(SampleID, taxon) %>% 
  summarize(total = sum(eff, na.rm = TRUE)) %>% ungroup()

data_wide<-as.data.frame(spread(test1, taxon, total, fill = 0))

## VALUES ARE THE SAME USING THE NE_R_FIN OR NE DATASETS
all.equal(data_wide[,2], EUCLID_Rstot[,2]) 
all.equal(specnumber(EUCLID_Rstot[,-1]), specnumber(data_wide[,-1]), check.names = FALSE )



## 3. Add explanatory variables and remove unnecessary columns----------------------



## Save cleaned data----------------------
write.csv(NE, "Output/NatEnemiesDiv_clean.csv")

# 2. Predation of Lobesia eggs-------------------------------------------------------

## Load data--------------------------------------------
Pred <- read.csv("Output/PredationPest_raw.csv")

## 1. Calculate predation rates---------------------------------

## 2. Add explanatory variables---------------------------------



# Save cleaned data---------------------------------
write.csv("Output/PredationPest_clean.csv")


## Supplementary Analysis - Year effect