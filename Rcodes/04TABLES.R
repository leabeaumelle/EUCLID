## Create table to show model results



## Functions----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)

## Data----------------------------------------------------------------------------
# Load model results from scripts 03_


## Make Table ---------------------------------------------------------------------------

ColumnsT1 <- c("Variable", "Predictor", "Chi-sq", "df", "p-value")


## Store final table in Table folder----------------------------------------------------
write.csv(Table1, "Tables/Table1.csv")