## Load functions and packages for the analysis

library(tidyverse)
library(lme4)
library(MuMIn)
library(DHARMa)
library(vegan)
library(MASS)
library(sjPlot)
library(sjmisc)
library(cowplot)
library(arm)
library(languageR)
library(lmerTest)
library("blmeco")
library("sjstats")
library(viridisLite)  # palette viridis
library(gridExtra)
library(RColorBrewer)
library(ggeffects)

# install.packages(c("ggplot2", "lme4", "MuMIn", "DHARMa", "vegan", "MASS", "sjPlot", "sjmisc", "grid", "cowplot"))
# install.packages(c("arm", "languageR", "lmerTest", "blmeco", "sjstats", "viridisLite", "stringr", "gridExtra", "RColorBrewer", "tidyverse"))
# remove.packages(c("arm", "languageR", "lmerTest", "blmeco", "sjstats", "viridisLite", "stringr", "gridExtra", "RColorBrewer", "tidyverse"), lib="~/R/R-3.5.2/library")


# overdisp_fun -> residual dispersion estimation with a glmer.nb
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  (rdf <- nrow(model@frame)-model.df)
  rp <- residuals(model)
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE,log.p=TRUE)
  c(chisq=Pearson.chisq,ratio=prat,p=exp(pval),logp=pval)
}

