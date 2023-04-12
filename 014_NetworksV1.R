# Script for looking at microbiome networks
# Combined Data

set.seed(18)

library(tidyverse)
library(devtools)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggplot2)
library(phyloseq)
library(gridExtra)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(ggpubr)
#install_github("microbiome/microbiome")
library(microbiome)

#turn off scientific notation
options(scipen=999)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")


load("phy_combined.rdata")
metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")

#devtools::install_github('schuyler-smith/phylosmith')
library(phylosmith)



# Test the whole network across everything!
combined_phylum <- conglomerate_taxa(phy_combined, "Phylum")


# Combine them all
net_fig <- co_occurrence_network(combined_phylum, treatment = c("Location","tissue"),
                      classification = 'Phylum', buffer = 5) + 
  ggtitle("Co_occurance Network of Bacteria and Fungal Phylums across the US")

# The table we want: now what
co_df <- co_occurrence(combined_phylum, treatment = c('Location', 'tissue'), method = 'spearman',
                                cores = 4)


# Correlation to variables - 
#I want correlations of microbes to environmental factors within tissue
# can only use env variables I pulled from the weather API

#filtered_combined <- taxa_filter(phy_combined, frequency = 0.1)

vc_table_genus_TEMP <- variable_correlation(phy_combined, variables = "temp",
                             classification = "Genus", 
                           treatment = "tissue",
                             method = 'spearman',cores = 4)

# Maybe have a bunch of them showing how temp is higher correlated and has phylum with more relative abundance?
temp_corr <- variable_correlation_network(phy_combined, variables = "temp",
                     classification = "Phylum", 
                     treatment = "tissue",
                     method = 'spearman',
                     p_threshold = .05, rho_threshold = c(-0.01, 0.01)) +
  ggtitle(" Correlation between Precipitation and Phylum")

precip_corr <- variable_correlation_network(phy_combined, variables = "precip",
                             classification = "Phylum", 
                             treatment = "tissue",
                             method = 'spearman',
                             p_threshold = .05, rho_threshold = c(-0.01, 0.01)) + 
  ggtitle(" Correlation between Precipitation and Phylum")

precip_corr_heat <- variable_correlation_heatmap(phy_combined, variables = "precip",
                                                 classification = "Phylum", 
                                                 treatment = "tissue",
                                                 method = 'spearman', cores = 4, 
                                                 significance_color = 'black',
                                                 colors = c("#2C7BB6", "white", "#D7191C")) + 
  ggtitle(" Correlation between Precipitation and Phylum")

###### Final step, get a measure of each network, and do a MLM with a covariate table of env data
### maybe network properties, connectivity is what I really want. 

# NetCoMi would work really well but it sucks at running? Try to fix. 



# per sample %>% netConstruct %>% netAnalyze %>% net analyze %>% grab data

# Try her example


