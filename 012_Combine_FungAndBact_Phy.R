# Script for comparing alpha diversity (Observed, Simpson, and Shannon)
# 1. With filtered ITS data as a check. 
# permanova of alpha diversity values
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

#turn off scientific notation
options(scipen=999)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")

# Read in the data - Bacteria
load("phy_filtered_bacteria.rdata")

phy_filtered_bacteria

# have to remove the sloth fur sample
phy_filtered_b = subset_samples(phy_filtered_bacteria, tissue != "Sloth fur")


# Read in the data - Fungus
load("phy_filtered_fungi_genus.rdata")

phy_filtered_fungus_genus
# have to remove the sloth fur sample
phy_filtered_g = subset_samples(phy_filtered_fungus_genus, tissue != "Sloth fur")

# Metadata
metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")

combined_otu <- merge_phyloseq(otu_table(phy_filtered_g), otu_table(phy_filtered_b))
combined_taxa <- merge_phyloseq(tax_table(phy_filtered_g), tax_table(phy_filtered_b))
combined_meta <- merge_phyloseq(sample_data(phy_filtered_g), sample_data(phy_filtered_b))

phy_combined <- phyloseq(combined_otu,combined_taxa,combined_meta)
phy_combined

save(phy_combined, file = "phy_combined.rdata")


test_otu <- as.data.frame.matrix(otu_table(phy_combined))
test_otu <- as.matrix(sapply(test_otu, as.numeric))

length(unique(rownames(test_otu)))

mybinarymap <- heatmap(test_otu, Rowv=NA, Colv=NA, col = c("white","black"))


# get a better color pallet because there's a lot
manualcolors <- c('black','forestgreen', 'red2', 'orange', 'cornflowerblue', 
                  'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                  'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                  "tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
                  'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
                  "darkolivegreen1" ,"tan2" ,   "tomato3" , "gold3","gainsboro")


c_dist <- ordinate(phy_combined, "PCoA", "dpcoa")
com_bc <- plot_ordination(phy_combined, c_dist, type="samples", 
                          color = "Microbe") + geom_point(size=5) +
  ggtitle("Combined x Bray-Curtis") + scale_color_manual(values = manualcolors)
com_bc





