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

# Read in the data

metadata <- read.csv("QZA_Files/MetadataKey5.txt", header = TRUE, sep = "\t")

load("phy_filtered_bacteria.rdata")

phy_filtered_bacteria

# have to remove the sloth fur sample
phy_filtered = subset_samples(phy_filtered_bacteria, tissue != "Sloth fur")

# try with rareified data
phy_filtered <- rarefy_even_depth(phy_filtered, sample.size = 500, replace = FALSE)

soil <- subset_samples(phy_filtered, tissue=="Soil")
rootwash <- subset_samples(phy_filtered, tissue=="Root wash")
root <- subset_samples(phy_filtered, tissue=="Root")
leaf <- subset_samples(phy_filtered, tissue=="Leaf")


alpha_plot1 <- plot_richness(phy_filtered, x="tissue", measures=c("Observed", "Shannon", "Simpson"), 
                             color = "Location",
                            title = "Alpha Diversity: Location by Tissues") + geom_point(alpha = .05)
alpha_plot1


soil_alpha <- plot_richness(soil, x = "Location",measures=c("Observed", "Shannon", "Simpson"), 
                            color = "Location",
                            title = "Alpha Diversity: Bacteria - Soil") + theme(legend.position = "none") +
  geom_point(position = position_dodge(width = .5))
soil_alpha$layers <- soil_alpha$layers[-1]
soil_alpha

rootwash_alpha <- plot_richness(rootwash, x = "Location",measures=c("Observed", "Shannon", "Simpson"), 
                            color = "Location",
                            title = "Alpha Diversity: Bacteria - Root Wash") + theme(legend.position = "none") +
  geom_point(position = position_dodge(width = .5))
rootwash_alpha$layers <- rootwash_alpha$layers[-1]
rootwash_alpha

root_alpha <- plot_richness(root, x = "Location",measures=c("Observed", "Shannon", "Simpson"), 
                                color = "Location",
                                title = "Alpha Diversity: Bacteria - Roots") + theme(legend.position = "none") +
  geom_point(position = position_dodge(width = .5))
root_alpha$layers <- root_alpha$layers[-1]
root_alpha

leaf_alpha <- plot_richness(leaf, x = "Location",measures=c("Observed", "Shannon", "Simpson"), 
                            color = "Location",
                            title = "Alpha Diversity: Bacteria - Leaf") + theme(legend.position = "none") +
  geom_point(position = position_dodge(width = .5))
leaf_alpha$layers <- leaf_alpha$layers[-1]
leaf_alpha

alpha_diversity <- ggarrange(soil_alpha, rootwash_alpha, root_alpha, leaf_alpha, ncol = 1, nrow = 4, labels = c("A","B","C","D"))
alpha_diversity # Figure for finaplot

ggsave("BacteriaAlphaDiversity.png", plot = alpha_diversity, path = "Results_Figs_Tables/Quick_Figures", dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")

######## it looks like location may have an effect on root alpha diversity

### Stats
library(vegan)

# all samples
all_rich <- estimate_richness(phy_filtered ,measures = c("Observed", "Shannon", "Simpson") )
all_rich <- tibble::rownames_to_column(all_rich, "sampleid")
all_rich <- merge(all_rich,metadata)

adonis2(all_rich$Observed ~ tissue + location + plant_stage, data = all_rich,
        by = "margin",na.action = na.exclude)

adonis2(all_rich$Shannon ~ tissue +location + plant_stage, data = all_rich,
        by = "margin",na.action = na.exclude)

adonis2(all_rich$Simpson ~ tissue + location + plant_stage, data = all_rich,
        by = "margin",na.action = na.exclude)

# roots: location is significant for all rarefied and unrarified
root_rich <- estimate_richness(root ,measures = c("Observed", "Shannon", "Simpson") )
root_rich <- tibble::rownames_to_column(root_rich, "sampleid")
root_rich <- merge(root_rich,metadata)

adonis2(root_rich$Observed ~ location + plant_stage, data = root_rich,
        by = "margin",na.action = na.exclude)

adonis2(root_rich$Shannon ~ location + plant_stage, data = root_rich,
        by = "margin",na.action = na.exclude)

adonis2(root_rich$Simpson ~ location + plant_stage, data = root_rich,
        by = "margin",na.action = na.exclude)

# soil - only Observed is significant
soil_rich <- estimate_richness(soil,measures = c("Observed", "Shannon", "Simpson") )
soil_rich <- tibble::rownames_to_column(soil_rich, "sampleid")
soil_rich <- merge(soil_rich,metadata)

adonis2(soil_rich$Observed ~ location + plant_stage, data = soil_rich,
        by = "margin",na.action = na.exclude)

adonis2(soil_rich$Shannon ~ location + plant_stage, data = soil_rich,
        by = "margin",na.action = na.exclude)

adonis2(soil_rich$Simpson ~ location + plant_stage, data = soil_rich,
        by = "margin",na.action = na.exclude)

# rootwash - only observed significant
rootwash_rich <- estimate_richness(rootwash ,measures = c("Observed", "Shannon", "Simpson") )
rootwash_rich <- tibble::rownames_to_column(rootwash_rich, "sampleid")
rootwash_rich <- merge(rootwash_rich,metadata)

adonis2(rootwash_rich$Observed ~ location + plant_stage, data = rootwash_rich,
        by = "margin",na.action = na.exclude)

adonis2(rootwash_rich$Shannon ~ location + plant_stage, data = rootwash_rich,
        by = "margin",na.action = na.exclude)

adonis2(rootwash_rich$Simpson ~ location + plant_stage, data = rootwash_rich,
        by = "margin",na.action = na.exclude)

# leafs - all significant!!!
leaf_rich <- estimate_richness(leaf ,measures = c("Observed", "Shannon", "Simpson") )
leaf_rich <- tibble::rownames_to_column(leaf_rich, "sampleid")
leaf_rich <- merge(leaf_rich,metadata)

adonis2(leaf_rich$Observed ~ location + plant_stage, data = leaf_rich,
        by = "margin",na.action = na.exclude)

adonis2(leaf_rich$Shannon ~ location + plant_stage, data = leaf_rich,
        by = "margin",na.action = na.exclude)

adonis2(leaf_rich$Simpson ~ location + plant_stage, data = leaf_rich,
        by = "margin",na.action = na.exclude)

# leaf and root alpha values are significantly changed by location (rarefied and unrarified) while soil and rootwash are not!



