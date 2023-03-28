# Script for comparing Beta diversity WITHOUT PHYLOGENETIC METHODS AS 16S TREE IS IFFY
# 1. With filtered ITS data as a check jaccard and bray curtis
# permanova of beta diversity values
set.seed(18)

library(ggplot2)
library(phyloseq)
library(plotrix)
library(vegan)
library(tidyverse)
library(devtools)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggplot2)
library(phyloseq)
library(gridExtra)
library(tidyr)
#library(microbiome) #hmmmmm
library(ape)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(ggpubr)
library(rbiom)

#turn off scientific notation
options(scipen=999)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")

# Read in the data

metadata <- read.csv("QZA_Files/MetadataKey5.txt", header = TRUE, sep = "\t")

load("phy_filtered_bacteria.rdata")

phy_filtered_bacteria
# have to remove the sloth fur sample
phy_filtered = subset_samples(phy_filtered_bacteria, tissue != "Sloth fur")
phy_filtered = transform_sample_counts(phy_filtered, function(x) (x / sum(x) * 10000))

#phy_filtered <- rarefy_even_depth(phy_filtered, sample.size = 1000, replace = TRUE)

# Color by tissue type for a nice comparison
all_dist <- ordinate(phy_filtered, "PCoA", "wunifrac")

tissue_bc <- plot_ordination(phy_filtered, all_dist, type="samples", color = "tissue") + geom_point(size=5) +
  ggtitle("Tissue Type Weigted-Unifrac") + scale_color_manual(values = c("Soil" = "brown",
                                                                         "Root wash" = "cyan 3",
                                                                         "Root" = "gold3",
                                                                         "Leaf" = "green4")) 

gd <- (tissue_bc$data)
# get a better color pallet because there's a lot
manualcolors <- c('black','forestgreen', 'red2', 'orange', 'cornflowerblue', 
                'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                "tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
                'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
                "darkolivegreen1" ,"tan2" ,   "tomato3" , "gold3","gainsboro")



soil <- subset_samples(phy_filtered, tissue=="Soil")
rootwash <- subset_samples(phy_filtered, tissue=="Root wash")
root <- subset_samples(phy_filtered, tissue=="Root")
leaf <- subset_samples(phy_filtered, tissue=="Leaf")

soil_dist <- ordinate(soil, "PCoA", "wunifrac")
sbc <- plot_ordination(soil, soil_dist, type="samples", color = "location") + geom_point(size=5) +
  ggtitle("Soil x Weigted-Unifrac") + theme(legend.position = "none") + scale_color_manual(values = manualcolors)

rw_dist <- ordinate(rootwash, "PCoA", "wunifrac")
rwbc <- plot_ordination(rootwash, rw_dist, type="samples", color = "location") + geom_point(size=5) +
  ggtitle("RootWash x Weigted-Unifrac") + theme(legend.position = "none") + scale_color_manual(values = manualcolors)

root_dist <- ordinate(root, "PCoA", "wunifrac")
rootbc <- plot_ordination(root, root_dist, type="samples", color = "location") + geom_point(size=5) +
  ggtitle("Roots x Weigted-Unifrac") + theme(legend.position = "none") + scale_color_manual(values = manualcolors)

leaf_dist <- ordinate(leaf, "PCoA", "wunifrac")
leafbc <- plot_ordination(leaf, leaf_dist, type="samples", color = "location") + geom_point(size=5) +
  ggtitle("Leaf x Weigted-Unifrac") + theme(legend.position = "none") + scale_color_manual(values = manualcolors)


beta_diversity <- ggarrange(sbc, rwbc, rootbc, leafbc, ncol = 2, nrow = 2, labels = c("A","B","C","D"),
                            common.legend = TRUE, legend = c("right"))
beta_plot <- ggarrange(tissue_bc,beta_diversity,nrow = 2)
beta_plot

ggsave("Bact_wunifrac_relabund_BetaDiversity.png", plot = beta_plot, path = "Results_Figs_Tables/Quick_Figures", dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")

# Statistics 

Bray_dist_all <- phyloseq::distance(phy_filtered,"wunifrac")

# can only test for tissue and location as all weather and soil is nested within location
adonis2(Bray_dist_all ~ sample_data(phy_filtered)$tissue +
          sample_data(phy_filtered)$temp +
          sample_data(phy_filtered)$precip +
          sample_data(phy_filtered)$Texture +
          sample_data(phy_filtered)$cloudcover +
          sample_data(phy_filtered)$solarradiation +
          sample_data(phy_filtered)$Soil.Temp +
          sample_data(phy_filtered)$Soil.Moisture +
          sample_data(phy_filtered)$Soil.pH,
        by = "margin",na.action = na.omit)

adonis2(Bray_dist_all ~ sample_data(phy_filtered)$tissue +
          sample_data(phy_filtered)$location +
          sample_data(phy_filtered)$collection_date +
          sample_data(phy_filtered)$plant_stage,
        by = "margin")

