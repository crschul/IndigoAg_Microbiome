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

tissue_bacteria_wu <- plot_ordination(phy_filtered, all_dist, type="samples", color = "tissue") + geom_point(size=5) +
  ggtitle("Bacteria Beta Diversity - Weighted Unifrac") + scale_color_manual(values = c("Soil" = "brown",
                                                                         "Root wash" = "cyan 3",
                                                                         "Root" = "gold3",
                                                                         "Leaf" = "green4"),
                                                          labels=c('Leaf', 'Root', 'Rhizosphere', 'Soil')) +
  labs(color="Tissue") 
tissue_bacteria_wu


all_dist <- ordinate(phy_filtered, "PCoA", "bray")

tissue_bacteria_bc <- plot_ordination(phy_filtered, all_dist, type="samples", color = "tissue") + geom_point(size=5) +
  ggtitle("Bacteria Beta Diversity - Bray Curtis") + scale_color_manual(values = c("Soil" = "brown",
                                                                     "Root wash" = "cyan 3",
                                                                     "Root" = "gold3",
                                                                     "Leaf" = "green4"),
                                                          labels=c('Leaf', 'Root', 'Rhizosphere', 'Soil')) +
  labs(color="Tissue") 
tissue_bacteria_bc


#turn off scientific notation
options(scipen=999)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")

# Read in the data

metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")

load("phy_filtered_fungi_asv.rdata")

phy_filtered_fungus_asv
# have to remove the sloth fur sample
phy_filtered = subset_samples(phy_filtered_fungus_asv, tissue != "Sloth fur")
phy_filtered = transform_sample_counts(phy_filtered, function(x) (x / sum(x) * 10000))

# Color by tissue type for a nice comparison
all_dist <- ordinate(phy_filtered, "PCoA", "bray")
tissue_fung <- plot_ordination(phy_filtered, all_dist, type="samples", color = "tissue") + geom_point(size=5) +
  ggtitle("Fungal Beta Diversity - Bray Curtis") + scale_color_manual(values = c("Soil" = "brown",
                                                                     "Root wash" = "cyan 3",
                                                                     "Root" = "gold3",
                                                                     "Leaf" = "green4"),
                                                        labels=c('Leaf', 'Root', 'Rhizosphere', 'Soil')) +
  labs(color="Tissue")

tissue_fung



Beta_fig2 <- ggarrange(tissue_bacteria_bc,tissue_fung, nrow = 2,labels = c("A","B"),
                       common.legend = TRUE, legend = c("none"))
Beta_fig2

Beta_fig2 <- ggarrange(Beta_fig2,tissue_bacteria_wu, ncol = 2,labels = c("","C"),
                       common.legend = TRUE, legend = c("right"), widths = c(1,.5))
Beta_fig2


Beta_def <- ggarrange(tissue_bacteria_wu, tissue_fung, nrow = 1, ncol = 2, common.legend = TRUE, legend = c("right"))

ggsave("Fig2_BetaDiversity.png", plot = Beta_fig2, path = "Results_Figs_Tables/Figures_Draft1", dpi = 700, 
       width = 12, height = 8, units = c("in"), device = "png")

ggsave("Fig2_DEFENSE.png", plot = Beta_def, path = "Results_Figs_Tables/Figures_Draft1", dpi = 700, 
       width = 12, height = 6, units = c("in"), device = "png")


