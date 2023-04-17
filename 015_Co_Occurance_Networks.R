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
#devtools::install_github('schuyler-smith/phylosmith')
library(phylosmith)
library(DESeq2)

#turn off scientific notation
options(scipen=999)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")


load("phy_combined.rdata")
metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")


################### Normalization
# pseudo counts -> normalize -> correlate

# pseudo counts + relative abundance: 20,000 things (suspicious)
# loose relative abundance: CLR, VST (DESEQ)
#phy_pseudo <- transform_sample_counts(phy_combined, function(OTU) OTU +1)

# TSS is just relative abundance
phy_filtered = transform_sample_counts(phy_combined, function(x) (x / sum(x) * 10000))

phy_pseudo <- transform_sample_counts(phy_filtered, function(OTU) OTU +1)

#phy_clr <- microbiome::transform(phy_pseudo, "clr")

# Test the whole network across everything!
combined_phylum <- conglomerate_taxa(phy_pseudo, "Phylum")


##### Rare and Abundant Co_Occurrence Networks


#####################################################################
# Look at networks across ALL tissues

# Correlation Table of all phylum, with a rho cutoff of .85
co_table <- co_occurrence(combined_phylum, treatment = c('Location'), method = 'spearman',
                          cores = 4, rho = 0.8)
# The actual network layout
network_layout_obj <- network_layout_ps(combined_phylum,treatment = c("Location"),
                                        co_occurrence_table = co_table)

# filter abundant and rare taxa through col in network layout obj

layout_abundant <- filter(network_layout_obj, `Mean Relative Abundance` != "(0,0.03]")
layout_rare <- filter(network_layout_obj, `Mean Relative Abundance` == "(0,0.03]")

# All

all_all <- co_occurrence_network(combined_phylum,treatment = c("Location"),
                                 co_occurrence_table = co_table,
                                 layout = network_layout_obj,
                                 negative_positive_colors = c('tomato3','grey22'), #Need to change if only positive or negative
                                 classification = 'Phylum') + 
  ggtitle("All Tissue Abundant Taxa")
all_all

ggsave("CO_Network_All_AllTissue.png", plot = all_all, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")

all_abundant <- co_occurrence_network(combined_phylum,treatment = c("Location"),
                                      co_occurrence_table = co_table,
                                      layout = layout_abundant,
                                      negative_positive_colors = c('grey22'), #Need to change if only positive or negative
                                      classification = 'Phylum') + 
  ggtitle("All Tissue Abundant Taxa")

# Remove labels for now, can always mess with them another time
# all_abundant <- all_abundant + geom_label(aes(x = all_abundant$data$x,y=all_abundant$data$y,
#                                               label = all_abundant$data$Phylum),
#                                           size = 6,nudge_y = .05,
#                                           fontface="bold",
#                                           vjust="inward",hjust="inward")
all_abundant

ggsave("CO_Network_Abundant_AllTissue.png", plot = all_abundant, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")


# Rare

all_rare <- co_occurrence_network(combined_phylum,treatment = c("Location","tissue"),
                                      co_occurrence_table = co_table,
                                      layout = layout_rare,
                                      negative_positive_colors = c('tomato3','grey22'), #Need to change if only positive or negative
                                      classification = 'Phylum') + 
  ggtitle("All Tissue Rare Taxa")

# Too many to label
# all_rare <- all_rare + geom_label(aes(x = all_rare$data$x,y=all_rare$data$y,
#                                               label = all_rare$data$Phylum),
#                                           size = 2,nudge_y = .05,
#                                           fontface="bold", check_overlap = TRUE,
#                                           vjust="inward",hjust="inward")
all_rare

ggsave("CO_Network_Rare_AllTissue.png", plot = all_rare, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")




#####################################################################
# Look at networks across soil

soil_d <- subset_samples(combined_phylum, tissue == "Soil")

# Correlation Table of all phylum, with a rho cutoff of .85
co_soil <- co_occurrence(soil_d, treatment = c('Location'), method = 'spearman',
                          cores = 4, rho = 0.8)
# The actual network layout
network_layout_soil <- network_layout_ps(soil_d,treatment = c("Location"),
                                        co_occurrence_table = co_soil,
                                        algorithm = 'circle')

# filter abundant and rare taxa through col in network layout obj

layout_abundant <- filter(network_layout_soil, `Mean Relative Abundance` != "(0,0.03]")
layout_rare <- filter(network_layout_soil, `Mean Relative Abundance` == "(0,0.03]")


# all
soil_all <- co_occurrence_network(soil_d,treatment = c("Location"),
                                       co_occurrence_table = co_soil,
                                       layout = network_layout_soil,
                                       negative_positive_colors = c('tomato3','grey22'), #Need to change if only positive or negative
                                       classification = 'Phylum') + 
  ggtitle("Soil All Taxa")

soil_all

ggsave("CO_Network_All_Soil.png", plot = soil_all, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")

# Abundant
soil_abundant <- co_occurrence_network(soil_d,treatment = c("Location"),
                                      co_occurrence_table = co_soil,
                                      layout = layout_abundant,
                                      negative_positive_colors = c('grey22'), #Need to change if only positive or negative
                                      classification = 'Phylum') + 
  ggtitle("Soil Abundant Taxa")

soil_abundant

ggsave("CO_Network_Abundant_Soil.png", plot = soil_abundant, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")


# Rare

soil_rare <- co_occurrence_network(soil_d,treatment = c("Location"),
                                  co_occurrence_table = co_soil,
                                  layout = layout_rare,
                                  negative_positive_colors = c('tomato3','grey22'), #Need to change if only positive or negative
                                  classification = 'Phylum') + 
  ggtitle("Soil Rare Taxa")

soil_rare

ggsave("CO_Network_Rare_Soil.png", plot = soil_rare, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")





#####################################################################
# Look at networks across Root Wash

rw_d <- subset_samples(combined_phylum, tissue == "Root wash")

# Correlation Table of all phylum, with a rho cutoff of .85
co_rw <- co_occurrence(rw_d, treatment = c('Location'), method = 'spearman',
                         cores = 4, rho = 0.8)
# The actual network layout
network_layout_rw <- network_layout_ps(rw_d,treatment = c("Location"),
                                         co_occurrence_table = co_rw,
                                         algorithm = 'circle')

# filter abundant and rare taxa through col in network layout obj

layout_abundant <- filter(network_layout_rw, `Mean Relative Abundance` != "(0,0.03]")
layout_rare <- filter(network_layout_rw, `Mean Relative Abundance` == "(0,0.03]")

# all
rw_all <- co_occurrence_network(rw_d,treatment = c("Location"),
                                     co_occurrence_table = co_rw,
                                     layout = network_layout_rw,
                                     negative_positive_colors = c('tomato3','grey22'), #Need to change if only positive or negative
                                     classification = 'Phylum') + 
  ggtitle("Root Wash All Taxa")

rw_all

ggsave("CO_Network_All_RW.png", plot = rw_all, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")


# Abundant
rw_abundant <- co_occurrence_network(rw_d,treatment = c("Location"),
                                       co_occurrence_table = co_rw,
                                       layout = layout_abundant,
                                       negative_positive_colors = c('grey22'), #Need to change if only positive or negative
                                       classification = 'Phylum') + 
  ggtitle("Root Wash Abundant Taxa")

rw_abundant

ggsave("CO_Network_Abundant_RW.png", plot = rw_abundant, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")


# Rare

rw_rare <- co_occurrence_network(rw_d,treatment = c("Location"),
                                   co_occurrence_table = co_rw,
                                   layout = layout_rare,
                                   negative_positive_colors = c('tomato3','grey22'), #Need to change if only positive or negative
                                   classification = 'Phylum') + 
  ggtitle("Root Wash Rare Taxa")

rw_rare

ggsave("CO_Network_Rare_RW.png", plot = rw_rare, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")




#####################################################################
# Look at networks across Roots

root_d <- subset_samples(combined_phylum, tissue == "Root")

# Correlation Table of all phylum, with a rho cutoff of .85
co_root <- co_occurrence(root_d, treatment = c('Location'), method = 'spearman',
                       cores = 4, rho = 0.8)
# The actual network layout
network_layout_root <- network_layout_ps(root_d,treatment = c("Location"),
                                       co_occurrence_table = co_root,
                                       algorithm = "circle")

# filter abundant and rare taxa through col in network layout obj

layout_abundant <- filter(network_layout_root, `Mean Relative Abundance` != "(0,0.04]")
layout_rare <- filter(network_layout_root, `Mean Relative Abundance` == "(0,0.04]")


# All
root_all <- co_occurrence_network(root_d,treatment = c("Location"),
                                       co_occurrence_table = co_root,
                                       layout = network_layout_root,
                                       negative_positive_colors = c('tomato3','grey22'), #Need to change if only positive or negative
                                       classification = 'Phylum') + 
  ggtitle("Root All Taxa")

root_all

ggsave("CO_Network_All_Root.png", plot = root_all, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")


# Abundant
root_abundant <- co_occurrence_network(root_d,treatment = c("Location"),
                                     co_occurrence_table = co_root,
                                     layout = layout_abundant,
                                     negative_positive_colors = c('grey22'), #Need to change if only positive or negative
                                     classification = 'Phylum') + 
  ggtitle("Root Abundant Taxa")

root_abundant

ggsave("CO_Network_Abundant_Root.png", plot = root_abundant, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")


# Rare

root_rare <- co_occurrence_network(root_d,treatment = c("Location"),
                                 co_occurrence_table = co_root,
                                 layout = layout_rare,
                                 negative_positive_colors = c('tomato3','grey22'), #Need to change if only positive or negative
                                 classification = 'Phylum') + 
  ggtitle("Root Rare Taxa")

root_rare

ggsave("CO_Network_Rare_Root.png", plot = root_rare, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")



#####################################################################
# Look at networks across Leafs

leaf_d <- subset_samples(combined_phylum, tissue == "Leaf")

# Correlation Table of all phylum, with a rho cutoff of .85
co_leaf <- co_occurrence(leaf_d, treatment = c('Location'), method = 'pearson',
                         cores = 4, rho = 0.2)
# The actual network layout
network_layout_leaf <- network_layout_ps(leaf_d,treatment = c("Location"),
                                         co_occurrence_table = co_leaf,
                                         algorithm = 'circle')

# filter abundant and rare taxa through col in network layout obj

layout_abundant <- filter(network_layout_leaf, `Mean Relative Abundance` != "(0,0.04]")
layout_rare <- filter(network_layout_leaf, `Mean Relative Abundance` == "(0,0.04]")

# all
leaf_all <- co_occurrence_network(leaf_d,treatment = c("Location"),
                                       co_occurrence_table = co_leaf,
                                       layout = network_layout_leaf,
                                       negative_positive_colors = c('tomato3','grey22'), #Need to change if only positive or negative
                                       classification = 'Phylum') + 
  ggtitle("Leaf All Taxa")

leaf_all

ggsave("CO_Network_All_Leaf.png", plot = leaf_all, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")




# Abundant
leaf_abundant <- co_occurrence_network(leaf_d,treatment = c("Location"),
                                       co_occurrence_table = co_leaf,
                                       layout = layout_abundant,
                                       negative_positive_colors = c('grey22'), #Need to change if only positive or negative
                                       classification = 'Phylum') + 
  ggtitle("Leaf Abundant Taxa")

leaf_abundant

ggsave("CO_Network_Abundant_Leaf.png", plot = leaf_abundant, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")


# Rare

leaf_rare <- co_occurrence_network(leaf_d,treatment = c("Location"),
                                   co_occurrence_table = co_leaf,
                                   layout = layout_rare,
                                   negative_positive_colors = c('tomato3','gray22'), #Need to change if only positive or negative
                                   classification = 'Phylum') + 
  ggtitle("Leaf Rare Taxa")

leaf_rare

ggsave("CO_Network_Rare_Leaf.png", plot = leaf_rare, 
       path = "Results_Figs_Tables/Quick_Figures/Quick_CO_Networks"
       , dpi = 700, 
       width = 10, height = 10, units = c("in"), device = "png")



############# I think there is a huge problem breaking up the rare and abundant taxa



