# Script for looking at core microbiome
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
library(microbiome)1


#turn off scientific notation
options(scipen=999)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")


load("phy_combined.rdata")
metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")

# make compositional
phy_combined_comp <- microbiome::transform(phy_combined, "compositional")

fung_data <- subset_samples(phy_combined_comp, Microbe == "fung")
bact_data <- subset_samples(phy_combined_comp, Microbe == "bacteria")

bact_data_g <- tax_glom(bact_data, taxrank="Genus")

# subset all the data
soil_f <- subset_samples(fung_data, tissue == "Soil")
rootwash_f <- subset_samples(fung_data, tissue == "Root wash")
root_f <- subset_samples(fung_data, tissue == "Root")
leaf_f <- subset_samples(fung_data, tissue == "Leaf")

soil_b <- subset_samples(bact_data_g, tissue == "Soil")
rootwash_b <- subset_samples(bact_data_g, tissue == "Root wash")
root_b <- subset_samples(bact_data_g, tissue == "Root")
leaf_b <- subset_samples(bact_data_g, tissue == "Leaf")


# core tables
core_soil_f <- as.data.frame(tax_table(soil_f))[rownames(as.data.frame(tax_table(soil_f))) %in% 
                                                  microbiome::core_members(soil_f, detection = 0, prevalence = .9),]
core_rootwash_f <- as.data.frame(tax_table(rootwash_f))[rownames(as.data.frame(tax_table(rootwash_f))) %in% 
                                                  microbiome::core_members(rootwash_f, detection = 0, prevalence = .9),]
core_root_f <- as.data.frame(tax_table(root_f))[rownames(as.data.frame(tax_table(root_f))) %in% 
                                                  microbiome::core_members(root_f, detection = 0, prevalence = .9),]
core_leaf_f <- as.data.frame(tax_table(leaf_f))[rownames(as.data.frame(tax_table(leaf_f))) %in% 
                                                  microbiome::core_members(leaf_f, detection = 0, prevalence = .9),]


# no core ASVs for bacteria
core_soil_b <- as.data.frame(tax_table(soil_b))[rownames(as.data.frame(tax_table(soil_b))) %in% 
                                                  microbiome::core_members(soil_b, detection = 0, prevalence = .9),]
core_rootwash_b <- as.data.frame(tax_table(rootwash_b))[rownames(as.data.frame(tax_table(rootwash_b))) %in% 
                                                  microbiome::core_members(rootwash_b, detection = 0, prevalence = .9),]
core_root_b <- as.data.frame(tax_table(root_b))[rownames(as.data.frame(tax_table(root_b))) %in% 
                                                  microbiome::core_members(root_b, detection = 0, prevalence = .9),]
core_leaf_b <- as.data.frame(tax_table(leaf_b))[rownames(as.data.frame(tax_table(leaf_b))) %in% 
                                                  microbiome::core_members(leaf_b, detection = 0, prevalence = .9),]
# n genus
fung_data # 
bact_data_g #

# more core fungi than bacteria, what kind of figure should I use? or just a table?
# core picrust next? 



