# Script for looking at environment correlation heatmaps
# Combined Data
# Subset for samples that actually have full metadata
#Benjamini & Hochberg (1995) correction through p.adjust

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

table7 <- read.csv("Results_Figs_Tables/Fig7Table.csv")

f7b <- ggplot(table7, aes(x = Tissue.1, y = Tissue.2, fill = Bacteria.KEGG.Terms)) + 
  geom_tile(color="black") + 
  scale_fill_gradient2(low = "red", mid = "red",high = "blue") + theme_bw() +
  geom_text(aes(label = Bacteria.KEGG.Terms))

f7b


f7f <- ggplot(table7, aes(x = Tissue.2, y = Tissue.1, fill = Fungi.EC.terms)) + 
  geom_tile(color="black") + 
  scale_fill_gradient2(low = "red", mid = "red",high = "blue") + theme_bw() +
  geom_text(aes(label = Fungi.EC.terms))

f7f
