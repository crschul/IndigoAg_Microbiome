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
#library(DESeq2)
library(NetCoMi)

#turn off scientific notation
options(scipen=999)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")


load("phy_combined.rdata")
metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")

phylum_combined <- tax_glom(phy_combined,"Phylum")

Soil_phy <- subset_samples(phylum_combined, tissue == "Soil")
RW_phy <- subset_samples(phylum_combined, tissue == "Root wash")
Root_phy <- subset_samples(phylum_combined, tissue == "Root")
Leaf_phy <- subset_samples(phylum_combined, tissue == "Leaf")


# cclasso, tss, pseudo, .8 look good <- USE ME
# sparcc, clr, pseudo, .8 look good

?plot.microNetProps # Try to find the answers to making better figs

# Get Colors
kingvec <- as.data.frame(tax_table(phylum_combined))$Kingdom
names(kingvec) <- as.data.frame(tax_table(phylum_combined))$Phylum


#### Combined data
combined_net <- netConstruct(phylum_combined,
                             measure = "cclasso",
                             normMethod = "TSS",
                             zeroMethod = "pseudo",
                             thresh = 0.8,
                             verbose = 2,
                             taxRank = "Phylum")

combined_anal <- netAnalyze(combined_net, clustMethod = "cluster_fast_greedy")

#summary(net_anal)


combined_plot <- plot(combined_anal,
                      title1 = "Network on Phylum level: All Tissues",
                      showTitle = TRUE,
                      cexTitle = 1.5,
                      rmSingles = TRUE,
                      labelScale = FALSE,
                      labelFont = 2,
                      cexLabels = .85,
                      posCol = "blue", 
                      negCol = "orange",
                      highlightHubs=FALSE,
                      nodeSize = "normCounts",
                      hubBorderCol = "darkgrey",
                      cexNodes = 5,
                      edgeWidth = 1,
                      edgeTranspLow = 0,
                      edgeTranspHigh = 0,
                      cut = 0,
                      nodeColor = "feature", 
                      featVecCol = kingvec,
                      colorVec = c("Fungi" = "tan","d__Bacteria" = "lightgreen", "d__Archaea" = "tan" )) # why does it call my fungi archea?




##################################################### Soil data
combined_net <- netConstruct(Soil_phy,
                             measure = "cclasso",
                             normMethod = "TSS",
                             zeroMethod = "pseudo",
                             thresh = 0.8,
                             verbose = 2,
                             taxRank = "Phylum")

combined_anal <- netAnalyze(combined_net, clustMethod = "cluster_fast_greedy")

#summary(net_anal)

combined_plot <- plot(combined_anal,
                      showTitle = TRUE,
                      cexTitle = 1.5,
                      rmSingles = TRUE,
                      labels = FALSE,
                      cexLabels = .85,
                      posCol = "blue", 
                      negCol = "red",
                      highlightHubs=TRUE,
                      nodeSize = "normCounts",
                      hubBorderCol = "black",
                      hubBorderWidth = 4,
                      cexNodes = 5,
                      edgeWidth = 1,
                      edgeTranspLow = 0,
                      edgeTranspHigh = 0,
                      cut = 0,
                      nodeColor = "feature",
                      nodeTransp = 1,
                      featVecCol = kingvec,
                      colorVec = c("Fungi" = "purple","d__Bacteria" = "lightgreen", "d__Archaea" = "tan")) # why does it call my fungi archea?



############################################# RootWash data
combined_net <- netConstruct(RW_phy,
                             measure = "cclasso",
                             normMethod = "TSS",
                             zeroMethod = "pseudo",
                             thresh = 0.8,
                             verbose = 2,
                             taxRank = "Phylum")

combined_anal <- netAnalyze(combined_net, clustMethod = "cluster_fast_greedy")

#summary(net_anal)
combined_plot <- plot(combined_anal,
                      showTitle = TRUE,
                      cexTitle = 1.5,
                      rmSingles = TRUE,
                      labels = FALSE,
                      cexLabels = .85,
                      posCol = "blue", 
                      negCol = "red",
                      highlightHubs=TRUE,
                      nodeSize = "normCounts",
                      hubBorderCol = "black",
                      hubBorderWidth = 4,
                      cexNodes = 5,
                      edgeWidth = 1,
                      edgeTranspLow = 0,
                      edgeTranspHigh = 0,
                      cut = 0,
                      nodeColor = "feature",
                      nodeTransp = 1,
                      featVecCol = kingvec,
                      colorVec = c("Fungi" = "purple","d__Bacteria" = "lightgreen", "d__Archaea" = "tan" )) # why does it call my fungi archea?



#### Root data
combined_net <- netConstruct(Root_phy,
                             measure = "cclasso",
                             normMethod = "TSS",
                             zeroMethod = "pseudo",
                             thresh = 0.8,
                             verbose = 2,
                             taxRank = "Phylum")

combined_anal <- netAnalyze(combined_net, clustMethod = "cluster_fast_greedy")

#summary(net_anal)

#summary(net_anal)
combined_plot <- plot(combined_anal,
                      showTitle = TRUE,
                      cexTitle = 1.5,
                      rmSingles = TRUE,
                      labels = FALSE,
                      cexLabels = .85,
                      posCol = "blue", 
                      negCol = "red",
                      highlightHubs=TRUE,
                      nodeSize = "normCounts",
                      hubBorderCol = "black",
                      hubBorderWidth = 4,
                      cexNodes = 5,
                      edgeWidth = 1,
                      edgeTranspLow = 0,
                      edgeTranspHigh = 0,
                      cut = 0,
                      nodeColor = "feature",
                      nodeTransp = 1,
                      featVecCol = kingvec,
                      colorVec = c("Fungi" = "purple","d__Bacteria" = "lightgreen", "d__Archaea" = "tan" )) # why does it call my fungi archea?



#### Leaf data
combined_net <- netConstruct(Leaf_phy,
                             measure = "cclasso",
                             normMethod = "TSS",
                             zeroMethod = "pseudo",
                             thresh = 0.8,
                             verbose = 2,
                             taxRank = "Phylum")

combined_anal <- netAnalyze(combined_net, clustMethod = "cluster_fast_greedy")

#summary(net_anal)

combined_plot <- plot(combined_anal,
                      showTitle = TRUE,
                      cexTitle = 1.5,
                      rmSingles = TRUE,
                      labels = FALSE,
                      cexLabels = .85,
                      posCol = "blue", 
                      negCol = "red",
                      highlightHubs=TRUE,
                      nodeSize = "normCounts",
                      hubBorderCol = "black",
                      hubBorderWidth = 4,
                      cexNodes = 5,
                      edgeWidth = 1,
                      edgeTranspLow = 0,
                      edgeTranspHigh = 0,
                      cut = 0,
                      nodeColor = "feature",
                      nodeTransp = 1,
                      featVecCol = kingvec,
                      colorVec = c("Fungi" = "tan","d__Bacteria" = "lightgreen", "d__Archaea" = "tan" )) # why does it call my fungi archea?

legend(.4, 1.1, cex = 1.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("blue","red"), 
       bty = "n", horiz = TRUE)
legend(.4, .9, cex = 1.2, pt.cex = 2.5, title = "Kingdom:", 
       legend=c("Fungi","Bacteria","Archaea"), col = c("tan","lightgreen","purple"), bty = "n", pch = 16)


