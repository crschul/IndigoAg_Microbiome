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

legend(.3, 1.1, cex = 1.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("blue","darkorange"), 
       bty = "n", horiz = TRUE)
legend(.6, .9, cex = 1.2, pt.cex = 2.5, title = "Kingdom:", 
       legend=c("Fungi","Bacteria"), col = c("tan","lightgreen"), bty = "n", pch = 16) 



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
                      title1 = "Network on Phylum level: Soil",
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

legend(.55, 1.1, cex = 1.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("blue","darkorange"), 
       bty = "n", horiz = TRUE)
legend(.8, .9, cex = 1.2, pt.cex = 2.5, title = "Kingdom:", 
       legend=c("Fungi","Bacteria"), col = c("tan","lightgreen"), bty = "n", pch = 16)



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
                      title1 = "Network on Phylum level: Root Wash",
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
                      colorVec = c("Fungi" = "purple","d__Bacteria" = "lightgreen", "d__Archaea" = "tan" )) # why does it call my fungi archea?

legend(.5, 1.1, cex = 1.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("blue","darkorange"), 
       bty = "n", horiz = TRUE)
legend(.9, .9, cex = 1.0, pt.cex = 2.5, title = "Kingdom:", 
       legend=c("Fungi","Bacteria","Archaea"), col = c("tan","lightgreen","purple"), bty = "n", pch = 16)


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
                      title1 = "Network on Phylum level: Roots",
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
                      colorVec = c("Fungi" = "purple","d__Bacteria" = "lightgreen", "d__Archaea" = "tan" )) # why does it call my fungi archea?

legend(.5, 1.1, cex = 1.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("blue","darkorange"), 
       bty = "n", horiz = TRUE)
legend(.7, -.5, cex = 1.2, pt.cex = 2.5, title = "Kingdom:", 
       legend=c("Fungi","Bacteria"), col = c("tan","lightgreen"), bty = "n", pch = 16)




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
                      title1 = "Network on Phylum level: Leaf",
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

legend(.4, 1.1, cex = 1.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("blue","darkorange"), 
       bty = "n", horiz = TRUE)
legend(.4, .9, cex = 1.2, pt.cex = 2.5, title = "Kingdom:", 
       legend=c("Fungi","Bacteria"), col = c("tan","lightgreen"), bty = "n", pch = 16) 















###### Testing below: ignore!!!
################################################### compare both methods!
#### Leaf data - NetCoMi
combined_net <- netConstruct(Leaf_phy,
                             measure = "spearman",
                             normMethod = "TSS",
                             zeroMethod = "pseudo",
                             thresh = 0.8,
                             verbose = 2,
                             taxRank = "Phylum")

combined_anal <- netAnalyze(combined_net, clustMethod = "cluster_fast_greedy")

#summary(net_anal)

combined_plot <- plot(combined_anal,
                      nodeColor = "grey",
                      title1 = "Network on Phylum level: Leaf",
                      showTitle = TRUE,
                      rmSingles = TRUE,
                      labelScale = FALSE,
                      cexLabels = .85,
                      posCol = "blue", 
                      negCol = "darkorange",
                      highlightHubs=TRUE,
                      nodeSize = "normCounts")

legend(.6, 1.1, cex = 1.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("blue","darkorange"), 
       bty = "n", horiz = TRUE)


##### Leaf Data PhyloSmith
# pseudo counts + relative abundance: 20,000 things (suspicious)
# loose relative abundance: CLR, VST (DESEQ)
#phy_pseudo <- transform_sample_counts(phy_combined, function(OTU) OTU +1)

# TSS is just relative abundance
phy_filtered = transform_sample_counts(phy_combined, function(x) (x / sum(x) * 10000))

phy_pseudo <- transform_sample_counts(phy_filtered, function(OTU) OTU +1)

#phy_clr <- microbiome::transform(phy_pseudo, "clr")

# Test the whole network across everything!
combined_phylum <- conglomerate_taxa(phy_pseudo, "Phylum")



# Look at networks across Leafs
leaf_d <- subset_samples(combined_phylum, tissue == "Leaf")

# Correlation Table of all phylum, with a rho cutoff of .85
co_leaf <- co_occurrence(leaf_d, treatment = c('Location'), method = 'pearson',
                         cores = 4, rho = 0.2) # Drop because it just wont run
# The actual network layout
network_layout_leaf <- network_layout_ps(leaf_d,treatment = c("Location"),
                                         co_occurrence_table = co_leaf,
                                         algorithm = 'circle')

# filter abundant and rare taxa through col in network layout obj

# Abundant
leaf_abundant <- co_occurrence_network(leaf_d,treatment = c("Location"),
                                       co_occurrence_table = co_leaf,
                                       layout = network_layout_leaf,
                                       negative_positive_colors = c('tomato3','grey22'), #Need to change if only positive or negative
                                       classification = 'Phylum') + 
  ggtitle("Leaf Network")

leaf_abundant

####################### Figure out the color and shape with their data

# Load data sets from American Gut Project (from SpiecEasi package)
data("amgut1.filt")

# Network construction
amgut_net <- netConstruct(amgut1.filt, measure = "pearson",
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 50),
                          zeroMethod = "pseudoZO", normMethod = "clr",
                          sparsMethod = "threshold", thresh = 0.3)

# Network analysis
amgut_props <- netAnalyze(amgut_net)

### Network plots ###
# Clusters are used for node coloring: 
plot(amgut_props, 
     nodeColor = "cluster")

# Remove singletons
plot(amgut_props, 
     nodeColor = "cluster", 
     rmSingles = TRUE)

# A higher repulsion places nodes with high edge weight closer together
plot(amgut_props, 
     nodeColor = "cluster", 
     rmSingles = TRUE, 
     repulsion = 1.2)



kingvec <- as.data.frame(tax_table(Leaf_phy))$Kingdom
names(kingvec) <- as.data.frame(tax_table(Leaf_phy))$Phylum




# A feature vector is used for node coloring
# (this could be a vector with phylum names of the ASVs)
set.seed(123456)
featVec <- sample(1:5, nrow(amgut1.filt), replace = TRUE)

# Names must be equal to ASV names
names(featVec) <- colnames(amgut1.filt)

featVec
plot(amgut_props, 
     rmSingles = TRUE, 
     nodeColor = "feature", 
     featVecCol = featVec,
     colorVec = heat.colors(5))

# Use a further feature vector for node shapes
shapeVec <- sample(1:3, ncol(amgut1.filt), replace = TRUE)
names(shapeVec) <- colnames(amgut1.filt)

plot(amgut_props,
     rmSingles = TRUE,
     nodeColor = "feature",
     featVecCol = featVec,
     colorVec = heat.colors(5), 
     nodeShape = c("circle", "square", "diamond"),
     featVecShape = shapeVec, 
     highlightHubs = FALSE)




############################## Try to convert this to ggnet
library(GGally)
library(network)
library(sna)

#### Leaf data - NetCoMi
leaf_net_test <- netConstruct(Leaf_phy,
                             measure = "cclasso",
                             normMethod = "TSS",
                             zeroMethod = "pseudo",
                             thresh = 0.8,
                             verbose = 2,
                             taxRank = "Phylum")

leaf_edge <- leaf_net_test$edgelist1
leaf_assoc <- leaf_net_test$assoMat1
leaf_adj <- leaf_net_test$adjaMat1


ggnet2(leaf_edge, label = TRUE)

leaf_net <- network(leaf_edge,
                    ignore.eval = FALSE)



taxtab <- as(tax_table(phylum_combined), "matrix")
phyla <- as.factor(gsub("p__", "", taxtab[, "Phylum"]))
names(phyla) <- taxtab[, "Kingdom"]


col = c("Ascomycota" = "tan", "Proteobacteria" = "orange","Firmicutes" = "darkgreen", "Basdiomycota" = "purple",
        "Bacteroidota" = "darkgreen", "Verrucomicrobiota" = "darkgreen", "Acidobacteriota" = "darkgreen")

set.edge.attribute(leaf_net, "color", ifelse(leaf_net %e% "asso" > 0, "black", "red"))

set.edge.value(leaf_net, "weight", value = leaf_edge$asso)

ggnet2(leaf_net, label = TRUE, edge.color = "color",
       edge.size = 5, node.color = col)





############# Theirs
# weighted adjacency matrix
net = rgraph(10, mode = "graph", tprob = 0.5)
net = network(net, directed = FALSE)

# vertex names
network.vertex.names(net) = letters[1:10]

net %v% "phono" = ifelse(letters[1:10] %in% c("a", "e", "i"), "vowel", "consonant")

ggnet2(net, color = "phono")





