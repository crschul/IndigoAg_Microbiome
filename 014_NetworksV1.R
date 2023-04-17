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
                      classification = 'Phylum', buffer = 20) + 
  ggtitle("Co_occurance Network of Bacteria and Fungal Phylums across the US")

net_fig_str <- ggplot_build(net_fig)
  

# Try piecing together our own

# Tomorrow - Relative abundance of phyla drop anything smaller than .2

co_table <- co_occurrence(combined_phylum, treatment = c('Location', 'tissue'), method = 'spearman',
                          cores = 4, rho = 0.85)

network_layout_obj <- network_layout_ps(combined_phylum,treatment = c("Location","tissue"),
                                        co_occurrence_table = co_table,
                                        algorithm = 'circle')
# filter rare taxa through col in network layout obj

layout_abundant <- filter(network_layout_obj, `Mean Relative Abundance` != "(0,0.04]")

all_abundant <- co_occurrence_network(combined_phylum,treatment = c("Location","tissue"),
                      co_occurrence_table = co_table,
                      layout = layout_abundant,
                      negative_positive_colors = c('gray22','tomato3'),
                      classification = 'Phylum') + 
  ggtitle("All Tissue Abundant Taxa")

all_abundant <- all_abundant + geom_label(aes(x = all_abundant$data$x,y=all_abundant$data$y,
                                             label = all_abundant$data$Phylum),
                                         size = 6,nudge_y = .05,
                                         fontface="bold",
                                         vjust="inward",hjust="inward")
all_abundant


layout_rare <- filter(network_layout_obj, `Mean Relative Abundance` == "(0,0.04]")

co_occurrence_network(combined_phylum,treatment = c("Location","tissue"),
                      co_occurrence_table = co_table,
                      layout = layout_rare,
                      classification = 'Phylum')

ggsave("BacteriaSankeyRelAbun_POSTER.png", plot = flow_and_relabund, path = "Results_Figs_Tables/Quick_Figures", dpi = 700, 
       width = 30, height = 10, units = c("in"), device = "png")

library(igraph)
g <- graph(network_layout_obj)

l <- layout_in_circle(network_layout_obj)
plot(network_layout_obj,layout = l)



write.csv(igraph::as_edgelist(network_layout_obj, names = T),file = "igraph_net.csv",sep = ",")


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

### all this with picrust 


###### Final steps, get a measure of each network, and do a MLM with a covariate table of env data
### maybe network properties, connectivity is what I really want. 

# Dan - correlated pulled soil data with existing soil data see if I can trust and extrapolate

# NetCoMi would work really well but it sucks at running? Try to fix. 



# per sample %>% netConstruct %>% netAnalyze %>% net analyze %>% grab data

# Try her example
library(NetCoMi)
data("soilrep")

soil_warm_yes <- phyloseq::subset_samples(soilrep, warmed == "yes")
soil_warm_no  <- phyloseq::subset_samples(soilrep, warmed == "no")

net_seas_p <- netConstruct(soil_warm_yes, soil_warm_no,
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 500),
                           zeroMethod = "pseudo",
                           normMethod = "clr",
                           measure = "pearson",
                           verbose = 0)


netprops1 <- netAnalyze(net_seas_p, clustMethod = "cluster_fast_greedy")

summary(netprops1)

nclust <- as.numeric(max(names(table(netprops1$clustering$clust1))))
col <- c(topo.colors(nclust), rainbow(6))



# my data
fung_phy <- phyloseq::subset_samples(phy_combined, Microbe == "fung")
bac_phy  <- phyloseq::subset_samples(phy_combined, Microbe == "bacteria")

phylum_combined <- tax_glom(phy_combined,"Phylum")

combined_net <- netConstruct(phylum_combined,
                             measure = "pearson",
                             zeroMethod = "pseudo",
                             normMethod = "clr",
                             verbose = 2,
                             taxRank = "Phylum")



net_anal <- netAnalyze(combined_net, clustMethod = "cluster_fast_greedy")

summary(net_anal)


plot(net_anal,
     layout = "circle",
     nodeColor = "gray")


plot(nat_anal,
     sameLayout = TRUE, 
     layoutGroup = "union", 
     colorVec = col,
     borderCol = "gray40", 
     nodeSize = "degree", 
     cexNodes = 0.9, 
     nodeSizeSpread = 3, 
     edgeTranspLow = 80, 
     edgeTranspHigh = 50, 
     showTitle = TRUE, 
     cexTitle = 2.8,
     mar = c(1,1,3,1), 
     repulsion = 0.9, 
     labels = FALSE, 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 10, 
     nodeTransp = 50, 
     hubTransp = 30)



# plot(netprops1, 
#      sameLayout = TRUE, 
#      layoutGroup = "union", 
#      colorVec = col,
#      borderCol = "gray40", 
#      nodeSize = "degree", 
#      cexNodes = 0.9, 
#      nodeSizeSpread = 3, 
#      edgeTranspLow = 80, 
#      edgeTranspHigh = 50,
#      groupNames = c("Warming", "Non-warming"), 
#      showTitle = TRUE, 
#      cexTitle = 2.8,
#      mar = c(1,1,3,1), 
#      repulsion = 0.9, 
#      labels = FALSE, 
#      rmSingles = "inboth",
#      nodeFilter = "clustMin", 
#      nodeFilterPar = 10, 
#      nodeTransp = 50, 
#      hubTransp = 30)

