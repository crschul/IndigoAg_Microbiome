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


load("phy_combined.rdata")
metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")

#devtools::install_github('schuyler-smith/phylosmith')
library(phylosmith)

# TSS is just relative abundance
phy_filtered = transform_sample_counts(phy_combined, function(x) (x / sum(x) * 10000))

phy_pseudo <- transform_sample_counts(phy_filtered, function(OTU) OTU +1)



#################### Change metadata columns for ploting
oldnames <- c("precip","temp","humidity",
              "cloudcover","Soil.pH","Organic.Matter",
              "Nitrate.N.ppm.N", "Potassium.ppm.K","Sulfate.S.ppm.S",
              "Calcium.ppm.Ca","Magnesium.ppm.Mg","Sodium.ppm.Na")

# if it has Soil.pH, it has our other soil variables!
phy_pseudo <- subset_samples(phy_pseudo, !is.na(Soil.pH))

colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "precip"] = "Precipitation"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "temp"] = "Temperature"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "humidity"] = "Humidity"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "cloudcover"] = "Cloud.Cover"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "Soil.pH"] = "Soil.pH"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "Soil.Moisture"] = "Soil.Moisture"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "Soil.Temp"] = "Soil Temperature"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "Organic.Matter"] = "Organic.Matter"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "Nitrate.N.ppm.N"] = "Nitrate.ppm"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "Potassium.ppm.K"] = "Potassium.ppm"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "Sulfate.S.ppm.S"] = "Sulfate.ppm"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "Calcium.ppm.Ca"] = "Calcium.ppm"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "Magnesium.ppm.Mg"] = "Magnesium.ppm"
colnames(sample_data(phy_pseudo))[colnames(sample_data(phy_pseudo)) == "Sodium.ppm.Na"] = "Sodium.ppm"



# Test the whole network across everything!
combined_phylum <- conglomerate_taxa(phy_pseudo, "Phylum")


#### API doesn't have 2017 soil measures :( so we have to clean our data
# subset for just the samples that have soil data

# how many samples do we lose if we just drop NAs across the board? 734 / 1093

sample_data(combined_phylum) <- sample_data(combined_phylum)[, -which(names(sample_data(combined_phylum)) %in%
                                                                        c("Soil.EC","UV.Light","Co2.ppm"))]

clean_phy <- combined_phylum
foo <- as.data.frame(sample_data(clean_phy))
#######################################################

# do bacteria and fungal graphs right next to each other!!!

Bacteria_phy <- subset_samples(clean_phy, Microbe == "bacteria")
Fungi_phy <- subset_samples(clean_phy, Microbe == "fung")

#### Fungus
fungi_corr_heat <- variable_correlation(Fungi_phy, 
                      variables = c("Precipitation",
                                    "Temperature",
                                    "Humidity",
                                    "Cloud.Cover",
                                    "Soil.pH",
                                    "Soil.Moisture",
                                    "Soil.Temperature",
                                    "Organic.Matter",
                                    "Nitrate.ppm",
                                    "Potassium.ppm",
                                    "Sulfate.ppm",
                                    "Calcium.ppm",
                                    "Magnesium.ppm", 
                                    "Sodium.ppm"),
                          classification = "Phylum", 
                          treatment = "tissue",
                          method = 'spearman', cores = 4)

fungi_adj <- fungi_corr_heat
fungi_adj$p <- p.adjust(fungi_adj$p, method = "BH")

fungi_adj <- filter(fungi_adj, p < .05)

corr_df_f <- fungi_adj

gf <- ggplot(corr_df_f, aes(x = Y, y = X, fill = rho)) + 
  geom_tile(color="black") + facet_wrap(~Treatment, nrow = 4) +
  scale_fill_gradient2(low = "red", mid = "white",high = "blue") + theme_bw() +
  ggtitle(("Effect of Environment on Fungal Phylum"))+ 
  theme(panel.grid.major.x = element_blank()) + 
  xlab(element_blank()) + ylab(element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) + scale_y_discrete(limits=rev)

gf


# Bacteria
bact_corr_heat <- variable_correlation(Bacteria_phy, 
                                        variables = c("Precipitation",
                                                      "Temperature",
                                                      "Humidity",
                                                      "Cloud.Cover",
                                                      "Soil.pH",
                                                      "Soil.Moisture",
                                                      "Soil.Temperature",
                                                      "Organic.Matter",
                                                      "Nitrate.ppm",
                                                      "Potassium.ppm",
                                                      "Sulfate.ppm",
                                                      "Calcium.ppm",
                                                      "Magnesium.ppm", 
                                                      "Sodium.ppm"),
                                        classification = "Phylum", 
                                        treatment = "tissue",
                                        method = 'spearman', cores = 4)

bact_adj <- bact_corr_heat
bact_adj$p <- p.adjust(bact_adj$p, method = "BH")

bact_adj <- filter(bact_adj, p < .05)

corr_df_b <- bact_adj
gb <- ggplot(corr_df_b, aes(x = Y, y = X, fill = rho)) + 
  geom_tile(color="black") + facet_wrap(~Treatment, nrow = 4) +
  scale_fill_gradient2(low = "red", mid = "white",high = "blue") + theme_bw()  +
  ggtitle(("Effect of Environment on Bacteria Phylum")) + 
  theme(panel.grid.major.x = element_blank()) + 
  xlab(element_blank()) + ylab(element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) + scale_y_discrete(limits=rev)

gb

the_heatmap <- ggarrange(gb,gf, ncol = 2, common.legend = TRUE, legend = "right",
                         labels = c("A","B"))


ggsave("Environmental_Corr_Heatmap.png", plot = the_heatmap, path = "Results_Figs_Tables/Quick_Figures", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

