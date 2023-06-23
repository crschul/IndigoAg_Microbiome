# Script for looking at microbiome networks MLM
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
#library(phylosmith)
#library(DESeq2)
library(NetCoMi)

#turn off scientific notation
options(scipen=999)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")


load("phy_combined.rdata")
metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")


# subset our data for samples that have meta data"
phy_clean <- subset_samples(phy_combined, !is.na(Soil.pH))


# see how dissimilarity changes at different taxa levels
# phylum_data <- tax_glom(phy_clean, "Phylum")    # 48
# class_data <- tax_glom(phy_clean, "Class")      # 143
# order_data <- tax_glom(phy_clean, "Order")      # 323
# family_data <- tax_glom(phy_clean, "Family")    # 722
# genus_data <- tax_glom(phy_clean, "Genus")      # 1516
asv_data <- phy_clean                           # 6578

new_meta <- as.data.frame(sample_data(phy_clean))

our_locations <- unique(new_meta$Location)
our_tissues <- unique(new_meta$tissue)


#### The loop to generate the data is in file 16, it takes forever so lets read it in
# you broke all the files you saved manually fix them at somepoint!
connectivity_asv <- read.csv(file = "Results_Figs_Tables/ASV_connectivity.csv", 
                             header = TRUE, sep = ",")


Soil <- filter(connectivity_asv, Tissue == "Soil")
Root <- filter(connectivity_asv, Tissue == "Root")
RootWash <- filter(connectivity_asv, Tissue == "Root wash")
Leaf <- filter(connectivity_asv, Tissue == "Leaf")

my_comparisons <- list( c("Soil", "Root wash"), c("Soil", "Leaf"), c("Soil", "Root"),
                        c("Root wash", "Leaf"), c("Root wash", "Root"), c("Leaf", "Root"))

######### Significance statistics!
anova <- aov(Connectivity ~ Tissue, data = connectivity_asv)
tukey <- TukeyHSD(anova)
tukey
tukey.table <- as.data.frame(tukey$Tissue)

# Change connectivity
connectivity_asv[connectivity_asv == "Root wash"] <- ("Rhizosphere")


plot_conn_tissue <- ggplot(connectivity_asv, aes(x = factor(Tissue, level=c('Soil', 'Rhizosphere', 'Root','Leaf')),
                                                 y = Connectivity)) + 
  geom_violin(aes(fill = Tissue)) + 
  geom_point() +  theme_bw()  +
  ggtitle(("Natural Connectivity across Tissues - ASV")) + 
  theme(panel.grid.major.x = element_blank()) +
  labs(y = "Natural Connectivity", x = element_blank()) + 
  theme(axis.title.y = element_text(face = "bold",size = 12)) +
  theme(axis.text.x = element_text(face = "bold",size = 12)) +
  theme(axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 16, face = "bold")) + 
  scale_fill_manual(values = c("Soil" = "brown",
                               "Rhizosphere" = "cyan 3",
                               "Root" = "gold3",
                               "Leaf" = "green4")) +
  scale_y_continuous(limits=c(.15,.48)) +
  theme(legend.position = "none") + 
  geom_label(label = "a",x = 1, y = .4, label.size = 0, size = 8, color = "black") +
  geom_label(label = "a",x = 2, y = .32, label.size = 0, size = 8, color = "black") +
  geom_label(label = "b",x = 3, y = .47, label.size = 0, size = 8, color = "black") +
  geom_label(label = "c",x = 4, y = .425, label.size = 0, size = 8, color = "black")

plot_conn_tissue


ggsave("Fig4_Connectivity_D1.png", plot = plot_conn_tissue, 
       path = "Results_Figs_Tables/Figures_Draft1", dpi = 700, 
       width = 8, height = 8, units = c("in"), device = "png")

write.table(tukey.table,file = "Results_Figs_Tables/Tissue_Connect_Tukeypvalues.csv", 
            sep = ",", col.names = TRUE)


