# Sankey Flow Figure
# Bacteria only 

options(scipen=999)
set.seed(18)

library(ggplot2)
library(phyloseq)
library(plotrix)
library(vegan)
library(tidyverse)
library(devtools)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(gridExtra)
library(tidyr)
#library(microbiome)
library(ape)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(ggpubr)
library(rbiom)
devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(ggthemes)
library(RColorBrewer)
library(pals)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")

# Read in the data

# Read in the data

metadata <- read.csv("QZA_Files/MetadataKey5.txt", header = TRUE, sep = "\t")

load("phy_filtered_bacteria.rdata")

phy_filtered_bacteria
# have to remove the sloth fur sample
phy_filtered = subset_samples(phy_filtered_bacteria, tissue != "Sloth fur")

### make data frames from phyloseq object
f_otu = as.data.frame(otu_table(phy_filtered,taxa_are_rows = TRUE))
f_taxa = as.data.frame(tax_table(phy_filtered))

### sum based on tissue
soilsum <- as.data.frame(taxa_sums(subset_samples(phy_filtered, tissue=="Soil")))
rootwashsum <- as.data.frame(taxa_sums(subset_samples(phy_filtered, tissue=="Root wash")))
rootsum <- as.data.frame(taxa_sums(subset_samples(phy_filtered, tissue=="Root")))
leafsum <- as.data.frame(taxa_sums(subset_samples(phy_filtered, tissue=="Leaf")))

tiss_tax <- do.call("cbind", list(soilsum,rootwashsum,rootsum,leafsum))
colnames(tiss_tax) <- c("Soil","Rhizosphere","Root","Leaf")

tiss_tax2 <- cbind(tiss_tax,f_taxa$Phylum)
colnames(tiss_tax2) <- c("Soil","Rhizosphere","Root","Leaf","Phylum")
tiss_tax2$Soil <- ifelse(tiss_tax2$Soil > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2$Rhizosphere <- ifelse(tiss_tax2$Rhizosphere > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2$Root <- ifelse(tiss_tax2$Root > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2$Leaf <- ifelse(tiss_tax2$Leaf > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2

taxa_long <- tiss_tax2 %>% make_long(Soil,Rhizosphere,Root,Leaf)
taxa_filt <- filter(taxa_long, node != 0)

# change to rhizosphere
#taxa_filt[taxa_filt == "RootWash"] <- as.factor("Rhizosphere")


unique(taxa_filt$node)
my_labels <- c("Verrucomicrobia","Proteobacteria",
               "Chloroflexi","Actinobacteriota","Acidobacteriota",
               "Bacteroidota", "Fusobacteriota", "Firmicutes","Methylomirabilota","Myxococcota",
               "Planctomycetota", "Cyanobacteria","Gemmatimonadota",
               "Nitrospirota","Verrucomicrobiota")

taxa_filt$labels <- ifelse((taxa_filt$node %in% my_labels & taxa_filt$x=='Soil'), taxa_filt$node, NA)


# Fig with geom_sankey
gs <- ggplot(taxa_filt, aes(x = x
                            , next_x = next_x
                            , node = node
                            , next_node = next_node
                            , fill = (factor(node))
                            , label = labels)
) + geom_sankey(flow.alpha = 0.5
                , node.color = "black"
                ,show.legend = TRUE) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 18)) + theme(axis.title.x = element_blank()
                                                       , axis.title.y = element_blank()
                                                       , axis.text.y = element_blank()
                                                       , axis.ticks = element_blank()  
                                                       , panel.grid = element_blank(),
                                                       panel.background = element_blank()) + 
  guides(fill="none") + geom_sankey_label(size = 4,
                                          color = "black",
                                          fill= "white") + ylab("Flow of ASVs") + scale_fill_manual(values=as.vector(polychrome()))


phy_filtered
phy.c <- transform_sample_counts(phy_filtered, function(x) x/sum(x))
phy.melt <- psmelt(phy.c)

phy.melt[phy.melt == 'Root wash'] <- 'Rhizosphere'

relbar <- ggplot(data = phy.melt, mapping = aes_string(x = "tissue",y = "Abundance")) +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill") 

relbar <- relbar + theme(legend.position = "none",
                         axis.title.x = element_blank(),
                         axis.title.y = element_text(size = 16,face = "bold"),
                         axis.text.x = element_text(size = 18, angle = 90, vjust = .7),
                         axis.text.y = element_text(size = 12),
                         panel.background = element_rect(fill = "white",color="white"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_blank() ) + ylab("Relative Abundance of Reads") +
  scale_fill_manual(values=as.vector(polychrome())) + 
  scale_color_manual(values=as.vector(polychrome())) + scale_y_reverse() + 
  scale_x_discrete(limits = c("Soil", "Rhizosphere", "Root","Leaf"))

relbar

flow_and_relabund_bs <- ggarrange(gs,relbar, widths = c(1,.5))


# Fungus


load("phy_filtered.rdata")

phy_filtered
# have to remove the sloth fur sample
phy_filtered = subset_samples(phy_filtered, tissue != "Sloth fur")


### make data frames from phyloseq object
f_otu = as.data.frame(otu_table(phy_filtered,taxa_are_rows = TRUE))
f_taxa = as.data.frame(tax_table(phy_filtered))

### sum based on tissue
soilsum <- as.data.frame(taxa_sums(subset_samples(phy_filtered, tissue=="Soil")))
rootwashsum <- as.data.frame(taxa_sums(subset_samples(phy_filtered, tissue=="Root wash")))
rootsum <- as.data.frame(taxa_sums(subset_samples(phy_filtered, tissue=="Root")))
leafsum <- as.data.frame(taxa_sums(subset_samples(phy_filtered, tissue=="Leaf")))

tiss_tax <- do.call("cbind", list(soilsum,rootwashsum,rootsum,leafsum))
colnames(tiss_tax) <- c("Soil","Rhizosphere","Root","Leaf")

tiss_tax2 <- cbind(tiss_tax,f_taxa$Phylum)
colnames(tiss_tax2) <- c("Soil","Rhizosphere","Root","Leaf","Phylum")
tiss_tax2$Soil <- ifelse(tiss_tax2$Soil > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2$Rhizosphere <- ifelse(tiss_tax2$Rhizosphere > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2$Root <- ifelse(tiss_tax2$Root > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2$Leaf <- ifelse(tiss_tax2$Leaf > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2

taxa_long <- tiss_tax2 %>% make_long(Soil,Rhizosphere,Root,Leaf)
taxa_filt <- filter(taxa_long, node != 0)

unique(taxa_filt$node)
my_labels <- c("Ascomycota","Basidiomycota","Chytridiomycota","Kickxellomycota","Mucoromycota")

taxa_filt$labels <- ifelse((taxa_filt$node %in% my_labels & taxa_filt$x=='Soil'), taxa_filt$node, NA)


# Fig with geom_sankey
gs <- ggplot(taxa_filt, aes(x = x
                            , next_x = next_x
                            , node = node
                            , next_node = next_node
                            , fill = (factor(node))
                            , label = labels)
) + geom_sankey(flow.alpha = 0.5
                , node.color = "black"
                ,show.legend = TRUE) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 18)) + theme(axis.title.x = element_blank()
                                                       , axis.title.y = element_blank()
                                                       , axis.text.y = element_blank()
                                                       , axis.ticks = element_blank()  
                                                       , panel.grid = element_blank(),
                                                       panel.background = element_blank()) + 
  guides(fill="none") + geom_sankey_label(size = 4,
                                          color = "black",
                                          fill= "white") + ylab("Flow of ASVs") + 
  scale_fill_manual(values=as.vector(polychrome()))



phy_filtered
phy.c <- transform_sample_counts(phy_filtered, function(x) x/sum(x))
phy.melt <- psmelt(phy.c)
phy.melt[phy.melt == 'Root wash'] <- 'Rhizosphere'

relbar <- ggplot(data = phy.melt, mapping = aes_string(x = "tissue",y = "Abundance")) +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill") 

relbar <- relbar + theme(legend.position = "none",
                         axis.title.x = element_blank(),
                         axis.title.y = element_text(size = 16,face = "bold"),
                         axis.text.x = element_text(size = 18, angle = 90, vjust = .7),
                         axis.text.y = element_text(size = 12),
                         panel.background = element_rect(fill = "white",color="white"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_blank() ) + ylab("Relative Abundance of Reads") +
  scale_fill_manual(values=as.vector(polychrome())) + 
  scale_color_manual(values=as.vector(polychrome()))  + scale_y_reverse()  + 
  scale_x_discrete(limits = c("Soil", "Rhizosphere", "Root","Leaf"))
relbar



flow_and_relabund_fung <- ggarrange(gs,relbar, widths = c(1,.5))



Figure_1 <- ggarrange(flow_and_relabund_bs,flow_and_relabund_fung, ncol = 2, labels = c("A","B"))

ggsave("Fig1_SankeyCombined.png", plot = Figure_1, path = "Results_Figs_Tables/Figures_Draft1", dpi = 700, 
       width = 16, height = 8, units = c("in"), device = "png")

