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
phylum_data <- tax_glom(phy_clean, "Phylum")    # 48
class_data <- tax_glom(phy_clean, "Class")      # 143
order_data <- tax_glom(phy_clean, "Order")      # 323
family_data <- tax_glom(phy_clean, "Family")    # 722
genus_data <- tax_glom(phy_clean, "Genus")      # 1516
asv_data <- phy_clean                           # 6578

new_meta <- as.data.frame(sample_data(phy_clean))

our_locations <- unique(new_meta$Location)
our_tissues <- unique(new_meta$tissue)


# combined_net <- netConstruct(family_data,
#                              measure = "cclasso",
#                              normMethod = "TSS",
#                              zeroMethod = "pseudo",
#                              thresh = 0.3,
#                              verbose = 2)
# 
# combined_anal <- netAnalyze(combined_net, clustMethod = "cluster_fast_greedy")
# 
# sum_net <- summary(combined_anal)
# globalstats <- as.data.frame(sum_net[1])
# globalstats[6,]

# Natural connectivity tissue - location
# cclasso - super low connectivity, the anova still stands though with soil higher than leaf
# spearman - high connectivity
# phylum  - 0.02307 | 0.77475 - 0.66241 - 41.7 seconds
# class   - 0.01986 | 0.80539 - 0.7142  - 8.895147 mins
# order   - 0.06356 | 0.80272 - 0.74099 - 31 mins
# family  - 0.04915 | .82,.81 2 for 12 mins
# genus - Inf         Inf - 2 for 1.2 hours. 

# Natural Connectivity: 
# https://iopscience.iop.org/article/10.1088/0256-307X/27/7/078902#:~:text=The%20natural%20connectivity%20has%20a,addition%20or%20deletion%20of%20edges.


#Empty DF

Connectivity_df <- data.frame(matrix(ncol=4,nrow=0,
                                     dimnames=list(NULL, c("Location", "Tissue", "Connectivity","nTaxa"))))

# loop for all locations
start_time <- Sys.time()

for(l in 1:length(our_locations)){ 
  print(our_locations[l])
  
  for(t in 1:length(our_tissues)){
   loc <- our_locations[l]
   tiss <- our_tissues[t]
   
   location_phy <- subset_samples(phylum_data, Location == loc & 
                                    tissue == tiss)
   location_phy <- prune_taxa(taxa_sums(location_phy) > 0, location_phy) # drop zeros!
  
   sums <- taxa_sums(location_phy) > 0
   n <- length(sums[sums==TRUE])

   location_net <- netConstruct(location_phy,
                                 measure = "spearman",
                                 normMethod = "clr",
                                 sparsMethod = "threshold",
                                 thresh = 0.5,
                                 verbose = 0)
   
    location_anal <- netAnalyze(location_net, clustMethod = "cluster_fast_greedy",
                                graphlet = FALSE,normDeg = FALSE, weightDeg = TRUE)
  
   sum_net <- location_anal$globalProps
   globalstats <- as.data.frame(sum_net[15])
   
   Connectivity_df <- rbind(Connectivity_df,
                            data.frame(Location = loc, Tissue = tiss, 
                                       Connectivity = globalstats, nTaxa = n ))
  }
  print("________________________________")
}
end_time <- Sys.time()
print(end_time - start_time)


names(connect_phy)[3] ="Connectivity"

# drop sample specific columns and clean up to get unique rows for each 
new_meta1 <- as.data.frame(new_meta[ , -which(names(new_meta) %in% 
                                  c("tissue","count","crude_16s","crude_its","Microbe",
                                    "collection_date", "latitude","longitude","Loc.YR",
                                    "Soil.EC","UV.Light", "Co2.ppm","Texture"))])
rownames(new_meta1) <- NULL
new_meta2 <- as.data.frame(as.matrix(new_meta1))
distinct_meta <- new_meta2 %>% distinct()

mlm_df <- merge(connect_phy,distinct_meta, by = "Location") ### change taxa level df here

# center and scale 
mlm_df[,8:50] <- sapply(mlm_df[,8:50], as.numeric)
mlm_df[,8:50] <- scale(mlm_df[,8:50])

################## statistics
# library(broom)
library(car)
# library(data.table)
# library(lme4)
# library(lmerTest)

# will need to center and scale the covariates

# Linear Models ####################################################
# Tissue and Location are fixed
linmod <- lm(Connectivity ~ Tissue + Location + nTaxa, data = mlm_df)
Anova(linmod, type = c("II"))

linmod <- lm(Connectivity ~ nTaxa + Location + Tissue, data = mlm_df)
Anova(linmod, type = c("II"))

## do something similar to endophyte ANOVA type II!!!

plot_conn_tissue <- ggplot(mlm_df, aes(x = Tissue, y = Connectivity)) + 
  geom_violin(aes(fill = Tissue)) + geom_point()

plot_conn_tissue

plot_conn_location <- ggplot(mlm_df, aes(x = Location, y = Connectivity)) + 
  geom_violin(aes(fill = Location)) + geom_point()

plot_conn_location


# We see massive differences in tissue due to differences in # of taxa

## do something similar to endophyte ANOVA type II!!!
############ break it up and see what has the most effect in each compartment?

Soil <- filter(mlm_df, Tissue == "Soil")
Root <- filter(mlm_df, Tissue == "Root")
RootWash <- filter(mlm_df, Tissue == "Root wash")
Leaf <- filter(mlm_df, Tissue == "Leaf")



Environmental_Vars <- c("temp","precip", "Organic.Matter",
                        "Nitrate.N.ppm.N", "Potassium.ppm.K",
                        "Sulfate.S.ppm.S", "Calcium.ppm.Ca", 
                        "Magnesium.ppm.Mg", "Sodium.ppm.Na",
                        "tempmax", "tempmin", 
                        "temp", "dew", 
                        "humidity", "snow", 
                        "windgust","windspeed",
                        "sealevelpressure", "cloudcover", 
                        "solarradiation", "solarenergy",
                        "Soil.Temp", "Soil.Moisture")

# Soil
Signif_Soil <- data.frame(matrix(ncol=3,nrow=0,))
for (e in 1:length(Environmental_Vars)) {
  formula <- paste("Connectivity ~ nTaxa +",Environmental_Vars[[e]], sep = "")
  linmod <- lm(formula, data = Soil)
  anvtab <- as.data.frame(Anova(linmod, type = c("II")))
  print(anvtab)
  Signif_Soil <- rbind(Signif_Soil, as.data.frame(anvtab[2,]))
}
# Root Wash
Signif_RootWash <- data.frame(matrix(ncol=3,nrow=0,))
for (e in 1:length(Environmental_Vars)) {
  formula <- paste("Connectivity ~ nTaxa +",Environmental_Vars[[e]], sep = "")
  linmod <- lm(formula, data = RootWash)
  anvtab <- as.data.frame(Anova(linmod, type = c("II")))
  #print(anvtab)
  Signif_RootWash <- rbind(Signif_RootWash, as.data.frame(anvtab[2,]))
}
# Root
Signif_Root <- data.frame(matrix(ncol=3,nrow=0,))
for (e in 1:length(Environmental_Vars)) {
  formula <- paste("Connectivity ~ nTaxa +",Environmental_Vars[[e]], sep = "")
  linmod <- lm(formula, data = Root)
  anvtab <- as.data.frame(Anova(linmod, type = c("II")))
  Signif_Root <- rbind(Signif_Root, as.data.frame(anvtab[2,]))
}
# Leaf
Signif_Leaf <- data.frame(matrix(ncol=3,nrow=0,))
for (e in 1:length(Environmental_Vars)) {
  formula <- paste("Connectivity ~ nTaxa +",Environmental_Vars[[e]], sep = "")
  linmod <- lm(formula, data = RootWash)
  anvtab <- as.data.frame(Anova(linmod, type = c("II")))
  Signif_Leaf <- rbind(Signif_Leaf, as.data.frame(anvtab[2,]))
}


# cant use Location as you run out of dfs


# get a correlation of phylum and family level ###########################



#Empty DF

Connectivity_df <- data.frame(matrix(ncol=4,nrow=0,
                                     dimnames=list(NULL, c("Location", "Tissue", "Connectivity","nTaxa"))))

# loop for all locations
start_time <- Sys.time()

for(l in 1:length(our_locations)){ 
  print(our_locations[l])
  
  for(t in 1:length(our_tissues)){
    loc <- our_locations[l]
    tiss <- our_tissues[t]
    
    location_phy <- subset_samples(phylum_data, Location == loc & 
                                     tissue == tiss)
    location_phy <- prune_taxa(taxa_sums(location_phy) > 0, location_phy) # drop zeros!
    
    sums <- taxa_sums(location_phy) > 0
    n <- length(sums[sums==TRUE])
    
    location_net <- netConstruct(location_phy,
                                 measure = "spearman",
                                 normMethod = "clr",
                                 sparsMethod = "threshold",
                                 thresh = 0.5,
                                 verbose = 0)
    
    location_anal <- netAnalyze(location_net, clustMethod = "cluster_fast_greedy",
                                graphlet = FALSE,normDeg = FALSE, weightDeg = TRUE)
    
    sum_net <- location_anal$globalProps
    globalstats <- as.data.frame(sum_net[15])
    
    Connectivity_df <- rbind(Connectivity_df,
                             data.frame(Location = loc, Tissue = tiss, 
                                        Connectivity = globalstats, nTaxa = n ))
  }
  print("________________________________")
}
end_time <- Sys.time()
print(end_time - start_time)

connect_phy <- Connectivity_df




#Empty DF

Connectivity_df <- data.frame(matrix(ncol=4,nrow=0,
                                     dimnames=list(NULL, c("Location", "Tissue", "Connectivity","nTaxa"))))

# loop for all locations
start_time <- Sys.time()

for(l in 1:length(our_locations)){ 
  print(our_locations[l])
  
  for(t in 1:length(our_tissues)){
    loc <- our_locations[l]
    tiss <- our_tissues[t]
    
    location_phy <- subset_samples(family_data, Location == loc & 
                                     tissue == tiss)
    location_phy <- prune_taxa(taxa_sums(location_phy) > 0, location_phy) # drop zeros!
    
    sums <- taxa_sums(location_phy) > 0
    n <- length(sums[sums==TRUE])
    
    location_net <- netConstruct(location_phy,
                                 measure = "spearman",
                                 normMethod = "clr",
                                 sparsMethod = "threshold",
                                 thresh = 0.5,
                                 verbose = 0)
    
    location_anal <- netAnalyze(location_net, clustMethod = "cluster_fast_greedy",
                                graphlet = FALSE,normDeg = FALSE, weightDeg = TRUE)
    
    sum_net <- location_anal$globalProps
    globalstats <- as.data.frame(sum_net[15])
    
    Connectivity_df <- rbind(Connectivity_df,
                             data.frame(Location = loc, Tissue = tiss, 
                                        Connectivity = globalstats, nTaxa = n ))
  }
  print("________________________________")
}
end_time <- Sys.time()
print(end_time - start_time)

connect_family <- Connectivity_df





#Empty DF

Connectivity_df <- data.frame(matrix(ncol=4,nrow=0,
                                     dimnames=list(NULL, c("Location", "Tissue", "Connectivity","nTaxa"))))

# loop for all locations
start_time <- Sys.time()

for(l in 1:length(our_locations)){ 
  print(our_locations[l])
  
  for(t in 1:length(our_tissues)){
    loc <- our_locations[l]
    tiss <- our_tissues[t]
    
    location_phy <- subset_samples(phy_clean, Location == loc & 
                                     tissue == tiss)
    location_phy <- prune_taxa(taxa_sums(location_phy) > 0, location_phy) # drop zeros!
    
    sums <- taxa_sums(location_phy) > 0
    n <- length(sums[sums==TRUE])
    
    location_net <- netConstruct(location_phy,
                                 measure = "spearman",
                                 normMethod = "clr",
                                 sparsMethod = "threshold",
                                 thresh = 0.5,
                                 verbose = 0)
    
    location_anal <- netAnalyze(location_net, clustMethod = "cluster_fast_greedy",
                                graphlet = FALSE,normDeg = FALSE, weightDeg = TRUE)
    
    sum_net <- location_anal$globalProps
    globalstats <- as.data.frame(sum_net[15])
    
    Connectivity_df <- rbind(Connectivity_df,
                             data.frame(Location = loc, Tissue = tiss, 
                                        Connectivity = globalstats, nTaxa = n ))
  }
  print("________________________________")
}
end_time <- Sys.time()
print(end_time - start_time)

connect_asv <- Connectivity_df




# check correlation of different taxa levels:
# phylum_conn <- Connectivity_df
# class_conn <- Connectivity_df
# 
cor.test(connect_phy$Connectivity,connect_asv$Connectivity) # .614 p < .000001

###################### Make a Network Connectivity Figure for the paper


plot_conn_tissue <- ggplot(mlm_df, aes(x = Tissue, y = Connectivity)) + 
  geom_violin(aes(fill = Tissue)) + 
  geom_point() +  theme_bw()  +
  ggtitle(("Natural Connectivity across Tissues - Phylum")) + 
  theme(panel.grid.major.x = element_blank()) +
  labs(y = "Natural Connectivity", x = element_blank()) + 
  theme(axis.title.y = element_text(face = "bold",size = 12)) +
  theme(axis.text.x = element_text(face = "bold",size = 12)) +
  theme(axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 16, face = "bold")) + 
  scale_fill_manual(values = c("Soil" = "brown",
                                "Root wash" = "cyan 3",
                                "Root" = "gold3",
                                "Leaf" = "green4")) +
  theme(legend.position = "none")       


plot_conn_tissue

ggsave("NaturalConnectivity_Tissues_Phylum.png", plot = plot_conn_tissue, 
       path = "Results_Figs_Tables/Quick_Figures", dpi = 700, 
       width = 8, height = 8, units = c("in"), device = "png")


write.table(connect_asv,file = "ASV_connectivity.csv", 
            sep = "\t", col.names = TRUE)
write.table(connect_family,file = "Family_connectivity.csv", 
            sep = "\t", col.names = TRUE)
write.table(connect_phy,file = "Phylum_connectivity.csv", 
            sep = "\t", col.names = TRUE)

write.table(Signif_Soil,file = "Soil_Cov_Signif.csv", 
            sep = "\t", col.names = TRUE)
write.table(Signif_RootWash,file = "RootWash_Cov_Signif.csv", 
            sep = "\t", col.names = TRUE)
write.table(Signif_Root,file = "Root_Cov_Signif.csv", 
            sep = "\t", col.names = TRUE)
write.table(Signif_Leaf,file = "Leaf_Cov_Signif.csv", 
            sep = "\t", col.names = TRUE)



