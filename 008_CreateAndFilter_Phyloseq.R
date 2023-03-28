# Take pre proccessed qiime objects, filter with phyloseq, and use for the rest of the analysis

library(tidyverse)
library(devtools)
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggplot2)
library(phyloseq)
library(gridExtra)

#turn off scientific notation
options(scipen=999)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")

# Read in the data

metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")



SVs <- read_qza("QZA_Files/Fungal_dada_table.qza")

taxonomy <- read_qza("QZA_Files/Fungal_taxonomy.qza")
taxtable <-taxonomy$data %>% as.data.frame() 
taxtable <- taxtable %>% separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") #convert the table into a tabular split version

phy <- qza_to_phyloseq(features="QZA_Files/Fungal_dada_table.qza", taxonomy="QZA_Files/Fungal_taxonomy.qza", metadata="Metadata_Indigo_Clean.tsv")
phy



##### FILTERING

##    Taxonomic Filtering: Supervised
# create taxonomy table: number of features for each phyla
table(tax_table(phy)[, "Phylum"], exclude = NULL)

# Remove phyla for which only one feature was observed. And ambiguous calls


#Filter out Chlorophlasts and Mitochondria
phy
psm = subset_taxa(phy,Family!="Mitochondria")
psm
psu = subset_taxa(psm, Phylum != "unidentified")
psu
psc = subset_taxa(psu,Order!="Chloroplast")
psc

# remove problem sample
psc = subset_samples(psc, sample_names(psc) !="sample5324")

phy_f_genus <- tax_glom(psc, taxrank = "Genus", NArm = FALSE)



#### Compute filtering options

# Filter any Genus that is in less than 2 samples and has less than 14 reads

filter_compare_reads <- function(phy_obj){
  
  prevdf = apply(X = otu_table(phy_obj),
                 MARGIN = ifelse(taxa_are_rows(phy_obj), yes = 1, no = 0),
                 FUN = function(x){sum(x > 0)})
  
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(phy_obj),
                      tax_table(phy_obj))
  
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(phy_obj, "Phylum"))
  
  # Total taxa sums
  TotalReads <- sum(prevdf1$TotalAbundance)
  readThreshold <- 2
  
  prevalenceThreshold = 2 #We have such a wide variety of locations and tissues that we don't want to lose alot 0.01 * nsamples(psc)
  prevalenceThreshold
  
  # taxa must be in at least two samples and make up .0001% of total reads
  keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold & prevdf1$TotalAbundance >= readThreshold)]
  phy_filt = prune_taxa(keepTaxa, phy_obj)
  print(sum(sample_sums(phy_filt)))
  return(phy_filt)

}

sum(sample_sums(phy_f_genus)) # 14,214,797 to start with

# Filter any taxa that is in less than 2 samples and has less than 14 reads
phy_asv <- filter_compare_reads(psc) # 10,044,327
phy_genus <- filter_compare_reads(phy_f_genus) # 14,210,687


### How many total reads do we have?
sample_sums(phy_genus)

phy_filtered_fungus_asv <- prune_samples(sample_sums(phy_asv) >= 500, phy_f_genus)
phy_filtered_fungus_genus <- prune_samples(sample_sums(phy_genus) >= 500, phy_f_genus)
                              
min(sample_sums(phy_filtered))

# save the phyloseq object so we can load it later!
save(phy_filtered_fungus_asv, file = "phy_filtered_fungi_asv.rdata")
save(phy_filtered_fungus_genus, file = "phy_filtered_fungi_genus.rdata")


