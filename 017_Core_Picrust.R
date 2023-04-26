# Script for looking at core picrust

# We need a combined (bacteria and fungal) phyloseq obj at the asv level to export to picrust 
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
library(phylosmith)

#turn off scientific notation
options(scipen=999)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")


load("phy_combined.rdata")
metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")


############### make a phyloseq object with bacteria and fungal ASVs!


# fungus
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

# Filter any ASV that is in less than 2 samples and has less than 2 reads

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

# Filter any taxa that is in less than 2 samples and has less than 14 reads
phy_asv <- filter_compare_reads(psc) # 10,044,327

phy_filtered_fungus_asv <- prune_samples(sample_sums(phy_asv) >= 500, phy_f_genus)




### bacteria

SVs <- read_qza("QZA_Files/table-sequence_runs_16s_2017.qza")

taxonomy <- read_qza("QZA_Files/taxonomy.qza")
btreey <- read_qza("QZA_Files/tree.qza")
taxtable <-taxonomy$data %>% as.data.frame() 
taxtable <- taxtable %>% separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") #convert the table into a tabular split version

phy <- qza_to_phyloseq(features="QZA_Files/table-sequence_runs_16s_2017.qza", taxonomy="QZA_Files/taxonomy.qza", tree="QZA_Files/tree.qza", metadata="Metadata_Indigo_Clean.tsv")
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
#psc = subset_samples(psc, sample_names(psc) !="sample5324")

#phy_f_genus <- tax_glom(psc, taxrank = "Genus", NArm = FALSE)



#### Compute filtering options

# Filter any Genus that is in less than 1 samples and has less than 2 reads

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



# Filter any taxa that is in less than 2 samples and has less than 2 reads
phy_asv <- filter_compare_reads(psc) # 10,044,327

phy_filtered_bacteria <- prune_samples(sample_sums(phy_asv) >= 500, phy_asv)



################# Combine them!

# have to remove the sloth fur sample
phy_filtered_b = subset_samples(phy_filtered_bacteria, tissue != "Sloth fur")

phy_filtered_g = subset_samples(phy_filtered_fungus_asv, tissue != "Sloth fur")

combined_otu <- merge_phyloseq(otu_table(phy_filtered_g), otu_table(phy_filtered_b))
combined_taxa <- merge_phyloseq(tax_table(phy_filtered_g), tax_table(phy_filtered_b))
combined_meta <- merge_phyloseq(sample_data(phy_filtered_g), sample_data(phy_filtered_b))

phy_combined_asv <- phyloseq(combined_otu,combined_taxa,combined_meta)
phy_combined_asv

save(phy_combined_asv, file = "phy_combined_asv.rdata")

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome/PICRUST_folder") 

phyloseq2qiime2<-function(physeq){
  #take a phyloseq object,check for individual parts, write to files ready for qiime2 upload
  library(phyloseq)
  library(biomformat)
  library(ape)
  library(Biostrings)
  library(dada2)
  if(packageVersion("biomformat") < "1.7") {
    stop("This will only work with biomformat version > 1.7")
  }
  ps_name <-deparse(substitute(physeq))
  taxa_are_rows_logical<-taxa_are_rows(physeq)
  #write OTU table to biom file
  if(is.null(access(physeq,"otu_table"))==FALSE){
    if(taxa_are_rows_logical==TRUE) {
      otu<-as(otu_table(physeq),"matrix")
      otu_biom<-make_biom(data=otu)
      write_biom(otu_biom,biom_file=paste0(ps_name,"_features-table.biom"))
      print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
    } else if (taxa_are_rows_logical==FALSE) {
      otu<-t(as(otu_table(physeq),"matrix"))
      otu_biom<-make_biom(data=otu)
      write_biom(otu_biom,biom_file=paste0(ps_name,"_feature-table.biom"))
      print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
    }
  }
  #write sample data (metadata) to tsv
  if(is.null(access(physeq,"sam_data"))==FALSE){
    write.table(sample_data(physeq),file=paste0(ps_name,"_sample-metadata.txt"),
                sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    print(paste0("Writing sample metadata to ",ps_name,"_sample-metadata.txt"))
  }
  #write taxonomy table to qiime2 formatted taxonomy
  if(is.null(access(physeq,"tax_table"))==FALSE){
    tax<-as(tax_table(physeq),"matrix")
    tax_cols <- colnames(tax)
    tax<-as.data.frame(tax)
    tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
    for(co in tax_cols) tax[co]<-NULL
    write.table(tax, file=paste0(ps_name,"_tax.txt"), quote=FALSE, col.names=FALSE, sep="\t")
    print(paste0("Writing taxonomy table to ",ps_name,"_tax.txt"))
  }
  #write phylogenetic tree to newick formwat
  if(is.null(access(physeq,"phy_tree"))==FALSE){
    if(is.rooted(phy_tree(physeq))==TRUE) {
      ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-rooted.newick"))
      print(paste0("Writing rooted tree to ",ps_name,"_tree-rooted.newick"))
    } else if (is.rooted(phy_tree(physeq))==FALSE) {
      ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-unrooted.newick"))
      print(paste0("Writing unrooted tree to ",ps_name,"_tree-unrooted.newick"))
    }
  }
  #write representative sequences to fasta format
  if(is.null(access(physeq,"refseq"))==FALSE){
    writeXStringSet(refseq(physeq),filepath=paste0(ps_name,"_ref-seqs.fasta"))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
  } else if (taxa_are_rows_logical==FALSE && unique(grepl("[^ATCG]",colnames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
    uniquesToFasta(t(otu), fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
  } else if (taxa_are_rows_logical==TRUE && unique(grepl("[^ATCG]",rownames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
    uniquesToFasta(otu, fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
  }
}

# Export to qiime qza and also a biome file
phyloseq2qiime2(phy_combined_asv)

# picrust wont run a combined analysis so export them seperately 

phyloseq2qiime2(phy_filtered_b)
phyloseq2qiime2(phy_filtered_g)




# biome file = phyCmbFiltClean_features-table.biom

# fasta table = dna-sequences.fasta from /home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_rep_seqs.qza

# metadata is metadata



### Command line -> conda activate picrust2
### picrust2_pipeline.py -s dna-sequences.fasta -i phyCmbFiltClean_features-table.biom -o picrust2_out_pipeline -p 4

# Load Raw Descriptors and Agglomerated Pathways - Shared Taxa
setwd("D:/Manual_Backup/May_2022/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Picrust2/picrust2_phy2noblank_out/KO_metagenome_out")

KO_table <- read.csv("ph2noblank_KOtable_descript.tsv", header=TRUE, sep="\t", row.names=1)
head(KO_table)[1:10]

# Create a phyloseq object out of the otu table and the metadata. 
#head(EC_table)[1:10]
# E_C table column names have . instead of - 
colnames(KO_table) <- gsub(x = colnames(KO_table), pattern = "\\.", replacement = "-")
Meta_EC <- metadata
row.names(Meta_EC) <- Meta_EC$SampleID

KO_table$description <- gsub("\\[.*\\]","",as.character(KO_table$description))
KO_table <- KO_table[!duplicated(KO_table$description), ]
dtable <- KO_table#[,-(1:2)]
rownames(dtable) <- KO_table[,1]
dtable <- dtable[,-1]

# Drop the descriptors
#descriptions = KO_table[,0:2]

Raw_phy <- phyloseq(otu_table(dtable, taxa_are_rows = TRUE), sample_data(Meta_EC))
Raw_phy # this is our functional phyloseq object

#####
### try to categorize by function similar to picrust1
setwd("D:/Manual_Backup/May_2022/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Picrust2")
#KO_table <- read.csv("KO_pred_metagenome_unstrat_descript.tsv", header=TRUE, sep="\t", row.names=1)
head(KO_table)[1:10]

kegg_brite_map <- read.table("picrust1_KO_BRITE_map.tsv",
                             header=TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, comment.char="", row.names=1)
#head(kegg_brite_map)
dim(kegg_brite_map)

# remove rows in brite map related to human diseases and a few other things
# Look at unique groups
unique(unlist(strsplit(as.character(kegg_brite_map$metadata_KEGG_Pathways), ";")))

brite_map_trim <- kegg_brite_map[!grepl("Human Diseases", kegg_brite_map$metadata_KEGG_Pathways),]
dim(brite_map_trim)



categorize_by_function_l3 <- function(in_ko, kegg_brite_mapping) {
  # Function to create identical output as categorize_by_function.py script,
  # but with R objects instead of BIOM objects in Python.
  # Input KO table is assumed to have rownames as KOs and sample names as columns.
  
  out_pathway <- data.frame(matrix(NA, nrow=0, ncol=(ncol(in_ko) + 1)))
  
  colnames(out_pathway) <- c("pathway", colnames(in_ko))
  
  for(ko in rownames(in_ko)) {
    
    # Skip KO if not in KEGG BRITE mapping df
    # (this occurs with newer KOs that weren't present in PICRUSt1).
    if(! ko %in% rownames(kegg_brite_mapping)) {
      next
    }
    
    pathway_list <- strsplit(kegg_brite_mapping[ko, "metadata_KEGG_Pathways"], "\\|")[[1]]
    
    for(pathway in pathway_list) {
      
      pathway <- strsplit(pathway, ";")[[1]][3]
      
      new_row <- data.frame(matrix(c(NA, as.numeric(in_ko[ko,])), nrow=1, ncol=ncol(out_pathway)))
      colnames(new_row) <- colnames(out_pathway)
      new_row$pathway <- pathway
      out_pathway <- rbind(out_pathway, new_row)
    }
    
  }
  
  out_pathway = data.frame(aggregate(. ~ pathway, data = out_pathway[,-2], FUN=sum))
  
  rownames(out_pathway) <- out_pathway$pathway
  
  out_pathway <- out_pathway[, -which(colnames(out_pathway) == "pathway")]
  
  if(length(which(rowSums(out_pathway) == 0)) > 0) {
    out_pathway <- out_pathway[-which(rowSums(out_pathway) == 0), ]
  }
  
  return(out_pathway)
  
}

### categorize by function similar to picrust1

### Run function to categorize all KOs by level 3 in BRITE hierarchy.
table_ko_L3 <- categorize_by_function_l3(KO_table, brite_map_trim)
# test_ko_L3_sorted <- test_ko_L3[rownames(orig_ko_L3), ]
head(table_ko_L3)[1:10]
ko_l3 <- table_ko_L3
# Now do the phyloseq object and diff abundance. 
colnames(ko_l3) <- gsub(x = colnames(ko_l3), pattern = "\\.", replacement = "-")
#ko_l3 <- cbind(rownames(ko_l3), data.frame(ko_l3, row.names = NULL))
#ko_l3 <- subset(ko_l3, select = -c(description))
descriptions = row.names(ko_l3)
ko_l3 = ko_l3[, -c(1)]

head(ko_l3)
colnames(ko_l3)
Meta_EC <- metadata
row.names(Meta_EC) <- Meta_EC$SampleID
Aglom_phy <- phyloseq(otu_table(ko_l3, taxa_are_rows = TRUE), sample_data(Meta_EC))
Aglom_phy # this is our functional phyloseq object - call it ec even though its kegg for simplicity

gplots::venn(list(metadata=rownames(Meta_EC), featuretable=colnames(ko_l3)))

###### Data Sets
Aglom_phy
Raw_phy


# Deseq
library("DESeq2")
library("ggplot2")

stalks_Ag <- subset_samples(Aglom_phy, Sample_Type_Blanks_differentiated=="Stalk")
rhizos_Ag <- subset_samples(Aglom_phy, Sample_Type_Blanks_differentiated=="Rhizosphere")
roots_Ag <- subset_samples(Aglom_phy, Sample_Type_Blanks_differentiated=="Root")
stalks_Raw <- subset_samples(Raw_phy, Sample_Type_Blanks_differentiated=="Stalk")
rhizos_Raw <- subset_samples(Raw_phy, Sample_Type_Blanks_differentiated=="Rhizosphere")
roots_Raw <- subset_samples(Raw_phy, Sample_Type_Blanks_differentiated=="Root")