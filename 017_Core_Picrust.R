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
library(microbiome)

#turn off scientific notation
options(scipen=999)

setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome")


#load("phy_combined.rdata")
metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")
metadata <- filter(metadata, tissue != "Rhizosphere")


############### make a phyloseq object with bacteria and fungal ASVs!
# 
# 
# # fungus
# SVs <- read_qza("QZA_Files/Fungal_dada_table.qza")
# 
# taxonomy <- read_qza("QZA_Files/Fungal_taxonomy.qza")
# taxtable <-taxonomy$data %>% as.data.frame() 
# taxtable <- taxtable %>% separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") #convert the table into a tabular split version
# 
# phy <- qza_to_phyloseq(features="QZA_Files/Fungal_dada_table.qza", taxonomy="QZA_Files/Fungal_taxonomy.qza", metadata="Metadata_Indigo_Clean.tsv")
# phy
# 
# 
# 
# ##### FILTERING
# 
# ##    Taxonomic Filtering: Supervised
# # create taxonomy table: number of features for each phyla
# table(tax_table(phy)[, "Phylum"], exclude = NULL)
# 
# # Remove phyla for which only one feature was observed. And ambiguous calls
# 
# 
# #Filter out Chlorophlasts and Mitochondria
# phy
# psm = subset_taxa(phy,Family!="Mitochondria")
# psm
# psu = subset_taxa(psm, Phylum != "unidentified")
# psu
# psc = subset_taxa(psu,Order!="Chloroplast")
# psc
# 
# # remove problem sample
# psc = subset_samples(psc, sample_names(psc) !="sample5324")
# 
# phy_f_genus <- tax_glom(psc, taxrank = "Genus", NArm = FALSE)
# 
# 
# 
# #### Compute filtering options
# 
# # Filter any ASV that is in less than 2 samples and has less than 2 reads
# 
# filter_compare_reads <- function(phy_obj){
#   
#   prevdf = apply(X = otu_table(phy_obj),
#                  MARGIN = ifelse(taxa_are_rows(phy_obj), yes = 1, no = 0),
#                  FUN = function(x){sum(x > 0)})
#   
#   # Add taxonomy and total read counts to this data.frame
#   prevdf = data.frame(Prevalence = prevdf,
#                       TotalAbundance = taxa_sums(phy_obj),
#                       tax_table(phy_obj))
#   
#   prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(phy_obj, "Phylum"))
#   
#   # Total taxa sums
#   TotalReads <- sum(prevdf1$TotalAbundance)
#   readThreshold <- 2
#   
#   prevalenceThreshold = 2 #We have such a wide variety of locations and tissues that we don't want to lose alot 0.01 * nsamples(psc)
#   prevalenceThreshold
#   
#   # taxa must be in at least two samples and make up .0001% of total reads
#   keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold & prevdf1$TotalAbundance >= readThreshold)]
#   phy_filt = prune_taxa(keepTaxa, phy_obj)
#   print(sum(sample_sums(phy_filt)))
#   return(phy_filt)
#   
# }
# 
# # Filter any taxa that is in less than 2 samples and has less than 14 reads
# phy_asv <- filter_compare_reads(psc) # 10,044,327
# 
# phy_filtered_fungus_asv <- prune_samples(sample_sums(phy_asv) >= 500, phy_f_genus)
# 
# 
# 
# 
# ### bacteria
# 
# SVs <- read_qza("QZA_Files/table-sequence_runs_16s_2017.qza")
# 
# taxonomy <- read_qza("QZA_Files/taxonomy.qza")
# btreey <- read_qza("QZA_Files/tree.qza")
# taxtable <-taxonomy$data %>% as.data.frame() 
# taxtable <- taxtable %>% separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") #convert the table into a tabular split version
# 
# phy <- qza_to_phyloseq(features="QZA_Files/table-sequence_runs_16s_2017.qza", taxonomy="QZA_Files/taxonomy.qza", tree="QZA_Files/tree.qza", metadata="Metadata_Indigo_Clean.tsv")
# phy
# 
# 
# 
# ##### FILTERING
# 
# ##    Taxonomic Filtering: Supervised
# # create taxonomy table: number of features for each phyla
# table(tax_table(phy)[, "Phylum"], exclude = NULL)
# 
# # Remove phyla for which only one feature was observed. And ambiguous calls
# 
# 
# #Filter out Chlorophlasts and Mitochondria
# phy
# psm = subset_taxa(phy,Family!="Mitochondria")
# psm
# psu = subset_taxa(psm, Phylum != "unidentified")
# psu
# psc = subset_taxa(psu,Order!="Chloroplast")
# psc
# 
# # remove problem sample
# #psc = subset_samples(psc, sample_names(psc) !="sample5324")
# 
# #phy_f_genus <- tax_glom(psc, taxrank = "Genus", NArm = FALSE)
# 
# 
# 
# #### Compute filtering options
# 
# # Filter any Genus that is in less than 1 samples and has less than 2 reads
# 
# filter_compare_reads <- function(phy_obj){
#   
#   prevdf = apply(X = otu_table(phy_obj),
#                  MARGIN = ifelse(taxa_are_rows(phy_obj), yes = 1, no = 0),
#                  FUN = function(x){sum(x > 0)})
#   
#   # Add taxonomy and total read counts to this data.frame
#   prevdf = data.frame(Prevalence = prevdf,
#                       TotalAbundance = taxa_sums(phy_obj),
#                       tax_table(phy_obj))
#   
#   prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(phy_obj, "Phylum"))
#   
#   # Total taxa sums
#   TotalReads <- sum(prevdf1$TotalAbundance)
#   readThreshold <- 2
#   
#   prevalenceThreshold = 2 #We have such a wide variety of locations and tissues that we don't want to lose alot 0.01 * nsamples(psc)
#   prevalenceThreshold
#   
#   # taxa must be in at least two samples and make up .0001% of total reads
#   keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold & prevdf1$TotalAbundance >= readThreshold)]
#   phy_filt = prune_taxa(keepTaxa, phy_obj)
#   print(sum(sample_sums(phy_filt)))
#   return(phy_filt)
#   
# }
# 
# 
# 
# # Filter any taxa that is in less than 2 samples and has less than 2 reads
# phy_asv <- filter_compare_reads(psc) # 10,044,327
# 
# phy_filtered_bacteria <- prune_samples(sample_sums(phy_asv) >= 500, phy_asv)
# 
# 
# 
# ################# Combine them!
# 
# # have to remove the sloth fur sample
# phy_filtered_b = subset_samples(phy_filtered_bacteria, tissue != "Sloth fur")
# 
# phy_filtered_g = subset_samples(phy_filtered_fungus_asv, tissue != "Sloth fur")
# 
# combined_otu <- merge_phyloseq(otu_table(phy_filtered_g), otu_table(phy_filtered_b))
# combined_taxa <- merge_phyloseq(tax_table(phy_filtered_g), tax_table(phy_filtered_b))
# combined_meta <- merge_phyloseq(sample_data(phy_filtered_g), sample_data(phy_filtered_b))
# 
# phy_combined_asv <- phyloseq(combined_otu,combined_taxa,combined_meta)
# phy_combined_asv
# 
# save(phy_combined_asv, file = "phy_combined_asv.rdata")
# 
# setwd("//wsl.localhost/Ubuntu-18.04/home/crschul/IndigoAgMicrobiome/PICRUST_folder") 
# 
# phyloseq2qiime2<-function(physeq){
#   #take a phyloseq object,check for individual parts, write to files ready for qiime2 upload
#   library(phyloseq)
#   library(biomformat)
#   library(ape)
#   library(Biostrings)
#   library(dada2)
#   if(packageVersion("biomformat") < "1.7") {
#     stop("This will only work with biomformat version > 1.7")
#   }
#   ps_name <-deparse(substitute(physeq))
#   taxa_are_rows_logical<-taxa_are_rows(physeq)
#   #write OTU table to biom file
#   if(is.null(access(physeq,"otu_table"))==FALSE){
#     if(taxa_are_rows_logical==TRUE) {
#       otu<-as(otu_table(physeq),"matrix")
#       otu_biom<-make_biom(data=otu)
#       write_biom(otu_biom,biom_file=paste0(ps_name,"_features-table.biom"))
#       print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
#     } else if (taxa_are_rows_logical==FALSE) {
#       otu<-t(as(otu_table(physeq),"matrix"))
#       otu_biom<-make_biom(data=otu)
#       write_biom(otu_biom,biom_file=paste0(ps_name,"_feature-table.biom"))
#       print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
#     }
#   }
#   #write sample data (metadata) to tsv
#   if(is.null(access(physeq,"sam_data"))==FALSE){
#     write.table(sample_data(physeq),file=paste0(ps_name,"_sample-metadata.txt"),
#                 sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#     print(paste0("Writing sample metadata to ",ps_name,"_sample-metadata.txt"))
#   }
#   #write taxonomy table to qiime2 formatted taxonomy
#   if(is.null(access(physeq,"tax_table"))==FALSE){
#     tax<-as(tax_table(physeq),"matrix")
#     tax_cols <- colnames(tax)
#     tax<-as.data.frame(tax)
#     tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
#     for(co in tax_cols) tax[co]<-NULL
#     write.table(tax, file=paste0(ps_name,"_tax.txt"), quote=FALSE, col.names=FALSE, sep="\t")
#     print(paste0("Writing taxonomy table to ",ps_name,"_tax.txt"))
#   }
#   #write phylogenetic tree to newick formwat
#   if(is.null(access(physeq,"phy_tree"))==FALSE){
#     if(is.rooted(phy_tree(physeq))==TRUE) {
#       ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-rooted.newick"))
#       print(paste0("Writing rooted tree to ",ps_name,"_tree-rooted.newick"))
#     } else if (is.rooted(phy_tree(physeq))==FALSE) {
#       ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-unrooted.newick"))
#       print(paste0("Writing unrooted tree to ",ps_name,"_tree-unrooted.newick"))
#     }
#   }
#   #write representative sequences to fasta format
#   if(is.null(access(physeq,"refseq"))==FALSE){
#     writeXStringSet(refseq(physeq),filepath=paste0(ps_name,"_ref-seqs.fasta"))
#     print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
#   } else if (taxa_are_rows_logical==FALSE && unique(grepl("[^ATCG]",colnames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
#     uniquesToFasta(t(otu), fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
#     print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
#   } else if (taxa_are_rows_logical==TRUE && unique(grepl("[^ATCG]",rownames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
#     uniquesToFasta(otu, fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
#     print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
#   }
# }
# 
# # Export to qiime qza and also a biome file
# phyloseq2qiime2(phy_combined_asv)
# 
# # picrust wont run a combined analysis so export them seperately 
# 
# phyloseq2qiime2(phy_filtered_b)
# phyloseq2qiime2(phy_filtered_g)
# 



# biome file = phyCmbFiltClean_features-table.biom

# fasta table = dna-sequences.fasta from /home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_rep_seqs.qza

# metadata is metadata


###############################################################################################################
# Bacteria
###############################################################################################################


### Command line -> conda activate picrust2
### picrust2_pipeline.py -s dna-sequences.fasta -i phyCmbFiltClean_features-table.biom -o picrust2_out_pipeline -p 4

# Load Raw Descriptors and Agglomerated Pathways - Shared Taxa

picrust_descript <- "PICRUST_folder/picrust_bacteria_out/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv"
KO_table <- read.table(picrust_descript, sep = "\t", header = TRUE)

# Create a phyloseq object out of the otu table and the metadata.
#head(EC_table)[1:10]
# E_C table column names have . instead of -
Meta_KO <- metadata
row.names(Meta_KO) <- Meta_KO$Sample_ID

KO_table$description <- gsub("\\[.*\\]","",as.character(KO_table$description))
KO_table <- KO_table[!duplicated(KO_table$description), ]
dtable <- KO_table#[,-(1:2)]
rownames(dtable) <- KO_table[,2]
dtable <- dtable[,-1:-2]

# Drop the descriptors
descriptions = KO_table[,0:2]

Raw_phy <- phyloseq(otu_table(dtable, taxa_are_rows = TRUE), sample_data(Meta_KO))
Raw_phy # this is our functional phyloseq object

#####
### try to categorize by function similar to picrust1

kegg_brite_map <- read.table("PICRUST_folder/picrust1_KO_BRITE_map.tsv",
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
picrust_descript <- "PICRUST_folder/picrust_bacteria_out/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv"
KO_table <- read.table(picrust_descript, sep = "\t", header = TRUE)
rownames(KO_table) <- KO_table[,1]
KO_table <- KO_table[,-1:-2]

### Run function to categorize all KOs by level 3 in BRITE hierarchy.
table_ko_L3 <- categorize_by_function_l3(KO_table, brite_map_trim)
# test_ko_L3_sorted <- test_ko_L3[rownames(orig_ko_L3), ]
head(table_ko_L3)[1:10]
ko_l3 <- table_ko_L3
# Now do the phyloseq object and diff abundance. 

Aglom_phy <- phyloseq(otu_table(ko_l3, taxa_are_rows = TRUE), sample_data(Meta_KO))
Aglom_phy # this is our functional phyloseq object - call it ec even though its kegg for simplicity


###### Data Sets
Aglom_phy
Raw_phy

# Deseq
library("DESeq2")
library("ggplot2")

# first you need to make just core

Soil_pi <- subset_samples(Aglom_phy, tissue == "Soil")
RW_pi <- subset_samples(Aglom_phy, tissue == "Root wash")
Root_pi <- subset_samples(Aglom_phy, tissue == "Root")
Leaf_pi <- subset_samples(Aglom_phy, tissue == "Leaf")

core_soil <- prune_taxa(rownames(as.data.frame(otu_table(Soil_pi))) %in% 
                          microbiome::core_members(Soil_pi, detection = 0,
                                                   prevalence = .9), Soil_pi)
core_rw <- prune_taxa(rownames(as.data.frame(otu_table(RW_pi))) %in% 
                          microbiome::core_members(RW_pi, detection = 0,
                                                   prevalence = .9), RW_pi)
core_root <- prune_taxa(rownames(as.data.frame(otu_table(Root_pi))) %in% 
                          microbiome::core_members(Root_pi, detection = 0,
                                                   prevalence = .9), Root_pi)
core_leaf <- prune_taxa(rownames(as.data.frame(otu_table(Leaf_pi))) %in% 
                          microbiome::core_members(Leaf_pi, detection = 0,
                                                   prevalence = .9), Leaf_pi)

core_phyloseq <- merge_phyloseq(core_soil,core_rw,core_root,core_leaf)

# filter functional groups that just don't make sense
core_phyloseq <- prune_taxa(!rownames(as.data.frame(otu_table(core_phyloseq))) %in% c("Function Uknown","Meiosis - yeast",
                                                                               "Others","General function prediction only",
                                                                               "Function uknown","Apoptosis"
                                                                               ), core_phyloseq)



alpha = .001

# Deseq of all 4 groups

compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~ tissue)
deseq_obj = DESeq2::DESeq(phydesq, test = "Wald", fitType = "parametric")

res = DESeq2::results(deseq_obj, cooksCutoff = FALSE,contrast = c("tissue","Leaf", "Root"))

vst <- DESeq2::varianceStabilizingTransformation(deseq_obj, blind = FALSE)

# Check distance matrix and PCA/PERMANOVA this gets all the data
DESeq2::plotPCA(vst, intgroup="tissue")

# tplot <- DESeq2::plotPCA(vst, intgroup="tissue")
# lplot <- DESeq2::plotPCA(vst, intgroup="Location")
# pca_df <- tplot$data
# pca_df$Location <- lplot$data$Location
# 
# pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, shape = tissue, color = Location)) + geom_point() +
#   facet_grid(. ~ tissue)
# pca_plot
# 
# 
# Bray_dist_all <- phyloseq::distance(core_phyloseq,"bray")
# 
# adonis2(Bray_dist_all ~ sample_data(core_phyloseq)$tissue +
#           sample_data(core_phyloseq)$Location,
#         by = "margin")




dists <- dist(t(assay(vst)))

library("pheatmap")

select <- order(rowMeans(counts(deseq_obj,normalized=TRUE)),
                decreasing=TRUE)[1:50] # need to filter so you can actually read it
 
kegg_stuff <- assay(vst)[select,]

tiss <- as.data.frame(vst$tissue)
rownames(tiss) <- colnames(vst)
colnames(tiss)[1] ="Tissue"

tiss_sorted <- tiss
tiss_sorted$Tissue <- factor(tiss_sorted$Tissue, 
                                   levels=c("Soil","Root wash", "Root", "Leaf"))

tiss_move <- as.data.frame(tiss_sorted[order(tiss_sorted$Tissue),, drop = FALSE])

new_order <- rownames(tiss_move)

kegg_sorted <- kegg_stuff[, new_order]

annot_colors=list(Tissue=c(Soil = "brown",
                           `Root wash` = "cyan 3",
                           Root = "gold3",
                           Leaf = "green4"))

heat <- pheatmap(kegg_sorted, cluster_rows=FALSE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=TRUE,annotation_col=tiss_move, 
         main = "Bacteria Predicted Functional Genomics", 
         annotation_colors = annot_colors, legend = FALSE, treeheight_col = 0)

# ggsave("Results_Figs_Tables/Bacteria_picrust_heat.png", plot = last_plot(), device = "png",
#        dpi = 700, width = 10, height = 10, units = c("in"))

#### generate and save tables because there is too many differentially abundant things

Diff_table_func_tissue <- function(phyloseq_obj, cat1, cat2, alpha_value){
  compartment_t = transform_sample_counts(phyloseq_obj, function(OTU) OTU +1)
  phydesq = phyloseq_to_deseq2(compartment_t, design = ~tissue)
  deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
  
  res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("tissue", cat1, cat2))
  alpha = alpha_value
  sigtab = res[which(res$padj < alpha), ]
  Functional_Group <- (c(rownames(sigtab)))
  sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
  
  return((sigtab))
}

Soil_vs_RW <- Diff_table_func_tissue(core_phyloseq,"Soil","Root wash",.01)
Soil_vs_Root <- Diff_table_func_tissue(core_phyloseq,"Soil","Root",.01)
Soil_vs_Leaf <- Diff_table_func_tissue(core_phyloseq,"Soil","Leaf",.01)
RW_vs_Root <- Diff_table_func_tissue(core_phyloseq,"Root wash","Root",.01)
RW_vs_Leaf <- Diff_table_func_tissue(core_phyloseq,"Root wash","Leaf",.01)
Root_vs_Leaf <- Diff_table_func_tissue(core_phyloseq,"Root","Leaf",.01)

write.csv(Soil_vs_RW,"Results_Figs_Tables/Picrust_Tables/Soil_vs_RW_bact.csv")
write.csv(Soil_vs_Root,"Results_Figs_Tables/Picrust_Tables/Soil_vs_Root_bact.csv")
write.csv(Soil_vs_Leaf,"Results_Figs_Tables/Picrust_Tables/Soil_vs_Leaf_bact.csv")
write.csv(RW_vs_Root,"Results_Figs_Tables/Picrust_Tables/RW_vs_Root._bact.csv")
write.csv(RW_vs_Leaf,"Results_Figs_Tables/Picrust_Tables/RW_vs_Leaf._bact.csv")
write.csv(Root_vs_Leaf,"Results_Figs_Tables/Picrust_Tables/Root_vs_Leaf_bact.csv")



# compare 1 vs all!!!
core_meta <- sample_data(core_phyloseq)

sample_data(core_phyloseq)$Soil_Label <- core_meta$tissue
sample_data(core_phyloseq)$RW_Label <- core_meta$tissue
sample_data(core_phyloseq)$Root_Label <- core_meta$tissue
sample_data(core_phyloseq)$Leaf_Label <- core_meta$tissue
sample_data(core_phyloseq)$In_or_Out <- core_meta$tissue

sample_data(core_phyloseq)$Soil_Label[sample_data(core_phyloseq)$Soil_Label != "Soil"] <- "Other"
sample_data(core_phyloseq)$RW_Label[sample_data(core_phyloseq)$RW_Label != "Root wash"] <- "Other"
sample_data(core_phyloseq)$Root_Label[sample_data(core_phyloseq)$Root_Label != "Root"] <- "Other"
sample_data(core_phyloseq)$Leaf_Label[sample_data(core_phyloseq)$Leaf_Label != "Leaf"] <- "Other"

sample_data(core_phyloseq)$In_or_Out[sample_data(core_phyloseq)$In_or_Out == "Soil"] <- "Outside"
sample_data(core_phyloseq)$In_or_Out[sample_data(core_phyloseq)$In_or_Out == "Root wash"] <- "Outside"
sample_data(core_phyloseq)$In_or_Out[sample_data(core_phyloseq)$In_or_Out == "Root"] <- "Inside"
sample_data(core_phyloseq)$In_or_Out[sample_data(core_phyloseq)$In_or_Out == "Leaf"] <- "Inside"

core_meta <- sample_data(core_phyloseq)


# Soil
compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~Soil_Label)
deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("Soil_Label", "Soil", "Other"))
alpha = alpha
sigtab = res[which(res$padj < alpha), ]
Functional_Group <- (c(rownames(sigtab)))
sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)

soil_v_all <- sigtab
write.csv(soil_v_all,"Results_Figs_Tables/Picrust_Tables/Soil_v_all_bact.csv")


# Root wash
compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~RW_Label)
deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("RW_Label", "Root wash", "Other"))
alpha = alpha
sigtab = res[which(res$padj < alpha), ]
Functional_Group <- (c(rownames(sigtab)))
sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)

rw_v_all <- sigtab
write.csv(rw_v_all,"Results_Figs_Tables/Picrust_Tables/RW_v_all_bact.csv")


# Root 
compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~Root_Label)
deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("Root_Label", "Root", "Other"))
alpha = alpha
sigtab = res[which(res$padj < alpha), ]
Functional_Group <- (c(rownames(sigtab)))
sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)

root_v_all <- sigtab
write.csv(rw_v_all,"Results_Figs_Tables/Picrust_Tables/Root_v_all_bact.csv")


# Leaf 
compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~Leaf_Label)
deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("Leaf_Label", "Leaf", "Other"))
alpha = alpha
sigtab = res[which(res$padj < alpha), ]
Functional_Group <- (c(rownames(sigtab)))
sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)

leaf_v_all <- sigtab
write.csv(rw_v_all,"Results_Figs_Tables/Picrust_Tables/Leaf_v_all_bact.csv")




# In or Out
compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~In_or_Out)
deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("In_or_Out", "Inside", "Outside"))
alpha = alpha
sigtab = res[which(res$padj < alpha), ]
Functional_Group <- (c(rownames(sigtab)))
sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)

In_v_Out <- sigtab
write.csv(In_v_Out,"Results_Figs_Tables/Picrust_Tables/In_v_out_bact.csv")











###############################################################################################################
# Fungi
###############################################################################################################
metadata <- read.csv("Metadata_Indigo_Clean.tsv", header = TRUE, sep = "\t")
metadata <- filter(metadata, tissue != "Rhizosphere")


### Command line -> conda activate picrust2
### picrust2_pipeline.py -s dna-sequences.fasta -i phyCmbFiltClean_features-table.biom -o picrust2_out_pipeline -p 4

# Load Raw Descriptors and Agglomerated Pathways - Shared Taxa

picrust_fungi <- "PICRUST_folder/picrust_fungi_out/ec_ITS_counts.txt_metagenome_out/pred_metagenome_unstrat_descrip.tsv"
EC_table <- read.table(picrust_fungi, sep = "\t", header = TRUE)

# Create a phyloseq object out of the otu table and the metadata.
#head(EC_table)[1:10]
# E_C table column names have . instead of -
Meta_KO <- metadata
row.names(Meta_KO) <- Meta_KO$Sample_ID

EC_table$description <- gsub("\\[.*\\]","",as.character(EC_table$description))
EC_table <- EC_table[!duplicated(EC_table$description), ]
dtable <- EC_table#[,-(1:2)]
rownames(dtable) <- EC_table[,2]
dtable <- dtable[,-1:-2]

# Drop the descriptors
descriptions = EC_table[,0:2]

Raw_phy <- phyloseq(otu_table(dtable, taxa_are_rows = TRUE), sample_data(Meta_KO))
Raw_phy # this is our functional phyloseq object

##### cant cantegorize EC data like you can KO and KEGG :(

###### Data Sets
Raw_phy

# Deseq
library("DESeq2")
library("ggplot2")

# first you need to make just core

Soil_pi <- subset_samples(Raw_phy, tissue == "Soil")
RW_pi <- subset_samples(Raw_phy, tissue == "Root wash")
Root_pi <- subset_samples(Raw_phy, tissue == "Root")
Leaf_pi <- subset_samples(Raw_phy, tissue == "Leaf")

core_soil <- prune_taxa(rownames(as.data.frame(otu_table(Soil_pi))) %in% 
                          microbiome::core_members(Soil_pi, detection = 0,
                                                   prevalence = .9), Soil_pi)
core_rw <- prune_taxa(rownames(as.data.frame(otu_table(RW_pi))) %in% 
                        microbiome::core_members(RW_pi, detection = 0,
                                                 prevalence = .9), RW_pi)
core_root <- prune_taxa(rownames(as.data.frame(otu_table(Root_pi))) %in% 
                          microbiome::core_members(Root_pi, detection = 0,
                                                   prevalence = .9), Root_pi)
core_leaf <- prune_taxa(rownames(as.data.frame(otu_table(Leaf_pi))) %in% 
                          microbiome::core_members(Leaf_pi, detection = 0,
                                                   prevalence = .9), Leaf_pi)

core_phyloseq <- merge_phyloseq(core_soil,core_rw,core_root,core_leaf)

# filter functional groups that just don't make sense
core_phyloseq <- prune_taxa(!rownames(as.data.frame(otu_table(core_phyloseq))) %in% c("not_found"), core_phyloseq)



alpha = .001

# Deseq of all 4 groups

compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~ tissue)
deseq_obj = DESeq2::DESeq(phydesq, test = "Wald", fitType = "parametric")

res = DESeq2::results(deseq_obj, cooksCutoff = FALSE,contrast = c("tissue","Leaf", "Root"))

vst <- DESeq2::varianceStabilizingTransformation(deseq_obj, blind = FALSE)
DESeq2::plotPCA(vst, intgroup="tissue")

library("pheatmap")
library(plyr)

select <- order(rowMeans(DESeq2::counts(deseq_obj,normalized=TRUE)),
                decreasing=TRUE) # need to filter so you can actually read it

kegg_stuff <- assay(vst)[select,]

tiss <- as.data.frame(vst$tissue)
rownames(tiss) <- colnames(vst)
colnames(tiss)[1] ="Tissue"

tiss_sorted <- tiss
tiss_sorted$Tissue <- factor(tiss_sorted$Tissue, 
                             levels=c("Soil","Root wash", "Root", "Leaf"))

tiss_move <- as.data.frame(tiss_sorted[order(tiss_sorted$Tissue),, drop = FALSE])

tiss_move <- transform(tiss_move,
                       Tissue = plyr::revalue(Tissue,c("Root wash" = "Rhizosphere")))

new_order <- rownames(tiss_move)

kegg_sorted <- kegg_stuff[, new_order]

annot_colors=list(Tissue=c(Soil = "brown",
                           Rhizosphere = "cyan 3",
                           Root = "gold3",
                           Leaf = "green4"))

heat <- pheatmap(kegg_sorted, cluster_rows=FALSE, show_rownames=FALSE, show_colnames = FALSE,
                 cluster_cols=TRUE,annotation_col=tiss_move, 
                 main = "Fungal Predicted Functional Genomics", 
                 annotation_colors = annot_colors,legend = FALSE, treeheight_col = 0)

# ggsave("Results_Figs_Tables/Bacteria_picrust_heat.png", plot = last_plot(), device = "png",
#        dpi = 700, width = 10, height = 10, units = c("in"))

#### generate and save tables because there is too many differentially abundant things

Diff_table_func_tissue <- function(phyloseq_obj, cat1, cat2, alpha_value){
  compartment_t = transform_sample_counts(phyloseq_obj, function(OTU) OTU +1)
  phydesq = phyloseq_to_deseq2(compartment_t, design = ~tissue)
  deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
  
  res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("tissue", cat1, cat2))
  alpha = alpha_value
  sigtab = res[which(res$padj < alpha), ]
  Functional_Group <- (c(rownames(sigtab)))
  sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
  
  return((sigtab))
}

Soil_vs_RW <- Diff_table_func_tissue(core_phyloseq,"Soil","Root wash",.01)
Soil_vs_Root <- Diff_table_func_tissue(core_phyloseq,"Soil","Root",.01)
Soil_vs_Leaf <- Diff_table_func_tissue(core_phyloseq,"Soil","Leaf",.01)
RW_vs_Root <- Diff_table_func_tissue(core_phyloseq,"Root wash","Root",.01)
RW_vs_Leaf <- Diff_table_func_tissue(core_phyloseq,"Root wash","Leaf",.01)
Root_vs_Leaf <- Diff_table_func_tissue(core_phyloseq,"Root","Leaf",.01)

write.csv(Soil_vs_RW,"Results_Figs_Tables/Picrust_Tables/Soil_vs_RW_Fungi.csv")
write.csv(Soil_vs_Root,"Results_Figs_Tables/Picrust_Tables/Soil_vs_Root_Fungi.csv")
write.csv(Soil_vs_Leaf,"Results_Figs_Tables/Picrust_Tables/Soil_vs_Leaf_Fungi.csv")
write.csv(RW_vs_Root,"Results_Figs_Tables/Picrust_Tables/RW_vs_Root_Fungi.csv")
write.csv(RW_vs_Leaf,"Results_Figs_Tables/Picrust_Tables/RW_vs_Leaf_Fungi.csv")
write.csv(Root_vs_Leaf,"Results_Figs_Tables/Picrust_Tables/Root_vs_Leaf_Fungi.csv")







# compare 1 vs all!!!
core_meta <- sample_data(core_phyloseq)

sample_data(core_phyloseq)$Soil_Label <- core_meta$tissue
sample_data(core_phyloseq)$RW_Label <- core_meta$tissue
sample_data(core_phyloseq)$Root_Label <- core_meta$tissue
sample_data(core_phyloseq)$Leaf_Label <- core_meta$tissue
sample_data(core_phyloseq)$In_or_Out <- core_meta$tissue

sample_data(core_phyloseq)$Soil_Label[sample_data(core_phyloseq)$Soil_Label != "Soil"] <- "Other"
sample_data(core_phyloseq)$RW_Label[sample_data(core_phyloseq)$RW_Label != "Root wash"] <- "Other"
sample_data(core_phyloseq)$Root_Label[sample_data(core_phyloseq)$Root_Label != "Root"] <- "Other"
sample_data(core_phyloseq)$Leaf_Label[sample_data(core_phyloseq)$Leaf_Label != "Leaf"] <- "Other"

sample_data(core_phyloseq)$In_or_Out[sample_data(core_phyloseq)$In_or_Out == "Soil"] <- "Outside"
sample_data(core_phyloseq)$In_or_Out[sample_data(core_phyloseq)$In_or_Out == "Root wash"] <- "Outside"
sample_data(core_phyloseq)$In_or_Out[sample_data(core_phyloseq)$In_or_Out == "Root"] <- "Inside"
sample_data(core_phyloseq)$In_or_Out[sample_data(core_phyloseq)$In_or_Out == "Leaf"] <- "Inside"
core_meta <- sample_data(core_phyloseq)


# Soil
compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~Soil_Label)
deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("Soil_Label", "Soil", "Other"))
alpha = alpha
sigtab = res[which(res$padj < alpha), ]
Functional_Group <- (c(rownames(sigtab)))
sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)

soil_v_all <- sigtab
write.csv(soil_v_all,"Results_Figs_Tables/Picrust_Tables/Soil_v_all_fung.csv")


# Root wash
compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~RW_Label)
deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("RW_Label", "Root wash", "Other"))
alpha = alpha
sigtab = res[which(res$padj < alpha), ]
Functional_Group <- (c(rownames(sigtab)))
sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)

rw_v_all <- sigtab
write.csv(rw_v_all,"Results_Figs_Tables/Picrust_Tables/RW_v_all_fung.csv")


# Root 
compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~Root_Label)
deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("Root_Label", "Root", "Other"))
alpha = alpha
sigtab = res[which(res$padj < alpha), ]
Functional_Group <- (c(rownames(sigtab)))
sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)

root_v_all <- sigtab
write.csv(rw_v_all,"Results_Figs_Tables/Picrust_Tables/Root_v_all_fung.csv")


# Leaf 
compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~Leaf_Label)
deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("Leaf_Label", "Leaf", "Other"))
alpha = alpha
sigtab = res[which(res$padj < alpha), ]
Functional_Group <- (c(rownames(sigtab)))
sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)

leaf_v_all <- sigtab
write.csv(rw_v_all,"Results_Figs_Tables/Picrust_Tables/Leaf_v_all_fung.csv")


# In or Out
compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~In_or_Out)
deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("In_or_Out", "Inside", "Outside"))
alpha = alpha
sigtab = res[which(res$padj < alpha), ]
Functional_Group <- (c(rownames(sigtab)))
sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)

In_v_Out <- sigtab
write.csv(In_v_Out,"Results_Figs_Tables/Picrust_Tables/In_v_out_fung.csv")


























###############################################################################################################
# Bacteria - EC RAW
###############################################################################################################


### Command line -> conda activate picrust2
### picrust2_pipeline.py -s dna-sequences.fasta -i phyCmbFiltClean_features-table.biom -o picrust2_out_pipeline -p 4

# Load Raw Descriptors and Agglomerated Pathways - Shared Taxa

picrust_bact_ec <- "PICRUST_folder/picrust_bacteria_out/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv"
EC_table <- read.table(picrust_bact_ec, sep = "\t", header = TRUE)

# Create a phyloseq object out of the otu table and the metadata.
#head(EC_table)[1:10]
# E_C table column names have . instead of -
Meta_KO <- metadata
row.names(Meta_KO) <- Meta_KO$Sample_ID

EC_table$description <- gsub("\\[.*\\]","",as.character(EC_table$description))
EC_table <- EC_table[!duplicated(EC_table$description), ]
dtable <- EC_table#[,-(1:2)]
rownames(dtable) <- EC_table[,2]
dtable <- dtable[,-1:-2]

dtable <- dtable[complete.cases(dtable),]

# Drop the descriptors
descriptions = EC_table[,0:2]

Raw_phy <- phyloseq(otu_table(dtable, taxa_are_rows = TRUE), sample_data(Meta_KO))
Raw_phy # this is our functional phyloseq object

##### cant cantegorize EC data like you can KO and KEGG :(


###### Data Sets

# Deseq
library("DESeq2")
library("ggplot2")

# first you need to make just core

Soil_pi <- subset_samples(Raw_phy, tissue == "Soil")
RW_pi <- subset_samples(Raw_phy, tissue == "Root wash")
Root_pi <- subset_samples(Raw_phy, tissue == "Root")
Leaf_pi <- subset_samples(Raw_phy, tissue == "Leaf")

core_soil <- prune_taxa(rownames(as.data.frame(otu_table(Soil_pi))) %in% 
                          microbiome::core_members(Soil_pi, detection = 0,
                                                   prevalence = .1), Soil_pi)
core_rw <- prune_taxa(rownames(as.data.frame(otu_table(RW_pi))) %in% 
                        microbiome::core_members(RW_pi, detection = 0,
                                                 prevalence = .9), RW_pi)
core_root <- prune_taxa(rownames(as.data.frame(otu_table(Root_pi))) %in% 
                          microbiome::core_members(Root_pi, detection = 0,
                                                   prevalence = .9), Root_pi)
core_leaf <- prune_taxa(rownames(as.data.frame(otu_table(Leaf_pi))) %in% 
                          microbiome::core_members(Leaf_pi, detection = 0,
                                                   prevalence = .9), Leaf_pi)

core_phyloseq <- merge_phyloseq(core_soil,core_rw,core_root,core_leaf)

# filter functional groups that just don't make sense
core_phyloseq <- prune_taxa(!rownames(as.data.frame(otu_table(core_phyloseq))) %in% c("not_found"), core_phyloseq)



alpha = .01

# Deseq of all 4 groups

compartment_t = transform_sample_counts(core_phyloseq, function(OTU) OTU +1)
phydesq = phyloseq_to_deseq2(compartment_t, design = ~ tissue)
deseq_obj = DESeq2::DESeq(phydesq, test = "Wald", fitType = "parametric")

res = DESeq2::results(deseq_obj, cooksCutoff = FALSE,contrast = c("tissue","Leaf", "Root"))

vst <- DESeq2::varianceStabilizingTransformation(deseq_obj, blind = FALSE)
DESeq2::plotPCA(vst, intgroup="tissue")

library("pheatmap")

select <- order(rowMeans(DESeq2::counts(deseq_obj,normalized=TRUE)),
                decreasing=TRUE) # need to filter so you can actually read it

kegg_stuff <- assay(vst)[select,]

tiss <- as.data.frame(vst$tissue)
rownames(tiss) <- colnames(vst)
colnames(tiss)[1] ="Tissue"

tiss_sorted <- tiss
tiss_sorted$Tissue <- factor(tiss_sorted$Tissue, 
                             levels=c("Soil","Root wash", "Root", "Leaf"))

tiss_move <- as.data.frame(tiss_sorted[order(tiss_sorted$Tissue),, drop = FALSE])

new_order <- rownames(tiss_move)

kegg_sorted <- kegg_stuff[, new_order]

annot_colors=list(Tissue=c(Soil = "brown",
                           `Root wash` = "cyan 3",
                           Root = "gold3",
                           Leaf = "green4"))

heat <- pheatmap(kegg_sorted, cluster_rows=FALSE, show_rownames=FALSE, show_colnames = FALSE,
                 cluster_cols=TRUE,annotation_col=tiss_move, 
                 main = "Bacteria Predicted Functional Genomics - EC", 
                 annotation_colors = annot_colors)










