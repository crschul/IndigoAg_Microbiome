#! /bin/bash
# This trims out adapters and trims the ends based on quality of the reads to get rid of crappy sequenceing. 

source activate qiime2-2022.11

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences /home/crschul/IndigoAgMicrobiome/QZA_Files/Fungal_seqs.qza \
  --p-cores 4 \
  --p-adapter-r TCGCATAGTGAATCATCGAATC \
  --p-adapter-f TCCTCCGCTTATTGATATGC \
  --p-quality-cutoff-5end 20 \
  --p-quality-cutoff-3end 20 \
  --o-trimmed-sequences /home/crschul/IndigoAgMicrobiome/QZA_Files/Fungal_Trim_2.qza


