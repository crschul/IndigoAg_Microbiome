# dada 2 trims, denoises, and joins your reads. 
# you can use metadata tabulate to look at how much of each sample is retained:
#ie not chimeras, or trimmed etc. 

source activate qiime2-2022.11

wrkdir=/home/crschul/IndigoAgMicrobiome/QZA_Files

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $wrkdir/Fungal_Trim_2.qza \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250 \
  --o-representative-sequences $wrkdir/Fungal_rep_seqs.qza \
  --o-table $wrkdir/Fungal_dada_table.qza \
  --o-denoising-stats $wrkdir/Fungal_dada_stats.qza
