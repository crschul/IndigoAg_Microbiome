# use the UNITE database we just build to assign fungal taxonomy

source activate qiime2-2022.11

wrkdir=/home/crschul/IndigoAgMicrobiome/

qiime feature-classifier classify-sklearn \
  --i-classifier $wrkdir/UNITE/unite-ver9-99-CLASSIFIER-16.10.2022.qza \
  --i-reads $wrkdir/QZA_Files/Fungal_rep_seqs.qza \
  --o-classification  $wrkdir/QZA_Files/Fungal_taxonomy.qza