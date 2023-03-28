# Dowloaded UNITE database from 
# https://doi.plutof.ut.ee/doi/10.15156/BIO/2483915
# Followed this tutorial to scrub it and load below:
# https://forum.qiime2.org/t/fungal-its-analysis-tutorial/7351

source activate qiime2-2022.11

wrkdir=/home/crschul/IndigoAgMicrobiome/

qiime tools import \
 --type FeatureData[Sequence] \
 --input-path $wrkdir/UNITE/developer/sh_refs_qiime_ver9_99_16.10.2022_dev_uppercase.fasta \
 --output-path $wrkdir/UNITE/unite-ver9_99_16.10.2022.qza

 qiime tools import \
 --type FeatureData[Taxonomy] \
 --input-path $wrkdir/UNITE/developer/sh_taxonomy_qiime_ver9_99_16.10.2022_dev.txt \
 --output-path $wrkdir/UNITE/unite-ver9-99-tax-16.10.2022.qza \
 --input-format HeaderlessTSVTaxonomyFormat

qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads $wrkdir/UNITE/unite-ver9_99_16.10.2022.qza \
 --i-reference-taxonomy $wrkdir/UNITE/unite-ver9-99-tax-16.10.2022.qza \
 --o-classifier $wrkdir/UNITE/unite-ver9-99-CLASSIFIER-16.10.2022.qza