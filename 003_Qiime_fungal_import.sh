#! /bin/bash
source activate qiime2-2022.11

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /home/crschul/IndigoAgMicrobiome/Indigo_Fungal_Manifest.tsv \
  --output-path /home/crschul/IndigoAgMicrobiome/QZA_Files/Fungal_seqs.qza\
  --input-format PairedEndFastqManifestPhred33V2