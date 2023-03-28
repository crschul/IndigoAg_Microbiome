#! /bin/bash

mkdir /home/crschul/IndigoAgMicrobiome/FastQC_Check1

wrkdir=/home/crschul/IndigoAgMicrobiome/Indigo_RawData/Fungal_Sequences

for File in $(ls $wrkdir/*R*);do

fastqc --extract $File --outdir=/home/crschul/IndigoAgMicrobiome/FastQC_Check1

done


#cd into this directory and run multiqc . to generate a report for EVERY sample
