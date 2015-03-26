#!/bin/bash
#this file is to run fastqc

#run paralell 
#-j indicates run as many jobs as possible
#+0 indicates add 0 job to cpu core(s)
#FASTQCDIR='/home/wjidea/Files/Turf_bac_rawseq/20130608_DNASeq_PE/fastqc_files'
#echo $FASTQCDIR
ls *.gz | time parallel -j+0 --eta 'fastqc {} --outdir=/home/sclpea/FastQC/'

# for FILE in ./*.gz; 
#   do fastqc $FILE --outdir=/home/wjidea/Files/Turf_bac_rawseq/20130608_DNASeq_PE/fastqc_files;
# done;
# fastqc can define the outdir as below:
# fastqc s1_* s2_* --outdir=/root/Dropbox/fastqc
