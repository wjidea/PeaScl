#!/bin/bash

for file in /home/peascl/scl_de/STARout_scl_unstr/ab_initio/all_bam/*.bam; do
	cufflinks -p 4 -g ./sclerotinia_sclerotiorum_2_transcripts.gtf -b ./Sclsc1_AssemblyScaffolds.fasta -o ${file%.Aligned*} $file;
done; 
