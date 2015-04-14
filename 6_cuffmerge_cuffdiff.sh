#Creating assmblies.txt for cuffmerge in clout_genes folder
ls */transcripts.gtf > assemblies.txt

#Running cuffmerge in clout_genes folder
~/Apps/cufflinks-2.2.1.Linux_x86_64/cuffmerge -g ../ref_scl/sclerotinia_sclerotiorum_2_transcripts.gtf -s ../ref_scl/Sclsc1_AssemblyScaffolds.fasta -p 8 assemblies.txt 

#Running cuffquant
~/Apps/cufflinks-2.2.1.Linux_x86_64/cuffquant -o ./cuffquant_out -p 8 -b ../ref_scl/Sclsc1_AssemblyScaffolds.fasta ./merged_asm/merged.gtf /home/peascl/scl_de/thout_genes_scl/Lifter-Scl-12-1/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Lifter-Scl-12-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Lifter-Scl-24-1/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Lifter-Scl-24-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Lifter-Scl-48-1/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Lifter-Scl-48-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/PI240515-Scl-12-1/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/PI240515-Scl-12-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/PI240515-Scl-24-1/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/PI240515-Scl-24-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/PI240515-Scl-48-1/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/PI240515-Scl-48-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Scl0205-12hr-1/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Scl0205-12hr-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Scl0205-24-1/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Scl0205-24-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Scl0205-48-1/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Scl0205-48-BR/accepted_hits.bam 

#Running cuffdiff 
## OLD cufflinks
cuffdiff -o diff_out -b ../ref_scl/Sclsc1_AssemblyScaffolds.fasta -p 8 -L Lifter-Scl-12,Lifter-Scl-24,Lifter-Scl-48,PI240515-Scl-12,PI240515-Scl-24,PI240515-Scl-48,Scl0205-12,Scl0205-24,Scl0205-48h -u ./merged_asm/merged.gtf /home/peascl/scl_de/thout_genes_scl/Lifter-Scl-12-1/accepted_hits.bam,/home/peascl/scl_de/thout_genes_scl/Lifter-Scl-12-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Lifter-Scl-24-1/accepted_hits.bam,/home/peascl/scl_de/thout_genes_scl/Lifter-Scl-24-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Lifter-Scl-48-1/accepted_hits.bam,/home/peascl/scl_de/thout_genes_scl/Lifter-Scl-48-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/PI240515-Scl-12-1/accepted_hits.bam,/home/peascl/scl_de/thout_genes_scl/PI240515-Scl-12-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/PI240515-Scl-24-1/accepted_hits.bam,/home/peascl/scl_de/thout_genes_scl/PI240515-Scl-24-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/PI240515-Scl-48-1/accepted_hits.bam,/home/peascl/scl_de/thout_genes_scl/PI240515-Scl-48-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Scl0205-12hr-1/accepted_hits.bam,/home/peascl/scl_de/thout_genes_scl/Scl0205-12hr-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Scl0205-24-1/accepted_hits.bam,/home/peascl/scl_de/thout_genes_scl/Scl0205-24-BR/accepted_hits.bam /home/peascl/scl_de/thout_genes_scl/Scl0205-48-1/accepted_hits.bam,/home/peascl/scl_de/thout_genes_scl/Scl0205-48-BR/accepted_hits.bam

##Running cuffdiff
## NEW cufflinks (v2.2.1)
~/Apps/cufflinks-2.2.1.Linux_x86_64/cuffdiff -o diff_out -b ../ref_scl/Sclsc1_AssemblyScaffolds.fasta -p 8 -L Lifter-Scl-12,Lifter-Scl-24,Lifter-Scl-48,PI240515-Scl-12,PI240515-Scl-24,PI240515-Scl-48,Scl0205-12,Scl0205-24,Scl0205-48h -u ./merged_asm/merged.gtf /home/peascl/scl_de/clout_gene/cuffquant_out/Lifter-Scl-12-1/abundances.cxb,/home/peascl/scl_de/clout_gene/cuffquant_out/Lifter-Scl-12-BR/abundances.cxb /home/peascl/scl_de/clout_gene/cuffquant_out/Lifter-Scl-24-1/abundances.cxb,/home/peascl/scl_de/clout_gene/cuffquant_out/Lifter-Scl-24-BR/abundances.cxb /home/peascl/scl_de/clout_gene/cuffquant_out/Lifter-Scl-48-1/abundances.cxb,/home/peascl/scl_de/clout_gene/cuffquant_out/Lifter-Scl-48-BR/abundances.cxb /home/peascl/scl_de/clout_gene/cuffquant_out/PI240515-Scl-12-1/abundances.cxb,/home/peascl/scl_de/clout_gene/cuffquant_out/PI240515-Scl-12-BR/abundances.cxb /home/peascl/scl_de/clout_gene/cuffquant_out/PI240515-Scl-24-1/abundances.cxb,/home/peascl/scl_de/clout_gene/cuffquant_out/PI240515-Scl-24-BR/abundances.cxb /home/peascl/scl_de/clout_gene/cuffquant_out/PI240515-Scl-48-1/abundances.cxb,/home/peascl/scl_de/clout_gene/cuffquant_out/PI240515-Scl-48-BR/abundances.cxb /home/peascl/scl_de/clout_gene/cuffquant_out/Scl0205-12hr-1/abundances.cxb,/home/peascl/scl_de/clout_gene/cuffquant_out/Scl0205-12hr-BR/abundances.cxb /home/peascl/scl_de/clout_gene/cuffquant_out/Scl0205-24-1/abundances.cxb,/home/peascl/scl_de/clout_gene/cuffquant_out/Scl0205-24-BR/abundances.cxb /home/peascl/scl_de/clout_gene/cuffquant_out/Scl0205-48-1/abundances.cxb,/home/peascl/scl_de/clout_gene/cuffquant_out/Scl0205-48-BR/abundances.cxb