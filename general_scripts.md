#General scripts to manipulate files during analyses

Recover number of mapped reads from the Log file in the STAR output
```
for i in ./*/*.Log.final.out; 
	do sed -rn 's/ {1,}Uniquely mapped reads number \|\t([0-9]{2,})/\1/p' $i; print $i;
 done
```

Create a file with names for bam files to use with htseq-count script using parallel
```
ls *.bam | sed -rn 's/(.*)\.Aligned\.sortedByCoord\.out\.bam/\1/p' > BAM_files.txt
```
