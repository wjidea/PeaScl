# Trimmomatic.jar
#java -jar <path to trimmomatic.jar> PE [-threads <threads] [-phred33 | -phred64] [-trimlog 
#<logFile>] >] [-basein <inputBase> | <input 1> <input 2>] [-baseout <outputBase> |
#<unpaired output 1> <paired output 2> <unpaired output 2> <step 1>

TRIM_PATH='/home/wjidea/apps/Trimmomatic-0.32/trimmomatic-0.32.jar'
RAW_SEQ='/home/wjidea/Files/Turf_bac_rawseq/20130608_DNASeq_PE/raw_seqs'
ADAPTER='/home/wjidea/apps/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa'
FWPE=${1%%.fastq.gz}.PE.fq.gz
FWSE=${1%%.fastq.gz}.SE.fq.gz
RVPE=${2%%.fastq.gz}.PE.fq.gz
RVSE=${2%%.fastq.gz}.SE.fq.gz

java -jar $TRIM_PATH PE -trimlog trimlog.txt $1 $2 $FWPE $FWSE $RVPE $RVSE ILLUMINACLIP:$ADAPTER:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#echo $FWPE
