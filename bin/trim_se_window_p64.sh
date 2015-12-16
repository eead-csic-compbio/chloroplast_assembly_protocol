#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
# 
# Parameters:
# $1: basename of fastq file, including directories (no extension)
# $2: number of bases
# $3: quality threshold
# $4: if "yes", fastqc_compact report is generated
BINPATH="`dirname \"$0\"`"

if [ -f $1.gz ]
then
    inf=$1.gz
else
    inf=$1.fastq
fi

java -jar ${BINPATH}/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 4 -phred64 $inf $1.wind$2_$3.fastq SLIDINGWINDOW:$2:$3

if [ "$4" == "yes" ];
then `./fastqc.sh $1.wind$2_$3`;
fi;

