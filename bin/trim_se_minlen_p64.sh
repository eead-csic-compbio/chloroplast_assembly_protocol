#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
# Parameters:
# $1: basename of fastq file, including directories (no extension)
# $2: minlen value
# $3: if "yes", generate fastqc_compact report
BINPATH="`dirname \"$0\"`"

if [ -f $1.gz ]
then
    inf=$1.gz
else
    inf=$1.fastq
fi

java -jar ${BINPATH}/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 4 -phred64 $inf $1.mlen$2.fastq MINLEN:$2

if [ "$3" == "yes" ];
then `./fastqc_compact.sh $1.mlen$2`;
fi;

