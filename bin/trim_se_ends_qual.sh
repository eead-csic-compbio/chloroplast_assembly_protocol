#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
# Parameters:
# $1: basename of file, including directories but no extension
# $2: leading quality
# $3: trailing quality
# $4: if "yes", fastqc_compact report is generated
BINPATH="`dirname \"$0\"`"

java -jar ${BINPATH}/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 4 -phred33 $1.fastq $1.$2qtrim$3.fastq LEADING:$2 TRAILING:$3 

if [ $4 == "yes" ];
then `./fastqc_compact.sh $1.$2qtrim$3`;
fi;
