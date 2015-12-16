#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
# Parameters:
# $1: basename of fastq file, with directories
# $2: leading bases
# $3: tail bases
# $4: if "yes", fastqc_compact report is generated.
BINPATH="`dirname \"$0\"`"

if [ -f $1.gz ]
then
    inf=$1.gz
else
    inf=$1.fastq
fi

java -jar ${BINPATH}/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 4 -phred33 $inf $1.$2crop$3.fastq CROP:$3 HEADCROP:$2

if [ "$4" == "yes" ];
then `./fastqc_compact.sh $1.$2crop$3`;
fi;
