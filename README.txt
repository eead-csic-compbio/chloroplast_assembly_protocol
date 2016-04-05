# Carlos P Cantalapiedra (1), Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)
# 1) Estacion Experimental de Aula Dei-CSIC, Zaragoza, Spain
# 2) Escuela Politecnica Superior de Huesca, U.Zaragoza, Spain
# 3) Fundacion ARAID, Zaragoza, Spain

# This file explains to use the attached scripts for assemblying chloroplast genomes 
out of whole-genome reads.

0. Software dependencies
========================

This protocol has been tested on Linux x86_64 systems, although it should also work on Mac-OSX 
settings. It requires perl, which should be installed on all Linux environments, plus 
some third-party programs listed below, which are provided pre-compiled. In order to check 
whether they work on your machine, or to re-compile them otherwise, for instance on 
Mac-OSX, please run on the terminal:

perl ./install.pl

The script will tell which programs are set up (OK) and which require installing a compiler 
(gcc or g++) or java in order to run.

This is the full list of software, located in bin/, required for the protocol:
 
DUK March032011 http://duk.sourceforge.net
Trimmomatic 0.32  http://www.usadellab.org/cms/?page=trimmomatic
FastQC  V0.10.1 http://www.bioinformatics.babraham.ac.uk/projects/fastqc
musket  1.0.6 http://musket.sourceforge.net
BWA 0.7.6a  http://bio-bwa.sourceforge.net/
velvet  1.2.08  https://github.com/dzerbino/velvet/tree/master
SSPACE  2.0 http://www.baseclear.com/bioinformatics-tools
GapFiller v1-11 http://www.baseclear.com/bioinformatics-tools
seqtk NA  https://github.com/lh3/seqtk
split_pairs v0.5  [this site, uses seqtk code]

Optionally samtools is used in 'HOWTOcheck_assembly.txt':

samtools 0.1.19 [http://samtools.sourceforge.net]

1. Protocol at a glance
=======================
# Create directory with your samples ("SAMPLESDIR")
./0_get_cp_reads.pl SAMPLESDIR ASSEMBLYDIR
./1_cleanreads.pl ASSEMBLYDIR FASTA_REF_GENOME
# Create config file ("ASSEMBLYCONF")
./2_assemble_reads.pl ASSEMBLYDIR ASSEMBLYCONF --ref FASTA_REF_GENOME

2. Test run
===========
2.1. Using a reference genome
=============================

# The input files for this test are within 
# test/ directory

# Fish cp reads from WGS library
./0_get_cp_reads.pl test test_cp poaceae.fna

# Clean and trim reads to remove poor quality segments
./1_cleanreads.pl test_cp reference.fna

# Create config file test_cp/assembly_pe
cp test_cp/cleanreads.txt test_cp/assembly_pe

# and edit test_cp/assembly_pe file leaving only one row:

#1 testPE cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz FR 221 Sanger

# Finally, assemble cp genome
./2_assemble_reads.pl test_cp assembly_pe --ref ./reference.fna

2.2. Combining PE + MP libraries
================================

# Perform the steps from the previous section (2.1) up to ./1_cleanreads command (included).

# Then, create a different config file to use both read libraries.

cp test_cp/cleanreads.txt test_cp/assembly_mp

# and edit the file reordering rows so that testPE is number #1:

#1 testPE cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz FR 221 Sanger
#2 testMP cp-testMP.wind15_28.3crop70.mlen60.corr.12.fq.gz RF 4295 Sanger

# Finally, assemble combining PE + MP libraries
./2_assemble_reads.pl test_cp assembly_mp --ref reference.fna

2.3. de-novo, without reference genome
======================================

# The input files for this test are within 
# test/ directory

## Fish cp reads from WGS library
./0_get_cp_reads.pl test test_cp_noref poaceae.fna

# clean and trim reads to remove poor quality segments
./1_cleanreads.pl test_cp_noref noref

# Create config file test_cp_noref/assembly_pe

cp test_cp_noref/cleanreads.txt test_cp_noref/assembly_pe

# and edit test_cp/assembly_pe file leaving only one or two rows (as in sections 2.1 and 2.2).
# In this case we are leaving a single PE read library, and note that we have to provide orientation
# and insert size (FR and 221 in this example; see next section):

#1 testPE cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz FR 221 Sanger

# Finally, assemble cp genome
./2_assemble_reads.pl test_cp_noref assembly_pe --ref noref

3. Reference
============
3.1. 0_get_cp_reads.pl
======================

This script can be used to extract chloroplast reads from different WGS libraries.

./0_get_cp_reads.pl INPUT_DIR OUTPUT_DIR FASTA_CP_GENOMES

- INPUT_DIR: a directory containing the FASTQ files with reads to be used as input.
- OUTPUT_DIR: a path to a directory (which will be created if does not exist) to store output files.
- FASTA_CP_GENOMES: fasta file with sequences of chloroplast genomes from related species.

"poaceae.fna" file is included as an example of a file which could be used
when working with Poaceae chloroplast genomes.

3.2. 1_cleanreads.pl
====================

A script to improve quality of reads and estimate insert size and reads orientation.

./1_cleanreads.pl WORKING_DIR FASTA_REF_GENOME [-s]

- WORKING_DIR: a path to a directory with input files (the one created in the previous script).
- FASTA_REF_GENOME: to perform de-novo assembly this parameter should be "noref".
Otherwise, the name of the FASTA file to be used as reference.
- -s: skip Musket error correction. default: perform Musket error correction.

Other parameters could be changed by editing the script 1_cleanreads.pl.
- Parameters related to running time:

my $CPUTHREADS = 4; # number of CPU threads to be used for parallele jobs
my $READSAMPLESIZE = 10_000; # reads mapped to reference to estimate insert size

- These correspond to Trimmomatic trimming and filtering procedures:

my $TRIM5 = 3; # 5' bases to remove, adjust as needed after inspection of FastQC reports
my $TRIM3 = 3; # 3' bases to remove
my $MINREADLENGTH = 60; # this will allow assemblying with kmers up to this size, decrease if required
my $MINSURVIVALRATE = 50; # a warning will be used if less than these %reads survive trimmomatic

3.3. 2_assemble_reads.pl
====================

A script to assemble the reads into new chloroplast contig or contigs.

./2_assemble_reads.pl WORKING_DIR ASSEMBLY_NAME --ref FASTA_REF_GENOME
[--threads FLOAT] [--sample INTEGER] [--kmer INTEGER] [--outdir OUTPUT_DIR]

- WORKING_DIR: a path to a directory with input files (the one created in the previous script).
It will be used as output directory also, if --outdir option is not used (see below).
- --ref FASTA_REF_GENOME: to perform de-novo assembly this parameter should be "noref".
Otherwise, the name of the FASTA file to be used as reference.
- --threads FLOAT: number of CPU threads to use.
- --sample INTEGER: number of reads to be used for the assembly.
- --kmer INTEGER: k-mer size for Velvet assembler.
- --outdir WORKING_DIR: a path to a directory to store output files.


Other parameters could be changed by editing the script 1_cleanreads.pl.

- Parameters to be used by the assembler

my $VELVETGPARAMS = '-cov_cutoff 50 -min_contig_lgth 150 -exp_cov auto -unused_reads yes -scaffolding no';

4. Notes about reference genome, reads orientation and insert size
==================================================================

Assembling using a reference (reference.fna in examples in sections 2.1 and 2.2)
allows the scripts to estimate (thru BWA mappings) mean insert sizes and orientations
for each input sample.

Such information is included by 1_cleanreads.pl script in a config file (see section A.1)
and it is required to perform the assembly step (2_assemble_reads.pl).

However, when assembling de-novo (as in the example in section 2.3), insert size and 
orientation of reads can not be estimated by the scripts. The config file generated
has "nd" values in the orientation and insert size fields:

#1 testMP cp-testMP.wind15_28.3crop70.mlen60.corr.12.fq.gz nd nd Sanger
#2 testPE cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz nd nd Sanger

Those fields must be edited previous to the assembly step.
Orientation is FR for most standard PE runs, whereas MP runs are often RF.
Insert size should be estimated somehow, likely from DNA size selection
in wet-lab previous to sequencing. A valid config file would look like:

#1 testPE cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz FR 221 Sanger
#2 testMP cp-testMP.wind15_28.3crop70.mlen60.corr.12.fq.gz RF 4295 Sanger

Note also that the last field corresponds to base quality encoding
(see File Formats section below), which is set to "Sanger" by default but
should be changed according to the FASTQ encoding of the files being used as input.

5. Post-assembly inspection
===========================

The attached tutorial 'HOWTOcheck_assembly.txt' contains recipes to analyze and validate your 
assemblies.


6. Propagation of annotated features to assemblies
==================================================

If a good quality master annotation in GenBank format is at hand, it can be used to propagate 
its features to several assembled chloroplast genomes with accompanying script _annot_fasta_from_gbk.pl .
The way to invoke is as follows, assuming strain XYZ was assembled:

perl _annot_fasta_from_gbk.pl reference.gbk assembly.fa XYZ_chloroplast.gbk XYZ

NOTE: This script requires setting a valid path to BLAST+ binaries installed on your system.


7. Generating tracks from aligned assemblies
============================================

If several chloroplast assemblies are aligned they can then be further used to produce tracks for software CIRCOS <http://circos.ca>.
Script _check_matrix.pl does just that; it requires setting a reference strain 
and a window size, please check the code. An example call would be:

perl _check_matrix.pl chloroplast_alignment.fna

A. File formats
===============
A.1 Format of assembly config files (cleanreads.txt)
====================================================

Fields (columns) are separated by a blank space:
1. File number: the file "#1" will be used as main PE file.
	A second file "#2" will be used as auxiliary MP file.
2. Name: a short name identifying each sample.
3. Input filename: file with clean reads to be used as input 
	for the assembly step.
4. Orientation: relative orientation of the paired and mate-pair reads.
	"nd" (non-defined): it should be changed to another value, before assembly step.
	"FR" (forward-reverse): typical of Illumina PE reads.
	"RF" (reverse-forward): found in some MP libraries, for example.
5. Insert size: mean fragment size in the library.
	"nd" (non-defined): it should be changed to a valid value before assembly step.
6. Base quality encoding: quality codes change from a sequencing platform to another.
	"Sanger": valid for Sanger, Illumina 1.9 and Phred-33, for example.
	"1.5": valid for Illumina 1.5 and Phred+64, for example.

##
