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

1. Test run
===========

# fish cp reads from whole genome library
./0_get_cp_reads.pl test/ test_cp/ 

# clean and trim reads to remove poor quality segments
# output contains mean insert sizes and orientations
./1_cleanreads.pl test_cp reference.fna 

# assemble cp genome from a single PE library
cd test_cp
../2_assemble_reads.pl --PEfile cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz \
  --PEinsert 221 --ref ../reference.fna
  
# assemble cp genome combining PE + MP libraries
../2_assemble_reads.pl --PEfile cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz \
  --PEinsert 221 --ref ../reference.fna \
  --MPFile cp-testMP.wind15_28.3crop70.mlen60.corr.12.fq.gz --MPinsert 4295 
  
  
  
  
  
2. Post-assembly inspection
===========================

The attached tutorial 'HOWTOcheck_assembly.txt' contains recipes to analyze and validate your 
assemblies.



3. Propagation of annotated features to assemblies
==================================================

If a good quality master annotation in GenBank format is at hand, it can be used to propagate 
its features to several assembled chloroplast genomes with accompanying script _annot_fasta_from_gbk.pl .
The way to invoke is as follows, assuming strain XYZ was assembled:

perl _annot_fasta_from_gbk.pl reference.gbk assembly.fa XYZ_chloroplast.gbk XYZ

NOTE: This script requires setting a valid path to BLAST+ binaries installed on your system.



4. Generating tracks from aligned assemblies
============================================

If several chloroplast assemblies are aligned they can then be further used to produce tracks for software CIRCOS <http://circos.ca>. Script _check_matrix.pl does just that; it requires setting a reference strain 
and a window size, please check the code. An example call would be:

perl _check_matrix.pl chloroplast_alignment.fna
