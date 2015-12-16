# Carlos P Cantalapiedra (1), Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)
# 1) Estacion Experimental de Aula Dei-CSIC, Zaragoza, Spain
# 2) Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
# 3) Fundacion ARAID, Zaragoza, Spain

# This file explains to use the attached scripts for assemblying chloroplast genomes out of whole-genome reads.

0. Software dependencies
========================

Name  Shipped version Source
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

The attached tutorial 'HOWTOcheck_assembly.txt' contains recipes to analyze and validate your assemblies.


