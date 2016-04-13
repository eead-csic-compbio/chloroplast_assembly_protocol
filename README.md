# Chloroplast assembly protocol

A set of scripts for the assembly of chloroplast genomes out of whole-genome sequencing reads

**Authors**
Carlos P Cantalapiedra (1), Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)

1. Estación Experimental de Aula Dei-CSIC, Zaragoza, Spain
2. Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
3. Fundación ARAID, Zaragoza, Spain

## Software dependencies

This protocol has been tested on Linux x86_64 systems, although it should also work on Mac-OSX settings.
It requires perl, which should be installed on all Linux environments, plus some third-party programs listed below,
which are provided pre-compiled. In order to check whether they work on your machine, or to re-compile them otherwise,
for instance on Mac-OSX, please run on the terminal:
```{shell}  
perl ./install.pl
```
The script will tell which programs are set up (OK) and which require installing a compiler
(gcc or g++) or java in order to run.

This is the full list of software, located in bin/, required for the protocol:

| Name | shipped version | Source |
|:-----|:---------------:|:-------|
| DUK | March 03, 2011 | <http://duk.sourceforge.net> |
| Trimmomatic | 0.32 | <http://www.usadellab.org/cms/?page=trimmomatic> |
| FastQC | 0.10.1 | <http://www.bioinformatics.babraham.ac.uk/projects/fastqc> |
| musket | 1.0.6 | <http://musket.sourceforge.net> |
| BWA | 0.7.6a | <http://bio-bwa.sourceforge.net> |
| velvet | 1.2.08 | <https://github.com/dzerbino/velvet/tree/master> |
| SSPACE | 2.0 | <http://www.baseclear.com/bioinformatics-tools> |
| GapFiller | v1-11 | <http://www.baseclear.com/bioinformatics-tools> |
| seqtk | | <https://github.com/lh3/seqtk> |
| split_pairs | 0.5 | <https://github.com/eead-csic-compbio/split_pairs> |

Optionally [samtools](http://samtools.sourceforge.net) is used in [HOWTOcheck_assembly.txt](HOWTOcheck_assembly.txt).

## Input reads

This pipeline has been tested with paired-end (PE) and mate-pairs (MP) read libraries,
with the following FASTQ quality encodings: Illumina 1.5 (Phred+64) and Sanger/Illumina (Phred+33). 
The code assumes that PE libraries have FR orientation, meaning that read pairs are "--><--". 
MP libraries are expected to be RF but can be set to FR as well.

## Reference sequences

* [poaceae.fna](poaceae.fna) is a FASTA file containing chloroplast (cp) genomes of several representative Poaceae species which provide DNA k-mers that help fishing out cp reads out of whole-genome read sets. It is used by [0_get_cp_reads.pl](0_get_cp_reads.pl). This file should be adapted to the particular phylogenetic group under study.

* [reference.fna](reference.fna) is required by this protocol, and it should contain the cp genome of a related species.

## Examples

Input reads are within "test" directory.

##### Reference-guided assembly

* Fish cp reads from whole genome library, using provided test reads (see [flowchart](./pics/0_get_cp_reads_1_cleanreads.png)):
```{shell}
./0_get_cp_reads.pl test test_cp poaceae.fna
```    

* Clean and trim reads to remove poor quality segments; output includes mean insert sizes and orientations (see [flowchart](./pics/0_get_cp_reads_1_cleanreads.png)):
```{shell}  
./1_cleanreads.pl -folder test_cp -ref reference.fna 
```

######- using a single PE library

* Create a config file "test_cp/assembly_pe" `cp test_cp/cleanreads.txt test_cp/assembly_pe`
and edit it leaving only one row:

    > 1 testPE cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz FR 221 1.5

* Assemble cp genome from a single PE library (see [flowchart-1](./pics/2_assemble_reads-1.png) and [flowchart-2](./pics/2_assemble_reads-2.png)):
```{shell}
./2_assemble_reads.pl test_cp assembly_pe -ref ./reference.fna
```
######- using both PE and MP libraries

* Then, create a different config file to use both read libraries `cp test_cp/cleanreads.txt test_cp/assembly_mp` 
and edit it reordering rows so that testPE is number 1:

    > 1 testPE cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz FR 221 Sanger
    
    > 2 testMP cp-testMP.wind15_28.3crop70.mlen60.corr.12.fq.gz RF 4295 Sanger

* Assemble cp genome combining PE + MP libraries (see [flowchart-1](./pics/2_assemble_reads-1.png) and [flowchart-2](./pics/2_assemble_reads-2.png)):
```{shell}
./2_assemble_reads.pl test_cp assembly_mp -ref reference.fna
```

##### Example of de-novo assembly

* Fish cp reads from WGS library

`./0_get_cp_reads.pl test test_cp_noref poaceae.fna`

* clean and trim reads to remove poor quality segments

`./1_cleanreads.pl -folder test_cp_noref `

* Create config file test_cp_noref/assembly_pe `cp test_cp_noref/cleanreads.txt test_cp_noref/assembly_pe`
and edit test_cp/assembly_pe file leaving one (PE reads) or two rows (PE + MP reads; see previous examples).

    > In this case we are leaving a single PE library, and note that we have to provide orientation
    > and insert size (FR and 221 in this example):
    
    > 1 testPE cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz FR 221 Sanger

* Finally, assemble cp genome

`./2_assemble_reads.pl test_cp_noref assembly_pe`

###### For more info and parameters of the scripts see [README.txt](README.txt)

In addition, examples yielding assembled data could be run following the steps at:

- [HOWTO_brachy.txt](HOWTO_brachy.txt): assembly with Brachypodium Bd21 data (PE).
- [HOWTO_barley.txt](HOWTO_barley.txt): assembly with barley cultivar Morex data (PE + MP).

## Post-assembly inspection

The tutorial [HOWTOcheck_assembly.txt](HOWTOcheck_assembly.txt) contains recipes to analyze and validate your assemblies.


## Propagation of annotated features to assemblies

If a good quality master annotation in GenBank format is at hand, it can be used to propagate 
its features to several assembled chloroplast genomes with accompanying script [_annot_fasta_from_gbk.pl](_annot_fasta_from_gbk.pl).
The way to invoke it is as follows, assuming strain XYZ was assembled:
```{shell}
perl _annot_fasta_from_gbk.pl reference.gbk assembly.fa XYZ_chloroplast.gbk XYZ
```

__NOTE__: This script requires setting a valid path to BLAST+ binaries installed on your system, plus [Bioperl](http://www.bioperl.org/wiki/Main_Page) installed on your system.


## Generating tracks from aligned assemblies

If several chloroplast assemblies are aligned they can then be further used to produce tracks for software
[CIRCOS](http://circos.ca). Script [_check_matrix.pl](_check_matrix.pl) does just that; it requires setting a reference strain 
and a window size, please check the code. An example call would be:

```{shell}
perl _check_matrix.pl chloroplast_alignment.fna
```
