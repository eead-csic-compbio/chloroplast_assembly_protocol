# chloroplast_assembly_protocol

A set of scripts for the assembly of chloroplast genomes out of whole-genome sequencing reads

**Authors**
Carlos P Cantalapiedra (1), Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)

1) Estación Experimental de Aula Dei-CSIC, Zaragoza, Spain
2) Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
3) Fundación ARAID, Zaragoza, Spain

## Software dependencies

The following software packages, pre-compiled with the scripts and located in bin/, are required for the protocol. They should work out-of-the-box in Linux x86_64 systems, but should be re-compiled otherwise. They haven't been tested under Windows:

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
| split_pairs | 0.5 | this site, uses seqtk code |

## Examples

* Fish chloroplast (cp) reads from whole genome library, using provided test reads:
```{shell}
./0_get_cp_reads.pl test/ test_cp/ 
```    

* Clean and trim reads to remove poor quality segments; output includes mean insert sizes and orientations:
```{shell}  
./1_cleanreads.pl test_cp reference.fna 
```

* Assemble cp genome from a single PE library:
```{shell}
cd test_cp
../2_assemble_reads.pl --PEfile cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz \
  --PEinsert 221 --ref ../reference.fna
```

* Assemble cp genome combining PE + MP libraries:
```{shell}
../2_assemble_reads.pl --PEfile cp-testPE.wind15_28.3crop70.mlen60.corr.12.fq.gz \
  --PEinsert 221 --ref ../reference.fna \
  --MPFile cp-testMP.wind15_28.3crop70.mlen60.corr.12.fq.gz --MPinsert 4295 
```

## Post-assembly inspection

The tutorial [HOWTOcheck_assembly.txt](HOWTOcheck_assembly.txt) contains recipes to analyze and validate your assemblies.
