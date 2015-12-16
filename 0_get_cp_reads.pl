#!/usr/bin/perl 

# Script that fishes chloroplast reads from whole-genome read sets

#Carlos P Cantalapiedra (1), Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)
#1) Estacion Experimental de Aula Dei-CSIC, Zaragoza, Spain
#2) Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
#3) Fundacion ARAID, Zaragoza, Spain

use strict;
use FindBin '$Bin';

my $DUKEXE = $Bin.'/bin/duk/duk';

my $refcpFASTA = $Bin.'/poaceae.fna';
# contains several cps from Poaceae species:
#>NC_015820 Acidosasa purpurea
#>NC_008591 Agrostis stolonifera
#>NC_014062 Anomochloa marantoidea
#>NC_015830 Bambusa emeiensis
#>NC_012927 Bambusa oldhamii
#>NC_011032 Brachypodium distachyon
#>NC_013273 Coix lacryma-jobi
#>NC_013088 Dendrocalamus latiflorus
#>NC_015831 Ferrocalamus rimosivaginus
#>NC_011713 Festuca arundinacea
#>EF115541 Hordeum vulgare subsp. vulgare cultivar Morex
#>NC_015803 Indocalamus longiauritus
#>NC_009950 Lolium perenne
#>AY522329 Oryza sativa (indica cultivar-group) isolate 93-11
#>AY522331 Oryza sativa (japonica cultivar-group) isolate PA64S
#>NC_015990 Panicum virgatum
#>NC_015817 Phyllostachys edulis
#>NC_015826 Phyllostachys nigra var. henonis 
#>NC_006084 Saccharum hybrid cultivar NCo 310 
#>NC_008602 Sorghum bicolor 
#>NC_002762 Triticum aestivum 
#>X86563 Zea mays complete 


if(!$ARGV[1] || !-d $ARGV[0]){ die "# usage: $_ <folder with all-read files> <output folder with cp-read files>\n"; }

my ($inpDIR,$outDIR) = (@ARGV);

mkdir($outDIR) if(!-d $outDIR);


opendir(READS,$inpDIR);
my @readfiles = grep {!/^\./} grep {/.fastq/ || /.fq/} readdir(READS);
closedir(READS); 

foreach my $readf (@readfiles)
{ 
  print "# fishing cp reads from $inpDIR/$readf ...\n";
  system("$DUKEXE -k 24 -c 2 -m $outDIR/cp-$readf -o $outDIR/cp-$readf.duk.log $refcpFASTA $inpDIR/$readf");
  system("gzip $outDIR/cp-$readf");
}
