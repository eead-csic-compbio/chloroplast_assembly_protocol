#!/usr/bin/perl

# Script that assembles a chloroplast circular genomes by calling velvet, SSPACE & GapFiller
# It expects a pair-end (PE) read library and optionally a mate-pair (MP) read library

#Carlos P Cantalapiedra (1), Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)
#1) Estacion Experimental de Aula Dei-CSIC, Zaragoza, Spain
#2) Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
#3) Fundacion ARAID, Zaragoza, Spain

use strict;
use FindBin '$Bin';
use Getopt::Long;
use File::Basename;

my $BINPATH = $Bin.'/bin';

my $SPLITEXE   = $BINPATH.'/split_pairs_v0.5/split_pairs.pl';
my $BWAEXE     = $BINPATH.'/bwa-0.7.6a/bwa';
my $VELVETH    = $BINPATH.'/velvet_1.2.08/velveth'; # compiled with make 'LONGSEQUENCES=1'
my $VELVETG    = $BINPATH.'/velvet_1.2.08/velvetg';
my $SSPACEXE   = $BINPATH."/SSPACE-BASIC-2.0_linux-x86_64/SSPACE_Basic_v2.0.pl -n 32 ";
my $GAPFILLEXE = $BINPATH."/GapFiller_v1-11_linux-x86_64/GapFiller.pl -o 5 ";
my $SEQTKEXE   = $BINPATH.'/seqtk/seqtk seq';

# tested on our cp assemblies
my $VELVETGPARAMS = '-cov_cutoff 50 -min_contig_lgth 150 -exp_cov auto -unused_reads yes -scaffolding no';

# defaults
my $SAMPLESIZE = 500_000; # number of reads to be used; occasionally 1M yielded better assemblies
my $CPUTHREADS = 4;       # number of CPU threads to be used for parallele jobs
my $KMER       = 47;      # optimal for B.distachyon, 59, 73, 81 were also useful occasionally
my $PEfile     = '';
my $PEinsert   = 0;
my $PEencoding = 'Sanger';# encoding of quality values in FASTQ [Illumina 1.5 -> Phred+64, Sanger/Illumina 1.9  -> Phred+33] 
my $MPfile     = '';
my $MPinsert   = 0;
my $MPorient   = 'RF';    # default orientation of MP reads
my $MPencoding = 'Sanger';
my $outDIR     = ''; 
my $refFASTA   = '';
my $help;
  
usage() if(@ARGV < 1 || 
  !GetOptions(
    'help|h'     => \$help, 
    'ref=s'      => \$refFASTA, 
    'PEfile=s'   => \$PEfile, 
    'PEinsert=i' => \$PEinsert,
    'PEenc=s'    => \$PEencoding,
    'MPfile=s'   => \$MPfile, 
    'MPinsert=i' => \$MPinsert,
    'MPenc=s'    => \$MPencoding,
    'MPorient=s' => \$MPorient,
    'threads=i'  => \$CPUTHREADS,
    'sample=i'   => \$SAMPLESIZE,
    'kmer=i'     => \$KMER,
    'outdir=s'   => \$outDIR) || defined($help) );
  
sub usage
{
  print   "\n[options]: \n";
  print   "-h this message\n";
  print   "--ref      reference genome in FASTA format             (required)\n"; 
  print   "--PEfile   input pair-end file in FASTQ format          (required)\n"; 
  print   "--PEinsert mean insert size of PE read pairs [integer]  (required)\n";
  print   "--PEenc    quality encoding of PE reads [Sanger|1.5]    (optional, default=$PEencoding)\n";
  print   "--MPfile   input mate-pair file in FASTQ format         (optional)\n"; 
  print   "--MPinsert mean insert size of MP read pairs [integer]  (optional, required with --MPfile)\n";
  print   "--MPenc    quality encoding of MP reads [Sanger|1.5]    (optional, default=$MPencoding)\n";
  print   "--MPorient orientation of MP read pairs [RF|FR]         (optional, default=$MPorient)\n";
  print   "--threads  number of CPU threads to use                 (optional, default=$CPUTHREADS)\n"; 
  print   "--sample   number of reads to be assembled              (optional, default=$SAMPLESIZE)\n";
  print   "--kmer     k-mer size for Velvet assembler              (optional, default=$KMER)\n"; 
  print   "--outdir   folder to store results                      (optional, default=assembly_kmer$KMER\_sample$SAMPLESIZE)\n\n";
  exit(-1);
}

if(!$refFASTA || !-s $refFASTA)
{ die "\n# $0 : need a valid --ref FASTA file, exit\n"; }

if(!$PEfile || !-s $PEfile || !$PEinsert || $PEinsert < 1)
{ die "\n# $0 : need --PEfile and --PEinsert, exit\n"; }
elsif($PEfile =~ /mlen(\d+)/)
{
  die "\n# $0 : $PEfile might contain reads < $1 b, please choose a smaller kmer\n" if($1 < $KMER);
}  

if($MPfile && ( !-s $MPfile || !$MPinsert || $PEinsert < 1) )
{ die "\n# $0 : --MPfile and/or --MPinsert are not valid, exit\n"; }

if(($PEencoding ne 'Sanger' && $PEencoding ne '1.5') || 
  ($MPencoding ne 'Sanger' && $MPencoding ne '1.5'))
{ die "\n# $0 : valid encodings are: Sanger|1.5, see [https://en.wikipedia.org/wiki/FASTQ_format]\n"; }

if($MPorient ne 'RF' && $MPorient ne 'FR')
{ die "\n# $0 : valid orientations are: FR|RF, exit\n"; }

if(!$outDIR){ $outDIR = "assembly_kmer$KMER\_sample$SAMPLESIZE" }
if(!-s $outDIR){ mkdir($outDIR) }
else{ print "# re-using existing output folder '$outDIR'\n\n"; }

printf("\n# %s --ref %s --PEfile %s --PEinsert %d --PEenc %s \\\n".
      "#  --MPfile %s --MPinsert %d --MPenc %s --MPorient %s \\\n".
      "#  --threads %d --sample %d --kmer %d --outdir %s\n\n",
	$0,$refFASTA,$PEfile,$PEinsert,$PEencoding,
  $MPfile,$MPinsert,$MPencoding,$MPorient,
  $CPUTHREADS,$SAMPLESIZE,$KMER,$outDIR); 

open(COMMAND,">$outDIR/command.txt");
printf COMMAND ("%s --ref %s --PEfile %s --PEinsert %d --PEenc %s \\\n".
      "  --MPfile %s --MPinsert %d --MPenc %s --MPorient %s \\\n".
      "  --threads %d --sample %d --kmer %d --outdir %s\n\n",
	$0,$refFASTA,$PEfile,$PEinsert,$PEencoding,
  $MPfile,$MPinsert,$MPencoding,$MPorient,
  $CPUTHREADS,$SAMPLESIZE,$KMER,$outDIR);
close(COMMAND);
      
#############################################################  
  
my ($nlines,$enc,$encMP,$root,$rootMP,$gapsOK,$velvet_cmd,@trash);
my ($intlfileMP,$pair1fileMP,$pair2fileMP,$samplefileMP,$samfileMP,$logfileMP);
my ($pair1file,$pair2file,$spaceparamfile,$gapparamfile,$finalfile);
my ($infile,$samplefile,$samfile,$logfile,$tmpfile,$intlfile);

## 0) check main input file and params
$infile = $PEfile;
$root = basename($infile);
$root =~ s/\.fq\.gz//;
  
if($PEencoding eq '1.5'){ $enc = 'Phred+64' }
else{ $enc = 'Phred+33' } 

if($MPencoding eq '1.5'){ $encMP = 'Phred+64' }
else{ $encMP = 'Phred+33' } 

## 1) extract sample read pairs and correct encoding and orientation if required

# 1.1) compulsory PE reads
$nlines = $SAMPLESIZE * 8; # 4 for F read and 4 for R

$tmpfile = $root.'.tmp';
$intlfile = $root .'.pairs.fq';
$pair1file = $root .'.1.fq';
$pair2file = $root .'.2.fq';
$samplefile = $root .'.sample.fq';

system("zcat $infile | head -$nlines > $tmpfile");
system("$SPLITEXE -n -i $tmpfile -l $intlfile");
unlink($tmpfile);
system("$SPLITEXE -n -i $intlfile -1 $pair1file -2 $pair2file");

# correct encoding to Phred+33 (Sanger) for BWA mem
$nlines = $SAMPLESIZE * 4;
if($enc eq 'Phred+64')
{
	print "# converting qualities to Sanger Phred+33 encoding...\n";
	system("head -$nlines $intlfile > $tmpfile");

	open(SAMPLE,">$samplefile") || die "# cannot create $samplefile\n";
	open(TMP,$tmpfile) || die "# cannot read $tmpfile\n";
	while(<TMP>)
	{
		chomp;
		#http://seqanswers.com/forums/showthread.php?t=5210
		if($.%4==0){ tr/\x40-\xff\x00-\x3f/\x21-\xe0\x21/ }
		print SAMPLE "$_\n";
	}	
	close(TMP);
	close(SAMPLE);
}
else{ unlink($samplefile); system("ln -s $intlfile $samplefile"); }

push(@trash,$intlfile,$pair1file,$pair2file,$samplefile);

# 1.2) optional MP reads
if($MPfile)
{
  $rootMP = basename($MPfile);
  $rootMP =~ s/\.fq\.gz//;

  $intlfileMP = $rootMP .'.pairs.fq';
  $pair1fileMP = $rootMP .'.1.fq';
  $pair2fileMP = $rootMP .'.2.fq';
  $samplefileMP = $rootMP .'.sample.fq';

  system("zcat $MPfile | head -$nlines > $tmpfile");
  system("$SPLITEXE -n -i $tmpfile -l $intlfileMP");
  unlink($tmpfile);
  system("$SPLITEXE -n -i $intlfileMP -1 $pair1fileMP -2 $pair2fileMP");

  # correct encoding to Phred+33 (Sanger) for BWA mem
  $nlines = $SAMPLESIZE * 4;
  if($encMP eq 'Phred+64')
  {
	  print "# converting qualities to Sanger Phred+33 encoding (MP)...\n";
	  system("head -$nlines $intlfileMP > $tmpfile");

	  open(SAMPLE,">$samplefileMP") || die "# cannot create $samplefileMP\n";
	  open(TMP,$tmpfile) || die "# cannot read $tmpfile\n";
	  while(<TMP>)
	  {
		  chomp;
		  if($.%4==0){ tr/\x40-\xff\x00-\x3f/\x21-\xe0\x21/ }
		  print SAMPLE "$_\n";
	  }	
	  close(TMP);
	  close(SAMPLE);
  }
  else{ unlink($samplefileMP); system("ln -s $intlfileMP $samplefileMP"); }

  push(@trash,$intlfileMP,$pair1fileMP,$pair2fileMP,$samplefileMP);

  if($MPorient eq 'RF') # get rc of mate pairs
  {
    print "# getting rc of MP reads ...\n";
	  system("$SEQTKEXE -r $pair1fileMP > $rootMP.rc.1.fq");
	  system("$SEQTKEXE -r $pair2fileMP > $rootMP.rc.2.fq");
	  unlink($pair1fileMP,$pair2fileMP);
    $pair1fileMP = $rootMP .'.rc.1.fq';
    $pair2fileMP = $rootMP .'.rc.2.fq';
    push(@trash,$pair1fileMP,$pair2fileMP);
  }
}


### 2) map reads to reference prior to assembly
$samfile = $root . '.sam';
$logfile = $root . '.bwa.log';

system("$BWAEXE index $refFASTA"); 
  
print "# BWA command:\n$BWAEXE mem -t $CPUTHREADS -p $refFASTA $samplefile > $samfile 2> $logfile\n";
open(BWA,"$BWAEXE mem -t $CPUTHREADS -p $refFASTA $samplefile > $samfile 2> $logfile |")
	|| die "# cannot run $BWAEXE mem -t $CPUTHREADS -p $refFASTA $samplefile > $samfile 2> $logfile\n";
while(<BWA>){ }
close(BWA); 

$ENV{'LC_ALL'} = 'POSIX';
system("sort $samfile > $tmpfile");
rename($tmpfile,$samfile);
push(@trash,$samfile,$logfile);

if($MPfile)
{
  $samfileMP = $rootMP . '.sam';
  $logfileMP = $rootMP . '.bwa.log';

  print "# BWA command:\n$BWAEXE mem -t $CPUTHREADS $refFASTA $pair1fileMP $pair2fileMP > $samfileMP 2> $logfileMP\n";
  open(BWA,"$BWAEXE mem -t $CPUTHREADS $refFASTA $pair1fileMP $pair2fileMP > $samfileMP 2> $logfileMP |")
	|| die "# cannot run $BWAEXE mem -t $CPUTHREADS $refFASTA $pair1fileMP $pair2fileMP> $samfileMP 2> $logfileMP\n";
  while(<BWA>){ }
  close(BWA); 

  system("sort $samfileMP > $tmpfile");
  rename($tmpfile,$samfileMP);
  push(@trash,$samfileMP,$logfileMP);
}



## 3) do assemble with velvet in two steps

# 3.1) make kmer hash table with reads mapped to split chloroplast reference
$velvet_cmd = "$VELVETH $outDIR $KMER -reference $refFASTA -shortPaired -sam $samfile";
if($MPfile)
{
  $velvet_cmd .= " -shortPaired2 -sam $samfileMP";
}

print "# velvetg command:\n$velvet_cmd\n";
open(VELVETH,"$velvet_cmd |") || die "# cannot run $velvet_cmd\n";
while(<VELVETH>){ }
close(VELVETH);

# 3.2) assemble hashed kmers with proper insert size
$velvet_cmd = "$VELVETG $outDIR -ins_length $PEinsert $VELVETGPARAMS";
#Final graph has 19 nodes and n50 of 22447, max 29336, total 131522, using 476527/483841 reads

if($MPfile)
{
  $velvet_cmd .= " -ins_length2 $MPinsert -shortMatePaired2 yes";
  #Final graph has 14 nodes and n50 of 75334, max 75334, total 134474, using 925984/984681 reads
}

print "# velveth command:\n$velvet_cmd\n";
open(VELVETG,"$velvet_cmd |") || die "# cannot run $velvet_cmd\n";
while(<VELVETG>)
{ 
	print if(/^Final graph/);
}
close(VELVETG);#system("head $outDIR/stats.txt"); 

$finalfile = "$outDIR/contigs.fa";
if(!-s $finalfile)
{
  unlink(@trash);
  die "\n# Velvet could not assembly any contig ($finalfile), exit\n\n";
}  

## 4) make scaffolds and fill gaps
$spaceparamfile = $root .'.sspace.params';
$gapparamfile = $root .'.gapfiller.params';

# 4.1) create scaffolding parameter file
# http://www.vcru.wisc.edu/simonlab/bioinformatics/programs/sspace/F132-01%20SSPACE_Basic_User_Manual_v2.0.pdf
# lib1 libscaf1.fq libscaf2.fq 4300 0.11 RF
open(PARAMS,">$spaceparamfile") || die "# cannot create $spaceparamfile\n";
print PARAMS "lib1 $pair1file $pair2file $PEinsert 0.2 FR\n";
if($MPfile)
{
  print PARAMS "lib2 $pair1fileMP $pair2fileMP $MPinsert 0.2 FR\n";
}
close(PARAMS);

# 4.2) merge contigs into scaffolds
print "# SSPACE command:\n$SSPACEXE -T $CPUTHREADS -l $spaceparamfile -s $outDIR/contigs.fa -b $root.sspace\n";
open(SSPACE,"$SSPACEXE -T $CPUTHREADS -l $spaceparamfile -s $outDIR/contigs.fa -b $root.sspace |\n") 
|| die "# cannot run $SSPACEXE -T $CPUTHREADS -l $spaceparamfile -s $outDIR/contigs.fa -b $root.sspace\n";
while(<SSPACE>){ }
close(SSPACE);

$gapsOK = 0;
open(EVIDENCE,"$root.sspace.final.evidence") if(-s "$root.sspace.final.evidence");
while(<EVIDENCE>)
{
	if(/gaps\d+/){ $gapsOK++ }
}
close(EVIDENCE);
system("rm -rf bowtieoutput pairinfo reads intermediate_results"); 


# 4.3) fill gaps if required
if($gapsOK)
{
	open(PARAMS,">$gapparamfile") || die "# cannot create $gapparamfile\n";
	print PARAMS "lib1 bowtie $pair1file $pair2file $PEinsert 0.2 FR\n";
  if($MPfile)
  {
    print PARAMS "lib2 bowtie $pair1fileMP $pair2fileMP $MPinsert 0.2 FR\n"; 
  }
	close(PARAMS);

	print "# GAPFILLER command:\n$GAPFILLEXE -T $CPUTHREADS -l $gapparamfile -s $root.sspace.final.scaffolds.fasta -b $root.gapfiller\n";
	open(GAPFILL,"$GAPFILLEXE -T $CPUTHREADS -l $gapparamfile -s $root.sspace.final.scaffolds.fasta -b $root.gapfiller |\n")
  	|| die "# cannot run $GAPFILLEXE -T $CPUTHREADS -l $gapparamfile -s $root.sspace.final.scaffolds.fasta -b $root.gapfiller\n";
	while(<GAPFILL>){ }
	close(GAPFILL);

  system("rm -rf $root.gapfiller/alignoutput $root.gapfiller/reads $root.gapfiller/intermediate_results");

  $finalfile = "$root.gapfiller/$root.gapfiller.gapfilled.final.fa";
}
else{ $finalfile = "$root.sspace.final.scaffolds.fasta" }

print "\n# final assembly: $finalfile\n\n";

unlink(@trash);


