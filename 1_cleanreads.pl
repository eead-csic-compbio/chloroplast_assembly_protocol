#!/usr/bin/perl 

# Script that cleans reads to be used for chloroplast genome assembly
# Uses a eeference genome; should be from a species as close as possible to the target species 

#Carlos P Cantalapiedra (1), Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)
#1) Estacion Experimental de Aula Dei-CSIC, Zaragoza, Spain
#2) Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
#3) Fundacion ARAID, Zaragoza, Spain

use strict;
use FindBin '$Bin';

my $CPUTHREADS = 4; # number of CPU threads to be used for parallele jobs

my $READSAMPLESIZE = 10_000; # reads mapped to reference to estimate insert size

my $TRIM5 = 3; # 5' bases to remove, adjust as needed after inspection of FastQC reports
my $TRIM3 = 3; # 3' bases to remove

my $MINREADLENGTH = 60; # this will allow assemblying with kmers up to this size, decrease if required

my $MINSURVIVALRATE = 50; # a warning will be used if less than these %reads survive trimmomatic

my $BINPATH = $Bin.'/bin';

my $FASTQCEXE     = $BINPATH.'/FastQC/fastqc';
my $SPLITPAIRSEXE = $BINPATH.'/split_pairs_v0.5/split_pairs.pl';
my $MUSKETEXE     = $BINPATH.'/musket-1.0.6/musket';
my $BWAEXE        = $BINPATH.'/bwa-0.7.6a/bwa';

my ($inpDIR,$refFASTA);

if(!$ARGV[1]){ die "# usage: $_ <folder with reads, results will be added there> <reference FASTA genome>\n"; }
else
{ 
  ($inpDIR,$refFASTA) = (@ARGV);
  print "# input_folder=$inpDIR reference=$refFASTA\n";
  print "# TRIM5=$TRIM5 TRIM3=$TRIM3 MINREADLENGTH=$MINREADLENGTH MINSURVIVALRATE=$MINSURVIVALRATE\n\n"; 
}

my ($inpDIR) = (@ARGV);
my ($encoding,$readlength,$N,$trimtype,$trimover,$minlength,@trash,$insert);
my ($readf,$gzfile,$outFQfile,$outFQsumfile,$outFQzipfile,$root,$orient);
my ($outFQfile2,$outFQsumfile2,$outFQzipfile2,$FQreportsDIR,$trimmedfile);

$FQreportsDIR = $inpDIR.'/reports/';
mkdir($FQreportsDIR) if(!-d $FQreportsDIR);

# process all files in input folder
opendir(READS,$inpDIR);
my @readfiles = sort grep { !-d "$inpDIR/$_" } grep {!/^\./} grep {!/.corr/} grep {/f[ast]*q.gz/} readdir(READS);
closedir(READS); 

foreach $readf (@readfiles)
{
	$gzfile = $readf; $gzfile =~ s/\.fastq//; 
	$root = $gzfile; $root =~ s/\.gz//;
	$outFQfile = $FQreportsDIR."/$gzfile"; $outFQfile =~ s/\.gz$/_fastqc\/fastqc_data.txt/;
	$outFQsumfile = $FQreportsDIR."/$gzfile"; $outFQsumfile =~ s/\.gz$/_fastqc\/summary.txt/;
	$outFQzipfile = $FQreportsDIR."/$gzfile"; $outFQzipfile =~ s/\.gz$/_fastqc.zip/;
	($encoding,$readlength,$N,$trimtype,$orient) = ('',0,0,'_p64','');
	@trash = ();

	print "> $root\n"; #print $inpDIR.'/'.$root.'*.corr.fq.gz'."\n";

  # make symb link
	system("ln -s $inpDIR/$readf $gzfile") if(!-e $gzfile);
	push(@trash,$gzfile);

	# run fastqc
	if(!-s $outFQfile)
	{
		open(FQ,"$FASTQCEXE $gzfile -o $FQreportsDIR -t $CPUTHREADS |") 
			|| die "# cannot run $FASTQCEXE $gzfile -o $FQreportsDIR -t $CPUTHREADS \n";
		while(<FQ>){ }
		close(FQ);
		
		if(-e $outFQzipfile)
		{
			system("cat $outFQsumfile");
			unlink($outFQzipfile);
		}
	}

	# find out encoding and length of these raw reads 
	#Encoding       Illumina 1.5    
	#Total Sequences        22604614        
	#Filtered Sequences     0       
	#Sequence length        95-100  
	open(FQDATA,$outFQfile) || die "# cannot read $outFQfile\n";
	while(<FQDATA>)
	{
		if(/^Encoding\t(.*)/){ $encoding = $1; } 
		elsif(/^Sequence length\t(\S+)/ || /^Sequence length\t\S+?-(\S+)/){ $readlength = $1; }
		elsif(/^Total Sequences\t(\S+)/){ $N = $1; }
	}	
	close(FQDATA); 

	if($encoding =~ "Illumina 1.8" || $encoding =~ "Sanger"){ $trimtype = '' }

	print "# Encoding=$encoding($trimtype) Total Sequences=$N\n"; 
     

	# clean with consecutive trimmomatic ops
	$trimover = $readlength-$TRIM5-$TRIM3;
	$minlength = $MINREADLENGTH;
	$trimmedfile = $root.".wind15_28.$TRIM5"."crop$trimover".".mlen$minlength";	
	$outFQfile2 = $FQreportsDIR."/$trimmedfile\_fastqc\/fastqc_data.txt";
  $outFQsumfile2 = $FQreportsDIR."/$trimmedfile\_fastqc\/summary.txt";
  $outFQzipfile2 = $FQreportsDIR."/$trimmedfile\_fastqc.zip";

	if(!-s "$trimmedfile.fastq")
	{
		open(TRIM1,"$BINPATH/trim_se_window$trimtype.sh $root 15 28 2>&1 |") 
			|| "# cannot run $BINPATH/trim_se_window$trimtype.sh $root 15 28\n";
		while(<TRIM1>)
    { 
      if(/Surviving: \S+ \((\S+)?%/)
      {
        print;
        if($1 < $MINSURVIVALRATE)
        {
          print "# WARNING: less than $MINSURVIVALRATE % reads survided, please check report at $FQreportsDIR\n";
          exit;
        }
      }
    }  
		close(TRIM1); 
	
		$root = $root . '.wind15_28'; push(@trash,$root.'.fastq');
		open(TRIM2,"$BINPATH/trim_se_ends_bases$trimtype.sh $root $TRIM5 $trimover 2>&1 |") 
			|| "# cannot run $BINPATH/trim_se_ends_bases$trimtype.sh $root $TRIM5 $trimover\n";
	        while(<TRIM2>)
          {
            print;
            if(/Surviving: \S+ \((\S+)?%/)
            {
              if($1 < $MINSURVIVALRATE)
              {
                print "# WARNING: less than $MINSURVIVALRATE % reads survided, please check report at $FQreportsDIR\n";
                exit;
              }
            }
          }
	        close(TRIM2);
	
		$root = $root . ".$TRIM5"."crop$trimover"; push(@trash,$root.'.fastq');
		open(TRIM3,"$BINPATH/trim_se_minlen$trimtype.sh $root $minlength 2>&1 |") 
			|| "# cannot run $BINPATH/trim_se_minlen$trimtype.sh $root $minlength\n";
	        while(<TRIM3>)
          {
            if(/Surviving: \S+ \((\S+)?%/)
            {
              print;
              if($1 < $MINSURVIVALRATE)
              {
                print "# WARNING: less than $MINSURVIVALRATE % reads survided, please check report at $FQreportsDIR\n";
                exit;
              }
            }
          }
	        close(TRIM3);
		# output file should be $trimmedfile.fastq
	} 	

	# produce final fastqc report
	if(!-s $outFQfile2)
  {
    open(FQ,"$FASTQCEXE $trimmedfile.fastq -o $FQreportsDIR -t $CPUTHREADS |") 
			|| die "# cannot run $FASTQCEXE $trimmedfile.fastq -o $FQreportsDIR -t $CPUTHREADS\n";
    while(<FQ>){ }
    close(FQ);

		if(-e $outFQzipfile2)
		{
			system("cat $outFQsumfile2");
			unlink($outFQzipfile2);
		}
  }
  
	# find read pairs and correct headers if required
	if(!-s "$trimmedfile.12.fq")
	{
      my $command = "$SPLITPAIRSEXE -i $trimmedfile.fastq -l $trimmedfile.12.fq -s $trimmedfile.s.fq ";
      if($encoding =~ 'Illumina 1.9'){ $command .= ' -R "(^.*?)\s([12]).*" -S "\1/\2"' } 

    	open(SPLIT,"$command 2>&1 |")|| die "# cannot run $command \n";
	    while(<SPLIT>)
	    {
	        print;
          if(/total pairs of reads = (\d+)/)
          {
            if($1 < 250_000)
            {
              print "\n# WARNING: too few clean pairs of reads ($1)\n";
              print "# WARNING: at least 0.25M are required to assemble chloroplasts\n\n";
            }
          }
	    }
	    close(SPLIT);
	    push(@trash,$trimmedfile.'.s.fq',"$trimmedfile.fastq");
	}
  
  # correct bona fide sequencing errors
  if(!-s "$trimmedfile.corr.fq")
  {
      open(MUSKET,"$MUSKETEXE $trimmedfile.12.fq -o $trimmedfile.corr.12.fq -lowercase -inorder -p $CPUTHREADS|")
          || die "# cannot run $MUSKETEXE $trimmedfile.12.fq -o $trimmedfile.corr.12.fq -lowercase -inorder -p $CPUTHREADS\n";
      while(<MUSKET>)
      {
          print;
      }
      close(MUSKET);
      push(@trash,"$trimmedfile.12.fq");
  }

	# extract sample pairs (for instance 100K)
  my $headlines = sprintf("-%d",4*$READSAMPLESIZE);
	system("head $headlines $trimmedfile.corr.12.fq > $trimmedfile.12.100K.fq") if(!-s "$trimmedfile.12.100K.fq"); 

	# estimate insert size and orientation by mapping against close reference
  # index reference if required
  system("$BWAEXE index $refFASTA"); 
  
	open(BWA,"$BWAEXE mem -p -t 4 $refFASTA $trimmedfile.12.100K.fq > $trimmedfile.bwa.sam 2> $trimmedfile.bwa.log |")
		|| die "# cannot run $BWAEXE mem -p -t 4 $refFASTA $trimmedfile.12.100K.fq > $trimmedfile.bwa.sam 2> $trimmedfile.bwa.log\n";
	while(<BWA>){ }
	close(BWA);
	push(@trash,"$trimmedfile.12.fq","$trimmedfile.12.100K.fq","$trimmedfile.bwa.sam");

	my $orientOK = 0;
	open(LOG,"$trimmedfile.bwa.log") || die "# cannot read $trimmedfile.bwa.log\n";
	while(<LOG>)
	{
		if(/# candidate unique pairs for \((\w+), (\w+), (\w+), (\w+)\): \((\d+), (\d+), (\d+), (\d+)\)/)
		{ 
			print;
			my @or = ($1,$2,$3,$4);
			my @n  = ($5,$6,$7,$8);
			$orient = $or[maxindex(\@n)]; #print "mira $orient $1 $2 $3 $4 @n\n";
		}
		elsif(/analyzing insert size distribution for orientation (\w{2})../)
		{ 
			print;
			if($1 eq $orient){ $orientOK = 1 }
		}
		elsif(/\] mean and std\.dev: \((\S+?), (\S+?)\)/)
		{ 
			print;
			if($orientOK){ $orient = "\n# most frequent orientation=$orient insert size=".int($1)."\n\n"; }
		}
	}
	close(LOG);
	push(@trash,"$trimmedfile.bwa.log"); 

	if($orient eq ''){ "\n# most frequent orientation=nd insert size=0\n\n"; }
	print $orient; 

	# compress and store final trimmed reads
	system("gzip $trimmedfile.corr.12.fq");
	system("mv $trimmedfile.corr.12.fq.gz $inpDIR");

    print "# fixed file: $inpDIR/$trimmedfile.corr.12.fq.gz\n";

	# remove unwanted files if required
	unlink(@trash);  
}

###########################

sub maxindex 
{
	my ($array_ref) = @_;
	my $max = -1;
	my $maxindex = -1;
	foreach my $i (0 .. $#{$array_ref}) 
	{
        	if($array_ref->[$i] > $max)
		{ 
			$max = $array_ref->[$i];
            		$maxindex = $i;  
        	}
	}
	return $maxindex;
}
