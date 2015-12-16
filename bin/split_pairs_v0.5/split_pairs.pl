#!/usr/bin/env perl 

use strict;
use Getopt::Std;
use FindBin '$Bin';
use File::Temp qw( tempdir );

my $VERSION = 0.5;

my $KSEQREADEXE = "$Bin/kseqread";
my $PIGZEXE     = "$Bin/pigz";
my $DEFPAIRSEP  = "/"; 

$ENV{'LC_ALL'} = 'POSIX'; # make sure shell sort behaves well

my ($INP_infile,$INP_intlvfile,$INP_pairfile1,$INP_pairfile2,$INP_singlefile); 
my ($INP_pairsep,$INP_nosort,$INP_regex,$INP_subspattern) = ($DEFPAIRSEP,0);
my ($INP_sanger,$INP_tmpdir,$INP_RAMbuffer,$INP_compress) = (0,'/tmp','1024M');

my ($tr_command,$sort_command,$check_command,$command,$tmpdir,%opts,@outfiles) = 
	("$KSEQREADEXE -L",$KSEQREADEXE,$KSEQREADEXE);

my $USAGE=<<EOU; 

USAGE INSTRUCTIONS

-h this message
-i input FASTA/FASTQ filename            (accepts GZIP compressed files)
-l output interleaved filename           (optional, read pairs are printed by default)
-1 output pair1 filename                 (optional, overrides -l)
-2 output pair2 filename                 (optional, overrides -l)
-s output singles filename               (optional, singles are discarded otherwise)
-p char separating pair identifiers      (optional, default: -p "$INP_pairsep", as in HWUSI:4:1101:3600:19#ATCA/1 )
-R Perl regex to shorten headers         (optional, example: -R "(^.*?)#\\w+?([12])" )
-S header substitution pattern           (optional, requires -R and must match -p, example: -S "\\1/\\2")
-H show first header and exit            (optional)
-B RAM buffer (Mb) while sorting         (optional, default: -B $INP_RAMbuffer )
-T path for temporary files              (optional, default: -T $INP_tmpdir )
-n do not sort input reads               (optional, reads are sorted with shell sort by default)
-z GZIP-compress outfiles with N threads (optional, requires -1 & -2 or -l, example: -z 8 )

EXAMPLES

\$ perl split_pairs.pl -i unsorted.fastq.gz -B 1G

\$ perl split_pairs.pl -i sorted.fastq -n -l sorted.12.fq -R "(^.*?)#\\w+?(/[12])" -S "\\1/\\2" 

\$ perl split_pairs.pl -i sorted.fastq -n -1 reads1.fq -2 reads2.fq -s singles.fq -z 8

\$ perl split_pairs.pl -i sorted1.9.fq.gz -p " "

CREDITS

v$VERSION 
Authors: Bruno Contreras-Moreira, Carlos P Cantalapiedra, MaJesus Garcia Pereira, Ruben Sancho Cohen     
EEAD-CSIC 2013-14 (www.eead.csic.es/compbio)

This software is based on http://lh3lh3.users.sourceforge.net/parsefastq.shtml and relies 
on kseq.h by Attractive Chaos <attractor\@live.co.uk> and pigz (http://zlib.net/pigz)
EOU

#-C convert to Sanger/Illumina1.9         (optional)

######################################################################

## 1) take care of arguments	
getopts('hHCnp:i:1:2:l:s:R:S:z:T:B:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) 
{ 
	print $USAGE; exit(-1);
}

if(defined($opts{'n'})){ $INP_nosort = 1 }

if(defined($opts{'i'}))
{
	if(-s $opts{'i'})
	{ 
		$INP_infile = $opts{'i'}; 
		$check_command .= " -i $INP_infile";
		if($INP_nosort){ $sort_command .= " -i $INP_infile"; }
		else{ $tr_command .= " -i $INP_infile"; }
	}
	else
	{
		warn "# ERROR: cannot open input file $INP_infile\n"; 
		exit(-2);
	}
}

if($opts{'R'} && $opts{'S'})
{ 
	$INP_regex = $opts{'R'};
	$INP_subspattern = $opts{'S'};
	$sort_command .= " -R \"$INP_regex\" -S \"$INP_subspattern\"";
	$check_command .= " -R \"$INP_regex\" -S \"$INP_subspattern\"";
}

if(defined($opts{'p'}))
{ 
	$INP_pairsep = $opts{'p'};
	$sort_command .= " -p \"$INP_pairsep\"";
	$check_command .= " -p \"$INP_pairsep\"";
}

open(KSEQREAD,"$check_command -H 2>&1 |") || die "# ERROR: $check_command -H failed\n";
while(<KSEQREAD>)
{
	if($opts{'H'})
	{
		print;
	} 
	elsif(/WARNING/)
	{
		print;
		exit(-3);
	}
	elsif(/first quality encoding/)
	{
		print;
	}
}
close(KSEQREAD); 

if(defined($opts{'H'})){ exit(0) }

if($opts{'s'})
{ 
	$INP_singlefile = $opts{'s'};
	$sort_command .= " -s $INP_singlefile"; 
	
	@outfiles = ( $INP_singlefile );
}

if($opts{'1'} || $opts{'2'})
{
	if($opts{'1'})
	{ 
		$INP_pairfile1 = $opts{'1'}; 
		$sort_command .= " -1 $INP_pairfile1";
	}
	else
	{
		warn "# ERROR: need a valid -1 pair filename\n";
		exit(-4);
	}
   
	if($opts{'2'})
	{ 
		$INP_pairfile2 = $opts{'2'};
		$sort_command .= " -2 $INP_pairfile2";
	}
	else
	{
		warn "# ERROR: need a valid -2 pair filename\n";
		exit(-5);
	}
	
	if(defined($opts{'z'}) && $opts{'z'}>=0){ $INP_compress = $opts{'z'}; }
	
	push(@outfiles,$INP_pairfile1,$INP_pairfile2);
}
elsif($opts{'l'})
{ 
	$INP_intlvfile = $opts{'l'}; 
	$sort_command .= " -l $INP_intlvfile"; 
	
	if(defined($opts{'z'}) && $opts{'z'}>=0){ $INP_compress = $opts{'z'} }
	
	push(@outfiles,$INP_intlvfile);
}
else
{
	$INP_intlvfile = 'stdout'; 
	$sort_command .= " -l $INP_intlvfile"; 
} 

if(defined($opts{'T'}))
{ 
	if(!-d $opts{'T'})
	{
		warn "# ERROR: need a valid directory for temporary files, using default $INP_tmpdir\n";
	}
	else{ $INP_tmpdir = $opts{'T'} }
}

if(defined($opts{'B'})){ $INP_RAMbuffer = $opts{'B'} }

#if(defined($opts{'C'})){ $INP_sanger = 1 }

	 
	 
## 2) split read pairs with selected arguments

if($INP_nosort)
{
	$command = $sort_command;
}
else
{
	$tmpdir = tempdir( '_split_pairs_XXXXXXXX',DIR=>$INP_tmpdir,CLEANUP=>1);
	$command = "$tr_command | sort -k 1,1 -t'\t' -u -S $INP_RAMbuffer -T $tmpdir | ". 
		'tr "\t" "\n" | '. $sort_command;
} #die "$command\n\n";

open(KSEQREAD,$command.' |'); # || die "# ERROR: $command failed, probably sort ran out of space at $INP_tmpdir\n";
while(<KSEQREAD>)
{
	print;
}
close(KSEQREAD);

if(defined($INP_compress))
{
	print "# compressing outfiles...\n";
	if($INP_compress){ $INP_compress = "-p $INP_compress" }
	else{ $INP_compress = '' }

	foreach my $outfile (@outfiles)
	{
		print "# $outfile\n";
		next if(!-e $outfile);
		system("$PIGZEXE $INP_compress $outfile");
	}
	print "# done\n\n";
}
