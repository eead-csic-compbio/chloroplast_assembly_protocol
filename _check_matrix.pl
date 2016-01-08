#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

# Script to analyze DNA polymorphisms along pre-aligned cp genomes 
# Produces data files to be used as tracks with circos

# Bruno Contreras Moreira (1,2)
#1) Estacion Experimental de Aula Dei-CSIC, Zaragoza, Spain
#2) Fundacion ARAID, Zaragoza, Spain

die "# usage: $0 <file with multiple alignment of DNA sequences in FASTA format> " if(!$ARGV[0]);

my $MATRIXFILE = $ARGV[0];

my $REFSTRAIN  = 'ABR6'; # name of reference strain, must be contained in multiple alignment
my $WINDOW     = 100;    # window size for calculating properties along the sequence

my @tracks = qw( N SNPs parsimony_info indels );

my ($feat,$gene,$CDScoords,$strand,$featseq,$exonseq,$strains,$length);
my ($present,$missing,$Ns,$protseq,$IR,$CDSok,$fnalign,$faalign);
my ($gaps,$SNPs,$pars,$align_length,$missense,$ref_tracks);

print "# params: MATRIXFILE=$MATRIXFILE REFSTRAIN=$REFSTRAIN WINDOW=$WINDOW\n";

# calc MSAs and print stats
($Ns,$gaps,$SNPs,$pars,$ref_tracks,$strains,$align_length) = 
  check_DNA_msa_ref($MATRIXFILE,$REFSTRAIN,$WINDOW);

print "\n#strains\tlength\tN\tindels\tSNPs\tparsimony-info\n"; 
print "$strains\t$align_length\t$Ns\t$gaps\t$SNPs\t$pars\n\n";

# save tracks to txt files
my @sorted_bins = sort {$a<=>$b} keys(%{$ref_tracks->{'N'}});
foreach my $track (@tracks)
{
  print "# saving track '$track'\n";
  open(TRACK,">_track$track.txt");
  foreach my $bin (@sorted_bins) 
  {
    #1 0 9999 0.00567
    #1 10000 19999 0.00578
    printf TRACK ("1 %d %d %1.5f\n",
      $bin,$bin+$WINDOW-1,$ref_tracks->{$track}{$bin}/$ref_tracks->{"t$track"} || 0.00000);
  }
  close(TRACK);
}

#################################################
#################################################

sub read_FASTA_file_matrix
{
	# in aligned FASTA format
	# returns a ref to a matrix representing the input aligned sequences, with column 0 containing the headers
	# first valid index (first sequence) is '0'
	
   my ( $infile, $user_strain ) = @_;
   my (@FASTA,$name,$seq,$nt,$letter,$n_of_sequences);
	 my $length = 0;
   my $user_strain_nb = -1;
	 
   $n_of_sequences = -1;
   open(FASTA,"<$infile") || die "# read_FASTA_sequence_matrix: cannot read $infile $!:\n";
   while(<FASTA>)
   {
   	 next if(/^$/);
     next if(/^#/);
     if(/^\>(.*?)[\n\r]/)
     {
      	$n_of_sequences++; # first sequence ID is 0
			  $name = $1; 
			  $nt = 0;
        $FASTA[$n_of_sequences][$nt++] = $name;
        if($name eq $user_strain){ $user_strain_nb = $n_of_sequences }
     }
     elsif($n_of_sequences>-1)
     {
      	$_ =~ s/[\s|\n]//g;
			  $_ =~ s/\·/-/g;
			
			  foreach $letter (split(//,$_))
			  {
       	  $FASTA[$n_of_sequences][$nt++] = uc($letter); 
			  }		
		 }
	}
	close(FASTA); 
	 
	my ($l,%length,); 
	foreach $seq (0 .. $#FASTA)
	{
		$l = scalar(@{$FASTA[$seq]}); #print "$l\n";
		$length{$l} = 1;
	}
	 
  if(scalar(keys(%length)) > 1)
	{
		die "# read_FASTA_sequence_matrix: input aligned sequences have different length\n";
	}
	 
	return (\@FASTA,$n_of_sequences+1,$l-1,$user_strain_nb);
}

# requires sequences in input file to be aligned
# calculates conservation of aligned DNA sequences, ignoring Ns,
# and returns tracks of values of $window size with $refstrain coords
sub check_DNA_msa_ref
{
	my ($infileFastaMSA,$refstrain,$window,$verbose) = @_;
	
	my ($fasta_ref,$N,$L,$refstrain_nb) = read_FASTA_file_matrix($infileFastaMSA,$refstrain);
	$N--;
  if(!defined($fasta_ref)){ die "# check_DNA_msa: failed parsing file $infileFastaMSA\n" }
  else{ print "# check_DNA_msa_ref: $refstrain -> $refstrain_nb\n\n"; }
	
	my ($cons_segment_length,$bin) = (0,0);
	my ($r,$c,$letter,$pars,$obsbases,%freq,%tracks);
	my ($refc,$tgaps,$tSNPs,$tpars,$tNs) = (-1,0,0,0,0);

  $tracks{'N'}{$bin} = $tracks{'SNPs'}{$bin} = $tracks{'parsimony_info'}{$bin} = $tracks{'indels'}{$bin} = 0;

  foreach $c (1 .. $L)
	{
		$pars=0;
		$obsbases='';
		$freq{'A'} = $freq{'C'} = $freq{'G'} = $freq{'T'} = $freq{'-'} = $freq{'N'} = 0; 
		
		foreach $r (0 .. $N) 
		{
			$letter = $fasta_ref->[$r][$c];
			$freq{$letter}++; 

      if($r == $refstrain_nb && grep(/$letter/,('A','C','G','T')))
      {
        $refc++; #print "$refc\n";
      }
		} 
		
		foreach $letter ('A','C','G','T')
		{
			next if($freq{$letter} == 0);
			if($freq{$letter} > 1){ $pars++ }
      $obsbases .= $letter;
		}
		
    # Ns
    $tNs += $freq{'N'};
    $tracks{'N'}{$bin} += $freq{'N'};
    $tracks{'tN'} += $freq{'N'};
    
    # SNPs
    if(length($obsbases)>1)
    { 
      $tSNPs++; 
      $tracks{'SNPs'}{$bin}++;
      $tracks{'tSNPs'}++;
    }
   
		# parsimony positions (al menos 2 letras aparecen al menos 2 veces)
    # http://www.megasoftware.net/mega4/WebHelp/glossary/rh_parsimony_informative_site.htm
		if($pars > 1)
    { 
      $tpars++; 
      $tracks{'parsimony_info'}{$bin}++;
      $tracks{'tparsimony_info'}++;
    }
		
		# columns with gaps
		if($freq{'-'} > 1)
    { 
      $tgaps++; 
      $tracks{'indels'}{$bin}++;
      $tracks{'tindels'}++;
    }
    
    # update bin
    if($refc > 0 && ($refc % $window) == 0)
    { 
      $bin += $window;
      $tracks{'N'}{$bin} = $tracks{'SNPs'}{$bin} = $tracks{'parsimony_info'}{$bin} = $tracks{'indels'}{$bin} = 0;
    }	
	}	
	
	return ($tNs,$tgaps,$tSNPs,$tpars,\%tracks,$N+1,$L);
}  
