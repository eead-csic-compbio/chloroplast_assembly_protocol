#!/usr/bin/perl -w
use strict;
use autodie;
use File::Basename; 
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Tools::CodonTable;

# Script to transfer features annotated on a reference GenBank file to another sequence (in FASTA format)
# which is output in GenBank formatr as well.
# Example of use:
# perl _annot_fasta_from_gbk.pl toyreference.gbk Gaz-8_complete.fa Gaz-8_complete.gbk Gaz-8 2> Gaz-8_complete.log

# Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)
#1) Estacion Experimental de Aula Dei-CSIC, Zaragoza, Spain
#2) Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
#3) Fundacion ARAID, Zaragoza, Spain

my $VERBOSE = 0;

my $NCBIBLASTPATH  = '/path/to/ncbi-blast-2.2.XX+/bin'; # please set appropriate path to blast+ binaries

my $BLASTNEXE      = "$NCBIBLASTPATH/blastn -task megablast -word_size 16 -outfmt \"6 std qlen qseq sseq\" ";
my $SHORTBLASTNEXE = "$NCBIBLASTPATH/blastn -task blastn-short -outfmt \"6 std qlen qseq sseq\" ";
my $MAKEBLASTDBEXE = "$NCBIBLASTPATH/makeblastdb -dbtype nucl";

my $TMPPATH        = (-d '/dev/shm/') ? '/dev/shm/' : '/tmp/';
my $MINQLENGHTPERC = 90; # to warn of partially matched features
my $MINQLENGTHIT   = 50; # % to filter out really short blast hits 
my @FEATURES       = qw ( CDS gene mRNA rRNA tRNA );

# Modify as required
my $GBKHEADER = <<'ENDHEADER';
LOCUS       EEEEEEE   LLLLLLL bp    DNA     circular   
DEFINITION  Genus species ecotype EEEEEEE chloroplast, complete genome.
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      chloroplast Genus species (common name)
  ORGANISM  Genus species
            Eukaryota; Viridiplantae; ...
REFERENCE   1  (bases 1 to LLLLLLL)
  AUTHORS   ...
  JOURNAL   Unpublished
REFERENCE   2  (bases 1 to LLLLLLL)
  AUTHORS   ...
  TITLE     Direct Submission
  JOURNAL   Submitted () ...
ENDHEADER

my ($gbkfile,$fnafile,$outgbkfile,$ecotype) = @ARGV;

die "# usage: $0 <reference GenBank file> <FASTA file to be annotated> <out GenBank file> <ecotype>\n" if(!$outgbkfile);

if(!$ecotype){ $ecotype = '_unknown_ecotype' }

my (@rubbish,$feat,$segment,$blastfile,$tmpblastfile,$qlength,$alnlength,$lengthdiff);
my ($value,$n_of_segments,$start,$end,$refstart,$refend,$missing5,$missing3);
my ($coordist,$fh,$tag,$tagref,$header,$dnalength,$subloc,$strand,$frame,$myCodonTable);


## 1) read input GenBank and store selected features
my $ref_array_features = extract_features_from_genbank($gbkfile,@FEATURES);

printf(STDERR "# $0 : total features in $gbkfile : %d\n",scalar(@$ref_array_features));


## 2) map reference features in FASTA/FNA file and get equivalent coordinates

# 2.1) format FNA file for searches
if($VERBOSE){ system("$MAKEBLASTDBEXE -in $fnafile") }
else{ system("$MAKEBLASTDBEXE -in $fnafile 2>&1 > /dev/null")}

push(@rubbish,$fnafile.'.nhr',$fnafile.'.nin',$fnafile.'.nsq');

# 2.2) put sequence in FASTA file in memory
my $outgbkseq = Bio::SeqIO->newFh(-file=>$fnafile)->getline();
$dnalength = $outgbkseq->length();

my $outgbksource = Bio::SeqFeature::Generic->new(
	-start =>1,
	-end   =>$dnalength,
   -primary_tag => 'source',
	-tag => { 
		organism => 'Brachypodium distachyon',
      organelle=> 'plastid:chloroplast',
		mol_type => 'genomic DNA',
      note     => 'type: DNA'}
	); $outgbkseq->add_SeqFeature($outgbksource);


# 2.3) iteratively map gbk features in input fasta file and write to out gbk file
$blastfile    = basename($gbkfile).'.'.basename($fnafile).'.blast';
$tmpblastfile = $TMPPATH.'_tmpBlastOutput.blast';

# 11. The Bacterial, Archaeal and Plant Plastid Code (transl_table=11) 
# http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11
# http://doc.bioperl.org/bioperl-live/Bio/Tools/CodonTable.html    
$myCodonTable  = Bio::Tools::CodonTable->new( -id => 11 ); 

open(BLASTOUT,">$blastfile");

printf(STDERR "# $0 : saving BLAST result in %s\n",$blastfile);

FEAT: foreach $feat (@$ref_array_features)
{
	#print "$feat->{'type'}\n";
	
	# create location object to then add locations, one per segment
	my $location = Bio::Location::Split->new(-splittype=>'join',-strand=>1);
	$n_of_segments	= 0;
	foreach $segment (@{$feat->{'segments'}})
	{
		my ($mindist,%blasthit) = (9999999999999999999); 
		
		#print "$segment->[0]\n";
		#next if($segment->[0] ne '64579'); # debugging
		#next if($segment->[0] ne '104761'); print "$feat->{'type'}\n";
		#next if($segment->[0] ne '125517');
		#next if($segment->[0] ne '71801');
	
		# print segment to be blasted for reference
		print BLASTOUT "$gbkfile\t$segment->[2]:$segment->[0]-$segment->[1]\t$segment->[3]\n";
		
		# try blastn by default
		#die "echo $segment->[3] | $BLASTNEXE -db $fnafile -outfmt \"6 std qlen qseq sseq\" |";
		system("echo $segment->[3] | $BLASTNEXE -db $fnafile > $tmpblastfile");
		
		if(!-s $tmpblastfile) # try blastn-short if that failed, as can happen with very short exons
		{
			system("echo $segment->[3] | $SHORTBLASTNEXE -db $fnafile  > $tmpblastfile");
		}	
		
		open(BLASTN,$tmpblastfile) || die "# cannot read $tmpblastfile\n";
		while(<BLASTN>)
		{
			print BLASTOUT;
			#Query_1	ABR2.2263.4.1841|scaffold1|size135171	100.00	40	0	0	1	40	51226	51265	9e-17	75.0	40	CGTGTAAACGAGTTGCTCTACCGAACTGAGCTATAGCCCT	CGTGTAAACGAGTTGCTCTACCGAACTGAGCTATAGCCCT
			my @data = split(/\t/,$_);
			
			($alnlength,$refstart,$refend,$start,$end,$qlength) = @data[3,6,7,8,9,12];
			
			# check 5' & 3' missing unaligned bases
			$missing5 = $refstart - 1;
			$missing3 = $qlength-$refend;
			
			if($start>$end)
			{ 
				$start += $missing5;
				$end -= $missing3;
				($start,$end) = ($end,$start);
			}
			else
			{
				$start -= $missing5;
				$end += $missing3;
			}
			
			# estimate distance in coordinate space
			$coordist = abs($segment->[0]-$start);
			
			# calculate difference between query length and alignment length to detect only partially matched features
			$lengthdiff = 100*$alnlength/$qlength; 
			if($lengthdiff>100){ $lengthdiff = 100 }
			
			next if($lengthdiff < $MINQLENGTHIT);
			
			#print "$refstart,$refend,$start,$end,$lengthdiff,$missing5,$missing3,$coordist\n"; 
			if($blasthit{$coordist})
			{
				printf(STDERR "# bad luck, two hits with same $coordist : $segment->[2]:$segment->[0]-$segment->[1]\n");
			}
			else
			{ 
				$blasthit{$coordist} = [$refstart,$refend,$start,$end,$lengthdiff,$missing5,$missing3]; 
				if($coordist < $mindist){ $mindist = $coordist }
			}
		}
		close(BLASTN);#print "mindist $mindist\n";
		
		if(!%blasthit)
		{
			printf(STDERR "# no hits, skipt it : $segment->[2]:$segment->[0]-$segment->[1]\n");
			next FEAT;
		}
		
		unlink($tmpblastfile);
		
		
		# work with hit with min distance in coordinate space
		($refstart,$refend,$start,$end,$lengthdiff,$missing5,$missing3) = @{$blasthit{$mindist}}; 
		
		if($lengthdiff < $MINQLENGHTPERC)
		{
			printf(STDERR "# partially matched feature (%1.1f), skip it : $segment->[2]:$segment->[0]-$segment->[1]\n",$lengthdiff);
		}
		elsif($missing5 || $missing3)
		{
			printf(STDERR "# unaligned bases 5'($missing5 ) 3'($missing3), skipt it : $segment->[2]:$segment->[0]-$segment->[1]\n");
		}
		else 
		{
			print "$segment->[0],$segment->[1] $start,$end $lengthdiff\n" if($VERBOSE);
			
			# add location of this segment conserving original order
			$location->add_sub_Location( Bio::Location::Simple->new(
				'-start'  => $start,
				'-end'    => $end,
				'-strand' => $segment->[2])
			); 
			
			$n_of_segments++;
		}
	}
	
	next if(!$n_of_segments);
	
	# actually add this feature to $seq object
	my $outgbkfeature = Bio::SeqFeature::Generic->new(-primary_tag =>$feat->{'type'}); 
		
	# add location to the feature	
	$outgbkfeature->location($location); 
	
	# debug: check gbk location string conserves order of sublocations
	#printf("%s\n",$location->to_FTstring(-no_sort =>1));
	#for my $splitloc ($location->sub_Location(0)){	printf("%d:%s..%s\n",$splitloc->strand,$splitloc->start,$splitloc->end);}	
			
	# add tags to this feature
	foreach $tagref (@{$feat->{'tags'}})
	{
		($tag,$value) = ($tagref->[0],$tagref->[1]);
		
		if($tag eq 'exception' && $value eq 'RNA editing')
		{
			$tag = 'note'; 
			$value = 'RNA editing of first codon annotated on ecotype Bd21, accession EU325680.1';
		}
		
		if($tag eq 'translation')
		{
			## cut DNA sequence of this feature by appending exons in the proper orientation
			my ($featseq,$exonseq,$splitloc,%strands);
			for my $splitloc ($location->sub_Location(0))
			{
				$strand = $splitloc->strand();
				$strands{$strand}++;
			}
			
			# two splicing scenarios:
			# A: cut and append exons first, finally take rc if requested
			if(scalar(keys(%strands)) == 1)
			{
				for $splitloc ($location->sub_Location(0))
				{
					#printf("%d:%s..%s\n",$splitloc->strand,$splitloc->start,$splitloc->end);
					$exonseq = $outgbkseq->subseq($splitloc->start,$splitloc->end);
					$featseq .= $exonseq;
				}	
    			
				if($strand == -1)
				{
					$featseq =~ tr/acgtnACGTN/tgcanTGCAN/;
    				$featseq = reverse($featseq);
				}
			}
			else # B: cut and take rc of each exon and then append
			{		
				#for $splitloc (sort { $a->start <=> $b->start } ($location->sub_Location()))
				for $splitloc ($location->sub_Location(0))
				{
					#printf("%d:%s..%s\n",$splitloc->strand,$splitloc->start,$splitloc->end);
					$exonseq = $outgbkseq->subseq($splitloc->start,$splitloc->end);
    				if($splitloc->strand == -1)
					{
						$exonseq =~ tr/acgtnACGTN/tgcanTGCAN/;
    					$exonseq = reverse($exonseq);
					}	
					$featseq .= $exonseq;
 				}
			}
	
			my $protseq = $myCodonTable->translate($featseq);
			if($myCodonTable->is_start_codon(substr($featseq,0,3)))
			{
				$protseq = 'M'. substr($protseq,1);
			}#print substr($featseq,0,3).' '.substr($protseq,0,1)."\n";
			
			if($protseq =~ /\*\S/)
			{
				my $locstring;
				#foreach $subloc (sort { $a->start <=> $b->start } ($location->sub_Location()))
				foreach $subloc ($location->sub_Location(0))
				{
	 				$locstring .= sprintf("%d:%s..%s",$subloc->strand,$subloc->start,$subloc->end).',';
            }
				printf(STDERR "# translated protein contains internal STOP codon : %s\n",$locstring);
			}
			
			# update $value with newly translated sequence
			$value = $protseq;
		}
			
		$outgbkfeature->add_tag_value($tag,$value);			
	}			
	$outgbkseq->add_SeqFeature($outgbkfeature);
}

close(BLASTOUT);


## 3) write output GenBank file
my $tmpgbkfile = $fnafile;
$tmpgbkfile =~ s/\.\S+$/.tmp.gbk/;

$outgbkfile = $fnafile;
$outgbkfile =~ s/\.\S+$/.gbk/;

# 3.1) produce gbk format straight from bioperl classes
open($fh,">$tmpgbkfile");
my $outgbk = Bio::SeqIO->new(-fh=>$fh,-format=>"genbank");
$outgbk->write_seq($outgbkseq);
close($fh);
push(@rubbish,$tmpgbkfile);

# 3.2) set header with our medatata
my $featOK = 0;

open(GBK,">$outgbkfile");

$header = $GBKHEADER;
$header =~ s/EEEEEEE/$ecotype/g;
$header =~ s/LLLLLLL/$dnalength/g;
print GBK $header;

open(TMPGBK,$tmpgbkfile);
while(<TMPGBK>)
{
	if(/^FEATURES/){ $featOK = 1 }
	print GBK if($featOK);
}
close(TMPGBK);
close(GBK);




## 4) clean
unlink(@rubbish);









sub extract_features_from_genbank
{
	my ($infile,@features2read) = (@_);

	# preferred order of tags for some features
	my %tagOrder = ( 
		'CDS'  => ['gene','codon_start','transl_table','product','translation'],
		'tRNA' => ['gene','product','note'] );
				
	my (@features,$start,$end,$strand,$featype,$tag,$value,$coords,$sequence);

	my $in = new Bio::SeqIO(-file => $infile, -format => 'genbank' );	
	#foreach my $it (keys(%$in)){ print "$it $in->{$it}\n"; }
	
	while( my $seq = $in->next_seq()) # scan next sequence, this is an iterator  
	{    
	   foreach my $f ($seq->get_SeqFeatures) 
	   {
			$featype = $f->primary_tag();
			
			next if(!grep(/^$featype$/,@features2read));

			my (%feature,@seqsegments);        
			
			# declare feature and save type ie CDS
			$feature{'type'} = $featype; #print "$featype\n";  
			
			if($f->location->isa('Bio::Location::SplitLocationI'))
         {
      		for my $location ( $f->location->sub_Location(0) ) 
				{
					$start = $location->start(); # 1..n naturals
            	$end   = $location->end();
            	$strand= $location->strand(); 
					$sequence = $seq->subseq($location->start,$location->end);
					#print "$strand:$start-$end:$sequence\n";
					push(@seqsegments,[$start,$end,$strand,$sequence]);
      		}
			}
			else
			{
				$start = $f->start(); # 1..n naturals
            $end   = $f->end();
            $strand= $f->location()->strand(); 
				$coords = $f->spliced_seq(); 
				$sequence = $coords->seq();#print "$strand:$start-$end:$sequence\n";
				push(@seqsegments,[$start,$end,$strand,$sequence]);
			}

			# store sequence segments linked to this feature
			$feature{'segments'} = \@seqsegments;

			# save tags & values associated to this feature
			#my %saved_tags;
			## save order of relevant tags if possible (no sirve de nada porque luego se imprime sin control nuestro)
			#if($tagOrder{$featype}){
			#	foreach $tag (@{$tagOrder{$featype}}){
			#		if($f->has_tag($tag)){  
			#			for $value ($f->get_tag_values($tag)){ push(@{$feature{'tags'}},[$tag,$value]) }
			#			$saved_tags{$tag} = 1;
			#		}
			#	}		  
			#}
			
			for $tag ($f->get_all_tags) 
			{           
				#next if($saved_tags{$tag}); 
      		for $value ($f->get_tag_values($tag)) 
				{              
					push(@{$feature{'tags'}},[$tag,$value]); #print "$tag $value\n";
      		}          
			}
			
			#foreach $tag (@{$feature{'tags'}}){ print "$featype $tag->[0] $tag->[1]\n";}
			
			# save this feature
			push(@features,\%feature);
		}
    }    
       	
	return (\@features);
}


#http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11
#11. The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
#
#  AAs  FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
#Starts ---M---------------M------------MMMM---------------M------------
#Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
