#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(max);
use File::Basename;
use Getopt::Long;

## Authors: Roy Ronen <rroen@ucsd.edu> and Christina Boucher <cboucher@eng.ucsd.edu> 
## University of San Diego, California.

## This script prepares an input directory for the SEQuel.jar program given 
## a contigs.fa file from assembly and paired-end reads from sequencing.

## Copyright 2012, The Regents of the University of California. 
## All Rights Reserved.

## Permission to use, copy, modify and distribute any part of this
## program for educational, research and non-profit purposes, without fee,
## and without a written agreement is hereby granted, provided that the
## above copyright notice, this paragraph and the following three paragraphs
## appear in all copies.

## Those desiring to incorporate this work into commercial
## products or use for commercial purposes should contact the Technology
## Transfer & Intellectual Property Services, University of California,
## San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910,
## Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu.

## IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
## FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
## INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN
## IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY
## OF SUCH DAMAGE.

## THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY
## OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
## ENHANCEMENTS, OR MODIFICATIONS. THE UNIVERSITY OF CALIFORNIA MAKES NO
## REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR
## EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
## MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF## THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.

## this script alignes 1/2 fasta files to given reference, and outputs alignment stats to STDOUT

###############################################################################
###############################################################################
###############################################################################

# random prefix to tmp file (to avoid clash in multiple runs)
my $random_number = int(rand(1000001)); 
my $tmp_pref = "/tmp/$random_number";

my ($name, $dirs, $suff) = fileparse($0);
my $usage = "\n\tusage: $name -a1 file1.fa [-a2 file2.fa] -r ref.fa\n\n";
my @to_aln;
my ($to_align_f1, $to_align_f2, $ref_seq, $aln, $blat);
my $result = GetOptions ("a1=s" => \$to_align_f1,
	                     "a2=s" => \$to_align_f2,
						 "r=s"  => \$ref_seq,
						 "b=s"    => \$blat,
						);

die $usage if(not defined $ref_seq or not defined $to_align_f1);
$blat = "blat" if(not defined $blat); # default

# add first file to align
my ($f1_short, $dirs1, $suff1) = fileparse($to_align_f1);
$aln = "$tmp_pref.$f1_short.psl";
push @to_aln, $to_align_f1;

# add second file to align if necessary
if(defined $to_align_f2){
	my ($f2_short, $dirs2, $suff2) = fileparse($to_align_f2);
	$aln = $tmp_pref.".".$f1_short."_".$f2_short.".psl";
	push @to_aln, $to_align_f2;
}

######################### align ###########################
foreach my $to_align_f (@to_aln){
	
	if(@ARGV == 3){
		`$blat -fine $ref_seq $to_align_f $aln 1>/dev/null`;
	}
	else{
		`$blat -fine -extendThroughN $ref_seq $to_align_f $aln 1>/dev/null`;
	}
	
	# read PSL output, parse, and print shorthand to screen
	my $curr_score;
	my $name = "";
	my @scores;
	open(PSL, "<", $aln) or die "could not open PSL output. Quitting...\n";
	while(my $l = <PSL>){
		if ($l =~/^psL/ or $l eq "\n"){ next;} # don't write  
		if($l =~ /^---/){ print $l; next; } # write separator line
		my @l_arr = split(/\t/,$l);
		if($l !~ /^[0-9]+/){
			print_blat_line(\@l_arr);
			next;
		}
		else{
			$name = print_blat_line(\@l_arr);
			$curr_score = get_aln_score(\@l_arr);
			push @scores, $curr_score;
		}
	}
	close PSL;
	
	if(@scores == 0){
		print "\nBest Score: no alignment\n\n";	
	}
	else{
		print "\nBest Score $name: " . max(@scores) ."\n\n";
	}
}

system("rm $aln");

###########################################################
#################### Subroutines ##########################
###########################################################

sub print_blat_line{
	my ($line) = @_;	
	for my $i (0 .. (scalar @{$line} - 1)){
			print $line->[$i] if ($i != 13 and $i != 14);
			print "\t" if ($i < (scalar @{$line} - 1) and $i != 13 and $i != 14);
	}
	return $line->[9];	
}

sub get_aln_score{
    my ($bline_r) = @_;
	my $score = 0;
	if(scalar @{$bline_r} < 20){ return 0; }

	# get alignment score
	my $go = 3; # gap open penalty
	my $gsize = max($bline_r->[5],$bline_r->[7]);
	my $gap_size_penalty = 0;
	if($gsize > 0){ $gap_size_penalty = sprintf("%.0f", $go * log($gsize)); } # very mild penalty for gap sizes
	
	$score = $bline_r->[0] + ($bline_r->[2]>>1) - $bline_r->[1] - $go*$bline_r->[4] - $go*$bline_r->[6] - $gap_size_penalty;
	return $score;
}

sub get_aln_score_old{
	my ($line) = @_;
	my $score = 0;
	if(scalar @{$line} < 20){ return 0; }
	else{
		# get alignment score 
		$score = $line->[0] + ($line->[2]>>1) - $line->[1] - $line->[4] -$line->[6];
	}
	return $score;
}

