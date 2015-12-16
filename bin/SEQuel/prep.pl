#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use File::Spec;
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
## MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
## THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.

#####################################################################
############################# PREP ##################################
#####################################################################
my $usage = usage();

# define inputs & set options
my ($reads1, $reads2, $contigs, $sam);
my ($num_threads, $min_len, $outdir, $bwa, $do_aln, $f) = (3, 0, "prep", "bwa", 1, 0); # defaults
my $h;

GetOptions(
	# INPUT
	'r1=s'    => \$reads1,
	'r2=s'    => \$reads2,
	'c=s'     => \$contigs,
	's=s'     => \$sam,
	# OPTIONS
	'o=s'     => \$outdir,
	'l=i'     => \$min_len,
	't=i'     => \$num_threads,
	'b=s'     => \$bwa,
	'h'       => \$h,
	'f=i'     => \$f, # add f N's at the beginning of each contig prior to alignment
);

die $usage if(defined $h);
die $usage unless( defined $contigs );
die $usage unless( (defined $reads1 and defined $reads2) or (defined $sam) );

$do_aln = 0 if(defined $sam);

# internal parameter settings
my ($small_contigs, $no_reads_contigs, $logfile)  = ("contigsSmall", "contigsNoReads", "prep.log");
my $contigs_flank_name = "contigs.flank.fa";
my ($tot_pe_reads, $tot_pe_aln) = (0, 0);
my ($noX0, , $flt_mult_pair_r1, $flt_mult_pair_r2, $flt_mult_pair_both, $flt_pair_from_r1_dest, $flt_pair_from_r2_dest) = (0, 0, 0, 0, 0, 0);
my %files; my %n_reads_c; # files & stats for paired reads
my %files_single; my %n_reads_c_single; # files & stats for single reads
my %contigs_h;
my ($tot_contigs, $no_reads, $too_small) = (0, 0, 0);
my @isizes;
my ($bwa_aln1_cmd, $bwa_aln2_cmd);
my ($path, $suf, $contigs_full, $contigs_name);

if($do_aln){
	$reads1 = File::Spec->rel2abs($reads1);
    $reads2 = File::Spec->rel2abs($reads2);	
	$contigs_full = File::Spec->rel2abs($contigs);
	($contigs_name, $path, $suf) = fileparse($contigs);
	
	$contigs_name = $contigs_flank_name if($f > 0);
	
	# BWA alignment commands (used to have: -e $max_ge, -o max_go, -L and -d 100)
	my ($max_ge, $max_go, $k, $n) = (100, 3, 3, 0.08);
	$bwa_aln1_cmd  = "$bwa aln $contigs_name $reads1 -t $num_threads -O 7 -E 2 -k $k -n $n -q 15";
	$bwa_aln2_cmd  = "$bwa aln $contigs_name $reads2 -t $num_threads -O 7 -E 2 -k $k -n $n -q 15";
}
else{ 
	$sam = File::Spec->rel2abs($sam);
}

# in case this didn't happen (or was overwritten) above
$contigs_full = File::Spec->rel2abs($contigs);
($contigs_name, $path, $suf) = fileparse($contigs);

mkdir($outdir) or die "[prep.pl]::Error:: could not mkdir $outdir, $!. Quitting...\n";
chdir($outdir) or die "[prep.pl]::Error:: could not chdir $outdir, $!. Quitting...\n";
system("cp $contigs_full .");

#####################################################################
############################# MAIN ##################################
#####################################################################
write_log_startup();
my $contigs_flank = flank_for_aln($contigs_name);
$sam = aln_reads_to_contigs($reads1, $reads2, $contigs_flank) if($do_aln); # contigs with flanking N's (default length f=0)
sep_contigs($contigs_name); # contigs without flanking N's (different if f>0)
disperse_reads_to_SAMs($sam);
close_SAMs_write_stats();
clean_short_and_detect_empty();
write_summary();
rename_contigs_file();

#####################################################################
######################### SUBROUTINES ###############################
#####################################################################
sub usage{
	my ($name, $path, $suf) = fileparse($0);
	my $usage  = "\n";
	   $usage .= "\tUSAGE1: $name -r1 1.fq -r2 2.fq -c contigs.fa [-o OUT-DIR] [-l MIN-LEN] [-t THREADS] [-b full-path-to-bwa]\n";
	   $usage .= "\tOR\n";
	   $usage .= "\tUSAGE2: $name -s aln.sam -c contigs.fa [-o OUT-DIR] [-l MIN-LEN]\n";
	   $usage .= "\n\n";
	   $usage .= "\tNOTES:\n";
	   $usage .= "\t1) USAGE1 requires the BWA v0.6.* executable (either in user's PATH, or passed via the -b argument).\n";
	   $usage .= "\t2) USAGE2 saves time on the 2nd (or higher) run, using the existing SAM file from the 1st run.\n";
	   $usage .= "\t   Also, if alternative alignment methods (other than BWA) are desired, this option can be\n";
	   $usage .= "\t   used with any SAM file of all reads aligned to all contigs.\n";
	   $usage .= "\t3) prep.pl creates a SAM file per contig. Due to runtime considerations, these files are kept\n";
	   $usage .= "\t   open simultaneously. Therefore, the user's maximum open-file-descriptor limit (check with\n";
	   $usage .= "\t   'ulimit -n') must exceed the number of contigs in the assembly.\n";
	   $usage .= "\n";
	   $usage .= "\tFor more details see the SEQuel manual at: http://bix.ucsd.edu/SEQuel/man.html\n";
	   $usage .= "\n";
			
	return $usage;
}

#####################################################################
sub flank_for_aln{
	
	my ($contigs) = @_;

	if($f > 0){
		print "Adding $f flanking N's to each contig...\n";

		# create flank string
		my $flank   = "N" x $f;
		   $flank .= "\n";
		
		# make new contigs file where each sequence is surrounded with the flank string
		open(OLD, "<", $contigs) or die "[prep.pl]::Error:: can not open $contigs, $!. Quitting...\n";
		open(NEW, ">", $contigs_flank_name) or die "[prep.pl]::Error:: can not open $contigs_flank_name, $!. Quitting...\n";
		
		while(my $l = <OLD>){
			if($l =~ /^>/){
				print NEW $l, $flank;
			}
			else{
				print NEW $l;
			}
		}
		print NEW $flank;
		
		close OLD;
		close NEW;

		return $contigs_flank_name; 
	}

	return $contigs;
}

#####################################################################
sub rename_contigs_file{
	# rename contigs-file if conforms to convention of prep for SEQuel
	if($contigs_name =~ /^([0-9]+)\.fa$/){ 
		my $new_name = "contigs_$1.fa";
		system("mv $contigs_name $new_name");
		$contigs_name = $new_name;
	}
	elsif($contigs_name =~ /^([0-9]+)\.fasta$/){
		my $new_name = "contigs_$1.fasta";
		system("mv $contigs_name $new_name");
		$contigs_name = $new_name;
	}
}

#####################################################################
sub write_log_startup{
	my $time = localtime;
	
	open(LOG, ">", $logfile) or die "[prep.pl]::Error:: could not open $logfile, $!. Quitting...\n";

	print LOG "\n";
    print LOG "$time --- Starting prep for SEQuel...\n\n";

	if($do_aln){
		print LOG "\tWill align reads:\n";
		print LOG 	"\t\t$reads1\n";
		print LOG 	"\t\t$reads2\n";
		print LOG "\tTo contigs:\n";
		print LOG 	"\t\t$contigs_full\n";
		print LOG "\n";
		print LOG "\tBWA parameters:\n";
		print LOG "\t\tbwa-aln-1: $bwa_aln1_cmd\n";
		print LOG "\t\tbwa-aln-2: $bwa_aln2_cmd\n";
		print LOG "\n";
	}

	print LOG "\tMinimum contig length: $min_len\n";
	close LOG;

	print "\n$time --- Starting SEQuel prep... \n\n";
}

#####################################################################
sub clean_short_and_detect_empty{
	# assumes contig and sam file names are of the form "ID.fa"
	mkdir($small_contigs) or die "[prep.pl]::Error:: could not mkdir $small_contigs, $!. Quitting...\n";
	mkdir($no_reads_contigs) or die "[prep.pl]::Error:: could not mkdir $no_reads_contigs, $!. Quitting...\n";
	
	# if alignment done, clean left overs to aln dir
    if($do_aln){
		# clean up
		mkdir("aln") or die "[prep.pl]::Error:: could not mkdir aln, $!. Quitting...\n";
		my $contigs_to_mv = $contigs_name;
		$contigs_to_mv = $contigs_flank_name if($f > 0);
		system("mv $contigs_to_mv.* *.sai aln/");
	}
	
	# cleanup short contigs
	my @dir_files = <*>;
	foreach my $file (@dir_files){
		if($file =~ /^([0-9]+)\.fa$/){
			my $id = $1;
			if(exists $contigs_h{$id} ){
				# long enough (check if sam files exist, if not move contig to no-reads directory)
				if(not -e "$id.sam" and not -e "$id.pair.sam" and not -e "$id.single.sam"){
					system("mv $file $no_reads_contigs/");
					$no_reads++;
				}
			}
			else{
				# too short (move contig & sam to small-contigs directory)
				system("mv $file $small_contigs/");
				if(-e "$id.sam"){ system("mv $id.sam $small_contigs/"); $too_small++; }
			}
		}
	}

	# cleanup remaining SAM files
	@dir_files = <*>;
	foreach my $file (@dir_files){
		my $id = sam_to_num($file);
		if($id ne ""){
			if(not exists $contigs_h{$id}){
				# sam file of (too) short contig
				system("mv $file $small_contigs/");
				$too_small++;
			}
		}

	}
}

#####################################################################
sub sam_to_num{
	my ($file) = @_;
	if($file =~ /^([0-9]+)\.sam$/){
		return $1;
	}
	elsif($file =~ /^([0-9]+)\.pair.sam$/){
		return $1;
	}
	elsif($file =~ /^([0-9]+)\.single.sam$/){
		return $1;
	}
	else{
		return "";
	}
}

#####################################################################
sub write_summary{
	
	# get mean insert size from alignment
	my $avg = 0;
	for my $is (@isizes){
		$avg += $is;
	}
	$avg = $avg / ( scalar @isizes );
	$avg = int($avg+0.5);

	# write to log file
	open(LOG, ">>", $logfile) or die "[prep.pl]::Error:: could not open $logfile, $!. Quitting...\n";
	my $time = localtime;
	print LOG "\n$time --- Done! (total contigs: $tot_contigs, 0-reads-aligned contigs: $no_reads, too-short contigs: $too_small)\n\n";
	print LOG "\nSummary:\n";
	print LOG "\t$tot_pe_reads read-pairs\n";
	print LOG "\t$tot_pe_aln aligned read-pairs\n";
	print LOG "\t", $tot_pe_aln - ($flt_mult_pair_r1 + $flt_mult_pair_r2 + $flt_mult_pair_both), " permissively aligned read-pairs\n";
	print LOG "\t\t$flt_mult_pair_r1\tfiltered-pairs on mult-map (only read1 had X0, and X0 > 1)\n";
	print LOG "\t\t$flt_mult_pair_r2\tfiltered-pairs on mult-map (only read2 had X0, and X0 > 1)\n";
	print LOG "\t\t$flt_mult_pair_both\tfiltered-pairs on mult-map (both reads had X0 > 1)\n";
	#print LOG "\t\t$flt_mult_r1\tfiltered read1 on mult-map (both reads had X0, only read1-X0 > 1)\n";
	#print LOG "\t\t$flt_mult_r2\tfiltered read2 on mult-map (both reads had X0, only read2-X0 > 1)\n";
	print LOG "\t\t$noX0\tkept-pairs with no info (both reads missing X0)\n";
	print LOG "\n";
	print LOG "\tAverage external insert-size: $avg\n";
	print LOG "\n";
	close LOG;
	
	print "\n$time --- Done! see $outdir/$logfile for details\n\n";
}

#####################################################################
sub close_SAMs_write_stats{
	open(STATS, ">", "stats.txt") or die "[prep.pl]::Error:: could not open stats.txt, $!. Quitting...\n";
	print STATS "#contig\t#-reads-permissively-aligned (pairs)\t#-reads-permissively-aligned (singles)\n";
	foreach my $c (sort keys %files){
		close $files{$c};
		print STATS "$c\t$n_reads_c{$c}";
		if(defined $n_reads_c_single{$c}){
			print STATS "\t$n_reads_c_single{$c}\n";
		}
		else{
			print STATS "\t0\n";
		}
	}
	close STATS;
}

#####################################################################
sub disperse_reads_to_SAMs{
	my ($sam) = @_;
	open(SAM,"<", $sam) or die "[prep.pl]::Error:: could not open $sam, $!. Quitting...\n";
	print "Processing reads into separate (contig) SAM files...\n\n";
	
	# read SAM header lines
	my $rewind = 0;
	while(my $l = <SAM>){
		$rewind = length $l;
		last if($l !~ /^@/);
	}
	seek(SAM, -1*$rewind , 1);

	# process alignments
	while(my $l1 = <SAM>){
		$tot_pe_reads++;
		print "\t$tot_pe_reads read-pairs processed\n" if($tot_pe_reads % 1_000_000 == 0);
		
		# get 2nd read and split
		my $l2 = <SAM>;
		my @a1 = split(/\t/, $l1);
		my @a2 = split(/\t/, $l2);

		# The following check ensures this SAM is of paired-end reads. 
		# When using single-reads from read-EC is implemented, this check
		# should be removed, or a duplicate "disperse_reads_to_SAMs_single"
		# function should be written to handle single-read-SAMs.

		# sanity check - verify reads are read1 & read2 of a pair
		my ($rname1, $rname2) = ($a1[0], $a2[0]);
		die "\n[prep.pl]::Error: epxecting paired-reads in SAM. See: $rname1 and $rname2\n\n" if($rname1 ne $rname2);
		
		my ($cname1, $cname2) = ($a1[2], $a2[2]);
		my ($bit1,$bit2) = ($a1[1], $a2[1]);
		
		#      Follwoing is permissive-alignment, s.t. more reads would be included.
		#      Only reads explicitely aligned to multiple positions are excluded.
		#      This should include reads that "dangle" off the edges of contigs, and reads where the mate aligned to another contig.
		#      Aditionally, if a read is asigned to a contig the pair will be as well (unless it aligned to multiple positions).
		
		next if( ($cname1 eq "*") and ($cname2 eq "*") );

		$tot_pe_aln++;

		# get X0 field in both reads
		$l1 =~ /X0:i:([0-9]+)/;
		my $X0_1 = $1;
		$l2 =~ /X0:i:([0-9]+)/;
		my $X0_2 = $1;
		
		my ($w_p, $w_s) = (1, 0); # write paired-read, write single-read
		my ($c1, $c2)   = (1, 1); # write to cname1-sam, write to cname2-sam
		
		#      Following is set up for using single (mate-less) reads. This
		#	   happens in error-correction of reads prior to assembly, when an 
		#      entire read is found erroneus and discarded. The mate-pair of such 
		#      a read goes into the assembly process as a single read, and should
		#      be used in SEQuel as well. 
		
		my ($single_cname, $single_samline); # single-read to write - obselete for now!
		
		# determine whether to write read-pair, single-read, or none (and where)
		if(defined $X0_1 and not defined $X0_2){
			# only read1 has X0 info, write pair to c1+c2 unless X0>1
			if($X0_1 > 1){
				$w_p = 0;
				$flt_mult_pair_r1++;
			}
		}
		elsif(defined $X0_2 and not defined $X0_1){
			# only read2 has X0 info, write pair to c1+c2 unless X0>1
			if($X0_2 > 1){
				$w_p = 0;
				$flt_mult_pair_r2++;
			}
		}
		elsif(defined $X0_1 and defined $X0_2){
			# both read1 & read2 have X0 info
			if($X0_1 > 1 and $X0_2 > 1){
				# do not write pair
				$w_p = 0;
				$flt_mult_pair_both++;
			}
			elsif($X0_1 > 1){
				# write pair only to c2
				$c1 = 0;
				$flt_pair_from_r1_dest++;
			}
			elsif($X0_2 > 1){
				# write pair only to c1
				$c2 = 0;
				$flt_pair_from_r2_dest++;
			}
		}
		else{
			# neither read1 or read2 have X0 info, write pair
			$noX0++;
		}

		# write single or pair read
		if($w_p){
			my $F;
			if($cname1 eq $cname2){
				# reads aligned to same contig
				if($cname1 ne "*"){
					# write pair to same contig SAM
					$F = get_file_handle_pair($cname1);
					print $F $l1, $l2;
				}
			}
			else{
				# reads aligned to different contigs
				if($cname1 ne "*" and $c1){
					# write pair to SAM of contig where read1 aligned
					$F = get_file_handle_pair($cname1);
					print $F $l1, $l2;
				}
				if($cname2 ne "*" and $c2){
					# write pair to SAM of contig where read2 aligned
					$F = get_file_handle_pair($cname2);
					print $F $l1, $l2;
				}
			}	
		}
		elsif($w_s){
			# obselete for now!
			my $F = get_file_handle_single($single_cname);
			print $F $single_samline;
		}
		
	}
	close SAM;
	print "\nDone.\n";
}

#####################################################################
sub get_file_handle_single{
	my ($name) = @_;
	my $F;
	if(exists $files_single{$name}){
		# already opened file handle for this contig
		$F = $files_single{$name};
		$n_reads_c_single{$name} += 1;
	}
	else{
		# first time, open it
		my $fname = $name . ".single.sam"; # for Euler contig-names
		my $num = get_num_from_name($name);
		if($num ne ""){
			$fname = $num . ".single.sam"; # for SPAdes/Velvet/SOAP/ALLPATHS contig-names
		}
		open($F, ">", $fname) or die "[prep.pl]::Error:: could not open $fname, $!. Quitting...\n";
		$files_single{$name} = $F;
		$n_reads_c_single{$name} = 1;
	}
	return $F;
}

#####################################################################
sub get_file_handle_pair{	
	my ($name) = @_;
	my $F;
	if(exists $files{$name}){
		# already opened file handle for this contig
		$F = $files{$name};
		$n_reads_c{$name} += 2;
	}
	else{
		# first time, open it
		my $fname = $name . ".pair.sam"; # for Euler contig-names
		my $num = get_num_from_name($name);
		if($num ne ""){
			$fname = $num . ".pair.sam"; # for SPAdes/Velvet/SOAP/ALLPATHS contig-names
		}
		open($F, ">", $fname) or die "[prep.pl]::Error:: could not open $fname, $!. Quitting...\n";
		$files{$name} = $F;
		$n_reads_c{$name} = 2;
	}
	return $F;
}

#####################################################################
sub get_num_from_name{
	my ($name) = @_;
	# extracts & returns contig-number from contig name
	# if none found, returns empty string
	if($name =~ /NODE_([0-9]+)_/){
		# SPAdes/Velvet contig name
		return $1;
	}
	elsif($name =~ /Contig_([0-9]+)/){
		# SOAP contig name
		return $1;
	}
	elsif($name =~ /contig_([0-9]+)/){
		# ALLPATHHS contig name
		return $1;
	} 
	return "";
}

#####################################################################
sub aln_reads_to_contigs{
	my ($reads1, $reads2, $contigs) = @_;
	print "\nGenerating read alignments...\n\n";
	
	# align with BWA
	my ($sai1, $sai2, $sam) = ("reads1_aln.sai", "reads2_aln.sai", "reads_aln.sam");
	my $cmd;
	$cmd = "$bwa index -a is $contigs 2>/dev/null";
	system($cmd);
	
	$cmd = $bwa_aln1_cmd . " > $sai1";
	system($cmd);
	$cmd = $bwa_aln2_cmd . " > $sai2";
	system($cmd);
	
	print "\nGenerating paired-end alignments & estimating insert-size...\n";
	my $is = "isize.txt";
	my $bwa_sampe_cmd = "$bwa sampe $contigs $sai1 $sai2 $reads1 $reads2 > $sam 2>$is";
	system($bwa_sampe_cmd);
	
	open(ISIZE, "<", $is) or die "\n\tCould not open $is, $!. Quitting...\n";
	while (my $line = <ISIZE>){
		chomp $line;
		if($line =~ /inferred external isize/){
			$line =~ /pairs: (\d+\.?\d+)/; #inferred external isize from 99779 pairs: 283.355 +/- 31.672
			push @isizes, $1;
		}
	}
	close ISIZE;
	system("rm $is"); # clean up
	print "\nDone!\n\n";
	return $sam;
}

#####################################################################
sub sep_contigs{
        print "Extracting contigs >= $min_len bp...";
        my ($f) = @_;
        open(CONTS, "<", $f) or die "[prep.pl]::Error:: could not open $f, $!. Quitting...\n";
        while(my $l = <CONTS>){
                if($l =~ /^>/){
						$tot_contigs++;

						# get header
                        my $header = $l;
                        chomp $header;
                        $header =~ s/^>//;
						
						# read contig sequence
						my $hold_sep = $/;
						$/ = ">";
						my $seq = <CONTS>; # read sequence						
						chop $seq;
						my $cont_len = length($seq);
						$/ = $hold_sep;
						seek(CONTS, -1, 1);

						# get contig ID
						my $cont_id;
			            if($header =~ /NODE_/){
							# Velvet
							my @line = split(/_/, $header);
							$cont_id = $line[1];
						}
						elsif($header =~ /Contig_/){
							# SOAP
							my @line = split(/_/, $header);
							$cont_id = $line[1];
						}
						elsif($header =~ /contig_/){
							# ALLPATHS
							my @line = split(/_/, $header);
							$cont_id = $line[1];
						}
						else{
							# Euler
							my @line = split(/\s+/, $header);
							$cont_id = $line[0];
						}
						
                        if($cont_len >= $min_len){
							# long enough (do not put in short-contig directory)
							$contigs_h{$cont_id} = $header;
						}
                        my $file_name = $cont_id.".fa";
                        my $temp = ">".$cont_id."\n"; 
                        
                        open(FA, ">", $file_name) or die "[prep.pl]::Error:: could not open $file_name, $!. Quitting...\n";
						$l .= "\n" if($l !~ /\n$/); # append new line to header if not there (should not happen!!)
						print FA $l, $seq;
                        close FA;
                }
        }
        close CONTS;
        print " done.\n\n";
}

