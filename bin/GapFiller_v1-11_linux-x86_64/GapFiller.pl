#!/usr/bin/perl --
#GapFiller v1.11 Marten Boetzer - Walter Pirovano, July 2012

#AUTHORS
# Marten Boetzer and Walter Pirovano (c) 2011-2012

#CITATION
#If you use GapFiller in a scientific publication, please cite:
# Boetzer, M. and Pirovano, W., Toward almost closed genomes with GapFiller, Genome Biology, 13(6), 2012
# walter.pirovano@baseclear.com

#LICENSE
#   GapFiller Copyright (c) 2011-2012 BaseClear B.V. The Netherlands. All rights reserved.
#   GapFiller can be freely used by academic institutes or non-profit organizations.
#   Commercial parties need to acquire a license. For more information about commercial licenses look at our website or email info@baseclear.com.

#TERMS
#   This software comes without any warranty. BaseClear cannot guarantee that this software will fulfill any of your particular purposes or needs.
#   BaseClear will not be liable for data loss or any damages related to the software.
#   This software and eventual license herein granted shall not be copied, shared or offered for re-sale without the explicit approval of BaseClear.

#DOCUMENTATION
#   README, MANUAL and TUTORIAL distributed with this software @ www.baseclear.com
#   Boetzer, M. and Pirovano, W., Toward almost closed genomes with GapFiller, Genome Biology, 13(6), 2012
#   http://www.baseclear.com/sequencing/data-analysis/bioinformatics-tools/
#   We hope this code is useful to you -- Please send comments & suggestions to Walter.Pirovano@baseclear.com
#   If you use either the GapFiller code or ideas, please cite our work appropriately and accurately


###CHANGES in version 1.10:
#  Fixed a bug where internal sequences were reverse complemented

###CHANGES in version 1.9:
#  Included the option -S to skip reading of the files if it was already done in a previous analysis
#  Hashing of the reads is now done during filling the gaps, instead of during the mapping stage. This saves memory and time.

  use strict;
  use File::Basename;
  use Storable;
  use File::Path;
  use File::Copy;
  use FindBin qw($Bin);
  use threads;
  use Text::Wrap;

  use threads::shared;
  my $totalReadsProcessed :shared;
  my $totalReadFiles :shared;

  $Text::Wrap::columns = 61;
   require "getopts.pl";

  use vars qw($opt_m $opt_o $opt_v $opt_p $opt_k $opt_a $opt_z $opt_s $opt_b $opt_n $opt_l $opt_x $opt_u $opt_t $opt_T $opt_g $opt_r $opt_d $opt_S $opt_i);
  &Getopts('m:o:v:p:k:a:z:s:b:n:l:x:u:t:T:g:r:d:S:i:');
   my ($keep, $trim, $gaps, $threads, $difference,$base_overlap,$min_overlap,$min_base_ratio,$base_name, $min_tig_overlap, $numiteration)= (0, 10,1,1, 50, 2,29, 0.7, "standard_output", 10,10);

  my $seplines = ("-" x 60)."\n";
  my $version = "[GapFiller_v1-11_Final]";
#-------------------------------------------------READ OPTIONS

  if(!($opt_l) || !($opt_s)){
     print "ERROR: Parameter -l is required. Please insert a library file\n" if(!$opt_l);
     print "ERROR: Parameter -s is required. Please insert a scaffold fastA file\n" if(!$opt_s);
     print "\nUsage: $0 $version\n\n";

     print "============ General Parameters ============\n";
     print "-l  Library file containing two paired-read files with insert size, error and orientation indication.\n";
     print "-s  Fasta file containing scaffold sequences used for extension.\n";
     print "============ Extension Parameters ============\n";
     print "-m  Minimum number of overlapping bases with the edge of the gap (default -m $min_overlap)\n";
     print "-o  Minimum number of reads needed to call a base during an extension (default -o $base_overlap)\n";
     print "-r  Percentage of reads that should have a single nucleotide extension in order to close a gap in a scaffold (Default: $min_base_ratio)\n";
     print "-d  Maximum difference between the gapsize and the number of gapclosed nucleotides. Extension is stopped if it matches this parameter + gap size (default -d $difference, optional).\n";
     print "-n  Minimum overlap required between contigs to merge adjacent sequences in a scaffold (default -n $min_tig_overlap, optional)\n";
     print "-t  Number of reads to trim off the start and begin of the sequence (usually missambled/low-coverage reads) (default -t $trim, optional)\n";
     print "-i  Number of iterations to fill the gaps (default -i $numiteration, optional)\n";
     print "============ Bowtie Parameters ============\n";
     print "-g  Maximum number of allowed gaps during mapping with Bowtie. Corresponds to the -v option in Bowtie. (default -g $gaps, optional)\n";
     print "============ Additional Parameters ============\n";
     print "-T  Number of threads to run (default -T $threads)\n";
     print "-S  Skip reading of the input files again\n";
     die "-b Base name for your output files (optional)\n";

  }
  my $scaffold = $opt_s if($opt_s);
  my $libraryfile;
  $libraryfile = $opt_l if ($opt_l);
  $min_overlap = $opt_m if ($opt_m);
  $base_overlap = $opt_o if ($opt_o);
  $base_name = $opt_b if($opt_b);
  $min_tig_overlap = $opt_n if($opt_n);
  $min_base_ratio = $opt_r if ($opt_r);
  $threads = $opt_T if ($opt_T);
  $gaps = $opt_g if ($opt_g || $opt_g eq 0);
  $difference = $opt_d if($opt_d);
  $trim = $opt_t if($opt_t || $opt_t eq 0);
  $keep = $opt_S if($opt_S || $opt_S eq 0);
  $numiteration = $opt_i if($opt_i);
  mkpath("$base_name");
  mkpath("$base_name/reads");

  my $files = "";
  my $textline = "Your inserted inputs on $version at ".getDate().":\n\t-s $scaffold\n\t-l $libraryfile\n\t-b $base_name\n\t-o $base_overlap\n\t-m $min_overlap\n";
  $textline .= "\t-r $min_base_ratio\n\t-n $min_tig_overlap\n\t-T $threads\n\t-g $gaps\n\t-d $difference\n\t-t $trim\n\t-i $numiteration\n\n";

  my $summaryfile = "$base_name/$base_name.summaryfile.final.txt";

  open (SUMFILE, ">$summaryfile") || die "Can't open $summaryfile -- fatal\n";
  print SUMFILE "$textline";
  print "$textline";
  &printMessage("\n=>".getDate().": Reading and processing paired-read files\n");

  #process the library and read the input files with multiple threads
  open(FILELIB, "< $libraryfile") || die "Can't open $libraryfile -- fatal\n";
  my ($prevlib,$sequenceReadsInGaps,$filehash, $maxGapScaf, $gapPosInScafHash);
  my $ctlib = 0;
  my $totalfiles = 0;
  while(<FILELIB>){
    my ($lib, $aligner, $fileA, $fileB, $insert_size, $insert_stdev, $orient) = split(/\s+/, $_);
    next if($lib eq '');
    die "ERROR: Invalid aligner in library $lib: $aligner. Should be either 'bowtie', 'bwa' or 'bwasw' -- fatal\n" if($aligner ne "bwa" && $aligner ne "bwasw" && $aligner ne "bowtie");
    die "ERROR: Invalid file in library $lib: $fileA -- fatal\n" if(!(-e $fileA));
    die "ERROR: Invalid file in library $lib: $fileB -- fatal\n" if(!(-e $fileB));
    die "ERROR: Orientation must have length of 2 characters and should contain one of the following; FR, FF, FR or RF. Your library $lib has orientation of $orient ...Exiting.\n" if(!(length($orient) == 2) || !($orient =~ /[FR][FR]/));

    $ctlib++;
    $prevlib = $lib;
    if(!$keep){
      my $thr = threads->create(\&generateInputFiles, $ctlib, $fileA, $fileB, "$aligner.reads.lib$ctlib", $orient);
      if(!($ctlib % $threads)){
        foreach my $thr (threads->list()) {
          my @res = $thr->join();
          my ($libs,$reads) = split(":",$res[0]);
          my @readarray = split(",",$reads);
          foreach my $readf (@readarray){
            $filehash->{$libs}{$readf}++;
          }
        }
      }
    }else{
      my @readfiles = <$base_name/reads/$base_name.$aligner.reads.lib$ctlib*>;
      if($#readfiles < 0){
        die "ERROR: no readfiles found for library $lib line:\n$_\nYou have inserted option -S 1, change this to -S 0\n";
      }
      foreach my $readfile (@readfiles){
        $filehash->{$ctlib}{$readfile}++;
      }
    }
  }
  if(!$keep){
    foreach my $thr (threads->list()) {
      my @res = $thr->join();
      my ($libs,$reads) = split(":",$res[0]);
      my @readarray = split(",",$reads);
      foreach my $readf (@readarray){
        $filehash->{$libs}{$readf}++;
      }
    }
  }
  close FILELIB;

###get the total number of files of all libraries
  foreach my $libs(keys %$filehash){
    my $filelist = $filehash->{$libs};
    foreach my $readfiles(keys %$filelist){
      $totalReadFiles++;
    }
  }
  #exit;
  print "\n";
  my $prevscaffold = $scaffold;
  my $iteration = 0;
  my ($finalFile, $finalSummary, $finalClose, $finalFill);
  my ($prevFillFile,$prevCloseFile);
  mkpath("$base_name/intermediate_results");
  while(++$iteration <= $numiteration){
    $totalReadsProcessed=0;
    $scaffold = $prevscaffold;
    print SUMFILE "After iteration $iteration:\n$seplines\n\n";
    print "ITERATION $iteration:\n";
    my $closefile = "$base_name/intermediate_results/$base_name.closed.evidence.iteration$iteration.txt";
    my $fillFile = "$base_name/intermediate_results/$base_name.filled.iteration$iteration.txt";

    open (CLOSEFILE, ">$closefile") || die "Can't open $closefile -- fatal\n";
    $ctlib=0;
    open(FILELIB, "< $libraryfile") || die "Can't open $libraryfile -- fatal\n";
    &printMessage("\n=>".getDate().": Mapping reads to scaffolds, reading alignment output and storing reads\n");
    my $ctthreads=0;
  #read the libraries again and map the reads to the scaffolds

    system("rm -rf $base_name/alignoutput/*");  #remove previous folder, if exists
    mkpath("$base_name/alignoutput");
    $prevlib = "";
    my ($bwa,$bowtie) = (0,0);
    while(<FILELIB>){
      my ($lib, $aligner, $fileA, $fileB, $insert_size, $insert_stdev, $reverse) = split(/\s+/, $_);
      next if($lib eq '' || $prevlib eq $lib || $fileA eq "TAB");
      $ctlib++;
      $prevlib = $lib;
      my $min_allowed = ($insert_stdev * $insert_size);
      my $maxlib = int($insert_size + $min_allowed);
      my $minlib = int($insert_size - $min_allowed);
      my $contigFile = "$base_name/alignoutput/$base_name.gapclosure.fa";
      processScaffold($scaffold, $contigFile, $trim, $maxlib) if($ctlib == 1);
      my $library = "$lib.$ctlib";
      die "Scaffold file ($contigFile) not found. Exiting...\n" if(!(-e $contigFile));

#BWA!
      if($aligner eq "bwa" || $aligner eq "bwasw"){
        my @readfiles= <$base_name/reads/$base_name.$aligner.reads.lib$ctlib*>;
        my $bwapath = "$Bin"."/bwa/bwa";
        $bwapath =~ s/ /\\ /g;
        my $bwaout = $base_name . ".bwaIndex";
        my $filesize = -s "$contigFile";
        my $index = "bwtsw";
        $index = "is" if($filesize <= 10000000);
        open(STDERR, ">$base_name/alignoutput/tmpbwa_logfile");
        system("$bwapath index -a $index $contigFile -p $base_name/alignoutput/$bwaout") == 0 || die "\nBwa error; $?" if(!$bwa); # returns exit status values
        my $filenum=1;
        foreach my $readfile (@readfiles){
          if($aligner eq "bwa"){
            my $bwaoutputaln = "$base_name/alignoutput/$base_name.$library.$ctlib.$filenum.bwa";
            my $procline = "$bwapath samse $base_name/alignoutput/$bwaout $bwaoutputaln $readfile |";
            my $samseline = "$bwapath aln -i 0 $base_name/alignoutput/$bwaout $readfile > $bwaoutputaln";
            my $thr = threads->create(\&getUnmappedReads, $minlib, $maxlib, $min_allowed, $procline, $samseline);
          }else{
            my $procline = "$bwapath bwasw $base_name/alignoutput/$bwaout $readfile |";
            my $thr = threads->create(\&getUnmappedReads, $minlib, $maxlib, $min_allowed, $procline);
          }
          $filenum++;
          getUnmappedThreadResult() if(!(++$ctthreads % $threads));
        }
        $bwa=1;
      }
#bowtie!
      if($aligner eq "bowtie"){
        my @readfiles = <$base_name/reads/$base_name.$aligner.reads.lib$ctlib*>;
        my $bowtieout = "$base_name.bowtieIndex";
        my $bowtiepath = "$Bin"."/bowtie/bowtie";
        $bowtiepath =~ s/ /\\ /g;
        my $bowbuildpath = $bowtiepath."-build";
        system("$bowbuildpath $contigFile $base_name/alignoutput/$bowtieout --quiet --noref") == 0 || die "\nBowtie-build error; $?" if(!$bowtie); # returns exit status values
        foreach my $readfile (@readfiles){
          my $procline = "$bowtiepath -v $gaps $base_name/alignoutput/$bowtieout -f $readfile -S --sam-nohead --quiet |";
          my $thr = threads->create(\&getUnmappedReads, $minlib,$maxlib, $min_allowed, $procline);
          getUnmappedThreadResult() if(!(++$ctthreads % $threads));
        }
        $bowtie = 1;
      }
      $finalSummary = $summaryfile;
      $finalClose = $closefile;
      $finalFill = $fillFile;
    }
    close FILELIB;
    getUnmappedThreadResult();
    CounterPrint("                                                     ");
  #=cut till this part should be included with the 'steps'

    $finalFile = "$base_name/intermediate_results/$base_name.gapfilled.iteration$iteration.fa";
    my $preNcount = goThroughScaffold($scaffold, $finalFile, $fillFile, $prevFillFile,$prevCloseFile);

    my ($sumorig,$preNcountOrig) = writesummaryfiles($scaffold);
    my ($sumclosed,$afterNcount) = writesummaryfiles($finalFile);
    
    $preNcount = $preNcountOrig if($iteration > 1);
    my $nucClosed = ($preNcount - $afterNcount);

    print SUMFILE "\tClosed $nucClosed out of $preNcount nucleotides\n\n$seplines\n";
    print "\tClosed $nucClosed out of $preNcount nucleotides\n\n";

    print SUMFILE "\n$seplines Scaffold statistics:\n$seplines";
    print SUMFILE "\n\tOriginal:\n$sumorig";

    print SUMFILE "\tAfter gapclosing:\n$sumclosed$seplines\n";
    $prevscaffold = $finalFile;
    close CLOSEFILE;
    $trim = 0;

    if($preNcount == $nucClosed){
      print "\tAll gaps are closed. Exiting...\n\n";
      $iteration++;
      last;
    }
    if($nucClosed == 0){
      print "\tCould not close any more gaps. Exiting...\n\n";
      $iteration++;
      last;
    }
    $prevFillFile = $fillFile;
    $prevCloseFile = $closefile
  }
  --$iteration;
  copy("$base_name/intermediate_results/$base_name.closed.evidence.iteration$iteration.txt", "$base_name/$base_name.closed.evidence.final.txt");
  copy("$base_name/intermediate_results/$base_name.filled.iteration$iteration.txt", "$base_name/$base_name.filled.final.txt");
  copy("$base_name/intermediate_results/$base_name.gapfilled.iteration$iteration.fa", "$base_name/$base_name.gapfilled.final.fa");




  my $time = (time - $^T);
  my $minutes = int ($time / 60);
  $time = $time % 60;
  print "Process run succesfully on ".getDate()." in $minutes"." minutes and $time"." seconds\n\n\n";
  print SUMFILE "\nProcess run succesfully on ".getDate()." in $minutes"." minutes and $time"." seconds\n\n\n";
  close SUMFILE;
  close CLOSEFILE;

#####################################    END OF MAIN SCRIPT

#process the read files and cut it into subsets of 1 million pairs
sub generateInputFiles{
  my ($lib,$fileA, $fileB, $outfile, $orient, $fname) = @_;
  my ($name,$seq1,$seq2);
  my ($ctstep,$counterext, $Ncount, $countsinglet, $fastq, $step) = (1,0,0,0,0,1000000);

  my ($ori_1,$ori_2) = split(//, $orient);
  open (OUTSINGLEFILE, ">$base_name/reads/$base_name.$outfile.1") || die "Can't write to single file $base_name/reads/$base_name.$outfile.1 -- fatal\n";
  my $files="$lib:$base_name/reads/$base_name.$outfile.1";
  if ($fileA =~ /\.gz$/) {
    open(TEST, "gunzip -c  $fileA |") || die "can't open .gz $fileB with gunzip";
  }else{
    open(TEST, "< $fileA");
  }
  $name = <TEST>;
  close TEST;
  $fastq = 1 if ($name =~ /^[@]/);
  if ($fileA =~ /\.gz$/) {
    open(FILEA, "gunzip -c  $fileA |") || die "can't open .gz $fileA with gunzip";
  }else{
    open(FILEA, "< $fileA");
  }
  if ($fileB =~ /\.gz$/) {
    open(FILEB, "gunzip -c  $fileB |") || die "can't open .gz $fileB with gunzip";
  }else{
    open(FILEB, "< $fileB");
  }
  while(<FILEA>) {
    <FILEB>;
    $seq1 = uc(<FILEA>); chomp $seq1;
    $seq2 = uc(<FILEB>); chomp $seq2;
    #FASTQ FORMAT
    <FILEA>,<FILEA>,<FILEB>,<FILEB> if ($fastq);

    if(++$countsinglet == $step){
      CounterPrint("reading files of library line $lib @ $countsinglet");
      $step = $step + 1000000;
      close OUTSINGLEFILE;
      $ctstep++;
      open (OUTSINGLEFILE, ">$base_name/reads/$base_name.$outfile.$ctstep") || die "Can't write to single file $base_name/reads/$base_name.$outfile.$ctstep -- fatal\n";
      $files.=",$base_name/reads/$base_name.$outfile.$ctstep";
    }    
    #reverse sequences, if not in --> <-- orientation
    $seq1 = reverseComplement($seq1) if($ori_1 eq "R");
    $seq2 = reverseComplement($seq2) if($ori_2 eq "F");
    print OUTSINGLEFILE ">read$countsinglet.1\n$seq1\n>read$countsinglet.2\n$seq2\n";
  }
  CounterPrint("                                            ");
  close OUTSINGLEFILE;
  close FILEB;
  close FILEA;
  return $files;
}

#process the scaffold file, retreive only the edges of the sequences surrounding the gap.
sub processScaffold{
  my ($scaffile,$processedFile, $trim, $maxlib) = @_;
  my ($seq, $name, $scafct) = ("","",0,0);
  open(SCAF,"$scaffile") || die "Can't open $scaffile -- fatal\n";
  open (BOWIN, ">$processedFile") || die "Can't write to file BOWIN -- fatal\n";

  while(<SCAF>){
    chomp;
    my $line = $_;
    $seq.= uc($line) if(eof(SCAF));
    if (/\>(\S+)/ || eof(SCAF)){
      if($seq ne ''){
        $scafct++;
        my @seqgaps = split(/[N]{1,}/, $seq);
        if($#seqgaps > 0){
          my @countgapper;
          my $gapct=1;
          #go through each gap, and estimate the length and the position along the scaffold. Add extra N's if trimming (-t) is set
          while ($seq =~ /(N+)/g) {
            my $start = (pos($seq)-length($1))-$trim;
            my $end = pos($seq)+$trim;
            push @countgapper, length($1)+($trim*2);
            $gapPosInScafHash->{"scaf$scafct"}{$gapct}{'start'} =$start-$trim;
            $gapPosInScafHash->{"scaf$scafct"}{$gapct}{'end'}= $end+$trim;
            $gapct++;
          }
          my ($ctgcount,$fullseq,$start,$end)=(0,"",0,0);
          my $numgap = $#seqgaps;
          #go threagh each 'contig' of the scaffold, trim the bases and estimate the position of the contig along the scaffold. Of the trimmed contig, only take the edges, to save time during alignment
          foreach my $ctgseq (@seqgaps){
            $ctgseq = substr($ctgseq,$trim) if($ctgcount > 0);
            $ctgseq = substr($ctgseq,0,length($ctgseq)-$trim);
            my $gapsize = $countgapper[$ctgcount];
            $end = $start + length($ctgseq);
            $ctgcount++;
            my $len = length($ctgseq);
            if($len > (($maxlib * 2))){
              my $upper = (length($ctgseq) - ($maxlib));
              my $first = substr($ctgseq, 0, $maxlib);
              my $second = substr($ctgseq, $upper);
              my $newseq = $first.$second;
              $newseq = $first if($ctgcount > $numgap);
              print BOWIN ">scaf$scafct"."_$start"."_$end"."_sub\n";
              print BOWIN "$newseq\n";
            }else{
              $ctgseq = substr($ctgseq, 0, $maxlib) if($ctgcount > $numgap);
              print BOWIN ">scaf$scafct"."_$start"."_$end\n";
              print BOWIN "$ctgseq\n";
            }
            $start = $end + $gapsize;
           }
        }
      }
      $seq='';
      $name = $1;
    }else{
      $seq.=uc($line);
    }
  }
  close SCAF;
  close BOWIN;
}

#get reads that are not mapped within a given distance (insertsize + (insertsize+deviation))
#store all reads per gap in a hash for later extension
sub getUnmappedReads{
  my ($minlib,$maxlib, $min_allowed, $procline, $aligninput) = @_;
  system("$aligninput") if($aligninput ne "");
  open(IN, "$procline") || die "Can't open baw output: Process: $procline -- fatal\n";
  my $usereadhash;

  #go through the SAM file, and determine the actual position along the scaffold, since only edges of the contigs are used instead of whole scaffold.
  while(my $line1 = <IN>){
    next if($line1 =~ /^@/);
    my $line2 = <IN>;
    my @t1 = split(/\t/,$line1);
    my @t2 = split(/\t/,$line2);
    #if first read is aligned, calculate actual position, and in which Gap(s) the other read potentially can fall
    if($t1[2] ne "*"){
      my @p1 = split(/_/,$t1[2]);
      if($p1[3] eq "sub"){
         if($t1[3] >= $maxlib){
           $t1[3] = $p1[2]-($maxlib - ($t1[3]-$maxlib));
         }else{
           $t1[3] = $t1[3] + $p1[1];
         }
      }else{
        $t1[3] = $t1[3] + $p1[1];
      }
      #determine in which gap(s) the other read can potentially fall
      $usereadhash = estimateReadInGap($usereadhash,$p1[0],$t1[1],$t1[3],$t2[9], $minlib,$maxlib, $min_allowed,$t1[9]);
    }
    #if second read is aligned, calculate actual position, and in which Gap(s) the other read potentially can fall
    if($t2[2] ne "*"){
      my @p2 = split(/_/,$t2[2]);
      if($p2[3] eq "sub"){
         if($t2[3] >= $maxlib){
           $t2[3] = $p2[2]-($maxlib - ($t2[3]-$maxlib));
         }else{
           $t2[3] = $t2[3] + $p2[1];
         }
      }else{
        $t2[3] = $t2[3] + $p2[1];
      }
      $usereadhash = estimateReadInGap($usereadhash,$p2[0],$t2[1],$t2[3],$t1[9], $minlib,$maxlib, $min_allowed,$t2[9]);
    }
  }
  close IN;
  $totalReadsProcessed++;
  my $perc = sprintf("%.1f", ($totalReadsProcessed/$totalReadFiles)*100);
  CounterPrint("...Processed ".($totalReadsProcessed*1000000)." paired-reads (~$perc%)");
  return $usereadhash;
}

sub estimateReadInGap{
  my ($usereadhash, $scaf,$ori,$pos,$seq,$mindist,$maxdist,$min_allowed) = @_;
  my $ctglist = $gapPosInScafHash->{$scaf};
  my $strand = (($ori&0x0010))? 'R':'F';
  #if strand is -->, search for gaps right of the position
  if($strand eq "F"){
    my $regionstart = $pos + $mindist;
    my $regionend = $pos + $maxdist;
    foreach my $gapnum (sort {$a<=>$b} keys %$ctglist){
      next if($regionstart > ($ctglist->{$gapnum}{'end'}+$min_allowed));
      last if($regionend < ($ctglist->{$gapnum}{'start'}-$min_allowed));
      $usereadhash->{"$scaf.$gapnum"}{$seq}++;
    }
  }else{ #if strand is <--, search for gaps left of the position
    my $regionstart2 = $pos - $maxdist;
    my $regionend2 = $pos - $mindist;
  
    foreach my $gapnum (sort {$a<=>$b} keys %$ctglist){
      next if($regionstart2 > ($ctglist->{$gapnum}{'end'}+$min_allowed));
      last if($regionend2 < ($ctglist->{$gapnum}{'start'}-$min_allowed));
      $usereadhash->{"$scaf.$gapnum"}{$seq}++;
    }
  }
  return $usereadhash;
}

#get the unmapped reads of all the threads and store it in a new hash, per gap. Go through the collection, and save all reads into a file (OUT) per gap
sub getUnmappedThreadResult{
  foreach my $thr (threads->list()) {
    my @ret = $thr->join();
    my $unmaphash = $ret[0];
    foreach my $scaf (keys %$unmaphash){
      my $list=$unmaphash->{$scaf};
     # print "SCAF = $scaf\n";
      open (OUT, ">>$base_name/alignoutput/$scaf.txt");
      my $num = 0;
      foreach my $gapread (keys %$list){
        $num++;

        print OUT (">$scaf.read\n$gapread\n" x $unmaphash->{$scaf}{$gapread});
      }
      close OUT;
    }
  }
}

sub goThroughScaffold{
  my ($scaffile, $finalfile, $fillFile, $prevfillFile,$prevCloseFile) = @_;
  my ($seq, $prevhead, $ct) = ("","",0);
  open(SCAF,"$scaffile") || die "Can't open scaffile -- fatal\n";
  open (CLOSED, ">$finalFile") || die "Can't write to file -- fatal\n";
  open (FILL, ">$fillFile") || die "Can't write to file $fillFile-- fatal\n";
  #if there was already a closure, update the previous positions using the intermediate files
  open (PREVCLOSE, "$prevCloseFile") || die "Can't open file $prevCloseFile-- fatal\n" if($prevCloseFile ne "");
  open (PREVFILL, "$prevfillFile") || die "Can't open file $prevfillFile-- fatal\n" if($prevfillFile ne "");

  my ($totalgap, $totalclosed, $scafct, $totalgapnuc)= (0,0,0,0);
  &printMessage("\n=>".getDate().": Filling gaps\n");
  while(<SCAF>){
    chomp;
    $seq.= uc($_) if(eof(SCAF));
    if (/\>(\S+)/ || eof(SCAF)){
      my $head = $1;
      if($seq ne ""){
        $scafct++;
        my @gap = split(/N+/, $seq);
        my @countgapper;
        while ($seq =~ /(N+)/g) {
          my $start = (pos($seq)-length($1));
          my $end = pos($seq);
          my $lengap = length($1);
          my $gapstring = "$start|$end|$lengap";
          push @countgapper, $gapstring;
        }
        my ($finalseq, $closed) = ("", 0);
        #if there is a scaffold with at least one gap
        if($#gap >0){
          print FILL ">$prevhead:\n" if($prevCloseFile eq "");
          $totalgap = $totalgap + $#gap;
          my ($totgapIncl, $gapnum) = (0,0);
          # go through each gap, determine the position of the gap along the scaffold and try to fill the gaps.
          for (my $x=0;$x < $#gap;$x++){
            my ($alreadyfilled,$prevext1, $prevext2) = (0,0,0);
            $gapnum++;
            my $scafgap = "scaf$scafct.$gapnum";
            CounterPrint("Filling: $prevhead = scaffold$scafct"."|gap$gapnum      ");
            my $gapinfo = $countgapper[$x];
            my ($start,$end,$lengap) = split(/\|/,$gapinfo);
            $lengap = ($lengap+($trim*2)) if($prevCloseFile eq "");
            $totalgapnuc+=$lengap;
            $totgapIncl+=$trim;
            $start = $start-$totgapIncl;
            $end = $start+$lengap;

            my $gap1 = $gap[$x];
            $gap1 = substr($gap1,0,length($gap1)-$trim);
            $gap1 = $finalseq if($x >0);
            my $gap2 = $gap[($x+1)];
            my $subgap1 = substr($gap1,-300);
            $gap2 = substr($gap2,$trim);
            $gap2 = substr($gap2,0,length($gap2)-$trim) if($gapnum < $#gap);
            #use the information of previous extensions, only if multiple iterations are set
            if($prevCloseFile ne ""){
              while(<PREVFILL>){
                if(/closed=no/){
                  my (undef,$prevgapn,$prevstart,$prevend,$prevgapsize) = split("\t",$_);
                  $lengap = $1 if(/gapsize=(.*?)\s+/);
                  $alreadyfilled = $1 if(/filled=(.*?)\s+/);
                  $prevext1 = $1 if(/ext1=(.*?)\s+/);
                  $prevext2 = $1 if(/ext2=(.*?)\s+/);
                  print FILL "\t$prevgapn\t$prevstart\t$prevend\t$prevgapsize\t";
                  last;
                }else{
                  print FILL $_;
                }
              }
              while(<PREVCLOSE>){
                if(/=>/){
                  print CLOSEFILE $_ if(/extending with/ || /^>/);
                }else{
                  if(/sequences are NOT merged/){
                    print CLOSEFILE "\tIteration $iteration:\n";
                    last;
                  }  
                  if(/sequences are merged/){
                    print CLOSEFILE $_."\n";
                  }elsif(!/^\s*$/){
                    print CLOSEFILE $_;
                  }
                }
              }
            }else{
              print CLOSEFILE ">$prevhead"."|gap$gapnum\n";
              print FILL "\tgapnum=$gapnum\t";
              print FILL "start=$start\t";
              print FILL "end=$end\t";
              print FILL "gapsize=$lengap\t";
            }

            my ($gapseq, $closed) = fillGap($scafgap,$ct, $subgap1, $gap2, $lengap, $alreadyfilled,$prevext1, $prevext2);
            my $subgap3 = substr($gapseq,300);
            print FILL "closed=$closed\n";

            $finalseq = $gap1.$subgap3 if($x == 0);
            $finalseq = $finalseq.$subgap3 if($x >0);
            $closed++ if($closed ne "no");
            $totalclosed++ if($closed ne "no");

            $ct++;
          }
          print CLOSED ">$prevhead\n".wrap('', '', $finalseq)."\n";
        }else{
          print CLOSED ">$prevhead\n". wrap('', '', $seq)."\n";
        }
      }
      $prevhead = $head;
      $seq='';
    }else{
      $seq.= uc($_);
    }
  }
  CounterPrint("                                                  ");
  print SUMFILE "GAPCLOSING STATS:\n$seplines";
  print SUMFILE "\n\tClosed $totalclosed out of $totalgap gaps \n";
  print "\n\tClosed $totalclosed out of $totalgap gaps \n";
  close SCAF;
  close EVI;
  close FILL;
  close CLOSED;
  return $totalgapnuc;
}

sub fillGap{
  my ($gapnum, $ct1, $ctg1, $ctg2, $gapsize, $prevfilled, $prevext1, $prevext2) = @_;
  my $gapclosed = $prevfilled;
  my $readfile = "$base_name/alignoutput/$gapnum.txt";
  my $seqstig;
  my $numreads = 0;
  #collect the reads of a gap
  open(READS,"$readfile");
  while(<READS>){
    my $read = <READS>;
    $numreads++;
    chomp $read;
    my $subct = 0;
    while($subct < length($read)-$min_overlap){
      my $subseq = substr($read, $subct, $min_overlap+1);
      if(index($subseq,"N") == -1){
        if(defined $seqstig->{$subseq}){
          $seqstig->{$subseq}++;
        }else{
          $subseq = reverseComplement($subseq);
          $seqstig->{$subseq}++;
        }
      }
      $subct++;
    }
  }
  print FILL "reads=$numreads\t";
  my $seqstig_backup = $seqstig;
  my $prevlen1 = length($ctg1);
  my $prevlen2 = length($ctg2);
  #before extending, first try to merge the contigs
  my ($newseq, $exitstatus) = MergeTwoSequences($ctg1, $ctg2);
  my ($extend3output, $extend5output) = ("","");
  #if the gap can not be closed, either if the contigs are not merged, or if the difference between gapclosed and gapsize does not meet parameter -d
  if($exitstatus <=0 || !(($gapclosed) >= $gapsize-$difference && ($gapclosed) < $gapsize+$difference)){
    $exitstatus = 0;
    my ($extended3prime, $extended5prime, $iteration, $prevclosed) = (1,1,1,0);
    while(($extended3prime || $extended5prime) && $exitstatus <= 0){
      #extend one base at 3'
      ($seqstig, $newseq, $ctg1, $exitstatus,$gapclosed, $extended3prime, $extend3output) = gapFillExtension("3", $seqstig, $ctg1, $ctg2,$gapclosed,$gapsize);
      #if the gap is not closed already, extend one base at 5'
      if($exitstatus <=0){
        last if(($gapclosed-$min_tig_overlap)>($gapsize+$difference));
        $ctg2 = reverseComplement($ctg2);
        $ctg1 = reverseComplement($ctg1);
        ($seqstig, $newseq, $ctg2, $exitstatus,$gapclosed, $extended5prime, $extend5output) = gapFillExtension("5", $seqstig, $ctg2, $ctg1,$gapclosed,$gapsize);
        $ctg2 = reverseComplement($ctg2);
        $ctg1 = reverseComplement($ctg1);
        $newseq= reverseComplement($newseq);
      }
      last if(($gapclosed-$min_tig_overlap)>($gapsize+$difference));
    }
  }
  #try to merge the contigs with LocalAlignment (allowing mismatches/small indels) instead of perfect match alignment
  if($exitstatus == 0){
    print CLOSEFILE "$extend3output" if($extend3output ne "");
    print CLOSEFILE "$extend5output" if($extend5output ne "");
    ($newseq, $exitstatus) = LocalAlignment($ctg1, $ctg2, $gapclosed, $gapsize) ;
  }
  if($exitstatus <=0 || !(($gapclosed) >= $gapsize-$difference)){
    $exitstatus = 0;
    $newseq = "";
  }
  #determine statistics, like number of bases extended 5' and 3'
  my $extsize1 =  (length($ctg1)-$prevlen1);
  my $extsize2 =  (length($ctg2)-$prevlen2);
  $extsize1 = 0 if($extsize1 <0);
  $extsize2 = 0 if($extsize2 <0);
  $extsize1 += $prevext1;
  $extsize1 += $prevext2;

  my $filled = ($extsize1+$extsize2)-$exitstatus;
  $filled = 0 if($filled<0);
  print FILL "filled=$filled\text1=$extsize1\text2=$extsize2\tmerged=$exitstatus\t";
  my $Ngap = 1;
  $Ngap = ($gapsize-$filled) if($filled < $gapsize);
  my $finalseq = $ctg1.("N" x $Ngap).$ctg2;
  if($newseq ne ""){
    $finalseq = $newseq;
    print CLOSEFILE "\tsequences are merged\n\n";
    print FILL "remaining=0\t";
  }else{
    print FILL "remaining=$Ngap\t";
    print CLOSEFILE "\tsequences are NOT merged\n\n";
  }
  my $closed = "no";
  $closed = "yes" if($exitstatus>0);
  return $finalseq, $closed;
}

sub gapFillExtension{

   my ($direction, $bin, $seq, $seq2, $totalclosed, $gapsize) = @_;
   my $extended = 0;
   my $subseq = substr($seq, -$min_overlap);
   my $stopoutput = "";

   my $revseq = reverseComplement($subseq);
   my $overhang;
   #get number of occurences for extension of A,C,G,T
   $overhang->{'A'} = $bin->{$subseq."A"}+$bin->{"T$revseq"};
   $overhang->{'C'} = $bin->{$subseq."C"}+$bin->{"G$revseq"};
   $overhang->{'G'} = $bin->{$subseq."G"}+$bin->{"C$revseq"};
   $overhang->{'T'} = $bin->{$subseq."T"}+$bin->{"A$revseq"};

   #obtain total coverage
   my $coverage = $overhang->{'A'}+$overhang->{'C'}+$overhang->{'G'}+$overhang->{'T'};
   my $info = "\tclosed=$totalclosed/$gapsize\tdir$direction: Total:$coverage\tA:$overhang->{'A'}\tT:$overhang->{'T'}\tG:$overhang->{'G'}\tC:$overhang->{'C'}";
   if ($coverage < $base_overlap){
     $stopoutput =  "$info => coverage too low\n";
     return $bin, "", $seq, 0, $totalclosed, 0,$stopoutput;
   }
   my ($ct_dna, $previous_bz) = (0, "");
   my $extend_nuc = "";
   #determine most likely extension
   BASE:
   foreach my $bz (sort {$overhang->{$b}<=>$overhang->{$a}} keys %$overhang){
      if($ct_dna == 1){## the two most abundant bases at that position
        my $bestcoverage = $overhang->{$previous_bz} + $overhang->{$bz};
        if($overhang->{$previous_bz} < $base_overlap){
          $stopoutput = "$info => coverage too low\n";
          return $bin, "", $seq, 0, $totalclosed, 0,$stopoutput;
        }
        elsif(($overhang->{$previous_bz} / $bestcoverage) >= $min_base_ratio){### a simple consensus btw top 2
          $extend_nuc = "$previous_bz";
          print CLOSEFILE "$info => extending with $previous_bz\n";
          last BASE;
        }
        else{
           $stopoutput = "$info => below ratio\n";
           return $bin, "", $seq, 0, $totalclosed, 0,$stopoutput;
         }
      }
      $previous_bz = $bz;
      $ct_dna++;
   }
   my $checkseq = $seq . $extend_nuc;
   $totalclosed++;
   #check if two sequences can merge, if so, check if they correspond with the estimated gap distance
   my ($newseq, $size) = MergeTwoSequences($checkseq, $seq2);
   if($size > 0){
     my $gaprangemin = $gapsize-$difference;
     my $gaprangemax = $gapsize+$difference;
     if(($totalclosed) >= $gaprangemin){
       $extended = 0;
       return $bin, $newseq, $checkseq, $size,$totalclosed, 0,$stopoutput;
     }else{
       print CLOSEFILE "\t\tsequences can be merged, but difference between gap and total closed\n";
     }
   }
   deleteData($bin, $subseq."$extend_nuc");
   $seq = $checkseq;
   $extended = 1;
   if(($totalclosed-$min_tig_overlap)>($gapsize+$difference)){
     $extended = 0;
   }
   return $bin, "", $seq, 0, $totalclosed, $extended,$stopoutput;
}


#try to merge two sequences by at least -n basepairs using local alignment to remove single nucleotide differences at adges
sub LocalAlignment{
   my ($ctg1,$ctg2,$extended, $gapsize) = @_;
   my $max_overlap = ($min_tig_overlap+$min_overlap+$extended);
   my $seq1 = $ctg1;
   my $seq2 = substr($ctg2,0,$max_overlap);
   my $MATCH    =  1; # +1 for letters that match
   my $MISMATCH = -1; # -1 for letters that mismatch
   my $GAP      = -5; # -5 for any gap

  my $res = 0;
  my $resend = length($seq1);
  # initialization
  my @matrix;
  my ($i,$j) = (0,0);
  $matrix[0][$res]{score}   = 0;
  $matrix[0][$res]{pointer} = "none";
  while(++$j <= $resend){
      $matrix[0][$j]{score}   = 0;
      $matrix[0][$j]{pointer} = "none";
  }
  while(++$i <= length($seq2)){
      $matrix[$i][$res]{score}   = 0;
      $matrix[$i][$res]{pointer} = "none";
  }
  # fill
  my $align = "";
  my $max_i     = 0;
  my $max_j     = 0;
  my $max_score = 0;
  my $i = 0;
  while(++$i <= length($seq2)){
   my $j = $res;
   while(++$j <= $resend){
          my ($diagonal_score, $left_score, $up_score);

          # calculate match score
          my $letter1 = substr($seq1, $j-1, 1);
          my $letter2 = substr($seq2, $i-1, 1);
          if ($letter1 eq $letter2) {
              $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
          }
          else {
              $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
          }

          # calculate gap scores
          $up_score   = $matrix[$i-1][$j]{score} + $GAP;
          $left_score = $matrix[$i][$j-1]{score} + $GAP;

          if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
              $matrix[$i][$j]{score}   = 0;
              $matrix[$i][$j]{pointer} = "none";
              next; # terminate this iteration of the loop
          }

          # choose best score
          if ($diagonal_score >= $up_score) {
              if ($diagonal_score >= $left_score) {
                  $matrix[$i][$j]{score}   = $diagonal_score;
                  $matrix[$i][$j]{pointer} = "diagonal";
              }
              else {
                  $matrix[$i][$j]{score}   = $left_score;
                  $matrix[$i][$j]{pointer} = "left";
              }
          } else {
              if ($up_score >= $left_score) {
                  $matrix[$i][$j]{score}   = $up_score;
                  $matrix[$i][$j]{pointer} = "up";
              }
              else {
                  $matrix[$i][$j]{score}   = $left_score;
                  $matrix[$i][$j]{pointer} = "left";
              }
          }

          # set maximum score
          if ($matrix[$i][$j]{score} > $max_score) {
              $max_i     = $i;
              $max_j     = $j;
              $max_score = $matrix[$i][$j]{score};
          }
      }
  }

  # trace-back
  my $align1_1 = "";
  my $align2_1 = "";
  my $j = $max_j;
  my $i = $max_i;
  my ($start1,$end, $variantsfound,$start2) = (0,0,0,0);
  while (1) {
      if ($matrix[$i][$j]{pointer} eq "none" || $i == 0 || $j == 0){
        $start1 = ($j+1);
        $start2 = ($i+1);
        last;
      }

      if ($matrix[$i][$j]{pointer} eq "diagonal") {
          $end = ($i+1) if($end == 0);
          my $sub1 = substr($seq1, $j-1, 1);
          my $sub2 = substr($seq2, $i-1, 1);
            $align1_1 .= $sub1;
            $align2_1 .= $sub2;
          if($sub1 ne $sub2){
            my $diff1 = length($seq1)-$j;
            my $diff2 = $i;
            if($diff1 > $diff2){
              $align .= $sub1;
            }else{
              $align .= $sub2;
            }
          }else{
            $align .= $sub1;
          }
          $i--; $j--;
      }
      elsif ($matrix[$i][$j]{pointer} eq "left") {
          my $sub1 = substr($seq1, $j-1, 1);
            $align1_1 .= $sub1;
          $j--;
      }
      elsif ($matrix[$i][$j]{pointer} eq "up") {
          my $sub2 = substr($seq2, $i-1, 1);
            $align2_1 .= $sub2;
          $i--;
      }
  }
  $align = reverse($align);
  $align1_1 = reverse($align1_1);
  $align2_1 = reverse($align2_1);
  my $newseq = "";
  #if the alignment length is above parameter -n and the alignment is at the beginning of one of the contigs, the alignment is valid
  if(length($align) > $min_tig_overlap && $max_j >= length($seq1)-($trim+$difference) && $start2 <= ($trim+$difference)){
    my $index1 = rindex($ctg1, $align1_1);
    my $index2 = index($ctg2, $align2_1);

    my $sub1 = substr($ctg1,0,$index1);
    my $sub2 = substr($ctg2,$index2+length($align));
    my $trim1 = length($ctg1) - ($index1 + length($align1_1));
    my $trim2 = $index2;
    if($trim1 == 0 || $trim2 == 0){
      $newseq = $sub1.lc($align).$sub2;
      return ($newseq, length($align));
    }else{
      return ("", 0);
    }
  }
  return ("", 0);
}

#do a simple search for overlapping between the contigs, based on -n parameter
sub MergeTwoSequences{
  my ($ctg1, $ctg2) = @_;

  my ($max_overlap, $newseq) = ($min_tig_overlap+$min_overlap, "");
  while($max_overlap >= $min_tig_overlap){
    my $seq2 = $ctg2;
    my $seq1 = $ctg1;
    my $subseq2 = substr($ctg2,0,$max_overlap);
    my $subseq1 = substr($ctg1,-$max_overlap);
    if($subseq1 eq $subseq2){
      my $newctg1 = substr($ctg1,0,-$max_overlap);
      my $newctg2 = substr($ctg2,$max_overlap);
      $newseq = $newctg1.lc($subseq1).$newctg2;
      return ($newseq, length($subseq1));
    }
    $max_overlap--;
  }
  return ("", 0);
}

###DELETE READ DATA IF IT HAS BEEN USED FOR FILLING A GAP
sub deleteData {
   my ($bin, $sequence) = @_;
   my $comp_seq = reverseComplement($sequence);
   delete $bin->{$sequence};
   delete $bin->{$comp_seq};
   return $bin;
}

###WRITE SUMMARY STATISTICS FOR ALL CONTIGS OR SCAFFOLDS
sub writesummaryfiles{
  my ($input_file, $sumfile) = @_;

  open (INFILE, $input_file) || die "Can't open input file $input_file.\n";

  my ($seq, $name,$counter,$sum,$totalNcount,$totallen, $totalGC, $totalgap) = ("","",0,0,0,0,0,0);
  my (@line, @lengths);
  while (<INFILE>) {
    chomp;
    $seq.=$_ if(eof(INFILE));
    if ($_ =~ /^[>]/ || eof(INFILE)) { 
      if($seq ne ""){
        $counter++;
         push(@lengths, length($seq));
         my $len = length($seq);
         my $Gcount = () = $seq =~ /G/g;
         $totalGC = $totalGC + $Gcount;
         my $Ccount = () = $seq =~ /C/g;
         $totalGC = $totalGC + $Ccount;
         $sum+= $len;
         my @gap = split(/N+/, $seq);
         $totalgap = $totalgap + $#gap;
         my $Ncount = () = $seq =~ /[Nn]/g;
         $totalNcount += $Ncount;
         $name = "";
         $seq = "";
      }
  
      $name = $_;
    }
    else {
       $seq .= $_;
    }               
  }
  
  my $half_length = $sum/2;
  my $N25 = $half_length/2;
  my $N75 = $half_length/2+$half_length;
  
  my @lengths2 = reverse sort { $a <=> $b } @lengths;

  my ($sumN50, $foundN50, $foundN25, $foundN75) = (0,0,0,0);
  for(my $i = 0; $i <= $#lengths; $i++)
  {
    $sumN50 += @lengths2[$i];
    $foundN50 = @lengths2[$i] if($sumN50 >= $half_length && $foundN50 == 0);
    $foundN25 = @lengths2[$i] if($sumN50 >= $N25 && $foundN25 == 0);
    $foundN75 = @lengths2[$i] if($sumN50 >= $N75 && $foundN75 == 0);
  }
  my $GCcontent = sprintf("%.2f", (($totalGC/($sum-$totalNcount))*100));

  $sumfile .= "\t\tTotal number of scaffolds = $counter\n";
  $sumfile .= "\t\tSum (bp) = ". $sum. "\n";
  $sumfile .= "\t\t\tTotal number of N's = $totalNcount\n";
  $sumfile .= "\t\t\tSum (bp) no N's = ". ($sum-$totalNcount)."\n";
  $sumfile .= "\t\tGC Content = $GCcontent\%\n";
  $sumfile .= "\t\tMax scaffold size = ". @lengths2[0]."\n";
  $sumfile .= "\t\tMin scaffold size = ". @lengths2[$#lengths]."\n";
  $sumfile .= "\t\tAverage scaffold size = ".int($sum/$counter)."\n";
  $sumfile .= "\t\tN25 = ". $foundN25. "\n";
  $sumfile .= "\t\tN50 = ". $foundN50. "\n";
  $sumfile .= "\t\tN75 = ". $foundN75. "\n\n";
  
  close (INFILE);
  close OUTFILE;

  return $sumfile, $totalNcount;
}

###FUNCTION TO GET THE CURRENT DATE
sub getDate{
  my $date = scalar(localtime);
  return $date;
}

###FUNCTION TO PRINT MESSAGES TO THE SCREEN AND TO THE LOG FILE
sub printMessage{
  my $message = shift;
  print $message;
}

###FLUSHES THE SUMMARY AND LOG FILE
sub FlushFiles{
  select((select(SUMFILE), $| = 1)[0]);
  select((select(LOG), $| = 1)[0]);
  $|++;
}
###FUNCTION TO REVERSE COMPLEMENT A SEQUENCE
sub reverseComplement{
   $_ = shift;
   tr/ATGCatgc/TACGtacg/;
   return (reverse());
}

sub CounterPrint{
  my $countingMessager = shift;
  print "\r$countingMessager";
  $|++;
}
