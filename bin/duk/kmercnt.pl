#!/jgi/tools/bin/perl -w

use strict;

my $qryFile = $ARGV[0];

my $rdname;
my $rdseq;
my %kmercnt;
my $kmer;

my $mersize = 5;
my $step = 1;
my $offside;
my $rdlen;
open SEQ, $qryFile or die " Can not open $qryFile";
my $lc = 0;
while(<SEQ>) {
  chomp;
  if(/@/) {
   $lc = 1;
   next;
  }
  # seq line
  if ($lc == 1) {
    $rdseq = $_;
    $rdlen = length($rdseq);
    $offside = 0;
   while ( ($offside + $mersize)<= $rdlen) {
      $kmer = substr($rdseq, $offside, $mersize);
      if (exists $kmercnt{$kmer}) {
         $kmercnt{$kmer}++;
      } else {
         $kmercnt{$kmer} = 1;
      }
      $offside += $step;
   }

   $lc++;
  }

}
#output kmer count
print "Kmer\tCount\n";
foreach $kmer (keys %kmercnt) {
  if ($kmer !~ /N/) { 
    print "$kmer\t", $kmercnt{$kmer}, "\n";
  }
}
