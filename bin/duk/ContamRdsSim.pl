#!/jgi/tools/bin/perl -w
use strict;
use Getopt::Long;
# get input parameters
# input parameters
#
my $manual = <<ENDMANUAL;
 $0 [options] ref.fa contam.fa sim.fq 
 
   OPTIONS:
     --contam      contamimation rate (defualt 0.5)
     --error       error rate (default 0.005)
     --indel       indel rate (default 0.1)
     --intExt      insertion extension rate (default 0.2)
     --N           number of simulation reads (default 10000)
     --readLen     read length of the simulated reads (default 76)
     --contamLen   the mininum contaminated bases in a contaminated read (default 20)
     --help        help message

ENDMANUAL


my ($contamRate, $errorRate, $indelRate, $insertExtRate); 
my ($numRds, $readLenth, $mincontamLength); 
my ($refFile, $contamFile, $simFile);
my $help = 0;
# 
GetOptions("contam=f"=>\$contamRate,
            "error=f"=>\$errorRate,
            "indel=f"=>\$indelRate,
            "intExt=f"=>\$insertExtRate,
            "N=i",\$numRds,
            "readLen=i"=> \$readLenth,
            "contamLen=i"=> \$mincontamLength,
            "help" => \$help);
 
if ($help)  {
	print $manual, "\n";
	exit(0);
}           

$contamRate = 0.5 if (! defined($contamRate));
$errorRate = 0.005 if (! defined($errorRate));
$indelRate = 0.1 if (!defined($indelRate));
$numRds = 10000 if (!defined ($numRds));
$readLenth = 76 if (!defined($readLenth));
$mincontamLength = 20 if (!defined($mincontamLength));
$insertExtRate = 0.2 if (!defined($insertExtRate));

if ( $#ARGV < 2) {
	print $manual, "\n";
	exit(0);
}else {
  ($refFile, $contamFile, $simFile) = @ARGV;
}

print "contam   = ", $contamRate, "\n",
      "error    = ", $errorRate, "\n",
      "indel    = ", $indelRate, "\n",
      "intExt   = ", $insertExtRate, "\n",
      "N        = ", $numRds, "\n",
      "ReadLen  = ", $readLenth, "\n",
      "contamLen= ", $mincontamLength, "\n",
      "Referene file: ", $refFile, "\n",
      "Contam Tag file: ", $contamFile, "\n",
      "Simulated reads file: ", $simFile, "\n";
      
#print join "\t", $contamRate, $errorRate, $indelRate, $numRds, $readLenth, $mincontamLength, "\n";
#print join "\t", $refFile, $contamFile, $simFile, "\n";
my $refReads = "";
open REF, $refFile or die "Can not open $refFile \n";

while(<REF>) {
	chomp;
	last if (/^>/);
}

while(<REF>) {
	chomp;
	$refReads .= $_;
}

my $refLen = length($refReads);
close(REF);

my (@ContamRdNames, @ContamRdSeqs);
my $rdname = "";
my $rdSeq = ""; 
open CONTAM, $contamFile or die "Can not open $contamFile \n";
while( <CONTAM>) {
	chomp;
	if (/^>/) {
		if ($rdname) {
			push @ContamRdNames, substr($rdname,1);
			push @ContamRdSeqs, $rdSeq;		
		}
		$rdname = $_;
		$rdSeq = "";
	} else {
		$rdSeq .= $_;
	}
}
if ($rdname) {
	push @ContamRdNames, substr($rdname,1);
	push @ContamRdSeqs, $rdSeq;
}

close(CONTAM);
my $numContamRds = $#ContamRdNames + 1;

# print out
#for( my $i=0; $i < $numContamRds; $i++) {
#	print ">", $ContamRdNames[$i], "\n";
#	print $ContamRdSeqs[$i], "\n";
#}
my $cloneLength = 3*$readLenth;
my $refSmplLen = $refLen - $cloneLength;
my ($tagSub, $tagIndel, $genSub, $genIndel);
my $clone;
my @cloneRds;
my $tagLen;
my ($i, $j, $k, $ip);
my $simRdName;
my @simRdSeq;
my $startPos;
my %subBases = ( 'A' => ['G', 'C', 'T'],
	             'G' => ['A', 'C', 'T'],
	             'C' => ['A', 'G', 'T'],
	             'T' => ['A', 'G', 'C']);
my @bases = ('A','G', 'C', 'T' );
my $qval = chr(int(-10*log($errorRate)/log(10.0)) + 64);
my $simRdqual = "";
for ($i=0; $i <$readLenth; $i++) {
	$simRdqual .= $qval;
}
# open simulation read file for write
open SIMLF, ">$simFile" or die "Can not open $simFile for write \n";
# simulate the reads 
for ($i=0; $i < $numRds; $i ++ ) {
	#sample clone
	$startPos = int(rand($refSmplLen));
	$clone = substr($refReads, $startPos, $cloneLength);
	$simRdName = "P" . $startPos;
	# sample tags.
	$tagLen =0;
	if (rand() < $contamRate) {
		$j = int(rand($numContamRds)); # which tag
		$ip = length($ContamRdSeqs[$j]);
		$tagLen = int(rand($ip - $mincontamLength +1)) + $mincontamLength;
	    $clone = substr($ContamRdSeqs[$j], $ip - $tagLen, $tagLen) . $clone;
	    
	    $simRdName .= "&". $ContamRdNames[$j] . "&" .$tagLen;
	} else {
		$simRdName .= "&" . "NoTag&0";
	}
	# 
	#print "\@$simRdName\n",
	#      $clone, "\n";
	# simulate sequencing clone
	@cloneRds = split //, $clone;
	$j = 0;
	$k = 0;
	($tagSub, $tagIndel, $genSub, $genIndel) = (0, 0,0,0);
	while ( $j < $readLenth ) {
	    if( rand() < $errorRate ) { # sequencing error
	    	if (rand() < $indelRate) { #indel error
	    		if(rand() < 0.5) {  # deletion 
	    		   if ($k < $tagLen) {
	    		     $tagIndel++;
	    		   } else {
	    		    	$genIndel++;
	    		   }
	    		   #  
	    		   $k++;
	    		}else { #insertion
	    			$simRdSeq[$j] = $bases[int(rand(4))];
	    		    if ($k < $tagLen) {
	    		       $tagIndel++;
	    		    } else {
	    		    	$genIndel++;
	    		    }
	    		    $j++;
	    		    # insertion extension
	    		    while($j < $readLenth && rand() < $insertExtRate ) {
	    		    	$simRdSeq[$j] = $bases[int(rand(4))];
	    		        if ($k < $tagLen) {
	    		           $tagIndel++;
	    		        } else {
	    		    	    $genIndel++;
	    		        }
	    		        $j++;
	    		    }	    			
	    			
	    		}
	    	} else { #substitution error
	    		$ip = int(rand(3));
	    		$simRdSeq[$j] = $subBases{$cloneRds[$k]}[$ip];
	    		if ($k < $tagLen) {
	    		   $tagSub++;
	    		} else {
	    			$genSub++;
	    		} 
	    		$j++;
	    		$k++;
	    	}
	    } else {
	       $simRdSeq[$j] = $cloneRds[$k];
	       $j++;
	       $k++;	
	    }
	
	}
	# add the error to simulated read name
	$simRdName .= "&". $tagSub . ":". $tagIndel . ":" . $genSub . ":" . $genIndel;
	# write out the sequence.  
	print SIMLF '@', $simRdName, "\n",
	      join("", @simRdSeq), "\n",
	     "+\n",
	      $simRdqual, "\n";
	     
}

close(SIMLF);
