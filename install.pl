#!/usr/bin/env perl

# script that checks/compiles software required by the assembly protocol

#Carlos P Cantalapiedra (1), Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)
#1) Estacion Experimental de Aula Dei-CSIC, Zaragoza, Spain
#2) Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
#3) Fundacion ARAID, Zaragoza, Spain

use strict;
use Cwd;
use File::Basename;
use FindBin '$Bin';
use lib "$Bin/bin";

my $BINPATH = $Bin.'/bin';

my %binsoft = (
  'bowtie' => [ $BINPATH.'/SSPACE-BASIC-2.0_linux-x86_64/bowtie/bowtie ', 'Usage' ]
);

my %Csoft = ( 
  'bwa' => [ $BINPATH.'/bwa-0.7.6a/bwa', 'Program' ] ,
  'seqtk' => [  $BINPATH.'/seqtk/seqtk', 'Usage' ] ,
  'velvet' => [ $BINPATH.'/velvet_1.2.08/velvetg', 'Usage' ] 
);
  
my %CPPsoft = ( 
  'duk' => [ $BINPATH.'/duk/duk', 'Reference' ] ,
  'musket' => [ $BINPATH.'/musket-1.0.6/musket', 'Usage' ] ,
  'split_pairs' => [ $BINPATH.'/split_pairs_v0.5/split_pairs_bf', 'usage' ]  
);

my %javasoft = (
  'fastqc' => [ $BINPATH.'/FastQC/fastqc -v', 'FastQC v0.10.1' ],
  'trimmomatic' => [ "java -jar $BINPATH/Trimmomatic-0.32/trimmomatic-0.32.jar", 'Usage' ]
);

my ($soft,$path,$output);
my $cwd = getcwd();

print "\n## checking software provided in binary form: \n";

foreach $soft (keys(%binsoft))
{    
  $path = dirname($binsoft{$soft}[0]);
  print "# $soft ($path)\n";
  
  $output = `$binsoft{$soft}[0] 2>&1`;
  if($output !~ /$binsoft{$soft}[1]/)
  {
    die "<<fail, please obtain binaries suited to your system and place them in $path\n"; 
  }
  else{ print ">> OK\n"; }
}

print "\n## checking pre-compiled C software: \n";

foreach $soft (keys(%Csoft))
{    
  $path = dirname($Csoft{$soft}[0]);
  print "# $soft ($path)\n";
  
  $output = `$Csoft{$soft}[0] 2>&1`;
  if($output !~ /$Csoft{$soft}[1]/)
  {
    print "\n# compiling $soft in $path (requires gcc compiler) ...\n";
    chdir($path);
    if($soft eq 'velvet'){ system("make CATEGORIES=8 MAXKMERLENGTH=201 LONGSEQUENCES=1 2>&1") }
    else{ system("make 2>&1") }
    chdir($cwd);
    
    $output = `$Csoft{$soft}[0] 2>&1`;
    if($output !~ /$Csoft{$soft}[1]/)
    {
      die "<<fail auto-compiling, please install gcc compiler and/or resolve missing dependencies\n"; 
    }  
  }
  else{ print ">> OK\n"; }
}

print "\n## checking pre-compiled C++ software: \n";

foreach $soft (keys(%CPPsoft))
{    
  $path = dirname($CPPsoft{$soft}[0]);
  print "# $soft ($path)\n";
  
  $output = `$CPPsoft{$soft}[0] 2>&1`;
  if($output !~ /$CPPsoft{$soft}[1]/)
  {
    print "\n# compiling $soft in $path (requires g++ compiler) ...\n";
    chdir($path);
    system("make 2>&1");
    chdir($cwd);
    
    $output = `$CPPsoft{$soft}[0] 2>&1`;
    if($output !~ /$CPPsoft{$soft}[1]/)
    {
      die "<<fail auto-compiling, please install g++ compiler and/or resolve missing dependencies\n"; 
    }  
  }
  else{ print ">> OK\n"; }
}

print "\n## checking JAVA software: \n";

foreach $soft (keys(%javasoft))
{    
  $path = dirname($javasoft{$soft}[0]);
  print "# $soft ($path)\n";
  
  $output = `$javasoft{$soft}[0] 2>&1`;
  if($output !~ /$javasoft{$soft}[1]/)
  {
    die "<<fail auto-compiling, please install java run-time enviroment (jre)\n";   
  }
  else{ print ">> OK\n"; }
}
