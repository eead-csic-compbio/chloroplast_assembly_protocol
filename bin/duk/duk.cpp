#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <getopt.h>
#include "sys/time.h"
#include <unistd.h>
#include <stdlib.h>
#include "seqfiles.h"
#include "kmercoder.h"
#include "kmerhash.h"
#include "kmerhash.cpp"

using namespace std;

//define parameters
typedef  struct {
   string dbFile;
   string queryFile;
   string noMatchRdsFile;  
   string matchRdsFile;  
   string logFile;

   int merSize;
   int stepSize;
   int cutoff;

} Parameters;

int parseCmdOpts(Parameters& param, int argc, char ** argv);
double cdf_poisson(int x, float mu);

int main (int argc, char * argv[]) {
  Parameters toolParams;
  ostream * logOut;


  struct timeval tv1,tv2;
  time_t startClk = 0;
  time_t endClk = 0;
  

  
  parseCmdOpts(toolParams,  argc,  argv);
  // get starting time
  gettimeofday(&tv1,NULL);
  startClk = clock();
  // output parameters used
  if (!toolParams.logFile.empty() ) {
     logOut = new ofstream(toolParams.logFile.c_str(), ios_base::app);
     if (logOut->fail()) {
        delete logOut;
        logOut = &cout;
     }
  }else 
     logOut = &cout;
     
  *logOut << "\n##INPUT PARAMETERS##" <<endl
       << "#Reference file: " << toolParams.dbFile << endl
       << "#Query file: " << toolParams.queryFile << endl
       << "#Not matched reads file: " << toolParams.noMatchRdsFile << endl
       << "#Matched reads file: " << toolParams.matchRdsFile << endl
       << "#Output file: " << toolParams.logFile << endl
       << "#Mer size: " << toolParams.merSize << endl
       << "#Step size: " << toolParams.stepSize << endl
       << "#Cut off: " << toolParams.cutoff << endl;

// read the database file
 string rdname;
 int rdnumber = 0;
 const char* rdseq, *pseq;
 const char* rdqual;
 int dbstepSize = 1;
 long nKmers, nUkmers;
 long nRefRds, nRefBases;

 //
 nKmers = 0;
 nUkmers = 0;
 nRefRds = 0;
 nRefBases = 0;
// check the input sequence file typs
if ( seqfiletype( toolParams.dbFile) != 1) {
   *logOut << "\n\n# The reference file " << toolParams.dbFile << " must be in fasta format\n\n"<< endl;
    return 0;
}
int qryFileType = seqfiletype(toolParams.queryFile);
if ( qryFileType!= 1 && qryFileType != 2) {
   *logOut << "\n\n# The query file " << toolParams.queryFile << " must be in either fasta or fastq format\n\n"<< endl;
    return 0;
}
// build the database
 CFasta refSeqs(  toolParams.dbFile );
 if( refSeqs.open() != 1) {
   *logOut << " Can not open the database file " << toolParams.dbFile << endl;
   return 0;
 }
 // build up hash table
 CKmerCoder kmercoder( toolParams.merSize);
  unsigned long kmer = 0;
  unsigned long revkmer = 0;

  CKmerhash<int> Contamkmers;
 while ( refSeqs.getaread( rdname, rdnumber,  rdseq)== 1 ) {
   // update number of reads and bases;
   nRefRds++;
   nRefBases += strlen(rdseq);
   // get kmers
   kmer = 0;
   revkmer = 0;
   //first mer
   pseq = rdseq;
   if (kmercoder.getakmer( pseq,  kmer, toolParams.merSize,  toolParams.merSize) != 1)
      continue;
    //update hash and computer reverse kmer 
    nKmers += 2;
    Contamkmers[kmer]++ ;
    revkmer = kmer;
    kmercoder.reversekmer(revkmer);
    Contamkmers[revkmer]++;
   // get subsequent kmers 
   while( kmercoder.getakmer( pseq,  kmer,  dbstepSize,   toolParams.merSize) == 1) {
      nKmers += 2;
      Contamkmers[kmer]++ ;
       revkmer = kmer;
       kmercoder.reversekmer(revkmer);
        Contamkmers[revkmer]++;
   }
 }
 refSeqs.close(); 
 // output reference stat.
   nUkmers = Contamkmers.numKeys();
   *logOut << "\n##REFERECE STAT##" << endl
           << "#Total Reads: " << nRefRds << endl
	   << "#Total Bases: " << nRefBases << endl
	   << "#Total " << toolParams.merSize << "mers (FW+RV): " << nKmers << endl
	   << "#Total Unique " << toolParams.merSize << "mers(FW+RV): " << nUkmers << endl;

//Contamkmers.analysis();

  //read the query fastq file, classify whether each read is contamiated or not.
  CSeqfile * qrySeqs;
  if (qryFileType == 1 )
    qrySeqs = new CFasta(toolParams.queryFile);
  else
    qrySeqs = new CFastq(toolParams.queryFile);
  if(qrySeqs->open() != 1) {
     *logOut << " Can not open the query file " <<  toolParams.queryFile << endl;
     delete qrySeqs;
     return 0;
  }

  int bOutputGoodRds, bOutputBadRds;
  unsigned long totalRds, nBadRds;
  int badOcc =0;
  char strOcc[12];

 CSeqfile *goodSeqs,* badSeqs;
 if (!toolParams.noMatchRdsFile.empty() ) {
     bOutputGoodRds = 1;
     if (toolParams.noMatchRdsFile == "-") 
        if (qryFileType == 1)
           goodSeqs = new CFasta( "", 1);
	 else
           goodSeqs = new CFastq( "", 1);
     else 
        if (qryFileType == 1)
          goodSeqs = new CFasta( toolParams.noMatchRdsFile, 1);
        else
          goodSeqs = new CFastq( toolParams.noMatchRdsFile, 1);

     goodSeqs->open();
  } else
     bOutputGoodRds = 0;
  
  if (!toolParams.matchRdsFile.empty() ) {
     bOutputBadRds = 1;
     if (qryFileType == 1)
        badSeqs = new CFasta( toolParams.matchRdsFile, 1);
     else
        badSeqs = new CFastq( toolParams.matchRdsFile, 1);
     badSeqs->open();
  }else
     bOutputBadRds = 0;
// bad kmer occurrence histogram
  long histKmer[256];
  int histBin = 255;

  for (int i=0; i<=histBin; i++) 
     histKmer[i] = 0;
  // initialize the total number of reads and bad reads;
  totalRds =0;
  nBadRds = 0;
  //
  long nTotalKmers = 0;
// read through query file
 while ( qrySeqs->getaread( rdname, rdnumber,  rdseq, rdqual)== 1 ) {
  //  cout << ">" <<  rdname << "  " << rdnumber << endl << rdseq << endl;
   totalRds++;

   kmer = 0;
   badOcc = 0;
   //first mer
   pseq = rdseq;
   if (kmercoder.getakmer( pseq,  kmer, toolParams.merSize,  toolParams.merSize) != 1) {
     /*
      nBadRds++;
     if (bOutputBadRds) {
         sprintf( strOcc, "\t%d\t", badOcc);
	 rdname = rdname +   strOcc;
        badSeqs.writearead(rdname , rdnumber,  rdseq, rdqual); 
     }
      continue; 
     */
   } else {
     nTotalKmers++;
     if (Contamkmers.find(kmer))
       badOcc++;
   }
   while( kmercoder.getakmer( pseq,  kmer, toolParams.stepSize,   toolParams.merSize) == 1) {
     nTotalKmers++;
     if (Contamkmers.find(kmer))
       badOcc++;
   }
    //
   if( badOcc >= toolParams.cutoff ) {
      nBadRds++;
     if (bOutputBadRds) {
        // sprintf( strOcc, "\t%d\t", badOcc);
	 //rdname = rdname +   strOcc;
        badSeqs->writearead(rdname , rdnumber,  rdseq, rdqual); 
     }
     // update the histogram
     if (badOcc <  histBin)
        histKmer[badOcc]++;
     else
       histKmer[histBin]++;  

   }else {
     if (bOutputGoodRds)
        goodSeqs->writearead(rdname, rdnumber,  rdseq, rdqual); 
   }
   
 }
  qrySeqs->close();
  delete qrySeqs;
  if( bOutputBadRds){ 
     badSeqs->close();
     delete badSeqs;
  }
  if( bOutputGoodRds) {
     goodSeqs->close();
     delete goodSeqs;
  }
  // compute expected value;
  float avgKmers = nTotalKmers*1.0/totalRds ;
  unsigned long nTKs = 1;
  nTKs = nTKs << (2* toolParams.merSize);
  float mu = avgKmers*nUkmers/nTKs;
  float pValue = 1- cdf_poisson(toolParams.cutoff , mu);
  // get computer time
  gettimeofday(&tv2,NULL);
  endClk = clock();
//output the running time
 *logOut << "\n## ELAPSED TIME##\n"
      << "# Readl time: " << tv2.tv_sec - tv1.tv_sec << " seconds, " 
      << "# CPU+SYS time: " << float(endClk-startClk)/CLOCKS_PER_SEC << "seconds" << endl;
 //output the simple stat.
  *logOut << "\n##QUERY FILE STAT##\n"
       << "# Total number of reads:    " << totalRds << endl
       << "# Total number of matched reads: " << nBadRds << endl
       << "# Match ratio:       " << (float)nBadRds/totalRds << endl;
  //output the p-value
  *logOut <<"\n##P-VALUE##\n"
          << "#Avg number of Kmer for each read:  " << avgKmers << endl
          <<"# P value for the given threshold " << toolParams.cutoff << " is " << pValue << endl;
//output the histogram
  int lastBin=histBin;
  while(histKmer[lastBin] == 0 && lastBin > -1)
      lastBin--;
  *logOut << "\n## Histogram of kmer occurance for reads with at least one occurance ## \n"
       << "#NumOcc\tNumReads\tPercentage\n";
  for (int i=0; i<= lastBin; i++)
     *logOut << i << "\t" << histKmer[i] << "\t" << setprecision(4) << histKmer[i] * 1.0/nBadRds << endl;


if (logOut != &cout) {
    delete logOut;
 }
 return 0;


}
