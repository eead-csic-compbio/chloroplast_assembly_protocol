#include <iostream>
#include <string>
#include "seqfiles.h"
#include "kmercoder.h"
#include "kmerhash.h"
#include "kmerhash.cpp"

using namespace std;

int main(int argc, char* argv[]) {
  string qryFile;
  int mersize = 5;

  // get query file name
  if (argc <2) 
    qryFile = "";
  else 
    qryFile = argv[1];

  cerr<< "Input Fastq file " << qryFile << endl;

  // get input 
  CKmerCoder kmercoder( mersize);
  unsigned long kmer = 0;
  CKmerhash<int> kmerCounter;

  CFastq qrySeqs(  qryFile );
  if (qrySeqs.open() != 1) {
   cerr << "Can't open the file " << qryFile << endl;
   exit(0);
  }

  int stepsize = 1;
  string rdname;
  int rdnumber;
  const char* rdseq, *pseq;

  while ( qrySeqs.getaread( rdname, rdnumber,  rdseq)== 1 ) {
    kmer = 0;
    //first kmer
    pseq = rdseq;
    if (kmercoder.getakmer( pseq,  kmer, mersize,  mersize) != 1)
         continue;

    kmerCounter[kmer]++;
    // get subsequent kmers
    while( kmercoder.getakmer( pseq,  kmer,  stepsize,   mersize) == 1) {
      kmerCounter[kmer]++;
   }
 }

 qrySeqs.close();

 int  nUkmers = kmerCounter.numKeys();
 cerr <<"Number of unique kmers " << nUkmers << endl;

 // output kmers and count

   CKmerhash<int>::iterator it1 = kmerCounter.begin();
   CKmerhash<int>::iterator itEnd = kmerCounter.end();
 
 char kmerword[20];
 kmerword[mersize] = '\0';
 int cnt;
  cout << "Key\tCount\n";

   for ( CKmerhash<int>::iterator it = it1; it!= itEnd; it++) {
       kmer = *it;
       cnt = kmerCounter[kmer];
       kmercoder.decodeKmer(kmerword,  mersize,  kmer); 
       cout <<  kmerword << "\t" <<  cnt << endl;
   }


}
