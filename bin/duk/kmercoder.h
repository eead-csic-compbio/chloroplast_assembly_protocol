#ifndef KMERCODER_H
#define KMERCODER_H



#include <iostream>
using namespace std;

class CKmerCoder{
    int coder[255];
    char decoder[4];
   int mersize;
   unsigned long bitmask;
 public:

     CKmerCoder(int msize);
     ~CKmerCoder(){ };

     int getakmer(const char * & seq, unsigned long& kmer, int stepsize, int window);
     int reversekmer(unsigned long& kmer);

     char * decodeKmer(char * seq, int mersize, unsigned long kmer);
     void  printKmer(unsigned long kmer);
};

#endif

