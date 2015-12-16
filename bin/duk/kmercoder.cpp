#include <iostream>
#include "kmercoder.h"
using namespace std;

// initialize the code.
CKmerCoder::CKmerCoder(int msize): mersize( msize) {
  
  int i;

  for (i=0; i<255; i++) {
      coder[i] = 0;
  }
// code
   coder[67] = 1; //C
   coder[99] = 1; //c
   coder[71] = 2; //G
   coder[103] = 2; //g
   coder[84] = 3; //T
   coder[116] = 3; //t
// decode
   decoder[0] = 'A';
   decoder[1] = 'C';
   decoder[2] = 'G';
   decoder[3] = 'T';
   // bitmask;
   bitmask = 0;
   for(i=0; i<mersize; i++) {
      bitmask <<=2;
      bitmask |= 3;
   }
   //
  //   printKmer(bitmask);
}
// 
int CKmerCoder::getakmer(const char * & seq, unsigned long& kmer, int stepsize, int window)
{ 
  // shift kmer
  while( (*seq) != '\0' && stepsize > 0) {
   if( *seq == 'N') {
      seq++;
      stepsize = window;
      continue;
    }
    kmer <<=  2;
    kmer |= coder[*seq];
      
    seq++;
    stepsize--;
  }
  kmer &= bitmask;

  if (stepsize > 0)
    return 0;
  else
    return 1;
}
// given a kmer, output its sequence
char * CKmerCoder::decodeKmer(char * seq, int mersize, unsigned long kmer) {

   for ( int i=mersize-1; i> -1; i--) {
      seq[i] = decoder[kmer & 3];
      kmer >>= 2;
   }

   return seq;

}
//
// get the reverse kmer of a given kmer.
int CKmerCoder::reversekmer(unsigned long& kmer) {
  kmer = ~kmer; 
  
  kmer = (kmer & 0x3333333333333333) << 2 | (kmer & 0xCCCCCCCCCCCCCCCC) >> 2;
  kmer = (kmer & 0x0F0F0F0F0F0F0F0F) << 4 | (kmer & 0xF0F0F0F0F0F0F0F0) >> 4;
  kmer = (kmer & 0x00FF00FF00FF00FF) << 8 | (kmer & 0xFF00FF00FF00FF00) >> 8;
  kmer = (kmer & 0x0000FFFF0000FFFF) << 16 | (kmer  & 0xFFFF0000FFFF0000) >> 16;
  kmer = (kmer & 0x00000000FFFFFFFF) << 32 | (kmer  & 0xFFFFFFFF00000000) >> 32;

  kmer >>= (64-2*mersize);

  kmer &= bitmask;
  return 1;

}

// print a kmer in binary format, two bits for each base.
void CKmerCoder::printKmer(unsigned long kmer) 
{
   int i;
   int lz = sizeof(long)*4;
   unsigned long mask = 3 ;
    mask = mask << (8*sizeof(long)-2);
   for ( i=lz-1; i>=0; i--) {
      cout << ( (kmer & mask) >> (i*2)  )<< " " ;
      mask >>=  2;
   }
   cout << endl;
}


