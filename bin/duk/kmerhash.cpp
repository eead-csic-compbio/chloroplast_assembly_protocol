#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include "kmerhash.h"

using namespace std;

template<class ValueClass>
CKmerhash<ValueClass>::CKmerhash( int kmer): mersize(kmer), tablesize(0), numkeys(0) {
  //  mersize = 14;
  nbits = 2*mersize;
  tablesize = 1 << nbits;
  loadingLimit = (int)(tablesize * LOADINGFACTORTHRES);
 

  bitmask = 1;
  for (int i=1; i< nbits; i++)
     bitmask |=  1 << i;

   hashtable=NULL;
   hashValue=NULL;
  
   hashtable = new unsigned long[tablesize];
   hashValue = new ValueClass[tablesize];
  //
  for(int i=0; i<tablesize; i++)
     hashtable[i] = EMPTY;

}
//
template<class ValueClass>
CKmerhash<ValueClass>::~CKmerhash() {

  delete [] hashtable;
  delete[] hashValue;

}
// index -- find the index of key  
template<class ValueClass>
int CKmerhash<ValueClass>:: index(ulong key) {

int ind = ((key >> nbits) ^  key);
    ind += ~(ind << 9);
    ind = ind & bitmask;

while ( (hashtable[ind] != EMPTY) && (hashtable[ind] != key) ) {
    ind ++;
    // reach to the end of table,start search from beginiing
    if (ind >= tablesize) 
      ind = 0;
 }
 return ind;
}
//

template<class ValueClass>
void CKmerhash<ValueClass>:: doubletablesize() {
   int oldTablesize = tablesize;
   int ind;
   ulong *oldHashtable = hashtable;
   ValueClass *oldHashValue=hashValue;
    // doulbe tablesize
    nbits++ ;
    tablesize = tablesize << 1;
    bitmask = (bitmask << 1) | 1;
    loadingLimit =  loadingLimit << 1;
    // allocate the new memory
    hashtable = new ulong[tablesize];
    hashValue = new ValueClass[tablesize];
    //initialize the new table;
    for(int i=0; i<tablesize; i++)
        hashtable[i] = EMPTY;
    // loop throught the old table
    for(int i=0; i<oldTablesize; i++) {
      if (oldHashtable[i] != EMPTY) {
         ind = index(oldHashtable[i]);
	 hashtable[ind] = oldHashtable[i];
	 hashValue[ind] = oldHashValue[i];

//	 cout << i << "\t" << ind << "\t" << oldHashtable[i] << "\t" << hashtable[ind] << endl;
      }
    }
   // release the old table's memeor
   delete[] oldHashtable;
   delete[] oldHashValue;
}

//
template<class ValueClass>
int CKmerhash<ValueClass>:: insert(ulong key, const ValueClass& val) {
    // resize the hashing table if the table is loaded too much
    if ( numkeys >  loadingLimit ) {
      doubletablesize();
    }

    int ind = index(key);
    if ( hashtable[ind] != key) 
       numkeys++;

    hashtable[ind] = key; 
    hashValue[ind] = val;
    return ind;

}
//

template<class ValueClass>
int CKmerhash<ValueClass>:: insert(ulong key) {
    if ( numkeys >  loadingLimit ) {
         doubletablesize();
    }

    int ind = index(key);
    if ( hashtable[ind] != key)
         numkeys++;
 
   hashtable[ind] = key; 

   return ind;
}
//
template<class ValueClass>
int CKmerhash<ValueClass>:: find(ulong key, ValueClass& val) {
 int ind = index(key);

 if ( hashtable[ind] != EMPTY ) {
       val = hashValue[ind];
       return 1;
 } else
    return 0;

}
//

template<class ValueClass>
int CKmerhash<ValueClass>:: find(ulong key) {
 
 int ind = index(key);

 if ( hashtable[ind] != EMPTY ) {
       return 1;
 } else
    return 0;

}

//

template<class ValueClass>
ValueClass& CKmerhash<ValueClass>:: operator[] (ulong key) {

   if ( numkeys >  loadingLimit ) {
       doubletablesize();
   }

    int ind = index(key);
    if ( hashtable[ind] != key)
          numkeys++;

  hashtable[ind] = key;

  return hashValue[ind];
 
 } 

// analyze the hash performance
template<class ValueClass>
void  CKmerhash<ValueClass>::analysis () {
 int totalkey = 0;
 ulong totalProbe = 0;
 int maxProbe = 0;

 cout << "Hash table information " << endl;

 int i;
 long key, ind, diff;

 for(i=0; i< tablesize; i++) {
    if( hashtable[i] != EMPTY) {
      totalkey++;
      key = hashtable[i];

    ind = ((key >> nbits) ^  key);
    ind += ~(ind << 9);
    ind = ind & bitmask;
  
  if ( i >= ind ) 
        diff = i - ind ;
      else
        diff = i+ (tablesize -ind) + 1; 
      if (diff > maxProbe) {
         maxProbe = diff;
      }
      totalProbe += diff;

   //   cout << i << "  " << ind <<  "   " << key <<"\t"<<hashValue[i]<< "\t" << diff << endl;
    }

 }
 // output
 cout << "Loading Limit: " << loadingLimit << endl;
 cout << "total keys:  " << totalkey << endl
      << "total keys:  "<< numkeys << endl
      << "tablesize:  "<< tablesize << endl
      << "Max Probe:   " << maxProbe << endl
      << "total Probe:   " << totalProbe << endl
      << "Avg Probe:   " << (totalProbe * 1.0)/totalkey << endl;
}
