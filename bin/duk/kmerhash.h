#ifndef KMERHASH_H
#define KMERHASH_H

#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <time.h>

using namespace std;
//.
// constant
typedef unsigned long ulong;

const ulong EMPTY = 0xFFFFFFFFFFFFFFFF;
//const ulong EMPTY = -1;
const float LOADINGFACTORTHRES = 0.5;
template<class ValueClass>
class CKmerhash{
   int mersize;
   int tablesize;
   int numkeys;

   int nbits;
   ulong bitmask;
   int loadingLimit;

   ulong * hashtable;
   ValueClass * hashValue;

   int index( ulong key);
   void doubletablesize();
 public:
   CKmerhash(int kmer=10);
   ~CKmerhash();
   
   int insert(  ulong key, const ValueClass& val);
   int insert( ulong key);
   int find( ulong key, ValueClass& val );
   int find( ulong key );

   ValueClass& operator[](ulong key);
   //analyze the hash performance
   void analysis(); 
   int numKeys() { return numkeys; }

   // iterator class
   class iterator;
   friend class iterator;

   class iterator {
       CKmerhash& kmerhash;
       int index;
       void findkey() {
         while ( index <  kmerhash.tablesize && kmerhash.hashtable[index] == EMPTY)
	     index++;
       }
     public:
       iterator(  CKmerhash<ValueClass>& newHash) : kmerhash(newHash), index(0) { findkey(); }
       iterator( CKmerhash<ValueClass>& newHash, bool): kmerhash(newHash), index(newHash.tablesize)  {}
       iterator(const iterator& i1) : kmerhash(i1.kmerhash), index(i1.index) {}

       iterator& operator++ () {
        if ( index <  kmerhash.tablesize ) {
            index++;
	    findkey();
	 } 
	  return *this;
       }

       iterator& operator++ ( int) {
	  return operator++();
       }
       ulong operator * () const {
            return kmerhash.hashtable[index]; 
       }
       
       bool operator == (CKmerhash<ValueClass>::iterator& hash2) {
          return index == hash2.index;
       }

       bool operator != (CKmerhash<ValueClass>::iterator& hash2) {
          return index != hash2.index;
       }

   };

   iterator begin()  {
     return iterator(*this);
   }

   iterator end()  {
     return iterator(*this, true);
   }
};

#endif
