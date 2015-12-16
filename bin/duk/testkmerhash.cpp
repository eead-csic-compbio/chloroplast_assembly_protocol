#include <iostream>
#include "kmerhash.h"
#include "kmerhash.cpp"
using namespace std;

int main(int argc, char* argv[]) {
   
   CKmerhash<int> wordCnt;

   unsigned long key;

   int value;

   wordCnt.insert( 1, 4);
   wordCnt.insert(5, 2);
   wordCnt.insert(3, 3);
   wordCnt.insert(11, 10);

   
   wordCnt[3] = 13;
   CKmerhash<int>::iterator it1 = wordCnt.begin();
   CKmerhash<int>::iterator itEnd = wordCnt.end();
   /*
   while( it1 != itEnd) {
      cout << "key " <<  *it1 << "\t" << wordCnt[*it1] << endl;
     it1++;
   }
   */
   for ( CKmerhash<int>::iterator it = it1; it!= itEnd; it++) {
      cout << "key " <<  *it << "\t" << wordCnt[*it] << endl;
   }
}

