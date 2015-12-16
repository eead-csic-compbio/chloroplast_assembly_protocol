#include  <iostream>
#include <fstream>
#include <string>
#include "seqfiles.h"
using namespace std;
// CSeqfile class implementation
int CSeqfile::openfile() {
    // if seqFilename is empty, use std input or output according to openMode.
    if (seqFilename.empty() ) {
        if (openMode)
	    seqStrm = (fstream *) & cout;
        else
	    seqStrm = (fstream *) & cin;
        return 1;
    }
    // open the file
    if (openMode)
        seqStrm = new  fstream(seqFilename.c_str(), fstream::out);
    else
        seqStrm = new  fstream(seqFilename.c_str(), fstream::in);
    // check whether the file is open
    if (seqStrm->is_open()) 
       return 1;
    else
       return -1;
}
// close the file stream and release space
void CSeqfile::closefile() {
  // if seqStrm is not cin or cout or NULL, delete it.
  if ( !seqFilename.empty() && seqStrm ) {
       delete seqStrm;
       seqStrm = NULL;
   }
}
// parseReadName
// Given a read header, parse it and save read name and number into variables readname and readnumber;
// If success, return 1
//   otherwise, return 0
int CSeqfile::parseReadName(string& rdHeader) {

  // parse the read name
  if (rdHeader.empty())
     return 0;
  
  size_t ipos= rdHeader.find('/');
  if ( ipos == string::npos) {
     readnumber = 0;
     readname = rdHeader.substr(1);
  } else {
     readnumber = rdHeader[ipos+1] - '0';
     readname = rdHeader.substr(1, ipos-1);
 }

 return 1;
}


// CFasta class implementation
CFasta::~CFasta() {
  close();
}
// open the file stream
int CFasta::open() {
  int i = openfile();
  // return if the file stream is not open properly. 
  if ( i != 1 || openMode)
    return i;

  // get the read name and number.
  // if the first line is not namne line or end of file, return -1
  if (!getline(*seqStrm,  readname) || readname[0] != '>' ) 
     return 0;
  else 
     return  parseReadName( readname);
}
//
//
int CFasta::open(  const string& filename  ) {
   seqFilename = filename;
   return open();
}

int CFasta::open(  const string& filename, bool bOpenMode) {
    seqFilename = filename;
    openMode = bOpenMode;

    return open();
}

//
int CFasta::getaread(string& rdname, int& rdnumber, const char *& rdseq) { 
  string tempstr;
  //check whether the file is open
  if (!seqStrm || openMode) {
     return -1;
  }
  // asign the read name and number
  rdnumber = readnumber;
  rdname = readname;
 // get the read sequence 
  readseq.clear();
  while(getline(*seqStrm,  tempstr)) {
   if (tempstr[0] == '>')
       break;
   readseq.append(tempstr);
 }
// check whether sequence is empty
  if (readseq.empty()) {
     rdseq = NULL;
     return 0;
  }
  rdseq = readseq.c_str();
  // parse the next read name 
  parseReadName( tempstr);
 
  return 1;
}
//
 int CFasta::getaread(string& rdname, int& rdnumber, const char *& rdseq, const char*& rdqual) {

   return getaread(rdname, rdnumber, rdseq);
 }
//

int CFasta::writearead(string& rdname, int rdnumber, const char * rdseq){

  //check whether the file is open
  if (!seqStrm || !openMode) {
     return -1;
  }
 // write out
  if (rdnumber != 0)
   ( *seqStrm )<< '>' << rdname << '/' <<  rdnumber << endl
           << rdseq << endl;
   else
    (*seqStrm) << '>' << rdname << endl
            << rdseq << endl;
  
  return 1;
}
//
int CFasta::writearead(string& rdname, int rdnumber, const char * rdseq, const char* rdqual) {
  return writearead( rdname, rdnumber,  rdseq);
}
// implementation for CFastq  class
//deconstructor
CFastq::~CFastq()  {
  close();
}
// open the fastq file stream
int CFastq::open() {
    
    return openfile();
}
int CFastq::open( const string& filename, bool newOpenMode ) {
   // reset the file name and openmode 
    seqFilename = filename;
    openMode = newOpenMode;

    return openfile();
}
int CFastq::open( const string& filename ) {
   // reset the file name  
    seqFilename = filename;

    return openfile();
}

// get a read

int CFastq::getaread(string& rdname, int& rdnumber, const char *& rdseq, const char*& rdqual ) {
 int i = getaread(rdname, rdnumber, rdseq);

 if (i != 1) 
   return i;
 else {
   rdqual = readqual.c_str();
   return 1;
 }
   

}
// get a read 
int CFastq::getaread(string& rdname, int& rdnumber, const char *& rdseq) {

 string  tempstr;
   //check whether the file is open
 if (!seqStrm || openMode) {
     return -1;
 }
 // get the read name
  while(getline(*seqStrm,  readname)) {
     if (readname[0] == '@')
        break;
 }
 // parse the read name
  
  int is = parseReadName( readname);
  if ( is != 1)
     return is;
  // get read sequence
  if( ! getline(*seqStrm,  readseq)) {
       return 0;
  }
  // 
   //get '+' symbol
  if (! getline(*seqStrm,  tempstr))
     return 0; 
  
  //get quality sequence
  if (! getline(*seqStrm,  readqual))
     return 0;

  // assign the value to input papameters
  rdname = readname;
  rdnumber = readnumber;
  rdseq = readseq.c_str();

  return 1;
}
// write a read
int  CFastq::writearead(string& rdname, int rdnumber, const char * rdseq, const char* rdqual) {
   // check the file stream
 if (!seqStrm || !openMode) {
     return -1;
 }
 // write the read to the output file stream
   if (rdnumber != 0)
    ( *seqStrm )<< '@' << rdname  << '/' <<  rdnumber << endl 
	        << rdseq << endl << '+' << endl << rdqual << endl;
    else
     (*seqStrm) << '@' << rdname << endl
	        << rdseq << endl << '+' << endl << rdqual << endl;

   return 1;

}
//
// seqfiletype -- return the sequence file type, 
//    If the input file name string is empty, the standard input (cin) is used.
//    The subroutine peeks at the first character of the file, it is a fasta file if the file starts with '>', 
//        it is a fastq file if the file starts with '@'.
//  return value
//    -1 - file not exist
//    0  - not known
//    1  - fasta
//    2  - fastq
int seqfiletype (string & seqfilename) {
   fstream* fstrm;
   
   // if the input string is empty, use cin
   if (seqfilename.empty()) {
      fstrm = (fstream *) &cin;
   } else { 
      fstrm = new fstream(seqfilename.c_str(), fstream::in);
      // if file can not be opened, return -1. 
      if ( !fstrm->is_open()) {
          delete fstrm;
	  return -1;
      }
   }
   // peek at the first character.
   char c = fstrm->peek();
   if (!seqfilename.empty())
      delete fstrm;
   if ( c== '>' )
     return 1;
   else if (c == '@')
          return 2;
        else 
	  return 0;

}
