#ifndef SEQFILES_H
#define SEQFILES_H
//
#include  <iostream>
#include <fstream>
#include <string>
using namespace std;

class CSeqfile {
  protected:
    string seqFilename;
    fstream *seqStrm;
    bool openMode;

    string readname;
    int readnumber;
    string readseq;
    string readqual;

     int openfile();
     void closefile();
     int parseReadName(string& rdHeader);
  public:
     CSeqfile ( ):seqStrm(NULL), openMode(0), readnumber(0) {}
     CSeqfile ( bool fmode):seqStrm(NULL), openMode(fmode), readnumber(0) {}
     CSeqfile (const string& seqFile, bool fmode=0 ):seqFilename(seqFile), seqStrm(NULL), openMode(fmode), readnumber(0) {}
     virtual ~CSeqfile() {} 
    //set data member
    //  void setFilename(const string & fname) {seqFilename=fname; }
     // void setOpenmode(bool newMode) { openMode = newMode;}
     void printVars() { cout << seqFilename << endl << openMode << endl;}
    // operation
     //open and close the file stream    
     virtual int open(  const string& filename, bool newOpenMode  ) =0;
     virtual int open(  const string& filename  ) =0;
     virtual int open()=0; 
     virtual void close() {closefile();}
     // read and write operation
     virtual int getaread(string& rdname, int& rdnumber, const char *& rdseq) =0; 
     virtual int getaread(string& rdname, int& rdnumber, const char *& rdseq, const char*& rdqual) { return -1;} 

     virtual int writearead(string& rdname, int rdnumber, const char * rdseq, const char* rdqual)
       {return 1;} 
     virtual int writearead(string& rdname, int rdnumber, const char * rdseq)
       {return 1;} 

};
// fastq file
class CFastq : public CSeqfile {
  public:
     CFastq(bool fmode) :CSeqfile(fmode) { }
     CFastq():CSeqfile(0) {}
     CFastq( const string& seqFile, bool fmode=0):CSeqfile(seqFile, fmode) {}

     ~CFastq();
    // operation
    // open and close the file stream
     int open(  const string& filename, bool newOpenMode ); 
      int open(  const string& filename  ) ;
     int open();
     void close() { closefile();}
    // read and write
    int getaread(string& rdname, int& rdnumber, const char *& rdseq);
    int getaread(string& rdname, int& rdnumber, const char *& rdseq, const char*& rdqual); 
    int writearead(string& rdname, int rdnumber, const char * rdseq, const char* rdqual);
};

//fastq file
class CFasta : public CSeqfile {
   
   public:
   CFasta(bool fmode) :CSeqfile(fmode) { }
   CFasta():CSeqfile(0) {}
   CFasta( const string& seqFile, bool fmode=0):CSeqfile(seqFile, fmode) {}

   ~CFasta();
    int open(  const string& filename, bool bOpenMode) ; 
     int open(  const string& filename  ); 
    int open() ;
    void close() { closefile();}
    // read and write
    int getaread(string& rdname, int& rdnumber, const char *& rdseq) ;
    int writearead(string& rdname, int rdnumber, const char * rdseq);
    int getaread(string& rdname, int& rdnumber, const char *& rdseq, const char*& rdqual); 
    int writearead(string& rdname, int rdnumber, const char * rdseq, const char* rdqual);
};

int seqfiletype (string & seqfilename);
#endif

