#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include "sys/time.h"
using namespace std;

//help message
#define HELP "\nUsage: duk [options]  ref.fa query \n\
   \n\
        Options:  \n\
	    -o, -output  <file>       print log information to file, default is stdout. \n\
            -n, -nomatch <file>       output the not matched reads to file,  \n\
	                               the opton value - stands for standard output  \n\
	    -m, -match <file>      output matched reads to file. \n\
	    -k, -kmer                 the k mer size,  default is 16.  \n\
	    -s, -step                 the step size, default is 4. \n\
	    -c, -cutoff               the cut off threshold for matched reads, default is 1.\n\
            -h, -help                 print out the help information\n\
         \n\n\
	Identify the reads in query file whether they match to ref.fa reads. \n\
	The ref.fa must be in fasta format and the query can be in fastq or fasta format. \n\
	If there is no query file, the tool gets reads from standard input \n\n"

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
/* 
  parse the command auguments and save the results in the Parameters struct. 

 */
int parseCmdOpts(Parameters& param, int argc, char ** argv) {
   //initialize the defaults
   param.merSize = 16;
   param.stepSize = 4;
   param.cutoff = 1;
   param.noMatchRdsFile = "";
   param.matchRdsFile = "";
   param.logFile = "";

   int c;
   int i;
   static int help_flag = 0;
   while (1) {
     static struct option long_options[] = {
        {"output", required_argument, 0, 'o'},
        {"nomatch", required_argument, 0, 'n'},
        {"match", required_argument, 0, 'm'},
        {"step", required_argument, 0, 's'},
        {"cutoff", required_argument, 0, 'c'},
        {"kmer", required_argument, 0, 'k'},
        {"help", no_argument, &help_flag, 1},
        {0, 0, 0, 0}
     };

     int option_index = 0;
     c = getopt_long_only (argc, argv, "o:n:m:s:c:hk:",
                      long_options, &option_index);
     if (c == -1) 
        break;

     switch (c) {
        case 0:
	  break;

        case 'h':
	   help_flag = 1; 
	   break;
        
        case 'o':
           param.logFile = optarg;
	   break;
        
	case 'n':
           param.noMatchRdsFile = optarg;
	   break;

	case 'm':
           param.matchRdsFile = optarg;
	   break;

	case 's':
	  i = atoi(optarg);
           if (!i ) {
	      cout << "step size must be an integer " << endl;
	      exit(0);
	    }
           if (i < 1 ) {
	      cout << "step size must be larger than 0" << endl;
	      exit(0);
	    }
	    param.stepSize = i;
	   break;

	case 'c':
	  i = atoi(optarg);
           if (!i ) {
	      cout << "Cut off value  must be an integer " << endl;
	      exit(0);
	    }
           if (i < 1 ) {
	      cout << "Cut off value  must be larger than 0" << endl;
	      exit(0);
	    }
	    param.cutoff = i;
	   break;
	
	case 'k':
	  i = atoi(optarg);
           if (!i ) {
	      cout << "Kmer size  must be an integer " << endl;
	      exit(0);
	    }
           if (i < 1 ) {
	      cout << "Kmer size  must be larger than 0" << endl;
	      exit(0);
	    }
	    param.merSize = i;
	   break;

        case '?':
          cout << HELP << endl;
	  exit(0);
	  break;

	default:
	  exit(0);

     }
    
    if (help_flag) { 
      cout << HELP << endl;
      exit(0);
    }

   }
  
  // contamination file
  if (optind >= argc) {
     cout << "Reference file must be supplied at fasta format" << endl << endl;
     cout <<HELP << endl;
     exit(0);
  } else {
     param.dbFile = argv[optind++];
  }
  // query file
  if (optind >= argc) {
     param.queryFile = "";
  } else {
     param.queryFile = argv[optind++];
  }

 return 0;
}
