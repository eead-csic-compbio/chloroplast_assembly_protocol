/*
 * Option.cpp
 *
 *  Created on: Apr 23, 2012
 *      Author: yongchao
 */
#include "option.h"
#include <map>

using namespace std;

void PrintUsage(ProgramOptions& opt) {
	cerr
			<< endl
			<< "MUSKET (version "
			<< VERSION << ") is a parallel multi-stage k-mer based error corrector"
			<< endl;

	cerr << "Usage: musket [options] file [file1 ...]" << endl;
	cerr
			<< "Basic Options:"
			<< endl
			<< "\t-k <int uint> (specify two paramters: k-mer size and estimated total number of k-mers for this k-mer size)"
			<< endl
			<< "\t   (e.g. estimated number of k-mers: 67108864, 134217728, 268435456 and 536870912)"
			<< endl
			<< "\t-o <str> (the single output file name)"
			<< endl
			<< "\t-omulti <str> (prefix of output file names, one input corresponding one output)"
			<< endl
			<< "\t-p <int> (number of threads [>=2], default " << opt.nthreads
			<< ")" << endl
			<< "\t-maxtrim <int> (maximal number of bases that can be trimmed, default "
			<< opt.maxTrim << ")" << endl
			<< "\t-inorder (keep sequences outputed in the same order with the input)"
			<< endl
			<< "\t-lowercase (write corrected bases in lowercase, default "
			<< opt.lowercase << ")" <<endl;

	cerr
			<< "Advanced:"
			<< endl
			<< "\t-maxbuff <int> (capacity of message buffer for each worker, default "
			<< opt.maxKmerPerSlot
			<< ")"
			<< endl
			<< "\t-multik <bool> (enable the use of multiple k-mer sizes, default "
			<< opt.enableMK
			<< ")"
			<< endl
			<< "\t-maxerr <int> (maximal number of mutations in any region of length #k, default "
			<< opt.numErrors
			<< ")"
			<< endl
			<< "\t-maxiter <int> (maximal number of correcting iterations per k-mer size, default "
			<< opt.maxIters << ")" << endl;
}

void ParseOptions(int argc, char **argv, ProgramOptions &opt) {
	int value;
	size_t uvalue;
	string correctedFileNamePrefix = "";
	int index = 1;
	if (argc < 2) {
		PrintUsage(opt);
		exit(1);
	}

	/*remove repeated kmer-size*/
	map<int, size_t> kmerSizes;
	while (index < argc) {
		if (argv[index][0] != '-'){
			/*assume that these are input file names*/
			opt.files.push_back(argv[index]);
			++index;
			continue;
		}
		if (!strcmp(argv[index], "-k")) {
			index++;
			if (index + 2 >= argc) {
				cerr << "Not specify sufficient values for option "
						<< argv[index - 1] << endl;
				PrintUsage(opt);
				exit(1);
			}

			/*get the k-mer size*/
			if (argv[index][0] != '-') {
				sscanf(argv[index], "%d", &value);
				if (value >= MAX_KMER_SIZE) {
					value = MAX_KMER_SIZE - 1;
				}
				if (value < 13) {
					value = 13;
				}
			} else {
				cerr << "Two values for option -k" << endl;
				exit(1);
			}
			/*get the total number of k-mers*/
			index++;
			if (argv[index][0] != '-') {
				sscanf(argv[index], "%lu", &uvalue);
				if (uvalue < 1000) {
					uvalue = 1000;
				}
			} else {
				cerr << "Two values for option -k" << endl;
				exit(1);
			}
			index++;

			/*save the two values*/
			kmerSizes.insert(make_pair(value, uvalue));
		} else if (!strcmp(argv[index], "-multik")) {
			index++;
			if (index < argc) {
				sscanf(argv[index], "%d", &value);
				if (value != 0) {
					value = 1;
				}
				opt.enableMK = value;
			} else {
				cerr << "not specify value for option " << argv[index - 1]
						<< endl;
				PrintUsage(opt);
				exit(1);
			}
			index++;
		} else if (!strcmp(argv[index], "-p")) {
			index++;
			if (index < argc) {
				sscanf(argv[index], "%d", &value);
				if (value < 2) {
					value = 2;
				}
				opt.nthreads = value;
			} else {
				cerr << "not specify value for option " << argv[index - 1]
						<< endl;
				PrintUsage(opt);
				exit(1);
			}
			index++;
		} else if (!strcmp(argv[index], "-maxtrim")) {
			index++;
			if (index < argc) {
				sscanf(argv[index], "%d", &value);
				if (value < 0) {
					value = 0;
				}
				opt.maxTrim = value;
			} else {
				cerr << "not specify value for option " << argv[index - 1]
						<< endl;
				PrintUsage(opt);
				exit(1);
			}
			index++;
		} else if (!strcmp(argv[index], "-maxiter")) {
			index++;
			if (index < argc) {
				sscanf(argv[index], "%d", &value);
				if (value < 1) {
					value = 1;
				}
				opt.maxIters = value;
			} else {
				cerr << "not specify value for option " << argv[index - 1]
						<< endl;
				PrintUsage(opt);
				exit(1);
			}
			index++;
		} else if (!strcmp(argv[index], "-o")) {
			index++;
			if (index < argc) {
				stringstream ss(argv[index]);
				if ((ss >> opt.correctedFileName).fail()) {
					cerr << "corrected reads output file name failed" << endl;
					PrintUsage(opt);
					exit(1);
				}
			} else {
				cerr << "not specify value for option " << argv[index - 1]
						<< endl;
				PrintUsage(opt);
				exit(1);
			}
			index++;
    } else if (!strcmp(argv[index], "-omulti")) {
      index++;
      if (index < argc) {
        stringstream ss(argv[index]);
        if ((ss >> correctedFileNamePrefix).fail()) {
          cerr << "corrected reads file name prefix failed" << endl;
          PrintUsage(opt);
          exit(1);
        }
      } else {
        cerr << "not specify value for option " << argv[index - 1]
            << endl;
        PrintUsage(opt);
        exit(1);
      }
      index++;
		} else if (!strcmp(argv[index], "-maxbuff")) {
			index++;
			if (index < argc) {
				sscanf(argv[index], "%d", &value);
				if (value < 16)
					value = 16;
				opt.maxKmerPerSlot = value;
			} else {
				cerr << "not specify value for option " << argv[index - 1]
						<< endl;
				PrintUsage(opt);
				exit(1);
			}
			index++;
		} else if (!strcmp(argv[index], "-maxerr")) {
			index++;
			if (index < argc) {
				sscanf(argv[index], "%d", &value);
				if (value < 2)
					value = 2;
				opt.numErrors = value;
			} else {
				cerr << "not specify value for option " << argv[index - 1]
						<< endl;
				PrintUsage(opt);
				exit(1);
			}
			index++;
		} else if (!strcmp(argv[index], "-inorder") || !strcmp(argv[index], "-paired")) {
			opt.keepOrder = true;
			++index;
		}else if(!strcmp(argv[index], "-lowercase")){
			opt.lowercase = true;
			++index;
		} else {
			cerr << "Unknown option " << argv[index] << endl;
			PrintUsage(opt);
			exit(1);
		}
	}
	/*check the input file list*/
	if(opt.files.size() == 0){
		cerr << "Not specify any input file" << endl;
		PrintUsage(opt);
		exit(1);
	}

	/*if keeping the order*/
	if (opt.keepOrder) {
		if (opt.maxReadsPerBatch < 256 * 1024) {
			opt.maxReadsPerBatch = 256 * 1024;
		}
	}
	/*check the output file*/
	if (opt.correctedFileName.length() == 0 && correctedFileNamePrefix.length() == 0) {
		cerr << "Must specif the output using either option \"-o\" or \"-omulti\"" << endl;
		exit(0);
	}
	if(correctedFileNamePrefix.length() > 0 && opt.correctedFileName.length() > 0){
		cerr << "Can only use one of options \"-o\" or \"-omulti\"" << endl;
		exit(0);
	}
	/*generate the output file list*/
	if(correctedFileNamePrefix.length() > 0){
		char buffer[1024];
		opt.outfiles.reserve(opt.files.size());
		for(size_t i = 0; i < opt.files.size(); ++i){
			string name = correctedFileNamePrefix + ".";
			sprintf(buffer, "%ld", i);
			name += buffer;
			opt.outfiles.push_back(name);
		}
		correctedFileNamePrefix.clear();
	}

	/*check the availability of k-mer sizes*/
	for (map<int, size_t>::reverse_iterator iter = kmerSizes.rbegin();
			iter != kmerSizes.rend(); ++iter) {
		opt.klist.push_back(make_pair(iter->first, iter->second));
	}
	kmerSizes.clear();

	if (opt.klist.size() == 0) {
		cerr << "Not specify the k-mer size and will use the default ones"
				<< endl;
		if (opt.enableMK) {
			opt.klist.push_back(make_pair(21, 536870912));
			opt.klist.push_back(make_pair(23, 536870912));
		} else {
			opt.klist.push_back(make_pair(21, 536870912));
		}
	} else {
		if (!opt.enableMK) {
			if (opt.klist.size() > 1) {
				cerr
						<< "Only the first k-mer size is used since the multiple k-mer size mode is disabled"
						<< endl;
			}
			opt.klist.resize(1);
		}
	}
	opt.watershed.resize(opt.klist.size(), 0);

	cerr << "The kmer sizes that will be used are as follows: " << endl;
	for (size_t i = 0; i < opt.klist.size(); ++i) {
		cerr << "\tKmer size: " << opt.klist[i].first
				<< " and estimated number of k-mers: " << opt.klist[i].second
				<< endl;
	}

	/*check the availability of the input files*/
	struct stat stFileInfo;
	vector<string>::iterator it;
	int intStat;
	for (it = opt.files.begin(); it != opt.files.end(); ++it) {
		intStat = stat(it->c_str(), &stFileInfo);
		if (intStat != 0) {
			cerr << "Error: file not found, " << *it << endl;
			PrintUsage(opt);
			exit(1);
		}
	}

	//calculate the number of masters and workers
	opt.nmasters = opt.nthreads / opt.ratio;
	if (opt.nmasters < 1)
		opt.nmasters = 1;

	opt.nworkers = opt.nthreads - opt.nmasters;
	cerr << "#master threads: " << opt.nmasters << " #worker threads: "
			<< opt.nworkers << endl;

	/*estimate the total number of k-mers per worker*/
	for (size_t i = 0; i < opt.klist.size(); ++i) {
		opt.klist[i].second = (opt.klist[i].second + opt.nworkers - 1)
				/ opt.nworkers;
	}

	/*
	 for(int i = 0; i <= 127; ++i){
	 int ch = toupper(i);
	 switch(ch){
	 case 'A':
	 cerr << "0,";break;
	 case 'C':
	 cerr << "1,";break;
	 case 'G': cerr << "2, "; break;
	 case 'T': cerr << "3, "; break;
	 default: cerr << "4, "; break;
	 }
	 if(i && i % 32 == 0){
	 cerr << endl;
	 }
	 }*/
}

