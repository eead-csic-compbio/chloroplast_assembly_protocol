/*
 * Option.h
 *
 *  Created on: Apr 23, 2012
 *      Author: yongchao
 */

#ifndef OPTION_H_
#define OPTION_H_
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <functional>
#include "thread.h"
#include <omp.h>
#include <bitset>
#include <map>

using namespace std;

enum {
	BASE_NONE, BASE_UNIQUE,
/*introduced a unique substitution*/
};
#define VERSION "1.0.6"

#ifndef MAX_SEQ_LENGTH
#define MAX_SEQ_LENGTH					255
#endif

#define MAX_SEQ_BUFFER_SIZE			(MAX_SEQ_LENGTH + 1)

struct ProgramOptions {
	ProgramOptions() {
		klist.reserve(4);
		watershed.reserve(4);
		nthreads = 2;
		ratio = 4;
		nmasters = 1;
		nworkers = 0;
		numErrors = 4;
		maxTrim = 0;
		maxIters = 2;
		enableMK = false;
		keepOrder = false;
		lowercase = false;
		maxKmerPerSlot = 1024;
		maxReadsPerBatch = 4096;
		correctedFileName = "";
		files.clear();
		outfiles.clear();
	}
	/*the list of kmer size*/
	vector<pair<int, unsigned int> > klist;
	size_t nmasters; //the number of masters
	size_t nworkers; //the number of work thread
	size_t nthreads; //the number of threads
	int numErrors; //maximum number of mutaitons per k-mer
	int maxTrim; //maximal number of bases to be trimmed
	int maxIters; //maximum iterations for each k-mer
	size_t ratio;
	bool enableMK;
	bool keepOrder;
	bool lowercase;	//corrected bases are written in lowercase
	vector<size_t> watershed;
	size_t maxKmerPerSlot;
	size_t maxReadsPerBatch;
	string correctedFileName;
	vector<string> files;	/*input file names*/
	vector<string> outfiles;	/*output file names*/
};

void ParseOptions(int argc, char **argv, ProgramOptions &opt);

/*get the multiplicity of the k-mer*/
inline size_t getKmerMulti(hmap_t *kmap, Kmer& km) {
	hmap_t::iterator iter = kmap->find(km);
	if (iter != kmap->end()) {
		return iter->GetVal();
	}
	return 0;
}

inline double getSysTime() {
	double dtime;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	dtime = (double) tv.tv_sec;
	dtime = dtime + (double) (tv.tv_usec) / 1000000.0;

	return dtime;
}

#endif /* OPTION_H_ */
