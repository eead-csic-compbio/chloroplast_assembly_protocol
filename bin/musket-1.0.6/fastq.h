#ifndef FASTQ_H
#define FASTQ_H

#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include <vector>
#include <string>

#include "kseq.h"

// initialize kseq structures
//KSEQ_INIT(gzFile, gzread);

using namespace std;

class FastqFile {
public:
	FastqFile(const vector<string> fnames);

	~FastqFile();

	void reset(const string file);
	void reset(const vector<string> fnames);
	void close();
	void reopen();
	int read_next(char *read, size_t *read_len, char *seq, size_t *seq_len,
			unsigned int *file_id);
	int
	read_next_with_qual(char *read, size_t *read_len, char *seq,
			size_t *seq_len, char *qual, size_t* qual_length,
			unsigned int *file_id);
	int
	read_next_batch(vector<string>& batch, size_t maxReadsPerBatch, bool& done,
			unsigned int* file_id);
	int
	read_next_batch_with_name_and_qual(vector<string>& batch,
			vector<string>& names, vector<string>& qualityScores,
			size_t maxReadsPerBatch, bool& done, unsigned int* file_id);
private:
	vector<string>::const_iterator
	open_next();

	vector<string>::const_iterator fnit;
	unsigned int file_no;

	/*input file list*/
	vector<string> fnames;

	/*file pointer*/
	gzFile fp;

	/*sequence pointer*/
	kseq_t *kseq;
};

#endif
