#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "kseq.h"
#include "fastq.h"

FastqFile::FastqFile(const vector<string> files) {
	fnames = files;
	fnit = fnames.begin();
	file_no = 0;
	fp = gzopen(fnit->c_str(), "r");
	kseq = kseq_init(fp);
}

FastqFile::~FastqFile() {
	close();
}
void FastqFile::reset(const string file)
{
	reset(vector<string>(1, file));
}
void FastqFile::reset(const vector<string> files)
{
	/*close previous files*/
	close();

	/*release the old file name list*/
	for(vector<string>::iterator iter = fnames.begin(); iter != fnames.end(); ++iter){
		iter->clear();
	}
	fnames.clear();

	/*re-initialize the file list*/
	fnames = files;
	fnit = fnames.begin();
	file_no = 0;
	fp = gzopen(fnit->c_str(), "r");
	kseq = kseq_init(fp);
}
void FastqFile::reopen() {
	close();
	fnit = fnames.begin();
	fp = gzopen(fnit->c_str(), "r");
	kseq = kseq_init(fp);
}

// returns >=0 (length of seq), -1 end of last file, -2 truncated quality string
int FastqFile::read_next(char *read, size_t *read_len, char *seq,
		size_t *seq_len, unsigned int *file_id) {
	int r;
	r = kseq_read(kseq);
	if (r >= 0) {
		memcpy(read, kseq->name.s, kseq->name.l + 1); // 0-terminated string
		*read_len = kseq->name.l;
		memcpy(seq, kseq->seq.s, kseq->seq.l + 1); // 0-terminated string
		*seq_len = kseq->seq.l;
		if (file_id != NULL) {
			*file_id = file_no / 2;
		}
	} else if (r == -1) {
		open_next();
		if (fnit != fnames.end()) {
			return read_next(read, read_len, seq, seq_len, file_id);
		} else {
			return -1;
		}
	}
	return r;
}

// returns >=0 (length of seq), -1 end of last file, -2 truncated quality string
int FastqFile::read_next_with_qual(char *read, size_t *read_len, char *seq,
		size_t *seq_len, char *qual, size_t* qual_length,
		unsigned int *file_id) {
	int r;
	r = kseq_read(kseq);
	if (r >= 0) {
		memcpy(read, kseq->name.s, kseq->name.l + 1); // 0-terminated string
		*read_len = kseq->name.l;
		memcpy(seq, kseq->seq.s, kseq->seq.l + 1); // 0-terminated string
		*seq_len = kseq->seq.l;
		//check whether it is FASTQ format
		if (kseq->qual.l == kseq->seq.l) {
			//FASTQ format
			memcpy(qual, kseq->qual.s, kseq->qual.l + 1);
			*qual_length = kseq->qual.l;
		} else {
			//FASTA format
			*qual_length = 0;
		}
		if (file_id != NULL) {
			*file_id = file_no / 2;
		}
	} else if (r == -1) {
		open_next();
		if (fnit != fnames.end()) {
			return read_next_with_qual(read, read_len, seq, seq_len, qual,
					qual_length, file_id);
		} else {
			return -1;
		}
	}
	return r;
}

int FastqFile::read_next_batch(vector<string>& batch, size_t maxReadsPerBatch,
		bool& done, unsigned int* file_id) {
	char name[8192], seq[8192];
	size_t nameLen, seqLen;
	//clear the batch
	for(vector<string>::iterator iter = batch.begin(); iter != batch.end(); ++iter){
		iter->clear();
	}
	batch.clear();

	//read a batch of reads
	for (size_t numReads = 0; numReads < maxReadsPerBatch; numReads++) {
		if (read_next(name, &nameLen, seq, &seqLen, file_id) < 0) {
			done = true;
			break;
		}
		//save the read
		batch.push_back(string(seq, seqLen));
	}
	return batch.size();
}

int FastqFile::read_next_batch_with_name_and_qual(vector<string>& batch,
		vector<string>& names, vector<string>& qualityScores,
		size_t maxReadsPerBatch, bool& done, unsigned int* file_id) {
	char name[8192], seq[8192], qual[8192];
	size_t nameLen, seqLen, qualLen;

	//clear the batch
	for(vector<string>::iterator iter = batch.begin(); iter != batch.end(); ++iter){
		iter->clear();
	}
	batch.clear();
	for(vector<string>::iterator iter = names.begin(); iter != names.end(); ++iter){
		iter->clear();
	}
	names.clear();
	for(vector<string>::iterator iter = qualityScores.begin(); iter != qualityScores.end(); ++iter){
		iter->clear();
	}
	qualityScores.clear();

	//read a batch of reads
	for (size_t numReads = 0; numReads < maxReadsPerBatch; numReads++) {
		if (read_next_with_qual(name, &nameLen, seq, &seqLen, qual, &qualLen,
				file_id) < 0) {
			done = true;
			break;
		}
		//save the read
		batch.push_back(string(seq, seqLen));
		names.push_back(string(name, nameLen));
		if (qualLen > 0) {
			qualityScores.push_back(string(qual, seqLen));
		} else {
			qualityScores.push_back("");
		}
	}
	return batch.size();
}
vector<string>::const_iterator FastqFile::open_next() {
	if (fnit != fnames.end()) {
		// close current file
		kseq_destroy(kseq);
		gzclose(fp);
		kseq = NULL;

		// get next file
		++fnit;
		++file_no;
		if (fnit != fnames.end()) {
			fp = gzopen(fnit->c_str(), "r");
			kseq = kseq_init(fp);
		}
	}
	return fnit;
}

void FastqFile::close() {
	if (kseq != NULL) {
		kseq_destroy(kseq);
		kseq = NULL;
		gzclose(fp);
		fnit = fnames.end();
	}
}
