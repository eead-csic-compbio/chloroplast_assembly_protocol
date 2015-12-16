/*
 * paracorrect.h
 *
 *  Created on: Apr 23, 2012
 *      Author: yongchao
 */

#ifndef PARACORRECT_H_
#define PARACORRECT_H_

#include "option.h"
typedef bitset<MAX_SEQ_LENGTH> ReadBitSet;
struct Correction {
	Correction() {
		index = 0;
		base = 0;
	}
	Correction(short _index, short _base) {
		index = _index;
		base = _base;
	}
	short index;
	short base;
};
struct CorrectRegion {
	CorrectRegion() {
		leftKmer = rightKmer = -1;
	}
	CorrectRegion(int _leftKmer, int _rightKmer) {
		leftKmer = _leftKmer;
		rightKmer = _rightKmer;
	}
	inline int length() const {
		return rightKmer - leftKmer + 1;
	}
	int leftKmer;
	int rightKmer;
};

class CorrectRegionPingComp {
public:
	bool operator()(CorrectRegion& a, CorrectRegion &b) const {
		int res = a.length() - b.length();

		if (res < 0)
			return true;

		return false;
	}
};
class CorrectRegionPongComp {
public:
	bool operator()(CorrectRegion& a, CorrectRegion &b) const {
		int res = a.length() - b.length();

		if (res > 0)
			return true;

		return false;
	}
};

inline void printSeq(char* seq, int seqLen) {
	for (int i = 0; i < seqLen; ++i) {
		cerr << seq[i];
	}
	cerr << endl;
}

bool ParaCorrect(ProgramOptions& opt, vector<hmap_t*>&mykmap, int k,
		size_t watershed, char* seq, int seqLen, uint8_t votes[][4]);

#endif /* WITHOUTGAPS_H_ */
