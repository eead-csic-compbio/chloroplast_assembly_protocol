/*
 * Voting.h
 *
 *  Created on: Apr 23, 2012
 *      Author: yongchao
 */

#ifndef VOTING_H_
#define VOTING_H_
#include "option.h"
#include "paracorrect.h"

int voteErroneousRead(vector<hmap_t*>&mykmap, int k, size_t watershed,
		char* seq, int seqLen, uint8_t votes[][4]);
bool fixErroneousBase(int k, char* seq, int seqLen, uint8_t votes[][4],
		size_t maxVote, size_t minVote, vector<Correction>& corrections);

bool isTrimmable(vector<hmap_t*>&mykmap, int k, size_t wateshed, char* seq,
		int seqLen, int maxTrim, pair<int, int>& region);

#endif /* VOTING_H_ */
