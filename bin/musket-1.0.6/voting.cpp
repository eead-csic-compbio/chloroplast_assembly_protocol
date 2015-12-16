#include "voting.h"

int voteErroneousRead(vector<hmap_t*>&mykmap, int k, size_t watershed,
		char* seq, int seqLen, uint8_t votes[][4]) {
	/*check the sequence length*/
	if (seqLen < k) {
		return 0;
	}

	//initialize the votes vector
	for (int i = 0; i < seqLen; ++i) {
		for (int j = 0; j < 4; ++j) {
			votes[i][j] = 0;
		}
	}

	//cast votes for mutations
	bool errorFree = true;
	bool revcomp;
	int desWorker;
	int nworkers = mykmap.size();
	Kmer km(k, seq);
	for (int ipos = 0; ipos <= seqLen - k; ++ipos) {
		if (ipos > 0) {
			km = km.forwardBase(seq[ipos + k - 1]);
		}
		Kmer rep = km.twin();
		revcomp = true;
		if (km < rep) {
			rep = km;
			revcomp = false;
		}
		desWorker = rep.hash() % nworkers;

		/*if the k-mer is healthy*/
		if (getKmerMulti(mykmap[desWorker], rep) >= watershed) {
			continue;
		}
		errorFree = false;

		//for each offset
		for (int off = 0; off < k; ++off) {
			/*for each possible mutation*/
			int baseIndex = revcomp ? k - off - 1 : off;
			int originalBase = rep.get_base(baseIndex);

			/*get all possible alternatives*/
			for (int base = (originalBase + 1) & 3; base != originalBase;
					base = (base + 1) & 3) {
				Kmer alternativeKmer = rep.set_base(base, baseIndex);
				Kmer alternativeSearch = alternativeKmer.twin();
				if (alternativeKmer < alternativeSearch) {
					alternativeSearch = alternativeKmer;
				}
				desWorker = alternativeSearch.hash() % nworkers;
				/*if finding a trusted k-mer*/
				if (getKmerMulti(mykmap[desWorker], alternativeSearch)
						>= watershed) {
					++votes[ipos + off][revcomp ? base ^ 3 : base];
				}
			}
		}
	}
	/*error-free read*/
	if (errorFree) {
		return 0;
	}

	/*select the maximum vote*/
	int vote, maxVote = 0;
	for (int ipos = 0; ipos < seqLen; ++ipos) {
		for (int base = 0; base < 4; ++base) {
			vote = votes[ipos][base];
			if (vote > maxVote) {
				maxVote = vote;
			}
		}
	}
	return maxVote;
}

bool fixErroneousBase(int k, char* seq, int seqLen, uint8_t votes[][4],
		size_t maxVote, size_t minVote, vector<Correction>& corrections) {

	corrections.clear();
	if (maxVote < minVote || seqLen < k) {
		return false;
	}

	/*find the qualified base mutations*/
	int alternativeBase;
	for (int pos = 0; pos < seqLen; ++pos) {
		alternativeBase = -1;
		for (int base = 0; base < 4; ++base) {
			if (votes[pos][base] == maxVote) {
				if (alternativeBase == -1) {
					alternativeBase = base;
				} else {
					alternativeBase = -2;
				}
			}
		}
		if (alternativeBase >= 0) {
			corrections.push_back(Correction(pos, alternativeBase));
		}
	}
	return corrections.size() > 0;
}
bool isTrimmable(vector<hmap_t*>& mykmap, int k, size_t watershed, char* seq,
		int seqLen, int maxTrim, pair<int, int>& region) {
	bool trustedRegion = true;
	int leftKmer = -1, rightKmer = -1;
	int desWorker;
	int nworkers = mykmap.size();

	/*initialize the region*/
	region = make_pair(0, -1);

	/*check the trustiness of each k-mer*/
	Kmer km(k, seq);
	for (int ipos = 0; ipos <= seqLen - k; ++ipos) {
		if (ipos > 0) {
			km = km.forwardBase(seq[ipos + k - 1]);
		}
		Kmer rep = km.twin();
		if (km < rep) {
			rep = km;
		}
		desWorker = rep.hash() % nworkers;
		if (getKmerMulti(mykmap[desWorker], rep) >= watershed) {
			//found a healthy k-mer!
			if (trustedRegion) {
				trustedRegion = false;
				leftKmer = rightKmer = ipos;
			} else {
				++rightKmer;
			}
		} else {
			//save the trusted region
			if (leftKmer >= 0) {
				if (rightKmer < leftKmer) {
					cerr << "rightKmer: " << rightKmer << "leftKmer: "
							<< leftKmer << endl;
					exit(0);
				}
				if (rightKmer - leftKmer > region.second - region.first) {
					region = make_pair(leftKmer, rightKmer);
				}
				leftKmer = rightKmer = -1;
			}
			trustedRegion = true;
		}
	}
	if (trustedRegion == false && leftKmer >= 0) {
		if (rightKmer - leftKmer > region.second - region.first) {
			region = make_pair(leftKmer, rightKmer);
		}
	}
	if (region.second < region.first) {
		return false;
	}
	/*check the longest trusted region*/
	region.second += k;
	if (seqLen - (region.second - region.first) > maxTrim) {
		return false;
	}
	return true;
}
