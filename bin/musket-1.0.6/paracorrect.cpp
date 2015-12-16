/*n
 * paracorrect.cpp
 *
 *  Created on: Apr 23, 2012
 *      Author: yongchao
 */
#include "paracorrect.h"
#include "voting.h"

/*check the correctness of its successing k-mer*/
static bool successor(vector<hmap_t*>&mykmap, size_t watershed, char* seq,
		int seqLen, Kmer km, int dist) {
	int k = km.k;
	if (seqLen < k || dist <= 0) {
		return true;
	}

	int nworkers = mykmap.size();
	int endPos = min(seqLen - k, dist - 1);
	for (int ipos = 0; ipos <= endPos; ++ipos) {
		km = km.forwardBase(seq[ipos + k - 1]);
		Kmer rep = km.twin();
		if (km < rep) {
			rep = km;
		}
		int desWorker = rep.hash() % nworkers;
		if (getKmerMulti(mykmap[desWorker], rep) < watershed) {
			return false;
		}
	}
	return true;
}
/*check the correctness of its precessing k-mer*/
static bool predecessor(vector<hmap_t*>&mykmap, size_t watershed, char* seq,
		int seqLen, Kmer km, int dist) {
	if (seqLen <= 0 || dist <= 0) {
		return true;
	}
	int nworkers = mykmap.size();
	int startPos = max(0, seqLen - dist);
	for (int ipos = seqLen - 1; ipos >= startPos; --ipos) {
		km = km.backwardBase(seq[ipos]);
		Kmer rep = km.twin();
		if (km < rep) {
			rep = km;
		}
		int desWorker = rep.hash() % nworkers;
		if (getKmerMulti(mykmap[desWorker], rep) < watershed) {
			return false;
		}
	}
	return true;
}

static int selectUniqueMutation(vector<hmap_t*>& mykmap, int k,
		size_t watershed, Kmer km, bool revcomp, int index) {

	int alternativeBase = -1;
	int nworkers = mykmap.size();
	int baseIndex = revcomp ? k - 1 - index : index;
	int originalBase = km.get_base(baseIndex);
	assert(index >= 0 && index < k);

	/*get all possible alternatives*/
	for (int base = (originalBase + 1) & 3; base != originalBase;
			base = (base + 1) & 3) {
		Kmer alternativeKmer = km.set_base(base, baseIndex);
		Kmer alternativeSearch = alternativeKmer.twin();
		if (alternativeKmer < alternativeSearch) {
			alternativeSearch = alternativeKmer;
		}
		int desWorker = alternativeSearch.hash() % nworkers;
		if (getKmerMulti(mykmap[desWorker], alternativeSearch) >= watershed) {
			/*save the base*/
			if (alternativeBase != -1) {
				return -1;
			}
			alternativeBase = base;
		}
	}
	return revcomp ? alternativeBase ^ 3 : alternativeBase;
}
static int selectAllMutations(vector<hmap_t*>& mykmap, int k, size_t watershed,
		Kmer km, bool revcomp, int index, pair<int, int> bases[4]) {
	int nworkers = mykmap.size();
	/*set the number of possible mutations*/
	int numBases = 0;

	/*get the original base at index*/
	int baseIndex = revcomp ? k - 1 - index : index;
	int originalBase = km.get_base(baseIndex);
	assert(index >= 0 && index < k);

	//cerr << "revcomp: " << revcomp << " kmer: " << km.toString() << endl;
	//cerr << "baseIndex: " << baseIndex << " index: " << index << endl;

	/*get all possible alternatives*/
	for (int base = (originalBase + 1) & 3; base != originalBase;
			base = (base + 1) & 3) {
		Kmer alternativeKmer = km.set_base(base, baseIndex);
		Kmer alternativeSearch = alternativeKmer.twin();
		if (alternativeKmer < alternativeSearch) {
			alternativeSearch = alternativeKmer;
		}
		int desWorker = alternativeSearch.hash() % nworkers;
		size_t multi = getKmerMulti(mykmap[desWorker], alternativeSearch);
		if (multi >= watershed) {
			bases[numBases++] = make_pair(revcomp ? base ^ 3 : base, multi);
		}
	}
	return numBases;
}
/*error correcting that only allows <= 1 errors in any k-mer of a read*/
static bool paraCorrectCoreTwoSides(vector<hmap_t*>&mykmap, int k,
		size_t watershed, char* seq, int seqLen,
		vector<Correction>& corrections) {
	ReadBitSet solids;
	int numBases, lnumBases, ipos;
	int nworkers = mykmap.size();
	pair<int, int> bases[4], lbases[4];

	/*check the sequence length*/
	corrections.clear();
	if (seqLen < k) {
		return true; /*too short do not correct it any more*/
	}

	/*iterate each k-mer to check its solidity*/
	Kmer km(k, seq);
	solids.reset();
	for (ipos = 0; ipos <= seqLen - k; ++ipos) {
		if (ipos > 0) {
			km = km.forwardBase(seq[ipos + k - 1]);
		}
		Kmer rep = km.twin();
		if (km < rep) {
			rep = km;
		}
		int desWorker = rep.hash() % nworkers;
		if (getKmerMulti(mykmap[desWorker], rep) >= watershed) {
			for (int i = ipos; i < ipos + k; ++i) {
				solids[i] = 1;
			}
		}
	}
	/*if the read is error-free*/
	if (solids.count() == (size_t) seqLen) {
		return true;
	}

	/*fix the unique errors relying on k-mer neighboring information*/
	km.set_kmer(seq);
	Kmer lkm = km;

	/*for bases from index 0 to index seqLen - k*/
	for (ipos = 0; ipos <= seqLen - k; ++ipos) {
		if (ipos > 0) {
			km = km.forwardBase(seq[ipos + k - 1]);
		}
		if (ipos >= k) {
			lkm = lkm.forwardBase(seq[ipos]);
		}
		if (solids[ipos] == 1) {
			continue;
		}

		/*try to find all mutations for the rightmost k-mer*/
		Kmer rep = km.twin();
		bool revcomp = true;
		if (km < rep) {
			rep = km;
			revcomp = false;
		}
		numBases = selectAllMutations(mykmap, k, watershed, rep, revcomp, 0,
				bases);

		/*try to find mutations for the leftmost k-mer*/
		rep = lkm.twin();
		revcomp = true;
		if (lkm < rep) {
			rep = lkm;
			revcomp = false;
		}
		lnumBases = selectAllMutations(mykmap, k, watershed, rep, revcomp,
				ipos >= k ? k - 1 : ipos, lbases);

		/*if the two bases are the same, this base is modified*/
		for (int i = 0; i < numBases; ++i) {
			int base = bases[i].first;
			for (int j = 0; j < lnumBases; ++j) {
				int lbase = lbases[j].first;
				if (base == lbase) {
					corrections.push_back(Correction(ipos, base));
				}
			}
		}
		/*check if some correction is found*/
		if (corrections.size() == 1) {
			return false;
		} else if (corrections.size() > 1) {
			corrections.clear();
		}
	}

	/*for the remaining k-1 bases*/
	Kmer rep = km.twin();
	bool revcomp = true;
	if (km < rep) {
		rep = km;
		revcomp = false;
	}
	for (; ipos < seqLen; ++ipos) {
		lkm = lkm.forwardBase(seq[ipos]);
		if (solids[ipos] == 1) {
			continue;
		}
		/*try to fix using the current k-mer*/
		numBases = selectAllMutations(mykmap, k, watershed, rep, revcomp,
				ipos - (seqLen - k), bases);

		/*try to fix using the left k-mer*/
		Kmer lrep = lkm.twin();
		bool lrevcomp = true;
		if (lkm < lrep) {
			lrep = lkm;
			lrevcomp = false;
		}
		lnumBases = selectAllMutations(mykmap, k, watershed, lrep, lrevcomp,
				k - 1, lbases);

		/*if the two bases are the same, this base is modified*/
		for (int i = 0; i < numBases; ++i) {
			int base = bases[i].first;
			for (int j = 0; j < lnumBases; ++j) {
				int lbase = lbases[j].first;
				if (base == lbase) {
					corrections.push_back(Correction(ipos, base));
				}
			}
		}
		/*check if some correction is found*/
		if (corrections.size() == 1) {
			return false;
		} else if (corrections.size() > 1) {
			corrections.clear();
		}
	}

	/*not error-free*/
	return false;
}
/*core function of error correction without gaps*/
static bool paraCorrectCore(vector<hmap_t*>&mykmap, int k, size_t watershed,
		char* seq, int seqLen, int maxErrorPerKmer,
		vector<Correction>& outcorrections, int stride, int pingpong) {

	/*check the sequence length*/
	outcorrections.clear();
	if (seqLen < k) {
		return false;
	}

	vector<Correction> corrections;
	vector<CorrectRegion> regions;
	CorrectRegionPingComp pingcomp;
	CorrectRegionPongComp pongcomp;

	int numBases;
	pair<int, int> bases[4];
	int desWorker;
	int leftKmer = -1, rightKmer = -1;
	int nworkers = mykmap.size();
	bool solidRegion = false;

	ReadBitSet solids;

	Kmer km(k, seq);
	solids.reset();
	for (int ipos = 0; ipos <= seqLen - k; ++ipos) {
		if (ipos > 0) {
			km = km.forwardBase(seq[ipos + k - 1]);
		}
		Kmer rep = km.twin();
		if (km < rep) {
			rep = km;
		}
		desWorker = rep.hash() % nworkers;
		//cerr << getKmerMulti(mykmap[desWorker], rep) << " " << km.toString() << " " << km.twin().toString() << endl;
		if (getKmerMulti(mykmap[desWorker], rep) >= watershed) {
			//found a healthy k-mer!
			if (!solidRegion) {
				solidRegion = true;
				leftKmer = rightKmer = ipos;
			} else {
				++rightKmer;
			}
			for (int i = ipos; i < ipos + k; ++i) {
				solids[i] = 1;
			}
		} else {
			//save the trusted region
			if (leftKmer >= 0) {
				if (rightKmer < leftKmer) {
					cerr << "rightKmer: " << rightKmer << "leftKmer: "
							<< leftKmer << endl;
					exit(0);
				}
				regions.push_back(CorrectRegion(leftKmer, rightKmer));
				if(pingpong){
					push_heap(regions.begin(), regions.end(), pingcomp);
				}else{
					push_heap(regions.begin(), regions.end(), pongcomp);
				}
				leftKmer = rightKmer = -1;
			}
			solidRegion = false;
		}
	}
	if (solidRegion == true && leftKmer >= 0) {
		regions.push_back(CorrectRegion(leftKmer, rightKmer));
    if(pingpong){
    	push_heap(regions.begin(), regions.end(), pingcomp);
   	}else{
    	push_heap(regions.begin(), regions.end(), pongcomp);
   	}
	}

	/*calculate the minimal number of votes per base*/
	if (regions.size() == 0) {
		/*this read is error-free*/
		return true;
	}

#if 0
	cerr << solids.count() << "   -----------------------------------" << endl;
	for(size_t i = 0; i < regions.size(); ++i) {
		cerr << regions[i].leftKmer << " " << regions[i].rightKmer << " " << regions[i].rightKmer - regions[i].leftKmer + 1<< endl;
	}
	getchar();
#endif

	bool hasChanged = false;
	int lastPosition, numCorrectionsPerKmer;
	while (regions.size() > 0) {

		/*get the region with the highest priority*/
		CorrectRegion region = regions[0];
		if(pingpong){
			pop_heap(regions.begin(), regions.end(), pingcomp);
		}else{
			pop_heap(regions.begin(), regions.end(), pongcomp);
		}
		regions.pop_back();

		/*get the k-mer range*/
		leftKmer = region.leftKmer;
		rightKmer = region.rightKmer;
		//cerr <<"left: " << leftKmer << "right: " << rightKmer << endl;

		/*form the starting k-mer*/
		km.set_kmer(&seq[rightKmer]);
		lastPosition = -1;
		numCorrectionsPerKmer = 0;
		corrections.clear();
		for (int ipos = rightKmer + 1; ipos <= seqLen - k; ++ipos) {
			int targetPos = ipos + k - 1;
			/*check the solids of the k-mer*/
			if (solids[targetPos]) {
				/*if it reaches another trusted region*/
				break;
			}
			if (targetPos >= seqLen) {
				cerr << "Targetpos exceeds the maximal index" << targetPos
						<< " / " << seqLen << endl;
				exit(0);
			}
			km = km.forwardBase(seq[targetPos]);
			Kmer rep = km.twin();
			bool revcomp = true;
			if (km < rep) {
				rep = km;
				revcomp = false;
			}

			desWorker = rep.hash() % nworkers;
			//cerr << "left ipos: " << ipos << " " << km.toString() << " " << getKmerMulti(mykmap[desWorker], rep) << endl;
			if (getKmerMulti(mykmap[desWorker], rep) < watershed) {
				/*select all possible mutations*/
				numBases = selectAllMutations(mykmap, k, watershed, rep,
						revcomp, k - 1, bases);

				/*start correcting*/
				bool done = false;
				//cerr << "right numBases: " << numBases << km.toString() << endl;
				if (numBases == 1) {
					int base = bases[0].first;
					corrections.push_back(Correction(targetPos, base));
					/*set the last base*/
					km = km.set_base(base, k - 1);
					done = true;
#ifdef DEBUG
					cerr << "km: " << km.toString() << " tpos: " << targetPos << " numBases: " << numBases << " multi: " << bases[0].second << endl;
#endif
				} else {
					/*select the best substitution*/
					pair<int, int> best = make_pair(0, 0);
					for (int i = 0; i < numBases; ++i) {
						int base = bases[i].first;
						if (successor(mykmap, watershed, seq + ipos + 1,
								seqLen - ipos - 1, km.set_base(base, k - 1),
								stride)) {
							/*check the multiplicity*/
							if (best.second < bases[i].second) {
									best = bases[i];
							}
						}
					}
					/*if finding a best one*/
					if (best.second > 0) {
						corrections.push_back(
								Correction(targetPos, best.first));
						km = km.set_base(best.first, k - 1);
						done = true;
#ifdef DEBUG
						cerr << "km: " << km.toString() << " pos: " << targetPos << endl;
#endif
					}
				}
				/*if finding one correction*/
				if (done) {
					/*recording the position*/
					if (lastPosition < 0) {
						lastPosition = targetPos;
					}
					/*check the number of errors in any k-mer length*/
					if (targetPos - lastPosition < k) {
						/*increase the number of corrections*/
						++numCorrectionsPerKmer;
						if (numCorrectionsPerKmer > maxErrorPerKmer) {
							for (int i = 0; i < numCorrectionsPerKmer; ++i) {
								corrections.pop_back();
							}
							break;
						}
					} else {
						lastPosition = targetPos;
						numCorrectionsPerKmer = 0;
					}
					continue;
				}
				break;
			}
		} /*for i*/

		/*save the corrections in this region*/
		for (size_t i = 0; i < corrections.size(); ++i) {
			outcorrections.push_back(corrections[i]);
		}
		hasChanged |= corrections.size() > 0;

		/*towards the beginning of the sequence from this position*/
		lastPosition = corrections.size() > 0 ? corrections[0].index : -1;
		numCorrectionsPerKmer = 0;
		corrections.clear();
		if (leftKmer > 0) {
			int ipos;
			km.set_kmer(&seq[leftKmer]);
			//cerr << km.toString() << endl;
			for (ipos = leftKmer - 1; ipos >= 0; --ipos) {
				//if(ipos - k + 1 >= 0 && solids[ipos -k + 1]){
				if (solids[ipos]) {
					/*if it reaches another trusted region*/
					break;
				}
				km = km.backwardBase(seq[ipos]);
				Kmer rep = km.twin();
				bool revcomp = true;
				if (km < rep) {
					rep = km;
					revcomp = false;
				}
				//cerr << "left ipos: " << ipos << " " << km.toString() << endl;
				desWorker = rep.hash() % nworkers;
				if (getKmerMulti(mykmap[desWorker], rep) < watershed) {
					/*select all possible mutations*/
					numBases = selectAllMutations(mykmap, k, watershed, rep,
							revcomp, 0, bases);

					/*start correcting*/
					bool done = false;
					//cerr << "left numBases: " << numBases << km.toString() << endl;
					if (numBases == 1) {
						int base = bases[0].first;
						corrections.push_back(Correction(ipos, base));
						km = km.set_base(base, 0);
						done = true;
#ifdef DEBUG
						cerr << "km: " << km.toString() << " pos: " << ipos << endl;
#endif
					} else {
						pair<int, int> best = make_pair(0, 0);
						for (int i = 0; i < numBases; ++i) {
							int base = bases[i].first;
							if (predecessor(mykmap, watershed, seq, ipos - 1,
									km.set_base(base, 0), stride)) {
								if (best.second < bases[i].second) {
									best = bases[i];
								}
							}
						}
						if (best.second > 0) {
							corrections.push_back(
									Correction(ipos, best.first));
							km = km.set_base(best.first, 0);
							done = true;
#ifdef DEBUG
							cerr << "km: " << km.toString() << " ipos: " << ipos << " numBases: " << numBases << " multi: " << maxMulti << endl;
#endif
						}
					}
					if (done) {
						if (lastPosition < 0) {
							lastPosition = ipos;
						}
						if (lastPosition - ipos < k) {
							++numCorrectionsPerKmer;
							if (numCorrectionsPerKmer > maxErrorPerKmer) {
								for (int i = 0; i < numCorrectionsPerKmer;
										++i) {
									corrections.pop_back();
								}
								break;
							}
						} else {
							lastPosition = ipos;
							numCorrectionsPerKmer = 0;
						}
						continue;
					}
					break; /*check the next region*/
				}
			}
		}/*if leftKmer > 0*/

		/*save the corrections in this region*/
		for (size_t i = 0; i < corrections.size(); ++i) {
			outcorrections.push_back(corrections[i]);
		}
		hasChanged |= corrections.size() > 0;
	}
	return false;
}
bool ParaCorrect(ProgramOptions& opt, vector<hmap_t*>&mykmap, int k,
		size_t watershed, char* seq, int seqLen, uint8_t votes[][4]) {

	int maxVote;
	vector<Correction> corrections;

	/*stage 1: conservative correcting using two-sided correcting*/
	for (int iter = 0; iter < opt.maxIters; ++iter) {
		if (paraCorrectCoreTwoSides(mykmap, k, watershed, seq, seqLen,
				corrections)) {
			return true;
		}
		/*make changes to the input sequence*/
		if (corrections.size() == 0)
			break;
		for (size_t i = 0; i < corrections.size(); ++i) {
			Correction& corr = corrections[i];
			seq[corr.index] = opt.lowercase ? Kmer::decodeToLower[corr.base] : Kmer::decodeToUpper[corr.base];
		}
	}

	/*stage 2: aggressive correcting from one side*/
	for (int nerr = 1; nerr <= opt.numErrors; ++nerr) {
		int pingpong = 1;
		for (int iter = 0; iter < opt.maxIters; ++iter) {
			if (paraCorrectCore(mykmap, k, watershed, seq, seqLen, nerr,
					corrections, opt.numErrors - nerr + 1, (pingpong++) & 1)) {
				return true;
			}
			/*make changes to the input sequence*/
			if (corrections.size() == 0)
				break;
			for (size_t i = 0; i < corrections.size(); ++i) {
				Correction& corr = corrections[i];
				seq[corr.index] = opt.lowercase ? Kmer::decodeToLower[corr.base] : Kmer::decodeToUpper[corr.base];
			}
		}
		/*stage 3: voting and fix erroneous bases*/
		if ((maxVote = voteErroneousRead(mykmap, k, watershed, seq, seqLen,
				votes)) == 0) {
			/*error free*/
			return true;
		}
		fixErroneousBase(k, seq, seqLen, votes, maxVote, 3, corrections);

		/*make changes to the input sequence*/
		for (size_t i = 0; i < corrections.size(); ++i) {
			Correction& corr = corrections[i];
			seq[corr.index] = opt.lowercase ? Kmer::decodeToLower[corr.base] : Kmer::decodeToUpper[corr.base];
		}
	}

	return false;
}
