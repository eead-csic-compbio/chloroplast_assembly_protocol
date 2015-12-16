#ifndef _MY_KMERINTPAIR_H
#define _MY_KMERINTPAIR_H

#include "kmer.h"
#include "sparsehash/sparse_hash_map"

using namespace std;
using namespace google;

struct KmerIntPair {
	KmerIntPair() {
	}
	;
	KmerIntPair(const Kmer &km, unsigned int k);

	char v[sizeof(Kmer) + sizeof(char)];
	unsigned int
	GetVal() const;
	void
	SetVal(const unsigned int k);
	const Kmer&
	GetKey() const;
	void
	SetKey(const Kmer& km);

	static const size_t KmerOffset;
	static const size_t IntOffset;
};

struct SelectKmerKey {
	const Kmer&
	operator()(const KmerIntPair &p) const {
		return p.GetKey();
	}
	typedef Kmer result_type;
};

struct SetKmerKey {
	void operator()(KmerIntPair *value, const Kmer& km) const {
		memcpy(value + KmerIntPair::KmerOffset, &km, sizeof(Kmer));
	}
};

typedef sparse_hashtable<KmerIntPair, Kmer, KmerHash, SelectKmerKey, SetKmerKey,
		std::equal_to<Kmer>, std::allocator<KmerIntPair> > hmap_t;

#endif

