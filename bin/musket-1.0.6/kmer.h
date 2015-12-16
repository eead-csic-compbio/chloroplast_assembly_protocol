#ifndef INCLUDE_KMER_HPP
#define INCLUDE_KMER_HPP

#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 28
#endif

#include <stdio.h>
#include <stdint.h>
#include "hash.h"
#include <cassert>
#include <cstring>
#include <string>
using namespace std;

struct KmerConfig {
	KmerConfig(int _k) {
		set(_k);
	}
	inline void set(int _k) {
		assert(_k <= MAX_KMER_SIZE);
		assert(_k > 0);

		// we can only call this once
		k = _k;
		k_bytes = (k + 3) / 4;
		//  k_longs = (_k+15)/16;
		k_modmask = (1 << (2 * ((k % 4) ? k % 4 : 4))) - 1;
	}
	int k;
	int k_bytes;
	unsigned int k_modmask;
};

char *
int2bin(uint32_t a, char *buffer, int buf_size);

class Kmer {
public:
	Kmer();
	Kmer(int _k);
	Kmer(const Kmer& o);
	Kmer(int _k, const char *s);

	Kmer&
	operator=(const Kmer& o);

	void
	set_deleted();

	bool
	operator<(const Kmer& o) const;

	bool
	operator==(const Kmer& o) const;

	bool operator!=(const Kmer& o) const {
		return !(*this == o);
	}

	void
	set_kmer(const char *s);
	Kmer
	set_base(uint8_t base, int index);
	uint8_t
	get_base(int index);

	uint64_t
	hash() const;

	Kmer
	twin() const;

	Kmer
	forwardBase(const char b) const;

	Kmer
	backwardBase(const char b) const;

	unsigned int
	getLastBase() const;
	unsigned int
	getFirstBase() const;

	void
	printBinary() const;

	void
	toString(char * s) const;
	string toString() const;

	// static functions
	static void
	set_k(int _k);

	static const int MAX_K = MAX_KMER_SIZE;
	static const char encode[128];
	static const char decodeToLower[4];
	static const char decodeToUpper[4];

	uint8_t k;
private:

	// data fields
	uint8_t bytes[MAX_K / 4];

	// private functions
	void
	shiftRight(int shift);

	void
	shiftLeft(int shift);
}__attribute__((packed));

struct KmerHash {
	size_t operator()(const Kmer &km) const {
		return km.hash();
	}
};

#endif //KMER_HPP
