#include "kmer.h"

char *
int2bin(uint32_t a, char *buffer, int buf_size) {
	//buffer += (buf_size - 1);

	for (int i = 7; i >= 0; i--) {
		*buffer++ = (a & 1) + '0';
		a >>= 1;
	}

	return buffer;
}
#define NBASE 		0
const char Kmer::encode[128] = { NBASE, NBASE, NBASE, NBASE, NBASE, NBASE,
		NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE,
		NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE,
		NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE,
		NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE,
		NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE,
		NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, 0, NBASE,
		1, NBASE, NBASE, NBASE, 2, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE,
		NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, 3, NBASE, NBASE, NBASE, NBASE,
		NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, 0, NBASE, 1,
		NBASE, NBASE, NBASE, 2, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE, NBASE,
		NBASE, NBASE, NBASE, NBASE, NBASE, 3, NBASE, NBASE, NBASE, NBASE, NBASE,
		NBASE, NBASE, NBASE, NBASE, NBASE, NBASE };
const char Kmer::decodeToLower[4] = { 'a', 'c', 'g', 't' };
const char Kmer::decodeToUpper[4] = { 'A', 'C', 'G', 'T' };

static const uint8_t base_swap[256] = { 0x00, 0x40, 0x80, 0xc0, 0x10, 0x50,
		0x90, 0xd0, 0x20, 0x60, 0xa0, 0xe0, 0x30, 0x70, 0xb0, 0xf0, 0x04, 0x44,
		0x84, 0xc4, 0x14, 0x54, 0x94, 0xd4, 0x24, 0x64, 0xa4, 0xe4, 0x34, 0x74,
		0xb4, 0xf4, 0x08, 0x48, 0x88, 0xc8, 0x18, 0x58, 0x98, 0xd8, 0x28, 0x68,
		0xa8, 0xe8, 0x38, 0x78, 0xb8, 0xf8, 0x0c, 0x4c, 0x8c, 0xcc, 0x1c, 0x5c,
		0x9c, 0xdc, 0x2c, 0x6c, 0xac, 0xec, 0x3c, 0x7c, 0xbc, 0xfc, 0x01, 0x41,
		0x81, 0xc1, 0x11, 0x51, 0x91, 0xd1, 0x21, 0x61, 0xa1, 0xe1, 0x31, 0x71,
		0xb1, 0xf1, 0x05, 0x45, 0x85, 0xc5, 0x15, 0x55, 0x95, 0xd5, 0x25, 0x65,
		0xa5, 0xe5, 0x35, 0x75, 0xb5, 0xf5, 0x09, 0x49, 0x89, 0xc9, 0x19, 0x59,
		0x99, 0xd9, 0x29, 0x69, 0xa9, 0xe9, 0x39, 0x79, 0xb9, 0xf9, 0x0d, 0x4d,
		0x8d, 0xcd, 0x1d, 0x5d, 0x9d, 0xdd, 0x2d, 0x6d, 0xad, 0xed, 0x3d, 0x7d,
		0xbd, 0xfd, 0x02, 0x42, 0x82, 0xc2, 0x12, 0x52, 0x92, 0xd2, 0x22, 0x62,
		0xa2, 0xe2, 0x32, 0x72, 0xb2, 0xf2, 0x06, 0x46, 0x86, 0xc6, 0x16, 0x56,
		0x96, 0xd6, 0x26, 0x66, 0xa6, 0xe6, 0x36, 0x76, 0xb6, 0xf6, 0x0a, 0x4a,
		0x8a, 0xca, 0x1a, 0x5a, 0x9a, 0xda, 0x2a, 0x6a, 0xaa, 0xea, 0x3a, 0x7a,
		0xba, 0xfa, 0x0e, 0x4e, 0x8e, 0xce, 0x1e, 0x5e, 0x9e, 0xde, 0x2e, 0x6e,
		0xae, 0xee, 0x3e, 0x7e, 0xbe, 0xfe, 0x03, 0x43, 0x83, 0xc3, 0x13, 0x53,
		0x93, 0xd3, 0x23, 0x63, 0xa3, 0xe3, 0x33, 0x73, 0xb3, 0xf3, 0x07, 0x47,
		0x87, 0xc7, 0x17, 0x57, 0x97, 0xd7, 0x27, 0x67, 0xa7, 0xe7, 0x37, 0x77,
		0xb7, 0xf7, 0x0b, 0x4b, 0x8b, 0xcb, 0x1b, 0x5b, 0x9b, 0xdb, 0x2b, 0x6b,
		0xab, 0xeb, 0x3b, 0x7b, 0xbb, 0xfb, 0x0f, 0x4f, 0x8f, 0xcf, 0x1f, 0x5f,
		0x9f, 0xdf, 0x2f, 0x6f, 0xaf, 0xef, 0x3f, 0x7f, 0xbf, 0xff };

Kmer::Kmer() {
	k = 0;
	memset(bytes, 0, MAX_K / 4);
}
Kmer::Kmer(int _k) {
	k = _k;
	memset(bytes, 0, MAX_K / 4);
}

Kmer::Kmer(const Kmer& o) {
	k = o.k;
	memcpy(bytes, o.bytes, MAX_K / 4);
}

Kmer::Kmer(int _k, const char *s) {
	k = _k;
	set_kmer(s);
}

Kmer&
Kmer::operator=(const Kmer& o) {
	if (this != &o) {
		k = o.k;
		memcpy(bytes, o.bytes, MAX_K / 4);
	}
	return *this;
}

void Kmer::set_deleted() {
	memset(bytes, 0xff, MAX_K / 4);
}

bool Kmer::operator<(const Kmer& o) const {

	return memcmp(bytes, o.bytes, MAX_K / 4) < 0;
}

bool Kmer::operator==(const Kmer& o) const {

	return memcmp(bytes, o.bytes, MAX_K / 4) == 0;
}

void Kmer::set_kmer(const char *s) {
	int i, j, l;
	//fprintf(stderr,"%s\n",s);
	memset(bytes, 0, MAX_K / 4);

	//fprintf(stderr,"before   "); this->printBinary();

	for (i = 0; i < k; ++i) {
		j = i % 4;
		l = i / 4;
		assert(*s != '\0');
		switch (encode[(uint8_t) *s]) {
		case 0:
			break;
		case 1:
			bytes[l] |= (0x01 << (2 * j));
			break;
		case 2:
			bytes[l] |= (0x02 << (2 * j));
			break;
		case 3:
			bytes[l] |= (0x03 << (2 * j));
			break;
		}
		//fprintf(stderr,"step %2d %c",i,*s);  printBinary();

		s++;
	}
}
Kmer Kmer::set_base(uint8_t base, int index) {
	int j, l;

	Kmer km(*this);
	j = index % 4;
	l = index / 4;

	assert(index >= 0 && index < k);
	assert(base >= 0 && base < 4);

	km.bytes[l] &= ~(3 << (2 * j));
	km.bytes[l] |= (base << (2 * j));

	return km;
}

uint8_t Kmer::get_base(int index) {
	int j, l;
	j = index % 4;
	l = index / 4;
	assert(index >= 0 && index < k);

	uint8_t base = bytes[l] >> (2 * j);
	return base & 3;
}
uint64_t Kmer::hash() const {
	uint64_t ret;
	MurmurHash3_x64_64((const void*) bytes, (k + 3) >> 2, 0, &ret);
	return ret;
}

Kmer Kmer::twin() const {

	Kmer km(*this);
	int k_bytes = (k + 3) >> 2;
	uint32_t k_modmask = (1 << (2 * ((k & 3) ? k & 3 : 4))) - 1;

	for (int i = 0; i < k_bytes; i++) {
		km.bytes[i] = ~bytes[i];
	}

// fprintf(stderr,"~       "); km.printBinary();
	km.bytes[k_bytes - 1] ^= ~k_modmask;
// fprintf(stderr,"!       ");     km.printBinary();

// shift to byte alignment
	km.shiftLeft(8 * k_bytes - 2 * k);

//fprintf(stderr,"shift   ");    km.printBinary();

	uint8_t tmp;
	for (int i = 0; i < k_bytes / 2; ++i) {
		tmp = km.bytes[i];
		km.bytes[i] = base_swap[km.bytes[k_bytes - 1 - i]];
		km.bytes[k_bytes - 1 - i] = base_swap[tmp];
	}

	if ((k_bytes & 1) == 1) {
		km.bytes[k_bytes >> 1] = base_swap[km.bytes[k_bytes >> 1]];
	}

//  fprintf(stderr,"rotated ");    km.printBinary(); fprintf(stderr,"\n");

	return km;
}

Kmer Kmer::forwardBase(const char b) const {
	int s = 2 * ((k + 3) & 3);
	int k_bytes = (k + 3) >> 2;
	uint32_t k_modmask = (1 << (2 * ((k & 3) ? k & 3 : 4))) - 1;

	Kmer km(*this);

	km.shiftRight(2);
//fprintf(stderr,"shift  "); km.printBinary();
	km.bytes[k_bytes - 1] &= k_modmask;

//fprintf(stderr,"mask   "); km.printBinary();
	switch (encode[(uint8_t) b]) {
	case 0:
		km.bytes[k_bytes - 1] |= 0x00 << s;
		break;
	case 1:
		km.bytes[k_bytes - 1] |= 0x01 << s;
		break;
	case 2:
		km.bytes[k_bytes - 1] |= 0x02 << s;
		break;
	case 3:
		km.bytes[k_bytes - 1] |= 0x03 << s;
		break;
	}
//fprintf(stderr,"after  "); km.printBinary();
//km.toString(tmp);
//fprintf(stderr,"%s\n\n",tmp);
	return km;
}

Kmer Kmer::backwardBase(const char b) const {

	int k_bytes = (k + 3) >> 2;
	uint32_t k_modmask = (1 << (2 * ((k & 3) ? k & 3 : 4))) - 1;

	Kmer km(*this);
//fprintf(stderr,"before "); km.printBinary();

	km.shiftLeft(2);
// fprintf(stderr,"shift  "); km.printBinary();
	km.bytes[k_bytes - 1] &= k_modmask;
	if (k % 4 == 0 && k_bytes < MAX_K / 4) {
		km.bytes[k_bytes] = 0x00;
	}
//fprintf(stderr,"mask   "); km.printBinary();

	switch (encode[(uint8_t) b]) {
	case 0:
		km.bytes[0] |= 0x00;
		break;
	case 1:
		km.bytes[0] |= 0x01;
		break;
	case 2:
		km.bytes[0] |= 0x02;
		break;
	case 3:
		km.bytes[0] |= 0x03;
		break;
	}

	return km;
}

unsigned int Kmer::getLastBase() const {
	Kmer km(*this);
	int k_bytes = (k + 3) >> 2;
	uint32_t k_modmask = (1 << (2 * ((k & 3) ? k & 3 : 4))) - 1;
	int s = 2 * ((k + 3) & 3);

	km.bytes[k_bytes - 1] &= k_modmask;
	int base = km.bytes[k_bytes - 1] >> s;
	return (base & 0x03);
}

unsigned int Kmer::getFirstBase() const {
	Kmer km(*this);

	int base = km.bytes[0] & 0x03;
	return base;
}

void Kmer::printBinary() const {
	char buff[9];
	buff[8] = '\0';
	fprintf(stderr, "binary:");
	for (unsigned int i = 0; i < MAX_K / 4; i++) {
		int2bin(bytes[i], buff, 8);
		fprintf(stderr, "%s", buff);
	}
	fprintf(stderr, "\n");
}
string Kmer::toString() const {
	int i, j, l;
	string s;

	for (i = 0; i < k; i++) {
		j = i % 4;
		l = i / 4;
		switch (((bytes[l]) >> (2 * j)) & 0x03) {
		case 0x00:
			s.push_back('A');
			break;
		case 0x01:
			s.push_back('C');
			break;
		case 0x02:
			s.push_back('G');
			break;
		case 0x03:
			s.push_back('T');
			break;
		}
	}
	return s;
}

void Kmer::toString(char * s) const {
	int i, j, l;

	for (i = 0; i < k; i++) {
		j = i % 4;
		l = i / 4;
		switch (((bytes[l]) >> (2 * j)) & 0x03) {
		case 0x00:
			*s = 'A';
			++s;
			break;
		case 0x01:
			*s = 'C';
			++s;
			break;
		case 0x02:
			*s = 'G';
			++s;
			break;
		case 0x03:
			*s = 'T';
			++s;
			break;
		}
	}
	*s = '\0';

	/*
	 fprintf(stderr,"bytes: ");
	 for (i = 0; i < k_bytes; i++ ) {
	 fprintf(stderr,"%x ",bytes[i]);
	 }
	 fprintf(stderr,"\n");
	 */

}

void Kmer::shiftLeft(int shift) {
	if (shift > 0) {
		if (shift < 8) {
			int k_bytes = (k + 3) >> 2;
			for (int i = k_bytes - 1; i > 0; i--) {
				bytes[i] <<= shift;
				bytes[i] |= (uint8_t) (bytes[i - 1] >> (8 - shift));
			}
			bytes[0] <<= shift;
		} else {
			// we should never need this!
			assert(0);
		}
	}
}

void Kmer::shiftRight(int shift) {
	if (shift > 0) {
		if (shift < 8) {
			int k_bytes = (k + 3) >> 2;
			for (int i = 0; i < k_bytes - 1; i++) {
				bytes[i] >>= shift;
				bytes[i] |= (uint8_t) (bytes[i + 1] << (8 - shift));
			}
			bytes[k_bytes - 1] >>= shift;
		} else {
			// bad
			assert(0);
		}
	}
}

