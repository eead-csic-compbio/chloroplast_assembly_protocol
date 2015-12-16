#include "bloom_filter.h"

bloom_filter operator &(const bloom_filter& a, const bloom_filter& b) {
	bloom_filter result = a;
	result &= b;
	return result;
}

bloom_filter operator |(const bloom_filter& a, const bloom_filter& b) {
	bloom_filter result = a;
	result |= b;
	return result;
}

bloom_filter operator ^(const bloom_filter& a, const bloom_filter& b) {
	bloom_filter result = a;
	result ^= b;
	return result;
}
