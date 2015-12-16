#ifndef _MY_BARRIER_H
#define _MY_BARRIER_H
#include <pthread.h>

class MyBarrier {
public:
	MyBarrier(int nthreads) {
		if (nthreads < 1)
			nthreads = 1;
		pthread_barrier_init(&barrier, NULL, nthreads);
	}
	~MyBarrier() {
		pthread_barrier_destroy(&barrier);
	}
	void wait() {
		pthread_barrier_wait(&barrier);
	}
private:
	pthread_barrier_t barrier;
};
#endif

