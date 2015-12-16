#ifndef _MY_LOCK_H
#define _MY_LOCK_H
#include <pthread.h>

class MyLock {
public:
	MyLock(bool locking = true) {
		_locking = locking;
		pthread_mutex_init(&mutex, NULL);
	}
	~MyLock() {
		pthread_mutex_destroy(&mutex);
	}
	inline void enable(bool v) {
		_locking = v;
	}
	inline void lock() {
		if (_locking)
			pthread_mutex_lock(&mutex);
	}
	inline void unlock() {
		if (_locking)
			pthread_mutex_unlock(&mutex);
	}
private:
	bool _locking;
	pthread_mutex_t mutex;
};
#endif

