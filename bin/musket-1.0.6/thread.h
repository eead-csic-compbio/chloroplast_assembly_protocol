#ifndef _MY_THREAD_H
#define _MY_THREAD_H
#include <list>
#include <vector>
#include <iostream>
#include "lock.h"
#include "barrier.h"
#include "message.h"
#include "fastq.h"
#include "kmer.h"
#include "bloom_filter.h"
#include "hash.h"
#include "kmerintpair.h"

#define DEFAULT_MAX_MULTI			512
class MyThread {
public:
	MyThread(const vector<pair<int, unsigned int> >& predicted_nkmers
			, size_t bits_per_element, size_t randseed) :
			barrier(2) {
		/*set the thread running status*/
		threadRunning = false;

		/*create the Bloom filters and hash tables for each k-mer size*/
		size_t numKmerSizes = predicted_nkmers.size();

		numRecvKmers.resize(numKmerSizes);
		totalRecvCov.resize(numKmerSizes);
		numRecvDels.resize(numKmerSizes);
		maxMults.resize(numKmerSizes);

		BFs.resize(numKmerSizes);
		kmaps.resize(numKmerSizes);
		for (size_t i = 0; i < numKmerSizes; ++i) {
			BFs[i] = new bloom_filter(predicted_nkmers[i].second,
					bits_per_element, randseed);
			kmaps[i] = new hmap_t;
		}

		//create a thread
		create();

		/*create buffer for multiplicity information*/
		kmerFreqs.resize(numKmerSizes);
		for (size_t i = 0; i < numKmerSizes; ++i) {
			kmerFreqs[i].resize(DEFAULT_MAX_MULTI, 0);
		}
	}
	~MyThread() {
		//stop the thread
		stop();

		//destroy hash table
		for (size_t i = 0; i < kmaps.size(); ++i) {
			if (kmaps[i]) {
				delete kmaps[i];
			}
		}
		kmaps.clear();

		//destroy bloom filter
		for (size_t i = 0; i < BFs.size(); ++i) {
			if (BFs[i]) {
				delete BFs[i];
			}
		}
		BFs.clear();

		//destroy node counts
		for (size_t i = 0; i < kmerFreqs.size(); ++i) {
			kmerFreqs[i].clear();
		}
		kmerFreqs.clear();
	}
	//statistic information
	inline size_t getItsNumKmers(int index) {
		return numRecvKmers[index];
	}
	inline size_t getItsTotalCov(int index) {
		return totalRecvCov[index];
	}
	inline size_t getItsNumDels(int index) {
		return numRecvDels[index];
	}
	inline size_t getItsKmapSize(int index) {
		return kmaps[index]->size();
	}
	inline size_t*
	getKmerFreqs(int index) {
		return &kmerFreqs[index][0];
	}
	inline hmap_t*
	getKmap(int index) {
		return kmaps[index];
	}
	inline size_t getMaxMulti(int index) {
		return maxMults[index];
	}

	//thread
	inline void create() {
		//create a thread
		threadRunning = true;
		if (pthread_create(&tid, NULL, threadfunc, this) != 0) {
			cerr << "Error: failed to create thread worker" << endl;
		}
	}
	inline void stop() {
		if (threadRunning) {
			//destroy the thread
			sendMsg(new MyMessageCommand(MSG_COMMAND_THREAD_EXIT));
			//wait for the barrier
			barrierwait();
			//set the flag
			threadRunning = false;
		}
	}
	//BF
	inline void releaseBF(int index) {
		if (BFs[index]) {
			delete BFs[index];
			BFs[index] = NULL;
		}
	}
	inline void releaseBF() {
		for (size_t i = 0; i < BFs.size(); ++i) {
			releaseBF(i);
		}
	}

	//messages
	size_t
	getNumMsg();
	void
	sendMsg(MyMessage* msg);
	MyMessage*
	recvMsg();

	//jobs
	void
	allfillup();
	void
	filtration();
	void
	excludeunique();

	//barrier
	inline void barrierwait() {
		barrier.wait();
	}
	inline pthread_t getTid() {
		return tid;
	}
private:
	//thread id
	pthread_t tid;
	bool threadRunning;

	//message list
	list<MyMessage*> msgs;

	//lock for this thread
	MyLock mutex;

	//barrier for this thread 
	MyBarrier barrier;

	//k-mer hash table
	vector<hmap_t*> kmaps;

	//currently used bloom filter
	vector<bloom_filter*> BFs;

	//variables
	vector<size_t> numRecvKmers;
	vector<size_t> totalRecvCov;
	vector<size_t> numRecvDels;
	vector<size_t> maxMults;
	vector<vector<size_t> > kmerFreqs;

	//thread function
	static void*
	threadfunc(void*);
	static const size_t maxCapacity = 32;
};
#endif

