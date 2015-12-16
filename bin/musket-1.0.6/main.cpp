#include "option.h"
#include "voting.h"
#include "paracorrect.h"

void findWatershed(ProgramOptions& opt, vector<size_t>& kmerFreqs, int ikmer) {

	size_t valleyIndex, peakIndex;
	size_t minMulti, maxMulti;
	const size_t lowerBoundMulti = 100;

	/*reaching the plateau*/
	valleyIndex = 2;
	for (size_t i = valleyIndex + 1; i < kmerFreqs.size(); ++i) {
		if (kmerFreqs[i] > kmerFreqs[valleyIndex]) {
			break;
		}
		valleyIndex++;
	}
	valleyIndex++;

	/*find the peak*/
	peakIndex = valleyIndex;
	maxMulti = kmerFreqs[peakIndex];
	for (size_t i = peakIndex + 1; i < kmerFreqs.size(); ++i) {
		if (kmerFreqs[i] > maxMulti) {
			maxMulti = kmerFreqs[i];
			peakIndex = i;
		}
	}
	/*find the valley*/
	minMulti = kmerFreqs[valleyIndex];
	for (size_t i = valleyIndex + 1; i < peakIndex; ++i) {
		if (kmerFreqs[i] < lowerBoundMulti) {
			continue;
		}
		if (kmerFreqs[i] < minMulti) {
			minMulti = kmerFreqs[i];
			valleyIndex = i;
		}
	}

#if 0
	/*adjust the valley index by calculting the number of k-mers on both sides*/
	size_t count = 0, count2 = 0;
	for (size_t i = 0; i < valleyIndex; ++i) {
		count += kmerFreqs[i];
	}
	for (size_t i = valleyIndex; i < kmerFreqs.size(); ++i) {
		count2 += kmerFreqs[i];
	}
	/*shift the valley index to the origin*/
	while (count > count2) {
		if (valleyIndex <= 2)
			break;
		--valleyIndex;
		count -= kmerFreqs[valleyIndex];
		count2 += kmerFreqs[valleyIndex];
	}
#endif
	/*save the minimum multiplicity*/
	opt.watershed[ikmer] = valleyIndex;

#if 0
	/*Calculate the mean and the standard devision of the coverage*/
	double avgCov = 0;
	size_t numKmers = 0;

	for(size_t cov = valleyIndex; cov < kmerFreqs.size(); ++cov)
	{
		/*check the multiplictiy of the k-mer*/
		if(kmerFreqs[cov] < lowerBoundMulti) {
			continue;
		}
		numKmers += kmerFreqs[cov];
		avgCov += kmerFreqs[cov] * cov;
	}
	avgCov /= numKmers;

	double stdCov = 0;
	for(size_t cov = valleyIndex; cov < kmerFreqs.size(); ++cov)
	{
		if(kmerFreqs[cov] < lowerBoundMulti) {
			continue;
		}
		stdCov += kmerFreqs[cov] * (cov - avgCov) * (cov - avgCov);
	}
	stdCov /= numKmers;
	stdCov = sqrt(stdCov);
	cerr << "Average coverage: " << avgCov << " +/- " << stdCov << endl;
	cerr << "Peak index: " << peakIndex << endl;
#endif
}

//global variables
MyLock* globalMutexInput, *globalMutexOutput;
MyBarrier* globalBarrier;
ofstream correctedReadsOutput;

//using master+worker model
void ParaKmerEC(ProgramOptions &opt) {

	uint64_t n_read = 0;
	//open the input files
	FastqFile FQ(opt.files);

	//create global barriers and mutex
	globalMutexInput = new MyLock();
	globalMutexOutput = new MyLock();
	globalBarrier = new MyBarrier(1 + opt.nworkers); //1 master + #workers

	//create worker threads
	vector<MyThread*> workers(opt.nworkers, NULL);
	for (size_t i = 0; i < opt.nworkers; i++) {
		workers[i] = new MyThread(opt.klist, 1, 7919); /*using fixed random seed*/
	}

	//create message buffer for each thread of each master
	const size_t maxKmerPerSlot = opt.maxKmerPerSlot;
	vector<vector<vector<pair<int, Kmer> >*> > msgSlots;
	msgSlots.resize(opt.nmasters);
	for (size_t i = 0; i < msgSlots.size(); i++) {
		msgSlots[i].resize(opt.nworkers);
		//allocate space for the workers of each master
		for (size_t j = 0; j < msgSlots[i].size(); j++) {
			msgSlots[i][j] = new vector<pair<int, Kmer> >;
			msgSlots[i][j]->reserve(maxKmerPerSlot);
		}
	}

	//create read batches for masters
	const size_t maxReadsPerBatch = opt.maxReadsPerBatch;
	vector<vector<string> > batches;
	batches.resize(opt.nmasters);
	for (size_t i = 0; i < opt.nmasters; i++) {
		batches[i].reserve(maxReadsPerBatch);
	}

	//set the number of masters
	omp_set_num_threads(opt.nmasters);

	//enable or disable global locking
	globalMutexInput->enable(opt.nmasters > 1);

	size_t progress = 1000 * maxReadsPerBatch;
	double stime, etime;
	bool done;

	/*********************STAGE 1*******************************\
	 *Determine all non-unique k-mers
	 ***********************************************************/
	stime = getSysTime();
	cerr << "Determine all non-unique k-mers" << endl;
	/*random number generator*/
	srand48(19);

	//send STARTUP message
	for (size_t i = 0; i < opt.nworkers; i++) {
		workers[i]->sendMsg(new MyMessageCommand(MSG_COMMAND_ALLFILLUP));
	}

	/*re-open all input files*/
	FQ.reopen();

	/*read the sequences*/
	done = false;
	n_read = 0;
#pragma omp parallel default(shared)
	for (;;) {
		int master = omp_get_thread_num();
		//read a batch of reads
		globalMutexInput->lock();
		if (done) {
			globalMutexInput->unlock();
			break;
		}
		//read a batch of reads
		FQ.read_next_batch(batches[master], maxReadsPerBatch, done, NULL);
		n_read += batches[master].size();
		if (n_read % progress == 0) {
			cerr << "processed " << n_read << " reads" << endl;
		}
		globalMutexInput->unlock();

		/*for each k-mer size in the k-mer size list*/
		for (size_t ikmer = 0; ikmer < opt.klist.size(); ++ikmer) {
			/*set the global k-mer size*/
			size_t k = opt.klist[ikmer].first;

			/*process all k-mers in this read batch usign the specified k-mer size*/
			for (size_t readIdx = 0; readIdx < batches[master].size();
					readIdx++) {
				char* seq = (char*) batches[master][readIdx].data();
				size_t seqLen = batches[master][readIdx].length();

				/*check the sequence length*/
				if (seqLen < k || seqLen > MAX_SEQ_LENGTH) {
					continue;
				}

				/*insert the k-mers into the Bloom filter*/
				Kmer km(k, seq);
				for (size_t i = 0; i <= seqLen - k; ++i) {
					if (i > 0) {
						km = km.forwardBase(seq[i + k - 1]);
					}
					Kmer rep = km.twin();
					if (km < rep) {
						rep = km;
					}

					//compute the destination worker
					int worker = rep.hash() % opt.nworkers;
					//add the canonical k-mer into the message buffer
					msgSlots[master][worker]->push_back(make_pair(ikmer, rep));
					if (msgSlots[master][worker]->size() >= maxKmerPerSlot) {
						//send a kmer message to the worker
						workers[worker]->sendMsg(
								new MyMessageKmers(msgSlots[master][worker]));
						msgSlots[master][worker] = new vector<pair<int, Kmer> >;
						msgSlots[master][worker]->reserve(maxKmerPerSlot);
					}
				}
			}
		}
	}

	//send all slots
	for (size_t master = 0; master < opt.nmasters; master++) {
		for (size_t worker = 0; worker < opt.nworkers; worker++) {
			//send a kmer message to the worker
			if (msgSlots[master][worker]->size() > 0) {
				workers[worker]->sendMsg(
						new MyMessageKmers(msgSlots[master][worker]));
				msgSlots[master][worker] = new vector<pair<int, Kmer> >;
				msgSlots[master][worker]->reserve(maxKmerPerSlot);
			}
		}
	}
	//send DONE message
	for (size_t worker = 0; worker < opt.nworkers; worker++) {
		workers[worker]->sendMsg(
				new MyMessageCommand(MSG_COMMAND_ALLFILLUP_DONE));
	}
	globalBarrier->wait(); //global synchronization
	etime = getSysTime();
	cerr << "Time taken: " << etime - stime << " seconds" << endl;

	/************************STAGE 2**********************
	 Calculate the multiplicity of each k-mer
	 * ***************************************************/
	stime = getSysTime();
	cerr << "Calculate the multiplicity of each k-mer" << endl;

	/*random number generator*/
	srand48(19);

	//send STARTUP message
	for (size_t i = 0; i < opt.nworkers; i++) {
		workers[i]->sendMsg(new MyMessageCommand(MSG_COMMAND_FILTRATION));
	}

	/*re-open all input files*/
	FQ.reopen();

	/*read the sequences*/
	done = false;
	n_read = 0;
#pragma omp parallel default(shared)
	for (;;) {
		int master = omp_get_thread_num();
		//read a sequence
		globalMutexInput->lock();
		if (done) {
			globalMutexInput->unlock();
			break;
		}
		//read a batch of reads
		FQ.read_next_batch(batches[master], maxReadsPerBatch, done, NULL);
		n_read += batches[master].size();
		if (n_read % progress == 0) {
			cerr << "processed " << n_read << " reads" << endl;
		}
		globalMutexInput->unlock();

		/*for each k-mer size in the k-mer size list*/
		for (size_t ikmer = 0; ikmer < opt.klist.size(); ++ikmer) {
			/*set the global k-mer size*/
			size_t k = opt.klist[ikmer].first;

			/*process all k-mers in this read batch using the specified k-mer size*/
			for (size_t readIdx = 0; readIdx < batches[master].size();
					readIdx++) {
				char* seq = (char*) batches[master][readIdx].data();
				size_t seqLen = batches[master][readIdx].length();

				/*check the sequence length*/
				if (seqLen < k || seqLen > MAX_SEQ_LENGTH) {
					continue;
				}
				Kmer km(k, seq);
				for (size_t i = 0; i <= seqLen - k; ++i) {
					if (i > 0) {
						km = km.forwardBase(seq[i + k - 1]);
					}
					Kmer rep = km.twin();
					if (km < rep) {
						rep = km;
					}

					//compute the destination worker
					int worker = rep.hash() % opt.nworkers;
					//add the canonical k-mer into the message buffer
					msgSlots[master][worker]->push_back(make_pair(ikmer, rep));
					if (msgSlots[master][worker]->size() >= maxKmerPerSlot) {
						//send a kmer message to the worker
						workers[worker]->sendMsg(
								new MyMessageKmers(msgSlots[master][worker]));
						msgSlots[master][worker] = new vector<pair<int, Kmer> >;
						msgSlots[master][worker]->reserve(maxKmerPerSlot);
					}
				}
			}
		}
	}

	//send all slots
	for (size_t master = 0; master < opt.nmasters; master++) {
		for (size_t worker = 0; worker < opt.nworkers; worker++) {
			if (msgSlots[master][worker]->size() > 0) {
				//send a kmer message to the worker
				workers[worker]->sendMsg(
						new MyMessageKmers(msgSlots[master][worker]));
				msgSlots[master][worker] = new vector<pair<int, Kmer> >;
				msgSlots[master][worker]->reserve(maxKmerPerSlot);
			}
		}
	}
	//send DONE message
	for (size_t worker = 0; worker < opt.nworkers; worker++) {
		workers[worker]->sendMsg(
				new MyMessageCommand(MSG_COMMAND_FILTRATION_DONE));
	}
	globalBarrier->wait(); //global synchronization

	etime = getSysTime();
	cerr << "Time taken: " << etime - stime << " seconds" << endl;

	/******************STAGE 3**********************
	 * Remove unique k-mers in the hash tables
	 * *********************************************/
	stime = getSysTime();
	cerr << "Remove unique k-mers in the hash tables" << endl;

	//send STARTUP message
	for (size_t i = 0; i < opt.nworkers; i++) {
		workers[i]->sendMsg(new MyMessageCommand(MSG_COMMAND_EXCLUDEUNIQUE));
	}
	//send DONE message
	for (size_t worker = 0; worker < opt.nworkers; worker++) {
		workers[worker]->sendMsg(
				new MyMessageCommand(MSG_COMMAND_EXCLUDEUNIQUE_DONE));
	}
	globalBarrier->wait(); //global synchronization

	//stop the execution of worker threads and release the bloom filters
	for (size_t worker = 0; worker < opt.nworkers; worker++) {
		//stop the worker thread
		workers[worker]->stop();

		//release its bloom filter structures
		workers[worker]->releaseBF();
	}
	//release message slots
	for (size_t master = 0; master < opt.nmasters; master++) {
		for (size_t worker = 0; worker < opt.nworkers; worker++) {
			msgSlots[master][worker]->clear();
			delete msgSlots[master][worker];
		}
		msgSlots[master].clear();
	}
	msgSlots.clear();

	//release the batches
	for (size_t i = 0; i < opt.nmasters; i++) {
		batches[i].clear();
	}
	batches.clear();

	etime = getSysTime();
	cerr << "Time taken: " << etime - stime << " seconds" << endl;

	/******************STAGE 4**********************
	 * Estimate the minimum multiplicity threshold
	 * **********************************************/
	stime = getSysTime();
	cerr << "Estimate the minimum multiplicity threshold" << endl;

	//gather node counts and establish cutoff
	vector<vector<size_t> > kmerFreqs;
	kmerFreqs.resize(opt.klist.size());
	for (size_t i = 0; i < kmerFreqs.size(); ++i) {
		kmerFreqs[i].resize(DEFAULT_MAX_MULTI, 0);
	}

	/*for each k-mer size*/
	for (size_t ikmer = 0; ikmer < opt.klist.size(); ++ikmer) {
		uint64_t num_kmers = 0;
		uint64_t kmap_size = 0;
		uint64_t total_cov = 0;
		uint64_t n_del = 0;
		unsigned int maxMulti = 0;

		cerr << "For k-mer size " << opt.klist[ikmer].first << endl;
		/*get the k-mer information for each k-mer size*/
		for (size_t worker = 0; worker < opt.nworkers; worker++) {
			total_cov += workers[worker]->getItsTotalCov(ikmer);
			num_kmers += workers[worker]->getItsNumKmers(ikmer);
			kmap_size += workers[worker]->getItsKmapSize(ikmer);
			n_del += workers[worker]->getItsNumDels(ikmer);

			/*get the maximal multiplicity*/
			if (workers[worker]->getMaxMulti(ikmer) > maxMulti) {
				maxMulti = workers[worker]->getMaxMulti(ikmer);
			}

			/*get the multiplicity information to build histogram from multiplicity 2*/
			size_t* dirKmerFreqs = &kmerFreqs[ikmer][0];
			size_t* srcKmerFreqs = workers[worker]->getKmerFreqs(ikmer);
			for (size_t j = 2; j < kmerFreqs[ikmer].size() && j < maxMulti;
					j++) {
				dirKmerFreqs[j] += srcKmerFreqs[j];
			}
		}
		/*resize the counts buffer*/
		if (maxMulti < kmerFreqs[ikmer].size()) {
			kmerFreqs[ikmer].resize(maxMulti);
		}
		cerr << "Maximum multiplicity " << maxMulti << endl;

		/*output the k-mer information for each k-mer size*/
		total_cov -= n_del;
		/*cerr << "\tprocessed " << num_kmers << " k-mers in " << n_read << " reads"
		 << endl;
		 cerr << "\tfound " << kmap_size << " non-filtered k-mers, removed "
		 << n_del << "(depending on BF)" << endl;
		 uint64_t filtered_kmers = num_kmers - total_cov;

		 cerr << "\ttotal coverage " << total_cov
		 << ", estimated number of k-mers " << filtered_kmers << endl;
		 cerr << "\taverage coverage " << (total_cov / ((double) kmap_size))
		 << endl;
		 cerr << "\tnum kmers " << num_kmers << " filtered kmers: "<< filtered_kmers << "kmers in hash " << kmap_size
		 << endl;*/
		cerr << "\tCoverage: " << total_cov << endl;
		cerr << "\tAverage coverage: " << total_cov / ((double) kmap_size)
				<< endl;
		cerr << "\tNumber of k-mers " << kmap_size << endl;

		cerr << "\tK-mer frequency histogram (>=2): " << endl;
		size_t* dirKmerFreqs = &kmerFreqs[ikmer][0];
		for (size_t i = 2; i < kmerFreqs[ikmer].size(); ++i) {
			cerr << dirKmerFreqs[i] << ",";
		}
		cerr << endl;

		/*estimate the minimal multiplicity to differentiate the trusted and un-trusted k-mers*/
		findWatershed(opt, kmerFreqs[ikmer], ikmer);
		cerr << "Min multi:  " << opt.watershed[ikmer] << endl;
	}
	etime = getSysTime();
	cerr << "Time taken: " << etime - stime << " seconds" << endl;

	/******************Stage 5 ERROR CORRECTION********************
	 * match file name to output options or input file name...	TODO
	 * check whether to make it fasta/fastq	TODO
	 * ***********************************************************/
	cerr << "Perform error correction" << endl;

	/*gather hash-maps from all worker threads*/
	vector<vector<hmap_t*> > kmaps;
	kmaps.resize(opt.klist.size());
	for (size_t i = 0; i < kmaps.size(); i++) {
		kmaps[i].resize(opt.nworkers);
	}
	for (size_t ikmer = 0; ikmer < kmaps.size(); ikmer++) {
		for (size_t worker = 0; worker < opt.nworkers; worker++) {
			kmaps[ikmer][worker] = workers[worker]->getKmap(ikmer);
		}
	}

	/*create read batches for all thread workers*/
	vector<vector<string> > batchNames, batchScores;
	batches.resize(opt.nthreads);
	batchNames.resize(opt.nthreads);
	batchScores.resize(opt.nthreads);
	for (size_t i = 0; i < opt.nthreads; i++) {
		batches[i].reserve(maxReadsPerBatch);
		batchNames[i].reserve(maxReadsPerBatch);
		batchScores[i].reserve(maxReadsPerBatch);
	}

	/*prepare for error correction*/
	done = false;
	n_read = 0;
	globalMutexInput->enable(true); //always enable the mutex regardless of the number of threads
	globalMutexOutput->enable(true);

	/*reset the number of OPENMP threads.
	 * Here, all specified threads are used for error correction
	 */
	omp_set_num_threads(opt.nthreads);

	/*re-create barrier*/
	delete globalBarrier;
	globalBarrier = new MyBarrier(opt.nthreads);

	/*check the output mode*/
	const bool singleOutput = opt.outfiles.size() == 0;

	/*input and output file iterators*/
	vector<string>::const_iterator inFileIter = opt.files.begin();
	vector<string>::const_iterator outFileIter = opt.outfiles.begin();

	if(singleOutput){
		/*open the output file of corrected reads*/
		correctedReadsOutput.open(opt.correctedFileName.c_str(),
			ios::binary | ios::out | ios::trunc);

		/*reopen all input files to conduct error correction*/
		FQ.reopen();
	}else{
		/*open the output file of corrected reads*/
		correctedReadsOutput.open(outFileIter->c_str(),
      ios::binary | ios::out | ios::trunc);

		/*open the input file*/
		FQ.reset(inFileIter->c_str());

		/*increase the iterators*/
		++inFileIter;
		++outFileIter;
	}

	/*start the main parallel loop*/
	if (opt.keepOrder) { /*if keeping the sequence order*/
		while (!done) {
			//read a batch of reads
			FQ.read_next_batch_with_name_and_qual(batches[0], batchNames[0],
					batchScores[0], maxReadsPerBatch, done, NULL);

			//record the number of processed reads
			n_read += batches[0].size();
			if (n_read % progress == 0) {
				cerr << "processed " << n_read << " reads" << endl;
			}

			/*correct this read batch*/
			size_t readIdx;
#pragma omp parallel for default(shared) private(readIdx) schedule(static, 1)
			for (readIdx = 0; readIdx < batches[0].size(); readIdx++) {
				char* seq = (char*) batches[0][readIdx].data();
				size_t seqLen = batches[0][readIdx].length();

				if (seqLen > MAX_SEQ_LENGTH) {
					continue;
				}
				/*perform ungapped error correction*/
				uint8_t votes[MAX_SEQ_LENGTH][4];
				bool errorFree = false;
				size_t numKmers = opt.klist.size();
				for (size_t ikmer = 0; ikmer < numKmers; ++ikmer) {
					/*set the global k-mer size*/
					int k = opt.klist[ikmer].first;

					/*invoke the core function of error correction*/
					if ((errorFree = ParaCorrect(opt, kmaps[ikmer], k,
							opt.watershed[ikmer], seq, seqLen, votes))
							== true) {
						break;
					}
				}
				if (!errorFree && opt.maxTrim > 0) {
					int ikmer = 0;
					pair<int, int> longestRegion = make_pair(0, -1);
					/*use the longest k-mer to trim the sequence*/
					int k = opt.klist[ikmer].first;
					/*attempt to trim the sequence using the largest k-mer size*/
					if (isTrimmable(kmaps[ikmer], k, opt.watershed[ikmer], seq,
							seqLen, opt.maxTrim, longestRegion)) {
						seqLen = longestRegion.second - longestRegion.first;
						if (longestRegion.first > 0) {
							memmove(seq, seq + longestRegion.first, seqLen);
						}
						seq[seqLen] = '\0';
					}
				}
			}/*for readIdx parallel section*/

			/*Having finished the read batch, output corrected reads*/
			for (size_t readIdx = 0; readIdx < batches[0].size(); readIdx++) {
				/*check if the quality scores exist*/
				if (batchScores[0][readIdx].size() == 0) {
					/*FASTA format*/
					correctedReadsOutput << ">" << batchNames[0][readIdx]
							<< '\n' << batches[0][readIdx] << '\n';
				} else {
					/*FASTQ format*/
					correctedReadsOutput << "@" << batchNames[0][readIdx]
							<< '\n' << batches[0][readIdx] << "\n+\n"
							<< batchScores[0][readIdx] << '\n';
				}
			}

			/*check the output mode*/
			if(done && !singleOutput){
				/*close the output file*/
				correctedReadsOutput.close();

				/*close the input file*/
				FQ.close();

				/*open the next file*/
				if(inFileIter != opt.files.end()){
				  /*open the output file of corrected reads*/
    			correctedReadsOutput.open(outFileIter->c_str(),
      			ios::binary | ios::out | ios::trunc);
				
    			/*open the input file*/
    			FQ.reset(inFileIter->c_str());

    			/*increase the iterators*/
				  ++inFileIter;
				  ++outFileIter;

					/*reset the value of done*/
					done = false;
				}
			}
		}
	} else {
#pragma omp parallel default(shared)
		for (;;) {
			const int worker = omp_get_thread_num(); //get the ID of the worker

			/*read a batch of reads*/
			globalMutexInput->lock();
			if(singleOutput){
				/*release the lock if DONE*/
				if (done) {
					globalMutexInput->unlock();
					break;
				}
			}else if(done){
				/*release the lock*/
				globalMutexInput->unlock();
				
				/*synchronize all threads*/
				globalBarrier->wait();

				/*re-open the file*/
				if(worker == 0){
       		/*close the output file*/
       		correctedReadsOutput.close();

       		/*close the input file*/
       		FQ.close();

 	      	/*open the next file*/
  	     	if(inFileIter != opt.files.end()){
 	        	/*open the output file of corrected reads*/
  	       	correctedReadsOutput.open(outFileIter->c_str(),
    	       	ios::binary | ios::out | ios::trunc);

      	   	/*open the input file*/
        	 	FQ.reset(inFileIter->c_str());

         		/*increase the iterators*/
         		++inFileIter;
         		++outFileIter;

         		/*reset the value of done*/
         		done = false;
					}
				}
				globalBarrier->wait();

				/*chech the status of global variable*/
				if(done){
					break;
				}

				/*re-enter the critical section*/
				globalMutexInput->lock();
			}
			//read a batch of reads
			FQ.read_next_batch_with_name_and_qual(batches[worker],
					batchNames[worker], batchScores[worker], maxReadsPerBatch,
					done, NULL);

			//record the number of processed reads
			n_read += batches[worker].size();
			if (n_read % progress == 0) {
				cerr << "processed " << n_read << " reads" << endl;
			}
			globalMutexInput->unlock();

			/*correct this read batch*/
			uint8_t votes[MAX_SEQ_LENGTH][4];
			for (size_t readIdx = 0; readIdx < batches[worker].size();
					readIdx++) {
				char* seq = (char*) batches[worker][readIdx].data();
				size_t seqLen = batches[worker][readIdx].length();

				if (seqLen > MAX_SEQ_LENGTH) {
					continue;
				}
				/*perform ungapped error correction*/
				bool errorFree = false;
				size_t numKmers = opt.klist.size();
				for (size_t ikmer = 0; ikmer < numKmers; ++ikmer) {
					/*set the global k-mer size*/
					int k = opt.klist[ikmer].first;

					/*invoke the core function of error correction*/
					if ((errorFree = ParaCorrect(opt, kmaps[ikmer], k,
							opt.watershed[ikmer], seq, seqLen, votes))
							== true) {
						break;
					}
				}
				if (!errorFree && opt.maxTrim > 0) {
					int ikmer = 0;
					pair<int, int> longestRegion = make_pair(0, -1);
					/*use the longest k-mer to trim the sequence*/
					int k = opt.klist[ikmer].first;
					/*attempt to trim the sequence using the largest k-mer size*/
					if (isTrimmable(kmaps[ikmer], k, opt.watershed[ikmer], seq,
							seqLen, opt.maxTrim, longestRegion)) {
						seqLen = longestRegion.second - longestRegion.first;
						if (longestRegion.first > 0) {
							memmove(seq, seq + longestRegion.first, seqLen);
						}
						seq[seqLen] = '\0';
					}
				}
			}/*for readIdx*/

			//Having finished the read batch, output corrected reads
			globalMutexOutput->lock();
			for (size_t readIdx = 0; readIdx < batches[worker].size();
					readIdx++) {
				/*check if the quality scores exist*/
				if (batchScores[worker][readIdx].size() == 0) {
					/*FASTA format*/
					correctedReadsOutput << ">" << batchNames[worker][readIdx]
							<< '\n' << batches[worker][readIdx] << '\n';
				} else {
					/*FASTQ format*/
					correctedReadsOutput << "@" << batchNames[worker][readIdx]
							<< '\n' << batches[worker][readIdx] << "\n+\n"
							<< batchScores[worker][readIdx] << '\n';
				}
			}
			globalMutexOutput->unlock();
		}
	}
	if(singleOutput){
		/*close the output file for corrected reads*/
		correctedReadsOutput.close();

		/*close all input files*/
		FQ.close();
	}
	cerr << "Processed " << n_read << " reads" << endl;

	/*release the batches*/
	for (size_t i = 0; i < opt.nthreads; i++) {
		vector<string>& batch = batches[i];
		vector<string>& names = batchNames[i];
		vector<string>& qualityScores = batchScores[i];
  	//clear the batch
  	for(vector<string>::iterator iter = batch.begin(); iter != batch.end(); ++iter){
    	iter->clear();
  	}
  	batch.clear();
  	for(vector<string>::iterator iter = names.begin(); iter != names.end(); ++iter){
    	iter->clear();
  	}
  	names.clear();
	  for(vector<string>::iterator iter = qualityScores.begin(); iter != qualityScores.end(); ++iter){
    	iter->clear();
  	}
  	qualityScores.clear();
	}
	batches.clear();
	batchNames.clear();
	batchScores.clear();

	/******************STAGE 6********************/
	cerr << "Finish the error correction" << endl;

//release all workers
	for (size_t worker = 0; worker < opt.nworkers; worker++) {
		delete workers[worker];
	}
	workers.clear();

//release resources
	delete globalMutexInput;
	delete globalMutexOutput;
	delete globalBarrier;
}

int main(int argc, char **argv) {

	double stime = getSysTime();

//parse command line options
	ProgramOptions opt;
	ParseOptions(argc, argv, opt);

	/*start the kernel of error correction*/
	ParaKmerEC(opt);

//debugKmer(opt);
	double etime = getSysTime();
	cerr << "Overall execution time " << etime - stime << " seconds" << endl;
}
