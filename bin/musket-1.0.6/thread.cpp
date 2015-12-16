#include "thread.h"

extern MyBarrier* globalBarrier;

void*
MyThread::threadfunc(void* arg) {
	int type;
	bool done = false;
	MyMessage* msg;
	MyMessageCommand* cmd;
	MyThread* entity = static_cast<MyThread*>(arg);

	//wait for commands from the master
	while (!done) {
		//receive a request
		msg = entity->recvMsg();
		type = msg->getType();
		//check the message type
		if (type != MSG_COMMAND_PACKET) {
			cerr << "Unexpected message: " << type << endl;
			exit(1);
		}
		//get the command message
		cmd = dynamic_cast<MyMessageCommand*>(msg);
		switch (cmd->getCommand()) {
		case MSG_COMMAND_ALLFILLUP:
			entity->allfillup();
			globalBarrier->wait();
			break;
		case MSG_COMMAND_FILTRATION:
			entity->filtration();
			globalBarrier->wait();
			break;
		case MSG_COMMAND_EXCLUDEUNIQUE:
			entity->excludeunique();
			globalBarrier->wait();
			break;
		case MSG_COMMAND_THREAD_EXIT:
			//wait for barrier
			entity->barrierwait();
			done = true;
			break;
		}
		//release the message
		delete msg;
	}

	return 0;
}
void MyThread::allfillup() {
	int type;
	bool done = false;
	MyMessage* msg;
	MyMessageCommand* cmd;
	MyMessageKmers* data;
	vector<pair<int, Kmer> >* kmers;

	/*clear the statistical values*/
	for (size_t i = 0; i < numRecvKmers.size(); ++i) {
		numRecvKmers[i] = 0;
	}
	//receive message from the master
	while (!done) {
		//receive a message
		msg = recvMsg();
		type = msg->getType();
		if (type == MSG_COMMAND_PACKET) {
			cmd = dynamic_cast<MyMessageCommand*>(msg);
			if (cmd->getCommand() == MSG_COMMAND_ALLFILLUP_DONE) {
			} else {
				cerr << "ALLFILLUP receives unexpected command "
						<< cmd->getCommand() << endl;
			}
			done = true;
		} else if (type == MSG_DATA_PACKET) {
			data = dynamic_cast<MyMessageKmers*>(msg);
			kmers = data->getKmers();
			//insert the k-mers into the bloom filter and hash table
			for (vector<pair<int, Kmer> >::iterator iter = kmers->begin();
					iter != kmers->end(); iter++) {
				if (BFs[iter->first]->contains(iter->second)) {
					kmaps[iter->first]->insert(KmerIntPair(iter->second, 0));
				} else {
					BFs[iter->first]->insert(iter->second);
				}
				numRecvKmers[iter->first]++;
			}
		} else {
			cerr << "ALLFILLUP receives unexpected message: " << type << endl;
			done = true;
		}
		//release the message
		delete msg;
	}
}
void MyThread::filtration() {
	int type;
	bool done = false;
	MyMessage* msg;
	MyMessageCommand* cmd;
	MyMessageKmers* data;
	vector<pair<int, Kmer> >* kmers;
	hmap_t::iterator miter;

	/*clear the statistical values*/
	for (size_t i = 0; i < totalRecvCov.size(); ++i) {
		totalRecvCov[i] = 0;
	}
	//receive message from the master
	while (!done) {
		//receive a message
		msg = recvMsg();
		type = msg->getType();
		if (type == MSG_COMMAND_PACKET) {
			cmd = dynamic_cast<MyMessageCommand*>(msg);
			if (cmd->getCommand() == MSG_COMMAND_FILTRATION_DONE) {
			} else {
				cerr << "FILTRATION receives unexpected command "
						<< cmd->getCommand() << endl;
			}
			done = true;
		} else if (type == MSG_DATA_PACKET) {
			data = dynamic_cast<MyMessageKmers*>(msg);
			kmers = data->getKmers();
			//insert the k-mers into the bloom filter and hash table
			for (vector<pair<int, Kmer> >::iterator iter = kmers->begin();
					iter != kmers->end(); iter++) {
				miter = kmaps[iter->first]->find(iter->second);
				if (miter != kmaps[iter->first]->end()) {
					/*accumulate the multiplicity of the k-mer*/
					miter->SetVal(miter->GetVal() + 1);

					/*calculate the total number of k-mers*/
					totalRecvCov[iter->first]++;
				}
			}
		} else {
			cerr << "FILTRATION receives unexpected message: " << type << endl;
			done = true;
		}
		//release the message
		delete msg;
	}
}
void MyThread::excludeunique() {
	int type;
	bool done = false;
	MyMessage* msg;
	MyMessageCommand* cmd;

	//delete unique k-mers
	Kmer km_del;
	km_del.set_deleted();

	/*clear the statistical values*/
	/*for each k-mer size*/
	for (size_t ikmer = 0; ikmer < kmaps.size(); ++ikmer) {
		hmap_t* kmap = kmaps[ikmer];
		size_t* freqs = &kmerFreqs[ikmer][0];
		size_t numRecvDel = 0;
		size_t maxMulti = 1;
		unsigned int multi;

		/*set the invalid key for deletion*/
		kmap->set_deleted_key(km_del);

		for (hmap_t::iterator iter = kmap->begin(); iter != kmap->end();) {
			/*get the multiplicity of this k-mer*/
			multi = iter->GetVal();
			if (multi < 2) {
				// remove k-mer that got through the bloom filter
				kmap->erase(iter++);
				numRecvDel++;
			} else {
				if (maxMulti < multi) {
					maxMulti = multi;
				}
				if (multi >= DEFAULT_MAX_MULTI) {
					freqs[DEFAULT_MAX_MULTI - 1]++;
				} else {
					freqs[multi]++;
				}

				iter++;
			}
		}
		/*save the statistical value*/
		numRecvDels[ikmer] = numRecvDel;
		maxMults[ikmer] = maxMulti;
	}
	//receive message from the master
	while (!done) {
		//receive a message
		msg = recvMsg();
		type = msg->getType();
		if (type == MSG_COMMAND_PACKET) {
			cmd = dynamic_cast<MyMessageCommand*>(msg);
			if (cmd->getCommand() == MSG_COMMAND_EXCLUDEUNIQUE_DONE) {
			} else {
				cerr << "EXCLUDEUNIQUE receives unexpected command "
						<< cmd->getCommand() << endl;
			}
			done = true;
		} else {
			cerr << "EXCLUDEUNIQUE receives unexpected message: " << type
					<< endl;
			done = true;
		}
		//release the message
		delete msg;
	}
}

size_t MyThread::getNumMsg() {
	size_t size;
	mutex.lock();
	size = msgs.size();
	mutex.unlock();

	return size;
}
void MyThread::sendMsg(MyMessage* msg) {
	while (getNumMsg() >= maxCapacity / 2)
		;

	mutex.lock();
	msgs.push_front(msg);
	mutex.unlock();
}

MyMessage*
MyThread::recvMsg() {
	MyMessage* msg;
	while (getNumMsg() == 0)
		;

	mutex.lock();
	msg = msgs.back();
	msgs.pop_back();
	mutex.unlock();

	return msg;
}

