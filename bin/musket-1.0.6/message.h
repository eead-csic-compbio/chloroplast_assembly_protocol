#ifndef _MY_MESSAGE_H
#define _MY_MESSAGE_H

#include <vector>
#include "kmer.h"

using namespace std;
//message type

enum {
	//COMMAND packet
	MSG_COMMAND_PACKET,
	//DATA packet
	MSG_DATA_PACKET
};

enum {
	MSG_COMMAND_THREAD_EXIT, //thread exit command
	MSG_COMMAND_ALLFILLUP, //fill up all kmers command
	MSG_COMMAND_ALLFILLUP_DONE,
	MSG_COMMAND_FILTRATION,
	MSG_COMMAND_FILTRATION_DONE,
	MSG_COMMAND_EXCLUDEUNIQUE,
	MSG_COMMAND_EXCLUDEUNIQUE_DONE
};

class MyMessage {
public:
	MyMessage(int type);
	virtual
	~MyMessage() = 0;

	inline int getType() {
		return _type;
	}

protected:
	int _type;
};

class MyMessageCommand: public MyMessage {
public:
	MyMessageCommand(int command);
	~MyMessageCommand();

	inline int getCommand() {
		return _command;
	}

private:
	int _command;
};

class MyMessageKmers: public MyMessage {
public:
	MyMessageKmers(vector<pair<int, Kmer> >* kmers);
	~MyMessageKmers();

	inline vector<pair<int, Kmer> >*
	getKmers() {
		return _kmers;
	}
	;
private:
	vector<pair<int, Kmer> >* _kmers;
};
#endif

