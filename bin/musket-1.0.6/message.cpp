#include "message.h"

MyMessage::MyMessage(int type) {
	_type = type;
}
MyMessage::~MyMessage() {
}

MyMessageCommand::MyMessageCommand(int command) :
		MyMessage(MSG_COMMAND_PACKET) {
	_command = command;
}
MyMessageCommand::~MyMessageCommand() {
}

MyMessageKmers::MyMessageKmers(vector<pair<int, Kmer> >* kmers) :
		MyMessage(MSG_DATA_PACKET) {
	_kmers = kmers;
}
MyMessageKmers::~MyMessageKmers() {
	_kmers->clear();
	delete _kmers;
}
