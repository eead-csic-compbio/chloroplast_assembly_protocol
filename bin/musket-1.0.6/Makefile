CXX = g++
MACROS = -DMAX_SEQ_LENGTH=200 -DMAX_KMER_SIZE=28 #-DDEBUG
CXXFLAGS = -O3 -Wall -I. -fopenmp -funroll-loops $(MACROS)
LDLIBS  = -lm -lz -lpthread

OBJECTS =  kmer.o hash.o bloom_filter.o fastq.o thread.o message.o kmerintpair.o kseq.o voting.o paracorrect.o option.o

target=musket
all: main.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) main.o -o $(target) $(LDLIBS)
	strip $(target)

main.o: main.cpp
thread.o: thread.cpp message.h bloom_filter.h kmer.h fastq.h hash.h
message.o: message.cpp
kmerintpair.o: kmerintpair.cpp kmerintpair.h kmer.h
fastq.o: fastq.h fastq.cpp  kseq.h kseq.cpp
kseq.o: kseq.h kseq.cpp
kmer.o: kmer.h kmer.cpp
bloom_filter.o: bloom_filter.h bloom_filter.cpp
hash.o: hash.h hash.cpp
option.o:option.h option.cpp
voting.o: voting.h voting.cpp
paracorrect.o: paracorrect.h paracorrect.cpp

clean:
	rm -rf *.o
	rm -rf $(target)
