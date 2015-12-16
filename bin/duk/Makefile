cc = g++
SOURCES = duk.cpp  kmercoder.cpp  kmerhash.cpp seqfiles.cpp parseOpts.cpp poisson.cpp
HEADS =  kmercoder.h  kmerhash.h seqfiles.h 
PROGS = duk 

OBJECTS = $(SOURCES:.cpp=.o)

all: $(OBJECTS)
	$(cc) -o $(PROGS) $(OBJECTS)

%.o: %.cpp %.h
	$(cc) -c $<

clean:
	rm -f  $(OBJECTS)

