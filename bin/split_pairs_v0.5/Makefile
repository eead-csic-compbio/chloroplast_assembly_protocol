CFLAGS=-g -Wall -O3 -DFULLHEADER 

all: split_pairs_bf kseqread pigz

split_pairs_bf:
	cd kseq; g++ $(CFLAGS) split_pairs_bf.cc -o $@ -lz -lpcrecpp; mv split_pairs_bf ..

kseqread: 
	cd kseq; g++ $(CFLAGS) kseqread.cc -o $@ -lz -lpcrecpp; mv kseqread ..		

pigz: 
	cd pigz-2.3.1; make pigz; rm -f unpigz; mv pigz ..

clean:
	rm -f split_pairs_bf kseqread pigz; cd pigz-2.3.1; make clean
