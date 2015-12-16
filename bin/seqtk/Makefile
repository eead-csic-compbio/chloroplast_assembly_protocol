CFLAGS=-g -Wall -O3 -DFULLHEADER

all: seqtk

seqtk: seqtk.c khash.h kseq.h
		gcc $(CFLAGS) seqtk.c -o $@ -lz -lm
		
clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk 
