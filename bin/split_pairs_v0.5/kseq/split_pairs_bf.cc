// Bruno Contreras-Moreira, Carlos P Cantalapiedra, MaJesus Garcia Pereira EEAD-CSIC 2013-14

#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "kseq.h"
#include <string>
#include <vector>
#include <pcrecpp.h> //http://www.gammon.com.au/pcre/pcrecpp.html

using namespace std;

// declare the type of file handler and the read() function  
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
    gzFile fp=NULL,tmpsfp=NULL;
    FILE *fp1=NULL,*fp2=NULL,*fps=NULL;
    kseq_t *seq,*useq;
    string regex,subspattern,pairsep = "/";
    string infile,intlvfile,pairfile1,pairfile2,singlefile;
    bool nosort=false, interleaved=false, parallel=false;
    bool shorten=false, showheader = false, oneline=false;
    string header1,header2,name1,name2,seq1,seq2,qual1,qual2;
    unsigned pairs=0,singles=0,threads=0,sortcache = 1024; //Mb, applies to default serial sort
    char strcache[21]; sprintf(strcache,"%dM ",sortcache);
    size_t slash1,slash2;
    int  length,flag;
    bool use_previous_read = false, rmsingles = false;
    string tmpunsorted = "/unsortedXXXXXX";
    string tmpsorted = "/sortedXXXXXX";
    string pathtotmp = "/tmp";
    struct stat st;

    string version = "v0.5 ";
    string usage = "-h this message\n";
    usage += "-i input FASTA/FASTQ filename      (accepts GZIP compressed files)\n";
    usage += "-1 output pair1 filename\n";
    usage += "-2 output pair2 filename\n";
	usage += "-l output interleaved filename     (optional, overrides -1 and -2)\n";
    usage += "-s output singles filename         (optional, singles are discarded otherwise)\n";
    usage += "-R Perl regex to shorten headers   (optional, example: -R \"(^.*?)#\\w+?([12])\" )\n";
    usage += "-S header substitution pattern     (optional, requires -R, example: -S \"\\1/\\2\" )\n";
    usage += "-p char separating pair numbers    (optional, default: -p \"/\" , as in HWUSI:4:1101:3600:1982#ATCACGA/1 )\n";
    usage += "-n do not sort input reads         (optional, by default reads are sorted with shell sort)\n";
    // disabled, see comments below
    //usage += "-P sort with GNU/parallel          (optional, by default shell sort is used)\n";
    //usage += "-t max threads for GNU/Parallel    (optional, default: all available CPU cores, requires -P)\n";
    usage += "-H show first header and exit      (optional)\n";
    usage += "-L make one-line reads and exit    (optional, creates file actually used for sorting reads)\n";
    usage += "-B RAM buffer (Mb) while sorting   (optional, default: -B ";
    usage += strcache;
    usage += ")\n";
    usage += "-T path for temporary files        (optional, default: -T ";
    usage += pathtotmp;
    usage += " )\n\n";
    usage += version;
    usage += __DATE__;
    usage += " ";
    usage += __TIME__;
    usage += "\n";
    usage += "This software is based on http://lh3lh3.users.sourceforge.net/parsefastq.shtml and\n";
    usage += "relies on kseq.h by Attractive Chaos <attractor@live.co.uk>\n\n";
    usage += "Authors: Bruno Contreras-Moreira, Carlos P Cantalapiedra, MaJesus Garcia Pereira EEAD-CSIC 2013-14\n";

    if(argc == 1)
    {	
	    fprintf(stderr,"usage:\n%s\n",usage.c_str());	
        exit(-1);	
    }

    while((flag = getopt (argc, argv,"PLHnhi:l:1:2:s:p:R:S:t:T:B:")) != -1)
    {
	    switch(flag)
        {
        	case 'p': pairsep = optarg; break;
		    case 'P': parallel  = true; break;
    		case 't': threads = atoi(optarg); break;
           	case 'n': nosort  = true; break;
    		case 's': singlefile = optarg; break;
            case '1': pairfile1 = optarg; break;
		    case '2': pairfile2 = optarg; break;
    		case 'l': intlvfile = optarg; break;
            case 'i': infile = optarg; break;
            case 'R': regex = optarg; break;
	    	case 'S': subspattern = optarg; break;
    		case 'H': showheader = true; break;
    		case 'L': oneline = true; break;
    		case 'T': pathtotmp = optarg; break;
            case 'B': sortcache = atoi(optarg); break;
    		case 'h': fprintf (stderr,"usage:\n%s\n",usage.c_str()); exit(-1); break;
      	}
    } //printf("t: %d\n",threads); // change to show any args desired

    // 1) check input and create outfiles
    fp = gzopen(infile.c_str(),"r");
    if(fp == NULL)
    {
        fprintf(stderr, "# ERROR: cannot open input file %s\n",infile.c_str());
        exit(-2);
    }

    if(showheader == true) // show first header if required, and exit
    {
	    useq = kseq_init(fp);
        while( (length=kseq_read(useq)) >= 0 )
        {
		    printf("# first header:\n%s %s\n",useq->name.s,useq->comment.s);
		    break;
        }
        kseq_destroy(useq);
        gzclose(fp);
	    exit(0);
    }

    if(singlefile.empty())
    {
        rmsingles = true;
    }
    else
    {
        fps=fopen(singlefile.c_str(),"w");
        if(fps == NULL)
        {
    	    fprintf(stderr, "# ERROR: cannot create singles file %s\n",singlefile.c_str());
            exit(-3);
        }
    }    
    
    if(!intlvfile.empty())
    {
	    interleaved = true;
    	fp1=fopen(intlvfile.c_str(),"w");
	    if(fp1 == NULL)
        {
        	fprintf(stderr, "# ERROR: cannot create interleaved outfile %s\n",intlvfile.c_str());
        	exit(-5);
    	}
	    fp2=fp1;
    }    
    else
    {
        if(pairfile1.empty())
        {
		    fprintf(stderr, "# ERROR: need a valid pairs1 file name\n");
            exit(-6);
	    }
    	fp1=fopen(pairfile1.c_str(),"w");
    	if(fp1 == NULL)
    	{
        	fprintf(stderr, "# ERROR: cannot create pairs1 file %s\n",pairfile1.c_str());
        	exit(-7);
    	}

	    if(pairfile2.empty())
        {
                fprintf(stderr, "# ERROR: need a valid pairs2 file name\n");
                exit(-8);
        }
	    fp2=fopen(pairfile2.c_str(),"w");
	    if(fp2 == NULL)
    	{
        	fprintf(stderr, "# ERROR: cannot create pairs2 file %s\n",pairfile2.c_str());
        	exit(-9);
    	}
    }

    // 2) check regex
    pcrecpp::RE cregex(regex);
    if( (!regex.empty() && subspattern.empty()) || 
	(regex.empty() && !subspattern.empty()) )
    {
	    fprintf(stderr, "# ERROR: need both -R regular expression and -S substitution pattern to shorten headers\n");
        exit(-10);  	
    } 
    else if(!regex.empty() && !subspattern.empty())
    {
	    if(cregex.error().length() > 0)
	    {
        	fprintf(stderr, "# ERROR: cannot accept regular expression, please check syntax: %s\n",
			cregex.error().c_str());
        	exit(-11);
       	}
	    else shorten = true;
    }

    // 3) sort reads with calls to shell if required (this way we can sort them even if they don't fit in RAM)
    if(nosort == false)
    {
        int sort_binary_OK = system("which sort > /dev/null");
        int tr_binary_OK = system("which tr > /dev/null");
        if(sort_binary_OK > 0 || tr_binary_OK > 0)
        {
            fprintf(stderr, "# ERROR: cannot find required shell binaries 'sort' and 'tr'\n");
            exit(-12);
        }

	    if (stat(pathtotmp.c_str(), &st) == 0 && S_ISDIR(st.st_mode))
    	{
        	tmpunsorted = pathtotmp + tmpunsorted;
	        tmpsorted = pathtotmp + tmpsorted; 
    	}
    	else
    	{
        	fprintf(stderr, "# ERROR: %s is not a valid directory for temporary files\n",pathtotmp.c_str());
        	exit(-13);
    	}

        if(sortcache < 10)
        {
            fprintf(stderr, "# ERROR: %d is not a valid buffer size for sorting\n",sortcache);
            exit(-14);
        }
        else sprintf(strcache,"%dM ",sortcache);

    	string command;
    	char *tmpunsortedfilename = &tmpunsorted[0];
    	char *tmpsortedfilename = &tmpsorted[0];
    	int tmpsfd = mkstemp(tmpsortedfilename); close(tmpsfd); 
    	int tmpusfd = mkstemp(tmpunsortedfilename); close(tmpusfd);
    	FILE *tmpfup = fopen(tmpunsortedfilename,"w");
        if(tmpfup == NULL)
        {
                fprintf(stderr, "# ERROR: cannot create temporary file %s\n",tmpunsortedfilename);
                exit(-15);
        }//else fprintf(stderr, "# temporary files %s , %s\n",tmpunsortedfilename,tmpsortedfilename);

	    useq = kseq_init(fp);
    	while( (length=kseq_read(useq))!=-1 )
    	{
            if(length<-1 || useq->name.l == 0)
            {
                fprintf(stderr,"# WARNING: skip bad sequence: '%s' (%d)\n",useq->name.s,length);
                continue;
            }

		    #ifdef FULLHEADER
            name1 = useq->name.s;
		    if(useq->comment.l)
        	{
                	name1 += " ";
	                name1 += useq->comment.s;
		    }
            #else
            name1 = useq->name.s;
            #endif

		    if(shorten == true) 
		    {
	            header1 = name1;
			    if(cregex.Extract(subspattern.c_str(),header1.c_str(),&name1) == false)
        		{
                		fprintf(stderr, "# WARNING: regex does not match header: %s\n",header1.c_str());
                		exit(-16);
			    }
        	}

        	if(useq->qual.l) fprintf(tmpfup,"@%s\t%s\t+\t%s\n",name1.c_str(),useq->seq.s,useq->qual.s);
        	else fprintf(tmpfup,">%s\t%s\n",name1.c_str(),useq->seq.s);
	    }
    	kseq_destroy(useq);
	    fclose(tmpfup);

        if(oneline == true)
	    {
		    printf("# output file with one-line reads: %s\n",tmpunsortedfilename);
		    exit(0);
	    }

    	//make sure sort handles correctly punctuation signs
	    setenv( "LC_ALL", "POSIX", 1 );
    	if(parallel == true) // deprecated until we find a way to make it work reliably 2013
	    {
		    //http://www.gnu.org/software/parallel/man.html
    		/*14052013 en nuestras pruebas ordena bien con -N 1000000, pero es ~2x mas lento!
    		time cat /tmp/unsortedXXiMri1w | parallel -N 1000000 --pipe --files sort -S 500M -k 1,1 -u | \
    			parallel -Xj1 sort -S 500M -k 1,1 -u -m {} ';' rm {} > XXiMri1wPall
    		real	7m3.143s
    		user	2m59.357s
    		sys	3m58.683s
    		time sort -k 1,1 -u -S 500M /tmp/unsortedXXiMri1w > XXiMri1w1C
    		real	2m55.332s
    		user	1m38.425s
    		sys	1m3.667s 
    		diff XXiMri1w1C XXiMri1wPall > diffXXiMri1w
    		-rw-rw-r--. 1 contrera contrera 12060529132 May 14 14:43 XXiMri1w1C
    		-rw-rw-r--. 1 contrera contrera 12060529132 May 14 11:38 XXiMri1wPall
    		-rw-rw-r--. 1 contrera contrera           0 May 14 11:48 diffXXiMri1w
    
    		There is 1 pair of difference!!
                    ./split_pairs -i chloroplast.2084.4.1756.fastq.gz -R "(^.*?)#\w+?(/[12])" -S "\1\2" -l intlv.fastq -s singles.fastq 
    		# total pairs of reads = 17549424 total single reads = 473972
    		./split_pairs -i chloroplast.2084.4.1756.fastq.gz -R "(^.*?)#\w+?(/[12])" -S "\1\2" -l intlv.fastq -s singles.fastq -P 
    		# total pairs of reads = 17549423 total single reads = 473974
    		*/
				
	    	command = "cat ";
            command += tmpunsorted;
    		if(threads > 0)
    		{
              	command += "| parallel -N 1000000 -P "; //subjobs of 10E6 seqs, to avoid creating too many for merging
    			char strthreads[21]; sprintf(strthreads,"%d",threads);
    			command += strthreads;
    			command += " --pipe --files sort -t'\t' -k 1,1 ";
    			command += "| parallel -P ";
    			command += strthreads;
    			command += " -Xj1 sort -t'\t' -k 1,1 -u -m {} ';' rm {} "; // merge all subsorted files in one go
    		}
    		else 
    		{
	    		command += "| parallel --pipe --files sort -t'\t' -k 1,1 ";
		    	command += "| parallel -Xj1 sort -t'\t' -k 1,1 -u -m {} ';' rm {} "; 
    		}
    		command += "| tr \"\t\" \"\n\" > ";
            command += tmpsorted;
    	}
    	else 
    	{
    		command = "sort -t'\t' -k 1,1 -u -S ";
            command += strcache;
    		command += "-T ";
    		command += pathtotmp; 
    		command += " ";
    		command += tmpunsorted;
    		command += "| tr \"\t\" \"\n\" > ";
	    	command += tmpsorted;
    	}
    
    	int sorttroutcome = system(command.c_str()); //printf("# sort: %s\n",command.c_str()); exit(-20);
	    if(sorttroutcome > 0)
        {
            fprintf(stderr, "# ERROR: could not execute shell sort\n");
            exit(-17);    
        } 

    	stat(tmpsortedfilename, &st);
    	tmpsfp = gzopen(tmpsortedfilename,"r");
    	if(tmpsfp == NULL)
    	{
            fprintf(stderr, "# ERROR: cannot open temporary file %s, sort probably failed!\n",tmpsortedfilename);
            exit(-18);
    	}	
    	else if(st.st_size == 0)
    	{
    		fprintf(stderr, "# ERROR: temporary file %s is empty, sort probably failed!\n",tmpsortedfilename);
    		exit(-19);
    	}
	    else
    	{
    		seq = kseq_init(tmpsfp);
    		unlink(tmpunsortedfilename);
    	}
    
    	shorten = false; // avoid shortening twice
    }
    else seq = kseq_init(fp);	

    //4) proceed with split
    while( (length=kseq_read(seq)) ) 
    {
        if(length == -1) // empty stream
        {
            if(use_previous_read == true)
            {
                if(rmsingles == false)
                {
                    if(qual2.size()>0) fprintf(fps,"@%s\n%s\n+\n%s\n",name2.c_str(),seq2.c_str(),qual2.c_str());
                    else fprintf(fps,">%s\n%s\n",name2.c_str(),seq2.c_str());
                }    
                singles++;
            }
            break;
	    }

        if(length<-1 || seq->name.l==0)
        {
            fprintf(stderr,"# WARNING: skip bad sequence: '%s' (%d)\n",seq->name.s,length);
            continue;
        } 

    	#ifdef FULLHEADER
        name1 = seq->name.s;
    	if(seq->comment.l)
    	{
          	name1 += " ";
    		name1 += seq->comment.s;
    	}
        #else
        name1 = seq->name.s;
        #endif
    	
    	//printf("antes: %s\n",name1.c_str());
    	if(shorten == true)
    	{
    		header1 = name1;	
    		if(cregex.Extract(subspattern.c_str(),header1.c_str(),&name1) == false)
    		{
    			fprintf(stderr, "# WARNING: regex does not match header: %s\n",header1.c_str());
               	exit(-20);
    		}//printf("ahora: %s\n",name1.c_str());
    	}
    	seq1  = seq->seq.s;
        if (seq->qual.l) qual1 = seq->qual.s;
        else qual1 = ""; 
              
        if(use_previous_read == false) //take next sequence from stream
        {
            length=kseq_read(seq);
            if(length == -1)
            {   
                if(rmsingles == false)
                {
                    if(qual2.size()>0) fprintf(fps,"@%s\n%s\n+\n%s\n",name1.c_str(),seq1.c_str(),qual1.c_str());
                    else fprintf(fps,">%s\n%s\n",name1.c_str(),seq1.c_str());
                }    
                singles++;
                break;
            }

            if(length<-1 || seq->name.l==0)
            { 
                fprintf(stderr,"# WARNING: skip bad sequence: '%s' (%d)\n",seq->name.s,length);
                continue;
            }
    
            #ifdef FULLHEADER
            name2 = seq->name.s;
    	    if(seq->comment.l)
            {
                name2 += " ";
                name2 += seq->comment.s;
            }
            #else
            name2 = seq->name.s;
            #endif
    
            seq2  = seq->seq.s;
            if(seq->qual.l) qual2 = seq->qual.s;
            else qual2 = "";
        }    
        else // reuse previous sequence
        {
            name1.swap(name2);
            seq1.swap(seq2);
            qual1.swap(qual2);
        }
    
        slash1 = name1.find_last_of(pairsep);
        if(slash1 != string::npos) 
        {
            slash2 = name2.find_last_of(pairsep);
            if(slash2 != string::npos) 
            {
                if(name1.compare(0,slash1,name2,0,slash2) == 0) // names are identical, pair found
                {                        
                    //printf("> %s\n> %s\n\n",name1.c_str(),name2.c_str());
                    use_previous_read = false;
                    if(qual1.size()>0)
                        fprintf(fp1,"@%s\n%s\n+\n%s\n",name1.c_str(),seq1.c_str(),qual1.c_str());
                    else fprintf(fp1,">%s\n%s\n",name1.c_str(),seq1.c_str());
        
                    if(qual2.size()>0)
                         fprintf(fp2,"@%s\n%s\n+\n%s\n",name2.c_str(),seq2.c_str(),qual2.c_str());
                    else fprintf(fp2,">%s\n%s\n",name2.c_str(),seq2.c_str());
                    pairs++;
        	    }
                else // not a pair
                {   
                    //printf("> %s\n< %s\n\n",name1.c_str(),name2.c_str());
                    use_previous_read = true;
                    singles++;
                 
                    if(rmsingles == false)
                    {
                        if(qual1.size()>0)
                            fprintf(fps,"@%s\n%s\n+\n%s\n",name1.c_str(),seq1.c_str(),qual1.c_str());
                        else fprintf(fps,">%s\n%s\n",name1.c_str(),seq1.c_str());	                    
                    }
                }
            }
        }
	    else
    	{
    		fprintf(stderr, "# WARNING: header does not contain char separating pairs '%s': %s\n",
                    pairsep.c_str(),name1.c_str());
            exit(-21);
    	}
    }
	
    kseq_destroy(seq);

    gzclose(fp);
    if(nosort==false) gzclose(tmpsfp);
    fclose(fp1);    
    if(rmsingles==false) fclose(fps);
    if(!interleaved) fclose(fp2);

    if(nosort==false) unlink(tmpsorted.c_str()); //delete last temp file
    	
    fprintf(stderr,"# total pairs of reads = %d total single reads = %d\n",pairs,singles);
 
    exit(0); 
}

