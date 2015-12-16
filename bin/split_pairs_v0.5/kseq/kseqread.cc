// Bruno Contreras-Moreira, Carlos P Cantalapiedra EEAD-CSIC 2013-14

#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "kseq.h"
#include <string>
#include <vector>
#include <algorithm>
#include <pcrecpp.h> //http://www.gammon.com.au/pcre/pcrecpp.html

using namespace std;

// declare the type of file handler and the read() function  
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
    gzFile fp=NULL;
    FILE *ifp=NULL,*fp1=NULL,*fp2=NULL,*fps=NULL;
    kseq_t *seq,*useq;
    string regex,subspattern,pairsep = "/";
    string infile,intlvfile,pairfile1,pairfile2,singlefile;
    bool interleaved=false, shorten=false, showheader=false, oneline=false;
    bool use_previous_read=false, rmsingles=false, fromstdin=false;
    string header1,header2,name1,name2,seq1,seq2,qual1,qual2;
    unsigned pairs=0,singles=0;
    size_t slash1,slash2;
    int  length,flag;

    string version = "v0.5 ";
    string usage = "-h this message\n";
    usage += "-i input FASTA/FASTQ filename      (optional, dy default reads stdin, accepts GZIP compressed files)\n";
    usage += "-1 output pair1 filename\n";
    usage += "-2 output pair2 filename\n";
    usage += "-l output interleaved filename     (optional, overrides -1 and -2, example: -l stdout)\n";
    usage += "-s output singles filename         (optional, singles are discarded otherwise)\n";
    usage += "-R Perl regex to shorten headers   (optional, example: -R \"(^.*?)#\\w+?([12])\" )\n";
    usage += "-S header substitution pattern     (optional, requires -R, example: -S \"\\1/\\2\" )\n";
    usage += "-p char separating pair numbers    (optional, default: -p \"/\" , as in HWUSI:4:1101:3600:1982#ATCACGA/1 )\n";
    usage += "-L make one-line reads and exit    (optional, prints tab-separated FASTQ reads, overrides -l,1,-2)\n";
    usage += "-H show first header and exit      (optional)\n\n";
    usage += version;
    usage += __DATE__;
    usage += " ";
    usage += __TIME__;
    usage += "\n";
    usage += "This FASTQ reader is based on http://lh3lh3.users.sourceforge.net/parsefastq.shtml and\n";
    usage += "relies on kseq.h by Attractive Chaos <attractor@live.co.uk>\n\n";
    usage += "Authors: Bruno Contreras-Moreira, Carlos P Cantalapiedra EEAD-CSIC 2013-14\n";

    if(argc == 1)
    {	
	    fprintf(stderr,"usage:\n%s\n",usage.c_str());	
        exit(-1);	
    }

    while((flag = getopt (argc, argv,"HLhi:l:1:2:s:p:R:S:")) != -1)
    {
	    switch(flag)
        {
        	case 'p': pairsep = optarg; break;
    		case 's': singlefile = optarg; break;
            case '1': pairfile1 = optarg; break;
		    case '2': pairfile2 = optarg; break;
    		case 'l': intlvfile = optarg; break;
            case 'i': infile = optarg; break;
            case 'R': regex = optarg; break;
	    	case 'S': subspattern = optarg; break;
    		case 'H': showheader = true; break;
			case 'L': oneline = true; break;
    		case 'h': fprintf (stderr,"usage:\n%s\n",usage.c_str()); exit(-1); break;
        }
    } 

    // 1) check input arguments and create outfiles if required
    if(!infile.empty())
    {
        ifp = fopen(infile.c_str(),"r");
        if(ifp == NULL)
        {
            fprintf(stderr, "# ERROR: cannot open input file %s\n",infile.c_str());
            exit(-2);
        }
    }
    else
    {
        ifp = stdin;
        fromstdin = true;
    }   
    
    fp = gzdopen(fileno(ifp),"r"); //http://biowize.wordpress.com/2013/03/05/using-kseq-h-with-stdin
 
    pcrecpp::RE cregex(regex);
    if( (!regex.empty() && subspattern.empty()) || (regex.empty() && !subspattern.empty()))
    {
	    fprintf(stderr, "# ERROR: need both -R regular expression and -S substitution pattern to shorten headers\n");
        exit(-3);  	
    } 
    else if(!regex.empty() && !subspattern.empty())
    {
	    if(cregex.error().length() > 0)
	    {
        	fprintf(stderr, "# ERROR: cannot accept regular expression, please check syntax: %s\n",
			cregex.error().c_str());
        	exit(-4);
       	}
	    else shorten = true;
    }
	  
    if(showheader == true) // if required show first header and check regex and pairsep; then exit
    {
	    useq = kseq_init(fp);
        while( (length=kseq_read(useq)) >= 0 )
        {
        		if(useq->qual.l) qual1 = useq->qual.s;
        		else qual1 = ""; 
				
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

            printf("# first header: %s\n",name1.c_str());
				
			if(useq->qual.l)
			{
				vector<unsigned int> qualities;
				float median_qual;
				for(unsigned int p=0;p<useq->qual.l;p++)
				{
					qualities.push_back(useq->qual.s[p]);
				}
				
				size_t size = qualities.size();
				sort(qualities.begin(), qualities.end());
				if (size  % 2 == 0) median_qual = (qualities[size / 2 - 1] + qualities[size / 2]) / 2;
				else median_qual = qualities[size / 2];
				
				printf("# first quality encoding: %d , %1.0f , %d         [Sanger range: 33,73]\n",
					qualities[0],median_qual,qualities[size-1]);
			}
			
			if(shorten == true)
    	    {
    		    header1 = name1;	
    		    if(cregex.Extract(subspattern.c_str(),header1.c_str(),&name1) == false)
    		    {
    			    fprintf(stderr, "# WARNING: regex does not match header: %s : %s\n",
					    header1.c_str(),regex.c_str());
    		    } //else printf("mira %s\n",name1.c_str());
                
                slash1 = subspattern.find_last_of(pairsep);
                if(slash1 == string::npos)
                {
                    fprintf(stderr, "# WARNING: -S pattern does not contain char separating pairs '%s': %s\n",
                        pairsep.c_str(),subspattern.c_str());
                }
          
                break;
    	    }
            else
            {
                slash1 = name1.find_last_of(pairsep);
                if(slash1 == string::npos)
                {
                    fprintf(stderr, "# WARNING: header does not contain char separating pairs '%s': %s\n",
                	    pairsep.c_str(),name1.c_str());
                }
		        break;
            }    
        }
        kseq_destroy(useq);
        if(fromstdin == false) gzclose(fp);
	    exit(0);
    }
    else if(oneline == true)
    {
        int tr_binary_OK = system("which tr > /dev/null");
        if(tr_binary_OK > 0)
        {
            fprintf(stderr, "# ERROR: cannot find required shell binary 'tr'\n");
            exit(-3);
        }

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

        	if(useq->qual.l) printf("@%s\t%s\t+\t%s\n",name1.c_str(),useq->seq.s,useq->qual.s);
        	else printf(">%s\t%s\n",name1.c_str(),useq->seq.s);
	    }
        kseq_destroy(useq);
        if(fromstdin == false) gzclose(fp);
	    exit(0);
    }
    else seq = kseq_init(fp);

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
            exit(-4);
        }
    }    
    
    if(!intlvfile.empty())
    {
	    interleaved = true;
        if(intlvfile == "stdout") fp1 = stdout;
        else fp1=fopen(intlvfile.c_str(),"w");
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
	 
    // 2) proceed with split
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
    			fprintf(stderr, "# WARNING: regex does not match header: %s : %s\n",
					header1.c_str(),regex.c_str());
               exit(-12);
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
   
            if(shorten == true)
            {
                header2 = name2;
                if(cregex.Extract(subspattern.c_str(),header2.c_str(),&name2) == false)
                {
                    fprintf(stderr, "# WARNING: regex does not match header: %s : %s\n",
                        header2.c_str(),regex.c_str());
                    exit(-13);
                }
            }

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
                    //printf(">%s|\n<%s|\n\n",name1.c_str(),name2.c_str());
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
            exit(-14);
    	}
    }
	
    kseq_destroy(seq);

    if(fromstdin == false) gzclose(fp);
    fclose(fp1);    
    if(rmsingles==false) fclose(fps);
    if(!interleaved) fclose(fp2);

    fprintf(stderr,"# total pairs of reads = %d total single reads = %d\n",pairs,singles);
 
    exit(0); 
}

