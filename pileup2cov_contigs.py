#!/usr/bin/env python2.7
# -*- coding: utf-8 -*- 

# Adrien Villain - PFBAG - 12/02/14 - avillain@pasteur.fr

import os
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts a .mipleup (samtools) generated coverage file to a coverage by windows of length n file (suited for SynTView usage)')
    parser.add_argument("mpileup", help='Coverage input file in .mpileup format (see http://samtools.sourceforge.net/)')
    parser.add_argument("fasta", help='Fasta file of the reference genome (needed to get the true length of genome as .mpileup file doesn\'t include 0 coverage positions')
    parser.add_argument("ptt", help='Protein table file (ptt-like SynTView2 format), needed to get the order of contigs right')
    parser.add_argument("output", help='Ouptut file, recommended extension is .coverage. Tabulated format of every line is : start_of_window mean_coverage')
    parser.add_argument("-l", dest="length",
        help="Length n of windows to calculate mean coverage from", default=250)

    if not sys.argv[1:] :
       sys.argv.append('-h')

    prog_name = os.path.basename(sys.argv[0])
    args = parser.parse_args()

    if not args.mpileup :
	parser.error("\nMpileup file missing.")
    else :
	if not os.access(args.mpileup, os.R_OK):
	    parser.error("Unable to read mpileup file %s " %args.mpileup)

    if not args.fasta :
        parser.error("\nFasta file missing.")
    else :
        if not os.access(args.fasta, os.R_OK):
            parser.error("Unable to read fasta file %s " %args.fasta)

    if not args.ptt :
        parser.error("\nPtt file missing.")
    else :
        if not os.access(args.ptt, os.R_OK):
            parser.error("Unable to read ptt file %s " %args.ptt)

    if not args.output :
        parser.error("\nOutput file name is missing.")
    else:
	if os.access(args.output, os.F_OK) and not os.access(args.output, os.W_OK):
            parser.error("\nOutput file %s exists but is not writable." %args.output)
    
    try:
	args.length=int(args.length)
    	print "Using windows of size : %d" %args.length
    except:
	parser.error("\nLength of window %s is not castable to int." %args.length)

    # Extract sequence of genome from fasta file (in order to get seq length)
    with open(args.fasta,'r') as filfasta:
	#seqs=[]
	seqs={}
    	seq=''
    	for f in filfasta:
	    if f[:1]==">":
		if seq:
	            seqs[name]=seq
		    seq=""
		name=f[1:-1]
	    else:
	        seq+=f[:-1]
	seqs[name]=seq

    with open(args.ptt,'r') as filptt:
	contigs=[]
	for f in filptt:
		if f[:1]==">":
		    contigs.append(f[1:-1].split(';')[0])
    #
    from pprint import pprint
    #print len(seqs['NODE_43_length_71884_cov_79.293930-1stcontig'])
    #pprint(seqs.keys())
    #exit(0)

    # Read mpileup file and write to output file
    coverages={}
    done=[]
    with open(args.output,'w') as filout:
	with open(args.mpileup,'r') as filin:
	    for k in range(len(seqs.keys())):
		coverage=""
	    	cov=0 # mean coverage
	    	n=0 # number of windows
		if k==0:
	    	    fields=filin.readline().split('\t')
	    	name=fields[0]
	    	if not name:
			break
		#filout.write(">%s\n"%name)
		#filout.write("Position\tTotReadsSNP\n")
		coverage+=">%s\n"%name
		coverage+="Position\tTotReadsSNP\n"
	        for i in range(len(seqs[name])): # for the whole contig
		    try:
		        if i+1==int(fields[1]): # if the ith position is in the mpileup file (ie not 0 coverage)
		            cov+=int(fields[3]) # add coverage
     		    	    fields=filin.readline().split('\t') # read next line
		    except:
		        pass # nothing to do since coverage is 0
		    if (i+1)%args.length==0 and i!=0: # if we are at the end of the window
		        #filout.write("%d\t%.2f\n"%(args.length*n+1,(cov+0.0)/args.length)) # write the last window 
		        coverage+="%d\t%.2f\n"%(args.length*n+1,(cov+0.0)/args.length)
			cov=0 # reset mean coverage
		        n+=1 # 1 more window is done
		#if name=='NODE_43_length_71884_cov_79.293930-1stcontig':
			#print n, cov, len(seqs[name]), len(seqs[name])-n*args.length
	        #filout.write("%d\t%.2f\n"%(args.length*n+1,(cov+0.0)/(len(seqs[name])-n*args.length))) # last window
		coverage+="%d\t%.2f\n"%(args.length*n+1,(cov+0.0)/(len(seqs[name])-n*args.length))
		coverages[name]=coverage
		done.append(name)
	
	if len(done)!=len(contigs):
	    for name in contigs:
		if name not in done:
		    coverage=">%s\nPosition\tTotReadsSNP\n" %name
		    n=0
		    for i in range(len(seqs[name])):
			if (i+1)%args.length==0 and i!=0:
			    coverage+="%d\t%.2f\n"%(args.length*n+1,0.00)
			    n+=1
		    coverage+="%d\t%.2f\n"%(args.length*n+1,0.00)
		    coverages[name]=coverage
		    done.append(name)
	
	for i in contigs:
	    filout.write(coverages[i])
