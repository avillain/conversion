#!/mount/biotools/bin/python2.7
# -*- coding: utf-8 -*-

# Villain Adrien
# 2014     
# avillain@pasteur.fr

import os
import sys
import argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes a multiptt file from genbankconvertor.py and reformat it (suited for SynTView 2 usage)')
    parser.add_argument("ptt", help='Multiptt file')
    parser.add_argument("output", help='Ouptut file : multiptt-like suited for SynTView2 usage')

    if not sys.argv[1:] :
       sys.argv.append('-h')

    prog_name = os.path.basename(sys.argv[0])
    args = parser.parse_args()

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

    with open(args.ptt,'r') as filptt:
        relength=re.compile('^(.+) - (\d+)\.\.(\d+)$')
        renbprot=re.compile('^(\d+) proteins$')
	headline=filptt.readline()
	headr=relength.findall(headline)[0]
	start=int(headr[1])
	end=int(headr[2])
	length=end-start+1
	name=headr[0].split(' ')[0]
	protline=filptt.readline()
	nbprot=int(renbprot.findall(protline)[0])
	plusline=">%s;feats=%d;length=%d\n" %(name,nbprot,length)
	allfile=[headr[0]+" - %d.."%start,'']
	totallength=length
	totalfeats=nbprot
	localine=filptt.readline()
	allfile.append("Position\tStrand\tLength\tAccNum (key)\tGene Synonym\tLocusTag\tCOG\tProduct\n")
	allfile.append(plusline)
	for line in filptt:
	    #allfile.append(line)
	    if line[:8]=='Location':
		allfile.append(plusline)
		totallength+=length
		totalfeats+=nbprot
	    elif relength.findall(line):
	    	linum=relength.findall(line)[0]
		start=int(linum[1])
		end=int(linum[2])
		length=end-start+1
		name=linum[0].split(' ')[0]
	    elif renbprot.findall(line):
		nbprot=int(renbprot.findall(line)[0])
		plusline=">%s;feats=%d;length=%d\n" %(name,nbprot,length)
	    else:
		fields=line.split('\t')
		line="\t".join([fields[0],fields[1],fields[2],fields[5],fields[4],fields[6],fields[5],fields[7],fields[8]])
		allfile.append(line)
	allfile[0]+="%d\n"%totallength
	allfile[1]="%d proteins\n"%totalfeats
	with open(args.output,'w') as filout:
	    for i in allfile:
		filout.write(i)

	

