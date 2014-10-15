#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-


# prism2vcf.py
#
__desc__ = """Converting PRISM output to VCF"""
__author__ = "Adrien Villain"

import argparse

def parse_fasta(infile):
	with open(infile,'r') as fastafile:
		dic={}
		seq=''
		id=''
		for line in fastafile:
			if line[:1]=='>':
				if seq!='':
					dic[id]=seq
					id=line[1:-1]
					seq=''
				else:
					id=line[1:-1]
			else:
				seq+=line[:-1]
		dic[id]=seq
	return(dic)

def parse_prism(infile, outfile, fastdic):
	with open(outfile,'w') as out:
		with open(infile,'r') as inf:
			for line in inf:
				if line[:1]!="@":
					fields=line.split('\t')
					chr=fields[0]
					sta=fields[1]
					end=fields[2]
					siz=int(end)-int(sta)+1
					typ=fields[4]
					cov=fields[5]
					ide="."
					qua="."
					fil="PASS"
					inf="DP=%s;SVTYPE=%s;SVLEN=%s"%(cov,typ,siz)
					if typ=="INS":
						ref="."
						alt=fields[21].strip()
					else:
						ref=fastdic[chr][int(sta)-1:int(end)]
						alt="."
						#print ref
					outline="\t".join([chr,sta,ide,ref,alt,qua,fil,inf]) + "\n"
					out.write(outline)
					#print outline

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('prism', help='Prism input file')
	parser.add_argument('fasta', help='Reference genome fasta file')
	parser.add_argument('output', help='Output vcf')
	args = parser.parse_args()
	fas=parse_fasta(args.fasta)
	#import pprint
	#pprint.pprint(fas)
	parse_prism(args.prism, args.output, fas)

