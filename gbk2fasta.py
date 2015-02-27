#!/usr/bin/env python2.7
# -*- coding: utf-8 -*- 

from Bio import SeqIO
import argparse
import sys

if __name__ == '__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('genbank', help='Genbank input file')
        parser.add_argument('fasta', nargs='?', help='Output genome fasta file')
        args = parser.parse_args()
	try:
                sys.stdout=open(args.fasta,"w") if args.fasta else sys.__stdout__
        except:
                print("Error, could not create file %s\n" %args.fasta)
		raise
	try:
		with open(args.genbank, "r") as genfile:
			for seq_record in SeqIO.parse(genfile, "genbank") :
				print >> sys.stderr, "Dealing with GenBank record %s" % seq_record.id
    				print ">%s %s\n%s\n" % (seq_record.id, seq_record.description, seq_record.seq)

	except:
		print("Error, could not open file %s\n" %args.genbank)
		raise
	if args.fasta:
		sys.stdout.close()

